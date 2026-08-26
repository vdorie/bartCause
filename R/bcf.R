## Bayesian causal forest.
##
## bcf() and the internal getBCFResponseFit (R/responseFit.R) are two front
## doors over one engine-facing core, fitBCF, which is the only place in the
## package that touches a dbarts multi-forest sampler. It builds the treatment
## basis, installs the two column masks, drives the sampler with four calls and
## packages the per-draw channels as an object of class 'bartBCF'.
##
## The prognostic forest mu splits on the confounders and the propensity score;
## the treatment forest tau splits on the moderators only and reaches the model
## through a basis of treatment indicators, whose amplitudes (b.0, b.1) are what
## make the counterfactual free. No test surface is involved -- a multi-forest
## sampler refuses one at creation -- so the counterfactual is recovered by
## swapping the amplitude the treatment forest is multiplied by.

## The forest-indexed containers are keyed by these names, in the engine's own
## order: the prognostic forest first, the treatment forest second.
bcfForestNames <- c("mu", "tau")

## Normalize a subset to positive integer positions. All four kinds reach the
## bart path, but the treatment arrives as a bare vector carrying no row names;
## a character subset has to be resolved against the frame it addresses before
## it can index that vector.
normalizeBCFSubset <- function(subset, n, rowNames)
{
  if (is.null(subset)) return(seq_len(n))

  ## an NA-bearing logical subset is refused downstream by both paths
  ## ("response contains missing values"), so it is not guarded here
  if (is.logical(subset)) return(which(subset))

  if (is.character(subset)) {
    if (is.null(rowNames))
      stop("character 'subset' requires data with row names to match against")
    matched <- match(subset, rowNames)
    if (anyNA(matched))
      stop("'subset' names rows not present in the data: ",
           paste0(subset[is.na(matched)], collapse = ", "))
    return(matched)
  }

  subset <- as.integer(subset)
  if (any(subset < 0L)) seq_len(n)[subset] else subset
}

## dbarts' run channels are predictor/observation-major with the chain last; the
## rest of bartCause is chain-major with the observation last, and drops the
## chain margin at one chain. This is dbarts:::convertSamplesFromDbartsToBart's
## uncombined branch, which is not exported.
convertBCFSamples <- function(x)
  if (length(dim(x)) > 2L) aperm(x, c(3L, 2L, 1L)) else t(x)

## One forest's slab of the run's forestFits channel, n.obs x n.forests x
## n.samples x n.chains, or n.obs x n.forests x n.samples at one chain.
sliceBCFForest <- function(x, index)
  if (length(dim(x)) > 3L) x[, index, , ] else x[, index, ]

## One row of the run's glue channel, sum(q) x n.samples x n.chains, or
## sum(q) x n.samples at one chain. Row 1 is the prognostic forest's a, rows 2
## and 3 the treatment forest's (b.0, b.1).
sliceBCFGlue <- function(x, index)
  if (length(dim(x)) > 2L) x[index, , ] else x[index, ]

## sigma arrives as n.samples x n.chains, or as a plain vector at one chain.
## It is stored as an [n.chains, n.samples] matrix at every chain count: both
## bartCause readers of a flat sigma try to un-flatten it themselves, and a
## matrix means neither reshape fires.
asBCFSigmaMatrix <- function(x, n.chains)
  if (is.null(dim(x))) matrix(x, nrow = n.chains) else t(x)

## mu under a fixed treatment condition: the observed surface where the
## condition matches the observed treatment and the counterfactual where it does
## not. The array form of the same helper extract.bartcFit builds locally.
bcfObsCfToTrtCtl <- function(obs, cf, trt) {
  if (length(dim(obs)) > 2L)
    aperm(aperm(obs, c(3L, 1L, 2L)) * trt + aperm(cf, c(3L, 1L, 2L)) * (1 - trt), c(2L, 3L, 1L))
  else
    t(t(obs) * trt + t(cf) * (1 - trt))
}

## (b.1 - b.0), one value per (chain, sample). The packaged glue is chain-major
## with the amplitude last, so as.vector() of this is exactly the leading-margin
## layout of the [n.chains, n.samples, n.obs] surfaces it multiplies.
bcfGlueDifference <- function(glue)
  if (length(dim(glue)) > 2L) glue[, , "b.1"] - glue[, , "b.0"] else glue[, "b.1"] - glue[, "b.0"]

fitBCF <- function(dbartsDataCall, evalEnv, z, treatmentName,
                   p.scoreName = NULL, moderators = NULL,
                   n.trees = 200L, base = 0.95, power = 2.0,
                   n.trees.treatment = 50L,
                   treatment.base = 0.25, treatment.power = 3,
                   sd.control = NULL, sd.moderate = 1,
                   b.prior.variance = 0.5,
                   update.a = TRUE, update.b = TRUE,
                   mu.interactions = NULL, tau.interactions = NULL,
                   tau.blocks = NULL,
                   n.samples = 500L, n.burn = 500L, n.chains = 4L,
                   n.threads = dbarts::guessNumCores(), combineChains = TRUE,
                   keepSampler = FALSE, verbose = TRUE, seed = NA_integer_,
                   call = NULL, ...)
{
  n.samples <- coerceOrError(n.samples, "integer")[1L]
  n.burn    <- coerceOrError(n.burn,    "integer")[1L]
  n.chains  <- coerceOrError(n.chains,  "integer")[1L]
  n.threads <- coerceOrError(n.threads, "integer")[1L]

  ## which interface the data call names: it decides where the row names come
  ## from, where the response is read for the missing-value refusal, and whether
  ## the formula needs the fit environment bound into it
  formulaArgument <- dbartsDataCall[[2L]]
  usesFormula <- is.call(formulaArgument) && identical(formulaArgument[[1L]], quote(`~`))

  dataValue <- if (is.null(dbartsDataCall$data)) NULL else eval(dbartsDataCall$data, evalEnv)
  rowNames <-
    if (is.data.frame(dataValue)) rownames(dataValue)
    else if (!usesFormula) rownames(eval(formulaArgument, evalEnv))
    else NULL

  ## resolve the subset one line before the data object is built, as the bart
  ## fitter does, and normalize it to row positions: it indexes the bare
  ## treatment vector here, and dbarts reads it as an ordinary row subscript
  ## when it restricts the basis, which a character subset would not survive
  subsetExpr <- dbartsDataCall$subset
  subsetValue <-
    if (is.null(subsetExpr)) NULL
    else if (is.list(dataValue)) eval(subsetExpr, dataValue, evalEnv)
    else eval(subsetExpr, evalEnv)
  subset <- normalizeBCFSubset(subsetValue, length(z), rowNames)

  z <- as.vector(z)
  if (anyNA(z) || !all(z %in% c(0, 1)))
    stop("'", treatmentName, "' must be a binary treatment coded 0 and 1")
  if (length(unique(z[subset])) < 2L)
    stop("every observation left after subsetting has '", treatmentName, "' equal to ",
         unique(z[subset])[1L],
         "; a treatment arm with no observations is an all-zero basis column, which fits ",
         "with an amplitude nothing identifies")

  ## the basis is built from the FULL-LENGTH treatment and handed over at that
  ## length: both interfaces check a forest basis against the data before
  ## 'subset' and restrict it to the kept rows themselves, and one arriving
  ## already at the kept-row count is refused. The column order is (1 - z, z),
  ## which is model.matrix(~ factor(z) - 1)'s and the order the run's glue
  ## channel stacks (a, b.0, b.1) in
  basis <- cbind(1 - z, z)

  fitEnv <- new.env(parent = evalEnv)
  fitEnv[["bcf.bases"]] <- list(NULL, basis)
  fitEnv[["bcf.subset"]] <- subset
  dbartsDataCall$bases  <- quote(bcf.bases)
  dbartsDataCall$subset <- quote(bcf.subset)
  ## model.frame re-evaluates 'subset' in the data and then in the formula's own
  ## environment, so the formula has to see the bindings installed just above
  if (usesFormula) {
    environment(formulaArgument) <- fitEnv
    dbartsDataCall[[2L]] <- formulaArgument
  }

  ## refuse a missing response ahead of dbarts, whose own refusal names the
  ## symptom rather than the reason
  responseValue <-
    if (!usesFormula) dataValue
    else if (is.list(dataValue)) tryCatch(eval(formulaArgument[[2L]], dataValue, evalEnv),
                                          error = function(e) NULL)
    else tryCatch(eval(formulaArgument[[2L]], evalEnv), error = function(e) NULL)
  if (!is.null(responseValue) && anyNA(responseValue[subset]))
    stop("a Bayesian causal forest cannot be fit with missing response values; recovering them ",
         "rides a counterfactual test surface, which a multi-forest sampler refuses at creation, ",
         "so drop the incomplete rows")

  responseData <- eval(dbartsDataCall, envir = fitEnv)
  n.obs <- nrow(responseData@x)

  ## the moderator exclusion. The treatment reaches the model through the
  ## treatment forest's basis, so a prognostic forest free to split on it could
  ## absorb the effect; the propensity score belongs to the prognostic forest
  ## alone. The score's name is the one the caller resolved, never re-derived
  ## from the design.
  designNames <- colnames(responseData@x)
  if (treatmentName %not_in% designNames)
    stop("the treatment column '", treatmentName, "' is not among the design columns (",
         paste0(designNames, collapse = ", "),
         "); a treatment the prognostic forest cannot be masked out of would absorb the effect")
  if (!is.null(p.scoreName) && p.scoreName %not_in% designNames)
    stop("the propensity score column '", p.scoreName, "' is not among the design columns (",
         paste0(designNames, collapse = ", "),
         "); it has to be maskable by name to be kept out of the treatment forest")
  muVars  <- setdiff(designNames, treatmentName)
  tauVars <- setdiff(designNames, c(treatmentName, p.scoreName))
  if (!is.null(moderators)) {
    if (!is.character(moderators))
      stop("'moderators' must be a character vector of predictor names")
    unknownModerators <- setdiff(moderators, designNames)
    if (length(unknownModerators) > 0L)
      stop("'moderators' names columns not in the design: ",
           paste0(unknownModerators, collapse = ", "))
    excludedModerators <- intersect(moderators, c(treatmentName, p.scoreName))
    if (length(excludedModerators) > 0L)
      stop("'moderators' cannot include ", paste0("'", excludedModerators, "'", collapse = " or "),
           ": the treatment enters the treatment forest through its basis and the propensity ",
           "score is a prognostic-forest column, so neither is available to split on")
    tauVars <- moderators
  }

  ## the sampler: four calls, mirroring bart2's standard path. Every per-draw
  ## channel arrives batched from the second run, so there is no per-sweep loop
  extraArgs <- list(...)
  controlFormals <- names(formals(dbarts::dbartsControl))
  controlArgs <- list(n.chains = n.chains, n.threads = n.threads,
                      n.trees = coerceOrError(n.trees, "integer")[1L],
                      n.burn = n.burn, n.samples = n.samples,
                      verbose = as.logical(verbose)[1L],
                      updateState = FALSE, seed = coerceOrError(seed, "integer")[1L])
  extraControl <- extraArgs[names(extraArgs) %in% controlFormals &
                            names(extraArgs) %not_in% names(controlArgs)]
  control <- do.call(dbarts::dbartsControl, c(controlArgs, extraControl))

  treePrior <- quote(cgm(power, base))
  treePrior[[2L]] <- power
  treePrior[[3L]] <- base

  samplerEnv <- new.env(parent = fitEnv)
  samplerEnv[["responseData"]] <- responseData
  samplerEnv[["control"]] <- control
  samplerEnv[["forests"]] <- list(
    dbarts::forest(vars = muVars, sd = sd.control, update.amplitude = update.a,
                   interactions = mu.interactions),
    dbarts::forest(vars = tauVars, n.trees = coerceOrError(n.trees.treatment, "integer")[1L],
                   base = treatment.base, power = treatment.power,
                   sd = sd.moderate, amplitude.prior.variance = b.prior.variance,
                   update.amplitude = update.b, interactions = tau.interactions,
                   blocks = tau.blocks))

  samplerCall <- quote(dbarts::dbarts(responseData, control = control, forests = forests))
  samplerCall$tree.prior <- treePrior
  samplerFormals <- names(formals(dbarts::dbarts))
  for (argName in names(extraArgs)[names(extraArgs) %in% samplerFormals &
                                   names(extraArgs) %not_in% names(samplerCall)]) {
    samplerEnv[[paste0("bcf.arg.", argName)]] <- extraArgs[[argName]]
    samplerCall[[argName]] <- as.symbol(paste0("bcf.arg.", argName))
  }

  sampler <- eval(samplerCall, envir = samplerEnv)

  ## the response transform, read before the sampler can go out of scope; the
  ## rows are per chain and the linear map is shared, which is asserted here
  ## rather than assumed
  calibration <- sampler$getCalibration(1L)
  response.scale <- calibration[1L, "response.scale"]
  response.shift <- calibration[1L, "response.shift"]
  if (any(calibration[, "response.scale"] != response.scale) ||
      any(calibration[, "response.shift"] != response.shift))
    stop("chains disagree on the response transform; the counterfactual identity needs one linear map")

  family <- sampler$model@family
  if (family == "auto") family <- if (sampler$control@binary) "probit" else "gaussian"
  responseIsBinary <- sampler$control@binary

  sampler$sampleTreesFromPrior(updateState = FALSE)
  ## run(0L, 0L) returns NULL rather than an empty set of channels, so a
  ## zero-burn fit skips the burn run outright
  burn <- if (n.burn > 0L) sampler$run(0L, n.burn, updateState = FALSE) else NULL
  samples <- sampler$run(0L, n.samples, updateState = FALSE)

  trt <- as.vector(responseData@x[, treatmentName])

  ## the counterfactual is a swap of the amplitude the treatment forest is
  ## multiplied by; both the offset and the response shift cancel out of the
  ## difference, so only the scale survives:
  ##   mu.hat.cf = train + response.scale * (b_{1-z} - b_z) * tau
  ## and b_{1-z} - b_z is (1 - 2z) * (b.1 - b.0)
  tau.internal <- sliceBCFForest(samples$forestFits, 2L)
  b.0 <- sliceBCFGlue(samples$glue, 2L)
  b.1 <- sliceBCFGlue(samples$glue, 3L)
  mu.hat.obs <- samples$train
  mu.hat.cf  <- mu.hat.obs + response.scale * outer(1 - 2 * trt, b.1 - b.0) * tau.internal

  ## for a binary response the combination is on the latent scale, so the link
  ## is applied LAST and to both surfaces
  if (responseIsBinary) {
    link <- if (family == "logistic") plogis else pnorm
    mu.hat.obs <- link(mu.hat.obs)
    mu.hat.cf  <- link(mu.hat.cf)
  }

  forestFits <- lapply(seq_along(bcfForestNames),
                       function(index) convertBCFSamples(sliceBCFForest(samples$forestFits, index)))
  names(forestFits) <- bcfForestNames

  varcount <- lapply(seq_along(bcfForestNames), function(index) {
    counts <- convertBCFSamples(sliceBCFForest(samples$varcount, index))
    dimnames(counts) <- c(rep(list(NULL), length(dim(counts)) - 1L), list(designNames))
    counts
  })
  names(varcount) <- bcfForestNames

  glue <- convertBCFSamples(samples$glue)
  dimnames(glue) <- c(rep(list(NULL), length(dim(glue)) - 1L), list(c("a", "b.0", "b.1")))

  result <- namedList(
    forests = forestFits,
    glue,
    mu.hat.obs = convertBCFSamples(mu.hat.obs),
    mu.hat.cf  = convertBCFSamples(mu.hat.cf),
    varcount,
    y = responseData@y,
    trt,
    name.trt = treatmentName,
    name.p.score = p.scoreName,
    family,
    response.scale,
    response.shift,
    n.trees = c(mu = control@n.trees, tau = coerceOrError(n.trees.treatment, "integer")[1L]),
    n.chains = n.chains,
    n.samples = n.samples,
    n.burn = n.burn,
    n.obs = n.obs,
    combineChains = as.logical(combineChains)[1L],
    data = responseData,
    call = call)

  ## a binary fit has no residual standard deviation: run() reports 1s, and
  ## bartCause keys "is this model binary?" on the element's absence
  if (!responseIsBinary) {
    result$sigma <- asBCFSigmaMatrix(samples$sigma, n.chains)
    result$first.sigma <-
      if (is.null(burn)) matrix(numeric(0L), n.chains, 0L)
      else asBCFSigmaMatrix(burn$sigma, n.chains)
  }

  if (keepSampler) result$fit <- sampler

  class(result) <- "bartBCF"
  result
}

bcf <- function(formula, data, subset, weights, offset,
                treatment,
                moderators = NULL,
                p.score = NULL,
                n.trees = 200L, base = 0.95, power = 2.0,
                n.trees.treatment = 50L,
                treatment.base = 0.25, treatment.power = 3,
                sd.control = NULL, sd.moderate = 1,
                b.prior.variance = 0.5,
                update.a = TRUE, update.b = TRUE,
                mu.interactions = NULL, tau.interactions = NULL,
                tau.blocks = NULL,
                n.samples = 500L, n.burn = 500L, n.chains = 4L,
                n.threads = dbarts::guessNumCores(), combineChains = TRUE,
                keepSampler = FALSE, verbose = TRUE, seed = NA_integer_, ...)
{
  matchedCall <- match.call()
  callingEnv  <- parent.frame(1L)

  dataAreMissing    <- missing(data)
  subsetIsMissing   <- missing(subset)
  weightsAreMissing <- missing(weights)
  offsetIsMissing   <- missing(offset)

  if (missing(formula))   stop("'formula' must be specified")
  if (missing(treatment)) stop("'treatment' variable must be specified")

  if (!is.null(matchedCall[["mu.blocks"]]) || !is.null(matchedCall[["blocks"]]))
    stop("a prognostic-forest block partition is not supported: 'blocks' is resolved against the ",
         "full design and must cover every column, including the treatment column the prognostic ",
         "forest is masked out of. Use 'tau.blocks' to block the treatment forest; blocking the ",
         "prognostic forest against its own column mask is not yet available")
  for (argName in c("forests", "bases", "test", "offset.test", "x.test"))
    if (!is.null(matchedCall[[argName]]))
      stop("'", argName, "' is set by bcf() itself and cannot be supplied")

  dataValue <- if (dataAreMissing) NULL else data
  usesFormula <- is.formula(formula)

  if (usesFormula && !dataAreMissing && !is.list(dataValue))
    stop("for the formula interface 'data' must be a data frame or list: bcf() adds the ",
         "treatment (and propensity score) to it as columns")
  if (!usesFormula && dataAreMissing)
    stop("'data' must be the response vector when 'formula' is a matrix of predictors")

  ## as.vector: a one-column matrix would enter the model frame as a matrix
  ## column and be expanded to a differently named design column, which the
  ## masks below key on by name
  z <- as.vector(resolveBCFColumn(matchedCall$treatment, dataValue, callingEnv))
  pihat <- if (is.null(matchedCall$p.score)) NULL else
    as.vector(resolveBCFColumn(matchedCall$p.score, dataValue, callingEnv))

  ## the treatment and the propensity score have to reach the model matrix as
  ## columns, so each needs a name that does not collide with the design's own
  reservedNames <- if (usesFormula) c(all.vars(formula), if (is.list(dataValue)) names(dataValue))
                   else colnames(as.matrix(formula))
  name.trt <- resolveBCFName(matchedCall$treatment, "z", reservedNames)
  name.p.score <- if (is.null(pihat)) NULL else
    resolveBCFName(matchedCall$p.score, "ps", c(reservedNames, name.trt))

  fitEnv <- new.env(parent = callingEnv)
  dbartsDataCall <- quote(dbarts::dbartsData(formula = bcf.formula, data = bcf.data))

  if (usesFormula) {
    modelFormula <- formula
    for (columnName in c(name.trt, name.p.score)) {
      newRHS <- quote(a + b)
      newRHS[[2L]] <- modelFormula[[3L]]
      newRHS[[3L]] <- str2lang(columnName)
      modelFormula[[3L]] <- newRHS
    }
    if (is.list(dataValue)) {
      dataValue[[name.trt]] <- z
      if (!is.null(name.p.score)) dataValue[[name.p.score]] <- pihat
      fitEnv[["bcf.data"]] <- dataValue
    } else {
      ## no data argument: the model frame falls back to the formula's own
      ## environment, so the columns are bound there
      fitEnv[[name.trt]] <- z
      if (!is.null(name.p.score)) fitEnv[[name.p.score]] <- pihat
      dbartsDataCall$data <- NULL
    }
    environment(modelFormula) <- fitEnv
    dbartsDataCall$formula <- modelFormula
  } else {
    x <- as.matrix(formula)
    if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
    x <- cbind(x, z)
    colnames(x)[ncol(x)] <- name.trt
    if (!is.null(name.p.score)) {
      x <- cbind(x, pihat)
      colnames(x)[ncol(x)] <- name.p.score
    }
    fitEnv[["bcf.x"]] <- x
    fitEnv[["bcf.data"]] <- data
    dbartsDataCall$formula <- quote(bcf.x)
  }

  if (!subsetIsMissing)   dbartsDataCall$subset  <- matchedCall$subset
  if (!weightsAreMissing) dbartsDataCall$weights <- matchedCall$weights
  if (!offsetIsMissing)   dbartsDataCall$offset  <- matchedCall$offset

  fitArgs <- list(dbartsDataCall = dbartsDataCall, evalEnv = fitEnv, z = z,
                  treatmentName = name.trt, p.scoreName = name.p.score,
                  moderators = moderators,
                  n.trees = n.trees, base = base, power = power,
                  n.trees.treatment = n.trees.treatment,
                  treatment.base = treatment.base, treatment.power = treatment.power,
                  sd.control = sd.control, sd.moderate = sd.moderate,
                  b.prior.variance = b.prior.variance,
                  update.a = update.a, update.b = update.b,
                  mu.interactions = mu.interactions, tau.interactions = tau.interactions,
                  tau.blocks = tau.blocks,
                  n.samples = n.samples, n.burn = n.burn, n.chains = n.chains,
                  n.threads = n.threads, combineChains = combineChains,
                  keepSampler = keepSampler, verbose = verbose, seed = seed,
                  call = matchedCall)

  ## quote = TRUE: the data call and the matched call are language objects, and
  ## do.call would otherwise splice them into the constructed call as code
  do.call(fitBCF, c(fitArgs, list(...)), quote = TRUE)
}

## Evaluate a column expression against 'data' first and the caller's frame
## second, the resolution order the rest of the package uses.
resolveBCFColumn <- function(expr, data, callingEnv)
{
  if (is.list(data)) {
    value <- tryCatch(eval(expr, data, callingEnv), error = function(e) e)
    if (!inherits(value, "error")) return(value)
  }
  eval(expr, callingEnv)
}

## A column supplied as a bare symbol keeps its own name -- it may already be in
## the design, in which case adding it again is a no-op. Anything else gets the
## default name, uniquified away from the design.
resolveBCFName <- function(expr, defaultName, reservedNames)
{
  if (is.symbol(expr)) return(as.character(expr))
  name <- defaultName
  while (name %in% reservedNames) name <- paste0(name, defaultName)
  name
}

print.bartBCF <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if (!is.null(x$call))
    cat("Call:\n", paste0(deparse(x$call), collapse = "\n  "), "\n\n", sep = "")

  cat("Bayesian causal forest, family '", x$family, "'\n", sep = "")
  cat("  forests: ", length(x$forests), " (",
      paste0(names(x$n.trees), ": ", x$n.trees, " trees", collapse = ", "), ")\n", sep = "")
  cat("  treatment: ", x$name.trt,
      if (!is.null(x$name.p.score)) paste0(", propensity score: ", x$name.p.score) else "",
      "\n", sep = "")
  cat("  ", x$n.chains, " chain", if (x$n.chains > 1L) "s" else "",
      " x ", x$n.samples, " samples, ", x$n.obs, " observations\n", sep = "")
  cat("  sampler ", if (!is.null(x$fit)) "kept" else "not kept", "\n", sep = "")

  invisible(x)
}

extract.bartBCF <-
  function(object,
           type = c("mu.obs", "mu.cf", "mu.1", "mu.0", "icate", "mu", "tau",
                    "glue", "sigma", "varcount"),
           forest = NULL,
           combineChains = object[["combineChains"]],
           ...)
{
  issueWarningForUnknownArguments()

  if (!is.character(type) || type[1L] %not_in% eval(formals(extract.bartBCF)$type))
    stop("type must be in '", paste0(eval(formals(extract.bartBCF)$type), collapse = "', '"), "'")
  type <- type[1L]

  if (!is.null(forest) && type != "varcount")
    stop("'forest' selects an element of a forest-indexed container and applies only to ",
         "type = 'varcount'; the per-forest fits are addressed as type = 'mu' and type = 'tau'")

  if (type == "sigma" && is.null(object[["sigma"]]))
    stop("binary response model does not have a residual standard deviation parameter (sigma)")

  n.chains <- object$n.chains

  if (type == "varcount") {
    result <- object$varcount
    if (!is.null(forest)) {
      if (length(forest) != 1L || (is.character(forest) && forest %not_in% names(result)) ||
          (!is.character(forest) && (forest < 1L || forest > length(result))))
        stop("'forest' must be one of '", paste0(names(result), collapse = "', '"),
             "' or an index into them")
      return(if (combineChains) combineChains(result[[forest]], n.chains) else result[[forest]])
    }
    return(if (combineChains) lapply(result, function(x) combineChains(x, n.chains)) else result)
  }

  ## icate is the treatment forest's contribution on the response scale,
  ## response.scale * (b.1 - b.0) * tau; on a gaussian fit it equals mu.1 - mu.0
  ## exactly, and on a binary one it is the LATENT-scale contrast, while mu.1
  ## and mu.0 have had the link applied
  result <- with(object, switch(type,
    mu.obs   = mu.hat.obs,
    mu.cf    = mu.hat.cf,
    mu.1     = bcfObsCfToTrtCtl(mu.hat.obs, mu.hat.cf, trt),
    mu.0     = bcfObsCfToTrtCtl(mu.hat.obs, mu.hat.cf, 1 - trt),
    icate    = response.scale * forests$tau * as.vector(bcfGlueDifference(glue)),
    mu       = forests$mu,
    tau      = forests$tau,
    glue     = glue,
    sigma    = sigma))

  if (combineChains) combineChains(result, n.chains) else result
}

fitted.bartBCF <-
  function(object,
           type = c("mu.obs", "mu.cf", "mu.1", "mu.0", "icate", "mu", "tau",
                    "glue", "sigma", "varcount"),
           ...)
{
  if (!is.character(type) || type[1L] %not_in% eval(formals(fitted.bartBCF)$type))
    stop("type must be in '", paste0(eval(formals(fitted.bartBCF)$type), collapse = "', '"), "'")
  type <- type[1L]

  result <- extract(object, type = type, combineChains = FALSE, ...)

  if (type == "sigma") return(mean(result))
  if (is.list(result)) return(lapply(result, function(x) apply(x, length(dim(x)), mean)))
  apply(result, length(dim(result)), mean)
}

residuals.bartBCF <- function(object, ...)
{
  issueWarningForUnknownArguments()

  as.vector(object$y) - fitted(object, type = "mu.obs")
}

predict.bartBCF <- function(object, newdata, type = c("mu", "tau"), ...)
{
  stop("predict is not available for a 'bartBCF' fit: out-of-sample mu(x) and tau(x) need ",
       "per-forest saved-tree replay, which dbarts does not expose, and the blended test ",
       "surface a multi-forest sampler would have to combine them through is refused at ",
       "creation. Use extract() or fitted() for the training-sample surfaces")
}
