getCommonSupportSubset <- function(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows)
{
  if (!is.character(commonSup.rule) || commonSup.rule[1L] %not_in% c("none", "sd", "chisq"))
    stop("commonSup.rule must be one of 'none', 'sd', or 'chisq'")
  commonSup.rule <- commonSup.rule[1L]
  
  if (commonSup.rule == "none") return(rep_len(TRUE, length(sd.obs)))
  
  if (commonSup.rule == "sd") {
    if (!is.numeric(commonSup.cut) || !is.finite(commonSup.cut))
      stop("for common support rule 'sd', commonSup.cut must be a real number")
  } else {
    if (!is.numeric(commonSup.cut) || commonSup.cut[1L] <= 0 || commonSup.cut[1L] >= 1)
      stop("for common support rule 'chisq', commonSup.cut must be in (0, 1)")
  }
  commonSup.cut <- commonSup.cut[1L]
  
  stat <- getCommonSupportStatistic(sd.obs, sd.cf, commonSup.rule, commonSup.cut)
  cut  <- getCommonSupportCutoff(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows)
  
  if (commonSup.rule == "sd")
    stat <= ifelse(trt, cut["trt"], cut["ctl"])
  else
    stat <= cut
}

getCommonSupportStatistic <- function(sd.obs, sd.cf, commonSup.rule, commonSup.cut)
  switch(commonSup.rule,
         sd = sd.cf,
         chisq = (sd.cf / sd.obs)^2,
         none = rep_len(NA_real_, length(sd.obs)))

getCommonSupportCutoff <- function(sd.obs, sd.cf, commonSup.rule, commonSup.cut, trt, missingRows)
  switch(commonSup.rule,
         sd = c(trt = max(sd.obs[!missingRows & trt > 0]), ctl = max(sd.obs[!missingRows & trt <= 0])) + commonSup.cut * sd(sd.obs[!missingRows]),
         chisq = qchisq(1 - commonSup.cut, 1),
         none = rep_len(NA_real_, length(sd.obs)))
