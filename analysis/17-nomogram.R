# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Aug 13 03:59:03 2021
# @DESCRIPTION: nomogram.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# Load data ---------------------------------------------------------------

os_risk_group_s_s <- readr::read_rds(file = "data/rda/os_risk_group_s_s.rds.gz")
pfs_risk_group_s_s <- readr::read_rds(file = "data/rda/pfs_risk_group_s_s.rds.gz")

# Function ----------------------------------------------------------------



# Nomogram ----------------------------------------------------------------


# OS ----------------------------------------------------------------------
os_risk_group_s_s_s <- os_risk_group_s_s %>% 
  dplyr::select(duration, event, CA125, platelet_count, stage_group, residual_group, riskscore_group, age) %>% 
  dplyr::mutate(stage_group = plyr::revalue(x = stage_group, replace = c("E" = "Early", "L" = "Late"))) %>% 
  dplyr::mutate(residual_group = plyr::revalue(x = residual_group, replace = c("R0" = "Optimal", "non-R0" = "Suboptimal")))

ddist <- rms::datadist(os_risk_group_s_s_s)
options(datadist = "ddist")

var.labels <- c(
  duration = "duration", 
  event = "event",
  CA125 = "CA125",
  platelet_count = "Platelet count",
  stage_group = "FIGO stage",
  residual_group = "Debulking",
  riskscore_group = "Risk score",
  age = "Age"
)

Hmisc::label(os_risk_group_s_s_s) <- lapply(names(var.labels), function(x) {
  Hmisc::label(os_risk_group_s_s_s[, x]) <- var.labels[x]
})

mod.cox <- rms::cph(
  formula = survival::Surv(duration, event) ~  CA125 + platelet_count + stage_group + residual_group + riskscore_group, 
  data = os_risk_group_s_s_s,
  surv = TRUE,
  x = TRUE,
  y = TRUE
)

surv.cox <- rms::Survival(mod.cox)
nom.cox <- rms::nomogram(
  fit = mod.cox,
  fun = list(
    function(x) surv.cox(12 * 3, x),
    function(x) surv.cox(12 * 5, x)
  ),
  funlabel = c(
    "3-Year OS probability",
    "5-Year OS probability"
  ),
  lp = FALSE
)

pdf(file = "data/newoutput/OS-nomogram.pdf", width = 12, height = 6)
plot(nom.cox, xfrac = 0.3, total.points.label = "Sum of all points", cex.axis = 1.1, force.label = FALSE, tcl = 0.3, lmgp = 0.1, vnames = "labels")
dev.off()

# PFS ---------------------------------------------------------------------

pfs_risk_group_s_s_s <- pfs_risk_group_s_s %>% 
  dplyr::select(duration, event, CA125, platelet_count, stage_group, residual_group, riskscore_group, age) %>% 
  dplyr::mutate(stage_group = plyr::revalue(x = stage_group, replace = c("E" = "Early", "L" = "Late"))) %>% 
  dplyr::mutate(residual_group = plyr::revalue(x = residual_group, replace = c("R0" = "Optimal", "non-R0" = "Suboptimal")))

ddist <- rms::datadist(pfs_risk_group_s_s_s)
options(datadist = "ddist")

var.labels <- c(
  duration = "duration", 
  event = "event",
  CA125 = "CA125",
  platelet_count = "Platelet count",
  stage_group = "FIGO stage",
  residual_group = "Debulking",
  riskscore_group = "Risk score",
  age = "Age"
)

Hmisc::label(pfs_risk_group_s_s_s) <- lapply(names(var.labels), function(x) {
  Hmisc::label(pfs_risk_group_s_s_s[, x]) <- var.labels[x]
})

mod.cox <- rms::cph(
  formula = survival::Surv(duration, event) ~  CA125 + platelet_count + stage_group + residual_group + riskscore_group, 
  data = pfs_risk_group_s_s_s,
  surv = TRUE,
  x = TRUE,
  y = TRUE
)

surv.cox <- rms::Survival(mod.cox)
nom.cox <- rms::nomogram(
  fit = mod.cox,
  fun = list(
    function(x) surv.cox(12 * 1, x),
    function(x) surv.cox(12 * 2, x),
    function(x) surv.cox(12 * 3, x)
  ),
  funlabel = c(
    "1-Year PFS probability",
    "2-Year PFS probability",
    "3-Year PFS probability"
  ),
  lp = FALSE
)

pdf(file = "data/newoutput/PFS-nomogram.pdf", width = 12, height = 6)
plot(nom.cox, xfrac = 0.3, total.points.label = "Sum of all points", cex.axis = 1.1, force.label = FALSE, tcl = 0.3, lmgp = 0.1, vnames = "labels")
dev.off()


# Save image --------------------------------------------------------------

save.image(file = "data/rda/17-nomogram.rda")
