
# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(DESeq2)
library(survival)

# Load residual -----------------------------------------------------------

residual <- readxl::read_excel(path = "data/metadata/Residual.xlsx") %>% 
  dplyr::filter(Residual != "NA") %>% 
  dplyr::mutate(residual = ifelse(Residual == "R0", "R0", "non-R0")) %>% 
  dplyr::select(barcode, residual)

platinum.se <- readr::read_rds(file = "data/rda/total351.platinum.se.rds.gz")
platinum.se@colData %>% 
  as.data.frame() %>% 
  dplyr::select(barcode, platinum) ->
  platinum

os_risk_group <- readxl::read_excel(path = "data/newoutput/Riskgroup-os.xlsx") %>% 
  dplyr::left_join(residual, by = "barcode") %>% 
  dplyr::left_join(platinum, by = "barcode")
pfs_risk_group <- readxl::read_excel(path = "data/newoutput/Riskgroup-pfs.xlsx") %>%
  dplyr::left_join(residual, by = "barcode") %>% 
  dplyr::left_join(platinum, by = "barcode")



# Function ----------------------------------------------------------------

fn_plot_risk_distribution <- function(.d) {
  .d %>% dplyr::filter(oc %in% c("OC521", "OC44")) %>% dplyr::pull(riskscore) %>% median ->
    .riskscore_cutoff
  .d %>% 
    dplyr::mutate(rank = rank(riskscore)) %>% 
    dplyr::mutate(event = as.factor(event)) %>% 
    ggplot(aes(x = rank, y = riskscore, fill = event, color = event)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = .riskscore_cutoff, linetype = 2) +
    scale_fill_viridis_d(name = "Status", direction = -1, label = c("Alive", "Dead")) +
    scale_color_viridis_d(name = "Status", direction = -1, label = c("Alive", "Dead")) +
    scale_x_continuous(expand = c(0.01,0)) +
    scale_y_continuous(expand = c(0.02,0)) +
    labs(
      x = "Patients",
      y = "Risk score"
    ) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = c(0.7, 0.1),
      legend.direction = "horizontal",
      plot.title = element_text(hjust = 0.5)
    )
}

human_read <- function(.x){
  .sign = ifelse(.x < 0 , TRUE, FALSE)
  .x <- abs(.x)
  
  if (.x >= 0.1) {
    .x %>% signif(digits = 2) %>% toString() -> .xx
  } else if (.x < 0.1 && .x >= 0.001 ) {
    .x %>% signif(digits = 2) %>% toString() -> .xx
  } else if (.x < 0.001 && .x > 0) {
    .x %>% format(digits = 3, scientific = TRUE) -> .xx
  } else {
    .xx <- '0'
  }
  
  ifelse(.sign, paste0('-',.xx), .xx)
}

human_read_latex_pval <- function(.x, .s = NA) {
  
  if (is.na(.s)) {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      latex2exp::TeX(glue::glue("<<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>"))
    } else {
      latex2exp::TeX(glue::glue("<<.x>>$", .open = "<<", .close = ">>"))
    }
  } else {
    if (grepl(pattern = "e", x = .x)) {
      sub("-0", "-", strsplit(split = "e", x = .x, fixed = TRUE)[[1]]) -> .xx
      latex2exp::TeX(glue::glue("<<.s>>, <<.xx[1]>> \\times 10^{<<.xx[2]>>}$", .open = "<<", .close = ">>"))
    } else {
      latex2exp::TeX(glue::glue("<<.s>>, <<.x>>$", .open = "<<", .close = ">>"))
    }
  }
  
}

fn_plot_hr <- function(.d) {
  .d %>% 
  ggplot(aes(x = hazard_ratio, y = formalname)) +
    geom_point(size = 3, color = "red", fill = "red", shape = 23) +
    geom_vline(xintercept = 1, linetype = 5, color = "black", size = 0.5) +
    geom_errorbarh(aes(xmax = ci.high, xmin = ci.low, height = 0.2), size = 1) +
    geom_text(aes(x = -8, y = formalname, label = hr_label), size = 6, hjust = 0.5) +
    geom_text(aes(x = 35, y = formalname, label = pval_label), size = 6, hjust = 1) +
    scale_y_discrete(expand = c(0.2, 0)) +
    scale_x_continuous(expand = c(0.18, 0, 0.05, 0)) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 14, colour = "black"),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 20),
      plot.title = element_text(size = 20, hjust = 0.3)
    ) +
    labs(x = "Hazard ratio") +
    annotate(geom = "text", x = -8, y = 8, label = "HR (95% CI)", size = 6, vjust = 1) +
    annotate(geom = "text", x = 35, y = 8, label = "P value", size = 6, vjust = 1, hjust = 1)
}

fn_unicox <- function(.group, .data) {
  
  coxph(formula = Surv(duration, event) ~ rlang::eval_bare(rlang::sym(.group)), data = .data) %>% 
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
    dplyr::select(
      term,
      hazard_ratio = estimate,
      ci.low = conf.low,
      ci.high = conf.high,
      pval = p.value
    )
}

# Multi-Cox ---------------------------------------------------------------

# OS ----------------------------------------------------------------------

os_risk_group %>% 
  dplyr::select(barcode, oc, stage, CA125, age, platelet_count, platinum, riskscore, riskscore_group = group, residual, event, duration) ->
  os_risk_group_s

os_risk_group_s %>% 
  fn_plot_risk_distribution() +
  labs(title = "OS risk score distribution") ->
  os_risk_group_s_riskscore_dis_plot

ggsave(
  filename = "OS-riskscore-distribution.pdf",
  plot = os_risk_group_s_riskscore_dis_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 5,
  height = 3.5
)

os_risk_group_s %>%
  tidyr::drop_na() %>% 
  dplyr::mutate(stage_group = factor(stage, levels = c("E", "L"))) %>% 
  dplyr::mutate(ca125_group = factor(ifelse(CA125 > 35, "CA125>35", "CA125<=35"), levels = c("CA125<=35", "CA125>35"))) %>% 
  dplyr::mutate(age_group = factor(ifelse(age > 50, "age>50", "age<=50"), levels = c("age<=50", "age>50"))) %>% 
  dplyr::mutate(plc_group = factor(ifelse(platelet_count > 350, "plc>350", "plc<=350"), levels = c("plc<=350", "plc>350"))) %>% 
  dplyr::mutate(platinum_group = factor(platinum, levels = c("sensitive", "resistant"), ordered = FALSE)) %>% 
  dplyr::mutate(residual_group = factor(residual, levels = c("R0", "non-R0"))) %>% 
  dplyr::mutate(riskscore_group = factor(riskscore_group, levels = c("Low", "High")))->
  os_risk_group_s_s

# Unicox ------------------------------------------------------------------
unicox_df <- 
  tibble::tibble(group = c(
    "riskscore_group",
    "ca125_group",
    "age_group",
    "plc_group",
    "stage_group",
    "residual_group",
    "platinum_group"
  ), formalname = c(
    "Risk score (high vs low)",
    "CA125 (>35 vs <=35)",
    "Age (>50 vs <=50)",
    "PLC (>350 vs <= 350)",
    "Stage (late vs early)",
    "Debulking (suboptimal vs optimal)",
    "Platinum (resitant vs senstive)"
  ))

unicox_df %>% 
  dplyr::mutate(unicox = purrr::map(.x = group, .f = fn_unicox, .data = os_risk_group_s_s)) %>% 
  tidyr::unnest(cols = c(unicox)) %>% 
  dplyr::select(term = group, hazard_ratio, ci.low, ci.high, pval, formalname) %>% 
  dplyr::arrange(-hazard_ratio) %>% 
  dplyr::mutate(formalname = factor(formalname, levels = formalname)) ->
  os_unicox_reg;os_unicox_reg

os_unicox_reg %>% 
  dplyr::mutate(hr_label = glue::glue("{round(hazard_ratio, 2)}({round(ci.low, 2)}-{round(ci.high, 2)})")) %>% 
  dplyr::mutate(pval_label = signif(pval, 2)) %>% 
  fn_plot_hr() +
  labs(title = "OS UniCox HR") ->
  os_hr_ucox_plot;os_hr_ucox_plot

ggsave(
  filename = "OS-HR-UCox.pdf",
  plot = os_hr_ucox_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 10,
  height = 5
)

# Multicox ----------------------------------------------------------------

coxph(formula = Surv(duration, event) ~ 
        riskscore_group + 
        ca125_group +
        age_group + 
        plc_group + 
        stage_group +
        residual_group +
        platinum_group, 
      data = os_risk_group_s_s) ->
  os.coxph

os.coxph %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    ci.low = conf.low,
    ci.high = conf.high,
    pval = p.value
  ) %>% 
  dplyr::mutate(formalname = c(
    "Risk score (high vs low)",
    "CA125 (>35 vs <=35)",
    "Age (>50 vs <=50)",
    "PLC (>350 vs <= 350)",
    "Stage (late vs early)",
    "Debulking (suboptimal vs optimal)",
    "Platinum (resitant vs senstive)"
    )) %>% 
  dplyr::arrange(-hazard_ratio) %>% 
  dplyr::mutate(formalname = factor(formalname, levels = formalname)) ->
  os_multicox_reg;os_multicox_reg

writexl::write_xlsx(os_multicox_reg, path = "data/newoutput/OS-multicox-reg.xlsx")

os_multicox_reg %>% 
  dplyr::mutate(hr_label = glue::glue("{round(hazard_ratio, 2)}({round(ci.low, 2)}-{round(ci.high, 2)})")) %>% 
  dplyr::mutate(pval_label = signif(pval, 2)) %>% 
  fn_plot_hr() +
  labs(title = "OS MultiCox HR") ->
  os_hr_mcox_plot;os_hr_mcox_plot

ggsave(
  filename = "OS-HR-MCox.pdf",
  plot = os_hr_mcox_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 10,
  height = 5
)


# PFS ---------------------------------------------------------------------

pfs_risk_group %>% 
  dplyr::select(barcode, oc, stage, CA125, age, platelet_count, platinum, riskscore, riskscore_group = group, residual, event, duration) ->
  pfs_risk_group_s


pfs_risk_group_s %>% 
  fn_plot_risk_distribution() +
  labs(title = "PFS risk score distribution") ->
  pfs_risk_group_s_riskscore_dis_plot

ggsave(
  filename = "PFS-riskscore-distribution.pdf",
  plot = pfs_risk_group_s_riskscore_dis_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 5,
  height = 3.5
)


pfs_risk_group_s %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(stage_group = factor(stage, levels = c("E", "L"))) %>% 
  dplyr::mutate(ca125_group = factor(ifelse(CA125 > 35, "CA125>35", "CA125<=35"), levels = c("CA125<=35", "CA125>35"))) %>% 
  dplyr::mutate(age_group = factor(ifelse(age > 50, "age>50", "age<=50"), levels = c("age<=50", "age>50"))) %>% 
  dplyr::mutate(plc_group = factor(ifelse(platelet_count > 300, "plc>350", "plc<=350"), levels = c("plc<=350", "plc>350"))) %>% 
  dplyr::mutate(platinum_group = factor(platinum, levels = c("sensitive", "resistant"), ordered = FALSE)) %>% 
  dplyr::mutate(residual_group = factor(residual, levels = c("R0", "non-R0"))) %>% 
  dplyr::mutate(riskscore_group = factor(riskscore_group, levels = c("Low", "High")))->
  pfs_risk_group_s_s



# Unicox ------------------------------------------------------------------


unicox_df %>% 
  dplyr::mutate(unicox = purrr::map(.x = group, .f = fn_unicox, .data = pfs_risk_group_s_s)) %>% 
  tidyr::unnest(cols = c(unicox)) %>% 
  dplyr::select(term = group, hazard_ratio, ci.low, ci.high, pval, formalname) %>% 
  dplyr::arrange(-hazard_ratio) %>% 
  dplyr::mutate(formalname = factor(formalname, levels = formalname)) ->
  pfs_unicox_reg;pfs_unicox_reg

pfs_unicox_reg %>% 
  dplyr::mutate(hr_label = glue::glue("{round(hazard_ratio, 2)}({round(ci.low, 2)}-{round(ci.high, 2)})")) %>% 
  dplyr::mutate(pval_label = signif(pval, 2)) %>% 
  fn_plot_hr() +
  labs(title = "OS UniCox HR") ->
  pfs_hr_ucox_plot;pfs_hr_ucox_plot

ggsave(
  filename = "PFS-HR-UCox.pdf",
  plot = pfs_hr_ucox_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 10,
  height = 5
)

# Multicox ----------------------------------------------------------------



coxph(formula = Surv(duration, event) ~ 
        riskscore_group + 
        ca125_group +
        age_group + 
        plc_group + 
        stage_group +
        residual_group +
        platinum_group, 
      data = pfs_risk_group_s_s) ->
  pfs.coxph


pfs.coxph %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  dplyr::select(
    term,
    hazard_ratio = estimate,
    ci.low = conf.low,
    ci.high = conf.high,
    pval = p.value
  ) %>% 
  dplyr::mutate(formalname = c(
    "Risk score (high vs low)",
    "CA125 (>35 vs <=35)",
    "Age (>50 vs <=50)",
    "PLC (>350 vs <= 350)",
    "Stage (late vs early)",
    "Debulking (suboptimal vs optimal)",
    "Platinum (resitant vs senstive)"
  )) %>% 
  dplyr::arrange(-hazard_ratio) %>% 
  dplyr::mutate(formalname = factor(formalname, levels = formalname))  ->
  pfs_multicox_reg

writexl::write_xlsx(pfs_multicox_reg, path = "data/newoutput/PFS-multicox-reg.xlsx")

pfs_multicox_reg %>% 
  dplyr::mutate(hr_label = glue::glue("{round(hazard_ratio, 2)}({round(ci.low, 2)}-{round(ci.high, 2)})")) %>% 
  dplyr::mutate(pval_label = signif(pval, 2)) %>% 
  fn_plot_hr() +
  labs(title = "PFS MultiCox HR") ->
  pfs_hr_mcox_plot;pfs_hr_mcox_plot

ggsave(
  filename = "PFS-HR-MCox.pdf",
  plot = pfs_hr_mcox_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 10,
  height = 5
)



# Save image --------------------------------------------------------------

save.image('data/rda/16-multivariate.rda')