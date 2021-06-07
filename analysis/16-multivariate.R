
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
  dplyr::mutate(stage_group = factor(stage, levels = c("E", "L"))) %>% 
  dplyr::mutate(ca125_group = factor(ifelse(CA125 > 35, "CA125>35", "CA125<=35"), levels = c("CA125<=35", "CA125>35"))) %>% 
  dplyr::mutate(age_group = factor(ifelse(age > 50, "age>50", "age<=50"), levels = c("age<=50", "age>50"))) %>% 
  dplyr::mutate(plc_group = factor(ifelse(platelet_count > 350, "plc>350", "plc<=350"), levels = c("plc<=350", "plc>350"))) %>% 
  dplyr::mutate(platinum_group = factor(platinum, levels = c("sensitive", "resistant"), ordered = FALSE)) %>% 
  dplyr::mutate(residual_group = factor(residual, levels = c("non-R0", "R0"))) %>% 
  dplyr::mutate(riskscore_group = factor(riskscore_group, levels = c("Low", "High")))->
  os_risk_group_s_s

coxph(formula = Surv(duration, event) ~ 
        riskscore_group + 
        age_group + 
        plc_group + 
        stage_group +
        residual_group +
        platinum_group, 
      data = os_risk_group_s_s) %>% 
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
    "Age (>50 vs <=50)",
    "PLC (>350 vs <= 350)",
    "Stage (late vs early)",
    "Debulking (suboptimal vs optimal)",
    "Platinum (resitant vs senstive)"
    )) %>% 
  dplyr::arrange(-hazard_ratio) %>% 
  dplyr::mutate(formalname = factor(formalname, levels = formalname)) ->
  os_multicox_reg

writexl::write_xlsx(os_multicox_reg, path = "data/newoutput/OS-multicox-reg.xlsx")

os_multicox_reg %>% 
  ggplot(aes(x = hazard_ratio, y = formalname)) +
  geom_point()


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
  dplyr::mutate(stage_group = factor(stage, levels = c("E", "L"))) %>% 
  dplyr::mutate(ca125_group = factor(ifelse(CA125 > 35, "CA125>35", "CA125<=35"), levels = c("CA125<=35", "CA125>35"))) %>% 
  dplyr::mutate(age_group = factor(ifelse(age > 50, "age>50", "age<=50"), levels = c("age<=50", "age>50"))) %>% 
  dplyr::mutate(plc_group = factor(ifelse(platelet_count > 350, "plc>350", "plc<=350"), levels = c("plc<=350", "plc>350"))) %>% 
  dplyr::mutate(platinum_group = factor(platinum, levels = c("sensitive", "resistant"), ordered = FALSE)) %>% 
  dplyr::mutate(residual_group = factor(residual, levels = c("non-R0", "R0"))) %>% 
  dplyr::mutate(riskscore_group = factor(riskscore_group, levels = c("Low", "High")))->
  pfs_risk_group_s_s

coxph(formula = Surv(duration, event) ~ 
        riskscore_group + 
        age_group + 
        plc_group + 
        stage_group +
        residual_group +
        platinum_group, 
      data = pfs_risk_group_s_s) %>% 
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

