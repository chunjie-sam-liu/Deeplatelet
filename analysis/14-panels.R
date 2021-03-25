
# Library -----------------------------------------------------------------

library(magrittr)
library(ggvenn)


# Load data ---------------------------------------------------------------


os.panel <- readr::read_rds(file='data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)
pfs.panel <- readr::read_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)
platinum.panel <- readr::read_rds(file = 'data/rda/panel.rds.gz')


ggvenn(
  list(
    "OS Panel" = os.panel,
    "PFS panel" = pfs.panel,
    "Platinum panel" = platinum.panel
  ),
  fill_color = RColorBrewer::brewer.pal(3, "Set3"),
  stroke_size = 0.5, set_name_size = 4
) ->
  panel_venn_plot

ggsave(
  filename = "Panel-venn.pdf",
  plot = panel_venn_plot,
  device = "pdf",
  path = "data/newoutput",
  width = 7,
  height = 6
)


# Save image --------------------------------------------------------------

save.image(file = "data/rda/14-panels.rda")
