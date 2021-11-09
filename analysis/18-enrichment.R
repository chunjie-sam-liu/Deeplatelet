# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Nov  9 19:41:35 2021
# @DESCRIPTION: 18-enrichment.R

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)

# src ---------------------------------------------------------------------

source(file = "src/convertID.R")

# Load panel --------------------------------------------------------------

os.panel <- readr::read_rds(file='data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)
pfs.panel <- readr::read_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)


# Convert ID --------------------------------------------------------------

os.panel.df <- fn_convertId(ids = os.panel)
pfs.panel.df <- fn_convertId(ids = pfs.panel)

os.pfs.panel.df <- list(
  "OS gene panel" = os.panel.df,
  "PFS gene panel" = pfs.panel.df
)

writexl::write_xlsx(
  x = os.pfs.panel.df, 
  path = "data/newoutput/OS-PFS-Panel-Description.xlsx"
  )


# MSigDB ------------------------------------------------------------------



# OS ----------------------------------------------------------------------



msig_df <- msigdbr(species = "Homo sapiens")

os.panel.msig <- enricher(
  os.panel,
  TERM2GENE = msig_df %>%
    dplyr::mutate(gs_name = glue::glue("{gs_cat}#{gs_name}")) %>%
    dplyr::select(gs_name, ensembl_gene)
)
os.panel.msig %>%
  tibble::as_tibble() %>%
  tidyr::separate(col = Description, into = c("gs_cat", "Description"), sep="#") ->
  os.panel.msig.result

os.panel.msig.result %>%
  dplyr::filter(gs_cat %in% c("C6", "H"))->
  os.panel.msig.result.c6.h

os.panel.msig.result.c6.h %>%
  dplyr::select(Description, p.adjust, Count, gs_cat) %>%
  dplyr::mutate(qvalue = -log10(p.adjust)) %>%
  dplyr::arrange(dplyr::desc(gs_cat), -qvalue) ->
  os.panel.msig.result.c6.h_ready

os.panel.msig.result.c6.h_ready %>%
  ggplot(aes(x = Description)) +
  geom_col(aes( y = qvalue, fill = gs_cat)) +
  geom_line(aes(y = Count, group = 1), color = "#B22222") +
  scale_y_continuous(
    name = "-log10(p.adjust)",
    sec.axis = sec_axis(~., name="Count"),
    expand = c(0.02, 0)
  ) +
  scale_x_discrete(
    limit = os.panel.msig.result.c6.h_ready$Description,
    labels =stringr::str_wrap(
      stringr::str_replace_all(
        string = os.panel.msig.result.c6.h_ready$Description,
        pattern = "_",
        replacement = " "
      ), 
      width = 10)
  ) +
  scale_fill_viridis_d(
    label = c("H6: Oncogenic signature", "Hallmark"),
    name = "MSigDB"
  ) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_line(color = "#B22222"),
    axis.text.y.right = element_text(color = "#B22222"),
    axis.title.y.right = element_text(color = "#B22222"),
    axis.ticks.y.right = element_line(color = "#B22222"),
    #legend.position = c(0.7, 0.7)
    legend.position = "top"
  ) +
  labs(x = "", y = "-log10(qvalue)") ->
  os.panel.msig.result.c6.h_ready_plot;os.panel.msig.result.c6.h_ready_plot

ggsave(
  filename = "OS-Panel-MSigDB-enrichment.pdf",
  path = "data/newoutput",
  plot = os.panel.msig.result.c6.h_ready_plot,
  device = "pdf",
  width =6,
  height = 4
)


# PFS ---------------------------------------------------------------------


pfs.panel.msig <- enricher(
  pfs.panel,
  TERM2GENE = msig_df %>%
    dplyr::mutate(gs_name = glue::glue("{gs_cat}#{gs_name}")) %>%
    dplyr::select(gs_name, ensembl_gene)
)
pfs.panel.msig %>%
  tibble::as_tibble() %>%
  tidyr::separate(col = Description, into = c("gs_cat", "Description"), sep="#") ->
  pfs.panel.msig.result

pfs.panel.msig.result %>%
  dplyr::filter(gs_cat %in% c("C6", "H"))->
  pfs.panel.msig.result.c6.h

pfs.panel.msig.result.c6.h %>%
  dplyr::select(Description, p.adjust, Count, gs_cat) %>%
  dplyr::mutate(qvalue = -log10(p.adjust)) %>%
  dplyr::arrange(dplyr::desc(gs_cat), -qvalue) ->
  pfs.panel.msig.result.c6.h_ready

pfs.panel.msig.result.c6.h_ready %>%
  ggplot(aes(x = Description)) +
  geom_col(aes( y = qvalue, fill = gs_cat)) +
  geom_line(aes(y = Count, group = 1), color = "#B22222") +
  scale_y_continuous(
    name = "-log10(p.adjust)",
    sec.axis = sec_axis(~., name="Count"),
    expand = c(0.02, 0)
  ) +
  scale_x_discrete(
    limit = pfs.panel.msig.result.c6.h_ready$Description,
    labels =stringr::str_wrap(
      stringr::str_replace_all(
        string = pfs.panel.msig.result.c6.h_ready$Description,
        pattern = "_",
        replacement = " "
      ), 
      width = 10)
  ) +
  scale_fill_viridis_d(
    label = c("H6: Oncogenic signature", "Hallmark"),
    name = "MSigDB"
  ) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_line(color = "#B22222"),
    axis.text.y.right = element_text(color = "#B22222"),
    axis.title.y.right = element_text(color = "#B22222"),
    axis.ticks.y.right = element_line(color = "#B22222"),
    #legend.position = c(0.7, 0.7)
    legend.position = "top"
  ) +
  labs(x = "", y = "-log10(qvalue)") ->
  pfs.panel.msig.result.c6.h_ready_plot;pfs.panel.msig.result.c6.h_ready_plot

ggsave(
  filename = "PFS-Panel-MSigDB-enrichment.pdf",
  path = "data/newoutput",
  plot = pfs.panel.msig.result.c6.h_ready_plot,
  device = "pdf",
  width =6,
  height = 4
)


# save image --------------------------------------------------------------

save.image(file = "data/rda/18-enrichment.rda")
