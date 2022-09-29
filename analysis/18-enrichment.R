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
library(org.Hs.eg.db)
# src ---------------------------------------------------------------------

source(file = "src/convertID.R")

# Load panel --------------------------------------------------------------

os.panel <- readr::read_rds(file='data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)
pfs.panel <- readr::read_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::pull(ensid)
# platinum.panel <- readr::read_rds(file = 'data/rda/panel.rds.gz')

# Convert ID --------------------------------------------------------------

os.panel.df <- fn_convertId(ids = os.panel)
pfs.panel.df <- fn_convertId(ids = pfs.panel) %>% 
  dplyr::filter(entrezgene_id != 100874074 )
# platinum.panel.df <- fn_convertId(ids = platinum.panel)
# platinum.panel.df %>% 
#   dplyr::group_by(ensembl_gene_id) %>% 
#   dplyr::filter(dplyr::n() > 1)

total416.os.expr.coxph.hazard_ratio <- readr::read_rds(file='data/rda/total416.os.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::rename(ensembl_gene_id = ensid)
total434.pfs.expr.coxph.hazard_ratio <- readr::read_rds(file = 'data/rda/total434.pfs.expr.coxph.hazard_ratio.rds.gz') %>% 
  dplyr::rename(ensembl_gene_id = ensid)


os.pfs.platinum.panel.df <- list(
  "OS gene panel" = os.panel.df %>% 
    dplyr::inner_join(total416.os.expr.coxph.hazard_ratio, by = "ensembl_gene_id"),
  "PFS gene panel" = pfs.panel.df %>% 
    dplyr::inner_join(total434.pfs.expr.coxph.hazard_ratio, by = "ensembl_gene_id"),
  "Platinum gene panel" = platinum.panel.df
)

os.pfs.platinum.panel.df <- list(
  "OS gene panel" = total416.os.expr.coxph.hazard_ratio %>% 
    dplyr::select(-name) %>% 
    dplyr::left_join(os.panel.df, by = "ensembl_gene_id") %>% 
    dplyr::filter(ensembl_gene_id  != "ENSG00000259174") %>% 
    dplyr::filter(entrezgene_id != 107080638 |is.na(entrezgene_id) )
)

writexl::write_xlsx(
  x = os.pfs.platinum.panel.df, 
  path = "data/newoutput/Panel-Description.xlsx"
  )


# MSigDB ------------------------------------------------------------------



# OS ----------------------------------------------------------------------
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

# GO ----------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
# a <- groupGO(
#   gene = os.panel,
#   OrgDb = org.Hs.eg.db,
#   keyType = "ENSEMBL",
#   ont      = "BP",
#   readable = TRUE
# )
ego_all <- enrichGO(
  gene = os.panel.df$ensembl_gene_id, 
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont           = "ALL",
  # pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  # qvalueCutoff  = 0.05,
  minGSSize = 2,
  readable      = TRUE
)

p <- dotplot(ego_all, showCategory = 10)
ggsave(
  filename = "enrichgo.pdf",
  plot = p,
  device = "pdf",
  path = "data/review",
  width = 7,
  height = 5
)
as.data.frame(ego_all)

as.data.frame(ego_all) %>% 
  head()

ego_kegg <- enrichKEGG(
  gene = os.panel.df$entrezgene_id, 
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

dotplot(ego_kegg)


ego_wp <- enrichWP(os.panel.df$entrezgene_id, organism = "Homo sapiens") 
library(ReactomePA)
ego_re <- enrichPathway(os.panel.df$entrezgene_id, organism = "Homo sapiens", pvalueCutoff = 0.2)









# MSigDB ------------------------------------------------------------------


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

# os.panel.msig.result %>% 
#   dplyr::filter(gs_cat == "C2")

os.panel.msig.result %>% 
  dplyr::filter(gs_cat %in% c("C5")) %>% 
  dplyr::filter(grepl(pattern = "^GO", x = Description)) %>% 
  dplyr::mutate(type = substr(x = Description, start = 1, stop = 4)) %>% 
  dplyr::mutate(qvalue = -log10(p.adjust)) %>% 
  dplyr::arrange(dplyr::desc(type), -qvalue) ->
  c5

c5 %>% 
  ggplot(aes(y = Description)) +
  geom_col(aes( x = qvalue, fill = type)) +
  # geom_line(aes(x = Count, group = 1), color = "#B22222") +
  # scale_x_continuous(
  #   name = "-log10(p.adjust)",
  #   sec.axis = sec_axis(~., name="Count"),
  #   expand = c(0.02, 0)
  # ) 
  scale_y_discrete(
    limit = c5$Description,
  ) +
  scale_fill_viridis_d(
    label = c("BP", "CC", "MF"),
    name = "GO"
  ) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_line(color = "#B22222"),
    axis.text.y.right = element_text(color = "#B22222"),
    axis.title.y.right = element_text(color = "#B22222"),
    axis.ticks.y.right = element_line(color = "#B22222"),
    #legend.position = c(0.7, 0.7)
    legend.position = "top"
  ) +
  labs(y = "", x = "-log10(qvalue)") ->
  c5_go

ggsave(
  filename = "c5_go.pdf",
  plot = c5_go,
  device = "pdf",
  path = "data/review",
  width = 12,
  height = 9
)

os.panel.msig.result %>% 
  dplyr::filter(gs_cat %in% c("C2")) %>% 
  dplyr::filter(grepl(pattern = "^PID|^REACTOME", x = Description)) %>% 
  dplyr::mutate(type = substr(x = Description, start = 1, stop = 8)) %>% 
  dplyr::mutate(type = ifelse(grepl(pattern = "^PID", x = Description), "PID", type)) %>% 
  dplyr::mutate(qvalue = -log10(p.adjust)) %>% 
  dplyr::arrange(dplyr::desc(type), -qvalue) ->
  c2

c2 %>% 
  ggplot(aes(y = Description)) +
  geom_col(aes( x = qvalue, fill = type)) +
  scale_y_discrete(
    limit = c2$Description,
  ) +
  scale_fill_brewer(
    palette  = "Set2",
    label = c("PID", "REACTOME"),
    name = "Pathway"
  ) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_line(color = "#B22222"),
    axis.text.y.right = element_text(color = "#B22222"),
    axis.title.y.right = element_text(color = "#B22222"),
    axis.ticks.y.right = element_line(color = "#B22222"),
    #legend.position = c(0.7, 0.7)
    legend.position = "top"
  ) +
  labs(y = "", x = "-log10(qvalue)") ->
  c2_go
  
ggsave(
  filename = "c2_REACTOME.pdf",
  plot = c2_go,
  device = "pdf",
  path = "data/review",
  width = 12,
  height = 5
)


os.panel.msig.result %>%
  dplyr::filter(gs_cat %in% c("C6", "H", "C2", "C5"))->
  os.panel.msig.result.c6.h

# os.panel.msig.result.c6.h %>% 
  # dplyr::filter(grepl(pattern = 'platelet', x = Description, ignore.case = TRUE))

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
    # labels =stringr::str_wrap(
    #   stringr::str_replace_all(
    #     string = os.panel.msig.result.c6.h_ready$Description,
    #     pattern = "_",
    #     replacement = " "
    #   ), 
    #   width = 50)
  ) +
  scale_fill_viridis_d(
    label = c("C2: Canonical pathways", "H6: Oncogenic signature", "Hallmark"),
    name = "MSigDB"
  ) +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
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
