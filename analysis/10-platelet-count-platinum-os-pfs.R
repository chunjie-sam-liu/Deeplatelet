
# Library -----------------------------------------------------------------


library(magrittr)
library(ggplot2)

# Load data ---------------------------------------------------------------

total351.platinum.se <- readr::read_rds(file = 'data/rda/total351.platinum.se.rds.gz')

metadata <- as.data.frame(total351.platinum.se@colData) %>% tibble::as_tibble()



# Function ----------------------------------------------------------------


# Analysis ----------------------------------------------------------------


# Platinum platelet count -------------------------------------------------


metadata %>% 
  dplyr::mutate(platelet_count = as.numeric(platelet_count)) %>% 
  dplyr::filter(!is.na(platelet_count)) %>%  
  ggpubr::ggboxplot(
    x = 'platinum',
    y = 'platelet_count',
    color = 'platinum',
    xlab = 'Platinum',
    ylab = 'Platelet count (in 10^9/L)',
    title = "Platinum sensitivity correlation with platelet count"
  ) +
  ggpubr::stat_compare_means() +
  ggthemes::theme_few() ->
  plot_platinum_platelet_count
ggsave(
  filename = 'platinum-platelet-count-correlation.pdf',
  plot = plot_platinum_platelet_count,
  device = 'pdf',
  path = 'data/output/',
  width = 6.5,
  height = 4.2
)


# platinum age ------------------------------------------------------------


metadata %>% 
  ggpubr::ggboxplot(
    x = 'platinum',
    y = 'age',
    color = 'platinum',
    xlab = 'Platinum',
    ylab = 'Age',
    title = 'Platinum sensitivity correlation with age'
  ) +
  ggpubr::stat_compare_means() +
  ggthemes::theme_few() ->
  plot_platinum_age

ggsave(
  filename = 'platinum-age-correlation.pdf',
  plot = plot_platinum_age,
  device = 'pdf',
  path = 'data/output/',
  width = 6.5,
  height = 4.2
)





# Platelet count os and pfs -----------------------------------------------


