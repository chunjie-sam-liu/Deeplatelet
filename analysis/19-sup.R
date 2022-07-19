# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Jun  5 04:09:12 2022
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)

# Load data ---------------------------------------------------------------

npat <- readxl::read_excel(path = "data/supdata/real-data.xlsx", sheet = 1)
