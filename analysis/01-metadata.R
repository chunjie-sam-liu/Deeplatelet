
# Library -----------------------------------------------------------------

library(magrittr)


# Path --------------------------------------------------------------------


filepath <- 'data/metadata/血小板预后2020年11月V1chijh.xlsx'


# Load data ---------------------------------------------------------------

metadata <- readxl::read_excel(path = filepath, sheet = 6)


