
# Library -----------------------------------------------------------------

library(magrittr)


# Path --------------------------------------------------------------------


filepath <- 'data/metadata/血小板预后2020年11月V1chijh.xlsx'


# Load data ---------------------------------------------------------------

metadata <- readxl::read_excel(path = filepath, sheet = 6, skip = 1)

metadata %>% 
  dplyr::select(barcode, patient_ID, oc, time_on_diagnosis, end_of_chemotherapy, platinum_sensitivity, palindromia, palindromia2, time_of_palindromia, time_of_palindromia2, pfs, alive_type, alive_type2, time_of_die, time_of_followup, os)


