
# Library -----------------------------------------------------------------

library(magrittr)


# Load data ---------------------------------------------------------------

d1 <- readxl::read_xlsx(path = "data/rda/platinum-os-pfs-used-data.xlsx", sheet = 1)
d2 <- readxl::read_xlsx(path = "data/rda/platinum-os-pfs-used-data.xlsx", sheet = 2)
d3 <- readxl::read_xlsx(path = "data/rda/platinum-os-pfs-used-data.xlsx", sheet = 3)

dplyr::bind_rows(
  d1 %>% dplyr::select(barcode:mapping_rate),
  d2 %>% dplyr::select(barcode:mapping_rate),
  d3 %>% dplyr::select(barcode:mapping_rate)
) %>% 
  dplyr::distinct() ->
  d

d %>% 
  dplyr::left_join(
    d1 %>% dplyr::select(barcode, platinum),
    by = "barcode"
  ) %>% 
  dplyr::left_join(
    d2 %>% dplyr::select(barcode, os_event = event, os_duration = duration),
    by = "barcode"
  ) %>% 
  dplyr::left_join(
    d3 %>% dplyr::select(barcode, pfs_event = event, pfs_duration = duration),
    by = "barcode"
  ) ->
  dd

writexl::write_xlsx(x = dd, path = "data/rda/merge-platinum-os-pfs-used-data.xlsx")


dd %>% 
  dplyr::mutate(group = dplyr::case_when(
    oc %in% c("OC521", "OC44") ~ "TC",
    oc == "OC79" ~ "VC1",
    oc == "OC172" ~ "VC2"
  )) ->
  ddd

ddd %>% 
  dplyr::group_by(group) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(iqr = purrr::map(.x = data, .f = function(.x) {
    .age <- .x$age
    .age %>% quantile() -> .q
    .all <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    .q <- .age[.age<=45] %>% quantile()
    .less45 <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    .q <- .age[.age>45] %>% quantile()
    .greater45 <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    tibble::tibble(
      name = c('Age', '<=45', '>45'),
      iqr = c(.all, .less45, .greater45)
    )
  })) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(cols = iqr) %>% 
  tidyr::spread(key = group, value = iqr)

ddd %>% 
  dplyr::group_by(group) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(iqr = purrr::map(.x = data, .f = function(.x) {
    .q <- .x$platelet_count %>% 
      as.numeric() %>% 
      na.omit() %>% 
      as.numeric() %>% 
      quantile()
    .all <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    tibble::tibble(
      name = "Platelet Count",
      iqr = .all
    )
  })) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(cols = iqr) %>% 
  tidyr::spread(key = group, value = iqr)


ddd %>% 
  dplyr::group_by(group) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(iqr = purrr::map(.x = data, .f = function(.x) {
    .q <- .x$os_duration %>% 
      as.numeric() %>% 
      na.omit() %>% 
      as.numeric() %>% 
      quantile()
    .all <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    tibble::tibble(
      name = "OS",
      iqr = .all
    )
  })) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(cols = iqr) %>% 
  tidyr::spread(key = group, value = iqr)

ddd %>% 
  dplyr::group_by(group) %>% 
  tidyr::nest() %>% 
  dplyr::mutate(iqr = purrr::map(.x = data, .f = function(.x) {
    .q <- .x$pfs_duration %>% 
      as.numeric() %>% 
      na.omit() %>% 
      as.numeric() %>% 
      quantile()
    .all <- glue::glue('{.q[3]} ({.q[2]}-{.q[4]})')
    tibble::tibble(
      name = "OS",
      iqr = .all
    )
  })) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(cols = iqr) %>% 
  tidyr::spread(key = group, value = iqr)


# Save image --------------------------------------------------------------

save.image(file = "data/rda/data-merge.rda")
