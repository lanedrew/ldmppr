# This code generates the medium_example_data object from a 50 m x 50 m
# window in the full Snodgrass 2015 dataset:
# https://data.ess-dive.lbl.gov/view/doi:10.15485/2476543.

# Selected window
a_x <- 327308.19109940575
a_y <- 4311057.818756809
L <- 50
b_x <- a_x + L
b_y <- a_y + L

medium_example_data <- readr::read_csv("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-3120c69b6a46352-20240513T174234713", show_col_types = FALSE) %>%
  dplyr::select(XTOP, YTOP, CANVOL2015) %>%
  dplyr::filter(XTOP >= a_x & XTOP <= b_x & YTOP >= a_y & YTOP <= b_y) %>%
  dplyr::rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  dplyr::arrange(desc(size)) |>
  dplyr::mutate(
    x = x - a_x,
    y = y - a_y
  ) %>%
  as.data.frame()

usethis::use_data(medium_example_data, overwrite = TRUE)
