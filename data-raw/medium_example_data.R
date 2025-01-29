# This code generates a real example dataset from the Snodgrass dataset located at https://data.ess-dive.lbl.gov/view/doi:10.15485/2476543.

# Specify the spatial extent
a_x <- 326496
a_y <- 4311439
b_x <- 326596 - 50
b_y <- 4311539 - 50

medium_example_data <- readr::read_csv("https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/ess-dive-3120c69b6a46352-20240513T174234713", show_col_types = FALSE) %>%
  dplyr::filter(LCmajority == 1) %>%
  dplyr::select(XTOP, YTOP, CANVOL2015) %>%
  dplyr::filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  dplyr::rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  dplyr::arrange(desc(size)) |>
  dplyr::mutate(x = x - a_x,
                y = y - a_y) %>%
  as.data.frame()

usethis::use_data(medium_example_data, overwrite = TRUE)
