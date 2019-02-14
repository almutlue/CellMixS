## Code to generate color palettes used
library(jcolors)

# histogram colors
col_hist <- c(jcolors('pal6'),jcolors('pal8')[c(10,6,12)])
names(colors_j) <- c()

#group colors
col_group <-c(jcolors('pal6'),jcolors('pal8')[c(10,6,12)])[c(1,8,10,5,2:4,6,7,9,11)]
names(col_group) <- c()

#save as internal data
usethis::use_data(col_hist, col_group, internal = TRUE, overwrite = TRUE)
