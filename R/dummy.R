#
# A new R CMD CHECK Note was introduced when I upgraded to R 4.2.0
# 07.06.2024
# A dummy function to pass "checking dependencies in R code...
# Namespaces in Imports field not imported from...
#
dummy <- function() {
  GGally::brew_colors
  bookdown::clean_book
  grid::arrow
  kableExtra::footnote
#  mcmcplots::caterplot
  qpdf::pdf_subset
}
