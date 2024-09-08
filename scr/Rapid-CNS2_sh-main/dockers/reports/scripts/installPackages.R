for (package in c('optparse', 'rmarkdown','kableExtra','knitr', 'ggplot2', 'openxlsx')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos = "http://cran.us.r-project.org")
    library(package, character.only=T)
  }
}