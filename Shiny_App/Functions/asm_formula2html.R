formula2html <- function(formula) {
  gsub("([0-9]+)", "<sub>\\1</sub>", formula)
}