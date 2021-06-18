container <- E(combined_network)
graph <- combined_network
attribute <- "edge_type"

container[edge_type]
use_series(container, attribute)
edge_attr(combined_network, "edge_type_3")
is.na(use_series(container, attribute))
library(rlang)
expr <- stop("boom") %>%
    expr()
str(expr)
expr
eval(expr)
message("foo") %>%
    {
        message("bar")
        .
    }
expr <- rlang %>%
    library()
library(csdR)
data("sick_expression")
