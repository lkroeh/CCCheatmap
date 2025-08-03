
# CCCheatmap

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of CCCheatmap is to build custom heatmaps from the output of CellChat.
These are modified versions of the heatmaps produced by their command \code{netVisual_heatmap}.
It allows the user to choose specific populations as senders and receivers for
which to consider signaling.
The colors are the sum of the communication probabilities for a given pathway 
between two designated cell types.


## Installation

You can install the development version of CCCheatmap from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lkroeh/CCCheatmap")
```

## Example

This is a basic example which shows you how to make a heatmap:

``` r
library(CCCheatmap)
cellchatobject
hm <- makehm(cellchatobject)
hm

#example color mapping
col_fun = colorRamp2(c(0, 0.0001, 0.01, 0.035, 2.5), c("grey", "lightyellow", "lightblue1","lightblue2", "blue"))
#The function will also output the quantiles of the data used in the heatmap, which will help you map the colors
```

