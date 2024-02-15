# DualSimplex algorithm's R package

## Installation
We haven't yet defined, which functions should be exported
from the package and which ones are internal, so for now, it is:
```
devtools::load_all()
```

After we define that, it will be: 
```
devtools::install_github("artyomovlab/DualSimplex")
```

After the publication, it will be:
```
install.packages("DualSimplex")
```


## Usage
A minimal usage example will be provided here.

For now, see Rmd examples in `inst/examples`.


## Code structure & Guidelines

The following files in `R/` directory represent different stages
of DualSimplex deconvolution process:
```
0. simulation.R
1. annotation.R
2. filtering.R
3. sinkhorn.R
4. projection.R
5. initialization.R
6. optimization.R
7. post_analysis.R
8. benchmarking.R
```

Ideally, a main logic functions in a stage shouldn't use 
functions from another stage, and a downstream stage 
should only use the objects generated on the previous stage as its input. 

Then, either the user or `DualSimplexSolver` use the main
functions from those packages to implement the whole control flow.

This rule of thumb leads to linear code logic and low code coupling,
which makes it simple to debug and introduce changes. However, there are
a few exceptions even in the current code, so we'll see how that works out.

> Note: for now there is a little bit too much functions, and there is
some redundancy between DualSimplexSolver and functions. Later we'll finalize
DualSimplexSolver and create some examples on using functions only. 
After that we'll determine, which of the functions will be exported
from the package, and remove unnecessary and unused ones.
It will be very easy to do it all, the code is quite decoupled.


### Current exceptions
Functions from `utils.R` use data from multiple stages, so they are separate.

Functions from `projected_plots.R` combine multiple stages too, because
those are the main plots and it is more convenient to keep the code together.

The `figures.R` file contains code only used for figure creation for the paper.
