tastyfish: a package to account for a non-linear functional response
when estimating prey dynamics using predator diet data
================
Matthew Robertson
04/8/2021

# Overview

<img src="hexsticker/hexsticker_tastyfish.PNG" style="width:35.0%" />

# Installation

Install the
[devtools](https://cran.r-project.org/web/packages/devtools/index.html)
package and run:

``` r
devtools::install_github(MatthewRobertson2452\tastyfish)
library(tastyfish)
```

If you are having problems with installation, you can install the
package locally as a ZIP file by clicking the Code menu and “download
ZIP” from the [github
page](https://github.com/MatthewRobertson2452/tastyfish). You can then
extract the folder in a local directory while recording the directory
name (which I will reference as download\_dir). To install, then use

``` r
devtools::install_local(path=download_dir, dependencies=FALSE)
library(tastyfish)
```

# Example

This example will run through how an index was developed for Northern
sand lance using stomach contents data from American plaice.

For this example, I will be analyzing trawl, full stomach contents, and
called stomach contents data. They are all formatted similarly

``` r
head(trawl)
```

    ##   year pa
    ## 1 1996  1
    ## 2 1996  0
    ## 3 1996  0
    ## 4 1996  1
    ## 5 1996  1
    ## 6 1996  0

We then need to organize all three types of data into one dataframe to
develop the data format needed for analyses. This is done by creating a
vector that describes the length of each dataset and then inputting the
data and that vector (n) into `make_data()`.

``` r
n<-c(length(trawl$year), length(call_sto$year), length(full_sto$year))

fishy_dat<-make_data(pa=c(trawl$pa, call_sto$pa, full_sto$pa),
          year=c(trawl$year, call_sto$year, full_sto$year),
          n=n)
```

Since both called and full stomach contents data will be treated as
having the same functional response we will identify them with the same
ID for the model. This ID will be different from the trawl data. The
order of the numbers in the `id` vector need to match the order of the
dataframes input in `make_data()`.

``` r
id=c(0,1,1)
```

If we wanted the different types of stomach contents data to estimate
separate functional response shapes we would use:

``` r
id=c(0,1,2)
```

Now that we have that pre-processing done, we can input that information
with some starting values to create the tmb data and parameter lists
using `model_data` and `model_param`. If you are estimating multiple
functional response forms (e.g. from multiple predators), lbeta and lchi
can be treated as vectors where the order of starting values will match
the order that the data was input as.

``` r
tmb_data<-model_data(fishy_dat=fishy_dat, n=n, type="nonlinear")

param_list<- model_param(tmb_data, lbeta=log(2), lchi=log(0.5), type="nonlinear")
```

We can then run the TMB model to calculate the prey index of abundance
while also estimating and accouting for the functional response of the
predator by using the following code:

``` r
obj<- TMB::MakeADFun(data = c(model = "NLFPM", # which model to use
                                tmb_data),
                       parameters = param_list,
                       DLL = "tastyfish_TMBExports", 
                       random=c("iye")) # package's DLL

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)


rep<-obj$report()
sdrep<-TMB::sdreport(obj)
```

I have also included code to quickly visualize the outputs of the model,
for the functional response:

``` r
plot_curve(rep, tmb_data)
```

![Estimated functional response for American plaice. Full stomach
contents shown as open circles and called stomach contents shown as
closed circles.](README_files/figure-gfm/unnamed-chunk-10-1.png)

And time-series:

``` r
plot_ts(tmb_data, rep, sdrep, year=unique(spring_camp_surv$year), dat_names=c("Trawl","Called","Full"), type="nonlinear")
```

![Estimated sand lance index of abundance from models with American
plaice stomach contents data. The shaded grey area represents the
standard error around the estimated trend. Dashed lines represent the
estimated trends from each data
source.](README_files/figure-gfm/unnamed-chunk-11-1.png)

We can compare these results to what we would obtain if assuming a
linear functional response by running the LFPM. We simply remove the
nonlinear functional response shape parameters from the `param_list`
function and change the type to “linear”:

``` r
tmb_data<-model_data(fishy_dat=fishy_dat, n=n, type="linear")

param_list<- model_param(tmb_data, type="linear")
```

We then change the TMB model to “LFPM”

``` r
obj<- TMB::MakeADFun(data = c(model = "LFPM", # which model to use
                                tmb_data),
                       parameters = param_list,
                       DLL = "tastyfish_TMBExports", 
                       random=c("iye")) # package's DLL

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)


rep<-obj$report()
sdrep<-TMB::sdreport(obj)
```

Finally, we can plot our output be specifying that we used a linear
functional response model in `plot_ts` function.

``` r
plot_ts(tmb_data, rep, sdrep, year=unique(spring_camp_surv$year), dat_names=c("Trawl","Called","Full"), type="linear")
```

![Estimated sand lance index of abundance from models with American
plaice stomach contents data. The shaded grey area represents the
standard error around the estimated trend. Dashed lines represent the
estimated trends from each data
source.](README_files/figure-gfm/unnamed-chunk-14-1.png)
