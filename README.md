tastyfish: a package to develop an index of prey abundance from stomach
contents data while accounting for predator functional response
================
Matt Robertson
23/11/2020

# Overview

# Installation

Install the
[devtools](https://cran.r-project.org/web/packages/devtools/index.html)
package and run:

``` r
devtools::install_github(MatthewRobertson2452\tastyfish)
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

Now that we have that pre-processing done, we can input that information
with some starting values to create the tmb data and parameter lists
using `model_data` and `model_param`.

``` r
tmb_data<-model_data(fishy_dat=fishy_dat, k=0.5, chi=0.1, id=id, n=n)

param_list<- model_param(tmb_data, lbeta=log(5))
```

We can now input these lists into `sel_proc` to run the model selection
procedure that is described in the paper.

``` r
best_model<-sel_proc(tmb_data, param_list)
```

This will provide us with information about the best model that we can
run as its own model after updating the data list by inputting the k and
chi values from the model selection procedure. We can then run the TMB
model by using the following code:

``` r
tmb_data<-model_data(fishy_dat=fishy_dat, k=best_model$k, chi=best_model$chi, id=id, n=n)

obj<- TMB::MakeADFun(data = c(model = "matrix_model_new", # which model to use
                                tmb_data),
                       parameters = param_list,
                       DLL = "tastyfish_TMBExports", 
                       random=c("yr_tau_mat","iye")) # package's DLL

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
closed circles.](README_files/figure-gfm/unnamed-chunk-9-1.png)

And time-series:

``` r
plot_ts(rep, sdrep, unique(trawl$year), dat_names=c("Trawl","Call","Full"), ylim=c(0.3,0.6))
```

![Estimated sand lance index of abundance from models with American
plaice stomach contents data. The shaded grey area represents the
standard error around the estimated trend. Dashed lines represent the
estimated trends from each data
source.](README_files/figure-gfm/unnamed-chunk-10-1.png)
