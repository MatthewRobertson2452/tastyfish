// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_tastyfish_TMBExports
#include <TMB.hpp>
#include "three_spp_model.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "three_spp_model") {
    return three_spp_model(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
