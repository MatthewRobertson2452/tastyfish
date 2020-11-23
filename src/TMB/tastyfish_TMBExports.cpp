// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_tastyfish_TMBExports
#include <TMB.hpp>
#include "matrix_model_new.hpp"
#include "type1_model.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "matrix_model_new") {
    return matrix_model_new(this);
  } else if(model == "type1_model") {
    return type1_model(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}