
// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
  Type type1_model(objective_function<Type>* obj) {

  using namespace density;
  Type nll = 0.0;
  Type zero = 0.0;
  Type one = 1.0;
  
  //input data
  DATA_INTEGER(n);
  DATA_INTEGER(nyrs);
  DATA_INTEGER(ndex);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(idex);
  DATA_VECTOR(pa);
  
  PARAMETER_MATRIX(yr_tau);
  PARAMETER_VECTOR(iye);
  PARAMETER(log_sd_resid);
  
  matrix<Type> resid(nyrs,ndex);
  Type sd_resid = exp(log_sd_resid);
  matrix<Type> new_mu(nyrs,ndex);
  matrix<Type> mu(nyrs,ndex);
  
  
  
  //Observation Model
  int id, iy;
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    
    mu(iy,id) = exp(yr_tau(iy, id))/(1+exp(yr_tau(iy, id)));
    
    if(isNA(pa(i))==false){
      nll -= dbinom(pa(i), one, mu(iy,id), true); 
    }
  }
  
  
  //Process Model //can examine yr_tau to eval difference in accuracy bwn Type 1 and Type 2 functional response assumption
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    
    resid(iy,id) =  mu(iy,id) - iye(iy);} //for the trawl samples we assume they are proportional to pop size, so can just subtract values directly to get year effect
  
  //RW on these likelihoods results in better estimates
  vector<Type> del_iye=iye;
  nll -= dnorm(del_iye(0),Type(0.5),one, true);
  for(int i = 1;i < nyrs;++i){
    nll -= dnorm(del_iye(i),del_iye(i-1),one, true);
  }
  
  for(int j = 0;j < ndex;++j){
    for(int i = 0;i < nyrs;++i){
      nll -= dnorm(resid(i,j), zero, sd_resid, true); 
    }
  }
  
  
  REPORT(resid);
  REPORT(yr_tau);
  REPORT(mu);
  REPORT(iye);
  
  ADREPORT(iye);
  
  return nll;
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
