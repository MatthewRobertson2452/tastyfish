
// **DON'T** #include <TMB.hpp> as it is not include-guarded
//#include <NLFPM.hpp> 

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj



// name of function below **probST** match filename
// (here it's three_spp_model)
template <class Type>
Type NLFPM(objective_function<Type>* obj) {

  using namespace density;
  Type nll = 0.0;
  Type zero = 0.0;
  Type one = 1.0;
  Type two = 2.0;
  
  //READ IN THE INPUT DATA
  DATA_INTEGER(n);
  DATA_INTEGER(nyrs);
  DATA_INTEGER(ndex);
  DATA_IVECTOR(iyear);
  DATA_IVECTOR(idex);
  DATA_SCALAR(k);
  DATA_VECTOR(pa);
  
  //READ THE PARAMETERS
  PARAMETER_VECTOR(iye);
  PARAMETER_VECTOR(lbeta);
  PARAMETER_VECTOR(lchi);
  PARAMETER(logrw_var);
  
  //TRANSFORM PARAMETERS AND CREATE DERIVED PARAMETERS
  Type rw_var = exp(logrw_var);
  vector<Type> beta=exp(lbeta);
  vector<Type> chi=exp(lchi);
  
  vector<Type> mu(nyrs);
  matrix<Type> new_mu(nyrs,ndex);
  
  int id, iy;
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    
    mu(iy) = one-exp(-iye(iy));  //EXPONENTIAL FOR TRAWL DATA
    
    if(id==0){new_mu(iy,0)=mu(iy);} //TRAWL PROBABILITY UNBIASED
    if(id>0){new_mu(iy,id)=(k*pow(mu(iy),beta(id-1)))/(pow(chi(id-1),beta(id-1))+pow(mu(iy),beta(id-1)));} //FUNCTIONAL RESPONSE FOR EACH STOMACH CONTENT DATA SET
    
    if(isNA(pa(i))==false){
      if(id==0){nll -= dbinom(pa(i), one, new_mu(iy,0), true);} //BERNOULLI LIKELIHOOD FOR PRESENCE/ABSENCE DATA
      if(id>0){nll -= dbinom(pa(i), one, new_mu(iy,id), true);}
    }
    
  }
  
  //GAUSSIAN RW FOR THE ABUNDANCE INDEX
  vector<Type> del_iye=log(iye);
  nll -= dnorm(del_iye(0),Type(10.0),rw_var, true);
  for(int i = 1;i < nyrs;++i){
    nll -= dnorm(del_iye(i),del_iye(i-1),rw_var, true);
  }
  
  
  REPORT(rw_var);
  REPORT(mu);
  REPORT(iye);
  REPORT(beta);
  REPORT(chi);
  REPORT(new_mu);
  
  ADREPORT(iye);
  
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
