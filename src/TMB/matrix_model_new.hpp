
// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// name of function below **probST** match filename
// (here it's three_spp_model)
template <class Type>
Type matrix_model_new(objective_function<Type>* obj) {

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
  DATA_IVECTOR(ifxn);
  DATA_SCALAR(k);
  DATA_VECTOR(pa);
  
  //READ THE PARAMETERS
  PARAMETER_VECTOR(ny);
  PARAMETER_VECTOR(lbeta);
  PARAMETER_VECTOR(lchi);
  
  //TRANSFORM PARAMETERS AND CREATE DERIVED PARAMETERS
  vector<Type> beta=exp(lbeta);
  vector<Type> chi=exp(lchi);
  vector<Type> prob(nyrs);
  matrix<Type> new_prob(nyrs,ndex);
  
  int id, iy, iw;
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    iw = ifxn(i);
    
    prob(iy) = exp(ny(iy))/(1+exp(ny(iy))); //LOGIT TRANSFORM THE ABUNDANCE INDEX
    
    if(id==0){new_prob(iy,0)=prob(iy);} //TRAWL PROBABILITY UNBIASED
    if(id>0){new_prob(iy,id)=(k*pow(prob(iy),beta(iw)))/(pow(chi(iw),beta(iw))+pow(prob(iy),beta(iw)));} //FUNCTIONAL RESPONSE FOR EACH STOMACH CONTENT DATA SET
    
    if(isNA(pa(i))==false){
      nll -= dbinom(pa(i), one, new_prob(iy,id), true); //BERNOULLI LIKELIHOOD FOR PRESENCE/ABSENCE DATA
    }
    
  }
  
  //GAUSSIAN RW FOR THE ABUNDANCE INDEX
  vector<Type> del_ny=ny;
  nll -= dnorm(del_ny(0),Type(10.0),one, true);
  for(int i = 1;i < nyrs;++i){
    nll -= dnorm(del_ny(i),del_ny(i-1),one, true);
  }
  
  
  REPORT(prob);
  REPORT(ny);
  REPORT(beta);
  REPORT(chi);
  REPORT(new_prob);
  
  ADREPORT(ny);
  
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
