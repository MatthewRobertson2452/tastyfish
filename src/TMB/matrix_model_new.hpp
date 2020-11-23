
// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// name of function below **MUST** match filename
// (here it's three_spp_model)
template <class Type>
Type matrix_model_new(objective_function<Type>* obj) {

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
  DATA_SCALAR(k);
  DATA_SCALAR(chi);
  DATA_VECTOR(pa);
  DATA_IVECTOR(t2dex);
  
  PARAMETER_MATRIX(yr_tau_mat);
  PARAMETER_VECTOR(iye);
  PARAMETER(lbeta);
  
  
  Type beta=exp(lbeta);
  matrix<Type> resid(nyrs,ndex);
  matrix<Type> mu(nyrs, ndex);
  
  //Observation Model
  int id, iy, t2;
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    
    if(id==0){mu(iy,0) = exp(yr_tau_mat(iy,0))/(1+exp(yr_tau_mat(iy,0)));}
    if(id==1){mu(iy,1) = exp(yr_tau_mat(iy,1))/(1+exp(yr_tau_mat(iy,1)));}
    if(id==2){mu(iy,2) = exp(yr_tau_mat(iy,2))/(1+exp(yr_tau_mat(iy,2)));}
    
    if(isNA(pa(i))==false){
      if(id==0){nll -= dbinom(pa(i), one, mu(iy,0), true);}
      if(id==1){nll -= dbinom(pa(i), one, mu(iy,1), true);}
      if(id==2){nll -= dbinom(pa(i), one, mu(iy,2), true);}
    }
    
  }
  
  
  matrix<Type> del_yr_tau=yr_tau_mat;
  for(int j = 0;j < ndex;++j){
    nll -= dnorm(del_yr_tau(0,j),Type(0.5),one, true);
    for(int i = 1;i < nyrs;++i){
      nll -= dnorm(del_yr_tau(i,j),del_yr_tau(i-1,j),one, true);
    }}
  
  matrix<Type> new_mu(nyrs,ndex);
  
  
  //Process Model //can examine yr_tau to eval difference in accuracy bwn Type 1 and Type 2 functional response assumption
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    t2 = t2dex(i);
    
    new_mu(iy,0)=mu(iy,0);
    if(id==1){new_mu(iy,1)=pow((-((mu(iy,1)-k)*pow(chi,-beta))/mu(iy,1)),-1/beta);}
    if(id==2){new_mu(iy,2)=pow((-((mu(iy,2)-k)*pow(chi,-beta))/mu(iy,2)),-1/beta);}
    
    if(t2 == zero){
      
      resid(iy,id) =  new_mu(iy,0) - iye(iy);} //for the trawl samples we assume they are proportional to pop size, so can just subtract values directly to get year effect
    
    if(t2 == one){
      
      resid(iy,id) = new_mu(iy,id) - iye(iy);}
    
  }
  
  
  
  //RW on these likelihoods results in better estimates
  vector<Type> del_iye=iye;
  nll -= dnorm(del_iye(0),Type(0.5),one, true);
  for(int i = 1;i < nyrs;++i){
    nll -= dnorm(del_iye(i),del_iye(i-1),one, true);
  }
  
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    t2 = t2dex(i);
    nll -= dnorm(resid(iy,id), zero, one, true);
  }
  
  
  REPORT(resid);
  REPORT(yr_tau_mat);
  REPORT(mu);
  REPORT(iye);
  REPORT(beta);
  REPORT(new_mu);
  
  ADREPORT(iye);
  
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
