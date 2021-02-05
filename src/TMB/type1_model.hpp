
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
    Type two = 2.0;
    
    //input data
    DATA_INTEGER(n);
    DATA_INTEGER(nyrs);
    DATA_INTEGER(ndex);
    DATA_IVECTOR(iyear);
    DATA_IVECTOR(idex);
    DATA_VECTOR(pa);
    
    PARAMETER_VECTOR(iye);
    PARAMETER(lk);
    
    vector<Type> mu(nyrs);
    matrix<Type> new_mu(nyrs,ndex);
    
    //Observation Model
    int id, iy;
    for(int i = 0;i < n;++i){
      iy = iyear(i);
      id = idex(i);
      
      mu(iy) = exp(iye(iy))/(1+exp(iye(iy)));
      
      if(id==0){new_mu(iy,0)=mu(iy);}
      if(id>0){new_mu(iy,id)=mu(iy);}
      
      if(isNA(pa(i))==false){
        if(id==0){nll -= dbinom(pa(i), one, new_mu(iy,0), true);}
        if(id>0){nll -= dbinom(pa(i), one, new_mu(iy,id), true);}
      }
      
    }
    
    //RW on these likelihoods results in better estimates
    vector<Type> del_iye=iye;
    nll -= dnorm(del_iye(0),Type(10.0),one, true);
    for(int i = 1;i < nyrs;++i){
      nll -= dnorm(del_iye(i),del_iye(i-1),one, true);
    }
    
    
    REPORT(mu);
    REPORT(iye);
    REPORT(new_mu);
    
    ADREPORT(iye);
  
  return nll;
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
