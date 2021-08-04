
// **DON'T** #include <TMB.hpp> as it is not include-guarded
//#include <LFPM.hpp> 

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type> //only list this in whichever cpp file is alphabetically first since thats the first cpp file read in and the others aren't include guarded meaning it will tell you that it has been re-defined
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// name of function below **probST** match filename
template <class Type>
Type LFPM(objective_function<Type>* obj) {

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
    DATA_VECTOR(pa);
    
    //READ THE PARAMETERS
    PARAMETER_VECTOR(iye);
    PARAMETER(lk); //CHECK CHECK PLACEHOLDER PARAMETER BECAUSE CAN'T ONLY HAVE RANDOM EFFECTS
    
    //CREATE DERIVED PARAMETERS
    vector<Type> mu(nyrs);
    matrix<Type> new_mu(nyrs,ndex);
    
    int id, iy;
    for(int i = 0;i < n;++i){
      iy = iyear(i);
      id = idex(i);
      
      mu(iy) = one-exp(-iye(iy)); //EXPONENTIAL FOR TRAWL DATA
      
      new_mu(iy,id)=mu(iy); 
      
      
      if(isNA(pa(i))==false){
        if(id==0){nll -= dbinom(pa(i), one, new_mu(iy,0), true);} //BERNOULLI LIKELIHOOD FOR PRESENCE/ABSENCE DATA
        if(id>0){nll -= dbinom(pa(i), one, new_mu(iy,id), true);} 
      }
      
    }
    
    //GAUSSIAN RW FOR THE ABUNDANCE INDEX
    vector<Type> del_iye=log(iye);
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
