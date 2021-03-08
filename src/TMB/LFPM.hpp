
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
    PARAMETER_VECTOR(ny);
    PARAMETER(lk); //PLACEHOLDER PARAMETER BECAUSE CAN'T ONLY HAVE RANDOM EFFECTS
    
    //TRANSFORM PARAMETERS AND CREATE DERIVED PARAMETERS
    vector<Type> prob(nyrs);
    matrix<Type> new_prob(nyrs,ndex);
    
    int id, iy;
    for(int i = 0;i < n;++i){
      iy = iyear(i);
      id = idex(i);
      
      prob(iy) = exp(ny(iy))/(1+exp(ny(iy))); //LOGIT TRANSFORM THE ABUNDANCE INDEX
      
      new_prob(iy,id)=prob(iy); // ALL UNBIASED DATA
      
      
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
    REPORT(new_prob);
    
    ADREPORT(ny);
  
  return nll;
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
