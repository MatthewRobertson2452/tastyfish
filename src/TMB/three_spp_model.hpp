/// @file ModelA.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
// (here it's ModelA)
template <class Type>
Type three_spp_model(objective_function<Type>* obj) {

  using namespace density;
  Type nll = 0.0;
  Type zero = 0.0;
  Type one = 1.0;
  //Type two = 2.0;
  //Type three = 3.0;

  //input data
  DATA_INTEGER(n);
  DATA_INTEGER(nyrs);
  //DATA_INTEGER(nyrs_p1);
  //DATA_INTEGER(nyrs_p2);
  DATA_INTEGER(ndex);
  DATA_IVECTOR(iyear);
  //DATA_IVECTOR(iyear_p1);
  //DATA_IVECTOR(iyear_p2);
  DATA_IVECTOR(idex);
  DATA_SCALAR(k);
  DATA_SCALAR(chi);
  DATA_VECTOR(pa);
  DATA_IVECTOR(t2dex);

  PARAMETER_VECTOR(yr_tau);
  PARAMETER_VECTOR(yr_tau2);
  PARAMETER_VECTOR(yr_tau3);
  PARAMETER_VECTOR(iye);
  //PARAMETER_VECTOR(log_sd_resid);
  PARAMETER(lbeta);

  //vector<Type>intercept(ndex);
  Type beta=exp(lbeta);
  matrix<Type> resid(nyrs,ndex);
  //vector<Type> resid2(nyrs);
  //vector<Type> resid3(nyrs);
  //vector<Type> sd_resid = exp(log_sd_resid);
  vector<Type> mu(nyrs);
  vector<Type> mu2(nyrs);
  vector<Type> mu3(nyrs);
  //Type intercept;

  //resid = zero;

  //Observation Model
  int id, iy, t2;
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    //iy1 = iyear_p1(i);
    //iy2 = iyear_p2(i);
    id = idex(i);

    if(id==0){mu(iy) = exp(yr_tau(iy))/(1+exp(yr_tau(iy)));}
    if(id==1){mu2(iy) = exp(yr_tau2(iy))/(1+exp(yr_tau2(iy)));}
    if(id==2){mu3(iy) =  exp(yr_tau3(iy))/(1+exp(yr_tau3(iy)));}

    if(isNA(pa(i))==false){
      if(id==0){nll -= dbinom(pa(i), one, mu(iy), true);}
      if(id==1){nll -= dbinom(pa(i), one, mu2(iy), true);}
      if(id==2){ nll -= dbinom(pa(i), one, mu3(iy), true);}
    }

  }


  vector<Type> del_yr_tau=yr_tau;
  vector<Type> del_yr_tau2=yr_tau2;
  vector<Type> del_yr_tau3=yr_tau3;
  nll -= dnorm(del_yr_tau(0),Type(0.5),one, true);
  nll -= dnorm(del_yr_tau2(0),Type(0.5),one, true);
  nll -= dnorm(del_yr_tau3(0),Type(0.5),one, true);
  for(int i = 1;i < nyrs;++i){
    nll -= dnorm(del_yr_tau(i),del_yr_tau(i-1),one, true);
    nll -= dnorm(del_yr_tau2(i),del_yr_tau2(i-1),one, true);
    nll -= dnorm(del_yr_tau3(i),del_yr_tau3(i-1),one, true);
  }
  //second_mu=mu2


  //third_mu=new_mu.col(2);

  vector<Type> new_mu(nyrs);
  vector<Type> second_mu(nyrs);
  vector<Type> third_mu(nyrs);



  //Process Model //can examine yr_tau to eval difference in accuracy bwn Type 1 and Type 2 functional response assumption
  for(int i = 0;i < n;++i){
    iy = iyear(i);
    id = idex(i);
    t2 = t2dex(i);
    //iy1 = iyear_p1(i);
    //iy2 = iyear_p2(i);
    new_mu(iy)=mu(iy);
    second_mu(iy)=pow((-((mu2(iy)-k)*pow(chi,-beta))/mu2(iy)),-1/beta);
    third_mu(iy)=pow((-((mu3(iy)-k)*pow(chi,-beta))/mu3(iy)),-1/beta);

    if(t2 == zero){

      resid(iy,id) =  new_mu(iy) - iye(iy);} //for the trawl samples we assume they are proportional to pop size, so can just subtract values directly to get year effect

    if(t2 == one){
      if(id==1){

        //if(second_mu(iy)==zero){resid2(iy)==zero;}else{
        resid(iy,id) = second_mu(iy) - iye(iy);}//}
      if(id==2){

        //if(third_mu(iy)==zero){resid3(iy)==zero;}else{
        resid(iy,id) = third_mu(iy) - iye(iy);} //}
      //}//for the predation samples we know they follow type 2 so need to find the difference between type 1 and type 2 estimates and then can get year effect
    }
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
    //}
  }


  REPORT(resid);
  //REPORT(resid2);
  //REPORT(resid3);
  REPORT(yr_tau);
  REPORT(yr_tau2);
  REPORT(yr_tau3);
  REPORT(mu);
  REPORT(mu2);
  REPORT(mu3);
  REPORT(iye);
  REPORT(beta);
  REPORT(new_mu);
  REPORT(second_mu);
  REPORT(third_mu);

  ADREPORT(iye);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
