

#' Model selection procedure for real data
#'
#' @param tmb_data TMB data object from the model_data() function
#' @param param_list Parameter list object from the model_param() function
#' @param k_seq Sequence of all possible k values to be cycled through in the model selection procedure, default is seq(from=0.1,to=1,by=0.05)
#' @param chi_seq Sequence of all possible chi values to be cycled through in the model selection procedure, default is seq(from=0.1,to=1,by=0.05)
#'
#' @return A list for outputs from the selected model including, its AIC, as well as the k, chi, and beta parameter estimates
#' @export 
#'
#' @examples 
#' best_model<-sel_proc(tmb_data=tmb_data, param_list=param_list, k_seq=seq(from=0.1,to=1,by=0.05), chi_seq=seq(from=0.05,to=1,by=0.05))
sel_proc<-function(tmb_data, param_list,
                   k_seq=seq(from=0.1,to=1,by=0.05), chi_seq=seq(from=0.05,to=1,by=0.05)){

  
grad_list<-matrix(NA, ncol=length(chi_seq), nrow=length(k_seq))
type2_given_aic<-matrix(NA, ncol=length(chi_seq), nrow=length(k_seq))

for(i in 1:length(k_seq)){
  for(j in 1:length(chi_seq)){
    
    tmb_data$k<-k_seq[i]
    tmb_data$chi<-chi_seq[j]
    
    obj<- TMB::MakeADFun(data = c(model = "matrix_model_new", # which model to use
                                  tmb_data),
                         parameters = param_list, inner.control=list(maxit=2000,trace=F),
                         DLL = "tastyfish_TMBExports", 
                         random=c("yr_tau_mat","iye"), silent=TRUE) # package's DLL
    ###CHECK FOR ERRORS IF MODELS CREATE NA
    x<-NULL
    x<-tryCatch(expr=nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE),
                error = function(e) e)
    
    x<-tryCatch(expr=nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE),
                error = function(e) e)
    
    x<-tryCatch(expr=nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE),
                error = function(e) e)
    
    if(inherits(x, "simpleError")){ next 
    }else{
      
      opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
      
      grad_list[i,j]<-max(obj$gr())
      
      type2_given_aic[i,j]<-2*opt$objective+ 2*length(opt$par)#aic

    }
  }
}



big_grad<-which(grad_list>0.01|grad_list<0.0000000001, arr.ind=TRUE)

for(i in 1:length(big_grad[,1])){
  type2_given_aic[big_grad[i,1],big_grad[i,2]]<-NA
}

type2_given_aic[is.nan(grad_list)]<-NA

small_aic<-which(type2_given_aic<(min(type2_given_aic, na.rm=TRUE)+20), arr.ind=TRUE)


aic_vals<-matrix(NA, nrow=length(small_aic[,1]), ncol=1)
for(i in 1:length(small_aic[,1])){
  aic_vals[i,1]<-type2_given_aic[small_aic[i,1],small_aic[i,2]]
}

sel_k<-k_seq[small_aic[which.min(aic_vals),1]]
sel_chi<-chi_seq[small_aic[which.min(aic_vals),2]]


tmb_data$k<-sel_k
tmb_data$chi<-sel_chi

obj<- TMB::MakeADFun(data = c(model = "matrix_model_new", # which model to use
                              tmb_data),
                     parameters = param_list,inner.control=list(maxit=2000,trace=F),
                     DLL = "tastyfish_TMBExports", 
                     random=c("yr_tau_mat","iye"), silent=TRUE) # package's DLL

opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)

rep<-obj$report()

output_list<-list(
  aic = 2*opt$objective+ 2*length(opt$par),
  k = sel_k,
  chi = sel_chi,
  beta = rep$beta
)

return(output_list)
}