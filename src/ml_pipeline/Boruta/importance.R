# Importance sources

#' randomForest importance adapters
#'
#' Those function is intended to be given to a \code{getImp} argument of \code{\link{Boruta}} function to be called by the Boruta algorithm as an importance source.
#' \code{getImpLegacyRfZ} generates default, normalized permutation importance, \code{getImpLegacyRfRaw} raw permutation importance, finally \code{getImpLegacyRfGini} generates Gini index importance, all using \code{\link[randomForest]{randomForest}} as a Random Forest algorithm implementation.
#' @name getImpLegacyRf
#' @rdname getImpLegacyRf
#' @aliases getImpLegacyRfZ getImpLegacyRfGini getLegacyImpRfRaw
#' @note The \code{getImpLegacyRfZ} function was a default importance source in Boruta versions prior to 5.0; since then \code{\link{ranger}} Random Forest implementation is used instead of \code{\link[randomForest]{randomForest}}, for speed, memory conservation and an ability to utilise multithreading.
#' Both importance sources should generally lead to the same results, yet there are differences.
#'
#' Most notably, ranger by default treats factor attributes as ordered (and works very slow if instructed otherwise with \code{respect.unordered.factors=TRUE}); on the other hand it lifts 32 levels limit specific to \code{\link[randomForest]{randomForest}}.
#' To this end, Boruta decision for factor attributes may be different.
#'
#' Random Forest methods has two main parameters, number of attributes tried at each split and the number of trees in the forest; first one is called \code{mtry} in both implementations, but the second \code{ntree} in \code{\link[randomForest]{randomForest}} and \code{num.trees} in \code{\link{ranger}}.
#' To this end, to maintain compatibility, \code{getImpRf*} functions still accept \code{ntree} parameter relaying it into \code{num.trees}.
#' Still, both parameters take the same defaults in both implementations (square root of the number all all attributes and 500 respectively).
#'
#' Moreover, \code{\link{ranger}} brings some addition capabilities to Boruta, like analysis of survival problems or sticky variables which are always considered on splits.
#'
#' Finally, the results for the same PRNG seed will be different.
#' @param x data frame of predictors including shadows.
#' @param y response vector.
#' @param ... parameters passed to the underlying \code{\link[randomForest]{randomForest}} call; they are relayed from \code{...} of \code{\link{Boruta}}.
#' @examples
#' set.seed(777)
#' #Add some nonsense attributes to iris dataset by shuffling original attributes
#' iris.extended<-data.frame(iris,apply(iris[,-5],2,sample))
#' names(iris.extended)[6:9]<-paste("Nonsense",1:4,sep="")
#' #Run Boruta on this data
#' Boruta(Species~.,getImp=getImpLegacyRfZ,
#'  data=iris.extended,doTrace=2)->Boruta.iris.extended
#' #Nonsense attributes should be rejected
#' print(Boruta.iris.extended)
#' @export
getImpLegacyRfZ<-function(x,y,...){
  randomForest::randomForest(x,y,
                             importance=TRUE,keep.forest=FALSE,...)->rf
  randomForest::importance(rf,1,scale=TRUE)[,1]
}
comment(getImpLegacyRfZ)<-'randomForest normalized permutation importance'

#' @rdname getImpLegacyRf
#' @export
getImpLegacyRfRaw<-function(x,y,...){
  randomForest::randomForest(x,y,
                             importance=TRUE,keep.forest=FALSE,...)->rf
  randomForest::importance(rf,1,scale=FALSE)[,1]
}
comment(getImpLegacyRfRaw)<-'randomForest raw permutation importance'

#' @rdname getImpLegacyRf
#' @export
getImpLegacyRfGini<-function(x,y,...){
  randomForest::randomForest(x,y,
                             keep.forest=FALSE,...)->rf
  randomForest::importance(rf,2,scale=FALSE)[,1]
}
comment(getImpLegacyRfGini)<-'randomForest Gini index importance'


#' ranger Random Forest importance adapters
#'
#' Those function is intended to be given to a \code{getImp} argument of \code{\link{Boruta}} function to be called by the Boruta algorithm as an importance source.
#' \code{getImpRfZ} generates default, normalized permutation importance, \code{getImpRfRaw} raw permutation importance, finally \code{getImpRfGini} generates Gini index importance.
#' @name getImpRf
#' @rdname getImpRf
#' @aliases getImpRfZ getImpRfGini getImpRfRaw
#' @param x data frame of predictors including shadows.
#' @param y response vector.
#' @param ntree  Number of trees in the forest; copied into \code{\link{ranger}}'s native num.trees, put to retain transparent compatibility with randomForest.
#' @param num.trees  Number of trees in the forest, as according to \code{\link{ranger}}'s nomenclature. If not given, set to \code{ntree} value. If both are given, \code{num.trees} takes precedence.
#' @param ... parameters passed to the underlying \code{\link{ranger}} call; they are relayed from \code{...} of \code{\link{Boruta}}.
#' @note Prior to Boruta 5.0, \code{getImpLegacyRfZ} function was a default importance source in Boruta; see \link{getImpLegacyRf} for more details.
#' @export
getImpRfZ<-function(x,y,ntree=500,num.trees=ntree,...){
  if(inherits(y,"Surv")){
    x$shadow.Boruta.time<-y[,"time"]
    x$shadow.Boruta.status<-y[,"status"]
    return(ranger::ranger(data=x,
                          dependent.variable.name="shadow.Boruta.time",
                          status.variable.name="shadow.Boruta.status",
                          num.trees=num.trees,importance="permutation",
                          scale.permutation.importance=TRUE,
                          write.forest=FALSE,...)$variable.importance)
  }
  #Abusing the fact that Boruta disallows attributes with names
  # starting from "shadow"
  x$shadow.Boruta.decision<-y
  ranger::ranger(data=x,dependent.variable.name="shadow.Boruta.decision",
                 num.trees=num.trees,importance="permutation",
                 scale.permutation.importance=TRUE,
                 write.forest=FALSE,...)$variable.importance
}
comment(getImpRfZ)<-'ranger normalized permutation importance'

#' @rdname getImpRf
#' @export
getImpRfGini<-function(x,y,ntree=500,num.trees=ntree,...){
  if(inherits(y,"Surv"))
    stop("Ranger cannot produce Gini importance for survival problems.")
  x$shadow.Boruta.decision<-y
  ranger::ranger(data=x,dependent.variable.name="shadow.Boruta.decision",
                 num.trees=num.trees,importance="impurity",
                 scale.permutation.importance=FALSE,
                 write.forest=FALSE,...)$variable.importance
}
comment(getImpRfGini)<-'ranger Gini index importance'

#' @rdname getImpRf
#' @export
getImpRfRaw<-function(x,y,ntree=500,num.trees=ntree,...){
  if(inherits(y,"Surv")){
    x$shadow.Boruta.time<-y[,"time"]
    x$shadow.Boruta.status<-y[,"status"]
    return(ranger::ranger(data=x,
                          dependent.variable.name="shadow.Boruta.time",
                          status.variable.name="shadow.Boruta.status",
                          num.trees=num.trees,importance="permutation",
                          write.forest=FALSE,...)$variable.importance)
  }
  x$shadow.Boruta.decision<-y
  ranger::ranger(data=x,dependent.variable.name="shadow.Boruta.decision",
                 num.trees=num.trees,importance="permutation",
                 scale.permutation.importance=FALSE,
                 write.forest=FALSE,...)$variable.importance
}
comment(getImpRfRaw)<-'ranger raw permutation importance'

#' ranger Extra-trees importance adapters
#'
#' Those function is intended to be given to a \code{getImp} argument of \code{\link{Boruta}} function to be called by the Boruta algorithm as an importance source.
#' \code{getImpExtraZ} generates default, normalized permutation importance, \code{getImpExtraRaw} raw permutation importance, finally \code{getImpExtraGini} generates Gini impurity importance.
#' @name getImpExtra
#' @rdname getImpExtra
#' @aliases getImpExtraZ getImpExtraGini getImpExtraRaw
#' @param x data frame of predictors including shadows.
#' @param y response vector.
#' @param ntree  Number of trees in the forest; copied into \code{\link{ranger}}'s native num.trees, put to retain transparent compatibility with randomForest.
#' @param num.trees  Number of trees in the forest, as according to \code{\link{ranger}}'s nomenclature. If not given, set to \code{ntree} value. If both are given, \code{num.trees} takes precedence.
#' @param ... parameters passed to the underlying \code{\link{ranger}} call; they are relayed from \code{...} of \code{\link{Boruta}}. Note that these function work just by setting \code{splitrule} to \code{"extratrees"}.
#' @export
getImpExtraZ<-function(x,y,ntree=500,num.trees=ntree,...)
  getImpRfZ(x,y,ntree=ntree,splitrule="extratrees",...)
comment(getImpExtraZ)<-'ranger normalized permutation importance'

#' @rdname getImpExtra
#' @export
getImpExtraGini<-function(x,y,ntree=500,num.trees=ntree,...)
  getImpRfGini(x,y,ntree=ntree,splitrule="extratrees",...)
comment(getImpExtraGini)<-'ranger extra-trees Gini index importance'

#' @rdname getImpExtra
#' @export
getImpExtraRaw<-function(x,y,ntree=500,num.trees=ntree,...)
  getImpRfRaw(x,y,ntree=ntree,splitrule="extratrees",...)
comment(getImpExtraRaw)<-'ranger extra-trees raw permutation importance'


#' Random Ferns importance
#'
#' This function is intended to be given to a \code{getImp} argument of \code{\link{Boruta}} function to be called by the Boruta algorithm as an importance source.
#' @param x data frame of predictors including shadows.
#' @param y response vector.
#' @param ... parameters passed to the underlying \code{\link[rFerns]{rFerns}} call; they are relayed from \code{...} of \code{\link{Boruta}}.
#' @export
#' @note Random Ferns importance calculation should be much faster than using Random Forest; however, one must first optimize the value of the \code{depth} parameter and
#' it is quite likely that the number of ferns in the ensemble required for the importance to converge will be higher than the number of trees in case of Random Forest.
getImpFerns<-function(x,y,...){
  f<-rFerns::rFerns(x,y,
                    saveForest=FALSE,importance=TRUE,...)
  f$importance[,1]
}
comment(getImpFerns)<-'rFerns importance'

#' Xgboost importance
#'
#' This function is intended to be given to a \code{getImp} argument of \code{\link{Boruta}} function to be called by the Boruta algorithm as an importance source.
#' This functionality is inspired by the Python package BoostARoota by Chase DeHan.
#' In practice, due to the eager way XgBoost works, this adapter changes Boruta into minimal optimal method, hence I strongly recommend against using this.
#' @param x data frame of predictors including shadows.
#' @param y response vector.
#' @param nrounds Number of rounds; passed to the underlying \code{\link[xgboost]{xgboost}} call.
#' @param verbose Verbosity level of xgboost; either 0 (silent) or 1 (progress reports). Passed to the underlying \code{\link[xgboost]{xgboost}} call.
#' @param ... other parameters passed to the underlying \code{\link[xgboost]{xgboost}} call.
#' Similarly as \code{nrounds} and \code{verbose}, they are relayed from \code{...} of \code{\link{Boruta}}.
#' For convenience, this function sets \code{nrounds} to 5 and verbose to 0, but this can be overridden.
#' @note Only dense matrix interface is supported; all predictions given to \code{\link{Boruta}} call have to be numeric (not integer).
#' Categorical features should be split into indicator attributes.
#' @references \url{https://github.com/chasedehan/BoostARoota}
#' @export
getImpXgboost<-function(x,y,nrounds=5,verbose=0,...){
  for(e in 1:ncol(x)) x[,e]<-as.numeric(x[,e])
  xgboost::xgb.importance(
    model=xgboost::xgboost(
      data=as.matrix(x),
      label=y,
      nrounds=nrounds,
      verbose=verbose,
      ...
    )
  )->imp
  stats::setNames(rep(0,ncol(x)),colnames(x))->ans
  ans[imp$Feature]<-imp$Gain
  ans
}
comment(getImpXgboost)<-'xgboost gain importance'