NeuralGasClustering <-function(Data,ClusterNo,PlotIt=FALSE,...){
   
  if (!requireNamespace('cclust',quietly = TRUE)) {
    message(
      'Subordinate clustering package (cclust) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      list(
        Cls = rep(1, nrow(Data)),
        Object = "Subordinate clustering package (cclust) is missing.
                Please install the package which is defined in 'Suggests'."
      )
    )
  }
  res=cclust::cclust(x=Data,centers=ClusterNo,method='neuralgas',...)
  Cls=res$cluster
  if(PlotIt){
    ClusterPlotMDS(Data,Cls)
  }
  Cls=ClusterRename(Cls,Data)
  return(list(Cls=Cls,Object=res))
  }