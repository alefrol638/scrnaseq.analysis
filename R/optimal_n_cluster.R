nbclust<-setClass("nbclust", slots=list(full_list = "list",
                                        top3 = "character", plots = "list",new.meta="data.frame"))
#' @exportClass nbclust

#' @title determine the optimal number of clusters using Nbclust
#' find optimal #cluster using NbClust
#' @param x Seurat Object
#' @param method character, Clustering method, closest coming to Seurat clustering is kmeans
#' @param index  Character vector, select indeces to determine best number of clusters one of: "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw
#'  if nothing is specified, all of them are used
#' @param min.nc,max.nc integer, minimal and maximal number of clusters to calculate for
#' @param only.metaif will skip the NbClust step and only create meta data with the nbclust clustering information from an existing nbclust class
#' @param name.graph which neighbor graph should be used? (<name of reduction_snn OR nn>)
#' @param dims number of dimensions to use for clustering
#' @export optimal_n_cluster
optimal_n_cluster<-function(x,out=new("nbclust"),method="kmeans",index="all",min.nc=2,max.nc=15,only.meta=F,name.graph="SCT_snn",dims=1:30)
{
  if(only.meta==F)
  {
    if(index=="all")
    {
      index_nbclust<-c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew",
                       "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2",
                       "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma",
                       "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")

    }else{index_nbclust=index}

    out@full_list<-foreach::foreach(i = index_nbclust,.packages = "NbClust",.errorhandling = "pass") %dopar%
      {
        R.utils::withTimeout({out@full_list[[i]]<-NbClust::NbClust(x@reductions$pca@cell.embeddings[,dims],diss=x@graphs[[name.graph]],distance=NULL,
                                                    method=method,index=i,min.nc = min.nc,max.nc = max.nc)},
                    timeout = 180,cpu=F,onTimeout = "warning")
      }
    names(out@full_list)<-index_nbclust
    #graphical approach ... find elbow in left dindex plot and peak in right
    indices<-list()

    for(i in c("hubert","dindex")){
      if(!is.character(out@full_list[[i]])){
        indices[[i]]<-data.frame("All.index"=out@full_list[[i]]$All.index,"cluster"=names(out@full_list[[i]]$All.index))
        for(k in 1:length(indices[[i]]$All.index)){
          indices[[i]]$diff[k]<-indices[[i]]$All.index[k]-indices[[i]]$All.index[k+1]
        }
        for(k in 1:length(indices[[i]]$All.index)){
          indices[[i]]$diff2nd[k]<-indices[[i]]$diff[k]-indices[[i]]$diff[k+1]
        }
        indices[[i]]$cluster<-as.integer(indices[[i]]$cluster)

        out@plots[[i]]<-cowplot::plot_grid(ggplot(indices[[i]],aes(y=All.index,x=cluster,group=1))+geom_line()+geom_point()+theme_classic(),ggplot(indices[[i]],aes(y=diff2nd,x=cluster,group=1))+geom_line()+geom_point()+theme_classic())
      }
    }
    #extract optimal #clust from each approaches
    nc<-c()
    for(i in names(out@full_list)[names(out@full_list)!="hubert"&names(out@full_list)!="dindex"]){
      if(!is.character(out@full_list[[i]]))
      {
          try({nc[i]<-out@full_list[[i]][["Best.nc"]][["Number_clusters"]]})

       }

    }

    #14 and 17 cluster seems to be dominant
    nc_freq<-table(nc)
    nc_freq<-sort(nc_freq,decreasing = T)
    out@top3<-names(nc_freq[1:3])
    nc_top3<-nc[nc%in%out@top3&!is.na(nc)]

  }


  # extract clustering to compare to seurat clusters
  k=1
  for(i in names(nc)){
    attach(out@full_list[[i]])
    if(exists("Best.nc")){
      if( out@full_list[[i]][["Best.nc"]][["Number_clusters"]]%in%out@top3)
      {
        clusters<-grep(paste("^",name.graph,"_res\\..*$",sep =""),names(x[[]]),fixed=F,value = T)
        cluster_n<-c()
        for(j in 1:length(clusters)){
          cluster_n[j]<- length(levels(x@meta.data[,clusters[j]]))
        }
        names(cluster_n)<-clusters
        res<-names(cluster_n)[cluster_n==out@full_list[[i]][["Best.nc"]][["Number_clusters"]]]
        if(length(res)!=0){
          try({x[[paste(res[1],k,sep="")]]<-out@full_list[[i]][["Best.partition"]]})
        }
      }
    }
    detach()
    k<-k+1
  }
  out@new.meta<-x@meta.data

  return(out)
}
