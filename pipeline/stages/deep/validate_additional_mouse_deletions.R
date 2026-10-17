# Every new-cohort mouse deletion, whole-genome models; unchanged companion
# cohort remains in the same 1057-hypothesis correction family per engine.
library(limma);library(edgeR)
B=Sys.getenv("OSDR_DEEP_DIR")
map=read.csv(file.path(Sys.getenv("OSDR_CORE_DIR"), "tables/verified_panel_mapping.csv"),stringsAsFactors=FALSE)
ann=read.csv(file.path(B,'prepared/additional_gene_annotation.csv'),stringsAsFactors=FALSE)
md=read.csv(file.path(B,'tables/additional_mouse_registry.csv'),stringsAsFactors=FALSE)
baseline=read.csv(file.path(B,'tables/additional_panel_complete.csv'),stringsAsFactors=FALSE)
focus=c('AIM2','AIF1','PYCARD','SLC1A3','S100B')
outputs=list();complete=list()
for(sid in c(525,457)){
  s0=md[md$OSD==paste0('OSD-',sid),,drop=FALSE]
  X0=as.matrix(read.csv(gzfile(file.path(B,'prepared',paste0('OSD-',sid,'_mouse_counts.csv.gz'))),row.names=1,check.names=FALSE))
  for(drop in s0$donor){
    s=s0[s0$donor!=drop,,drop=FALSE];X=X0[,s$key,drop=FALSE]
    if(sid==525){
      s$group=factor(s$flight,levels=c('GC','FLT','VIV'));cs=cbind(c(-1,1,0),c(0,1,-1),c(1,0,-1));tests=c('flight_vs_habitat_ground','flight_vs_vivarium','habitat_ground_vs_vivarium')
    }else{
      s$group=factor(paste(s$flight,s$genotype,sep='_'),levels=c('GC_WT','FLT_WT','GC_Nrf2KO','FLT_Nrf2KO'));cs=cbind(c(-1,1,0,0),c(0,0,-1,1),c(1,-1,-1,1),c(-.5,.5,-.5,.5));tests=c('flight_WT','flight_Nrf2KO','flight_by_Nrf2KO','flight_balanced_genotype')
    }
    D=model.matrix(~0+group,s);d=DGEList(X);keep=filterByExpr(d,design=D);d=calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method='TMM')
    v=voom(d,D,plot=FALSE);fit=eBayes(contrasts.fit(lmFit(v,D),cs),robust=TRUE);d=estimateDisp(d,D,robust=TRUE);ql=glmQLFit(d,D,robust=TRUE)
    for(method in c('limma_voom','edgeR_QL')){
      a=baseline[baseline$method==method,,drop=FALSE]
      for(k in seq_along(tests)){
        res=if(method=='limma_voom')topTable(fit,coef=k,number=Inf,sort.by='none') else topTags(glmQLFTest(ql,contrast=cs[,k]),n=Inf,sort.by='none')$table
        if(method=='edgeR_QL')res$P.Value=res$PValue
        res$SYMBOL=ann$SYMBOL[match(rownames(res),ann$ENSEMBL)]
        for(i in seq_len(nrow(map))){
          hits=which(!is.na(res$SYMBOL)&res$SYMBOL==map$mouse_symbol_primary[i]);ok=length(hits)==1
          idx=which(a$OSD==paste0('OSD-',sid)&a$contrast==tests[k]&a$analyte==map$analyte[i]);stopifnot(length(idx)==1)
          a$pvalue[idx]=if(ok)res$P.Value[hits]else 1;a$log2FC[idx]=if(ok)res$logFC[hits]else NA
          a$status[idx]=if(ok)'tested'else 'unavailable_or_nonunique'
        }
      }
      a$q_full_1057_family=p.adjust(a$pvalue,'BH');a$deleted_OSD=paste0('OSD-',sid);a$deleted_mouse=drop
      complete[[length(complete)+1]]=a[,c('deleted_OSD','deleted_mouse','OSD','method','contrast','analyte','log2FC','pvalue','q_full_1057_family','status')]
      outputs[[length(outputs)+1]]=a[a$OSD==paste0('OSD-',sid)&a$analyte%in%focus,c('deleted_OSD','deleted_mouse','OSD','method','contrast','analyte','log2FC','pvalue','q_full_1057_family','status')]
    }
    cat('Deleted',sid,drop,'\n')
  }
}
write.csv(do.call(rbind,complete),gzfile(file.path(B,'tables/additional_mouse_deletion_complete.csv.gz')),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(B,'tables/additional_mouse_deletion_focus.csv'),row.names=FALSE)
