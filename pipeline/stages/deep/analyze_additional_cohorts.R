library(limma);library(edgeR)
B=Sys.getenv("OSDR_DEEP_DIR")
map=read.csv(file.path(Sys.getenv("OSDR_CORE_DIR"), "tables/verified_panel_mapping.csv"),stringsAsFactors=FALSE)
ann=read.csv(file.path(B,'prepared/additional_gene_annotation.csv'),stringsAsFactors=FALSE)
md=read.csv(file.path(B,'tables/additional_mouse_registry.csv'),stringsAsFactors=FALSE)
rows=list();qcs=list()
for(sid in c(525,457)){
  s=md[md$OSD==paste0('OSD-',sid),,drop=FALSE]
  X=as.matrix(read.csv(gzfile(file.path(B,'prepared',paste0('OSD-',sid,'_mouse_counts.csv.gz'))),row.names=1,check.names=FALSE));X=X[,s$key,drop=FALSE]
  if(sid==525){
    s$group=factor(s$flight,levels=c('GC','FLT','VIV'))
    cs=cbind(c(-1,1,0),c(0,1,-1),c(1,0,-1));tests=c('flight_vs_habitat_ground','flight_vs_vivarium','habitat_ground_vs_vivarium')
  } else {
    s$group=factor(paste(s$flight,s$genotype,sep='_'),levels=c('GC_WT','FLT_WT','GC_Nrf2KO','FLT_Nrf2KO'))
    cs=cbind(c(-1,1,0,0),c(0,0,-1,1),c(1,-1,-1,1),c(-.5,.5,-.5,.5));tests=c('flight_WT','flight_Nrf2KO','flight_by_Nrf2KO','flight_balanced_genotype')
  }
  colnames(cs)=tests;D=model.matrix(~0+group,s);colnames(D)=levels(s$group)
  # Standard whole-genome group-aware abundance gate, applied before panel extraction.
  d=DGEList(X);keep=filterByExpr(d,design=D);d=calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method='TMM')
  v=voom(d,D,plot=FALSE);fit=eBayes(contrasts.fit(lmFit(v,D),cs),robust=TRUE)
  d=estimateDisp(d,D,robust=TRUE);ql=glmQLFit(d,D,robust=TRUE)
  write.csv(v$E,gzfile(file.path(B,'prepared',paste0('OSD-',sid,'_voom_logCPM.csv.gz'))))
  qcs[[length(qcs)+1]]=data.frame(OSD=paste0('OSD-',sid),n_mice=nrow(s),filtered_genes=sum(keep),design_rank=qr(D)$rank,residual_df=nrow(D)-qr(D)$rank,min_TMM=min(d$samples$norm.factors),max_TMM=max(d$samples$norm.factors))
  for(k in seq_along(tests))for(method in c('limma_voom','edgeR_QL')){
    if(method=='limma_voom'){
      res=topTable(fit,coef=k,number=Inf,sort.by='none');se=fit$stdev.unscaled[,k]*sqrt(fit$s2.post)
      res$ci_low=res$logFC-qt(.975,fit$df.total)*se;res$ci_high=res$logFC+qt(.975,fit$df.total)*se
    }else{
      res=topTags(glmQLFTest(ql,contrast=cs[,k]),n=Inf,sort.by='none')$table;res$P.Value=res$PValue;res$ci_low=NA;res$ci_high=NA
    }
    res$SYMBOL=ann$SYMBOL[match(rownames(res),ann$ENSEMBL)]
    write.csv(res,gzfile(file.path(B,'tables',paste0('OSD-',sid,'_',method,'_',tests[k],'_genome.csv.gz'))))
    for(i in seq_len(nrow(map))){
      hits=which(!is.na(res$SYMBOL)&res$SYMBOL==map$mouse_symbol_primary[i]);ok=length(hits)==1
      rows[[length(rows)+1]]=data.frame(OSD=paste0('OSD-',sid),method=method,contrast=tests[k],analyte=map$analyte[i],mouse_symbol=map$mouse_symbol_primary[i],
       log2FC=if(ok)res$logFC[hits]else NA,ci_low=if(ok)res$ci_low[hits]else NA,ci_high=if(ok)res$ci_high[hits]else NA,pvalue=if(ok)res$P.Value[hits]else 1,
       n_mice=nrow(s),status=if(ok)'tested'else 'unavailable_or_nonunique')
    }
  }
  cat('OSD-',sid,' ',nrow(s),'mice, ',sum(keep),'filtered genes\n')
}
a=do.call(rbind,rows);a$pvalue[!is.finite(a$pvalue)]=1
a$q_full_1057_family=ave(a$pvalue,a$method,FUN=function(p)p.adjust(p,'BH'))
a$q_within_study_family=ave(a$pvalue,interaction(a$OSD,a$method),FUN=function(p)p.adjust(p,'BH'))
write.csv(a,file.path(B,'tables/additional_panel_complete.csv'),row.names=FALSE)
write.csv(do.call(rbind,qcs),file.path(B,'tables/additional_model_QC.csv'),row.names=FALSE)
print(with(a,table(OSD,method,q_full_1057_family<.05)))
