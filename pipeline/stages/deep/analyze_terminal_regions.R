library(limma);library(edgeR)
B=Sys.getenv("OSDR_DEEP_DIR");md=read.csv(file.path(B,'tables/terminal_region_mouse_registry.csv'),stringsAsFactors=FALSE)
ann=read.csv(file.path(B,'prepared/additional_gene_annotation.csv'),stringsAsFactors=FALSE)
map=read.csv(file.path(Sys.getenv("OSDR_CORE_DIR"), "tables/verified_panel_mapping.csv"),stringsAsFactors=FALSE)
rows=list();qcs=list()
for(sid in c(612,613)){
 s=md[md$OSD==paste0('OSD-',sid),,drop=FALSE];X=as.matrix(read.csv(gzfile(file.path(B,'prepared',paste0('OSD-',sid,'_cerebral_mouse_counts.csv.gz'))),row.names=1,check.names=FALSE));X=X[,s$key,drop=FALSE]
 if(sid==612){s$group=factor(s$flight,levels=c('GC','FLT'));cs=matrix(c(-1,1),ncol=1);tests='flight_vs_ground'}else{
  s$group=factor(paste(s$flight,s$age,sep='_'),levels=c('GC_12','GC_29','FLT_12','FLT_29'))
  cs=cbind(c(-1,0,1,0),c(0,-1,0,1),c(-.5,-.5,.5,.5),c(1,-1,-1,1));tests=c('flight_age12','flight_age29','flight_balanced_age','flight_by_age')
 }
 D=model.matrix(~0+group,s);stopifnot(qr(D)$rank==ncol(D));d=DGEList(X);keep=filterByExpr(d,design=D);d=calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method='TMM')
 v=voom(d,D,plot=FALSE);fit=eBayes(contrasts.fit(lmFit(v,D),cs),robust=TRUE);d=estimateDisp(d,D,robust=TRUE);ql=glmQLFit(d,D,robust=TRUE)
 write.csv(v$E,gzfile(file.path(B,'prepared',paste0('OSD-',sid,'_voom_logCPM.csv.gz'))))
 qcs[[length(qcs)+1]]=data.frame(OSD=paste0('OSD-',sid),n_mice=nrow(s),n_filtered_genes=sum(keep),design_rank=qr(D)$rank,residual_df=nrow(D)-qr(D)$rank)
 for(k in seq_along(tests))for(method in c('limma_voom','edgeR_QL')){
  if(method=='limma_voom'){
   res=topTable(fit,coef=k,number=Inf,sort.by='none');se=fit$stdev.unscaled[,k]*sqrt(fit$s2.post);res$ci_low=res$logFC-qt(.975,fit$df.total)*se;res$ci_high=res$logFC+qt(.975,fit$df.total)*se
  }else{res=topTags(glmQLFTest(ql,contrast=cs[,k]),n=Inf,sort.by='none')$table;res$P.Value=res$PValue;res$ci_low=NA;res$ci_high=NA}
  res$SYMBOL=ann$SYMBOL[match(rownames(res),ann$ENSEMBL)];write.csv(res,gzfile(file.path(B,'tables',paste0('OSD-',sid,'_',method,'_',tests[k],'_genome.csv.gz'))))
  for(i in seq_len(nrow(map))){hits=which(!is.na(res$SYMBOL)&res$SYMBOL==map$mouse_symbol_primary[i]);ok=length(hits)==1
   rows[[length(rows)+1]]=data.frame(OSD=paste0('OSD-',sid),method=method,contrast=tests[k],analyte=map$analyte[i],mouse_symbol=map$mouse_symbol_primary[i],log2FC=if(ok)res$logFC[hits]else NA,ci_low=if(ok)res$ci_low[hits]else NA,ci_high=if(ok)res$ci_high[hits]else NA,pvalue=if(ok)res$P.Value[hits]else 1,n_mice=nrow(s),status=if(ok)'tested'else 'unavailable_or_nonunique')
  }
 }
 cat(sid,nrow(s),'mice',sum(keep),'genes\n')
}
a=do.call(rbind,rows);a$q_full_755_family=ave(a$pvalue,a$method,FUN=function(p)p.adjust(p,'BH'))
write.csv(a,file.path(B,'tables/terminal_region_panel_complete.csv'),row.names=FALSE)
write.csv(do.call(rbind,qcs),file.path(B,'tables/terminal_region_model_QC.csv'),row.names=FALSE)
