# Independent implementations on the same samples: sensitivity, not replication.
library(limma);library(edgeR)
B=Sys.getenv("OSDR_EXTENDED_DIR");V1=Sys.getenv("OSDR_CORE_DIR")
map=read.csv(file.path(V1,'tables/verified_panel_mapping.csv'),stringsAsFactors=FALSE)
save_panel=function(res,ann,study,tag,method){
  idx=match(rownames(res),ann$ENSEMBL);res$symbol=ann$SYMBOL[idx];res$gene_id=rownames(res)
  rows=lapply(seq_len(nrow(map)),function(i){
    hits=which(!is.na(res$symbol)&res$symbol==map$mouse_symbol_primary[i]);ok=length(hits)==1
    data.frame(study=study,contrast=tag,method=method,analyte=map$analyte[i],mouse_symbol=map$mouse_symbol_primary[i],
               log2FC=if(ok)res$logFC[hits] else NA,pvalue=if(ok)res$P.Value[hits] else 1,
               gene_id=if(ok)res$gene_id[hits] else '',test_status=if(ok)'tested' else 'unavailable_or_nonunique')
  });do.call(rbind,rows)
}
run=function(counts,md,ann,study,core=FALSE){
  stopifnot(identical(colnames(counts),rownames(md)))
  keep=rowSums(counts>=10)>=3;counts=counts[keep,,drop=FALSE]
  md$group=factor(md$group);design=model.matrix(~0+group,md);colnames(design)=levels(md$group)
  stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
  vec=function(a,b)as.numeric(colnames(design)==a)-as.numeric(colnames(design)==b)
  contrasts=list()
  if(core){
    for(cohort in c('ISS_T','LAR'))for(age in c('12','29'))contrasts[[paste0('flight_',cohort,'_age',age)]]=vec(paste('FLT',cohort,age,sep='_'),paste('GC',cohort,age,sep='_'))
    for(f in c('GC','FLT')){
      for(cohort in c('ISS_T','LAR'))contrasts[[paste0('age_',f,'_',cohort)]]=vec(paste(f,cohort,'29',sep='_'),paste(f,cohort,'12',sep='_'))
      for(age in c('12','29'))contrasts[[paste0('cohort_',f,'_age',age)]]=vec(paste(f,'ISS_T',age,sep='_'),paste(f,'LAR',age,sep='_'))
    }
    for(age in c('12','29'))contrasts[[paste0('flight_by_cohort_age',age)]]=vec(paste('FLT','ISS_T',age,sep='_'),paste('GC','ISS_T',age,sep='_'))-vec(paste('FLT','LAR',age,sep='_'),paste('GC','LAR',age,sep='_'))
  }else contrasts[['flight_vs_ground']]=vec('FLT','GC')
  d=calcNormFactors(DGEList(counts=counts),method='TMM')
  v=voom(d,design,plot=FALSE);fit=eBayes(contrasts.fit(lmFit(v,design),do.call(cbind,contrasts)),robust=TRUE)
  d=estimateDisp(d,design,robust=TRUE);ql=glmQLFit(d,design,robust=TRUE)
  p=list();p2=list()
  for(k in seq_along(contrasts)){
    tag=names(contrasts)[k]
    res=topTable(fit,coef=k,number=Inf,sort.by='none',adjust.method='BH')
    write.csv(res,gzfile(file.path(B,'tables',paste0(study,'_voom_',tag,'_genome.csv.gz'))))
    p[[k]]=save_panel(res,ann,study,tag,'limma_voom_TMM')
    er=topTags(glmQLFTest(ql,contrast=contrasts[[k]]),n=Inf,sort.by='none')$table;er$P.Value=er$PValue
    write.csv(er,gzfile(file.path(B,'tables',paste0(study,'_edgeR_',tag,'_genome.csv.gz'))))
    p2[[k]]=save_panel(er,ann,study,tag,'edgeR_QL_TMM')
  }
  write.csv(data.frame(sample=colnames(counts),TMM_norm_factor=d$samples$norm.factors,library=d$samples$lib.size),file.path(B,'tables',paste0(study,'_TMM_factors.csv')),row.names=FALSE)
  write.csv(v$E,gzfile(file.path(B,'prepared',paste0(study,'_voom_logCPM.csv.gz'))))
  saveRDS(list(voom=v,fit=fit,design=design,metadata=md,annotation=ann,contrasts=contrasts),file.path(B,'prepared',paste0(study,'_voom_fit.rds')))
  cat(study,nrow(counts),ncol(counts),'genes/samples\n')
  rbind(do.call(rbind,p),do.call(rbind,p2))
}
core=list();external=list();cells=list()
for(sid in c(561,562,352,563,564)){
  old=sid %in% c(561,562);src=if(old)V1 else B
  md=read.csv(file.path(src,'tables',paste0('OSD-',sid,if(old)'_sample_metadata.csv' else '_metadata.csv')),row.names=1,check.names=FALSE)
  count=as.matrix(read.csv(gzfile(file.path(src,'prepared',paste0('OSD-',sid,'_STAR_endogenous.csv.gz'))),row.names=1,check.names=FALSE))
  count=count[,rownames(md),drop=FALSE]
  ann=read.csv(file.path(src,'prepared',paste0('OSD-',sid,'_annotation.csv')),stringsAsFactors=FALSE)
  res=run(count,md,ann,paste0('OSD-',sid),old)
  if(old)core[[length(core)+1]]=res else external[[length(external)+1]]=res
}
for(level in c('l2','l3')){
  count=as.matrix(read.csv(gzfile(file.path(B,'prepared',if(level=='l2')'OSD352_RNA_pseudobulk.csv.gz' else 'OSD352_RNA_subtype_pseudobulk.csv.gz')),row.names=1,check.names=FALSE))
  md=read.csv(file.path(B,'tables',paste0('OSD352_',level,'_registry.csv')),stringsAsFactors=FALSE)
  selected=if(level=='l2')c('CNS neurons','CNS glia','Immune cells') else c('Astroependymal cells','Oligodendrocytes')
  for(ct in selected){
    s=md[md$cell_type==ct,,drop=FALSE];stopifnot(nrow(s)==5,all(s$n_nuclei>=20));rownames(s)=s$column;s$group=ifelse(s$flight_status=='flight','FLT','GC')
    ann=read.csv(file.path(B,'prepared/OSD352_RNA_annotation.csv'),stringsAsFactors=FALSE);study=paste0('OSD352_',level,'_',gsub('[^A-Za-z0-9]+','_',ct))
    cells[[length(cells)+1]]=run(count[,s$column,drop=FALSE],s,ann,study)
  }
}
finish=function(items,name){
  a=do.call(rbind,items);a$pvalue[!is.finite(a$pvalue)]=1
  a$q_full_family=ave(a$pvalue,a$method,FUN=function(p)p.adjust(p,'BH'))
  write.csv(a,file.path(B,'tables',name),row.names=FALSE)
  print(with(a,table(method,q_full_family<.05)))
}
finish(core,'bulk_method_comparison_panel.csv');finish(external,'independent_bulk_panel.csv');finish(cells,'celltype_pseudobulk_panel.csv')
writeLines(capture.output(sessionInfo()),file.path(B,'logs/R_sessionInfo.txt'))
