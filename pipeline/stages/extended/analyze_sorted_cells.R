# External sorted-cell pseudobulk: age-adjusted and interaction contrasts.
library(limma);library(edgeR)
B=Sys.getenv("OSDR_EXTENDED_DIR");map=read.csv(file.path(Sys.getenv("OSDR_CORE_DIR"), "tables/verified_panel_mapping.csv"),stringsAsFactors=FALSE)
md=read.csv(file.path(B,'tables/GSE295428_sorted_cell_QC.csv'),stringsAsFactors=FALSE);ann=read.csv(file.path(B,'prepared/GSE295428_sorted_annotation.csv'),stringsAsFactors=FALSE)
rows=list();summary=list()
for(qc in c('published_gate','alternative_gate','marker_positive_gate')){
  filename=switch(qc,published_gate='GSE295428_sorted_donor_pseudobulk.csv.gz',alternative_gate='GSE295428_sorted_sensitivity_pseudobulk.csv.gz',marker_positive_gate='GSE295428_sorted_marker_positive_pseudobulk.csv.gz')
  counts=as.matrix(read.csv(gzfile(file.path(B,'prepared',filename)),row.names=1,check.names=FALSE))
  for(tissue in sort(unique(md$tissue))){
    s=md[md$tissue==tissue,,drop=FALSE];s$n_cells=switch(qc,published_gate=s$n_retained_cells,alternative_gate=s$n_sensitivity_cells,marker_positive_gate=s$n_marker_positive_cells)
    s=s[s$n_cells>=50,,drop=FALSE];s$group=factor(paste(s$flight,s$age,sep='_'),levels=c('GC_3','GC_8','FLT_3','FLT_8'));rownames(s)=s$sample
    available=nrow(s)>=6 && all(table(s$group)>0)
    tests=c('flight_age3','flight_age8','flight_balanced_age','flight_by_age');cs=cbind(c(-1,0,1,0),c(0,-1,0,1),c(-.5,-.5,.5,.5),c(1,-1,-1,1));colnames(cs)=tests
    if(available){
      X=counts[,s$sample,drop=FALSE];keep=rowSums(X>=10)>=3;X=X[keep,,drop=FALSE];D=model.matrix(~0+group,s);colnames(D)=levels(s$group)
      d=calcNormFactors(DGEList(X),method='TMM');v=voom(d,D,plot=FALSE);fit=eBayes(contrasts.fit(lmFit(v,D),cs),robust=TRUE);d=estimateDisp(d,D,robust=TRUE);ql=glmQLFit(d,D,robust=TRUE)
      write.csv(v$E,gzfile(file.path(B,'prepared',paste0('Sorted_',gsub('[^A-Za-z0-9]','_',tissue),'_',qc,'_voom_logCPM.csv.gz'))))
      out=list()
      for(k in seq_along(tests)){
        vv=topTable(fit,coef=k,number=Inf,sort.by='none');ee=topTags(glmQLFTest(ql,contrast=cs[,k]),n=Inf,sort.by='none')$table;ee$P.Value=ee$PValue
        out[[k]]=list(limma_voom=vv,edgeR_QL=ee)
      }
      summary[[length(summary)+1]]=data.frame(tissue=tissue,QC=qc,n_donors=nrow(s),n_flight=sum(s$flight=='FLT'),n_ground=sum(s$flight=='GC'),n_filtered_genes=nrow(X),n_retained_cells=sum(s$n_cells),status='model_fitted')
    }else summary[[length(summary)+1]]=data.frame(tissue=tissue,QC=qc,n_donors=nrow(s),n_flight=sum(s$flight=='FLT'),n_ground=sum(s$flight=='GC'),n_filtered_genes=0,n_retained_cells=sum(s$n_cells),status='insufficient_donor_or_design_cell_coverage')
    for(k in seq_along(tests))for(method in c('limma_voom','edgeR_QL')){
      if(available){res=out[[k]][[method]];res$SYMBOL=ann$SYMBOL[match(rownames(res),ann$ENSEMBL)];write.csv(res,gzfile(file.path(B,'tables',paste0('Sorted_',gsub('[^A-Za-z0-9]','_',tissue),'_',qc,'_',method,'_',tests[k],'_genome.csv.gz'))))}
      for(i in seq_len(nrow(map))){
        hits=if(available)which(!is.na(res$SYMBOL)&res$SYMBOL==map$mouse_symbol_primary[i]) else integer();ok=length(hits)==1
        rows[[length(rows)+1]]=data.frame(tissue=tissue,QC=qc,method=method,contrast=tests[k],analyte=map$analyte[i],mouse_symbol=map$mouse_symbol_primary[i],log2FC=if(ok)res$logFC[hits] else NA,pvalue=if(ok)res$P.Value[hits] else 1,n_donors=nrow(s),n_flight=sum(s$flight=='FLT'),n_ground=sum(s$flight=='GC'),status=if(ok)'tested' else if(!available)'insufficient_design_cell_coverage' else 'unavailable_or_nonunique')
      }
    }
    cat(tissue,qc,nrow(s),'donors, model available:',available,'\n')
  }
}
a=do.call(rbind,rows);a$pvalue[!is.finite(a$pvalue)]=1;a$q_full_3624_family=ave(a$pvalue,interaction(a$QC,a$method),FUN=function(p)p.adjust(p,'BH'));write.csv(a,file.path(B,'tables/sorted_cell_panel_complete.csv'),row.names=FALSE);write.csv(do.call(rbind,summary),file.path(B,'tables/sorted_cell_model_summary.csv'),row.names=FALSE)
print(with(a,table(QC,method,q_full_3624_family<.05)))
