# Inspect the installed CAMERA implementation and actual designs; no model refit.
library(limma)
B=Sys.getenv("OSDR_EXTENDED_DIR")
writeLines(c(paste('limma version',packageVersion('limma')),deparse(camera.default)),file.path(B,'logs/CAMERA_implementation_audit.txt'))
rows=list()
for(p in list.files(file.path(B,'prepared'),pattern='_voom_fit.rds$',full.names=TRUE)){
    x=readRDS(p);unit=sub('_voom_fit.rds','',basename(p));n=nrow(x$design);rank=qr(x$design)$rank;G=nrow(x$voom$E)
    stopifnot(rank==ncol(x$design))
    rows[[length(rows)+1]]=data.frame(unit=unit,n_samples=n,design_rank=rank,n_genes=G,fixed_correlation_df=G-2,estimated_correlation_df=min(n-rank,G-2))
}
md=read.csv(file.path(B,'tables/GSE295428_sorted_cell_QC.csv'))
counts=as.matrix(read.csv(gzfile(file.path(B,'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz')),row.names=1,check.names=FALSE))
for(tissue in sort(unique(md$tissue))){
    s=md[md$tissue==tissue & md$n_retained_cells>=50,,drop=FALSE];s$group=factor(paste(s$flight,s$age,sep='_'),levels=c('GC_3','GC_8','FLT_3','FLT_8'))
    X=counts[,s$sample,drop=FALSE];G=sum(rowSums(X>=10)>=3);D=model.matrix(~0+group,s);n=nrow(D);rank=qr(D)$rank
    stopifnot(rank==ncol(D))
    rows[[length(rows)+1]]=data.frame(unit=tissue,n_samples=n,design_rank=rank,n_genes=G,fixed_correlation_df=G-2,estimated_correlation_df=min(n-rank,G-2))
}
write.csv(do.call(rbind,rows),file.path(B,'tables/pathway_camera_degrees_of_freedom.csv'),row.names=FALSE)
cat('Audited',length(rows),'CAMERA study/population designs\n')
