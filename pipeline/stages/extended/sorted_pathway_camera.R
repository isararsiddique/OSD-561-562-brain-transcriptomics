# Full 50-set RR8 pathway check, primary QC, age-adjusted contrasts.
library(limma);library(edgeR)
B=Sys.getenv("OSDR_EXTENDED_DIR");md=read.csv(file.path(B,'tables/GSE295428_sorted_cell_QC.csv'));ann=read.csv(file.path(B,'prepared/GSE295428_sorted_annotation.csv'))
counts=as.matrix(read.csv(gzfile(file.path(B,'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz')),row.names=1,check.names=FALSE))
lines=strsplit(readLines(file.path(B,'sources/mh.all.v2025.1.Mm.symbols.gmt')),'\t');sets=lapply(lines,function(x)x[-c(1,2)]);names(sets)=vapply(lines,`[`,character(1),1)
tests=c('flight_age3','flight_age8','flight_balanced_age','flight_by_age');cs=cbind(c(-1,0,1,0),c(0,-1,0,1),c(-.5,-.5,.5,.5),c(1,-1,-1,1));colnames(cs)=tests
rows=list()
for(tissue in sort(unique(md$tissue))){
    s=md[md$tissue==tissue & md$n_retained_cells>=50,,drop=FALSE];s$group=factor(paste(s$flight,s$age,sep='_'),levels=c('GC_3','GC_8','FLT_3','FLT_8'));X=counts[,s$sample,drop=FALSE];X=X[rowSums(X>=10)>=3,,drop=FALSE]
    D=model.matrix(~0+group,s);d=calcNormFactors(DGEList(X),method='TMM');v=voom(d,D,plot=FALSE);symbols=ann$SYMBOL[match(rownames(X),ann$ENSEMBL)];indices=lapply(sets,function(g)which(symbols %in% g))
    for(k in seq_along(tests))for(mode in c('default_0.01','conservative_0.10','estimated')){
        cor=if(mode=='estimated')NA else if(mode=='default_0.01').01 else .10
        a=camera(v,indices,design=D,contrast=cs[,k],inter.gene.cor=cor,sort=FALSE);a$set=rownames(a);a$tissue=tissue;a$contrast=tests[k];a$correlation_model=mode;a$n_donors=nrow(s)
        if(!'Correlation' %in% names(a))a$Correlation=cor
        a$PValue[a$NGenes<10|!is.finite(a$PValue)]=1;rows[[length(rows)+1]]=a
    }
    cat(tissue,'complete\n')
}
a=do.call(rbind,rows);a$q_full_1200_family=ave(a$PValue,a$correlation_model,FUN=function(p)p.adjust(p,'BH'));write.csv(a,file.path(B,'tables/sorted_pathway_CAMERA_complete.csv'),row.names=FALSE)
print(with(a,table(correlation_model,q_full_1200_family<.05)))
