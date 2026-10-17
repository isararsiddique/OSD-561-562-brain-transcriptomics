# Competitive gene-set test with intergene correlation; no isolated gene cherry-picking.
library(limma)
B=Sys.getenv("OSDR_EXTENDED_DIR")
gmt=readLines(file.path(B,'sources/mh.all.v2025.1.Mm.symbols.gmt'))
sets=lapply(strsplit(gmt,'\t'),function(x)x[-c(1,2)]);names(sets)=vapply(strsplit(gmt,'\t'),`[`,character(1),1)
all=list()
files=list.files(file.path(B,'prepared'),pattern='_voom_fit.rds$',full.names=TRUE)
for(p in files){
  x=readRDS(p);sid=sub('_voom_fit.rds','',basename(p));symbols=x$annotation$SYMBOL[match(rownames(x$voom$E),x$annotation$ENSEMBL)]
  idx=lapply(sets,function(s)which(symbols %in% s))
  for(k in seq_along(x$contrasts))for(mode in c('default_0.01','conservative_0.10','estimated')){
    cor=if(mode=='estimated')NA else if(mode=='default_0.01').01 else .1
    z=camera(x$voom,index=idx,design=x$design,contrast=x$contrasts[[k]],inter.gene.cor=cor,sort=FALSE)
    z$set=rownames(z);z$study=sid;z$contrast=names(x$contrasts)[k];z$correlation_model=mode
    if(!'Correlation' %in% names(z))z$Correlation=cor
    z$coverage_status=ifelse(z$NGenes>=10,'at_least_10_genes','insufficient_gene_coverage')
    z$PValue[z$NGenes<10|!is.finite(z$PValue)]=1
    all[[length(all)+1]]=z
  }
}
a=do.call(rbind,all)
a$family=ifelse(a$study %in% c('OSD-561','OSD-562'),'core_50_sets_28_contrasts',ifelse(grepl('^OSD352_',a$study),'celltypes_50_sets_5_types','independent_bulk_50_sets_3_studies'))
a$q_full_family=ave(a$PValue,interaction(a$family,a$correlation_model),FUN=function(p)p.adjust(p,'BH'))
write.csv(a,file.path(B,'tables/pathway_CAMERA_complete.csv'),row.names=FALSE)
print(with(a,table(family,q_full_family<.05)))
