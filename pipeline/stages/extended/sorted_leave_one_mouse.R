# Actual count-model deletion refits; no new p-values from repeated analyses.
library(limma);library(edgeR)
B=Sys.getenv("OSDR_EXTENDED_DIR");md=read.csv(file.path(B,'tables/GSE295428_sorted_cell_QC.csv'));ann=read.csv(file.path(B,'prepared/GSE295428_sorted_annotation.csv'));panel=read.csv(file.path(Sys.getenv("OSDR_CORE_DIR"), "tables/verified_panel_mapping.csv"))
X=as.matrix(read.csv(gzfile(file.path(B,'prepared/GSE295428_sorted_donor_pseudobulk.csv.gz')),row.names=1,check.names=FALSE));rows=list();models=list();contrast=c(-.5,-.5,.5,.5)
for(tissue in sort(unique(md$tissue[grepl('^Myeloid_',md$tissue)]))){
    full=md[md$tissue==tissue,,drop=FALSE];C=X[,full$sample,drop=FALSE];C=C[rowSums(C>=10)>=3,,drop=FALSE]
    for(removed in full$donor){
        s=full[full$donor!=removed,,drop=FALSE];s$group=factor(paste(s$flight,s$age,sep='_'),levels=c('GC_3','GC_8','FLT_3','FLT_8'));ok=all(table(s$group)>0)
        if(ok){
            D=model.matrix(~0+group,s);counts=C[,s$sample,drop=FALSE];d=calcNormFactors(DGEList(counts),method='TMM');v=voom(d,D,plot=FALSE);f=contrasts.fit(lmFit(v,D),matrix(contrast,ncol=1));d=estimateDisp(d,D,robust=TRUE);q=glmQLFit(d,D,robust=TRUE)
            estimates=list(limma_voom=f$coefficients[,1],edgeR_QL=drop(q$coefficients %*% contrast)/log(2));symbols=ann$SYMBOL[match(rownames(counts),ann$ENSEMBL)]
        }
        models[[length(models)+1]]=data.frame(tissue=tissue,removed_donor=removed,n_remaining_donors=nrow(s),status=if(ok)'refitted' else 'missing_age_exposure_design_cell')
        for(method in c('limma_voom','edgeR_QL'))for(i in seq_len(nrow(panel))){
            hit=if(ok)which(!is.na(symbols)&symbols==panel$mouse_symbol_primary[i]) else integer();available=length(hit)==1
            rows[[length(rows)+1]]=data.frame(tissue=tissue,removed_donor=removed,method=method,analyte=panel$analyte[i],mouse_symbol=panel$mouse_symbol_primary[i],log2FC=if(available)estimates[[method]][hit] else NA,status=if(available)'refitted' else if(!ok)'missing_age_exposure_design_cell' else 'unavailable_or_nonunique')
        }
        cat(tissue,removed,ok,'\n')
    }
}
write.csv(do.call(rbind,rows),file.path(B,'tables/sorted_leave_one_mouse_effects.csv'),row.names=FALSE);write.csv(do.call(rbind,models),file.path(B,'tables/sorted_leave_one_mouse_models.csv'),row.names=FALSE)
