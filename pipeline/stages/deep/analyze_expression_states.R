# Post-result forebrain state analysis, full planned target/contrast families.
library(limma); library(edgeR)
B=Sys.getenv("OSDR_DEEP_DIR")
args=commandArgs(trailingOnly=TRUE);if(length(args)>0)B=args[1]
map=read.csv(file.path(Sys.getenv("OSDR_CORE_DIR"), "tables/verified_panel_mapping.csv"),stringsAsFactors=FALSE)
ann=read.csv(file.path(B,'prepared/state_annotation.csv'),stringsAsFactors=FALSE)
md=read.csv(file.path(B,'tables/state_mouse_QC.csv'),stringsAsFactors=FALSE)
tests=c('flight_age3','flight_age8','flight_balanced_age','flight_by_age')
cs=cbind(c(-1,0,1,0),c(0,-1,0,1),c(-.5,-.5,.5,.5),c(1,-1,-1,1)); colnames(cs)=tests
rows=list(); qc=list()
for(K in c(4,6)){
  counts=as.matrix(read.csv(gzfile(file.path(B,'prepared',paste0('state_K',K,'_pseudobulk.csv.gz'))),row.names=1,check.names=FALSE))
  for(state in seq_len(K)){
    s=md[md$K==K & md$state==state,,drop=FALSE]
    s$group=factor(paste(s$flight,s$age,sep='_'),levels=c('GC_3','GC_8','FLT_3','FLT_8'))
    ok=nrow(s)==8 && all(s$n_cells>=50) && all(table(s$group)==2)
    if(ok){
      X=counts[,s$key,drop=FALSE]; keep=rowSums(X>=10)>=3; X=X[keep,,drop=FALSE]
      D=model.matrix(~0+group,s); colnames(D)=levels(s$group)
      d=calcNormFactors(DGEList(X),method='TMM')
      v=voom(d,D,plot=FALSE); fit=eBayes(contrasts.fit(lmFit(v,D),cs),robust=TRUE)
      d=estimateDisp(d,D,robust=TRUE); ql=glmQLFit(d,D,robust=TRUE)
      write.csv(v$E,gzfile(file.path(B,'prepared',paste0('state_K',K,'_S',state,'_voom_logCPM.csv.gz'))))
      qc[[length(qc)+1]]=data.frame(K=K,state=state,model_available=TRUE,min_cells_per_mouse=min(s$n_cells),n_mice=8,n_filtered_genes=nrow(X),design_rank=qr(D)$rank,residual_df=nrow(D)-qr(D)$rank)
    } else qc[[length(qc)+1]]=data.frame(K=K,state=state,model_available=FALSE,min_cells_per_mouse=min(s$n_cells),n_mice=nrow(s),n_filtered_genes=0,design_rank=NA,residual_df=NA)
    for(k in seq_along(tests))for(method in c('limma_voom','edgeR_QL')){
      if(ok){
        if(method=='limma_voom'){
          res=topTable(fit,coef=k,number=Inf,sort.by='none')
          se=fit$stdev.unscaled[,k]*sqrt(fit$s2.post)
          res$ci_low=res$logFC-qt(.975,fit$df.total)*se
          res$ci_high=res$logFC+qt(.975,fit$df.total)*se
        }else{
          res=topTags(glmQLFTest(ql,contrast=cs[,k]),n=Inf,sort.by='none')$table
          res$P.Value=res$PValue;res$ci_low=NA;res$ci_high=NA
        }
        res$SYMBOL=ann$SYMBOL[match(rownames(res),ann$ENSEMBL)]
        write.csv(res,gzfile(file.path(B,'tables',paste0('state_K',K,'_S',state,'_',method,'_',tests[k],'_genome.csv.gz'))))
      }
      for(i in seq_len(nrow(map))){
        hits=if(ok)which(!is.na(res$SYMBOL)&res$SYMBOL==map$mouse_symbol_primary[i])else integer(); unique=length(hits)==1
        rows[[length(rows)+1]]=data.frame(K=K,state=state,method=method,contrast=tests[k],analyte=map$analyte[i],mouse_symbol=map$mouse_symbol_primary[i],
          log2FC=if(unique)res$logFC[hits] else NA,ci_low=if(unique)res$ci_low[hits]else NA,ci_high=if(unique)res$ci_high[hits]else NA,
          pvalue=if(unique)res$P.Value[hits]else 1,n_mice=nrow(s),min_cells_per_mouse=min(s$n_cells),
          status=if(unique)'tested' else if(!ok)'insufficient_all_mouse_cell_coverage' else 'unavailable_or_nonunique')
      }
    }
    cat('K=',K,' state=',state,' available=',ok,' minimum cells=',min(s$n_cells),'\n')
  }
}
a=do.call(rbind,rows);a$pvalue[!is.finite(a$pvalue)]=1
a$q_full_state_family=ave(a$pvalue,interaction(a$K,a$method),FUN=function(p)p.adjust(p,'BH'))
a$planned_family_size=151*a$K*4
write.csv(a,file.path(B,'tables/state_panel_complete.csv'),row.names=FALSE)
write.csv(do.call(rbind,qc),file.path(B,'tables/state_model_QC.csv'),row.names=FALSE)
print(with(a,table(K,method,q_full_state_family<.05)))
