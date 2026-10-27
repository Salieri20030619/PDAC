#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })
R2<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/02_results"
ff<-function(pat) list.files(R2, pattern=pat, recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
show<-function(f,n=15){cat("\n----",sub(R2,"",f),"----\n")
  x<-tryCatch(fread(f),error=function(e)NULL)
  if(is.null(x)){cat("(非表格/读不了)\n");return(invisible())}
  cat("列:",paste(names(x),collapse=", "),"\n"); print(head(x,n))}

cat("################ 0. 02_results 全部 tsv/csv/rds 清单 ################\n")
print(list.files(R2, pattern="\\.(tsv|csv|csv\\.gz|rds)$", recursive=TRUE))

cat("\n################ 1. Nomogram 自助校正 C-index ################\n")
for(nm in c(file.path(R2,"05_nomogram/nomogram_models.rds"),
            file.path(R2,"05_nomogram/nomogram_calibration.rds"))){
  if(file.exists(nm)){m<-tryCatch(readRDS(nm),error=function(e)NULL)
    if(!is.null(m)&&!is.null(m$optimism_corrected_c)) cat(basename(nm),": optimism_corrected_c =",round(m$optimism_corrected_c,3),"\n")}}
cat("(若上面没出值,看这些候选:)\n",paste(ff("nomogram.*rds"),collapse="\n"),"\n")

cat("\n################ 2. 免疫细胞-NeFRS 相关 (确认 myeloid/Treg/Tfh 负、Th2 正) ################\n")
for(f in ff("immune.*cor.*tsv|.*nefrs.*immune.*tsv|immune.*nefrs.*tsv")){
  x<-tryCatch(fread(f),error=function(e)NULL)
  if(!is.null(x)&&"rho"%in%names(x)){cat("\n[",sub(R2,"",f),"]\n")
    cat("--- 最负(|rho|) ---\n"); print(head(x[order(rho)],12))
    cat("--- 最正 ---\n"); print(head(x[order(-rho)],8))}}

cat("\n################ 3. GSEA 通路 (NES/padj) ################\n")
for(f in head(ff("gsea.*tsv|gsea.*csv|.*enrich.*tsv|.*fgsea.*tsv"),6)) show(f,20)

cat("\n################ 4. 药敏 ################\n")
for(f in head(ff("drug.*tsv|drug.*csv|.*sensitiv.*tsv|.*oncopredict.*tsv"),6)) show(f,20)

cat("\n################ 5. 伪时间关联 (NeFRS vs pseudotime) ################\n")
for(f in head(ff("pseudo.*tsv|pseudo.*csv|.*assoc.*tsv|tradeseq.*tsv|fitgam.*tsv"),6)) show(f,15)

cat("\n################ 6. 空间相关 ################\n")
for(f in head(ff("spatial.*tsv|spatial.*csv"),6)) show(f,15)

cat("\n################ 7. 突变/基因组 ################\n")
for(f in head(ff("mut.*tsv|mut.*csv|.*genom.*tsv|.*cnv.*tsv|.*maf.*tsv|.*driver.*tsv"),6)) show(f,15)

cat("\n################ 8. signature 对比 C-index ################\n")
for(f in head(ff("sig.*compar.*tsv|.*compar.*cindex.*tsv|.*compar.*tsv|.*benchmark.*tsv"),8)) show(f,25)

cat("\n=== 完成。把以上完整输出贴回。若某项 0 结果,我据第0节清单再定位。 ===\n")
