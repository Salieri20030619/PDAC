#!/usr/bin/env Rscript
# MF2 panel B 重做:per-SD Meta 森林(WYSIWYG 88x60mm)
# 方案B(默认): 仅 7 验证 + Pooled 菱形, 线性轴 1.0-2.5
# 方案A(SHOW_TRAIN<-TRUE): log 轴, 顶部加 PACA-AU 训练红点 + 分隔线
ROOT <- "/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"
source(file.path(ROOT,"01_scripts/00_style/fig_style.R"))
META <- file.path(ROOT,"02_results/20_meta")
FIGP <- file.path(ROOT,"03_figures/panels")

SHOW_TRAIN <- FALSE   # <- 改 TRUE 切方案A(log轴+保留训练集)

per  <- read.delim(file.path(META,"meta_per_cohort.tsv"), stringsAsFactors=FALSE, check.names=FALSE)
pool <- read.delim(file.path(META,"meta_pooled.tsv"),     stringsAsFactors=FALSE, check.names=FALSE)
cat("per cols :", paste(colnames(per),collapse=", "), "\n")
cat("pool cols:", paste(colnames(pool),collapse=", "), "\n")

# 容错取列
gc1 <- function(df,cands){ h<-intersect(tolower(cands),tolower(colnames(df)))[1]
  colnames(df)[tolower(colnames(df))==h] }
co_c <- gc1(per,c("cohort","Cohort","dataset")); hr_c <- gc1(per,c("HR","hr"))
lo_c <- gc1(per,c("lo","low","ci_lo","HR_lo")); hi_c <- gc1(per,c("hi","high","ci_hi","HR_hi"))

is_train <- grepl("PACA", per[[co_c]], ignore.case=TRUE)
val <- per[!is_train,]
val <- val[order(val[[hr_c]], decreasing=TRUE),]   # HR 大在上
trn <- per[is_train,]

pv <- pool[grepl("val", pool$set, ignore.case=TRUE),][1,]
nm <- function(x) gsub("_","-",x)

# ---- 坐标 ----
if(SHOW_TRAIN){ hrmin<-0.95; hrmax<-4.0; ticks<-c(1,2,3,4); uselog<-TRUE
} else { hrmin<-0.9; hrmax<-2.6; ticks<-c(1.0,1.5,2.0,2.5); uselog<-FALSE }
xL<-0.30; xR<-0.66
tf <- if(uselog) log else identity
fx <- function(hr) xL + (tf(hr)-tf(hrmin))/(tf(hrmax)-tf(hrmin))*(xR-xL)

draw <- function(){
  par(mai=c(0.02,0.02,0.02,0.02))
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")

  rows <- val[[co_c]]; rhr<-val[[hr_c]]; rlo<-val[[lo_c]]; rhi<-val[[hi_c]]
  rcol <- rep("black",nrow(val))
  if(SHOW_TRAIN && nrow(trn)>0){
    rows<-c(trn[[co_c]],rows); rhr<-c(trn[[hr_c]],rhr)
    rlo<-c(trn[[lo_c]],rlo);   rhi<-c(trn[[hi_c]],rhi)
    rcol<-c(COL_HIGH,rcol)
  }
  M <- length(rhr)
  yT<-0.88; yB<-0.42; step<-(yT-yB)/max(M-1,1)
  ry <- yT-(0:(M-1))*step
  py <- yB-1.35*step            # pooled 行
  ax <- py-1.25*step            # x 轴

  # 参考线 x=1
  segments(fx(1),ax,fx(1),yT+0.5*step,lty=2,col="grey55",lwd=1)

  # 各队列
  for(i in seq_len(M)){
    nmlab <- nm(rows[i]); if(SHOW_TRAIN&&i==1) nmlab<-paste0(nmlab," (train)")
    text(0.01,ry[i],nmlab,adj=0,cex=0.70,xpd=NA)
    segments(fx(max(rlo[i],hrmin)),ry[i],fx(min(rhi[i],hrmax)),ry[i],col=rcol[i],lwd=1.3,xpd=NA)
    points(fx(rhr[i]),ry[i],pch=15,col=rcol[i],cex=1.05,xpd=NA)
    text(0.99,ry[i],sprintf("%.2f (%.2f-%.2f)",rhr[i],rlo[i],rhi[i]),adj=1,cex=0.65,xpd=NA)
  }
  # 训练/验证分隔线(方案A)
  if(SHOW_TRAIN) segments(0,ry[1]-0.5*step,1,ry[1]-0.5*step,col="grey75",lwd=0.8,xpd=NA)

  # Pooled 菱形(val7)
  hd<-0.55*step
  polygon(c(fx(pv$lo),fx(pv$HR),fx(pv$hi),fx(pv$HR)),
          c(py,py+hd,py,py-hd),col=COL_LOW,border=COL_LOW,xpd=NA)
  text(0.01,py,"Pooled (7 val)",adj=0,cex=0.70,font=2,xpd=NA)
  text(0.99,py,sprintf("%.2f (%.2f-%.2f)",pv$HR,pv$lo,pv$hi),adj=1,cex=0.65,font=2,xpd=NA)
  # I2 标注:贴 Pooled 行,措辞明确
  text((xL+xR)/2, py+1.05*step,
       expression(paste("Random-effects (7 val): ",italic(I)^2,"=0%")),
       cex=0.62,xpd=NA)

  # x 轴
  segments(fx(hrmin),ax,fx(hrmax),ax,lwd=1,xpd=NA)
  for(t in ticks){ segments(fx(t),ax,fx(t),ax-0.012,lwd=1,xpd=NA)
    text(fx(t),ax-0.05,formatC(t,format="f",digits=1),cex=0.62,xpd=NA) }
  text((xL+xR)/2, ax-0.13,"HR per SD (95% CI)",cex=0.75,xpd=NA)
}

save_base(draw, file.path(FIGP,"MF2_forest"), 88, 60)
cat("done -> MF2_forest (SHOW_TRAIN=",SHOW_TRAIN,")\n")
