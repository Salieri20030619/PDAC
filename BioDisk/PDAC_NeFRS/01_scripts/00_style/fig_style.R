# ==============================================================
# PDAC NeFRS · 统一出图规范  fig_style.R  (v3: 10pt)
#   source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
# --------------------------------------------------------------
# 【格式】 每 panel 同出 PDF(矢量,cairo) + PNG + TIFF(600dpi, TIFF=LZW)
# 【标题】 panel 内统一无标题(title/main);轴标签保留
# 【文字】 全部 10pt(轴/图例/分面);geom 内文字用常量 GT;panel 字母 12pt 加粗
# 【配色】 高危 #C00000 低危 #2166AC;发散 RdBu;分类 PAL_CAT(10色)
# 【显著性】 ns/* <0.05/** <0.01/*** <0.001;p 写 "p<0.001" 不写 0
# 【尺寸】 按 w_mm/h_mm 指定;大图全宽 180mm,拼装统一缩放保持一致
# ==============================================================
suppressPackageStartupMessages({ library(ggplot2) })
FONT <- "sans"; BASE_PT <- 10
mm2in <- function(mm) mm/25.4
GT <- BASE_PT/.pt                 # geom_text/annotate 用的字号
W1U<-42; W2U<-87; W3U<-132; WFULL<-180

COL_HIGH<-"#C00000"; COL_LOW<-"#2166AC"
PAL_RISK<-c(Low=COL_LOW,High=COL_HIGH)
PAL_DIV <-c("#2166AC","#4393C3","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B")
PAL_CAT <-c("#4C72B0","#DD8452","#55A868","#C44E52","#8172B3","#937860","#DA8BC3","#8C8C8C","#CCB974","#64B5CD")

theme_pub <- function(base = BASE_PT){
  theme_classic(base_size=base, base_family=FONT) +
  theme(
    axis.line=element_line(linewidth=0.4,colour="black"),
    axis.ticks=element_line(linewidth=0.4,colour="black"),
    axis.text=element_text(size=base,colour="black"),
    axis.title=element_text(size=base,colour="black"),
    legend.title=element_text(size=base), legend.text=element_text(size=base),
    legend.key.size=unit(3.8,"mm"), legend.background=element_blank(),
    plot.title=element_blank(), plot.subtitle=element_blank(),
    plot.tag=element_text(size=12,face="bold"),
    plot.margin=margin(2,2,2,2,"mm"),
    strip.background=element_blank(), strip.text=element_text(size=base))
}
theme_set(theme_pub())

scale_risk_col <-function(...) scale_colour_manual(values=PAL_RISK,...)
scale_risk_fill<-function(...) scale_fill_manual(values=PAL_RISK,...)
scale_cat_col  <-function(...) scale_colour_manual(values=PAL_CAT,...)
scale_cat_fill <-function(...) scale_fill_manual(values=PAL_CAT,...)

save_panel <- function(plot, file, w_mm=90, h_mm=70, dpi=600){
  if(!is.null(plot$labels)){plot$labels$title<-NULL; plot$labels$subtitle<-NULL}
  stub<-sub("\\.[^.]+$","",file)
  ggsave(paste0(stub,".pdf"), plot, width=mm2in(w_mm), height=mm2in(h_mm), units="in", device=cairo_pdf)
  ggsave(paste0(stub,".png"), plot, width=mm2in(w_mm), height=mm2in(h_mm), units="in", dpi=dpi, bg="white")
  ggsave(paste0(stub,".tiff"),plot, width=mm2in(w_mm), height=mm2in(h_mm), units="in", dpi=dpi, bg="white",
         device="tiff", compression="lzw")
  message(sprintf("saved 3fmt: %s (%gx%g mm)",stub,w_mm,h_mm)); invisible(stub)
}
save_base <- function(draw, file, w_mm=90, h_mm=70, dpi=600){
  stub<-sub("\\.[^.]+$","",file)
  render<-function(openfun){openfun()
    par(mar=c(3.2,3.4,0.8,0.8),mgp=c(2,0.6,0),tcl=-0.25,lwd=0.7,las=1,cex.axis=1,cex.lab=1); draw(); dev.off()}
  render(function() cairo_pdf(paste0(stub,".pdf"), width=mm2in(w_mm), height=mm2in(h_mm), family=FONT, pointsize=BASE_PT))
  render(function() png(paste0(stub,".png"),  width=w_mm, height=h_mm, units="mm", res=dpi, family=FONT, pointsize=BASE_PT, type="cairo"))
  render(function() tiff(paste0(stub,".tiff"), width=w_mm, height=h_mm, units="mm", res=dpi, family=FONT, pointsize=BASE_PT, type="cairo", compression="lzw"))
  message("saved 3fmt: ",stub)
}
panel_tag <- function(lab) mtext(lab, side=3, adj=-0.16, line=0.2, font=2, cex=12/BASE_PT)
sig_star  <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns")))
fmt_p     <- function(p) ifelse(p<0.001,"p<0.001",sprintf("p=%.3f",p))
if(sys.nframe()==0) cat("fig_style.R v3 OK | ",BASE_PT,"pt | pdf+png+tiff | 无标题\n")
