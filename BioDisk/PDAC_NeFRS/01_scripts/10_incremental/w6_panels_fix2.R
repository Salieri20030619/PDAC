#!/usr/bin/env Rscript
# 修订: MF4 E/F/G(顶部留白+图例) + MF6J(VennDiagram 重画)
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; PFIG<-file.path(ROOT,"03_figures/panels")

## ---- MF4 E/F/G: 顶部留白 + 队列图例 ----
ctl<-fread(file.path(ROOT,"02_results/18_controls/controls_scores_persample.tsv"))
scat<-function(xv,xlab,file){d<-ctl[is.finite(get(xv))&is.finite(NeFRS)]
  ct<-suppressWarnings(cor.test(d[[xv]],d$NeFRS,method="spearman"))
  ann<-sprintf("Spearman rho = %.2f, %s",as.numeric(ct$estimate),fmt_p(ct$p.value))
  p<-ggplot(d,aes(.data[[xv]],NeFRS))+
    geom_point(aes(colour=cohort),size=1.1,alpha=0.45)+
    geom_smooth(method="lm",se=FALSE,colour="black",linewidth=0.7)+
    scale_colour_manual(values=PAL_CAT,name=NULL)+
    scale_y_continuous(expand=expansion(mult=c(0.05,0.18)))+
    annotate("text",x=-Inf,y=Inf,hjust=-0.05,vjust=1.3,label=ann,size=GT)+
    labs(x=xlab,y="NeFRS")+
    guides(colour=guide_legend(override.aes=list(size=2,alpha=1)))+
    theme(legend.position="right",legend.key.size=unit(3.2,"mm"))
  save_panel(p,file,w_mm=128,h_mm=84)}
scat("purity","Tumor purity (ESTIMATE, z)",file.path(PFIG,"MF4E_NeFRS_vs_purity"))
scat("prolif","Proliferation score (z)",file.path(PFIG,"MF4F_NeFRS_vs_prolif"))
scat("basal_prob","PurIST Basal probability",file.path(PFIG,"MF4G_NeFRS_vs_basal"))

## ---- MF6J: VennDiagram 重画 ----
suppressPackageStartupMessages({ library(VennDiagram); library(grid) })
flog.threshold(ERROR, name="VennDiagramLogger")   # 关掉日志
KOD<-file.path(ROOT,"02_results/10_knockout/Results_final")
gs<-function(f) unique(fread(file.path(KOD,f))$gene)
A<-gs("NRP1_sig_adjP005.csv");B<-gs("EGFR_sig_adjP005.csv");C<-gs("BMPR1B_sig_adjP005.csv")
nAB<-length(intersect(A,B));nBC<-length(intersect(B,C));nAC<-length(intersect(A,C));nABC<-length(Reduce(intersect,list(A,B,C)))
cols<-c(COL_HIGH,COL_LOW,"#55A868")
venn_one<-function(open){ open(); grid.newpage()
  draw.triple.venn(area1=length(A),area2=length(B),area3=length(C),
    n12=nAB,n23=nBC,n13=nAC,n123=nABC,
    category=c("NRP1 KO","EGFR KO","BMPR1B KO"),
    fill=adjustcolor(cols,0.40), col=cols, lwd=1.6,
    cex=1.1, cat.cex=1.2, cat.col=cols, cat.fontface=2,
    fontfamily="sans", cat.fontfamily="sans", margin=0.10)
  dev.off() }
stub<-file.path(PFIG,"MF6J_KO_convergence_venn")
venn_one(function() cairo_pdf(paste0(stub,".pdf"),width=mm2in(95),height=mm2in(92),family=FONT,pointsize=BASE_PT))
venn_one(function() png(paste0(stub,".png"),width=95,height=92,units="mm",res=600,family=FONT,pointsize=BASE_PT,type="cairo"))
venn_one(function() tiff(paste0(stub,".tiff"),width=95,height=92,units="mm",res=600,family=FONT,pointsize=BASE_PT,type="cairo",compression="lzw"))
cat("\n修订完成: MF4E/F/G(留白+图例), MF6J(VennDiagram)\n")
