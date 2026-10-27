#!/usr/bin/env Rscript
# W6 panel 统一重出 (10pt): MF2K/MF3H/MF3I/MF4EFG/MF6J/MF7J
source("/mnt/BioDisk/wangmingyuan/PDAC_NeFRS/01_scripts/00_style/fig_style.R")
suppressPackageStartupMessages({ library(survival); library(glmnet); library(data.table); library(ggplot2) })
ROOT<-"/mnt/BioDisk/wangmingyuan/PDAC_NeFRS"; PFIG<-file.path(ROOT,"03_figures/panels")
W2IN<-file.path(ROOT,"02_results/02_w2_inputs"); MIME<-file.path(ROOT,"02_results/04_mime")

## MF2K: 8 队列 C-index
cohort_data<-readRDS(file.path(W2IN,"cohort_data_8.rds"))
mi<-readRDS(file.path(MIME,"mime_input_list.rds")); rs0<-readRDS(file.path(MIME,"mime_results.rds"))
mg<-colnames(mi[[1]])[-(1:3)]; ml<-rs0$ml.res[["Enet[α=0.1]"]]
crs<-function(em){X<-matrix(0,ncol(em),length(mg),dimnames=list(colnames(em),mg))
  for(g in intersect(mg,rownames(em)))X[,g]<-em[g,]; setNames(as.numeric(predict(ml,newx=X,s="lambda.min")),colnames(em))}
ci<-data.frame()
for(cn in names(cohort_data)){co<-cohort_data[[cn]];cl<-co$clinical;r<-crs(co$expression)
  d<-data.frame(t=cl$os_time_month,s=cl$os_status,z=as.numeric(r[cl$sample_id]))
  d<-d[!is.na(d$t)&d$t>0&!is.na(d$s)&!is.na(d$z),]
  ci<-rbind(ci,data.frame(cohort=cn,C=concordance(coxph(Surv(t,s)~z,data=d))$concordance,
    set=ifelse(cn=="PACA_AU","Training","Validation")))}
ci$cohort<-factor(ci$cohort,levels=ci$cohort[order(-ci$C)])
p<-ggplot(ci,aes(cohort,C,fill=set))+geom_col(width=0.72)+
  geom_hline(yintercept=0.5,linetype=2,linewidth=0.4,colour="grey50")+
  geom_text(aes(label=sprintf("%.3f",C)),vjust=-0.5,size=GT)+
  scale_fill_manual(values=c(Training=COL_HIGH,Validation=COL_LOW))+
  coord_cartesian(ylim=c(0.45,0.82))+labs(x=NULL,y="C-index",fill=NULL)+
  theme(axis.text.x=element_text(angle=35,hjust=1),legend.position=c(0.88,0.9))
save_panel(p,file.path(PFIG,"MF2K_cindex_8cohort"),w_mm=120,h_mm=85)

## MF3H: 超越分期 C-index
bs<-fread(file.path(ROOT,"02_results/21_nri/beyond_stage_cindex.tsv"))
dh<-data.frame(model=factor(c("AJCC stage","Stage + NeFRS"),levels=c("AJCC stage","Stage + NeFRS")),C=c(bs$C_stage,bs$C_stage_NeFRS))
p<-ggplot(dh,aes(model,C,fill=model))+geom_col(width=0.6)+
  geom_hline(yintercept=0.5,linetype=2,linewidth=0.4,colour="grey50")+
  geom_text(aes(label=sprintf("%.3f",C)),vjust=-0.5,size=GT)+
  annotate("text",x=1.5,y=max(dh$C)+0.05,label=sprintf("ΔC = +%.3f",bs$dC),fontface=2,size=GT)+
  scale_fill_manual(values=c("AJCC stage"="grey70","Stage + NeFRS"=COL_HIGH),guide="none")+
  coord_cartesian(ylim=c(0.45,0.78))+labs(x=NULL,y="C-index")
save_panel(p,file.path(PFIG,"MF3H_beyondstage_C"),w_mm=72,h_mm=74)

## MF3I: IDI/NRI (短标题)
di<-data.frame(measure=factor(rep(c("IDI","Continuous NRI","Median improve"),2),levels=c("IDI","Continuous NRI","Median improve")),
  t0=factor(rep(c("12 mo","24 mo"),each=3),levels=c("12 mo","24 mo")),
  est=c(0.152,0.376,0.146,0.132,0.280,0.160),lo=c(0.107,0.291,0.106,0.094,0.206,0.104),hi=c(0.205,0.473,0.215,0.173,0.366,0.233))
p<-ggplot(di,aes(measure,est,colour=t0))+geom_hline(yintercept=0,linewidth=0.4,colour="grey60")+
  geom_point(position=position_dodge(0.5),size=2.2)+
  geom_errorbar(aes(ymin=lo,ymax=hi),position=position_dodge(0.5),width=0.18,linewidth=0.6)+
  scale_colour_manual(values=c("12 mo"=COL_LOW,"24 mo"=COL_HIGH))+scale_y_continuous(limits=c(0,NA))+
  labs(x=NULL,y="Improvement vs stage",colour=NULL)+
  theme(axis.text.x=element_text(angle=15,hjust=1),legend.position=c(0.85,0.92))
save_panel(p,file.path(PFIG,"MF3I_idi_nri"),w_mm=110,h_mm=82)

## MF4 E/F/G 散点
ctl<-fread(file.path(ROOT,"02_results/18_controls/controls_scores_persample.tsv"))
scat<-function(xv,xlab,file){d<-ctl[is.finite(get(xv))&is.finite(NeFRS)]
  ct<-suppressWarnings(cor.test(d[[xv]],d$NeFRS,method="spearman"))
  ann<-sprintf("Spearman rho = %.2f, %s",as.numeric(ct$estimate),fmt_p(ct$p.value))
  p<-ggplot(d,aes(.data[[xv]],NeFRS))+geom_point(aes(colour=cohort),size=1.1,alpha=0.45)+
    geom_smooth(method="lm",se=FALSE,colour="black",linewidth=0.7)+scale_colour_manual(values=PAL_CAT)+
    annotate("text",x=-Inf,y=Inf,hjust=-0.05,vjust=1.5,label=ann,size=GT)+
    labs(x=xlab,y="NeFRS")+theme(legend.position="none")
  save_panel(p,file,w_mm=95,h_mm=82)}
scat("purity","Tumor purity (ESTIMATE, z)",file.path(PFIG,"MF4E_NeFRS_vs_purity"))
scat("prolif","Proliferation score (z)",file.path(PFIG,"MF4F_NeFRS_vs_prolif"))
scat("basal_prob","PurIST Basal probability",file.path(PFIG,"MF4G_NeFRS_vs_basal"))

## MF6J: 三 KO 趋同 Venn
KOD<-file.path(ROOT,"02_results/10_knockout/Results_final")
gs<-function(f) unique(fread(file.path(KOD,f))$gene)
A<-gs("NRP1_sig_adjP005.csv");B<-gs("EGFR_sig_adjP005.csv");C<-gs("BMPR1B_sig_adjP005.csv")
nABC<-length(Reduce(intersect,list(A,B,C)));nAB<-length(intersect(A,B));nAC<-length(intersect(A,C));nBC<-length(intersect(B,C))
oA<-length(A)-nAB-nAC+nABC;oB<-length(B)-nAB-nBC+nABC;oC<-length(C)-nAC-nBC+nABC
oAB<-nAB-nABC;oAC<-nAC-nABC;oBC<-nBC-nABC
dc<-function(x,y,r,col){a<-seq(0,2*pi,length.out=200);polygon(x+r*cos(a),y+r*sin(a),col=adjustcolor(col,0.30),border=col,lwd=1.5)}
draw<-function(){par(mar=c(0.4,0.4,0.4,0.4));plot(NA,xlim=c(-1.9,1.9),ylim=c(-2,1.7),axes=FALSE,xlab="",ylab="",asp=1)
  dc(-0.5,0.35,0.95,COL_HIGH);dc(0.5,0.35,0.95,COL_LOW);dc(0,-0.55,0.95,"#55A868")
  tt<-function(x,y,n,cx=1)text(x,y,n,cex=cx)
  tt(-1,0.72,oA);tt(1,0.72,oB);tt(0,-1.25,oC);tt(0,0.78,oAB);tt(-0.64,-0.3,oAC);tt(0.64,-0.3,oBC);tt(0,0.03,nABC,1.25)
  text(-1.2,1.3,"NRP1 KO",font=2,col=COL_HIGH);text(1.2,1.3,"EGFR KO",font=2,col=COL_LOW);text(0,-1.6,"BMPR1B KO",font=2,col="#55A868")}
save_base(draw,file.path(PFIG,"MF6J_KO_convergence_venn"),w_mm=92,h_mm=92)

## MF7J: 免疫相关性热图 (缩写细胞名)
ic<-fread(file.path(ROOT,"02_results/14_immune/immune_nefrs_cor.tsv"))
ab<-c("Activated B cell"="Act. B","Activated CD4 T cell"="Act. CD4 T","Activated CD8 T cell"="Act. CD8 T",
"Activated dendritic cell"="Act. DC","CD56bright natural killer cell"="CD56bright NK","CD56dim natural killer cell"="CD56dim NK",
"Central memory CD4 T cell"="Tcm CD4","Central memory CD8 T cell"="Tcm CD8","Effector memory CD4 T cell"="Tem CD4",
"Effector memory CD8 T cell"="Tem CD8","Gamma delta T cell"="gd T","Immature B cell"="Imm. B","Immature dendritic cell"="Imm. DC",
"Mast cell"="Mast","Memory B cell"="Mem. B","Natural killer cell"="NK","Natural killer T cell"="NKT",
"Plasmacytoid dendritic cell"="pDC","Regulatory T cell"="Treg","T follicular helper cell"="Tfh",
"Type 1 T helper cell"="Th1","Type 17 T helper cell"="Th17","Type 2 T helper cell"="Th2")
ic[,cell2:=ifelse(cell %in% names(ab), ab[cell], cell)]
ic[,star:=sig_star(padj)]
ord<-ic[,.(m=mean(rho,na.rm=TRUE)),by=cell2][order(m),cell2]; ic[,cell2:=factor(cell2,levels=ord)]
lim<-max(abs(ic$rho),na.rm=TRUE)
p<-ggplot(ic,aes(cohort,cell2,fill=rho))+geom_tile(colour="white",linewidth=0.4)+
  geom_text(aes(label=star),size=GT,vjust=0.8)+
  scale_fill_gradient2(low=COL_LOW,mid="white",high=COL_HIGH,midpoint=0,limits=c(-lim,lim),name="Spearman\nrho")+
  labs(x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=35,hjust=1),axis.line=element_blank(),axis.ticks=element_blank())
save_panel(p,file.path(PFIG,"MF7J_immune_cor_heatmap"),w_mm=100,h_mm=148)
cat("\n全部 panel 已按 10pt 重出\n")
