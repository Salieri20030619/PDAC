#!/usr/bin/env python3
import subprocess, os
from PIL import Image, ImageDraw, ImageFont
FIG="03_figures"; TMP="/tmp/supp"; os.makedirs(TMP,exist_ok=True)
FNT=None
for p in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf","/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf","/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
    if os.path.exists(p): FNT=ImageFont.truetype(p,52); break
if FNT is None: FNT=ImageFont.load_default()
def pdf2png(pdf,name,dpi=300):
    o=os.path.join(TMP,name)
    subprocess.run(["pdftoppm","-png","-r",str(dpi),"-singlefile",os.path.join(FIG,pdf),o],check=True)
    return Image.open(o+".png").convert("RGB")
def build_supp(panels,ncol,stub,panel_w=1000,gap=44,pad=34,label_h=64):
    imgs=[]
    for pdf,lab in panels:
        im=pdf2png(pdf,stub+"_"+lab)
        h=int(im.size[1]*panel_w/im.size[0])
        imgs.append((im.resize((panel_w,h)),lab))
    nrow=(len(imgs)+ncol-1)//ncol
    row_h=[max(imgs[r*ncol+c][0].size[1] for c in range(ncol) if r*ncol+c<len(imgs)) for r in range(nrow)]
    W=pad*2+ncol*panel_w+(ncol-1)*gap
    H=pad*2+sum(row_h)+(nrow-1)*gap+nrow*label_h
    cv=Image.new("RGB",(W,H),"white"); dr=ImageDraw.Draw(cv); y=pad
    for r in range(nrow):
        x=pad
        for c in range(ncol):
            i=r*ncol+c
            if i>=len(imgs): break
            im,lab=imgs[i]
            dr.text((x,y),lab,fill="black",font=FNT)
            cv.paste(im,(x,y+label_h))
            x+=panel_w+gap
        y+=label_h+row_h[r]+gap
    cv.save(os.path.join(FIG,stub+".png"),dpi=(300,300))
    f=2200/cv.size[1]; cv.resize((int(cv.size[0]*f),2200)).save(os.path.join(FIG,stub+"_preview.png"))
    print(f"{stub}: {cv.size[0]}x{cv.size[1]}")

build_supp([
 ("Fig5_DCA/Fig5C_DCA_PACA_AU.pdf","A"),
 ("Fig5_DCA/Fig5C_DCA_E_MTAB_6134.pdf","B"),
 ("Fig5B_Calibration_PACA_AU.pdf","C"),
 ("Fig5B_Calibration_External.pdf","D"),
],2,"FigS3_clinical")

build_supp([
 ("Fig17_pancancer/Fig17_KM_HCC.pdf","A"),
 ("Fig17_pancancer/Fig17_KM_CRC.pdf","B"),
],2,"FigS4_pancancer")

build_supp([
 ("Fig6C_Werba_NeFRS_by_celltype.pdf","A"),
 ("Fig6B_NeFRS_vs_Hwang17_programs.pdf","B"),
 ("Fig8_pseudotime/Fig8A2_ductal_celltype.pdf","C"),
 ("Fig8_pseudotime/Fig8B2_ductal_NeFRS.pdf","D"),
 ("Fig8_pseudotime/Fig8C2_ductal_NeFRS_pt.pdf","E"),
 ("Fig8_pseudotime/Fig8E_nefrs_genes_smoothers.pdf","F"),
],3,"FigS5_singlecell")

build_supp([
 ("Fig10_TF/Fig10A_TF_heatmap.pdf","A"),
 ("Fig10_TF/Fig10B_keyTF_hilo.pdf","B"),
 ("Fig10_tf/Fig10B_v_quintile_heatmap.pdf","C"),
 ("Fig10_tf/Fig10B_v_quintile_line.pdf","D"),
 ("Fig10_tf/Fig10C_v_neural_dichotomy.pdf","E"),
 ("Fig11_cellchat/Fig11A_network_count.pdf","F"),
 ("Fig11_cellchat/Fig11B_neural_bubble.pdf","G"),
],3,"FigS6_mechanism")
