#!/usr/bin/env python3
import subprocess, os
from PIL import Image, ImageDraw, ImageFont
FIG="03_figures"; TMP="/tmp/supp2"; os.makedirs(TMP,exist_ok=True)
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
        im=pdf2png(pdf,stub+"_"+lab); h=int(im.size[1]*panel_w/im.size[0]); imgs.append((im.resize((panel_w,h)),lab))
    nrow=(len(imgs)+ncol-1)//ncol
    row_h=[max(imgs[r*ncol+c][0].size[1] for c in range(ncol) if r*ncol+c<len(imgs)) for r in range(nrow)]
    W=pad*2+ncol*panel_w+(ncol-1)*gap; H=pad*2+sum(row_h)+(nrow-1)*gap+nrow*label_h
    cv=Image.new("RGB",(W,H),"white"); dr=ImageDraw.Draw(cv); y=pad
    for r in range(nrow):
        x=pad
        for c in range(ncol):
            i=r*ncol+c
            if i>=len(imgs): break
            im,lab=imgs[i]; dr.text((x,y),lab,fill="black",font=FNT); cv.paste(im,(x,y+label_h)); x+=panel_w+gap
        y+=label_h+row_h[r]+gap
    cv.save(os.path.join(FIG,stub+".png"),dpi=(300,300))
    f=2200/cv.size[1]; cv.resize((int(cv.size[0]*f),2200)).save(os.path.join(FIG,stub+"_preview.png"))
    print(f"{stub}: {cv.size[0]}x{cv.size[1]}")
build_supp([
 ("Fig13_immune/Fig13C_adaptive_innate_box.pdf","A"),
 ("Fig13_immune/Fig13D_immune_boxplot.pdf","B"),
],2,"FigS7_immune")
build_supp([
 ("Fig7_spatial/Fig7C_microregion_coupling.pdf","A"),
 ("panels/S8_moran.pdf","B"),
 ("panels/S8_treat.pdf","C"),
],3,"FigS8_spatial")
