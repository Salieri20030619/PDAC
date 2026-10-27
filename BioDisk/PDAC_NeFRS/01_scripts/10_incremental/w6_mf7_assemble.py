#!/usr/bin/env python3
from PIL import Image, ImageDraw, ImageFont
import os
PAN="03_figures/panels"; OUT="03_figures"; SCALE=600/25.4
def mm(v): return int(round(v*SCALE))
def load(n): return Image.open(os.path.join(PAN,n+".png")).convert("RGB")
LET_PX=round(12/72*600); FNT=None
for p in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf","/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf","/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
    if os.path.exists(p): FNT=ImageFont.truetype(p,LET_PX); break
if FNT is None: FNT=ImageFont.load_default()
def build(rows,stub,MARG=0.9,CGAP=1.2,RGAP=0.8,LBAND=4.3):
    imgs={n:load(n) for r in rows for n,_ in r}
    wmm=lambda n: imgs[n].size[0]/SCALE; hmm=lambda n: imgs[n].size[1]/SCALE
    row_w=[sum(wmm(n) for n,_ in r)+CGAP*(len(r)-1) for r in rows]; content_w=max(row_w)
    cw=mm(content_w+2*MARG)
    ch=mm(sum(max(hmm(n) for n,_ in r)+LBAND for r in rows)+RGAP*(len(rows)-1)+2*MARG)
    cv=Image.new("RGB",(cw,ch),"white"); dr=ImageDraw.Draw(cv); y=MARG
    for r in rows:
        x=MARG; ytop=y+LBAND
        for n,lab in r:
            cv.paste(imgs[n],(mm(x),mm(ytop)))
            if lab: dr.text((mm(x),mm(y)),lab,fill="black",font=FNT)
            x+=wmm(n)+CGAP
        y+=max(hmm(n) for n,_ in r)+LBAND+RGAP
    cv.save(os.path.join(OUT,stub+".png"),dpi=(600,600))
    print(f"{stub}: {cv.size} = {cv.size[0]/SCALE:.1f} x {cv.size[1]/SCALE:.1f} mm")
    f=2200/cv.size[1]; cv.resize((int(cv.size[0]*f),2200)).save(os.path.join(OUT,stub+"_preview.png"))
build([
 [("MF7_A_cor","A"),("MF7_B_lollipop","B"),("MF7_C_desert","C")],
 [("MF7_D_gsea","D")],
 [("MF7_E_curveUp","E"),("MF7_F_curveDn","F"),("MF7_G_kegg","G")],
], "MF7_main")
