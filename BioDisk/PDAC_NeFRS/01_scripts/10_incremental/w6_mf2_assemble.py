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
def build(rows,stub,MARG=0.8,CGAP=1.2,RGAP=0.4,LBAND=4.3,NOLAB=0.4,RGAP_LAST=None):
    if RGAP_LAST is None: RGAP_LAST=RGAP
    imgs={n:load(n) for r in rows for n,_ in r}
    wmm=lambda n: imgs[n].size[0]/SCALE; hmm=lambda n: imgs[n].size[1]/SCALE
    bands=[LBAND if any(lab for _,lab in r) else NOLAB for r in rows]
    gaps=[RGAP]*(len(rows)-1)
    if gaps: gaps[-1]=RGAP_LAST              # 最后一个行间(进入 B/C)单独控制
    row_w=[sum(wmm(n) for n,_ in r)+CGAP*(len(r)-1) for r in rows]; content_w=max(row_w)
    cw=mm(content_w+2*MARG)
    ch=mm(sum(max(hmm(n) for n,_ in r)+bands[k] for k,r in enumerate(rows))+sum(gaps)+2*MARG)
    cv=Image.new("RGB",(cw,ch),"white"); dr=ImageDraw.Draw(cv); y=MARG
    for k,r in enumerate(rows):
        x=MARG; ytop=y+bands[k]
        for n,lab in r:
            cv.paste(imgs[n],(mm(x),mm(ytop)))
            if lab: dr.text((mm(x),mm(y)),lab,fill="black",font=FNT)
            x+=wmm(n)+CGAP
        y+=max(hmm(n) for n,_ in r)+bands[k]+(gaps[k] if k<len(rows)-1 else 0)
    cv.save(os.path.join(OUT,stub+".png"),dpi=(600,600))
    print(f"{stub}: {cv.size} = {cv.size[0]/SCALE:.1f} x {cv.size[1]/SCALE:.1f} mm")
    f=2200/cv.size[1]; cv.resize((int(cv.size[0]*f),2200)).save(os.path.join(OUT,stub+"_preview.png"))
build([
 [("MF2_KMc_PACA_AU","A"),("MF2_KMc_TCGA_PAAD",""),("MF2_KMc_E_MTAB_6134",""),("MF2_KMc_GSE71729","")],
 [("MF2_KMc_GSE85916",""),("MF2_KMc_GSE62452",""),("MF2_KMc_GSE57495",""),("MF2_KMc_GSE21501","")],
 [("MF2_forest","B"),("MF2_cindex","C")],
], "MF2_main", RGAP_LAST=-1.2)
