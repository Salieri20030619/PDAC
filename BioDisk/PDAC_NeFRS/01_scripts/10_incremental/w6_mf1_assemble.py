#!/usr/bin/env python3
# MF1 主图(7张,3行) + 补充S1(H/I)。服务器端 1:1 拼，180mm 宽
from PIL import Image, ImageDraw, ImageFont
import os
PAN="03_figures/panels"; OUT="03_figures"
SCALE=600/25.4
def mm(v): return int(round(v*SCALE))
def load(n): return Image.open(os.path.join(PAN,n+".png")).convert("RGB")
FNT=None
for p in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
          "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
    if os.path.exists(p): FNT=ImageFont.truetype(p,mm(5.0)); break
if FNT is None: FNT=ImageFont.load_default()

def build(rows, stub, MARG=2, CGAP=2, RGAP=2, LBAND=3):
    imgs={n:load(n) for r in rows for n,_ in r}
    wmm=lambda n: imgs[n].size[0]/SCALE
    hmm=lambda n: imgs[n].size[1]/SCALE
    row_w=[sum(wmm(n) for n,_ in r)+CGAP*(len(r)-1) for r in rows]
    content_w=max(row_w)
    cw=mm(content_w+2*MARG)
    row_h=[max(hmm(n) for n,_ in r)+LBAND for r in rows]
    ch=mm(sum(row_h)+RGAP*(len(rows)-1)+2*MARG)
    cv=Image.new("RGB",(cw,ch),"white"); dr=ImageDraw.Draw(cv)
    y=MARG
    for ri,r in enumerate(rows):
        x=MARG+(content_w-row_w[ri])/2.0; ytop=y+LBAND
        for n,lab in r:
            cv.paste(imgs[n],(mm(x),mm(ytop)))
            dr.text((max(mm(x)-2,2),mm(y)),lab,fill="black",font=FNT)
            x+=wmm(n)+CGAP
        y+=row_h[ri]+RGAP
    cv.save(os.path.join(OUT,stub+".png"),dpi=(600,600))
    print(f"{stub}: {cv.size} = {cv.size[0]/SCALE:.1f} x {cv.size[1]/SCALE:.1f} mm")
    f=2200/cv.size[1]; pv=cv.resize((int(cv.size[0]*f),2200))
    pv.save(os.path.join(OUT,stub+"_preview.png"))

# 主图：A B C / D E / F G
build([
 [("MF1_venn","A"),("MF1_funnel","B"),("MF1_Bopt2_cvcurve","C")],
 [("MF1_E_mime","D"),("MF1_coef_top30","E")],
 [("MF1_I_neural_heat","F"),("MF1_J_ferro_heat","G")],
], "MF1_main")

# 补充 S1：H 8队列分布 + I 方向一致性（一行）
build([
 [("MF1_score_dist","A"),("MF1_D_consist","B")],
], "FigS1")
