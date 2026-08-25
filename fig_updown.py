import csv, glob, os, random
import numpy as np, openslide
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

P = 224
allsvs = glob.glob(r"E:\KIRC\slides_tumor\**\*.svs", recursive=True)
def find(case):
    h=[f for f in allsvs if os.path.basename(f).startswith(case)]
    return h[0] if h else None

def tissue_positions(s):
    W,H=s.dimensions; ds=32
    gray=np.array(s.get_thumbnail((W//ds,H//ds)).convert("L")); tis=gray<210
    sx,sy=gray.shape[1]/W, gray.shape[0]/H
    pos=[]
    for y in range(0,H-P,P):
        for x in range(0,W-P,P):
            if tis[min(int(y*sy),gray.shape[0]-1),min(int(x*sx),gray.shape[1]-1)]: pos.append((x,y))
    return pos

def good_patches(case, k=2):
    # aday patch'lerden en doku-dolu/keskin k tanesini sec
    path=find(case); s=openslide.OpenSlide(path)
    pos=tissue_positions(s); random.seed(1)
    cand=random.sample(pos, min(14,len(pos)))
    scored=[]
    for p in cand:
        im=np.array(s.read_region(p,0,(P,P)).convert("RGB"))
        g=im.mean(2); score=(g<210).mean()*g.std()
        scored.append((score,im))
    scored.sort(key=lambda t:-t[0])
    s.close()
    return [im for _,im in scored[:k]]

rows=list(csv.DictReader(open(r"E:\KIRC\gene_morph_updown.csv")))
seen=[]; [seen.append(r["satir"]) for r in rows if r["satir"] not in seen]
gruplar={x:[r for r in rows if r["satir"]==x] for x in seen}
etik={x:gruplar[x][0]["etiket"] for x in seen}

# 4 satir x 4 sutun (2 hasta x 2 patch). Hasta ciftlerini ince bosluk ile ayir.
fig=plt.figure(figsize=(12,12.5))
gs=gridspec.GridSpec(4,4, hspace=0.32, wspace=0.08)

for ri,x in enumerate(seen):
    cs=gruplar[x][:2]
    col=0
    for ci,c in enumerate(cs):
        pats=good_patches(c["submitter_id"],2)
        for pj,im in enumerate(pats):
            a=fig.add_subplot(gs[ri,col]); a.imshow(im); a.set_xticks([]); a.set_yticks([])
            if col==0: a.set_ylabel(etik[x],fontsize=10)
            if pj==0: a.set_title(c["submitter_id"],fontsize=8.5, loc="left")
            col+=1

fig.suptitle("Morphology by up/down signature-gene expression (best image-predicted genes)",
             fontsize=14, weight="bold", y=0.925)
fig.text(0.5,0.045,"2 patients x 2 patches per group   |   DLX4: up-regulated/risky   ·   AJAP1: down-regulated/protective",
         fontsize=10, ha="center", color="#5B7682", style="italic")
plt.savefig(r"E:\KIRC\Fig_updown_morphology.png",dpi=170,bbox_inches="tight")
print("SAVED")
