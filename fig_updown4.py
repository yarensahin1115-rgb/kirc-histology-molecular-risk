import csv, glob, os, random
import numpy as np, openslide
from scipy.ndimage import laplace
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

P=224
allsvs=glob.glob(r"E:\KIRC\slides_tumor\**\*.svs",recursive=True)
def find(c):
    h=[f for f in allsvs if os.path.basename(f).startswith(c)]; return h[0] if h else None
def tissue_pos(s):
    W,H=s.dimensions; ds=32
    g=np.array(s.get_thumbnail((W//ds,H//ds)).convert("L")); tis=g<210
    sx,sy=g.shape[1]/W,g.shape[0]/H; pos=[]
    for y in range(0,H-P,P):
        for x in range(0,W-P,P):
            if tis[min(int(y*sy),g.shape[0]-1),min(int(x*sx),g.shape[1]-1)]: pos.append((x,y))
    return pos
def sharpness(im):
    g=im.mean(2)
    return laplace(g).var()        # netlik: dusukse bulanik

def scored_patches(case):
    path=find(case)
    if not path: return -1,[]
    s=openslide.OpenSlide(path); pos=tissue_pos(s)
    if len(pos)<4: s.close(); return -1,[]
    random.seed(1); cand=random.sample(pos,min(20,len(pos)))
    sc=[]
    for p in cand:
        im=np.array(s.read_region(p,0,(P,P)).convert("RGB"))
        g=im.mean(2); tissue=(g<210).mean()
        if tissue<0.4: continue           # cok bos patch'i at
        sharp=sharpness(im)
        sc.append((sharp, tissue, im))
    if len(sc)<2: s.close(); return -1,[]
    sc.sort(key=lambda t:-t[0])           # en keskin patch'ler once
    s.close()
    quality=np.median([x[0] for x in sc[:4]])   # hasta netlik skoru
    return quality,[im for _,_,im in sc[:2]]

rows=list(csv.DictReader(open(r"E:\KIRC\gene_morph_updown.csv")))
seen=[]; [seen.append(r["satir"]) for r in rows if r["satir"] not in seen]
gruplar={x:[r for r in rows if r["satir"]==x] for x in seen}
etik={x:gruplar[x][0]["etiket"] for x in seen}

chosen={}
for x in seen:
    scored=[]
    for r in gruplar[x]:
        q,pats=scored_patches(r["submitter_id"])
        if pats: scored.append((q,r["submitter_id"],pats))
    scored.sort(key=lambda t:-t[0])       # en keskin hastalar
    chosen[x]=scored[:2]
    print(x,"->",[(c[1],round(c[0],1)) for c in chosen[x]])

fig=plt.figure(figsize=(12,13))
gs=gridspec.GridSpec(4,4,hspace=0.34,wspace=0.08)
colt=["Patch 1","Patch 2","Patch 1","Patch 2"]
for ri,x in enumerate(seen):
    col=0
    for (q,case,pats) in chosen[x]:
        for pj,im in enumerate(pats):
            a=fig.add_subplot(gs[ri,col]); a.imshow(im); a.set_xticks([]); a.set_yticks([])
            if col==0: a.set_ylabel(etik[x],fontsize=10)
            if ri==0: a.text(0.5,1.14,colt[col],transform=a.transAxes,fontsize=9,ha="center",color="#555")
            if pj==0: a.text(0.0,1.03,case,transform=a.transAxes,fontsize=8.5,ha="left",color="#222")
            col+=1
fig.suptitle("Morphology by up/down signature-gene expression (best image-predicted genes)",
             fontsize=14,weight="bold",y=0.94)
fig.text(0.5,0.06,"2 patients x 2 patches per group   |   DLX4: up-regulated/risky   ·   AJAP1: down-regulated/protective",
         fontsize=10,ha="center",color="#5B7682",style="italic")
plt.savefig(r"E:\KIRC\Fig_updown_morphology.png",dpi=170,bbox_inches="tight")
print("SAVED")
