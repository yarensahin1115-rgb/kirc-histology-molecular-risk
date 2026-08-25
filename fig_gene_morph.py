import csv, glob, os, random
import numpy as np, openslide
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = 224
cases = list(csv.DictReader(open(r"E:\KIRC\gene_morphology_cases.csv")))
gen = cases[0]["gen"]
allsvs = glob.glob(r"E:\KIRC\slides_tumor\**\*.svs", recursive=True)

def patches_of(case, k=3):
    hit = [f for f in allsvs if os.path.basename(f).startswith(case)]
    if not hit: return None
    s = openslide.OpenSlide(hit[0]); W,H = s.dimensions
    ds=32; gray=np.array(s.get_thumbnail((W//ds,H//ds)).convert("L")); tis=gray<210
    sx,sy=gray.shape[1]/W, gray.shape[0]/H
    pos=[]
    for y in range(0,H-P,P):
        for x in range(0,W-P,P):
            if tis[min(int(y*sy),gray.shape[0]-1),min(int(x*sx),gray.shape[1]-1)]: pos.append((x,y))
    random.seed(0)
    pk = random.sample(pos, min(k,len(pos)))
    out=[s.read_region(p,0,(P,P)).convert("RGB") for p in pk]
    s.close(); return out

low = [c["submitter_id"] for c in cases if c["grup"]=="low"]
high= [c["submitter_id"] for c in cases if c["grup"]=="high"]

# 2 satir (low/high) x 6 sutun (2 hasta x 3 patch) -> sade tutmak icin 3 hasta x 1 temsili patch yerine
# 2 satir x 3 hasta, her hasta 1 temsili patch
fig, ax = plt.subplots(2, 3, figsize=(11, 7.5))
for col, case in enumerate(low):
    pp = patches_of(case, 1)
    ax[0,col].imshow(pp[0]); ax[0,col].axis("off")
    ax[0,col].set_title(case, fontsize=10)
for col, case in enumerate(high):
    pp = patches_of(case, 1)
    ax[1,col].imshow(pp[0]); ax[1,col].axis("off")
    ax[1,col].set_title(case, fontsize=10)
ax[0,0].set_ylabel(f"Low {gen}", fontsize=13); ax[0,0].axis("on")
ax[0,0].set_xticks([]); ax[0,0].set_yticks([])
ax[1,0].set_ylabel(f"High {gen}", fontsize=13); ax[1,0].axis("on")
ax[1,0].set_xticks([]); ax[1,0].set_yticks([])
plt.suptitle(f"Morphology vs {gen} expression (best image-predicted signature gene, r=0.62)", fontsize=13)
plt.tight_layout()
plt.savefig(r"E:\KIRC\Fig_gene_morphology.png", dpi=180, bbox_inches="tight")
print("SAVED: E:\\KIRC\\Fig_gene_morphology.png")
