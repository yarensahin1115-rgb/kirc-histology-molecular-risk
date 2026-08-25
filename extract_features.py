import os, glob, random, time, csv
import numpy as np
import torch
import openslide
from transformers import AutoImageProcessor, ViTModel

DEVICE  = "cuda" if torch.cuda.is_available() else "cpu"
PATCHES = 200          # slayt basina ornek patch
BATCH   = 32           # OOM olursa 16 yap
P       = 224
OUT     = r"E:\KIRC\slide_features.csv"
print("Cihaz:", DEVICE)

processor = AutoImageProcessor.from_pretrained("owkin/phikon")
model = ViTModel.from_pretrained("owkin/phikon", add_pooling_layer=False).to(DEVICE).eval()

def tissue_positions(slide):
    W, H = slide.dimensions
    ds = 32
    thumb = slide.get_thumbnail((max(1, W // ds), max(1, H // ds))).convert("L")
    tissue = np.array(thumb) < 210
    sx, sy = thumb.size[0] / W, thumb.size[1] / H
    pos = []
    for y in range(0, H - P, P):
        ty0, ty1 = int(y * sy), int((y + P) * sy)
        for x in range(0, W - P, P):
            tx0, tx1 = int(x * sx), int((x + P) * sx)
            b = tissue[ty0:max(ty1, ty0 + 1), tx0:max(tx1, tx0 + 1)]
            if b.size and b.mean() > 0.5:
                pos.append((x, y))
    return pos

def embed(imgs):
    inp = processor(images=imgs, return_tensors="pt")["pixel_values"].to(DEVICE)
    with torch.inference_mode():
        out = model(pixel_values=inp).last_hidden_state[:, 0, :]
    return out.cpu().numpy()

def slide_vector(path):
    slide = openslide.OpenSlide(path)
    pos = tissue_positions(slide)
    if not pos:
        slide.close(); return None
    random.seed(0)
    sample = random.sample(pos, min(PATCHES, len(pos)))
    feats, batch = [], []
    for (x, y) in sample:
        batch.append(slide.read_region((x, y), 0, (P, P)).convert("RGB"))
        if len(batch) == BATCH:
            feats.append(embed(batch)); batch = []
    if batch:
        feats.append(embed(batch))
    slide.close()
    return np.concatenate(feats, axis=0).mean(axis=0)

# Islenecek tum slaytlar (tumor + normal)
jobs = []
for grp, folder in [("Tumor", r"E:\KIRC\slides_tumor"), ("Normal", r"E:\KIRC\slides_normal")]:
    for f in glob.glob(os.path.join(folder, "**", "*.svs"), recursive=True):
        jobs.append((grp, f))
print("Toplam slayt:", len(jobs))

# Devam: bitmis slaytlari atla
done = set()
if os.path.exists(OUT):
    with open(OUT) as fh:
        for row in csv.reader(fh):
            if row: done.add(row[0])

new = not os.path.exists(OUT)
fh = open(OUT, "a", newline="")
w = csv.writer(fh)
if new:
    w.writerow(["slide_id", "case_barcode", "group"] + [f"f{i}" for i in range(768)])

t0 = time.time()
for i, (grp, path) in enumerate(jobs, 1):
    sid = os.path.basename(os.path.dirname(path))      # UUID klasoru = benzersiz
    if sid in done:
        continue
    case = "-".join(os.path.basename(path).split("-")[:3])
    try:
        v = slide_vector(path)
        if v is None:
            print(i, "DOKU YOK:", case); continue
        w.writerow([sid, case, grp] + [round(float(x), 5) for x in v])
        fh.flush()
    except Exception as e:
        print(i, "HATA:", case, str(e)[:80]); continue
    if i % 20 == 0:
        dt = time.time() - t0
        print(f"{i}/{len(jobs)}  ({dt/60:.1f} dk, ~{dt/max(i,1):.1f} sn/slayt)")
fh.close()
print("Bitti. Cikti:", OUT)
