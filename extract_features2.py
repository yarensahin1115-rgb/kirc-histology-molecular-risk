import os, glob, time, csv
import numpy as np
import multiprocessing as mp

PATCHES, BATCH, P = 200, 32, 224
OUT = r"E:\KIRC\slide_features.csv"
TIMEOUT = 1800   # bir slayt 30 dk'yi gecerse logla ve atla (sadece gercek sonsuz takilmayi yakalar)

def worker_loop(in_q, out_q):
    import torch, openslide, random
    from transformers import AutoImageProcessor, ViTModel
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    proc = AutoImageProcessor.from_pretrained("owkin/phikon")
    mdl = ViTModel.from_pretrained("owkin/phikon", add_pooling_layer=False).to(dev).eval()
    def embed(imgs):
        inp = proc(images=imgs, return_tensors="pt")["pixel_values"].to(dev)
        with torch.inference_mode():
            return mdl(pixel_values=inp).last_hidden_state[:, 0, :].cpu().numpy()
    out_q.put("ready")
    while True:
        path = in_q.get()
        if path is None: break
        try:
            slide = openslide.OpenSlide(path)
            W, H = slide.dimensions
            ds = 32
            thumb = slide.get_thumbnail((max(1, W//ds), max(1, H//ds))).convert("L")
            tissue = np.array(thumb) < 210
            sx, sy = thumb.size[0]/W, thumb.size[1]/H
            pos = []
            for y in range(0, H-P, P):
                ty0, ty1 = int(y*sy), int((y+P)*sy)
                for x in range(0, W-P, P):
                    tx0, tx1 = int(x*sx), int((x+P)*sx)
                    b = tissue[ty0:max(ty1,ty0+1), tx0:max(tx1,tx0+1)]
                    if b.size and b.mean() > 0.5: pos.append((x, y))
            if not pos:
                slide.close(); out_q.put(None); continue
            random.seed(0)
            sample = random.sample(pos, min(PATCHES, len(pos)))
            feats, batch = [], []
            for (x, y) in sample:
                batch.append(slide.read_region((x, y), 0, (P, P)).convert("RGB"))
                if len(batch) == BATCH:
                    feats.append(embed(batch)); batch = []
            if batch: feats.append(embed(batch))
            slide.close()
            out_q.put(np.concatenate(feats, axis=0).mean(axis=0))
        except Exception as e:
            out_q.put(("err", str(e)[:80]))

def main():
    jobs = []
    for grp, folder in [("Tumor", r"E:\KIRC\slides_tumor"), ("Normal", r"E:\KIRC\slides_normal")]:
        for f in glob.glob(os.path.join(folder, "**", "*.svs"), recursive=True):
            jobs.append((grp, f))
    print("Toplam slayt:", len(jobs))

    done = set()
    if os.path.exists(OUT):
        with open(OUT) as fh:
            for row in csv.reader(fh):
                if row: done.add(row[0])
    print("Zaten biten:", len(done))

    new = not os.path.exists(OUT)
    fh = open(OUT, "a", newline="")
    w = csv.writer(fh)
    if new:
        w.writerow(["slide_id", "case_barcode", "group"] + [f"f{i}" for i in range(768)])

    ctx = mp.get_context("spawn")
    in_q, out_q = ctx.Queue(), ctx.Queue()
    p = ctx.Process(target=worker_loop, args=(in_q, out_q), daemon=True); p.start()
    out_q.get(); print("Worker hazir.")

    n = 0
    for grp, path in jobs:
        sid = os.path.basename(os.path.dirname(path))
        if sid in done: continue
        case = "-".join(os.path.basename(path).split("-")[:3])
        n += 1; ts = time.time()
        in_q.put(path)
        try:
            res = out_q.get(timeout=TIMEOUT)
        except Exception:
            print(f"  DONDU, atlandi (>{TIMEOUT}s): {case}")
            p.terminate(); p.join()
            p = ctx.Process(target=worker_loop, args=(in_q, out_q), daemon=True); p.start()
            out_q.get(); continue
        if res is None:
            print(f"  doku yok: {case}"); continue
        if isinstance(res, tuple):
            print(f"  hata: {case} -> {res[1]}"); continue
        w.writerow([sid, case, grp] + [round(float(x), 5) for x in res]); fh.flush()
        print(f"{n}: {case} ({time.time()-ts:.1f}s)")
    in_q.put(None); fh.close()
    print("Bitti. Cikti:", OUT)

if __name__ == "__main__":
    main()
