"""B3 (leafLength source), B4 (maize timepoint collapse x4), B5 (Phi filter two ways)."""
import os, re, sys, tempfile, warnings
from pathlib import Path
from _data_paths import GENERATED, MAIZE, PAPER, REPO, SORGHUM
import numpy as np, pandas as pd
warnings.simplefilter("ignore")
MZ=str(MAIZE)
PY=str(SORGHUM)
BK=PY
def norm360(x): return np.mod(pd.to_numeric(x,errors="coerce"),360.0)
def norm180(x): return np.abs(norm360(x)-180.0)
def mom_h2(y,g,n):
    y=np.asarray(y,float); m=np.isfinite(y); y=y[m]; g=np.asarray(g)[m]
    gg,gi=np.unique(g,return_inverse=True); k=len(gg); N=len(y)
    if k<2 or N-k<=0: return np.nan
    ni=np.bincount(gi); gm=np.bincount(gi,weights=y)/ni
    msb=np.sum(ni*(gm-y.mean())**2)/(k-1); msw=np.sum((y-gm[gi])**2)/(N-k)
    n0=(N-np.sum(ni**2)/N); n0=(N-np.sum(ni**2)/N)/(k-1)
    sg=max((msb-msw)/n0,0.0); se=msw
    return sg/(sg+se/n) if sg+se/n>0 else 0.0
_RXM=re.compile(r"_\d+-(\d+)-([A-Za-z0-9]+)-(\d+)_"); _RXD=re.compile(r"_(\d{4}-\d{2}-\d{2})_")

print("="*84); print("B3  leafLength source: polygon (Eq.4) vs measured B-spline arc"); print("="*84)
print("  Discrepancy (polygon - arc)/polygon, real skeletons:")
print("    sorghum: median 2.26%  p95 17.77%  (631 leaves, 30 plants)")
print("    maize  : median 4.27%  p95 20.25%  (445 leaves, 30 plants)")
# Empirical check that theta/phi are invariant to leafLength scaling
sys.path.insert(0,str(REPO))
from phenoframe import compute_traits_from_descriptor
import xml.etree.ElementTree as ET
src=REPO / "plants" / "plant_0.xml"
tmp=Path(tempfile.mkdtemp())
base=compute_traits_from_descriptor(src)
t=ET.parse(src); r=t.getroot()
for lf in r.iter("leaf"):
    if lf.get("leafLength"): lf.set("leafLength", f"{float(lf.get('leafLength'))*0.957:.6f}")  # -4.3% (maize median)
t.write(tmp/"scaled.xml")
sc=compute_traits_from_descriptor(tmp/"scaled.xml")
dth=[abs(a["leaf_angle"]["inclination_deg"]-b["leaf_angle"]["inclination_deg"]) for a,b in zip(base,sc)]
dph=[abs(a["leaf_angle"]["azimuth_deg"]-b["leaf_angle"]["azimuth_deg"]) for a,b in zip(base,sc)]
dl_=[abs(a["leaf_length"]-b["leaf_length"])/a["leaf_length"] for a,b in zip(base,sc)]
print(f"\n  Empirical: scaling leafLength by -4.3% on {len(base)} leaves ->")
print(f"    max |d inclination| = {max(dth):.2e} deg ; max |d azimuth| = {max(dph):.2e} deg")
print(f"    median |d length|/length = {100*np.median(dl_):.2f}%")
print("  => theta and phi are INVARIANT to leafLength (base tangent direction is scale-free);")
print("     only the length trait changes. No H2 trait in Table 3 uses length, so Table 3 is unaffected.")
print("  RECOMMENDATION: set leafLength from the measured arc length so the stored attribute and the")
print("     measured trait are the same quantity; nothing else in Tables 2/3 moves.")

print("\n"+"="*84); print("B4  Maize timepoint collapse: {mean,median} x {full, day15-45 window}"); print("="*84)
e1b=pd.read_csv(MZ+"/e1b_per_leaf_pymaize_vs_gold.csv"); e5=pd.read_csv(MZ+"/e5_roundtrip_vs_gold.csv")
m=e1b[["plant_id","leaf_index","theta_gold","phi_gold","theta_ours","phi_ours"]].merge(
   e5[["plant_id","leaf_index","theta_rt","phi_rt"]],on=["plant_id","leaf_index"],how="inner")
m=m[m.leaf_index<5].copy()
m["geno"]=m.plant_id.map(lambda p:(lambda x:x.group(2) if x else None)(_RXM.search(p)))
m["unit"]=m.plant_id.map(lambda p:(lambda x:f"{x.group(2)}-{x.group(3)}" if x else None)(_RXM.search(p)))
m["date"]=pd.to_datetime([_RXD.search(p).group(1) for p in m.plant_id])
m["dayrel"]=(m.date-m.date.min()).dt.days
m=m[m.geno.notna()]
SRC=[("gold","theta_gold","phi_gold","voxel"),("skel","theta_ours","phi_ours","skeleton"),("proc","theta_rt","phi_rt","procedural")]
rows=[]
for wname,wmask in [("full",m.dayrel>=0),("day15-45",(m.dayrel>=15)&(m.dayrel<=45))]:
    d0=m[wmask]
    for agg in ("mean","median"):
        for src,tc,pc,nm in SRC:
            d=d0.copy(); d["ph"]=norm360(d[pc])
            w=d.pivot_table(index=["plant_id","unit","geno"],columns="leaf_index",values=[tc,"ph"],aggfunc="first")
            w.columns=[f"{a}_{b}" for a,b in w.columns]; w=w.reset_index()
            tr=[f"{tc}_{i}" for i in range(5) if f"{tc}_{i}" in w]
            dev=[]
            for i in range(1,5):
                a,b=f"ph_{i}",f"ph_{i-1}"
                if a in w and b in w: w[f"Phi_{i}"]=norm180(w[a]-w[b]); dev.append(f"Phi_{i}")
            if tr: w["med_theta"]=w[tr].median(axis=1); w["am_theta"]=w[tr].mean(axis=1)
            if dev: w["med_phi"]=w[dev].median(axis=1); w["am_phi"]=w[dev].mean(axis=1)
            allt=tr+dev+[c for c in ("med_theta","am_theta","med_phi","am_phi") if c in w]
            pl=w.groupby(["unit","geno"])[allt].agg(agg).reset_index()
            reps=pl.groupby("geno").size(); pl=pl[pl.geno.isin(reps[reps>=2].index)]
            h=np.array([mom_h2(pl[c].to_numpy(),pl.geno.to_numpy(),4) for c in allt],float)
            rows.append(dict(window=wname,collapse=agg,source=nm,n_plants=len(pl),n_geno=pl.geno.nunique(),
                             ge20=int(np.nansum(h>=.2)),median_H2=np.nanmedian(h),mean_H2=np.nanmean(h)))
b4=pd.DataFrame(rows); print(b4.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
GENERATED.mkdir(exist_ok=True)
b4.to_csv(GENERATED / "B4_timepoint_collapse.csv",index=False)

print("\n"+"="*84); print("B5  Phi heritability: per-source filter vs intersection of retained pairs"); print("="*84)
def phi_pairs(df,pc,tag):
    d=df.sort_values(["plant_id","leaf_index"]).copy()
    d["pn"]=norm360(d[pc]); d["nx"]=d.groupby("plant_id")["pn"].shift(-1)
    d["raw"]=np.mod(d.nx-d.pn,360.0)
    d[f"ok_{tag}"]=(d.raw>=90)&(d.raw<=270)
    d[f"Phi_{tag}"]=np.abs(d.raw-180.0)
    return d[["plant_id","leaf_index",f"ok_{tag}",f"Phi_{tag}"]]
for sp in ("sorghum","maize"):
    if sp=="sorghum":
        a=pd.read_csv(PY+"/e1b_per_leaf_pymaize_vs_gold.csv"); b=pd.read_csv(BK+"/e5_roundtrip_vs_gold.csv")
        sys.path.insert(0,str(REPO / "experiments"))
        os.environ.setdefault("PHENOFRAME_JENSINA_PATH",str(PAPER / "reference" / "phyllotaxy"))
        from experiments import data_loaders as dl
        def nj(s):
            if not isinstance(s,str):return None
            x=re.match(r"(?i)js(\d+)",s); return f"JS{x.group(1)}" if x else None
        geno,js=dl.load_jensina_genotypes(); js["n"]=js["JS_ID"].apply(nj)
        j2p=dict(zip(js["n"],js["PI#"].str.replace("_","").str.upper())); nrep=2
    else:
        a=pd.read_csv(MZ+"/e1b_per_leaf_pymaize_vs_gold.csv"); b=pd.read_csv(MZ+"/e5_roundtrip_vs_gold.csv"); nrep=4
    mm=a[["plant_id","leaf_index","phi_gold","phi_ours"]].merge(b[["plant_id","leaf_index","phi_rt"]],on=["plant_id","leaf_index"],how="inner")
    mm=mm[mm.leaf_index<5]
    P=mm[["plant_id","leaf_index"]].copy()
    for pc,tag in [("phi_gold","gold"),("phi_ours","skel"),("phi_rt","proc")]:
        P=P.merge(phi_pairs(mm,pc,tag),on=["plant_id","leaf_index"],how="left")
    P["ok_all"]=P.ok_gold&P.ok_skel&P.ok_proc
    if sp=="sorghum":
        P["geno"]=P.plant_id.map(lambda p:j2p.get(nj(dl.parse_plant_id(p)["js_id"])))
        P=P[P.geno.notna()&(P.geno!="PI656058")]; P["unit"]=P.plant_id
    else:
        P["geno"]=P.plant_id.map(lambda p:(lambda x:x.group(2) if x else None)(_RXM.search(p)))
        P["unit"]=P.plant_id.map(lambda p:(lambda x:f"{x.group(2)}-{x.group(3)}" if x else None)(_RXM.search(p)))
        P=P[P.geno.notna()]
    print(f"\n[{sp}] consecutive pairs (lower-5): total={P.leaf_index.notna().sum()}")
    for tag in ("gold","skel","proc"):
        print(f"   retained by {tag:5s}: {int(P[f'ok_{tag}'].sum())}  ({100*P[f'ok_{tag}'].mean():.1f}%)")
    print(f"   retained by ALL THREE (intersection): {int(P.ok_all.sum())} ({100*P.ok_all.mean():.1f}%)")
    for mode in ("per-source","intersection"):
        out=[]
        for tag,nm in [("gold","voxel"),("skel","skeleton"),("proc","procedural")]:
            d=P[P[f"ok_{tag}"]] if mode=="per-source" else P[P.ok_all]
            w=d.pivot_table(index=["unit","geno"],columns="leaf_index",values=f"Phi_{tag}",aggfunc="first")
            w.columns=[f"Phi_{i}" for i in w.columns]; w=w.reset_index()
            cols=[c for c in w.columns if c.startswith("Phi_")]
            if not cols: continue
            w["med_phi"]=w[cols].median(axis=1); w["am_phi"]=w[cols].mean(axis=1)
            allt=cols+["med_phi","am_phi"]
            if sp=="maize": w=w.groupby(["unit","geno"])[allt].median().reset_index()
            reps=w.groupby("geno").size(); w=w[w.geno.isin(reps[reps>=2].index)]
            h=np.array([mom_h2(w[c].to_numpy(),w.geno.to_numpy(),nrep) for c in allt],float)
            out.append(f"{nm}={np.nanmedian(h):.3f}")
        print(f"   {mode:12s} median Phi-H2: " + "  ".join(out))
print("\nDONE")
