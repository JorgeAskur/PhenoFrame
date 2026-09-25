"""R1-R3 definitive heritability. Correct order of operations:
   derive per-SCAN traits -> apply filters -> collapse scans (maize) -> REML H2.
Conventions: leaf set both-keep; Phi on common-pair intersection; maize collapse = median,
full window; sorghum n=2, maize n=4; estimator REML (statsmodels MixedLM)."""
import os, re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
warnings.simplefilter("ignore")
PAPER = Path(__file__).resolve().parents[1]
MZ = str(PAPER / "source_data" / "maize")
PY = str(PAPER / "source_data" / "sorghum")
BK = PY
_RXM=re.compile(r"_\d+-(\d+)-([A-Za-z0-9]+)-(\d+)_")
def n360(x): return np.mod(pd.to_numeric(x,errors="coerce"),360.0)
def n180(x): return np.abs(x-180.0)

def reml_h2(df,trait,n):
    s=df[["geno",trait]].dropna().rename(columns={trait:"y"})
    c=s.groupby("geno").size(); s=s[s.geno.isin(c[c>1].index)]
    if s.geno.nunique()<2: return np.nan,0,0
    sg=se=None
    for meth in ("bfgs","powell","lbfgs"):
        try:
            f=smf.mixedlm("y ~ 1",s,groups=s["geno"]).fit(method=meth,reml=True)
            if f.converged: sg=float(f.cov_re.iloc[0,0]); se=float(f.scale); break
        except Exception: continue
    if sg is None: return np.nan,len(s),s.geno.nunique()
    d=sg+se/n
    return (sg/d if d>0 else 0.0), len(s), s.geno.nunique()

def sorghum_geno(pids):
    os.environ.setdefault("PHENOFRAME_JENSINA_PATH", str(PAPER / "reference" / "phyllotaxy"))
    sys.path.insert(0, str(PAPER.parent / "experiments"))
    from experiments import data_loaders as dl
    def nj(s):
        if not isinstance(s,str): return None
        m=re.match(r"(?i)js(\d+)",s); return f"JS{m.group(1)}" if m else None
    _,js=dl.load_jensina_genotypes(); js["n"]=js["JS_ID"].apply(nj)
    j2p=dict(zip(js["n"],js["PI#"].str.replace("_","").str.upper()))
    return {p:j2p.get(nj(dl.parse_plant_id(p)["js_id"])) for p in pids}

def build(species):
    if species=="sorghum":
        a=pd.read_csv(PY+"/e1b_per_leaf_pymaize_vs_gold.csv"); b=pd.read_csv(BK+"/e5_roundtrip_vs_gold.csv")
    else:
        a=pd.read_csv(MZ+"/e1b_per_leaf_pymaize_vs_gold.csv"); b=pd.read_csv(MZ+"/e5_roundtrip_vs_gold.csv")
    m=a[["plant_id","leaf_index","theta_gold","phi_gold","theta_ours","phi_ours","keep"]].merge(
        b[["plant_id","leaf_index","theta_rt","phi_rt","keep"]].rename(columns={"keep":"keep_p"}),
        on=["plant_id","leaf_index"],how="inner")
    m=m[m.leaf_index<5].copy()
    if species=="sorghum":
        gm=sorghum_geno(m.plant_id.unique()); m["geno"]=m.plant_id.map(gm)
        m=m[m.geno.notna()&(m.geno!="PI656058")]
        counts=m[["plant_id","geno"]].drop_duplicates().groupby("geno").size()
        m=m[m.geno.isin(counts[counts>=2].index)]
        m["unit"]=m.plant_id
    else:
        m["geno"]=m.plant_id.map(lambda p:(lambda x:x.group(2) if x else None)(_RXM.search(p)))
        m["unit"]=m.plant_id.map(lambda p:(lambda x:f"{x.group(2)}-{x.group(3)}" if x else None)(_RXM.search(p)))
        m=m[m.geno.notna()]
    return m

SRC=[("gold","theta_gold","phi_gold","voxel"),("skel","theta_ours","phi_ours","skeleton"),
     ("proc","theta_rt","phi_rt","procedural")]

def leafmask(m,basis,src):
    if basis=="own":  return (m.keep==True) if src in("gold","skel") else (m.keep_p==True)
    if basis=="both": return (m.keep==True)&(m.keep_p==True)
    return pd.Series(True,index=m.index)

def per_scan_traits(m, basis, src, tc, pc, phi_mode):
    """Derive traits PER SCAN (correct order), then collapse."""
    d=m[leafmask(m,basis,src)].copy()
    d["ph"]=n360(d[pc])
    # wide per scan
    w=d.pivot_table(index=["plant_id","unit","geno"],columns="leaf_index",values=[tc,"ph"],aggfunc="first")
    w.columns=[f"{a}_{b}" for a,b in w.columns]; w=w.reset_index()
    th=[f"{tc}_{i}" for i in range(5) if f"{tc}_{i}" in w]
    # raw divergence per scan for EVERY source (needed for the intersection mask)
    raws={}
    for s2,_,p2,_ in SRC:
        dd=m.copy(); dd["ph2"]=n360(dd[p2])
        ww=dd.pivot_table(index=["plant_id"],columns="leaf_index",values="ph2",aggfunc="first")
        ww.columns=[f"p_{c}" for c in ww.columns]
        for i in range(1,5):
            a,b=f"p_{i}",f"p_{i-1}"
            raws[(s2,i)]=np.mod(ww[a]-ww[b],360.0) if (a in ww and b in ww) else None
        raws[(s2,"idx")]=ww.index
    dev=[]
    for i in range(1,5):
        a,b=f"ph_{i}",f"ph_{i-1}"
        if a not in w or b not in w: continue
        raw=np.mod(w[a]-w[b],360.0)
        ok=(raw>=90)&(raw<=270)
        if phi_mode=="intersection":
            idx=raws[(SRC[0][0],"idx")]
            okall=pd.Series(True,index=w.plant_id)
            for s2,_,_,_ in SRC:
                r2=raws[(s2,i)]
                if r2 is None: continue
                v=r2.reindex(w.plant_id).to_numpy()
                okall &= pd.Series((v>=90)&(v<=270),index=w.plant_id)
            ok=ok.to_numpy()&okall.to_numpy()
        c=f"Phi_{i}"; w[c]=np.where(ok,n180(raw),np.nan); dev.append(c)
    if th: w["med_theta"]=w[th].median(axis=1); w["am_theta"]=w[th].mean(axis=1)
    if dev: w["med_phi"]=w[dev].median(axis=1); w["am_phi"]=w[dev].mean(axis=1)
    allt=th+dev+[c for c in ("med_theta","am_theta","med_phi","am_phi") if c in w]
    if species_is_maize: w=w.groupby(["unit","geno"])[allt].median().reset_index()
    else: w=w[["unit","geno"]+allt]
    return w, th, dev, allt

print("="*96)
print("R1  Maize trait-count conflict — diagnosis")
print("="*96)
print("  Round-3 B1 pivoted with index=[unit,geno], aggfunc='mean': scans were collapsed AT THE PIVOT")
print("  by arithmetic MEAN, and phi (a circular variable) was averaged linearly BEFORE Phi was derived.")
print("  -> B1's label 'median' was wrong, and its Phi traits were corrupted (mean of 359 and 1 = 180).")
print("  -> B4 (derive per-scan, then collapse) is correct. B1's maize Table-3 rows are SUPERSEDED.")
print()

results={}
for species in ("sorghum","maize"):
    species_is_maize = (species=="maize")
    m=build(species); nrep = 4 if species_is_maize else 2
    for basis in ("own","both","all"):
        for phi_mode in ("per-source","intersection"):
            for src,tc,pc,nm in SRC:
                w,th,dev,allt=per_scan_traits(m,basis,src,tc,pc,phi_mode)
                hs={}
                for c in allt:
                    h,npl,ng=reml_h2(w,c,nrep); hs[c]=h
                arr=np.array([hs[c] for c in allt],float)
                results[(species,basis,phi_mode,nm)]=dict(
                    h=hs,allt=allt,th=th,dev=dev,arr=arr,n_plants=len(w),n_geno=w.geno.nunique())

def row(species,basis,phi_mode,nm):
    r=results[(species,basis,phi_mode,nm)]; a=r["arr"]
    return f"{np.nansum(a>=0.20):>2.0f}/{len(a)}   {np.nanmedian(a):.3f}   {np.nanmean(a):.3f}"

print("="*96)
print("R1/R2  Table 3 in REML — all three leaf bases (Phi per-source, to compare with published)")
print("="*96)
print(f"{'species':8} {'basis':6} {'source':11} {'>=.20':>7} {'median':>8} {'mean':>8}   n_plants/n_geno")
for species in ("sorghum","maize"):
    for basis,lab in (("own","own"),("both","both"),("all","all")):
        for _,_,_,nm in SRC:
            r=results[(species,basis,"per-source",nm)]
            print(f"{species:8} {lab:6} {nm:11} {row(species,basis,'per-source',nm)}   {r['n_plants']}/{r['n_geno']}")
print("\n  Published sorghum Table 3 = 0.32 / 0.29 / 0.48 (median H2, V/S/P) -> compare rows above.")
print("  Published maize   Table 3 = 0.77 / 0.82 / 0.79 at n=2 (paper); here maize is at n=4.")

print("\n"+"="*96)
print("R2  DEFINITIVE Table 3 — both-keep, Phi on common-pair intersection, REML")
print("="*96)
print(f"{'species':8} {'source':11} {'>=.20':>7} {'median':>8} {'mean':>8}   n_plants/n_geno")
for species in ("sorghum","maize"):
    for _,_,_,nm in SRC:
        r=results[(species,"both","intersection",nm)]
        print(f"{species:8} {nm:11} {row(species,'both','intersection',nm)}   {r['n_plants']}/{r['n_geno']}")

print("\n"+"="*96)
print("R3  Per-trait Phi heritabilities (REML), per-source vs intersection, both-keep")
print("="*96)
for species in ("sorghum","maize"):
    print(f"\n[{species}]  trait      per-source (V/S/P)        intersection (V/S/P)")
    dev=results[(species,"both","per-source","voxel")]["dev"]
    for c in dev+["med_phi","am_phi"]:
        a=[results[(species,"both","per-source",nm)]["h"].get(c,np.nan) for _,_,_,nm in SRC]
        b=[results[(species,"both","intersection",nm)]["h"].get(c,np.nan) for _,_,_,nm in SRC]
        print(f"   {c:10}  {a[0]:.3f} / {a[1]:.3f} / {a[2]:.3f}      {b[0]:.3f} / {b[1]:.3f} / {b[2]:.3f}")
    # claim checks
    ag=[results[(species,"both","intersection",nm)]["h"] for _,_,_,nm in SRC]
    mp=[g.get("med_phi",np.nan) for g in ag]; ap=[g.get("am_phi",np.nan) for g in ag]
    print(f"   CLAIM aggregate Phi H2 >= 0.45 in all three (intersection): "
          f"median-Phi min={np.nanmin(mp):.3f}  mean-Phi min={np.nanmin(ap):.3f} "
          f"-> {'HOLDS' if min(np.nanmin(mp),np.nanmin(ap))>=0.45 else 'FAILS'}")

# trait names differ per source (theta_gold_0 / theta_ours_0 / theta_rt_0) -> compare POSITIONALLY
def arr_of(species,basis,mode,nm):
    r=results[(species,basis,mode,nm)]; return np.array([r["h"][c] for c in r["allt"]],float), r["allt"]
LBL=["theta_0","theta_1","theta_2","theta_3","theta_4","Phi_1","Phi_2","Phi_3","Phi_4",
     "med_theta","am_theta","med_phi","am_phi"]
av,_=arr_of("maize","both","intersection","voxel")
as_,_=arr_of("maize","both","intersection","skeleton")
ap,_=arr_of("maize","both","intersection","procedural")
d=np.abs(av-ap); d=d[np.isfinite(d)]
print(f"\n  CLAIM maize max |dH2(P<->V)| across 13 traits (both-keep, intersection, REML, n=4): {d.max():.3f}")
print("\n  DEFINITIVE maize per-trait H2 (both-keep, intersection, REML, n=4):")
print(f"   {'trait':12} {'voxel':>7} {'skeleton':>9} {'procedural':>11}")
for i,lab in enumerate(LBL[:len(av)]):
    print(f"   {lab:12} {av[i]:>7.3f} {as_[i]:>9.3f} {ap[i]:>11.3f}")

# ---- R1: compare against the PUBLISHED e8 file to identify its basis ----
print("\n"+"="*96); print("R1  Which basis reproduces the PUBLISHED sorghum Table 3?"); print("="*96)
pub=pd.read_csv(PY+"/e8_heritability.csv")
pm={"gold":"voxel","e1b":"skeleton","e5":"procedural"}
pub["src"]=pub.source.map(pm)
print("  published (e8_heritability.csv, REML, n=2): median H2 per source")
for nm in ("voxel","skeleton","procedural"):
    s=pub[pub.src==nm]
    print(f"    {nm:11} median={s.H2_n2.median():.3f}  mean={s.H2_n2.mean():.3f}  >=.20={int((s.H2_n2>=.2).sum())}/13"
          f"  n_plants={s.n_plants.iloc[0]}  n_geno={s.n_genotypes.iloc[0]}")
print("  my REML reproduction (per-source Phi), by basis:")
for basis in ("own","both","all"):
    vals=[]
    for nm in ("voxel","skeleton","procedural"):
        a,_=arr_of("sorghum",basis,"per-source",nm); vals.append(np.nanmedian(a))
    print(f"    {basis:5}: {vals[0]:.3f} / {vals[1]:.3f} / {vals[2]:.3f}")

out=[]
for species in ("sorghum","maize"):
    for nm in ("voxel","skeleton","procedural"):
        a,cols=arr_of(species,"both","intersection",nm)
        for i,c in enumerate(cols):
            out.append(dict(species=species,source=nm,trait=LBL[i] if i<len(LBL) else c,H2=a[i]))
(PAPER / "generated").mkdir(exist_ok=True)
pd.DataFrame(out).to_csv(PAPER / "generated" / "R123_table3_definitive.csv",index=False)
print("\nsaved R123_table3_definitive.csv")
print("DONE")
