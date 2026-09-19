"""B1 — Tables 2 and 3 regenerated under three leaf-set definitions.

Leaf sets (all restricted to the analysis window):
  (a) own-retained : each source scored on ITS OWN retained leaves  [current published basis]
                     skeleton -> e1b.keep ; procedural -> e5.keep
  (b) both-keep    : e1b.keep AND e5.keep  (intersection of quality filters)
  (c) all-common   : every leaf measured by all three sources, NO quality filter
Table 2 window = lower-4 ; Table 3 (H2) window = lower-5.
Chamfer/Hausdorff are WHOLE-PLANT metrics (e4) and do not depend on the leaf set — reported once.
"""
import re, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.simplefilter("ignore")
PAPER = Path(__file__).resolve().parents[1]
MZ = str(PAPER / "source_data" / "maize")
PY = str(PAPER / "source_data" / "sorghum")
BK = PY
LOWER5=list(range(5))

def R2(a,b):
    a=pd.to_numeric(a,errors="coerce").to_numpy(float); b=pd.to_numeric(b,errors="coerce").to_numpy(float)
    m=np.isfinite(a)&np.isfinite(b)
    if m.sum()<3: return np.nan,0
    r=np.corrcoef(a[m],b[m])[0,1]; return r*r,int(m.sum())
def rcc(a,b):
    a=np.radians(pd.to_numeric(a,errors="coerce")); b=np.radians(pd.to_numeric(b,errors="coerce"))
    m=np.isfinite(a)&np.isfinite(b); a,b=a[m],b[m]
    if len(a)<3: return np.nan
    A=a-np.arctan2(np.sin(a).sum(),np.cos(a).sum()); B=b-np.arctan2(np.sin(b).sum(),np.cos(b).sum())
    return float(np.sum(np.sin(A)*np.sin(B))/np.sqrt(np.sum(np.sin(A)**2)*np.sum(np.sin(B)**2)))
def medabs(x): return float(np.nanmedian(np.abs(pd.to_numeric(x,errors="coerce"))))
def wrap(x): return ((pd.to_numeric(x,errors="coerce")+180)%360)-180
def norm360(x): return np.mod(pd.to_numeric(x,errors="coerce"),360.0)
def norm180(x): return np.abs(norm360(x)-180.0)
def mom_h2(y,g,n):
    y=np.asarray(y,float); m=np.isfinite(y); y=y[m]; g=np.asarray(g)[m]
    gg,gi=np.unique(g,return_inverse=True); k=len(gg); N=len(y)
    if k<2 or N-k<=0: return np.nan
    ni=np.bincount(gi); gm=np.bincount(gi,weights=y)/ni
    msb=np.sum(ni*(gm-y.mean())**2)/(k-1); msw=np.sum((y-gm[gi])**2)/(N-k)
    n0=(N-np.sum(ni**2)/N)/(k-1); sg=max((msb-msw)/n0,0.0); se=msw
    return sg/(sg+se/n) if sg+se/n>0 else 0.0

_RXM=re.compile(r"_\d+-(\d+)-([A-Za-z0-9]+)-(\d+)_")
def mz_geno(p):
    m=_RXM.search(p); return m.group(2) if m else None
def mz_plant(p):
    m=_RXM.search(p); return f"{m.group(2)}-{m.group(3)}" if m else None

def load(species):
    if species=="sorghum":
        e1b=pd.read_csv(PY+"/e1b_per_leaf_pymaize_vs_gold.csv"); e5=pd.read_csv(BK+"/e5_roundtrip_vs_gold.csv")
    else:
        e1b=pd.read_csv(MZ+"/e1b_per_leaf_pymaize_vs_gold.csv"); e5=pd.read_csv(MZ+"/e5_roundtrip_vs_gold.csv")
    m=e1b[["plant_id","leaf_index","theta_gold","phi_gold","theta_ours","phi_ours","keep"]].merge(
        e5[["plant_id","leaf_index","theta_rt","phi_rt","keep"]].rename(columns={"keep":"keep_p"}),
        on=["plant_id","leaf_index"], how="inner")
    return m

def sorghum_geno_map(plant_ids):
    import sys, os
    os.environ.setdefault("PHENOFRAME_JENSINA_PATH", str(PAPER / "reference" / "phyllotaxy"))
    sys.path.insert(0, str(PAPER.parent / "experiments"))
    from experiments import data_loaders as dl
    def nj(s):
        if not isinstance(s,str): return None
        mm=re.match(r"(?i)js(\d+)",s); return f"JS{mm.group(1)}" if mm else None
    geno,js_pi=dl.load_jensina_genotypes(); js_pi["JS_norm"]=js_pi["JS_ID"].apply(nj)
    j2p=dict(zip(js_pi["JS_norm"],js_pi["PI#"].str.replace("_","").str.upper()))
    return {p:j2p.get(nj(dl.parse_plant_id(p)["js_id"])) for p in plant_ids}

def mask_for(m, which, src):
    if which=="a":   # own-retained
        return m.keep==True if src in ("skel","gold") else m.keep_p==True
    if which=="b":   # both-keep
        return (m.keep==True)&(m.keep_p==True)
    return pd.Series(True,index=m.index)   # c: all-common

# ---------------- Table 2 ----------------
print("="*92)
print("B1 / TABLE 2  — R2 and median |error|, lower-4, under three leaf sets")
print("="*92)
t2rows=[]
for sp in ("sorghum","maize"):
    m=load(sp); m4=m[m.leaf_index<4]
    for which,label in [("a","(a) own-retained [published]"),("b","(b) both-keep"),("c","(c) all-common")]:
        for src,tc,pc in [("skel","theta_ours","phi_ours"),("proc","theta_rt","phi_rt")]:
            d=m4[mask_for(m4,which,src)]
            r2t,nt=R2(d.theta_gold,d[tc]); r2p,_=R2(d.phi_gold,d[pc])
            t2rows.append(dict(species=sp,leafset=label,source=src,n=nt,
                               R2_theta=r2t,R2_phi=r2p,r_cc_phi=rcc(d.phi_gold,d[pc]),
                               med_dtheta=medabs(d.theta_gold-d[tc]),
                               med_dphi=medabs(wrap(d.phi_gold-d[pc]))))
t2=pd.DataFrame(t2rows)
print(t2.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
(PAPER / "generated").mkdir(exist_ok=True)
t2.to_csv(PAPER / "generated" / "B1_table2_by_leafset.csv",index=False)

print("\n  Whole-plant geometry (e4) — INDEPENDENT of leaf set:")
for sp,p in [("sorghum",BK),("maize",MZ)]:
    e4=pd.read_csv(p+"/e4_geometric_fidelity.csv")
    print(f"   {sp}: median Chamfer  S<->V={e4.vo_cd_sym_cm.median():.2f}  P<->V={e4.vp_cd_sym_cm.median():.2f} cm"
          f" | median Hausdorff S<->V={e4.vo_hd_sym_cm.median():.2f}  P<->V={e4.vp_hd_sym_cm.median():.2f} cm  (N={len(e4)})")

# ---------------- Table 3 ----------------
print("\n"+"="*92)
print("B1 / TABLE 3  — broad-sense H2 under three leaf sets (sorghum n=2, maize n=4)")
print("="*92)
SRC=[("gold","theta_gold","phi_gold","voxel"),("skel","theta_ours","phi_ours","skeleton"),("proc","theta_rt","phi_rt","procedural")]
t3rows=[]
for sp in ("sorghum","maize"):
    m=load(sp); m5=m[m.leaf_index<5].copy()
    if sp=="sorghum":
        gmap=sorghum_geno_map(m5.plant_id.unique()); m5["geno"]=m5.plant_id.map(gmap)
        m5=m5[m5.geno.notna()&(m5.geno!="PI656058")]; m5["unit"]=m5.plant_id; nrep=2
    else:
        m5["geno"]=m5.plant_id.map(mz_geno); m5["unit"]=m5.plant_id.map(mz_plant)
        m5=m5[m5.geno.notna()]; nrep=4
    for which,label in [("a","(a) own-retained [published]"),("b","(b) both-keep"),("c","(c) all-common")]:
        for src,tc,pc,nm in SRC:
            d=m5[mask_for(m5,which,src)].copy()
            d["ph"]=norm360(d[pc])
            w=d.pivot_table(index=["unit","geno"],columns="leaf_index",values=[tc,"ph"],aggfunc="mean")
            w.columns=[f"{a}_{b}" for a,b in w.columns]; w=w.reset_index()
            traits=[]
            for i in LOWER5:
                c=f"{tc}_{i}"
                if c in w: traits.append(c)
            dev=[]
            for i in range(1,5):
                a,b=f"ph_{i}",f"ph_{i-1}"
                if a in w and b in w:
                    w[f"Phi_{i}"]=norm180(w[a]-w[b]); dev.append(f"Phi_{i}")
            if traits: w["med_theta"]=w[traits].median(axis=1); w["am_theta"]=w[traits].mean(axis=1)
            if dev: w["med_phi"]=w[dev].median(axis=1); w["am_phi"]=w[dev].mean(axis=1)
            allt=traits+dev+[c for c in ("med_theta","am_theta","med_phi","am_phi") if c in w]
            if sp=="maize":   # collapse timepoints -> one row per physical plant
                w=w.groupby(["unit","geno"])[allt].median().reset_index()
            reps=w.groupby("geno").size(); w=w[w.geno.isin(reps[reps>=2].index)]
            h=np.array([mom_h2(w[c].to_numpy(),w.geno.to_numpy(),nrep) for c in allt],float)
            t3rows.append(dict(species=sp,leafset=label,source=nm,n_traits=len(allt),
                               n_plants=len(w),n_geno=w.geno.nunique(),
                               ge20=int(np.nansum(h>=0.20)),median_H2=np.nanmedian(h),mean_H2=np.nanmean(h)))
t3=pd.DataFrame(t3rows)
print(t3.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
t3.to_csv(PAPER / "generated" / "B1_table3_by_leafset.csv",index=False)
print("\nsaved B1_table2_by_leafset.csv, B1_table3_by_leafset.csv")
print("DONE")
