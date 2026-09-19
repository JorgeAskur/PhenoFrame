"""R8 — recompute on the both-keep basis (theta0/theta4, bootstrap, per-trait 95% CI, harmonic-n).
R9 — format the generator SI tables (no re-run)."""
import os, re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
warnings.simplefilter("ignore")
PAPER = Path(__file__).resolve().parents[1]
MZ = str(PAPER / "source_data" / "maize")
PY = str(PAPER / "source_data" / "sorghum")
BK = PY
_RXM=re.compile(r"_\d+-(\d+)-([A-Za-z0-9]+)-(\d+)_")
RNG=np.random.default_rng(11); NB=1000
def n360(x): return np.mod(pd.to_numeric(x,errors="coerce"),360.0)
def n180(x): return np.abs(x-180.0)
def mom(y,g,n):
    y=np.asarray(y,float); m=np.isfinite(y); y=y[m]; g=np.asarray(g)[m]
    gg,gi=np.unique(g,return_inverse=True); k=len(gg); N=len(y)
    if k<2 or N-k<=0: return np.nan
    ni=np.bincount(gi); gm=np.bincount(gi,weights=y)/ni
    msb=np.sum(ni*(gm-y.mean())**2)/(k-1); msw=np.sum((y-gm[gi])**2)/(N-k)
    n0=(N-np.sum(ni**2)/N)/(k-1); sg=max((msb-msw)/n0,0.0)
    d=sg+msw/n; return sg/d if d>0 else 0.0
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
SRC=[("gold","theta_gold","phi_gold","voxel"),("skel","theta_ours","phi_ours","skeleton"),("proc","theta_rt","phi_rt","procedural")]
LBL=["theta_0","theta_1","theta_2","theta_3","theta_4","Phi_1","Phi_2","Phi_3","Phi_4","med_theta","am_theta","med_phi","am_phi"]

def build(species):
    if species=="sorghum":
        a=pd.read_csv(PY+"/e1b_per_leaf_pymaize_vs_gold.csv"); b=pd.read_csv(BK+"/e5_roundtrip_vs_gold.csv")
    else:
        a=pd.read_csv(MZ+"/e1b_per_leaf_pymaize_vs_gold.csv"); b=pd.read_csv(MZ+"/e5_roundtrip_vs_gold.csv")
    m=a[["plant_id","leaf_index","theta_gold","phi_gold","theta_ours","phi_ours","keep"]].merge(
        b[["plant_id","leaf_index","theta_rt","phi_rt","keep"]].rename(columns={"keep":"keep_p"}),
        on=["plant_id","leaf_index"],how="inner")
    m=m[(m.leaf_index<5)&(m.keep==True)&(m.keep_p==True)].copy()   # BOTH-KEEP
    if species=="sorghum":
        gm=sorghum_geno(m.plant_id.unique()); m["geno"]=m.plant_id.map(gm)
        m=m[m.geno.notna()&(m.geno!="PI656058")]; m["unit"]=m.plant_id
    else:
        m["geno"]=m.plant_id.map(lambda p:(lambda x:x.group(2) if x else None)(_RXM.search(p)))
        m["unit"]=m.plant_id.map(lambda p:(lambda x:f"{x.group(2)}-{x.group(3)}" if x else None)(_RXM.search(p)))
        m=m[m.geno.notna()]
    return m

def traits(m,species,tc,pc):
    d=m.copy(); d["ph"]=n360(d[pc])
    # intersection mask for Phi
    ok={}
    for _,_,p2,_ in SRC:
        w2=m.pivot_table(index="plant_id",columns="leaf_index",values=p2,aggfunc="first")
        w2=w2.apply(lambda c:n360(c))
        for i in range(1,5):
            if i in w2.columns and (i-1) in w2.columns:
                r=np.mod(w2[i]-w2[i-1],360.0); ok[(p2,i)]=((r>=90)&(r<=270))
    w=d.pivot_table(index=["plant_id","unit","geno"],columns="leaf_index",values=[tc,"ph"],aggfunc="first")
    w.columns=[f"{a}_{b}" for a,b in w.columns]; w=w.reset_index()
    th=[f"{tc}_{i}" for i in range(5) if f"{tc}_{i}" in w]; dev=[]
    for i in range(1,5):
        a,b=f"ph_{i}",f"ph_{i-1}"
        if a not in w or b not in w: continue
        raw=np.mod(w[a]-w[b],360.0); good=((raw>=90)&(raw<=270)).to_numpy()
        for _,_,p2,_ in SRC:
            if (p2,i) in ok:
                good=good&ok[(p2,i)].reindex(w.plant_id).fillna(False).to_numpy()
        c=f"Phi_{i}"; w[c]=np.where(good,n180(raw),np.nan); dev.append(c)
    if th: w["med_theta"]=w[th].median(axis=1); w["am_theta"]=w[th].mean(axis=1)
    if dev: w["med_phi"]=w[dev].median(axis=1); w["am_phi"]=w[dev].mean(axis=1)
    allt=th+dev+[c for c in("med_theta","am_theta","med_phi","am_phi") if c in w]
    if species=="maize": w=w.groupby(["unit","geno"])[allt].median().reset_index()
    else: w=w[["unit","geno"]+allt]
    return w,allt

print("="*94); print("R8  Both-keep basis (Phi on intersection). Point=REML from R123; CI/bootstrap=MoM."); print("="*94)
defin=pd.read_csv(PAPER / "source_data" / "summary" / "R123_table3_definitive.csv")
print("\n  (a) theta_0 and theta_4 heritabilities (REML, both-keep, intersection):")
for sp in ("sorghum","maize"):
    for t in ("theta_0","theta_4"):
        r=defin[(defin.species==sp)&(defin.trait==t)]
        if len(r):
            v={x.source:x.H2 for x in r.itertuples()}
            print(f"    {sp:8} {t:8} voxel={v.get('voxel',np.nan):.3f}  skeleton={v.get('skeleton',np.nan):.3f}  procedural={v.get('procedural',np.nan):.3f}")

allres={}
for sp in ("sorghum","maize"):
    m=build(sp); nrep=4 if sp=="maize" else 2
    tabs={}
    for _,tc,pc,nm in SRC:
        w,allt=traits(m,sp,tc,pc); tabs[nm]=(w,allt)
    allres[sp]=(tabs,nrep)
    genos=tabs["voxel"][0].geno.unique()
    # (b) bootstrap median-H2 comparison
    by={nm:{g:tabs[nm][0][tabs[nm][0].geno==g] for g in genos} for nm in tabs}
    beats=0; boot={nm:np.empty(NB) for nm in tabs}
    for b in range(NB):
        pick=RNG.choice(genos,genos.size,replace=True)
        med={}
        for nm,(w,allt) in tabs.items():
            parts=[]
            for j,g in enumerate(pick):
                p=by[nm][g].copy(); p["geno"]=f"{g}__{j}"; parts.append(p)
            rs=pd.concat(parts,ignore_index=True)
            med[nm]=np.nanmedian([mom(rs[c].to_numpy(),rs.geno.to_numpy(),nrep) for c in allt])
            boot[nm][b]=med[nm]
        if med["procedural"]>med["voxel"]: beats+=1
    print(f"\n  (b) {sp}: P(procedural median H2 > voxel median H2) = {beats/NB:.3f}  [both-keep, {NB} genotype resamples, MoM]")
    for nm in ("voxel","skeleton","procedural"):
        lo,hi=np.nanpercentile(boot[nm],[2.5,97.5])
        print(f"        {nm:11} median-H2 bootstrap 95% CI = [{lo:.3f}, {hi:.3f}]")
    # (c) per-trait 95% CI
    rows=[]
    for nm,(w,allt) in tabs.items():
        for k,c in enumerate(allt):
            pt=mom(w[c].to_numpy(),w.geno.to_numpy(),nrep)
            bs=np.empty(NB)
            for b in range(NB):
                pick=RNG.choice(genos,genos.size,replace=True)
                parts=[]
                for j,g in enumerate(pick):
                    p=by[nm][g].copy(); p["geno"]=f"{g}__{j}"; parts.append(p)
                rs=pd.concat(parts,ignore_index=True)
                bs[b]=mom(rs[c].to_numpy(),rs.geno.to_numpy(),nrep)
            lo,hi=np.nanpercentile(bs,[2.5,97.5])
            rows.append(dict(species=sp,source=nm,trait=LBL[k] if k<len(LBL) else c,H2_mom=pt,lo=lo,hi=hi,width=hi-lo))
    (PAPER / "generated").mkdir(exist_ok=True)
    df=pd.DataFrame(rows); df.to_csv(PAPER / "generated" / f"R8_per_trait_CI_{sp}.csv",index=False)
    print(f"  (c) {sp}: per-trait 95% CIs saved -> R8_per_trait_CI_{sp}.csv ; median CI width by source: "
          + ", ".join(f"{nm}={df[df.source==nm].width.median():.3f}" for nm in ("voxel","skeleton","procedural")))
    # (d) harmonic-mean n
    reps=tabs["voxel"][0].groupby("geno").size(); reps=reps[reps>=2]
    nh=len(reps)/np.sum(1.0/reps)
    print(f"  (d) {sp}: replicates/genotype arithmetic={reps.mean():.3f} harmonic={nh:.3f} (range {reps.min()}-{reps.max()})")
    for nm,(w,allt) in tabs.items():
        hn=np.nanmedian([mom(w[c].to_numpy(),w.geno.to_numpy(),nh) for c in allt])
        hd=np.nanmedian([mom(w[c].to_numpy(),w.geno.to_numpy(),nrep) for c in allt])
        print(f"        {nm:11} median H2: n={nrep} -> {hd:.3f} | n_harmonic={nh:.2f} -> {hn:.3f}")

print("\n"+"="*94); print("R9  Generator SI tables (format only)"); print("="*94)
syn=pd.read_csv(PAPER / "source_data" / "generation" / "sec2_5_synthetic_traits.csv"); real=pd.read_csv(PY+"/e1b_per_leaf_pymaize_vs_gold.csv")
def dev(df,c):
    d=df.sort_values(["plant_id","leaf_index"]).copy(); d["pn"]=n360(d[c])
    d["nx"]=d.groupby("plant_id")["pn"].shift(-1); d["Phi"]=n180(np.mod(d.nx-d.pn,360.0)); return d
rs=dev(real,"phi_gold"); ss=dev(syn,"phi")
print("\n  SI Table A — per-node distributional agreement, synthetic vs real (sorghum, lower-5)")
print(f"   {'node':>4} | {'theta KS':>9} {'theta W1 (deg)':>15} | {'Phi KS':>7} {'Phi W1 (deg)':>13}")
for i in range(5):
    rt=pd.to_numeric(real[real.leaf_index==i].theta_gold,errors="coerce").dropna()
    st=pd.to_numeric(syn[syn.leaf_index==i].theta,errors="coerce").dropna()
    rp=rs[rs.leaf_index==i].Phi.dropna(); sp_=ss[ss.leaf_index==i].Phi.dropna()
    print(f"   {i:>4} | {stats.ks_2samp(rt,st).statistic:>9.3f} {stats.wasserstein_distance(rt,st):>15.2f} | "
          f"{stats.ks_2samp(rp,sp_).statistic:>7.3f} {stats.wasserstein_distance(rp,sp_):>13.2f}")
def wide(df,c):
    w=df[df.leaf_index<5].pivot_table(index="plant_id",columns="leaf_index",values=c,aggfunc="first")
    return w.reindex(columns=range(5))
rw=wide(real,"theta_gold"); sw=wide(syn,"theta")
print("\n  SI Table B — between-node theta correlation, REAL (lower triangle) ")
print(np.round(rw.corr().values,3))
print("\n  SI Table B — between-node theta correlation, SYNTHETIC")
print(np.round(sw.corr().values,3))
tri=np.triu_indices(5,1)
print(f"\n   mean |off-diagonal|: real={np.abs(rw.corr().values[tri]).mean():.3f}  synthetic={np.abs(sw.corr().values[tri]).mean():.3f}")
print("DONE")
