"""R4 — refit the procedural model with polygon vs arc-length leafLength, 30 plants/species.
leaf_length is an INPUT to the fit (held fixed while theta, a_D are optimised), so changing its
source can move the fitted angles. Compares fitted theta, droopiness and residual."""
import importlib, os, sys, tempfile, warnings
from pathlib import Path
from _data_paths import GENERATED, MAIZE, REPO, SORGHUM
import numpy as np, pandas as pd
warnings.simplefilter("ignore")
sys.path.insert(0,str(REPO/"experiments")); sys.path.insert(0,str(REPO))
from scipy.optimize import minimize
import phenoframe.skeleton_to_descriptor as S
from phenoframe import compute_traits_from_descriptor

BOUNDS=[(0.5,89.5),(-300.0,300.0)]
def fit_one(obs_local, L, spline_points, stem_radius):
    obs_r=S._resample_polyline(obs_local, S._N_RESAMPLE)
    def cost(p):
        la,dr=p
        gen=S._forward_leaf_local(la,dr,L,spline_points,stem_radius,0.0)
        return float(np.sum((obs_r-S._resample_polyline(gen,S._N_RESAMPLE))**2))
    best=None
    for d0 in np.linspace(-250,250,S._N_STARTS):
        try:
            r=minimize(cost,x0=[30.0,d0],method="L-BFGS-B",bounds=BOUNDS)
            if best is None or r.fun<best.fun: best=r
        except Exception: pass
    return (best.x[0],best.x[1],float(best.fun)) if best is not None else (np.nan,np.nan,np.nan)

def run(species, voxel_root, ids, k=30):
    os.environ["PHENOFRAME_VOXEL_PATH"]=voxel_root
    from experiments import paths, data_loaders as dl
    importlib.reload(paths)
    importlib.reload(dl)
    tmp=Path(tempfile.mkdtemp()); rows=[]; done=0
    for pid in ids:
        if done>=k: break
        try:
            vox=dl.load_voxel_grid(pid)
            ctrl=S.skeleton_to_override_xml(vox, tmp/"o.xml")
            tr=compute_traits_from_descriptor(tmp/"o.xml")
            arc={t["leaf_index"]:t["leaf_length"] for t in tr}
            for i,c in enumerate(ctrl):
                c=np.asarray(c,float)
                if len(c)<4 or i not in arc or arc[i]<=0: continue
                Lp=float(np.sum(np.linalg.norm(np.diff(c,axis=0),axis=1)))   # polygon (Eq.4, current)
                La=float(arc[i])                                             # measured B-spline arc
                if Lp<=0: continue
                ap,dp,rp=fit_one(c,Lp,len(c),S.DEFAULT_STEM_RADIUS)
                aa,da,ra=fit_one(c,La,len(c),S.DEFAULT_STEM_RADIUS)
                rows.append(dict(species=species,plant=pid,leaf=i,L_poly=Lp,L_arc=La,
                                 rel=(Lp-La)/Lp,th_poly=ap,th_arc=aa,ad_poly=dp,ad_arc=da,
                                 res_poly=rp,res_arc=ra))
            done+=1
        except Exception: continue
    return pd.DataFrame(rows), done

out=[]
for sp,vr,src in [
    ("sorghum", os.environ.get("PHENOFRAME_SORGHUM_VOXEL_PATH"), SORGHUM / "e5_fitted_params.csv"),
    ("maize", os.environ.get("PHENOFRAME_MAIZE_VOXEL_PATH"), MAIZE / "e5_fitted_params.csv"),
]:
    if not vr:
        raise RuntimeError(f"Set PHENOFRAME_{sp.upper()}_VOXEL_PATH to rerun this raw-voxel sensitivity check")
    ids=pd.read_csv(src).plant_id.unique().tolist()
    df,n=run(sp,vr,ids,30)
    if len(df)==0:
        print(f"[{sp}] no leaves fitted"); continue
    df["dth"]=(df.th_arc-df.th_poly).abs(); df["dad"]=(df.ad_arc-df.ad_poly).abs()
    print(f"\n[{sp}] plants={n} leaves={len(df)}  median (L_poly-L_arc)/L_poly = {100*df.rel.median():.2f}%")
    print(f"   |d theta| : median={df.dth.median():.3f} deg  p95={df.dth.quantile(.95):.3f} deg  max={df.dth.max():.2f}")
    print(f"   |d a_D|   : median={df.dad.median():.3f} deg  p95={df.dad.quantile(.95):.3f} deg  max={df.dad.max():.2f}")
    print(f"   residual  : polygon median={df.res_poly.median():.3e}  arc median={df.res_arc.median():.3e}  "
          f"(arc better on {100*(df.res_arc<df.res_poly).mean():.1f}% of leaves)")
    print(f"   fits at a_D bound: polygon={100*(df.ad_poly.abs()>=299.5).mean():.1f}%  arc={100*(df.ad_arc.abs()>=299.5).mean():.1f}%")
    out.append(df)
if out:
    GENERATED.mkdir(exist_ok=True)
    pd.concat(out).to_csv(GENERATED / "R4_leaflength_refit.csv",index=False)
    print("\nsaved R4_leaflength_refit.csv")
print("DONE")
