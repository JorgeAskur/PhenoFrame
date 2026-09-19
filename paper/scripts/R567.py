"""R5 (self-consistency residual mechanism), R6 (flip-excluded r_cc on both-keep), R7 (maize |dtheta| by a_D bin)."""
import numpy as np, pandas as pd, warnings
from _data_paths import MAIZE, SORGHUM
warnings.simplefilter("ignore")
MZ=str(MAIZE)
PY=str(SORGHUM)
BK=PY
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
def dphi(g,o): return np.abs(((pd.to_numeric(g,errors="coerce")-pd.to_numeric(o,errors="coerce")+180)%360)-180)

print("="*92); print("R5  Self-consistency leaf-angle residual — mechanism search"); print("="*92)
sc=pd.read_csv(PY+"/descriptor_vs_measured_all_leaves.csv")
sc["res"]=pd.to_numeric(sc.py_inclination,errors="coerce")-pd.to_numeric(sc.xml_angle,errors="coerce")
sc["ares"]=sc.res.abs()
print(f"  n={len(sc)}  median|res|={sc.ares.median():.3f}  p99={sc.ares.quantile(.99):.2f}  min={sc.res.min():.2f}  max={sc.res.max():.2f}")
def show(col,bins,label,fmt="{:.0f}"):
    print(f"\n  by {label}:")
    print(f"   {'bin':>16} {'n':>6} {'median|res|':>12} {'p99':>8} {'max':>8}")
    s=sc.copy(); s["b"]=pd.cut(pd.to_numeric(s[col],errors="coerce"),bins)
    for b,g in s.groupby("b"):
        if len(g): print(f"   {str(b):>16} {len(g):>6} {g.ares.median():>12.3f} {g.ares.quantile(.99):>8.2f} {g.ares.max():>8.2f}")
show("xml_angle",[-1,5,15,30,45,60,75,85,91],"declared theta (deg)")
show("xml_splinePoints",[0,4,8,16,24,32,41],"n_i (splinePoints)")
show("xml_azimuth",[-181,-135,-90,-45,0,45,90,135,181],"declared azimuth (deg)")
show("leaf_id",[-1,0,1,2,3,4,6,8,20],"leaf rank")
# droopiness already ruled out; include for completeness
show("xml_droopiness",[-40,-30,-20,-10,0,10,20],"declared droopiness (deg)")
print("\n  20 worst leaves:")
w=sc.reindex(sc.ares.sort_values(ascending=False).index).head(20)
cols=["plant","leaf_id","xml_angle","xml_droopiness","xml_splinePoints","xml_azimuth","xml_length","py_inclination","res"]
print(w[cols].to_string(index=False,float_format=lambda v:f"{v:.3f}"))
# correlation summary
for c in ("xml_angle","xml_splinePoints","xml_droopiness","xml_azimuth","leaf_id","xml_length"):
    x=pd.to_numeric(sc[c],errors="coerce"); m=np.isfinite(x)&np.isfinite(sc.ares)
    print(f"  corr(|res|,{c:18}) = {np.corrcoef(x[m],sc.ares[m])[0,1]:+.3f}")

print("\n"+"="*92); print("R6  Flip-excluded r_cc(phi) and R2(phi) on BOTH-KEEP leaves"); print("="*92)
print(f"  {'species':8} {'mode':11} {'subset':22} {'n':>6} {'n_excl':>7} {'R2(phi)':>9} {'r_cc':>8}")
for sp,e1bp,e5p in [("sorghum",PY,BK),("maize",MZ,MZ)]:
    a=pd.read_csv(e1bp+"/e1b_per_leaf_pymaize_vs_gold.csv"); b=pd.read_csv(e5p+"/e5_roundtrip_vs_gold.csv")
    m=a[["plant_id","leaf_index","phi_gold","phi_ours","keep"]].merge(
        b[["plant_id","leaf_index","phi_rt","keep"]].rename(columns={"keep":"keep_p"}),
        on=["plant_id","leaf_index"],how="inner")
    m=m[(m.leaf_index<4)&(m.keep==True)&(m.keep_p==True)]
    for mode,pc in [("skeleton","phi_ours"),("procedural","phi_rt")]:
        d=m.copy(); fl=dphi(d.phi_gold,d[pc])>90
        r2a,na=R2(d.phi_gold,d[pc]); rca=rcc(d.phi_gold,d[pc])
        print(f"  {sp:8} {mode:11} {'all both-keep':22} {na:>6} {0:>7} {r2a:>9.3f} {rca:>8.3f}")
        d2=d[~fl]; r2b,nb=R2(d2.phi_gold,d2[pc]); rcb=rcc(d2.phi_gold,d2[pc])
        print(f"  {sp:8} {mode:11} {'flip-excluded':22} {nb:>6} {int(fl.sum()):>7} {r2b:>9.3f} {rcb:>8.3f}")

print("\n"+"="*92); print("R7  Maize procedural |dtheta| split by |a_D| bin (does the penalty sit at the bound?)"); print("="*92)
rt=pd.read_csv(MZ+"/e5_roundtrip_vs_gold.csv"); fp=pd.read_csv(MZ+"/e5_fitted_params.csv")
j=rt.merge(fp[["plant_id","leaf_index","droopiness","leaf_angle"]],on=["plant_id","leaf_index"],how="inner")
e1b=pd.read_csv(MZ+"/e1b_per_leaf_pymaize_vs_gold.csv")
j=j.merge(e1b[["plant_id","leaf_index","theta_ours","keep"]],on=["plant_id","leaf_index"],how="inner")
j=j[(j.leaf_index<4)&(j.keep==True)&(j.keep=="True" if False else True)]
if "keep" in rt.columns:
    j=j[j.keep_x==True] if "keep_x" in j.columns else j
j["ad"]=pd.to_numeric(j.droopiness,errors="coerce").abs()
j["dth_p"]=(pd.to_numeric(j.theta_gold,errors="coerce")-pd.to_numeric(j.theta_rt,errors="coerce")).abs()
j["dth_s"]=(pd.to_numeric(j.theta_gold,errors="coerce")-pd.to_numeric(j.theta_ours,errors="coerce")).abs()
print(f"  overall (lower-4, both-keep-ish n={len(j)}): procedural median|dtheta|={j.dth_p.median():.2f}  skeleton={j.dth_s.median():.2f}")
print(f"\n   {'|a_D| bin':>14} {'n':>7} {'frac':>6} {'proc med|dth|':>14} {'skel med|dth|':>14} {'penalty':>9}")
for lo,hi,lab in [(0,50,"0-50"),(50,100,"50-100"),(100,200,"100-200"),(200,299.5,"200-300"),(299.5,1e9,"at bound")]:
    s=j[(j.ad>=lo)&(j.ad<hi)]
    if len(s):
        print(f"   {lab:>14} {len(s):>7} {len(s)/len(j):>6.3f} {s.dth_p.median():>14.2f} {s.dth_s.median():>14.2f} {s.dth_p.median()-s.dth_s.median():>9.2f}")
sub=j[j.ad<299.5]
print(f"\n  excluding at-bound fits: procedural median|dtheta| = {sub.dth_p.median():.2f} (vs {j.dth_p.median():.2f} overall), n={len(sub)}")
print("DONE")
