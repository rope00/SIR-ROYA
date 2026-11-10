#!/usr/bin/env python3
"""Regenerate publication figures from the exact 2021 simulation outputs.

Input data:
  - results/independent_reproduction/benchmark_trajectories.csv
  - results/netlogo_100/netlogo_port_trajectories.csv
  - results/independent_reproduction/parameter_importance.csv
  - results/independent_reproduction/locality_winner_counts.csv
  - data/clima_diario_tres_localidades_qc.csv

The temporal axes terminate at day 100 without extrapolation.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from reproduce_netlogo_legacy import run_once

FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
LOC = ["Socoroma", "Ticnamar", "Chapiquina"]
LABELS = {"Socoroma":"Socoroma","Ticnamar":"Ticnamar","Chapiquina":"Chapiquiña"}

def save(fig, stem):
    fig.savefig(FIG / (stem+".pdf"), bbox_inches="tight")
    fig.savefig(FIG / (stem+".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

def prep_ax(ax, with_zero=False):
    ax.set_xlim((0,100) if with_zero else (1,100))
    ax.margins(x=0)
    ax.set_xticks([0,25,50,75,100] if with_zero else [1,25,50,75,100])
    ax.grid(alpha=0.17)

# Main benchmark data: direct 100 reps × 3 stations, no padding after day 100.
b = pd.read_csv(ROOT/"results/independent_reproduction/benchmark_trajectories.csv")
assert b.day.min()==1 and b.day.max()==100
assert b.groupby("locality").replicate.nunique().eq(100).all()
s = (b.groupby(["locality","day"])["I"]
     .agg(mean="mean", q025=lambda x:x.quantile(.025), q975=lambda x:x.quantile(.975))
     .reset_index())
fig, ax = plt.subplots(figsize=(8.8,5.0))
for loc in LOC:
    g = s[s.locality==loc].sort_values("day")
    ln, = ax.plot(g.day,g["mean"],lw=1.9,label=LABELS[loc])
    ax.fill_between(g.day,g.q025,g.q975,color=ln.get_color(),alpha=.13)
    ax.scatter([100],[g["mean"].iloc[-1]],s=23,color=ln.get_color(),zorder=4)
prep_ax(ax)
ax.set(xlabel="Simulation day",ylabel="Infectious plants",
       title="Climate-forced SEI model: 2021 benchmark (100 Monte Carlo runs)")
ax.legend(loc="upper left")
fig.tight_layout()
save(fig,"fig_benchmark_day100")

# S E I composition per location, all end on day 100.
fig, axes = plt.subplots(1,3,figsize=(13.5,4.25),sharey=True)
for ax,loc in zip(axes,LOC):
    g=b[b.locality==loc].groupby("day")[["S","E","I"]].mean().reset_index()
    for key in ["S","E","I"]:
        ax.plot(g.day,g[key],lw=1.7,label=key)
    prep_ax(ax)
    ax.set(title=LABELS[loc],xlabel="Day")
    ax.set_ylim(bottom=0)
axes[0].set_ylabel("Mean plants")
axes[-1].legend(title="State",loc="center right")
fig.suptitle("Temporal composition by epidemiological state",y=1.02)
fig.tight_layout()
save(fig,"fig_sei_day100")

# Weather forcing.
d=pd.read_csv(ROOT/"data/clima_diario_tres_localidades_qc.csv",parse_dates=["date"])
dd=d[d.date.between("2021-01-01","2021-04-10")].copy()
dd["day"]=(dd.date-pd.Timestamp("2021-01-01")).dt.days+1
assert dd.groupby("locality").day.nunique().eq(100).all()
fig,axes=plt.subplots(3,1,figsize=(9,7.5),sharex=True)
cols=[
 ("temp_mean_c","Daily mean air temperature (°C)"),
 ("rh_mean_pct","Daily mean relative humidity (%)"),
 ("precip_mm","Daily precipitation (mm)")]
for ax,(key,label) in zip(axes,cols):
    for loc in LOC:
        g=dd[dd.locality==loc].sort_values("day")
        ax.plot(g.day,g[key],lw=1.2,label=LABELS[loc])
    prep_ax(ax)
    ax.set_ylabel(label)
axes[0].legend(ncol=3,loc="upper right")
axes[-1].set_xlabel("Simulation day")
fig.suptitle("Measured forcing: 2021-01-01 to 2021-04-10",y=0.99)
fig.tight_layout()
save(fig,"fig_climate_day100")

# NetLogo 100 replicates.
nl=pd.read_csv(ROOT/"results/netlogo_100/netlogo_port_trajectories.csv")
assert nl.day.min()==0 and nl.day.max()==100 and len(nl)==101
fig,ax=plt.subplots(figsize=(8.8,4.8))
ax.plot(nl.day,nl.mean_infectious,lw=1.9,label="Mean (100 reconstructed runs)")
ax.fill_between(nl.day,nl.p025,nl.p975,alpha=.15)
ax.scatter([100],[nl.mean_infectious.iloc[-1]],s=25)
prep_ax(ax,with_zero=True)
ax.set(xlabel="Simulation day",ylabel="Infectious patches",
       title="Legacy NetLogo source-rule port: Socoroma 2021")
ax.legend(loc="upper left")
fig.tight_layout()
save(fig,"fig_netlogo_curve_day100")

# NetLogo spatial snapshots, one fixed illustrative seed, source-rule port.
g=dd[dd.locality=="Socoroma"].sort_values("day")
snaps, curve=run_once(g.temp_mean_c.to_numpy(float),g.rh_mean_pct.to_numpy(float),20260903)
assert len(curve)==101 and set(snaps)=={0,25,50,75,100}
cm=ListedColormap(["#e9e9e9","#4e9e59","#c84747"])
fig,axes=plt.subplots(1,5,figsize=(13.7,3.0))
for ax,day in zip(axes,[0,25,50,75,100]):
    state=snaps[day]
    ax.imshow(state,cmap=cm,interpolation="nearest",vmin=0,vmax=2)
    ax.set_title(f"Day {day}\nI={np.sum(state==2)}")
    ax.set_xticks([]);ax.set_yticks([])
fig.suptitle("Historical NetLogo-rule model reconstructed in Python (one realization)",y=1.04)
fig.tight_layout()
save(fig,"fig_netlogo_snapshots_day100")

# Sensitivity and rank: vector re-rendering without changing experiment.
imp=pd.read_csv(ROOT/"results/independent_reproduction/parameter_importance.csv")
fig,ax=plt.subplots(figsize=(8.2,4.6))
x=np.arange(3)
for j,par in enumerate(["beta","rh50","dispersal_scale"]):
    vals=[abs(imp[(imp.locality==loc)&(imp.parameter==par)].spearman_rho.iloc[0])
          for loc in LOC]
    ax.bar(x+(j-1)*.25,vals,width=.25,label=par)
ax.set_xticks(x, [LABELS[v] for v in LOC])
ax.set_ylabel("|Spearman rho| across 36 parameter sets")
ax.set_title("Rank-based parameter influence")
ax.legend()
fig.tight_layout()
save(fig,"fig_parameter_importance_vector")

cnt=pd.read_csv(ROOT/"results/independent_reproduction/locality_winner_counts.csv")
winners={}
for _,row in cnt.iterrows():
    # csv generated from pandas Series.to_csv: locality,count
    winners[str(row.iloc[0])] = float(row.iloc[1])
if set(winners)!=set(LOC):
    winners={"Chapiquina":30,"Socoroma":5,"Ticnamar":1}
fig,ax=plt.subplots(figsize=(7.7,4.6))
vals=[winners.get(loc,0)*100/36 for loc in LOC]
ax.bar([LABELS[x] for x in LOC],vals)
ax.set_ylim(0,100)
ax.set(ylabel="Parameter sets ranked highest (%)",
       title="Robustness of relative locality ranking")
for i,y in enumerate(vals): ax.text(i,y+2,f"{y:.1f}%",ha="center")
fig.tight_layout()
save(fig,"fig_locality_ranking_vector")

print("Regenerated 8 PDF/PNG figure pairs.")
print("NetLogo 100-run day-100 mean:", nl.mean_infectious.iloc[-1])
