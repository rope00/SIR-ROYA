#!/usr/bin/env python3
from pathlib import Path
import sys,time
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
sys.path.insert(0,str(Path(__file__).resolve().parent))
from simulate_v4 import RustModel,Parameters,weather_array
BASE=Path('results/convergence');BASE.mkdir(parents=True,exist_ok=True)
df=pd.read_csv('data/clima_diario_tres_localidades_qc.csv',parse_dates=['date'])
localities=['Socoroma','Ticnamar','Chapiquina']
climate={loc:weather_array(df,loc) for loc in localities}
cases=[('reference',0.05,80,1.25),('low_risk_near_tie',0.025,90,0.75),('near_tie_moderate',0.05,70,0.75),('near_tie_highRH',0.10,90,1.25),('high_intensity',0.15,70,2.00),('intermediate_chapiquina',0.10,80,2.00),('high_beta_high_RH',0.15,90,2.00)]
n_grid=[20,50,100,200,500];rows=[];started=time.monotonic()
for pi,(case,b,h,d) in enumerate(cases):
 par=Parameters(beta=b,rh50=h,dispersal_scale=d)
 for loc in localities:
  T,RH,PP=climate[loc]
  for rep in range(500):
   seed=20261110+pi*100000+rep
   result,_=RustModel(T,RH,PP,par,seed=seed).run(record=False)
   rows.append({'scenario':case,'beta':b,'rh50':h,'dispersal_scale':d,'locality':loc,'replicate':rep,'seed':seed,'I_final':result['I_final'],'cumulative_new_exposed':result['cumulative_new_exposed']})
 print(f'{case} done {time.monotonic()-started:.1f}s',flush=True)
all_runs=pd.DataFrame(rows);all_runs.to_csv(BASE/'mc_10500_realizations.csv',index=False)
summary=[]
for (case,loc),g in all_runs.groupby(['scenario','locality']):
 g=g.sort_values('replicate')
 for m in ['I_final','cumulative_new_exposed']:
  y=g[m].to_numpy(dtype=float);mu500=y.mean()
  for n in n_grid:
   a=y[:n];mu=a.mean();sd=a.std(ddof=1);se=sd/np.sqrt(n);half=student_t.ppf(.975,n-1)*se
   summary.append({'scenario':case,'locality':loc,'metric':m,'n':n,'mean':mu,'sd':sd,'mcse':se,'ci95_low':mu-half,'ci95_high':mu+half,'ci_halfwidth':half,'ci_halfwidth_relative_pct':100*half/mu if mu>0 else np.nan,'deviation_from_n500':mu-mu500,'abs_deviation_pct_of_n500':100*abs(mu-mu500)/mu500 if mu500 else np.nan})
summary_df=pd.DataFrame(summary);summary_df.to_csv(BASE/'convergence_prefix_summary.csv',index=False)
comparisons=[]
for case in all_runs.scenario.unique():
 X=all_runs.query('scenario==@case').pivot(index='replicate',columns='locality',values='cumulative_new_exposed')
 for n in n_grid:
  z=X.iloc[:n];means=z.mean();first=means.idxmax();runnerup=means.nlargest(2).index[-1]
  for a,b in [('Chapiquina','Socoroma'),('Chapiquina','Ticnamar'),('Socoroma','Ticnamar')]:
   da=z[a].to_numpy(float);db=z[b].to_numpy(float);diff=da-db
   paired_mean=diff.mean();paired_se=diff.std(ddof=1)/np.sqrt(n);paired_half=student_t.ppf(.975,n-1)*paired_se
   independent_se=np.sqrt(da.var(ddof=1)/n+db.var(ddof=1)/n);independent_half=1.96*independent_se
   comparisons.append({'scenario':case,'n':n,'a':a,'b':b,'leader':first,'runnerup':runnerup,'difference_mean':paired_mean,'paired_ci95_low':paired_mean-paired_half,'paired_ci95_high':paired_mean+paired_half,'independent_ci95_low':paired_mean-independent_half,'independent_ci95_high':paired_mean+independent_half,'paired_excludes_zero':(paired_mean-paired_half)>0 or (paired_mean+paired_half)<0,'independent_excludes_zero':(paired_mean-independent_half)>0 or (paired_mean+independent_half)<0})
pd.DataFrame(comparisons).to_csv(BASE/'pairwise_mc_convergence.csv',index=False)
B=3000;rng=np.random.default_rng(20261111);rank=[]
for case in all_runs.scenario.unique():
 X=all_runs.query('scenario==@case').pivot(index='replicate',columns='locality',values='cumulative_new_exposed')[localities]
 for n in n_grid:
  z=X.iloc[:n].to_numpy(float);means=z.mean(axis=0);leader_idx=int(np.argmax(means))
  idx=rng.integers(n,size=(B,n));boot_mean=z[idx].mean(axis=1)
  paired_prob=(np.argmax(boot_mean,axis=1)==leader_idx).mean()
  indep_mean=np.column_stack([z[rng.integers(n,size=(B,n)),j].mean(axis=1) for j in range(3)])
  indep_prob=(np.argmax(indep_mean,axis=1)==leader_idx).mean()
  rank.append({'scenario':case,'n':n,'leader':localities[leader_idx],'gap_top_two':np.sort(means)[-1]-np.sort(means)[-2],'paired_bootstrap_support':paired_prob,'independent_bootstrap_support':indep_prob})
rank_df=pd.DataFrame(rank);rank_df.to_csv(BASE/'winner_stability_by_n.csv',index=False)
print('RANKING n=500');print(rank_df[rank_df.n==500].to_string(index=False))
print('REFERENCE');print(summary_df.query("scenario=='reference' and metric=='cumulative_new_exposed'").to_string(index=False));print('Saved',BASE)
