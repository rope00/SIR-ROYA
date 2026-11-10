#!/usr/bin/env python3
"""Reproducible oregano-rust v4 scenario simulations.

This reproduces the S-E-I benchmark and factorial parameter sweep in the paper.

Example:
    python code/simulate_v4.py --input data/clima_diario_tres_localidades_qc.csv --output results/reproduced --benchmark
    python code/simulate_v4.py --input data/clima_diario_tres_localidades_qc.csv --output results/reproduced --full

No field parameter calibration or observed disease outcome is claimed.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import convolve
from scipy.stats import spearmanr

LOCATIONS=('Socoroma','Ticnamar','Chapiquina')
DATES=pd.date_range('2021-01-01',periods=100,freq='D')

@dataclass(frozen=True)
class Parameters:
    beta:float=.05
    rh50:float=80.
    dispersal_scale:float=1.25
    t_min:float=4.  # operational assumption, not measured threshold
    t_opt:float=20.
    t_max:float=27.
    rh_scale:float=7.5
    max_radius:int=3

def temp_suit(T:float,p:Parameters)->float:
    if T<=p.t_min or T>=p.t_max:return 0.
    return (T-p.t_min)/(p.t_opt-p.t_min) if T<=p.t_opt else (p.t_max-T)/(p.t_max-p.t_opt)

def latent_period(T:float)->float:
    limited=float(np.clip(T,5,22))
    return float(np.clip(253.291*limited**(-1.056),10,50))

def spatial_kernel(p:Parameters)->np.ndarray:
    r=p.max_radius
    mat=np.zeros((2*r+1,2*r+1),float)
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            d=np.hypot(x,y)
            if 0<d<=r:mat[x+r,y+r]=np.exp(-d/p.dispersal_scale)
    return mat

def weather_array(data:pd.DataFrame,locality:str):
    x=data[(data.locality==locality)&(data.date.isin(DATES))].sort_values('date')
    if len(x)!=100 or not x[['temp_mean_c','rh_mean_pct']].notna().all().all():
        raise ValueError(f'{locality}: missing/invalid benchmark weather')
    return (x.temp_mean_c.to_numpy(float),x.rh_mean_pct.to_numpy(float),x.precip_mm.to_numpy(float))

class RustModel:
    def __init__(self,T,RH,PP,parameters:Parameters,seed:int,
                 size:int=50,plant_occupancy:float=.8,initial_fraction:float=.01):
        self.T,self.RH,self.PP=T,RH,PP
        self.params=parameters
        self.rng=np.random.default_rng(seed)
        self.state=np.zeros((size,size),dtype=np.int8)
        self.state[self.rng.random((size,size))<plant_occupancy]=1
        candidates=np.argwhere(self.state==1)
        n_initial=max(1,int(round(len(candidates)*initial_fraction)))
        selected=candidates[self.rng.choice(len(candidates),n_initial,replace=False)]
        self.state[tuple(selected.T)]=3
        self.initial_infected=n_initial
        self.development=np.zeros_like(self.state,dtype=float)
        self.K=spatial_kernel(parameters)

    def run(self,record:bool=False):
        cumulative=0
        trajectory=[]
        for day,(T,RH,PP) in enumerate(zip(self.T,self.RH,self.PP),1):
            temperature=temp_suit(T,self.params)
            humidity=1/(1+np.exp(-(RH-self.params.rh50)/self.params.rh_scale))
            rainfall=0 if np.isnan(PP) else 1-np.exp(-max(PP,0)/2)
            env=temperature*max(humidity,rainfall)
            infectious_pressure=convolve((self.state==3).astype(float),self.K,mode='constant',cval=0)
            infection_prob=1-np.exp(-self.params.beta*env*infectious_pressure)
            new_E=(self.state==1)&(self.rng.random(self.state.shape)<infection_prob)
            self.state[new_E]=2
            self.development[new_E]=0
            cumulative+=int(new_E.sum())
            E=self.state==2
            self.development[E]+=1/latent_period(T)
            self.state[E&(self.development>=1)]=3
            if record:
                trajectory.append({'day':day,'S':int((self.state==1).sum()),
                    'E':int((self.state==2).sum()),'I':int((self.state==3).sum()),
                    'new_exposed':int(new_E.sum())})
        return {'I_final':int((self.state==3).sum()),
                'cumulative_new_exposed':cumulative,
                'affected_plants':self.initial_infected+cumulative},pd.DataFrame(trajectory)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',default='data/clima_diario_tres_localidades_qc.csv')
    p.add_argument('--output',default='results/reproduced')
    group=p.add_mutually_exclusive_group()
    group.add_argument('--benchmark',action='store_true',help='Benchmark: 100 replications per locality')
    group.add_argument('--full',action='store_true',help='Benchmark + 2160-run uncertainty sweep')
    group.add_argument('--smoke-test',action='store_true',help='2 replications per locality')
    args=p.parse_args()
    output=Path(args.output);output.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.input,parse_dates=['date'])
    climate={loc:weather_array(df,loc) for loc in LOCATIONS}
    n_rep=2 if args.smoke_test else 100
    benchmark=[];trajectories=[]
    for loc in LOCATIONS:
        T,RH,PP=climate[loc]
        for rep in range(n_rep):
            res,traj=RustModel(T,RH,PP,Parameters(),seed=20260903+rep).run(record=True)
            benchmark.append({'locality':loc,'replicate':rep,**res})
            traj['locality']=loc;traj['replicate']=rep;trajectories.append(traj)
    bb=pd.DataFrame(benchmark)
    bb.to_csv(output/'benchmark_replicates.csv',index=False)
    pd.concat(trajectories,ignore_index=True).to_csv(output/'benchmark_trajectories.csv',index=False)
    print('Benchmark mean infectious plants by locality:')
    print(bb.groupby('locality').I_final.mean())
    if not args.full:return

    grid=[(b,r,d) for b in [0.025,0.05,0.10,0.15]
          for r in [70.,80.,90.]
          for d in [0.75,1.25,2.0]]
    rows=[]
    for loc in LOCATIONS:
        T,RH,PP=climate[loc]
        for index,(b,r,d) in enumerate(grid):
            parameters=Parameters(beta=b,rh50=r,dispersal_scale=d)
            for rep in range(20):
                result,_=RustModel(T,RH,PP,parameters,seed=20260903+index*1000+rep).run()
                rows.append({'locality':loc,'beta':b,'rh50':r,'dispersal_scale':d,
                             'replicate':rep,**result})
    ss=pd.DataFrame(rows)
    ss.to_csv(output/'sensitivity_replicates.csv',index=False)
    avg=ss.groupby(['locality','beta','rh50','dispersal_scale'],as_index=False)['cumulative_new_exposed'].mean()
    avg.to_csv(output/'sensitivity_summary.csv',index=False)
    importance=[]
    for loc in LOCATIONS:
        z=avg[avg.locality==loc]
        for par in ['beta','rh50','dispersal_scale']:
            rho,pvalue=spearmanr(z[par],z.cumulative_new_exposed)
            importance.append({'locality':loc,'parameter':par,'spearman_rho':rho,'p_value':pvalue})
    pd.DataFrame(importance).to_csv(output/'parameter_importance.csv',index=False)
    rank=avg.pivot_table(index=['beta','rh50','dispersal_scale'],columns='locality',values='cumulative_new_exposed')
    winners=rank.idxmax(axis=1).value_counts()
    print('Locality with highest mean incidence by parameter set:')
    print(winners)
    winners.to_csv(output/'locality_winner_counts.csv')

if __name__=='__main__':main()
