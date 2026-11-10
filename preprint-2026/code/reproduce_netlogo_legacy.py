#!/usr/bin/env python3
"""Exploratory source-rule port of 2025 NetLogo `RoyaSocoromaModel.nlogo`.

This is NOT a native NetLogo run and not a validation of the new SEI model.
The original model loads `socoroma_clima_clean.csv`; this port instead uses the
independently quality-controlled Socoroma 2021 100-day climate, allowing comparison
of structural assumptions under common forcing. RNG streams are not identical
between NetLogo and NumPy. Source details: NetLogo 6.4, 33x33, density=78%,
3 initial infected, cardinal neighbors at radius=1, no mortality.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def old_infection(T,H):
    h=(H/12.0)**3
    return float(.5*np.exp(-((T-20.)**2)/(2*5.**2))*h/(1+h))


def run_once(T,RH,seed:int,frame_days=(0,25,50,75,100)):
    rng=np.random.default_rng(seed)
    state=np.zeros((33,33),np.int8)
    state[rng.random(state.shape)<.78]=1
    candidates=np.argwhere(state==1)
    initial=candidates[rng.choice(len(candidates),3,replace=False)]
    state[tuple(initial.T)]=2
    snaps={0:state.copy()}
    infected=[3]
    for day,(temp,hum) in enumerate(zip(T,RH),start=1):
        p=old_infection(temp,hum)
        # NetLogo's outer ask agentset is evaluated at tick start.
        infectious_at_tick_start=np.argwhere(state==2)
        infectious_at_tick_start=infectious_at_tick_start[rng.permutation(len(infectious_at_tick_start))]
        for x,y in infectious_at_tick_start:
            for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:
                nx,ny=x+dx,y+dy
                if 0<=nx<33 and 0<=ny<33 and state[nx,ny]==1 and rng.random()<p:
                    state[nx,ny]=2
        infected.append(int((state==2).sum()))
        if day in frame_days:snaps[day]=state.copy()
    return snaps,np.array(infected)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',default='data/clima_diario_tres_localidades_qc.csv')
    parser.add_argument('--output',default='results/netlogo_port')
    parser.add_argument('--replicates',type=int,default=100)
    args=parser.parse_args()
    data=pd.read_csv(args.input,parse_dates=['date'])
    data=data[(data.locality=='Socoroma')&data.date.between('2021-01-01','2021-04-10')].sort_values('date')
    if len(data)!=100 or not data[['temp_mean_c','rh_mean_pct']].notna().all().all():
        raise ValueError('Need 100 valid weather days for Socoroma, 2021')
    T=data.temp_mean_c.to_numpy(float)
    RH=data.rh_mean_pct.to_numpy(float)
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    runs=[]
    for k in range(args.replicates):
        snaps,curve=run_once(T,RH,20260903+k)
        runs.append(curve)
    curves=np.vstack(runs)
    df=pd.DataFrame({'day':np.arange(101),'mean_infectious':curves.mean(axis=0),
                     'p025':np.quantile(curves,.025,axis=0),
                     'p975':np.quantile(curves,.975,axis=0)})
    df.to_csv(out/'netlogo_port_trajectories.csv',index=False)

    cm=ListedColormap(['#f2f2f2','#439b57','#bf3a3a'])
    fig,axs=plt.subplots(1,5,figsize=(14,3.4),constrained_layout=True)
    for ax,(day,state) in zip(axs,snaps.items()):
        ax.imshow(state,cmap=cm,vmin=0,vmax=2,interpolation='nearest')
        ax.set_title(f'Day {day}\nI={int((state==2).sum())}')
        ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Exploratory NetLogo-rule port under Socoroma 2021 forcing')
    fig.savefig(out/'spatial_snapshots.png',dpi=230,bbox_inches='tight');plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.8,4.7))
    ax.plot(df.day,df.mean_infectious,lw=1.8)
    ax.fill_between(df.day,df.p025,df.p975,alpha=.15)
    ax.plot(100,df.mean_infectious.iloc[-1],'o')
    ax.set(xlim=(0,100),xlabel='Day (0--100)',ylabel='Infectious patches')
    ax.margins(x=0);ax.grid(alpha=.15);fig.tight_layout()
    fig.savefig(out/'disease_curve.png',dpi=230);plt.close(fig)
    print('Mean infectious patches at day 100:',df.mean_infectious.iloc[-1])

if __name__=='__main__':main()
