#!/usr/bin/env python3
"""Perform the documented, conservative QC of DMC RE1018 daily records.

No interpolation and no replacement weather values.
Use after code/download_re1018.py completes.
"""
import argparse
from pathlib import Path
import pandas as pd

def clean(inp, out):
    d=pd.read_csv(inp,parse_dates=['date'])
    required={'date','locality','temp_mean_c','rh_mean_pct','precip_mm'}
    if not required.issubset(d.columns):
        raise ValueError(f"Required columns missing: {required-set(d.columns)}")
    d['qc_flag']=''
    invalid=(d.locality=='Ticnamar') & (d.date.between('2020-02-22','2020-02-29')) & (d.temp_mean_c <= -20)
    if invalid.sum()!=8:
        raise ValueError(f"Expected eight documented Ticnamar anomalies, found {invalid.sum()}")
    d.loc[invalid,'qc_flag']='sentinel_temperature_around_minus40'
    d.loc[invalid,['temp_mean_c','rh_mean_pct']]=float('nan')
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(out,index=False,date_format='%Y-%m-%d')
    return len(d),int(invalid.sum())

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument('--input',default='clima_diario_tres_localidades.csv')
    p.add_argument('--output',default='data/clima_diario_tres_localidades_qc.csv')
    args=p.parse_args()
    print("Rows, flagged days:",clean(args.input,args.output))
