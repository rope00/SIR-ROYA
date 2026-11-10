#!/usr/bin/env python3
import calendar, time, re
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup

BASE = "https://climatologia.meteochile.gob.cl/application/mensual/climatMensualDatosDiarios/{station}/{year}/{month}"
STATIONS = {"180030":"Socoroma","180038":"Ticnamar","180046":"Chapiquina"}
YEARS = range(2020, 2025)
MONTHS = range(1, 13)

RAW = Path("re1018_raw_html")
RAW.mkdir(exist_ok=True)

session = requests.Session()
session.headers["User-Agent"] = "Mozilla/5.0"

def txt(x):
    return re.sub(r"\s+", " ", str(x)).strip()

def num(x):
    s = txt(x)
    if s in {"", ".", "-", "--"}:
        return np.nan
    try:
        return float(s.replace(",", ""))
    except Exception:
        return np.nan

def precip(x):
    s = txt(x).lower()
    if s in {"s/p", "sp"}:
        return 0.0
    return num(s)

def fetch(station, year, month):
    p = RAW / f"{station}_{year}_{month:02d}.html"
    url = BASE.format(station=station, year=year, month=month)
    if p.exists() and p.stat().st_size > 1000:
        return p.read_text(encoding="utf-8", errors="replace"), url

    r = session.get(url, timeout=60)
    r.raise_for_status()
    if len(r.text) < 1000:
        raise RuntimeError("Respuesta HTML demasiado corta")
    p.write_text(r.text, encoding="utf-8")
    return r.text, url

def parse_month(html, station, locality, year, month, url):
    soup = BeautifulSoup(html, "html.parser")
    table = None
    for t in soup.find_all("table"):
        s = txt(t.get_text(" ", strip=True)).lower()
        if "temperatura" in s and "humedad media" in s:
            table = t
            break
    if table is None:
        return pd.DataFrame()

    rows = []
    max_day = calendar.monthrange(year, month)[1]

    for tr in table.find_all("tr"):
        c = [txt(x.get_text(" ", strip=True)) for x in tr.find_all(["th","td"])]
        if len(c) != 10:
            continue
        try:
            day = int(c[0])
        except Exception:
            continue
        if not (1 <= day <= max_day):
            continue

        rows.append({
            "date": pd.Timestamp(year=year, month=month, day=day),
            "station_code": station,
            "locality": locality,
            "temp_min_c": num(c[1]),
            "temp_max_c": num(c[2]),
            "temp_mean_c": num(c[3]),
            "rh_mean_pct": num(c[4]),
            "precip_mm": precip(c[5]),
            "evapotranspiration_mm_day": num(c[6]),
            "radiation_total_daily_raw": num(c[7]),
            "radiation_max_w_m2": num(c[8]),
            "radiation_max_time_local": c[9] if c[9] not in {"", "."} else None,
            "source": "DMC RE1018",
            "source_url": url,
        })
    return pd.DataFrame(rows)

frames = []
failures = []

total = len(STATIONS) * len(list(YEARS)) * len(list(MONTHS))
i = 0

print(f"Consultas programadas: {total}")

for station, locality in STATIONS.items():
    for year in YEARS:
        for month in MONTHS:
            i += 1
            label = f"{locality} {year}-{month:02d}"
            print(f"[{i:03d}/{total}] {label}", end=" ... ")

            try:
                html, url = fetch(station, year, month)
                df = parse_month(html, station, locality, year, month, url)

                if df.empty:
                    print("SIN DATOS")
                    failures.append([station, locality, year, month, "sin datos"])
                else:
                    frames.append(df)
                    nv_t = df["temp_mean_c"].notna().sum()
                    nv_h = df["rh_mean_pct"].notna().sum()
                    print(f"OK ({len(df)} dias; T={nv_t}; HR={nv_h})")
            except Exception as e:
                print(f"ERROR: {e}")
                failures.append([station, locality, year, month, str(e)])

            time.sleep(0.35)

if not frames:
    raise SystemExit("No se recuperaron datos.")

full = pd.concat(frames, ignore_index=True)
full = full.sort_values(["station_code","date"]).reset_index(drop=True)

if full.duplicated(["station_code","date"]).any():
    raise SystemExit("Hay fechas duplicadas. Se detuvo por seguridad.")

full.to_csv("clima_diario_tres_localidades_completo.csv", index=False, date_format="%Y-%m-%d")

model = full[[
    "date","locality","temp_mean_c","rh_mean_pct",
    "precip_mm","evapotranspiration_mm_day","source"
]].copy()

model.to_csv("clima_diario_tres_localidades.csv", index=False, date_format="%Y-%m-%d")

expected = pd.date_range("2020-01-01","2024-12-31",freq="D")
report = []

for station, locality in STATIONS.items():
    x = full[full["station_code"] == station]
    valid = x.dropna(subset=["temp_mean_c","rh_mean_pct"])
    dates = pd.DatetimeIndex(valid["date"].unique())

    report.append({
        "station_code": station,
        "locality": locality,
        "valid_days": len(dates),
        "expected_days": len(expected),
        "coverage_pct": round(100*len(dates.intersection(expected))/len(expected),2),
        "missing_days": len(expected.difference(dates)),
        "mean_temp_c": round(valid["temp_mean_c"].mean(),3) if len(valid) else np.nan,
        "mean_rh_pct": round(valid["rh_mean_pct"].mean(),3) if len(valid) else np.nan,
    })

report = pd.DataFrame(report)
report.to_csv("reporte_cobertura_re1018.csv", index=False)

if failures:
    pd.DataFrame(
        failures,
        columns=["station_code","locality","year","month","reason"]
    ).to_csv("re1018_fallos.csv", index=False)

print("\n=== COBERTURA ===")
print(report.to_string(index=False))
print("\nGenerados:")
print("clima_diario_tres_localidades.csv")
print("clima_diario_tres_localidades_completo.csv")
print("reporte_cobertura_re1018.csv")
if failures:
    print("re1018_fallos.csv")
