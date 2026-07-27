#!/usr/bin/env python3
"""Observed vs simulated streamflow scatter, one panel per engine mode, with skill metrics.

Top row daily, bottom row monthly. Each panel: observed on x, simulated on y, 1:1
line, and NSE / KGE / PBIAS for that mode.

usage: plot_obs_scatter.py <runroot> <obs_csv> <out.png>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, xarray as xr

CFS = 0.0283168
SER_C, HRU_C, WAV_C = "#2a78d6", "#1baf7a", "#eb6834"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d8d7d2"

_BASE = ["flo", "sed", "orgn", "sedp", "no3", "solp", "chla", "nh3", "no2",
         "cbod", "dox", "san", "sil", "cla", "sag", "lag", "grv", "temp"]
SD = (["area", "precip", "evap", "seep"] + [f"{c}_stor" for c in _BASE]
      + [f"{c}_in" for c in _BASE] + [f"{c}_out" for c in _BASE])

MODES = [("serial", "Serial (1 core)", SER_C),
         ("rtser", "HRU-parallel · serial routing (shipped)", HRU_C),
         ("wave", "Full wavefront · parallel routing", WAV_C)]


def flow(run_dir):
    with xr.open_dataset(Path(run_dir) / "channel_sd_day.nc") as ds:
        idx = pd.to_datetime(dict(year=ds["yrc"].values.astype(int),
                                  month=ds["mo"].values.astype(int),
                                  day=ds["day_mo"].values.astype(int)))
        i = SD.index("flo_out")
        a = np.asarray(ds[f"v{i+1}"].values, dtype=np.float64)
    return pd.Series(a[:, 0] if a.ndim == 2 else a, index=idx)


def style(ax):
    ax.set_facecolor("#fcfcfb")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=8, length=3)
    for lb in ax.get_xticklabels() + ax.get_yticklabels():
        lb.set_color(MUTED)


def skill(o, s):
    j = pd.concat([o, s], axis=1, join="inner").dropna()
    ov, sv = j.iloc[:, 0].values, j.iloc[:, 1].values
    nse = 1 - ((ov - sv) ** 2).sum() / ((ov - ov.mean()) ** 2).sum()
    pb = 100 * (ov - sv).sum() / ov.sum()
    r = np.corrcoef(ov, sv)[0, 1]
    kge = 1 - np.sqrt((r - 1) ** 2 + (sv.std() / ov.std() - 1) ** 2
                      + (sv.mean() / ov.mean() - 1) ** 2)
    return ov, sv, dict(n=len(j), nse=nse, kge=kge, pbias=pb, r2=r ** 2)


def main():
    root, obs_csv, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    sims = {k: flow(root / k) for k, _, _ in MODES}

    o = pd.read_csv(obs_csv, parse_dates=["date"])
    o.loc[o["streamflow"] <= -90, "streamflow"] = np.nan
    obs_d = (o.set_index("date")["streamflow"] * CFS).dropna()
    g = obs_d.resample("MS")
    obs_m = g.mean()[(g.count() / g.size()) >= 0.90].dropna()

    fig, axes = plt.subplots(2, 3, figsize=(14.0, 9.4), facecolor="#fcfcfb")
    fig.subplots_adjust(hspace=0.34, wspace=0.24, top=0.93, bottom=0.08,
                        left=0.065, right=0.985)

    for row, (scale, obs, lab) in enumerate(
            [("daily", obs_d, "Daily"), ("monthly", obs_m, "Monthly mean")]):
        for col, (key, name, colour) in enumerate(MODES):
            ax = axes[row, col]; style(ax)
            s = sims[key] if scale == "daily" else sims[key].resample("MS").mean()
            ov, sv, m = skill(obs, s)
            pos = (ov > 0) & (sv > 0)
            # clip the lower bound to a low percentile: a single near-zero day
            # otherwise stretches four empty decades and shrinks the cloud
            floor = min(np.percentile(ov[pos], 0.5), np.percentile(sv[pos], 0.5))
            lim = [floor * 0.6, max(ov.max(), sv.max()) * 1.4]
            ax.plot(lim, lim, color=MUTED, lw=1.2, zorder=2)
            ax.scatter(ov[pos], sv[pos], s=13 if scale == "daily" else 30,
                       color=colour, alpha=0.42 if scale == "daily" else 0.8,
                       lw=0, zorder=4)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(lim); ax.set_ylim(lim)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel(f"Observed  (m³/s)", color=INK, fontsize=9)
            if col == 0:
                ax.set_ylabel(f"{lab} simulated  (m³/s)", color=INK, fontsize=9)
            ax.text(0.04, 0.96,
                    f"NSE  {m['nse']:.3f}\nKGE  {m['kge']:.3f}\n"
                    f"PBIAS {m['pbias']:+.1f}%\nR²   {m['r2']:.3f}\nn = {m['n']}",
                    transform=ax.transAxes, va="top", ha="left",
                    color=INK, fontsize=8.5, linespacing=1.5,
                    bbox=dict(boxstyle="round,pad=0.45", fc="#fcfcfb",
                              ec=GRID, lw=0.8))
            ax.text(0.0, 1.045, f"{'abcdef'[row*3+col]} · {lab} — {name}",
                    transform=ax.transAxes, color=INK, fontsize=8.6,
                    fontweight="bold")

    fig.savefig(out, dpi=200, facecolor="#fcfcfb", bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
