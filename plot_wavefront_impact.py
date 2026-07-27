#!/usr/bin/env python3
"""Sim-vs-obs and sim-vs-sim panels for serial (1 core) vs full-wavefront routing.

Makes the measured 4-9 % constituent deviation visible: whether it is a uniform
offset or concentrated in events, and what it does to the hydrograph against
observations.

usage: plot_wavefront_impact.py <runroot> <obs_csv> <out.png>
"""
import sys
from pathlib import Path

import matplotlib
import matplotlib.dates
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, xarray as xr
from matplotlib.ticker import FuncFormatter

CFS = 0.0283168
OBS_C, SER_C, WAV_C = "#52514e", "#2a78d6", "#eb6834"   # validated: gray / blue / orange
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d8d7d2"

_BASE = ["flo", "sed", "orgn", "sedp", "no3", "solp", "chla", "nh3", "no2",
         "cbod", "dox", "san", "sil", "cla", "sag", "lag", "grv", "temp"]
SD = (["area", "precip", "evap", "seep"] + [f"{c}_stor" for c in _BASE]
      + [f"{c}_in" for c in _BASE] + [f"{c}_out" for c in _BASE])


def load(run_dir):
    with xr.open_dataset(Path(run_dir) / "channel_sd_day.nc") as ds:
        idx = pd.to_datetime(dict(year=ds["yrc"].values.astype(int),
                                  month=ds["mo"].values.astype(int),
                                  day=ds["day_mo"].values.astype(int)))
        out = {}
        for v in ds.variables:
            if not v.startswith("v"):
                continue
            i = int(v[1:]) - 1
            if i < len(SD):
                a = np.asarray(ds[v].values, dtype=np.float64)
                out[SD[i]] = pd.Series(a[:, 0] if a.ndim == 2 else a, index=idx)
    return out


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


def main():
    root, obs_csv, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    S, W = load(root / "serial"), load(root / "wave")

    o = pd.read_csv(obs_csv, parse_dates=["date"])
    o.loc[o["streamflow"] <= -90, "streamflow"] = np.nan
    oq = o.set_index("date")["streamflow"] * CFS
    g = oq.resample("MS")
    obs_m = g.mean()[(g.count() / g.size()) >= 0.90].dropna()

    fig, axes = plt.subplots(3, 2, figsize=(13.5, 11.0), facecolor="#fcfcfb")
    fig.subplots_adjust(hspace=0.42, wspace=0.22, top=0.955, bottom=0.075,
                        left=0.075, right=0.985)

    # ---- (a) monthly flow vs observations --------------------------------
    ax = axes[0, 0]; style(ax)
    sm, wm = S["flo_out"].resample("MS").mean(), W["flo_out"].resample("MS").mean()
    win = obs_m.index.intersection(sm.index)
    ax.plot(win, obs_m.loc[win], color=OBS_C, lw=2.4, label="Observed (USGS)", zorder=3)
    ax.plot(win, sm.loc[win], color=SER_C, lw=2.0, label="Serial (1 core)", zorder=4)
    ax.plot(win, wm.loc[win], color=WAV_C, lw=2.0, ls=(0, (4, 2)),
            label="Full wavefront (8 threads)", zorder=5)
    ax.set_ylabel("Monthly mean flow  (m³/s)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="upper left")
    ax.text(0.0, 1.06, "a · Streamflow against observations — all three overlap",
            transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    # ---- (b) daily flow around the largest event, fully inside the run ----
    ax = axes[0, 1]; style(ax)
    flo = S["flo_out"]
    peak = flo.idxmax()
    lo = max(flo.index.min(), peak - pd.Timedelta(days=55))
    hi = min(flo.index.max(), peak + pd.Timedelta(days=55))
    sl = slice(lo, hi)
    d_obs = oq.loc[sl]
    ax.plot(d_obs.index, d_obs.values, color=OBS_C, lw=1.7, label="Observed", zorder=3)
    ax.plot(flo.loc[sl].index, flo.loc[sl].values, color=SER_C, lw=1.7,
            label="Serial", zorder=4)
    ax.plot(W["flo_out"].loc[sl].index, W["flo_out"].loc[sl].values, color=WAV_C, lw=1.7,
            ls=(0, (4, 2)), label="Full wavefront", zorder=5)
    ax.set_ylabel("Daily flow  (m³/s)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="upper left")
    ax.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
    ax.text(0.0, 1.06, f"b · Daily hydrograph around the largest event "
                       f"({lo:%b %Y} – {hi:%b %Y})",
            transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    # ---- (c-e) the three particulate constituents ------------------------
    panels = [("orgn_out", "Organic N", axes[1, 0], "c"),
              ("sedp_out", "Sediment P", axes[1, 1], "d"),
              ("nh3_out", "Ammonia", axes[2, 0], "e")]
    for key, label, ax, tag in panels:
        style(ax)
        s, w = S[key].resample("MS").sum(), W[key].resample("MS").sum()
        ax.plot(s.index, s.values, color=SER_C, lw=2.0, label="Serial (1 core)", zorder=4)
        ax.plot(w.index, w.values, color=WAV_C, lw=2.0, ls=(0, (4, 2)),
                label="Full wavefront", zorder=5)
        ax.fill_between(s.index, s.values, w.values, color=WAV_C, alpha=0.16,
                        lw=0, zorder=2)
        tot = 100 * (w.sum() - s.sum()) / s.sum()
        ax.set_ylabel(f"{label} load  (kg/month)", color=INK, fontsize=9)
        ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="upper left")
        ax.text(0.0, 1.06, f"{tag} · {label} — wavefront {tot:+.1f}% over the run",
                transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    # ---- (f) per-month relative difference --------------------------------
    ax = axes[2, 1]; style(ax)
    ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
    series = [("flo_out", "Flow", "#2a78d6", "-"), ("orgn_out", "Organic N", "#eb6834", "-"),
              ("sedp_out", "Sediment P", "#1baf7a", (0, (4, 2))),
              ("nh3_out", "Ammonia", "#eda100", (0, (1, 1.6)))]
    for key, label, c, ls in series:
        agg = "mean" if key == "flo_out" else "sum"
        s = getattr(S[key].resample("MS"), agg)()
        w = getattr(W[key].resample("MS"), agg)()
        m = np.abs(s.values) > 1e-9
        ax.plot(s.index[m], 100 * (w.values[m] - s.values[m]) / s.values[m],
                color=c, lw=1.9, ls=ls, label=label, zorder=4)
    ax.set_ylabel("Wavefront − serial  (% of serial)", color=INK, fontsize=9)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:+.0f}%"))
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, ncol=2, loc="lower left")
    ax.text(0.0, 1.06, "f · Month-by-month divergence — event-driven, not a constant offset",
            transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    fig.savefig(out, dpi=200, facecolor="#fcfcfb", bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
