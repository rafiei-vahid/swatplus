#!/usr/bin/env python3
"""Daily-resolution comparison: serial (1 core) vs full-wavefront routing, against observations.

Monthly aggregation hides the divergence; this shows it at the timestep the engine
actually computes on, plus the day-by-day scatter that reveals whether the error is
systematic or scattered.

usage: plot_daily_impact.py <runroot> <obs_csv> <out.png>
"""
import sys
from pathlib import Path

import matplotlib
import matplotlib.dates
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, xarray as xr

CFS = 0.0283168
OBS_C, SER_C, WAV_C = "#52514e", "#2a78d6", "#eb6834"
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


def nse(o, s):
    j = pd.concat([o, s], axis=1, join="inner").dropna()
    ov, sv = j.iloc[:, 0].values, j.iloc[:, 1].values
    return 1 - ((ov - sv) ** 2).sum() / ((ov - ov.mean()) ** 2).sum(), len(j)


def main():
    root, obs_csv, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    S, W = load(root / "serial"), load(root / "wave")
    o = pd.read_csv(obs_csv, parse_dates=["date"])
    o.loc[o["streamflow"] <= -90, "streamflow"] = np.nan
    obs = (o.set_index("date")["streamflow"] * CFS).dropna()

    fig, axes = plt.subplots(3, 2, figsize=(13.5, 11.0), facecolor="#fcfcfb")
    fig.subplots_adjust(hspace=0.45, wspace=0.22, top=0.955, bottom=0.075,
                        left=0.075, right=0.985)

    # ---- (a) full daily record, log scale ---------------------------------
    ax = axes[0, 0]; style(ax)
    ax.set_yscale("log")
    win = S["flo_out"].index
    ob = obs.reindex(win)
    ax.plot(win, ob.values, color=OBS_C, lw=1.0, label="Observed", zorder=3)
    ax.plot(win, S["flo_out"].values, color=SER_C, lw=1.0, label="Serial (1 core)", zorder=4)
    ax.plot(win, W["flo_out"].values, color=WAV_C, lw=1.0, alpha=0.85,
            label="Full wavefront", zorder=5)
    ns, n = nse(obs, S["flo_out"]); nw, _ = nse(obs, W["flo_out"])
    ax.set_ylabel("Daily flow  (m³/s, log)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="lower right", ncol=3)
    ax.text(0.0, 1.06, f"a · Full daily record — daily NSE serial {ns:.3f} vs wavefront {nw:.3f} "
                       f"(n={n})", transform=ax.transAxes, color=INK, fontsize=9.5,
            fontweight="bold")

    # ---- (b) the largest event, linear ------------------------------------
    ax = axes[0, 1]; style(ax)
    peak = S["flo_out"].idxmax()
    sl = slice(peak - pd.Timedelta(days=20), peak + pd.Timedelta(days=20))
    ax.plot(obs.loc[sl].index, obs.loc[sl].values, color=OBS_C, lw=2.0, marker="o",
            ms=3.5, label="Observed", zorder=3)
    ax.plot(S["flo_out"].loc[sl].index, S["flo_out"].loc[sl].values, color=SER_C, lw=2.0,
            marker="o", ms=3.5, label="Serial", zorder=4)
    ax.plot(W["flo_out"].loc[sl].index, W["flo_out"].loc[sl].values, color=WAV_C, lw=2.0,
            ls=(0, (4, 2)), marker="s", ms=3.5, label="Full wavefront", zorder=5)
    ax.set_ylabel("Daily flow  (m³/s)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="upper left")
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
    d = 100 * (W["flo_out"].loc[peak] - S["flo_out"].loc[peak]) / S["flo_out"].loc[peak]
    ax.text(0.0, 1.06, f"b · The largest event ±20 d — peak day differs by {d:+.1f}%",
            transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    # ---- (c,d) daily constituent series over one high-load year -----------
    yr = S["orgn_out"].resample("YS").sum().idxmax().year
    ysl = slice(f"{yr}-01-01", f"{yr}-12-31")
    for (key, label), ax, tag in ((("orgn_out", "Organic N"), axes[1, 0], "c"),
                                  (("sedp_out", "Sediment P"), axes[1, 1], "d")):
        style(ax)
        ax.plot(S[key].loc[ysl].index, S[key].loc[ysl].values, color=SER_C, lw=1.5,
                label="Serial (1 core)", zorder=4)
        ax.plot(W[key].loc[ysl].index, W[key].loc[ysl].values, color=WAV_C, lw=1.5,
                ls=(0, (4, 2)), label="Full wavefront", zorder=5)
        ax.set_ylabel(f"{label}  (kg/day)", color=INK, fontsize=9)
        ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="upper right")
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
        pk = 100 * (W[key].loc[ysl].max() - S[key].loc[ysl].max()) / S[key].loc[ysl].max()
        ax.text(0.0, 1.06, f"{tag} · {label} daily, {yr} — annual peak day {pk:+.1f}%",
                transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    # ---- (e,f) day-by-day scatter, wavefront vs serial ---------------------
    for (key, label, unit), ax, tag in ((("flo_out", "Flow", "m³/s"), axes[2, 0], "e"),
                                        (("orgn_out", "Organic N", "kg/d"), axes[2, 1], "f")):
        style(ax)
        s, w = S[key].values, W[key].values
        m = s > (np.nanmax(s) * 1e-4)
        lim = [np.nanmin(s[m]) * 0.8, np.nanmax(s) * 1.3]
        ax.plot(lim, lim, color=MUTED, lw=1.2, zorder=2)
        ax.scatter(s[m], w[m], s=9, color=WAV_C, alpha=0.45, lw=0, zorder=4)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(f"Serial  ({unit})", color=INK, fontsize=9)
        ax.set_ylabel(f"Full wavefront  ({unit})", color=INK, fontsize=9)
        off = 100 * np.median(np.abs(w[m] - s[m]) / s[m])
        ax.text(0.0, 1.06, f"{tag} · {label}, every day — median |error| {off:.1f}% "
                           f"(line is 1:1)", transform=ax.transAxes, color=INK,
                fontsize=9.5, fontweight="bold")

    fig.savefig(out, dpi=200, facecolor="#fcfcfb", bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
