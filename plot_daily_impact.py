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
# validated 4-series palette (dataviz six checks, light mode): gray / blue / aqua / orange
OBS_C, SER_C, HRU_C, WAV_C = "#52514e", "#2a78d6", "#1baf7a", "#eb6834"
L_SER = "Serial (1 core)"
L_HRU = "HRU-parallel · serial routing (shipped)"
L_WAV = "Full wavefront · parallel routing"
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
    S, H, W = load(root / "serial"), load(root / "rtser"), load(root / "wave")
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
    ax.plot(win, ob.values, color=OBS_C, lw=0.9, label="Observed", zorder=3)
    ax.plot(win, S["flo_out"].values, color=SER_C, lw=1.3, label=L_SER, zorder=4)
    ax.plot(win, H["flo_out"].values, color=HRU_C, lw=0.8, ls=(0, (2, 2)),
            label=L_HRU, zorder=6)
    ax.plot(win, W["flo_out"].values, color=WAV_C, lw=1.0, alpha=0.9,
            label=L_WAV, zorder=5)
    ns, n = nse(obs, S["flo_out"]); nw, _ = nse(obs, W["flo_out"])
    nh, _ = nse(obs, H["flo_out"])
    ax.set_ylabel("Daily flow  (m³/s, log)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=MUTED, loc="lower right", ncol=2)
    ax.text(0.0, 1.06, f"a · Daily NSE — serial {ns:.3f} · HRU-parallel {nh:.3f} · "
                       f"wavefront {nw:.3f}  (n={n})", transform=ax.transAxes, color=INK,
            fontsize=9.5, fontweight="bold")

    # ---- (b) the largest event, linear ------------------------------------
    ax = axes[0, 1]; style(ax)
    peak = S["flo_out"].idxmax()
    sl = slice(peak - pd.Timedelta(days=20), peak + pd.Timedelta(days=20))
    ax.plot(obs.loc[sl].index, obs.loc[sl].values, color=OBS_C, lw=2.0, marker="o",
            ms=3.5, label="Observed", zorder=3)
    ax.plot(S["flo_out"].loc[sl].index, S["flo_out"].loc[sl].values, color=SER_C, lw=2.4,
            marker="o", ms=3.8, label=L_SER, zorder=4)
    # Identity must READ as identity. A thin dashed line laid over a thick one looks like
    # "similar", and next to a second dashed line the eye groups the two dashed series
    # together — the exact misreading this panel caused. Spaced ring markers riding the
    # serial line cannot be mistaken for a separate trajectory.
    _h = H["flo_out"].loc[sl]
    ax.plot(_h.index[::2], _h.values[::2], color=HRU_C, lw=0, marker="o", ms=7,
            markerfacecolor="none", markeredgewidth=1.6, label=L_HRU + " — exactly on serial",
            zorder=6)
    ax.plot(W["flo_out"].loc[sl].index, W["flo_out"].loc[sl].values, color=WAV_C, lw=2.0,
            ls=(0, (5, 2)), marker="s", ms=3.5, label=L_WAV, zorder=5)
    ax.set_ylabel("Daily flow  (m³/s)", color=INK, fontsize=9)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=MUTED, loc="upper left")
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
        ax.plot(S[key].loc[ysl].index, S[key].loc[ysl].values, color=SER_C, lw=2.0,
                label=L_SER, zorder=4)
        ax.plot(H[key].loc[ysl].index, H[key].loc[ysl].values, color=HRU_C, lw=1.1,
                ls=(0, (2, 2)), label=L_HRU, zorder=6)
        ax.plot(W[key].loc[ysl].index, W[key].loc[ysl].values, color=WAV_C, lw=1.5,
                ls=(0, (5, 2)), label=L_WAV, zorder=5)
        ax.set_ylabel(f"{label}  (kg/day)", color=INK, fontsize=9)
        ax.legend(frameon=False, fontsize=7.5, labelcolor=MUTED, loc="upper right")
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
        pk = 100 * (W[key].loc[ysl].max() - S[key].loc[ysl].max()) / S[key].loc[ysl].max()
        ax.text(0.0, 1.06, f"{tag} · {label} daily, {yr} — annual peak day {pk:+.1f}%",
                transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    # ---- (e,f) day-by-day scatter, wavefront vs serial ---------------------
    for (key, label, unit), ax, tag in ((("flo_out", "Flow", "m³/s"), axes[2, 0], "e"),
                                        (("orgn_out", "Organic N", "kg/d"), axes[2, 1], "f")):
        style(ax)
        s, w, h = S[key].values, W[key].values, H[key].values
        m = s > (np.nanmax(s) * 1e-4)
        lim = [np.nanmin(s[m]) * 0.8, np.nanmax(s) * 1.3]
        ax.plot(lim, lim, color=MUTED, lw=1.2, zorder=2)
        ax.scatter(s[m], w[m], s=11, color=WAV_C, alpha=0.40, lw=0, zorder=4,
                   label=L_WAV)
        ax.scatter(s[m], h[m], s=7, color=HRU_C, alpha=0.75, lw=0, zorder=5,
                   label=L_HRU)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(f"Serial  ({unit})", color=INK, fontsize=9)
        ax.set_ylabel(f"Parallel mode  ({unit})", color=INK, fontsize=9)
        off = 100 * np.median(np.abs(w[m] - s[m]) / s[m])
        offh = 100 * np.median(np.abs(h[m] - s[m]) / s[m])
        ax.legend(frameon=False, fontsize=7.5, labelcolor=MUTED, loc="upper left",
                  markerscale=1.6)
        ax.text(0.0, 1.06, f"{tag} · {label}, every day — median |error| "
                           f"{offh:.2f}% vs {off:.1f}%  (line is 1:1)",
                transform=ax.transAxes, color=INK, fontsize=9.5, fontweight="bold")

    fig.savefig(out, dpi=200, facecolor="#fcfcfb", bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
