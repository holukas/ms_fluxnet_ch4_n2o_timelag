import logging
import warnings

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.font_manager as _fm
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import gaussian_kde

# ══════════════════════════════════════════════════════════════════════════════
#  CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
FILE1 = {
    "label": "CM [0–5 s, |cov| max]",
    "path": "OPENLAG-5s_2021_2_LGR_eddypro_CH-CHA_FR-20240725-121331_fluxnet_2024-07-25T121333_adv.csv",
}

FILE2 = {
    "label": "PWB (Vitale et al., 2024)",
    "paths": [
        "../../output/tlag_results_part1.csv",
        "../../output/tlag_results_part2.csv",
        "../../output/tlag_results_part3.csv",
        "../../output/tlag_results_part4.csv",
        "../../output/tlag_results_part5.csv",
        "../../output/tlag_results_part6.csv",
        "../../output/tlag_results_part7.csv",
    ],
    "ts_col": 0,
    "ts_pattern": r"(\d{8}-\d{4})",
    "ts_format": "%Y%m%d-%H%M",
    "cols": {
        "CH4_TLAG_ACTUAL": "ch4_tlag_sec",
        "N2O_TLAG_ACTUAL": "n2o_tlag_sec",
    },
}

SAVE_PATH = "timelags.png"
JITTER = 0.06
TILLAGE_DATE = "2021-08-20"
FERTILIZATION_DATE = "2021-07-29"
PRECIPITATION_DATE = "2021-08-16"
YLIM = (0, 5.1)
DPI = 300

COLOR_F1 = "#0072B2"  # Wong blue  — CM
COLOR_F2 = "#E05C2A"  # coral-orange — PWB
# ══════════════════════════════════════════════════════════════════════════════

# ── font ──────────────────────────────────────────────────────────────────────
_PREFERRED = ["Helvetica Neue", "Helvetica", "Arial", "Calibri",
              "Segoe UI", "Gill Sans MT", "Trebuchet MS",
              "Liberation Sans", "DejaVu Sans"]
_available = {f.name for f in _fm.fontManager.ttflist}
_FONT = next((f for f in _PREFERRED if f in _available), "sans-serif")

# ── rcParams ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": _FONT,
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "axes.edgecolor": "black",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#DDDDDD",
    "grid.linewidth": 0.4,
    "grid.linestyle": "--",
    "axes.grid.axis": "y",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.pad": 3,
    "ytick.major.pad": 3,
    "lines.linewidth": 1.0,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
})


# ── load ──────────────────────────────────────────────────────────────────────
def load_file1(cfg):
    df = pd.read_csv(cfg["path"], na_values=[-9999],
                     usecols=["TIMESTAMP_END", "CH4_TLAG_ACTUAL", "N2O_TLAG_ACTUAL"])
    df["timestamp"] = pd.to_datetime(df["TIMESTAMP_END"], format="%Y%m%d%H%M")
    return df.sort_values("timestamp").reset_index(drop=True)


def load_file2(cfg):
    frames = []
    for path in cfg["paths"]:
        df = pd.read_csv(path, na_values=[-9999])
        fn_col = df.columns[cfg["ts_col"]]
        df["timestamp"] = (
            df[fn_col]
            .str.extract(cfg["ts_pattern"])[0]
            .apply(lambda s: pd.to_datetime(s, format=cfg["ts_format"])
            if pd.notna(s) else pd.NaT)
        )
        df = df.rename(columns={v: k for k, v in cfg["cols"].items()})
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True).drop_duplicates(subset="timestamp")
    return merged.sort_values("timestamp").reset_index(drop=True)


ds1 = load_file1(FILE1)
ds2 = load_file2(FILE2)


# ── KDE helper ────────────────────────────────────────────────────────────────
def smooth_kde(series, n=500):
    lo, hi = series.min(), series.max()
    pad = (hi - lo) * 0.12
    y = np.linspace(lo - pad, hi + pad, n)
    return y, gaussian_kde(series, bw_method="scott")(y)


rng = np.random.default_rng(42)

# ── global x-limits ───────────────────────────────────────────────────────────
_all_ts = list(ds1["timestamp"].dropna()) + list(ds2["timestamp"].dropna())
x_min, x_max = min(_all_ts), max(_all_ts)

# ══════════════════════════════════════════════════════════════════════════════
#  LAYOUT  — 2 rows × 2 cols, each cell = [scatter | KDE]
#  Implemented as a 2×5 GridSpec: [scatter | kde | gap | scatter | kde]
# ══════════════════════════════════════════════════════════════════════════════
FIG_W = 28 / 2.54
FIG_H = 14 / 2.54

fig = plt.figure(figsize=(FIG_W, FIG_H))

gs = gridspec.GridSpec(
    2, 5, figure=fig,
    hspace=0.14, wspace=0.06,
    left=0.04, right=0.98,
    top=0.91, bottom=0.10,
    width_ratios=[5, 1, 0.25, 5, 1],
)

ax_master = fig.add_subplot(gs[0, 0])

axes = {
    ("CH4", "f1"): (ax_master, fig.add_subplot(gs[0, 1], sharey=ax_master)),
    ("CH4", "f2"): (fig.add_subplot(gs[0, 3], sharex=ax_master, sharey=ax_master),
                    fig.add_subplot(gs[0, 4], sharey=ax_master)),
    ("N2O", "f1"): (fig.add_subplot(gs[1, 0], sharex=ax_master, sharey=ax_master),
                    fig.add_subplot(gs[1, 1], sharey=ax_master)),
    ("N2O", "f2"): (fig.add_subplot(gs[1, 3], sharex=ax_master, sharey=ax_master),
                    fig.add_subplot(gs[1, 4], sharey=ax_master)),
}

# ── panel definitions ─────────────────────────────────────────────────────────
PANEL_DEFS = [
    ("CH4", "f1", ds1, COLOR_F1, "a", "CH$_4$ \u2014 " + FILE1["label"], True,  False, True),
    ("CH4", "f2", ds2, COLOR_F2, "b", "CH$_4$ \u2014 " + FILE2["label"], False, False, True),
    ("N2O", "f1", ds1, COLOR_F1, "c", "N$_2$O \u2014 " + FILE1["label"], True,  True,  False),
    ("N2O", "f2", ds2, COLOR_F2, "d", "N$_2$O \u2014 " + FILE2["label"], False, True,  False),
]

# ── event lines ───────────────────────────────────────────────────────────────
EVENTS = [
    (TILLAGE_DATE,       "Tillage",       "-",  "black", "right"),
    (FERTILIZATION_DATE, "Fertilization", ":",  "black", "right"),
    (PRECIPITATION_DATE, "Precipitation", "--", "black", "left"),
]


def draw_events(ax, show_labels):
    for date_str, label, ls, color, side in EVENTS:
        dt = pd.Timestamp(date_str)
        ax.axvline(dt, color=color, lw=1.0, ls=ls, zorder=6)
        if show_labels:
            offset = pd.Timedelta(days=1) if side == "right" else pd.Timedelta(days=-1)
            ha = "left" if side == "right" else "right"
            ax.text(dt + offset, 0.98, label,
                    fontsize=9, va="top", ha=ha, color=color,
                    transform=ax.get_xaxis_transform(),
                    rotation=90, clip_on=True)


# ── main plot loop ────────────────────────────────────────────────────────────
for gas, fkey, ds, color, letter, panel_label, show_ylabel, show_xlabel, show_event_lbl in PANEL_DEFS:
    col = "CH4_TLAG_ACTUAL" if gas == "CH4" else "N2O_TLAG_ACTUAL"
    ax_s, ax_k = axes[(gas, fkey)]

    series = ds[col].dropna()
    ts = ds.loc[series.index, "timestamp"]

    # scatter
    jittered = series.values + rng.uniform(-JITTER, JITTER, size=len(series))
    ax_s.scatter(ts, jittered, marker="o", color=color,
                 s=6, alpha=0.20, linewidths=0, zorder=3)

    # mode line (black dashed, same for both methods)
    if len(series) > 5:
        y_m, d_m = smooth_kde(series)
        mode_val = round(y_m[np.argmax(d_m)] / 0.05) * 0.05
        ax_s.axhline(mode_val, color="black", lw=1.2, ls="--", zorder=5)
        ax_k.axhline(mode_val, color="black", lw=1.0, ls="--", zorder=5)

    # KDE
    if len(series) > 5:
        y_g, dens = smooth_kde(series)
        dn = dens / dens.max()
        ax_k.fill_betweenx(y_g, dn, color=color, alpha=0.15)
        ax_k.plot(dn, y_g, color=color, lw=1.2, alpha=1.0)

    # event lines
    draw_events(ax_s, show_labels=show_event_lbl)

    # scatter axis styling
    ax_s.set_ylim(YLIM)
    ax_s.set_xlim(x_min, x_max)
    ax_s.xaxis.set_major_locator(mdates.MonthLocator())
    ax_s.xaxis.set_major_formatter(mdates.DateFormatter("%b"))

    if not show_xlabel:
        ax_s.tick_params(axis="x", labelbottom=False)
    else:
        fig.canvas.draw()
        for loc, lbl in zip(ax_s.xaxis.get_majorticklocs(),
                            ax_s.get_xticklabels()):
            if mdates.num2date(loc).year != x_min.year:
                lbl.set_visible(False)

    if show_ylabel:
        ax_s.set_ylabel("Time lag (s)", labelpad=4)
    else:
        ax_s.tick_params(axis="y", labelleft=False)

    # panel label
    ax_s.text(0.01, 1.10, "(" + letter + ") " + panel_label,
              transform=ax_s.transAxes,
              fontsize=9, fontweight="bold",
              va="top", ha="left", color="black")

    # KDE axis styling
    ax_k.spines["left"].set_visible(False)
    ax_k.tick_params(axis="y", left=False, labelleft=False)
    ax_k.set_xlim(0, 1.1)
    ax_k.xaxis.set_major_locator(ticker.FixedLocator([0, 0.5, 1.0]))
    ax_k.xaxis.set_major_formatter(ticker.FixedFormatter(["0", "0.5", "1"]))
    ax_k.set_xlabel("KDE", labelpad=3)
    ax_k.grid(axis="y", color="#DDDDDD", lw=0.4, ls="--")
    ax_k.grid(axis="x", visible=False)
    if not show_xlabel:
        ax_k.tick_params(axis="x", labelbottom=False)
    else:
        ax_k.tick_params(axis="x", rotation=0)
        for lbl in ax_k.get_xticklabels():
            lbl.set_ha("center")

# ── shared legend: CM (blue) | PWB (orange) | Mode (black dashed) ─────────────
shared_handles = [
    plt.scatter([], [], marker="o", color=COLOR_F1, s=20, linewidths=0,
                label=FILE1["label"]),
    plt.scatter([], [], marker="o", color=COLOR_F2, s=20, linewidths=0,
                label=FILE2["label"]),
    plt.Line2D([0], [0], color="black", lw=1.2, ls="--",
               label="Mode"),
]

fig.legend(
    handles=shared_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.00),
    ncol=len(shared_handles),
    frameon=False,
    fontsize=9,
    handletextpad=0.4, borderpad=0.5,
    columnspacing=1.5,
)

plt.savefig(SAVE_PATH, dpi=DPI, facecolor="white")
print(f"Saved → {SAVE_PATH}  (font: {_FONT})")
plt.show()