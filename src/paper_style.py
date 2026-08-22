# -*- coding: utf-8 -*-
"""
paper_style.py — a cohesive, publication-grade, ART-FIRST aesthetic for the "Self-falsifying quantum spectroscopy"
figures. Serif Computer-Modern-style typography (mathtext, no fragile usetex), a warm+cool palette matched to the
prior SQD paper, refined spines/grids, and hand-picked colormaps. Import and call apply().
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --- palette (warm-first, matched to the SQD paper) ---
INK      = '#1B1B1B'   # near-black ink
MUTE     = '#6E6E6E'   # muted grey
WARM     = '#E8770C'   # warm amber
WARMDEEP = '#7A1E10'   # deep oxblood
CRIMSON  = '#C0392B'
SLATE    = '#2471A3'   # cool slate blue
TEAL     = '#12876F'   # deep teal
PLUM     = '#7D3C98'   # plum / violet
GOLD     = '#D4A017'
GROUND   = '#FBFAF7'   # soft warm-white paper ground
PANEL    = '#FFFFFF'

# a warm sequential map (deep oxblood -> amber -> gold cream), for heat-style panels
WARMMAP = LinearSegmentedColormap.from_list('warmpaper', ['#160805', '#5A160B', '#B0360C', WARM, '#F6C453', '#FBEFC7'])
# a diverging warm/cool map for signed quantities
DIVMAP = LinearSegmentedColormap.from_list('divwc', [SLATE, '#9DC3DB', GROUND, '#F0B27A', WARMDEEP])

CYCLE = [SLATE, WARMDEEP, TEAL, WARM, PLUM, CRIMSON, GOLD]


def apply():
    mpl.rcParams.update({
        'figure.facecolor': GROUND, 'axes.facecolor': PANEL, 'savefig.facecolor': GROUND,
        'font.family': 'serif', 'font.serif': ['CMU Serif', 'DejaVu Serif', 'Times New Roman'],
        'mathtext.fontset': 'cm', 'mathtext.rm': 'serif',
        'font.size': 11.5, 'axes.titlesize': 12.5, 'axes.labelsize': 12,
        'axes.edgecolor': INK, 'axes.linewidth': 0.9, 'axes.labelcolor': INK,
        'axes.titlecolor': INK, 'axes.titleweight': 'normal', 'axes.titlepad': 9,
        'xtick.color': INK, 'ytick.color': INK, 'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5,
        'xtick.direction': 'out', 'ytick.direction': 'out', 'xtick.major.width': 0.9, 'ytick.major.width': 0.9,
        'axes.grid': True, 'grid.color': MUTE, 'grid.alpha': 0.14, 'grid.linewidth': 0.6,
        'axes.spines.top': False, 'axes.spines.right': False,
        'legend.frameon': True, 'legend.framealpha': 0.92, 'legend.edgecolor': '#DDD6CB',
        'legend.fontsize': 9.5, 'lines.linewidth': 2.2, 'lines.markersize': 6.5,
        'figure.dpi': 150, 'savefig.dpi': 200, 'savefig.bbox': 'tight',
        'axes.prop_cycle': mpl.cycler(color=CYCLE),
    })


def finish(ax):
    """subtle refinement per-axis: nudge spines inward-feeling, soften ticks."""
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK); ax.spines[s].set_linewidth(0.9)
    ax.tick_params(length=3.5, color=INK)
    return ax
