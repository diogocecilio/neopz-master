"""Common style and helpers of the figure scripts of the Modified Cam-Clay u-p examples.

Style of the figures of the article (scripts/figstyle.py of the Python code): sans-serif font (TeX Gyre Heros
when available), thin lines, light grid and the categorical palette blue, enggreen, magenta and violet.

The figure scripts (Projects/<Example>/plot_figures.py) read the CSV files written by the executables and are
run from the directory with these files:

    python3 <neopz>/Projects/<Example>/plot_figures.py [run directory] [-o output directory]

The figures are written as PDF and PNG in <run directory>/figures (or in the output directory).
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402

COMMON = os.path.dirname(os.path.abspath(__file__))
PROJECTS = os.path.dirname(COMMON)


def _heros(face):
    """File of the TeX Gyre Heros font (TeX installation or usual paths), None if not found."""
    try:
        p = subprocess.run(['kpsewhich', f'texgyreheros-{face}.otf'], capture_output=True, text=True).stdout.strip()
        if p:
            return p
    except (FileNotFoundError, OSError):
        pass
    c = glob.glob(f'/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyreheros-{face}.otf') + \
        glob.glob(f'/usr/share/fonts/**/texgyreheros-{face}.otf', recursive=True) + \
        glob.glob(f'/usr/local/texlive/*/texmf-dist/fonts/opentype/public/tex-gyre/texgyreheros-{face}.otf')
    return c[0] if c else None


def _installed(family):
    """True if matplotlib finds the font family without falling back to its default font."""
    try:
        font_manager.findfont(font_manager.FontProperties(family=family), fallback_to_default=False)
        return True
    except ValueError:
        return False


_found = False
for _face in ('regular', 'bold', 'italic', 'bolditalic'):
    _p = _heros(_face)
    if _p:
        font_manager.fontManager.addfont(_p)
        _found = True
# font of the article (TeX Gyre Heros); otherwise a font with the metrics of Helvetica, so that the legends and
# labels keep their size; DejaVu Sans (default of matplotlib) as the last option
FONT = 'TeX Gyre Heros' if _found else next(
    (f for f in ('Nimbus Sans', 'Nimbus Sans L', 'Liberation Sans', 'FreeSans', 'Arial', 'Helvetica') if _installed(f)),
    'DejaVu Sans')

INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#8f8e88'
ENGGREEN = '#1B7F5B'
C1, C2, C3, C4 = '#2a78d6', ENGGREEN, '#d9649a', '#4a3aa7'  # blue, enggreen, magenta, violet
GREEN_FACE, GREEN_EDGE = '#9fd2bc', (0.06, 0.30, 0.21, 0.55)
C5, C7 = '#e87ba4', '#4a3aa7'
SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
TEXTW = 6.54  # text width of the article (in), 166 mm

plt.rcParams.update({
    'font.family': FONT,
    'mathtext.fontset': 'dejavusans' if FONT == 'DejaVu Sans' else 'custom',
    'font.size': 8,
    'axes.labelsize': 8.5,
    'axes.titlesize': 8.5,
    'axes.titleweight': 'normal',
    'axes.titlelocation': 'left',
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'legend.fontsize': 7.3,
    'legend.frameon': False,
    'axes.edgecolor': AXIS,
    'axes.linewidth': 0.6,
    'axes.labelcolor': INK,
    'xtick.color': INK2,
    'ytick.color': INK2,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'axes.grid': True,
    'grid.color': GRID,
    'grid.linewidth': 0.5,
    'grid.linestyle': '-',
    'axes.axisbelow': True,
    'lines.linewidth': 1.3,
    'lines.markersize': 4.5,
    'figure.dpi': 150,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'text.color': INK,
})
if FONT != 'DejaVu Sans':
    plt.rcParams.update({'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic', 'mathtext.bf': FONT + ':bold'})


def panel_label(ax, s, x=-0.02, y=1.03):
    """Panel label (a), (b), ... above the axes, on the left."""
    ax.text(x, y, s, transform=ax.transAxes, ha='left', va='bottom', fontsize=8.5, color=INK)


def arguments(description):
    """Command line of the figure scripts: run directory (CSV files) and output directory."""
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('rundir', nargs='?', default='.', help='directory with the CSV files of the executable')
    ap.add_argument('-o', '--outdir', default=None, help='output directory (default: <rundir>/figures)')
    args = ap.parse_args()
    args.rundir = os.path.abspath(args.rundir)
    if not os.path.isdir(args.rundir):
        sys.exit(f'run directory not found: {args.rundir}')
    args.outdir = os.path.abspath(args.outdir or os.path.join(args.rundir, 'figures'))
    os.makedirs(args.outdir, exist_ok=True)
    return args


def read_csv(path):
    """Reads a CSV file written by mcc::WriteCSV: dict column name -> numpy array (float)."""
    if not os.path.exists(path):
        sys.exit(f'file not found: {path} (run the executable in this directory first)')
    with open(path, newline='') as f:
        rows = list(csv.reader(f))
    header = [h.strip() for h in rows[0]]
    data = np.array([[float(v) for v in r] for r in rows[1:] if r], dtype=float).reshape(-1, len(header))
    return {h: data[:, i] for i, h in enumerate(header)}


def reference(example_dir, name):
    """Loads a JSON file of digitized reference curves from <example>/reference."""
    with open(os.path.join(example_dir, 'reference', name)) as f:
        return json.load(f)


def save(fig, outdir, name):
    """Writes the figure as PDF and PNG."""
    fig.savefig(os.path.join(outdir, name + '.pdf'), metadata={'CreationDate': None})
    fig.savefig(os.path.join(outdir, name + '.png'), dpi=200)
    plt.close(fig)
    print('written', os.path.join(outdir, name + '.pdf'))
