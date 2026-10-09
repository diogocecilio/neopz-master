#!/usr/bin/env python3
"""Digitize Figs. 5, 8 and 9 of Ceron, Cecilio, Linn & Maghous (IJNAMG 2025,
doi:10.1002/nag.3993) into CSV reference data.

Two independent extraction paths are implemented:

(R) raster   450-dpi renders of the PDF pages (hi-09.png = Fig. 5, hi-13.png =
             Figs. 8 and 9).  The axes boxes, spines and tick marks are located
             in pixels (PIL + numpy), a least-squares pixel->data map is built
             from the tick marks (Fig. 8: log10 y axis from the decade ticks
             10^1, 10^2 and the minor ticks 3..9 x 10^k), and the curves are
             traced from dark pixel runs with coverage-weighted (sub-pixel)
             centres.  Curve separation: Fig. 5 by the band fill on one side of
             the run; Fig. 8 by continuity tracking from seed points picked on
             the figure (runs shared by two lines are discarded); Fig. 9 per
             column (the solid line is the lowest run, the dashes lie above it).
             Fig. 9 is a polyline with vertices every 5 deg: one straight centre
             line is fitted per interval (scanning across the band, with an edge
             fit that ignores the lopsided runs at the dash ends) and the node
             value is the weighted mean of the two adjacent lines.
(V) vector   Figs. 5 and 8 are vector graphics in the PDF (Fig. 9 is an embedded
             300-ppi bitmap).  `pdftocairo -svg` exposes the plotted geometry:
             Fig. 5's grey band is the matplotlib fill_between polygon whose
             vertices are the plotted data; Fig. 8's lines are stroke outlines
             whose long straight edges, shifted inwards by half the stroke
             width, lie exactly on the plotted polyline (vertices recovered by a
             continuous piecewise-linear fit).  Axes maps from the vector ticks.

Outputs (in --out-dir): paper_fig5.csv, paper_fig8.csv (vector values),
paper_fig9.csv (raster), paper_fig5_raster.csv / paper_fig8_raster.csv (raster
values of the vector figures, for validation), paper_fig8_vertices.csv (all
recovered polyline vertices of Fig. 8), digitize_check_fig{5,8,9}.png
(overlays) and digitize_meta.json (uncertainties, omitted/hidden points,
validation statistics).  With --no-vector the raster values go to the main
CSV names.

Usage:  python3 digitize_figures.py [--paper-dir DIR] [--out-dir DIR]
        DIR must contain paper.pdf; hi-09.png / hi-13.png (450-dpi renders of
        PDF pages 9 and 13) are rendered with pdftoppm when missing.
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image, ImageDraw

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "data"))
DEFAULT_PAPER_DIR = os.environ.get(
    "SLOPE_PAPER_DIR",
    "/tmp/claude-0/-home-user-neopz-master/339d491a-b062-561e-a71b-775db071f418/scratchpad/pdf")

DPI = 450.0
PX_PER_PT = DPI / 72.0          # 6.25 px per PDF point

# ----------------------------------------------------------------------------
# Tick label values, read from the page images (they are not machine readable)
# ----------------------------------------------------------------------------
FIG5_XTICKS = [15, 30, 45, 60, 75, 90]
FIG5_PANELS = [  # (alpha, y tick labels bottom->top), panel order: TL, TR, BL, BR
    (1, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
    (2, [0.00, 0.14, 0.28, 0.42, 0.56, 0.70]),
    (4, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]),
    (10, [0.00, 0.06, 0.12, 0.18, 0.24, 0.30]),
]
FIG5_BETAS = [15, 30, 45, 60, 75, 90]

FIG8_XTICKS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
FIG8_SOILS = ["London", "Israeli"]           # left, right panel
FIG8_BETAS = {"London": (30, 60), "Israeli": (35, 60)}
FIG8_HW = [0.0, 0.05] + [round(0.1 * k, 1) for k in range(1, 11)]

FIG9_XTICKS = [15, 30, 45, 60, 75, 90]
FIG9_YTICKS = [0, 1, 2, 3, 4, 5]
FIG9_ALPHAS = [1, 5, 10]                     # left, middle, right panel
FIG9_BETAS = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]


# ============================================================================
# Generic helpers
# ============================================================================
class LinMap:
    """value = a + b * pixel (value is log10(data) for a log axis)."""

    def __init__(self, pix, val, log=False):
        pix = np.asarray(pix, float)
        val = np.asarray(val, float)
        v = np.log10(val) if log else val
        A = np.vstack([np.ones_like(pix), pix]).T
        (self.a, self.b), *_ = np.linalg.lstsq(A, v, rcond=None)
        self.log = log
        self.resid_pix = (v - (self.a + self.b * pix)) / self.b
        self.n = len(pix)

    def to_val(self, p):
        v = self.a + self.b * np.asarray(p, float)
        return 10.0 ** v if self.log else v

    def to_pix(self, val):
        v = np.log10(val) if self.log else np.asarray(val, float)
        return (v - self.a) / self.b

    def lin(self, p):            # linear (log10 for log axes) coordinate
        return self.a + self.b * np.asarray(p, float)

    def rms(self):
        return float(np.sqrt(np.mean(self.resid_pix ** 2)))


def runs_1d(mask):
    """Contiguous True runs of a 1D boolean array -> list of (start, end) inclusive."""
    m = np.concatenate([[False], np.asarray(mask, bool), [False]])
    d = np.diff(m.astype(np.int8))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0] - 1
    return list(zip(starts, ends))


def robust_line(x, y, w=None, n_iter=4, k=3.0):
    """Weighted least-squares line y = a + b x with iterative outlier rejection."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    w = np.ones_like(x) if w is None else np.asarray(w, float)
    keep = np.ones(len(x), bool)
    a = b = np.nan
    for _ in range(n_iter):
        if keep.sum() < 2 or np.ptp(x[keep]) == 0:
            break
        W = w[keep]
        A = np.vstack([np.ones(keep.sum()), x[keep]]).T * np.sqrt(W)[:, None]
        (a, b), *_ = np.linalg.lstsq(A, y[keep] * np.sqrt(W), rcond=None)
        r = y - (a + b * x)
        s = 1.4826 * np.median(np.abs(r[keep])) + 1e-9
        new = np.abs(r) < max(k * s, 0.05)
        if np.array_equal(new, keep):
            break
        keep = new
    res = y[keep] - (a + b * x[keep]) if keep.sum() >= 2 else np.array([np.nan])
    return a, b, keep, float(np.sqrt(np.mean(res ** 2)))


# ============================================================================
# Vector path (pdftocairo -svg)
# ============================================================================
def svg_paths(pdf, page, workdir):
    """Convert one PDF page to SVG and return its drawn (non-glyph) paths.

    Each path: {'fill': (r,g,b) in [0,1] or None, 'subs': [ {'pts': (n,2), 'kind': [...]} ],
    'ncmd': number of drawing commands, 'hasC': bool}.  Coordinates are PDF points
    from the top-left page corner (y down), exactly as the 450-dpi render (x 6.25).
    """
    svg = os.path.join(workdir, "p%02d.svg" % page)
    subprocess.run(["pdftocairo", "-svg", "-f", str(page), "-l", str(page), pdf, svg],
                   check=True)
    s = open(svg).read()
    body = s[s.find("</defs>"):]
    out = []
    for m in re.finditer(r"<path ([^>]*)/>", body):
        a = m.group(1)
        if "stroke=" in a:              # stroked paths on these pages are text rules
            continue
        d = re.search(r'\bd="([^"]*)"', a).group(1)
        f = re.search(r'\bfill="rgb\(([^)]*)\)"', a)
        fill = tuple(float(v.strip().rstrip("%")) / 100.0 for v in f.group(1).split(",")) if f else None
        toks = re.findall(r"[MLCZ]|-?\d+\.?\d*(?:e-?\d+)?", d)
        subs, cur, j, hasC = [], None, 0, False
        while j < len(toks):
            t = toks[j]
            if t == "M":
                cur = {"pts": [(float(toks[j + 1]), float(toks[j + 2]))], "kind": []}
                subs.append(cur)
                j += 3
            elif t == "L":
                cur["pts"].append((float(toks[j + 1]), float(toks[j + 2])))
                cur["kind"].append("L")
                j += 3
            elif t == "C":
                hasC = True
                p0 = np.array(cur["pts"][-1])
                c = [float(v) for v in toks[j + 1:j + 7]]
                p1, p2, p3 = np.array(c[0:2]), np.array(c[2:4]), np.array(c[4:6])
                for tt in np.linspace(0, 1, 9)[1:]:
                    q = (1 - tt) ** 3 * p0 + 3 * (1 - tt) ** 2 * tt * p1 + 3 * (1 - tt) * tt ** 2 * p2 + tt ** 3 * p3
                    cur["pts"].append(tuple(q))
                    cur["kind"].append("C")
                j += 7
            elif t == "Z":
                j += 1
            else:
                raise ValueError("unexpected SVG token %r" % t)
        for sp in subs:
            sp["pts"] = np.array(sp["pts"], float)
        allp = np.vstack([sp["pts"] for sp in subs]) if subs else np.zeros((0, 2))
        out.append({"fill": fill, "subs": subs, "ncmd": len(re.findall("[MLC]", d)),
                    "hasC": hasC, "bbox": (allp[:, 0].min(), allp[:, 0].max(),
                                           allp[:, 1].min(), allp[:, 1].max()) if len(allp) else None})
    return out


def is_dark(fill):
    return fill is not None and max(fill) < 0.2


def vec_axes(paths, min_w=100.0):
    """Axes boxes of a page: white rectangles larger than min_w, refined by the
    dark spine rectangles.  Returns list of dict(L, R, T, B) in points."""
    boxes = []
    for p in paths:
        bb = p["bbox"]
        if p["fill"] and min(p["fill"]) > 0.99 and p["ncmd"] <= 6 and bb[1] - bb[0] > min_w \
                and bb[3] - bb[2] > 50:
            boxes.append({"L": bb[0], "R": bb[1], "T": bb[2], "B": bb[3]})
    boxes.sort(key=lambda b: (round(b["T"] / 50), b["L"]))
    return boxes


def vec_ticks(paths, box, tol=1.5):
    """Tick marks (small dark rectangles) outside the bottom and the left spine.

    Returns (xticks, yticks): xticks = sorted centre x, yticks = list of (centre y, length)."""
    xt, yt = [], []
    for p in paths:
        if not is_dark(p["fill"]) or p["ncmd"] > 6:
            continue
        x0, x1, y0, y1 = p["bbox"]
        w, h = x1 - x0, y1 - y0
        if w < 1.2 and 1.5 < h < 6 and abs(y0 - box["B"]) < tol and box["L"] - 2 < 0.5 * (x0 + x1) < box["R"] + 2:
            xt.append(0.5 * (x0 + x1))
        if h < 1.2 and 1.5 < w < 6 and abs(x1 - box["L"]) < tol and box["T"] - 2 < 0.5 * (y0 + y1) < box["B"] + 2:
            yt.append((0.5 * (y0 + y1), w))
    return sorted(xt), sorted(yt)


def polygon_centroid(P):
    """Area centroid of a closed polygon (vertex list, closing vertex optional)."""
    P = np.asarray(P, float)
    if np.hypot(*(P[-1] - P[0])) < 1e-9:
        P = P[:-1]
    x, y = P[:, 0], P[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    cr = x * yn - xn * y
    A = 0.5 * cr.sum()
    return np.array([np.sum((x + xn) * cr), np.sum((y + yn) * cr)]) / (6.0 * A)


def outline_centerline(sub, w, min_len=1.25):
    """Straight edges (longer than min_len) of a stroke outline polygon, shifted
    inwards by w/2: they lie on the stroked centreline.  Returns list of (p, q)."""
    pts, kind = sub["pts"], list(sub["kind"])
    if np.hypot(*(pts[-1] - pts[0])) > 1e-6:
        pts = np.vstack([pts, pts[:1]])
        kind = kind + ["L"]
    x, y = pts[:, 0], pts[:, 1]
    area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    sgn = 1.0 if area > 0 else -1.0
    segs = []
    for k in range(len(pts) - 1):
        if kind[k] != "L":
            continue
        p, q = pts[k], pts[k + 1]
        d = q - p
        L = np.hypot(*d)
        if L < min_len:
            continue
        t = d / L
        n = np.array([-t[1], t[0]]) * sgn
        segs.append((p + 0.5 * w * n, q + 0.5 * w * n))
    return segs


def outline_width(subs, min_len=1.25):
    """Stroke width = median distance between antiparallel long edges."""
    edges = []
    for sp in subs:
        pts = sp["pts"]
        if np.hypot(*(pts[-1] - pts[0])) > 1e-6:
            pts = np.vstack([pts, pts[:1]])
        for k in range(len(pts) - 1):
            if k < len(sp["kind"]) and sp["kind"][k] != "L":
                continue
            d = pts[k + 1] - pts[k]
            L = np.hypot(*d)
            if L >= min_len:
                edges.append((pts[k], d / L, L))
    ws = []
    for (p, t, L) in edges:
        m = p + 0.5 * L * t
        best = None
        for (p2, t2, L2) in edges:
            if np.dot(t, t2) > -0.99:
                continue
            s = np.dot(m - p2, t2)
            if s < -0.5 or s > L2 + 0.5:
                continue
            perp = abs((m - p2)[0] * t2[1] - (m - p2)[1] * t2[0])
            if best is None or perp < best:
                best = perp
        if best is not None:
            ws.append(best)
    return float(np.median(ws))


def polyline_lsq(pieces, points, nodes, lam=1e-9):
    """Continuous piecewise-linear least-squares fit y(x) with vertices at `nodes`.

    pieces: straight centreline pieces [(p, q)] (exact samples of the curve),
            each clipped to the node interval that contains its midpoint (so the
            few hundredths of a point by which an offset edge overshoots a
            vertex are discarded);
    points: isolated centreline points [(x, y)] (centres of the dots of a
            dash-dot line).
    A tiny second-difference penalty (lam) only matters for nodes without data.
    Returns y(nodes), rms residual of the samples, and the per-node weight
    of the data supporting it (0 -> value interpolated by the penalty)."""
    nodes = np.asarray(nodes, float)
    nn = len(nodes)
    X, Y, W = [], [], []
    for p, q in pieces:
        if p[0] > q[0]:
            p, q = q, p
        if q[0] - p[0] < 1e-9:
            continue
        k = np.searchsorted(nodes, 0.5 * (p[0] + q[0])) - 1
        if k < 0 or k >= nn - 1:
            continue
        a, b = max(p[0], nodes[k]), min(q[0], nodes[k + 1])
        if b - a < 1e-6:
            continue
        sl = (q[1] - p[1]) / (q[0] - p[0])
        for t in np.linspace(0, 1, 6):
            x = a + t * (b - a)
            X.append(x)
            Y.append(p[1] + sl * (x - p[0]))
            W.append((b - a) / 6.0)
    for x, y in points:
        if nodes[0] <= x <= nodes[-1]:
            X.append(x)
            Y.append(y)
            W.append(0.5 * (nodes[1] - nodes[0]) / 6.0)
    X, Y, W = np.array(X), np.array(Y), np.array(W)
    A = np.zeros((len(X), nn))
    k = np.clip(np.searchsorted(nodes, X, side="right") - 1, 0, nn - 2)
    t = (X - nodes[k]) / (nodes[k + 1] - nodes[k])
    A[np.arange(len(X)), k] = 1 - t
    A[np.arange(len(X)), k + 1] = t
    D = np.zeros((nn - 2, nn))
    for i in range(nn - 2):
        D[i, i:i + 3] = (1, -2, 1)
    sw = np.sqrt(W)
    M = np.vstack([A * sw[:, None], np.sqrt(lam * W.sum()) * D])
    rhs = np.concatenate([Y * sw, np.zeros(nn - 2)])
    yv, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    r = Y - A @ yv
    support = (A * W[:, None]).sum(0) / ((nodes[1] - nodes[0]) / 6.0)
    return yv, float(np.sqrt(np.sum(W * r ** 2) / W.sum())), float(np.max(np.abs(r))), support


# ============================================================================
# Raster helpers
# ============================================================================
def load_gray(path):
    return np.asarray(Image.open(path).convert("L"), dtype=float)


def find_lines(G, region, min_len, thr=128.0, horizontal=True):
    """Long straight dark lines (spines) inside region = (r0, r1, c0, c1).

    Returns list of dict(pos (sub-pixel centre), a, b (extent), thick)."""
    r0, r1, c0, c1 = region
    D = G[r0:r1, c0:c1] < thr
    if not horizontal:
        D = D.T
    found = []
    for i in range(D.shape[0]):
        for s, e in runs_1d(D[i]):
            if e - s + 1 >= min_len:
                found.append((i, s, e))
    groups = []
    for i, s, e in found:
        for g in groups:
            if i - g["rows"][-1] <= 1 and min(e, g["b"]) - max(s, g["a"]) > 0.5 * min_len:
                g["rows"].append(i)
                g["a"] = min(g["a"], s)
                g["b"] = max(g["b"], e)
                break
        else:
            groups.append({"rows": [i], "a": s, "b": e})
    out = []
    off_r, off_c = (r0, c0) if horizontal else (c0, r0)
    for g in groups:
        rows = np.arange(g["rows"][0] - 2, g["rows"][-1] + 3)
        a, b = g["a"] + 20, g["b"] - 20
        if b <= a:
            continue
        if horizontal:
            prof = np.array([np.mean(255.0 - G[off_r + r, off_c + a:off_c + b]) for r in rows])
        else:
            prof = np.array([np.mean(255.0 - G[off_c + a:off_c + b, off_r + r]) for r in rows])
        prof = np.clip(prof - prof.min(), 0, None)
        pos = off_r + np.sum(prof * rows) / np.sum(prof)
        out.append({"pos": pos, "a": off_c + g["a"], "b": off_c + g["b"],
                    "thick": float(np.sum(prof) / prof.max())})
    return out


def find_boxes(G, region, min_w, min_h):
    """Axes boxes made of 4 dark spines.  Returns list of dict(L, R, T, B, thick)."""
    H = find_lines(G, region, min_w, horizontal=True)
    V = find_lines(G, region, min_h, horizontal=False)
    boxes = []
    tol = 32   # extents include the outward tick marks (about 23 px long)
    for top in H:
        for bot in H:
            if bot["pos"] - top["pos"] < min_h or abs(bot["a"] - top["a"]) > tol or abs(bot["b"] - top["b"]) > tol:
                continue
            left = [v for v in V if abs(v["pos"] - top["a"]) < tol and abs(v["a"] - top["pos"]) < tol and abs(v["b"] - bot["pos"]) < tol]
            right = [v for v in V if abs(v["pos"] - top["b"]) < tol and abs(v["a"] - top["pos"]) < tol and abs(v["b"] - bot["pos"]) < tol]
            if left and right:
                boxes.append({"L": left[0]["pos"], "R": right[0]["pos"], "T": top["pos"], "B": bot["pos"],
                              "thick": float(np.mean([top["thick"], bot["thick"], left[0]["thick"], right[0]["thick"]]))})
    # remove duplicates (the same box found twice)
    uniq = []
    for b in sorted(boxes, key=lambda b: (round(b["T"] / 100), b["L"])):
        if not any(abs(b["L"] - u["L"]) < 5 and abs(b["T"] - u["T"]) < 5 for u in uniq):
            uniq.append(b)
    return uniq


def raster_ticks(G, box, out_len_px, thr=128.0):
    """Tick marks outside the bottom spine (x) and the left spine (y).

    Returns xt (list of centre columns) and yt (list of (centre row, length px))."""
    h = 0.5 * box["thick"]
    B, L = box["B"], box["L"]
    # x ticks: rows just below the bottom spine
    rows = np.arange(int(np.ceil(B + h + 2)), int(np.floor(B + out_len_px - 1)))
    cols = np.arange(int(box["L"] - 15), int(box["R"] + 16))
    prof = np.mean(255.0 - G[np.ix_(rows, cols)], axis=0)
    xt = []
    for s, e in runs_1d(prof > 0.4 * 255):
        c = cols[max(s - 2, 0):e + 3]
        w = prof[max(s - 2, 0):e + 3]
        xt.append(float(np.sum(c * w) / np.sum(w)))
    # y ticks: columns just left of the left spine
    cols = np.arange(int(np.floor(L - h - 6)), int(np.floor(L - h - 1)))
    rows = np.arange(int(box["T"] - 15), int(box["B"] + 16))
    prof = np.mean(255.0 - G[np.ix_(rows, cols)], axis=1)
    yt = []
    for s, e in runs_1d(prof > 0.4 * 255):
        r = rows[max(s - 2, 0):e + 3]
        w = prof[max(s - 2, 0):e + 3]
        rc = float(np.sum(r * w) / np.sum(w))
        # tick length: dark extent to the left of the spine on the centre row
        row = G[int(round(rc)), :int(L)]
        k = int(L - h - 2)
        while k > 0 and row[k] < thr:
            k -= 1
        yt.append((rc, float(L - k)))
    return xt, yt


def column_runs(gcol, thr, ink):
    """Dark runs of one pixel column with coverage-weighted sub-pixel centres.

    Coverage f = (bg - g) / (bg - ink), with the background taken separately
    above and below the run (white, grid grey or the Fig. 5 fill grey)."""
    n = len(gcol)
    out = []
    for s, e in runs_1d(gcol < thr):
        up = gcol[max(s - 4, 0):max(s - 1, 0)]
        dn = gcol[min(e + 2, n):min(e + 5, n)]
        bg_up = float(up.max()) if len(up) else 255.0
        bg_dn = float(dn.max()) if len(dn) else 255.0
        i0, i1 = max(s - 2, 0), min(e + 3, n)
        idx = np.arange(i0, i1)
        mid = 0.5 * (s + e)
        bg = np.where(idx < mid, bg_up, bg_dn)
        f = np.clip((bg - gcol[i0:i1]) / np.maximum(bg - ink, 1.0), 0.0, 1.0)
        if f.sum() <= 0:
            continue
        # number of band-fill grey pixels (Fig. 5) in 14-px windows above / below
        wu = gcol[max(s - 16, 0):max(s - 2, 0)]
        wd = gcol[min(e + 3, n):min(e + 17, n)]
        nfu = int(np.sum((wu >= 195) & (wu <= 224)))
        nfd = int(np.sum((wd >= 195) & (wd <= 224)))
        out.append({"c": float(np.sum(f * idx) / np.sum(f)), "s": int(s), "e": int(e),
                    "mass": float(f.sum()), "bg_up": bg_up, "bg_dn": bg_dn,
                    "nfill_up": nfu, "nfill_dn": nfd, "rows": idx.astype(float), "f": f})
    return out


def track_curve(runs, col0, y0, direction, gate=5.0, hist=16, short=5, max_gap=60, seed_tol=15.0, weak=None,
                row_lim=None):
    """Follow one curve through per-column dark runs starting near (col0, y0).

    runs: dict col -> list of runs (from column_runs, centres in image rows).
    Two predictions are made from the accepted centres: a linear fit over the
    last `hist` points (stable across dash gaps) and over the last `short`
    points (follows kinks).  A run is accepted when either prediction falls
    inside its extent (+- gate, the gate growing across gaps) and its centre
    is the closest one.  Partial runs at dash ends (much less ink than the
    recent runs) are skipped: their centres are biased.
    Tracking stops when the predicted position leaves row_lim = (rmin, rmax)
    (the curve exits the plotting area).  Returns list of (col, run)."""
    cols = sorted(runs)
    # seed: nearest run to y0 within +-12 columns of col0
    seed = None
    for dc in sorted(range(-12, 13), key=abs):
        for r in runs.get(col0 + dc, []):
            if (r["s"] - seed_tol <= y0 <= r["e"] + seed_tol) and (seed is None or abs(r["c"] - y0) < abs(seed[1]["c"] - y0)):
                seed = (col0 + dc, r)
        if seed is not None:
            break
    if seed is None:
        return []
    out = [seed]
    hx, hy, hm = [seed[0]], [seed[1]["c"]], [seed[1]["mass"]]
    gap = 0
    c = seed[0] + direction
    while cols[0] <= c <= cols[-1]:
        preds = []
        for n in (hist, short):
            xs, ys = hx[-n:], hy[-n:]
            if len(xs) >= 3 and np.ptp(xs) > 0:
                b, a = np.polyfit(xs, ys, 1)
                preds.append(a + b * c)
            else:
                preds.append(ys[-1])
        if row_lim is not None and all(p < row_lim[0] or p > row_lim[1] for p in preds):
            break
        g = gate + 0.3 * gap
        best, bd = None, None
        for r in runs.get(c, []):
            ok = any(r["s"] - g <= p <= r["e"] + g or abs(r["c"] - p) <= g for p in preds)
            if ok:
                d = min(abs(r["c"] - p) for p in preds)
                if best is None or d < bd:
                    best, bd = r, d
        mref = float(np.median(hm[-40:]))
        if best is not None and best["mass"] >= 0.55 * mref:
            out.append((c, best))
            hx.append(c)
            hy.append(best["c"])
            hm.append(best["mass"])
            gap = 0
        else:
            if best is not None and weak is not None:
                weak.append((c, best))          # partial run at a dash end

            gap += 1
            if gap > max_gap:
                break
        c += direction
    return out


def track_both(runs, col0, y0, **kw):
    """Track to both sides of the seed; with weak=list the partial dash-end runs are
    appended to that list."""
    left = track_curve(runs, col0, y0, -1, **kw)
    right = track_curve(runs, col0, y0, +1, **kw)
    d = {c: r for c, r in left}
    d.update({c: r for c, r in right})
    return [(c, d[c]) for c in sorted(d)]


def sample_two_sided(x, y, xt, half_w, excl, xmin, xmax):
    """Value at xt from line fits on [xt-half_w, xt-excl] and [xt+excl, xt+half_w]
    (averaged when both exist).  Exact at polyline vertices, insensitive to
    the bias of runs merged across a vertex."""
    x = np.asarray(x)
    y = np.asarray(y)
    vals = []
    for lo, hi in ((xt - half_w, xt - excl), (xt + excl, xt + half_w)):
        lo, hi = max(lo, xmin), min(hi, xmax)
        m = (x >= lo) & (x <= hi)
        if m.sum() >= 4 and np.ptp(x[m]) > 0.3 * half_w:
            a, b, keep, rms = robust_line(x[m], y[m])
            vals.append(a + b * xt)
    if not vals:
        return np.nan, np.nan
    return float(np.mean(vals)), (float(abs(vals[1] - vals[0])) if len(vals) == 2 else np.nan)


def local_poly(x, y, xt, half_w, deg=2, xmin=-np.inf, xmax=np.inf):
    """Robust local polynomial value at xt using points within xt +- half_w."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lo, hi = max(xt - half_w, xmin), min(xt + half_w, xmax)
    m = (x >= lo) & (x <= hi)
    if m.sum() < deg + 3:
        return np.nan, np.nan
    xs, ys = x[m] - xt, y[m]
    keep = np.ones(len(xs), bool)
    for _ in range(5):
        c = np.polyfit(xs[keep], ys[keep], deg)
        r = ys - np.polyval(c, xs)
        s = 1.4826 * np.median(np.abs(r[keep])) + 1e-9
        new = np.abs(r) < 3.5 * s + 1e-6
        if np.array_equal(new, keep) or new.sum() < deg + 3:
            break
        keep = new
    rms = float(np.sqrt(np.mean(r[keep] ** 2)))
    return float(np.polyval(c, 0.0)), rms


# ============================================================================
# Fig. 5
# ============================================================================
FILL_GREY_PT = 0.827   # fill colour of the band in Fig. 5 (rgb ~ 82.7 %)


def split_band_polygon(P):
    """Split a fill_between polygon into its two x-monotone chains."""
    P = np.asarray(P, float)
    if np.hypot(*(P[-1] - P[0])) < 1e-6:
        P = P[:-1]
    x = P[:, 0]
    xmin, xmax = x.min(), x.max()
    n = len(P)
    imin = [i for i in range(n) if abs(x[i] - xmin) < 1e-3]
    imax = [i for i in range(n) if abs(x[i] - xmax) < 1e-3]
    # rotate so that the polygon starts right after the last xmin vertex
    start = max(imin) if not (0 in imin and n - 1 in imin) else min(imin)
    Q = np.roll(P, -start, axis=0)
    xq = Q[:, 0]
    jmax = [i for i in range(n) if abs(xq[i] - xmax) < 1e-3]
    a = Q[:min(jmax) + 1]
    b = Q[max(jmax):]
    b = np.vstack([b, Q[:1]]) if abs(b[-1, 0] - xmin) > 1e-3 else b
    chains = []
    for c in (a, b):
        c = c[np.argsort(c[:, 0], kind="stable")]
        chains.append(c)
    # upper chain = smaller page y (larger value)
    chains.sort(key=lambda c: np.mean(c[:, 1]))
    return chains


def fig5_vector(pdf, workdir):
    paths = svg_paths(pdf, 9, workdir)
    boxes = vec_axes(paths)
    if len(boxes) != 4:
        raise RuntimeError("Fig. 5: expected 4 vector axes, found %d" % len(boxes))
    res = {}
    for box, (alpha, ylab) in zip(boxes, FIG5_PANELS):
        xt, yt = vec_ticks(paths, box)
        if len(xt) != 6 or len(yt) != 6:
            raise RuntimeError("Fig. 5 vector ticks: %s %s" % (xt, yt))
        xmap = LinMap(xt, FIG5_XTICKS)
        ymap = LinMap([y for y, _ in sorted(yt, reverse=True)], ylab)
        inside = lambda bb: bb[0] > box["L"] - 2 and bb[1] < box["R"] + 2 and bb[2] > box["T"] - 2 and bb[3] < box["B"] + 2
        band = [p for p in paths if p["fill"] and abs(p["fill"][0] - FILL_GREY_PT) < 0.01 and not p["hasC"]
                and p["ncmd"] > 50 and inside(p["bbox"])]
        if len(band) != 1:
            raise RuntimeError("Fig. 5 alpha=%d: %d band polygons" % (alpha, len(band)))
        up, lo = split_band_polygon(band[0]["subs"][0]["pts"])
        curves = {"J_FE": up, "minus_Jstar_opt": lo}
        # independent check: stroke outlines of the dashed (upper) and solid (lower) lines
        lines = [p for p in paths if is_dark(p["fill"]) and inside(p["bbox"]) and p["ncmd"] > 10
                 and p["bbox"][1] - p["bbox"][0] > 100]
        check = {}
        for p in lines:
            w = outline_width(p["subs"])
            pieces = []
            for sp in p["subs"]:
                pieces += outline_centerline(sp, w)
            name = "J_FE" if len(p["subs"]) > 3 else "minus_Jstar_opt"
            ch = curves[name]
            d = []
            for a_, b_ in pieces:
                for t in np.linspace(0.1, 0.9, 5):
                    q = a_ + t * (b_ - a_)
                    if ch[0, 0] < q[0] < ch[-1, 0]:
                        d.append(ymap.to_val(q[1]) - ymap.to_val(np.interp(q[0], ch[:, 0], ch[:, 1])))
            check[name] = (len(p["subs"]), w, float(np.max(np.abs(d))) if d else np.nan,
                           float(np.sqrt(np.mean(np.square(d)))) if d else np.nan)
        vals = {}
        for name, ch in curves.items():
            bx = xmap.to_val(ch[:, 0])
            by = ymap.to_val(ch[:, 1])
            vals[name] = {b: float(np.interp(b, bx, by)) for b in FIG5_BETAS}
            vals[name]["_dense"] = (bx, by)
        res[alpha] = {"vals": vals, "xmap": xmap, "ymap": ymap, "box": box, "check": check}
    return res


def fig5_raster(png):
    G = load_gray(png)
    boxes = find_boxes(G, (100, 1700, 600, 3200), min_w=600, min_h=400)
    if len(boxes) != 4:
        raise RuntimeError("Fig. 5 raster: expected 4 axes, found %d" % len(boxes))
    ink = float(np.percentile(G[int(boxes[0]["T"]):int(boxes[0]["B"]), int(boxes[0]["L"]):int(boxes[0]["R"])], 0.5))
    res = {}
    for box, (alpha, ylab) in zip(boxes, FIG5_PANELS):
        xt, yt = raster_ticks(G, box, out_len_px=3.7 * PX_PER_PT)
        if len(xt) != 6 or len(yt) != 6:
            raise RuntimeError("Fig. 5 raster ticks alpha=%d: %s %s" % (alpha, xt, yt))
        xmap = LinMap(xt, FIG5_XTICKS)
        ymap = LinMap([y for y, _ in sorted(yt, reverse=True)], ylab)
        h = 0.5 * box["thick"]
        r0, r1 = int(np.ceil(box["T"] + h + 2)), int(np.floor(box["B"] - h - 2))
        c0, c1 = int(np.ceil(box["L"] + h + 2)), int(np.floor(box["R"] - h - 2))
        pts = {"J_FE": [], "minus_Jstar_opt": []}
        for c in range(c0, c1 + 1):
            col = G[r0:r1 + 1, c]
            runs = column_runs(col, 120.0, ink)
            # upper curve: band fill below, none above; lower curve: the opposite
            ups = [r for r in runs if r["nfill_dn"] >= 3 and r["nfill_up"] == 0]
            los = [r for r in runs if r["nfill_up"] >= 3 and r["nfill_dn"] == 0]
            if ups:
                r = min(ups, key=lambda r: r["c"])
                pts["J_FE"].append((c, r0 + r["c"], r["mass"]))
            if los:
                r = max(los, key=lambda r: r["c"])
                pts["minus_Jstar_opt"].append((c, r0 + r["c"], r["mass"]))
        vals, dense = {}, {}
        for name, arr in pts.items():
            arr = np.array(arr)
            # reject runs merged with annotation arrows (too much ink in the column)
            m = arr[:, 2]
            med = np.array([np.median(m[max(0, i - 15):i + 16]) for i in range(len(m))])
            ok = m < 1.3 * med
            arr = arr[ok]
            bx = xmap.to_val(arr[:, 0])
            by = ymap.to_val(arr[:, 1])
            dense[name] = (bx, by)
            vals[name] = {}
            for b in FIG5_BETAS:
                v, rms = local_poly(bx, by, b, 3.0, deg=2, xmin=15.0, xmax=90.0)
                vals[name][b] = v
        res[alpha] = {"vals": vals, "dense": dense, "xmap": xmap, "ymap": ymap, "box": box}
    return res


# ============================================================================
# Fig. 8
# ============================================================================
def log_axis_map(ticks, major_vals=(100.0, 10.0)):
    """Log10 axis map from (centre, length) ticks: the long ticks are the decades
    (top first), the short ones the minor ticks k x 10^n assigned by proximity."""
    ticks = sorted(ticks)
    lens = np.array([l for _, l in ticks])
    thr = 0.5 * (lens.min() + lens.max())
    major = [p for p, l in ticks if l > thr]
    if len(major) != len(major_vals):
        raise RuntimeError("log axis: majors %s" % major)
    m0 = LinMap(major, major_vals, log=True)
    pix, vals = list(major), list(major_vals)
    for p, l in ticks:
        if l > thr:
            continue
        v = float(m0.to_val(p))
        e = np.floor(np.log10(v))
        k = round(v / 10 ** e)
        if 2 <= k <= 9 and abs(np.log10(v) - np.log10(k * 10 ** e)) < 0.02:
            pix.append(p)
            vals.append(k * 10 ** e)
    return LinMap(pix, vals, log=True)


def fig8_vector(pdf, workdir):
    paths = svg_paths(pdf, 13, workdir)
    boxes = [b for b in vec_axes(paths) if b["B"] < 200]          # Fig. 8 panels (top of page)
    if len(boxes) != 2:
        raise RuntimeError("Fig. 8: expected 2 vector axes, found %d" % len(boxes))
    nodes = np.round(np.arange(0, 41) * 0.025, 4)
    res = {}
    for box, soil in zip(sorted(boxes, key=lambda b: b["L"]), FIG8_SOILS):
        xt, yt = vec_ticks(paths, box)
        xmap = LinMap(xt, FIG8_XTICKS)
        ymap = log_axis_map(yt)
        ylim = (float(ymap.to_val(box["B"])), float(ymap.to_val(box["T"])))
        lines = [p for p in paths if is_dark(p["fill"]) and p["ncmd"] > 60
                 and p["bbox"][0] > box["L"] - 2 and p["bbox"][1] < box["R"] + 2
                 and p["bbox"][1] - p["bbox"][0] > 100]
        if len(lines) != 6:
            raise RuntimeError("Fig. 8 %s: %d curve paths" % (soil, len(lines)))
        by_style = {"Wu_rp025": [], "vopt": [], "FE": []}
        for p in lines:
            n = len(p["subs"])
            if n == 1:
                style = "Wu_rp025"                      # solid
            else:
                lens = []
                for sp in p["subs"]:
                    P = sp["pts"]
                    lens.append(np.sum(np.hypot(*np.diff(np.vstack([P, P[:1]]), axis=0).T)) / 2)
                lens = np.array(lens)
                # dash-dot: about half of the pieces are short dots; dashed: equal dashes
                style = "FE" if np.mean(lens < 0.4 * np.percentile(lens, 90)) > 0.3 else "vopt"
            by_style[style].append(p)
        curves = {}
        for style, ps in by_style.items():
            if len(ps) != 2:
                raise RuntimeError("Fig. 8 %s: style %s has %d curves" % (soil, style, len(ps)))
            ps.sort(key=lambda p: np.mean(np.vstack([sp["pts"] for sp in p["subs"]])[:, 1]))
            for beta, p in zip(FIG8_BETAS[soil], ps):        # upper curve = flatter slope
                w = outline_width(p["subs"])
                pieces, dots = [], []
                for sp in p["subs"]:
                    P = sp["pts"]
                    half_perimeter = np.sum(np.hypot(*np.diff(np.vstack([P, P[:1]]), axis=0).T)) / 2
                    if half_perimeter - w < 1.5:          # a dot of the dash-dot pattern
                        c = polygon_centroid(P)
                        dots.append((float(xmap.lin(c[0])), float(ymap.lin(c[1]))))
                        continue
                    for a_, b_ in outline_centerline(sp, w):
                        pieces.append((np.array([xmap.lin(a_[0]), ymap.lin(a_[1])]),
                                       np.array([xmap.lin(b_[0]), ymap.lin(b_[1])])))
                # coarsest vertex grid that reproduces the plotted polyline exactly
                # (sample residual at the level of the 1/256 pt coordinate rounding)
                for step in (0.1, 0.05, 0.025):
                    vn = np.round(np.arange(0, round(1 / step) + 1) * step, 4)
                    yv_s, rms, rmax, support = polyline_lsq(pieces, dots, vn)
                    if rms < 3e-5 and rmax < 2.5e-4:
                        break
                yv = np.interp(nodes, vn, yv_s)
                sl = np.diff(yv_s) / np.diff(vn)
                kinks = vn[1:-1][np.abs(np.diff(sl)) > 0.03 + 0.01 * np.abs(sl[1:])]
                curves[(beta, style)] = {"nodes": nodes, "log10H": yv, "rms": rms, "rmax": rmax,
                                         "step": step, "vnodes": vn, "vlog10H": yv_s,
                                         "support": support, "kinks": kinks, "width_pt": w,
                                         "nsub": len(p["subs"]), "pieces": pieces, "dots": dots}
        res[soil] = {"curves": curves, "xmap": xmap, "ymap": ymap, "box": box, "ylim": ylim}
    return res


# Raster Fig. 8: masks for in-axes text (axes fractions u = (col-L)/(R-L), v = (row-T)/(B-T)),
# read from zoomed crops of hi-13.png: legend box with the soil name and the beta labels.
FIG8_MASKS = {
    "London": [(0.52, 0.98, 0.04, 0.19), (0.68, 0.89, 0.34, 0.44), (0.09, 0.28, 0.84, 0.92)],
    "Israeli": [(0.57, 0.98, 0.04, 0.19), (0.68, 0.89, 0.34, 0.44), (0.04, 0.22, 0.89, 0.97)],
}
# Seed points (h_w/H, H_crit) picked on the figure where each curve is isolated, and the
# h_w/H range in which the curve can be followed unambiguously in the bitmap (outside it
# the line runs on top of / crosses another line within one line width).
FIG8_SEEDS = {
    "London": {(30, "Wu_rp025"): [(0.9, 19.5, 0.25, 1.0), (0.02, 92.0, 0.0, 0.045)],
               (30, "vopt"): [(0.9, 14.7, 0.25, 1.0)],
               (30, "FE"): [(0.9, 10.5, 0.07, 1.0)],
               (60, "vopt"): [(0.9, 8.0, 0.0, 1.0)],
               (60, "FE"): [(0.2, 9.9, 0.0, 0.38)],
               (60, "Wu_rp025"): [(0.2, 8.5, 0.0, 0.38)]},
    "Israeli": {(35, "Wu_rp025"): [(0.9, 12.1, 0.085, 1.0), (0.015, 132.0, 0.0, 0.026)],
                (35, "vopt"): [(0.9, 8.8, 0.0, 1.0)],
                (35, "FE"): [(0.9, 7.5, 0.0, 1.0)],
                (60, "vopt"): [(0.9, 4.8, 0.0, 1.0)],
                (60, "FE"): [(0.2, 6.2, 0.0, 0.30)],
                (60, "Wu_rp025"): [(0.2, 5.65, 0.0, 0.30)]},
}


def panel_runs(G, box, thr, ink, masks=(), col_margin=2.0, row_margin=2.0):
    """Per-column dark runs inside an axes box (image row coordinates)."""
    h = 0.5 * box["thick"]
    r0, r1 = int(np.ceil(box["T"] + h + row_margin)), int(np.floor(box["B"] - h - row_margin))
    c0, c1 = int(np.ceil(box["L"] + h + col_margin)), int(np.floor(box["R"] - h - col_margin))
    sub = G[r0:r1 + 1, :].copy()
    W, Hh = box["R"] - box["L"], box["B"] - box["T"]
    for (u0, u1, v0, v1) in masks:
        cc0, cc1 = int(box["L"] + u0 * W), int(box["L"] + u1 * W)
        rr0, rr1 = int(box["T"] + v0 * Hh) - r0, int(box["T"] + v1 * Hh) - r0
        sub[max(rr0, 0):max(rr1, 0), cc0:cc1] = 255.0
    runs = {}
    for c in range(c0, c1 + 1):
        rr = column_runs(sub[:, c], thr, ink)
        for r in rr:
            r["c"] += r0
            r["s"] += r0
            r["e"] += r0
            r["rows"] = r["rows"] + r0
        runs[c] = rr
    return runs, (r0, r1, c0, c1)


def fig8_raster(png):
    G = load_gray(png)
    boxes = [b for b in find_boxes(G, (100, 1300, 600, 3200), min_w=500, min_h=400) if b["B"] < 1300]
    if len(boxes) != 2:
        raise RuntimeError("Fig. 8 raster: expected 2 axes, found %d" % len(boxes))
    res = {}
    for box, soil in zip(sorted(boxes, key=lambda b: b["L"]), FIG8_SOILS):
        xt, yt = raster_ticks(G, box, out_len_px=3.7 * PX_PER_PT)
        xmap = LinMap(xt, FIG8_XTICKS)
        ymap = log_axis_map(yt)
        ink = float(np.percentile(G[int(box["T"]):int(box["B"]), int(box["L"]):int(box["R"])], 0.3))
        runs, (r0, r1, c0, c1) = panel_runs(G, box, 100.0, ink, FIG8_MASKS[soil])
        curves = {}
        ylim = (float(ymap.to_val(box["B"])), float(ymap.to_val(box["T"])))
        tracks, weak = {}, {}
        for key, seeds in FIG8_SEEDS[soil].items():
            tracks[key], weak[key] = [], []
            for hw0, H0, u0, u1 in seeds:
                sub = {c: rr for c, rr in runs.items() if u0 - 1e-9 <= xmap.to_val(c) <= u1 + 1e-9}
                col0 = int(round(xmap.to_pix(hw0)))
                tracks[key] += track_both(sub, col0, float(ymap.to_pix(H0)), gate=5.0, max_gap=60,
                                          row_lim=(r0, r1), weak=weak[key])
        # a run claimed by two tracks belongs to overlapping / crossing lines: its
        # centre is not attributable to either curve and is discarded
        owner = {}
        for key, tr in tracks.items():
            for c, r in tr:
                owner.setdefault((c, r["s"]), set()).add(key)
        for key, tr in tracks.items():
            keep = [(c, r) for c, r in tr if len(owner[(c, r["s"])]) == 1]
            cols = np.array([c for c, _ in keep], float)
            ys = np.array([r["c"] for _, r in keep])
            u = xmap.to_val(cols)
            lv = ymap.lin(ys)
            vals = {}
            for t in FIG8_HW:
                v, spread = sample_two_sided(u, lv, t, 0.03, 0.006, 0.0, 1.0)
                vals[t] = (10 ** v if np.isfinite(v) else np.nan, spread)
            cr = {c: r for c, r in keep}
            for c, r in weak[key]:
                if c not in cr and (c, r["s"]) not in owner:
                    cr[c] = r
            curves[key] = {"hw": u, "log10H": lv, "vals": vals, "n_merged": len(tr) - len(keep), "cruns": cr}
        res[soil] = {"curves": curves, "xmap": xmap, "ymap": ymap, "box": box, "ylim": ylim,
                     "G": G, "ink": ink, "rows": (r0, r1)}
    return res


def fig8_polyline_check(R8, V8, excl=0.008):
    """Accuracy of the Fig. 9 raster method (polyline_raster) measured on Fig. 8:
    the vertex grid of each curve is taken from the vector extraction (only the
    0.05 / 0.1 grids leave long enough segments), the raster node values are
    compared with the vector ones at the nodes inside the tracked range."""
    out = []
    for soil, d in R8.items():
        for key, c in d["curves"].items():
            vc = V8[soil]["curves"][key]
            if vc["step"] < 0.05 or len(c["cruns"]) < 50:
                continue
            nodes = vc["vnodes"]
            ms = np.array([r["mass"] for r in c["cruns"].values()])
            P = polyline_raster(d["G"], d["ink"], c["cruns"], d["xmap"], d["ymap"], nodes, excl,
                                float(np.percentile(ms, 60)), d["rows"][0] + 2, d["rows"][1] - 2, 100.0)
            lo, hi = c["hw"].min(), c["hw"].max()
            for x, v, sg, vv in zip(nodes, P["val"], P["sigma"], vc["vlog10H"]):
                if np.isfinite(v) and lo - 1e-6 <= x <= hi + 1e-6:
                    out.append((soil, key, float(x), float(10 ** vv), float(10 ** v), float(v - vv), float(sg)))
    return out


# ============================================================================
# Fig. 9 (embedded 300-ppi bitmap, raster only)
# ============================================================================
# the "alpha = .." label boxes (axes fractions u, v as for Fig. 8), read from hi-13.png
FIG9_MASKS = {1: [(0.08, 0.45, 0.81, 1.0)], 5: [(0.08, 0.45, 0.81, 1.0)], 10: [(0.08, 0.51, 0.81, 1.0)]}


def row_center(grow, x_pred, half, thr, ink):
    """Coverage-weighted centre and ink mass of the dark run nearest x_pred in one
    image row (window x_pred +- half); None when the run touches the window border."""
    n = len(grow)
    x0, x1 = int(np.floor(x_pred - half)), int(np.ceil(x_pred + half))
    if x0 < 0 or x1 >= n:
        return None
    seg = grow[x0:x1 + 1]
    rr = runs_1d(seg < thr)
    if not rr:
        return None
    s, e = min(rr, key=lambda se: abs(x0 + 0.5 * (se[0] + se[1]) - x_pred))
    if s < 3 or e > len(seg) - 4:
        return None
    lf, rt = seg[max(s - 4, 0):s - 1], seg[e + 2:e + 5]
    bl, br = float(lf.max()), float(rt.max())
    idx = np.arange(s - 2, e + 3)
    bg = np.where(idx < 0.5 * (s + e), bl, br)
    f = np.clip((bg - seg[s - 2:e + 3]) / np.maximum(bg - ink, 1.0), 0, 1)
    return x0 + float(np.sum(f * idx) / np.sum(f)), float(np.sum(f))


def band_center_line(x, c, m, m_ref=None):
    """Centre line of a straight band of constant thickness from scan-line runs.

    Each run gives two edges, c - m/2 and c + m/2 (exact for an anti-aliased box
    profile).  At the ends of the dashes the run is cut by the oblique butt cap:
    one of its edges is then the cap, always displaced towards the inside of the
    band.  The two edge lines (common slope) are therefore fitted with an
    asymmetric rejection: 'low' edge points lying inside the band (positive
    residual) and 'high' edge points inside the band (negative residual) are
    dropped.  Returns (a, b) of the centre line c = a + b x, rms (px), points used."""
    x, c, m = np.asarray(x, float), np.asarray(c, float), np.asarray(m, float)
    lo, hi = c - 0.5 * m, c + 0.5 * m
    # start from runs of the expected thickness (neither cut by a cap nor merged
    # with another line); a merged run has one edge outside the band, which the
    # loose outer bound (4 sigma + 0.5 px) rejects
    m_ref = float(np.median(m)) if m_ref is None else m_ref
    full = np.abs(m - m_ref) < 0.15 * m_ref
    if full.sum() < 3:
        full = m >= 0.85 * np.median(m)
    use_l, use_h = full.copy(), full.copy()
    for _ in range(8):
        nl, nh = use_l.sum(), use_h.sum()
        if nl + nh < 4 or np.ptp(np.concatenate([x[use_l], x[use_h]])) < 2:
            return None
        A = np.zeros((nl + nh, 3))
        A[:nl, 0] = 1
        A[nl:, 1] = 1
        A[:nl, 2] = x[use_l]
        A[nl:, 2] = x[use_h]
        rhs = np.concatenate([lo[use_l], hi[use_h]])
        (al, ah, b), *_ = np.linalg.lstsq(A, rhs, rcond=None)
        rl = lo - (al + b * x)
        rh = hi - (ah + b * x)
        res = np.concatenate([rl[use_l], rh[use_h]])
        sig = max(1.4826 * np.median(np.abs(res)), 0.15)
        new_l = (rl < 2.5 * sig + 0.3) & (rl > -4 * sig - 0.5)
        new_h = (rh > -2.5 * sig - 0.3) & (rh < 4 * sig + 0.5)
        if np.array_equal(new_l, use_l) and np.array_equal(new_h, use_h):
            break
        use_l, use_h = new_l, new_h
    rms = float(np.sqrt(np.mean(res ** 2)))
    xs = np.concatenate([x[use_l], x[use_h]])
    stats = {"n": len(xs), "xbar": float(xs.mean()), "sxx": float(np.sum((xs - xs.mean()) ** 2)), "rms": rms}
    return 0.5 * (al + ah), b, rms, int(use_l.sum() + use_h.sum()), stats


def segment_line_px(G, ink, cruns, c_lo, c_hi, w_px, r_min, r_max, thr):
    """Centre line (pixel coordinates) of one straight polyline segment between
    columns c_lo..c_hi from the runs of a traced curve.

    The scan direction is chosen across the band: columns for |slope| <= 1
    (the column runs of the trace), rows for steeper segments (row runs
    re-measured in the image).  The centre line comes from band_center_line,
    which is insensitive to the lopsided runs at the oblique dash ends.
    Returns (a, b) of row = a + b * col, the rms residual (px), samples used."""
    ent = []
    for c in sorted(cruns):
        if c_lo <= c <= c_hi:
            v = cruns[c]
            ent += [(c, r) for r in (v if isinstance(v, list) else [v])]
    if len(ent) < 3:
        return None
    xc = np.array([c for c, _ in ent], float)
    yc = np.array([r["c"] for _, r in ent])
    mc = np.array([r["mass"] for _, r in ent])
    if np.ptp(xc) < 2:
        return None
    a0 = np.polyfit(xc, yc, 1)
    if abs(a0[0]) <= 1.0:
        L = band_center_line(xc, yc, mc, w_px * np.sqrt(1 + a0[0] ** 2))
        if L is None:
            return None
        a, b, rms, n, st = L
        # standard error (rows) of the line at column cc; x2 for correlated samples
        se = lambda cc: 2.0 * max(st["rms"], 0.1) * np.sqrt(1.0 / st["n"] + (cc - st["xbar"]) ** 2 / max(st["sxx"], 1e-9))
        return a, b, rms, n, se
    # steep: scan rows
    th = np.arctan(abs(a0[0]))
    rows_lo = int(np.ceil(max(min(np.polyval(a0, c_lo), np.polyval(a0, c_hi)), r_min)))
    rows_hi = int(np.floor(min(max(np.polyval(a0, c_lo), np.polyval(a0, c_hi)), r_max)))
    half = 0.5 * w_px / np.sin(th) + 6
    R, X, M = [], [], []
    for r in range(rows_lo, rows_hi + 1):
        xp = (r - a0[1]) / a0[0]
        if xp < c_lo or xp > c_hi:
            continue
        rc = row_center(G[r], xp, half, thr, ink)
        if rc is not None:
            R.append(r)
            X.append(rc[0])
            M.append(rc[1])
    if len(R) < 4:
        return None
    L = band_center_line(R, X, M, w_px * np.sqrt(1 + 1 / a0[0] ** 2))     # col = a1 + b1 * row
    if L is None:
        return None
    a1, b1, rms, n, st = L

    def se(cc):
        rr = (cc - a1) / b1
        return 2.0 * max(st["rms"], 0.1) * np.sqrt(1.0 / st["n"] + (rr - st["xbar"]) ** 2 / max(st["sxx"], 1e-9)) / abs(b1)
    return -a1 / b1, 1.0 / b1, rms, n, se


def polyline_raster(G, ink, cruns, xmap, ymap, nodes, excl, w_px, r_min, r_max, thr):
    """Polyline node values from a traced curve: one centre line per node interval
    (segment_line_px on the columns farther than `excl` from both vertices).
    The node value is the inverse-variance weighted mean of the two adjacent
    interval lines evaluated at the node (data-linear units: log10 for a log
    axis).  Returns dict with val, sigma (statistical, from the line fits, and
    at least half the left/right disagreement), mismatch (right - left), the
    vertex abscissa where the adjacent lines intersect (polyline check) and the
    per-interval lines (a, b, rms_px, n) in data-linear units."""
    lines, ses = [], []
    for k in range(len(nodes) - 1):
        c_lo = int(np.ceil(xmap.to_pix(nodes[k] + excl)))
        c_hi = int(np.floor(xmap.to_pix(nodes[k + 1] - excl)))
        L = segment_line_px(G, ink, cruns, c_lo, c_hi, w_px, r_min, r_max, thr) if c_hi - c_lo >= 4 else None
        if L is None:
            lines.append(None)
            ses.append(None)
            continue
        a, b, rms, n, se = L
        p1, p2 = float(c_lo), float(c_hi)
        x1, x2 = float(xmap.lin(p1)), float(xmap.lin(p2))
        y1, y2 = float(ymap.lin(a + b * p1)), float(ymap.lin(a + b * p2))
        bd = (y2 - y1) / (x2 - x1)
        lines.append((y1 - bd * x1, bd, rms, n))
        ses.append(se)
    nn = len(nodes)
    val, sig, mis, xint = (np.full(nn, np.nan) for _ in range(4))
    for i, xn in enumerate(nodes):
        v, w = [], []
        for k in (i - 1, i):
            if 0 <= k < len(lines) and lines[k] is not None:
                v.append(lines[k][0] + lines[k][1] * xn)
                s_pix = ses[k](float(xmap.to_pix(xn)))
                w.append(1.0 / max(abs(ymap.b) * s_pix, 1e-6) ** 2)
        if not v:
            continue
        v, w = np.array(v), np.array(w)
        val[i] = float(np.sum(w * v) / np.sum(w))
        sig[i] = float(1.0 / np.sqrt(np.sum(w)))
        if len(v) == 2:
            mis[i] = v[1] - v[0]
            sig[i] = max(sig[i], 0.5 * abs(v[1] - v[0]))
            l0, l1 = lines[i - 1], lines[i]
            if abs(l0[1] - l1[1]) > 1e-9:
                xint[i] = (l1[0] - l0[0]) / (l0[1] - l1[1])
    return {"val": val, "sigma": sig, "mismatch": mis, "xint": xint, "lines": lines}


def fig9_raster(png, region=(1700, 2700, 600, 3200), px_per_pt=PX_PER_PT):
    G = load_gray(png)
    boxes = sorted(find_boxes(G, region, min_w=int(70 * px_per_pt), min_h=int(60 * px_per_pt)),
                   key=lambda b: b["L"])
    if len(boxes) != 3:
        raise RuntimeError("Fig. 9 raster: expected 3 axes, found %d" % len(boxes))
    nodes = np.array(FIG9_BETAS, float)
    res = {}
    for box, alpha in zip(boxes, FIG9_ALPHAS):
        xt, yt = raster_ticks(G, box, out_len_px=3.7 * px_per_pt)
        if len(xt) != 6 or len(yt) != 6:
            raise RuntimeError("Fig. 9 ticks alpha=%d: %s %s" % (alpha, xt, yt))
        xmap = LinMap(xt, FIG9_XTICKS)
        ymap = LinMap([y for y, _ in sorted(yt, reverse=True)], FIG9_YTICKS)
        ink = float(np.percentile(G[int(box["T"]):int(box["B"]), int(box["L"]):int(box["R"])], 0.3))
        runs, (r0, r1, c0, c1) = panel_runs(G, box, 110.0, ink, FIG9_MASKS[alpha])
        # Curve separation.  In every column the solid line (-grad u'_FE) is the
        # lowest dark run (it is continuous, and the dashed K^-1 v'_opt curve lies
        # above it over the whole beta range: checked on the page image); all runs
        # above it are dashes of the K^-1 v'_opt curve.  Runs cut by the top or
        # bottom of the plotting area are discarded.  Near beta = 90 deg (alpha = 1)
        # the last dashes overlap the solid line: such merged runs (too much ink for
        # one shallow line) are used through their outer edges only.
        fe_raw, vo_raw = {}, {}
        for c, rr in runs.items():
            if not rr:
                continue
            rr = sorted(rr, key=lambda r: r["c"])
            fe_raw[c] = rr[-1]
            vo_raw[c] = rr[:-1]
        fc = np.array(sorted(fe_raw), float)
        fm = np.array([fe_raw[c]["mass"] for c in sorted(fe_raw)])
        fy = np.array([fe_raw[c]["c"] for c in sorted(fe_raw)])
        vo_list = [(c, r) for c in sorted(vo_raw) for r in vo_raw[c]]
        vc = np.array([c for c, _ in vo_list], float)
        vm = np.array([r["mass"] for _, r in vo_list])
        cut = lambda r: not (r["s"] > r0 + 1 and r["e"] < r1 - 1)
        cr = {"FE": {}, "vopt": {}}
        n_merged = 0
        # line widths (ink thickness normal to the line) from the shallow parts
        sl_fe = np.gradient(fy, fc)
        sh = np.abs(sl_fe) < 0.3
        w_fe = float(np.median(fm[sh] / np.sqrt(1 + sl_fe[sh] ** 2)))
        w_vo = float(np.percentile(vm, 75)) if len(vm) else w_fe
        for i, c in enumerate(sorted(fe_raw)):
            r = fe_raw[c]
            if cut(r):
                continue
            near = np.abs(fc - c) <= 30
            sl = np.polyfit(fc[near], fy[near], 1)[0] if near.sum() > 3 else 0.0
            m_loc = w_fe * np.sqrt(1 + sl ** 2)
            if abs(sl) < 1.0 and r["mass"] > 1.3 * m_loc and not vo_raw[c]:
                vnear = np.abs(vc - c) <= 60
                if vnear.sum() == 0:
                    continue
                m_vo = min(w_vo, float(np.percentile(vm[vnear], 75)))
                rf, rv = dict(r), dict(r)
                rf["c"], rf["mass"] = r["c"] + 0.5 * r["mass"] - 0.5 * m_loc, m_loc
                rv["c"], rv["mass"] = r["c"] - 0.5 * r["mass"] + 0.5 * m_vo, m_vo
                cr["FE"].setdefault(c, []).append(rf)
                cr["vopt"].setdefault(c, []).append(rv)
                n_merged += 1
            else:
                cr["FE"].setdefault(c, []).append(r)
        for c, r in vo_list:
            if not cut(r):
                cr["vopt"].setdefault(c, []).append(r)
        curves = {}
        for key in ("vopt", "FE"):
            cruns = cr[key]
            cols = np.array([c for c in sorted(cruns) for _ in cruns[c]], float)
            ys = np.array([r["c"] for c in sorted(cruns) for r in cruns[c]])
            ms = np.array([r["mass"] for c in sorted(cruns) for r in cruns[c]])
            bx = xmap.to_val(cols)
            gy = ymap.to_val(ys)
            # line width: vertical ink thickness of the shallow parts (|slope| < 0.3)
            order = np.argsort(cols)
            slopes = np.abs(np.gradient(ys[order], cols[order] + 1e-3 * np.arange(len(cols))))
            shallow = slopes < 0.3
            # 75th percentile: the partial runs at the dash ends carry less ink
            w_px = float(np.percentile(ms[order][shallow] / np.sqrt(1 + slopes[shallow] ** 2), 75))
            P = polyline_raster(G, ink, cruns, xmap, ymap, nodes, 1.0, w_px, r0 + 2, r1 - 2, 110.0)
            val, mis, xint, lines = P["val"], P["mismatch"], P["xint"], P["lines"]
            visible = (val >= FIG9_YTICKS[0]) & (val <= FIG9_YTICKS[-1])
            curves[key] = {"beta": bx, "Gamma": gy, "nodes": nodes, "val": val, "mismatch": mis, "sigma": P["sigma"],
                           "xint": xint, "lines": lines, "visible": visible, "w_px": w_px, "cruns": cruns,
                           "n_shared": n_merged, "beta_range": (float(bx.min()), float(bx.max()))}
        res[alpha] = {"curves": curves, "xmap": xmap, "ymap": ymap, "box": box, "masks": FIG9_MASKS[alpha]}
    return res


# ============================================================================
# Output: CSV files, overlays, metadata
# ============================================================================
def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def fmt(v, nd):
    return ("%%.%df" % nd) % v


def draw_polyline(d, pts, color, width=2):
    pts = [(float(x), float(y)) for x, y in pts if np.isfinite(x) and np.isfinite(y)]
    if len(pts) > 1:
        d.line(pts, fill=color, width=width)


def draw_marker(d, x, y, color, r=9, width=3):
    d.ellipse((x - r, y - r, x + r, y + r), outline=color, width=width)


def ensure_render(paper_dir, page, png, dpi=450):
    """Render a PDF page at `dpi` with pdftoppm when the PNG is missing."""
    if os.path.exists(png):
        return png
    pdf = os.path.join(paper_dir, "paper.pdf")
    root = png[:-4] + "_tmp"
    subprocess.run(["pdftoppm", "-r", str(dpi), "-f", str(page), "-l", str(page), "-png", pdf, root], check=True)
    out = [f for f in os.listdir(os.path.dirname(root)) if f.startswith(os.path.basename(root))][0]
    shutil.move(os.path.join(os.path.dirname(root), out), png)
    return png


def run_fig5(pdf, png, out, work, meta):
    V = fig5_vector(pdf, work) if pdf else None
    R = fig5_raster(png)
    rows, rows_r = [], []
    unc = {}
    for alpha in [a for a, _ in FIG5_PANELS]:
        for name in ("J_FE", "minus_Jstar_opt"):
            for b in FIG5_BETAS:
                if V:
                    rows.append([alpha, b, name, fmt(V[alpha]["vals"][name][b], 5)])
                rows_r.append([alpha, b, name, fmt(R[alpha]["vals"][name][b], 5)])
        if V:
            dv = [R[alpha]["vals"][n][b] - V[alpha]["vals"][n][b] for n in ("J_FE", "minus_Jstar_opt") for b in FIG5_BETAS]
            unc[alpha] = {"vector_band_vs_line_outline_max": max(c[2] for c in V[alpha]["check"].values()),
                          "raster_minus_vector_max_abs": float(np.max(np.abs(dv))),
                          "raster_minus_vector_mean": float(np.mean(dv)),
                          "vector_axis_fit_rms_pt": [V[alpha]["xmap"].rms(), V[alpha]["ymap"].rms()],
                          "raster_axis_fit_rms_px": [R[alpha]["xmap"].rms(), R[alpha]["ymap"].rms()]}
    write_csv(os.path.join(out, "paper_fig5.csv"), ["alpha", "beta_deg", "curve", "value"], rows if V else rows_r)
    if V:
        write_csv(os.path.join(out, "paper_fig5_raster.csv"), ["alpha", "beta_deg", "curve", "value"], rows_r)
    # overlay on the 450-dpi render
    im = Image.open(png).convert("RGB")
    d = ImageDraw.Draw(im)
    col = {"J_FE": (220, 0, 0), "minus_Jstar_opt": (0, 90, 255)}
    for alpha, _ in FIG5_PANELS:
        xm, ym = R[alpha]["xmap"], R[alpha]["ymap"]
        for name in col:
            if V:
                bx, by = V[alpha]["vals"][name]["_dense"]
                draw_polyline(d, zip(xm.to_pix(bx), ym.to_pix(by)), col[name], 1)
            rx, ry = R[alpha]["dense"][name]
            for x, y in zip(xm.to_pix(rx[::4]), ym.to_pix(ry[::4])):
                d.point((x, y - 12), fill=(0, 160, 0))
            src = V[alpha]["vals"][name] if V else R[alpha]["vals"][name]
            for b in FIG5_BETAS:
                draw_marker(d, xm.to_pix(b), ym.to_pix(src[b]), (255, 140, 0))
    im.crop((700, 120, 3060, 1700)).save(os.path.join(out, "digitize_check_fig5.png"))
    meta["fig5"] = {
        "source": "vector fill_between polygon of the PDF (page 9)" if V else "raster hi-09.png",
        "curves": {"J_FE": "dashed upper boundary of the grey band, J(u'_FE)",
                   "minus_Jstar_opt": "solid lower boundary, -J*(v'_opt)"},
        "normalization": "k_h H^2 gamma_w^2 (paper text)",
        "uncertainty_per_alpha": unc,
        "uncertainty_note": "vector values: <= 0.0015 (alpha=1), <= 0.0009 (2), <= 0.0003 (4), <= 0.0002 (10) "
                            "= max difference between the band polygon and the independent line outline; "
                            "raster values: within 0.0014 of the vector values",
        "overlay": "red/blue thin lines = vector band edges, green dots (drawn 12 px above the line) = raster "
                   "trace, orange circles = CSV points",
    }
    return V, R


def run_fig8(pdf, png, out, work, meta):
    V = fig8_vector(pdf, work) if pdf else None
    R = fig8_raster(png)
    rows, rows_r, hidden, verts, flags = [], [], [], [], []
    for soil in FIG8_SOILS:
        for beta in FIG8_BETAS[soil]:
            for curve in ("vopt", "FE", "Wu_rp025"):
                key = (beta, curve)
                rc = R[soil]["curves"][key]
                ylo, yhi = R[soil]["ylim"]
                for t in FIG8_HW:
                    v = rc["vals"][t][0]
                    if np.isfinite(v) and ylo <= v <= yhi:
                        rows_r.append([soil, beta, curve, t, fmt(v, 3)])
                if not V:
                    continue
                vc = V[soil]["curves"][key]
                ylo, yhi = V[soil]["ylim"]
                i_node = {round(x, 4): i for i, x in enumerate(vc["nodes"])}
                for t in FIG8_HW:
                    i = i_node[round(t, 4)]
                    H = 10 ** vc["log10H"][i]
                    if ylo <= H <= yhi:
                        rows.append([soil, beta, curve, t, fmt(H, 3)])
                    else:
                        hidden.append({"soil": soil, "beta_deg": beta, "curve": curve, "hw_over_H": t,
                                       "Hcrit_m_from_clipped_vector_geometry": round(float(H), 3),
                                       "axis_range_m": [round(ylo, 2), round(yhi, 2)]})
                for x, lv, sup in zip(vc["vnodes"], vc["vlog10H"], vc["support"]):
                    H = 10 ** lv
                    verts.append([soil, beta, curve, round(float(x), 4), fmt(H, 3), int(ylo <= H <= yhi)])
                    if sup < 1.0:
                        flags.append({"soil": soil, "beta_deg": beta, "curve": curve, "hw_over_H": round(float(x), 4),
                                      "note": "no ink within this vertex' neighbourhood (dash gap); value "
                                              "interpolated / extrapolated by the polyline fit"})
    hdr = ["soil", "beta_deg", "curve", "hw_over_H", "Hcrit_m"]
    if V:
        write_csv(os.path.join(out, "paper_fig8.csv"), hdr, rows)
        write_csv(os.path.join(out, "paper_fig8_vertices.csv"), hdr + ["visible"], verts)
    write_csv(os.path.join(out, "paper_fig8_raster.csv" if V else "paper_fig8.csv"), hdr, rows_r)
    # overlay
    im = Image.open(png).convert("RGB")
    d = ImageDraw.Draw(im)
    col = {"vopt": (0, 150, 0), "FE": (220, 0, 0), "Wu_rp025": (0, 80, 255)}
    for soil in FIG8_SOILS:
        xm, ym = R[soil]["xmap"], R[soil]["ymap"]
        for key, rc in R[soil]["curves"].items():
            for x, y in zip(xm.to_pix(rc["hw"][::3]), ym.to_pix(10 ** rc["log10H"][::3])):
                d.point((x, y - 10), fill=col[key[1]])
            if V:
                vc = V[soil]["curves"][key]
                draw_polyline(d, zip(xm.to_pix(vc["vnodes"]), ym.to_pix(10 ** vc["vlog10H"])), col[key[1]], 1)
        for r in (rows if V else rows_r):
            if r[0] == soil:
                draw_marker(d, xm.to_pix(r[3]), ym.to_pix(float(r[4])), (255, 140, 0), r=7)
    im.crop((850, 150, 3010, 1150)).save(os.path.join(out, "digitize_check_fig8.png"))
    info = {}
    if V:
        dif = []
        for soil in FIG8_SOILS:
            for key, rc in R[soil]["curves"].items():
                vc = V[soil]["curves"][key]
                for t in FIG8_HW:
                    v = rc["vals"][t][0]
                    if np.isfinite(v):
                        dif.append(np.log10(v) - float(np.interp(t, vc["nodes"], vc["log10H"])))
        dif = np.array(dif)
        chk = fig8_polyline_check(R, V)
        dd = np.array([c[5] for c in chk])
        info = {
            "vector_polyline_fit": {"%s beta=%d %s" % (soil, k[0], k[1]): {
                "vertex_step_hw": V[soil]["curves"][k]["step"],
                "rms_log10": V[soil]["curves"][k]["rms"], "max_log10": V[soil]["curves"][k]["rmax"],
                "stroke_width_pt": V[soil]["curves"][k]["width_pt"]}
                for soil in FIG8_SOILS for k in sorted(V[soil]["curves"])},
            "vector_axis_fit_rms_pt": {s_: [V[s_]["xmap"].rms(), V[s_]["ymap"].rms()] for s_ in FIG8_SOILS},
            "raster_two_sided_minus_vector_log10": {"n": int(len(dif)), "rms": float(np.sqrt(np.mean(dif ** 2))),
                                                    "max_abs": float(np.max(np.abs(dif)))},
            "raster_polyline_method_minus_vector_log10": {"n": int(len(dd)), "rms": float(np.sqrt(np.mean(dd ** 2))),
                                                          "max_abs": float(np.max(np.abs(dd)))},
        }
    meta["fig8"] = {
        "source": "vector stroke outlines of the PDF (page 13)" if V else "raster hi-13.png",
        "curves": {"vopt": "dashed, seepage forces K^-1 . v'_opt", "FE": "dash-dot, seepage forces -grad u'_FE",
                   "Wu_rp025": "solid, Wu et al. r_p = 0.25 (digitized only)"},
        "y_axis": "log10, ticks 10^1, 10^2 and minor ticks; axis range %.2f .. %.1f m" % tuple(
            (V or R)["London"]["ylim"]),
        "omitted_clipped_points": hidden,
        "low_support_vertices": flags,
        "uncertainty_note": "vector values: polyline vertices recovered to <= 1.5e-4 in log10(H) (0.03 %), "
                            "axis map 2e-5; points flagged in low_support_vertices (end of a dash-dot line in a "
                            "gap) are interpolated / extrapolated over <= 0.025 h_w/H, about +-0.3 %. Values at "
                            "h_w/H between the data vertices are the plotted straight segments (log-linear "
                            "interpolation). Raster values (paper_fig8_raster.csv): only where the line is not "
                            "overlapped by another one; agreement with the vector values reported below.",
        "validation": info,
        "overlay": "thin coloured lines = vector polylines (green vopt, red FE, blue Wu), dots drawn 10 px above "
                   "= raster trace, orange circles = CSV points",
    }
    return V, R


def run_fig9(png, out, meta, extra_renders=()):
    R = fig9_raster(png)
    # resolution check: same pipeline on other renders of the page
    alt = []
    for path, dpi in extra_renders:
        ppp = dpi / 72.0
        reg = tuple(int(v * ppp / PX_PER_PT) for v in (1700, 2700, 600, 3200))
        alt.append((dpi, fig9_raster(path, region=reg, px_per_pt=ppp)))
    rows, hidden, unc = [], [], {}
    for alpha in FIG9_ALPHAS:
        for curve in ("FE", "vopt"):
            c = R[alpha]["curves"][curve]
            for i, b in enumerate(c["nodes"]):
                v = c["val"][i]
                if not np.isfinite(v):
                    continue
                spread = [a[1][alpha]["curves"][curve]["val"][i] for a in alt]
                spread = [x for x in spread if np.isfinite(x)]
                u = float(np.sqrt(c["sigma"][i] ** 2 + (0.5 * (max(spread + [v]) - min(spread + [v]))) ** 2))
                if c["visible"][i]:
                    rows.append([alpha, int(b), curve, fmt(v, 4)])
                    unc["alpha=%d %s beta=%d" % (alpha, curve, b)] = round(u, 4)
                else:
                    hidden.append({"alpha": alpha, "curve": curve, "beta_deg": int(b), "Gamma_extrapolated": round(float(v), 3),
                                   "note": "above the plotted range (Gamma > 5): straight-line extrapolation of the "
                                           "visible part of the adjacent polyline segment, not a visible point"})
    write_csv(os.path.join(out, "paper_fig9.csv"), ["alpha", "beta_deg", "curve", "Gamma"], rows)
    im = Image.open(png).convert("RGB")
    d = ImageDraw.Draw(im)
    col = {"vopt": (0, 150, 0), "FE": (220, 0, 0)}
    for alpha in FIG9_ALPHAS:
        xm, ym = R[alpha]["xmap"], R[alpha]["ymap"]
        box = R[alpha]["box"]
        for curve, c in R[alpha]["curves"].items():
            for x, y in zip(xm.to_pix(c["beta"][::2]), ym.to_pix(c["Gamma"][::2])):
                d.point((x, y - 14), fill=col[curve])
            ok = np.isfinite(c["val"])
            draw_polyline(d, zip(xm.to_pix(c["nodes"][ok]), ym.to_pix(c["val"][ok])), col[curve], 1)
        for r in rows:
            if r[0] == alpha:
                draw_marker(d, xm.to_pix(r[1]), ym.to_pix(float(r[3])), (255, 140, 0), r=8)
        W, Hh = box["R"] - box["L"], box["B"] - box["T"]
        for (u0, u1, v0, v1) in R[alpha]["masks"]:
            d.rectangle((box["L"] + u0 * W, box["T"] + v0 * Hh, box["L"] + u1 * W, box["T"] + v1 * Hh),
                        outline=(255, 160, 0), width=1)
    im.crop((850, 1820, 3010, 2620)).save(os.path.join(out, "digitize_check_fig9.png"))
    vtx = {"alpha=%d %s" % (a, k): [round(float(x - b), 3) for x, b, m in
                                     zip(R[a]["curves"][k]["xint"], R[a]["curves"][k]["nodes"],
                                         np.abs(np.diff(np.concatenate([[np.nan], [l[1] if l else np.nan for l in R[a]["curves"][k]["lines"]]]))))
                                     if np.isfinite(x) and np.isfinite(m) and m > 0.02]
           for a in FIG9_ALPHAS for k in ("FE", "vopt")}
    uv = np.array(list(unc.values()))
    meta["fig9"] = {
        "source": "raster only (Fig. 9 is an embedded 300-ppi bitmap): hi-13.png, polyline vertices at 5 deg",
        "curves": {"FE": "solid, seepage forces -grad u'_FE", "vopt": "dashed, seepage forces K^-1 . v'_opt"},
        "method": "curves separated per column (solid = lowest run); per 5-deg interval one straight centre line "
                  "(column or row scans across the band, asymmetric edge fit insensitive to dash ends); node = "
                  "weighted mean of the two adjacent lines",
        "vertex_check_deg": "intersection abscissa of adjacent segment lines minus node, at kinks with slope "
                            "change > 0.02/deg: %s" % json.dumps(vtx),
        "uncertainty_Gamma_per_point": unc,
        "uncertainty_note": "u = sqrt(sigma_fit^2 + (half spread over renders at %s dpi)^2); median %.4f, max %.4f. "
                            "Method validated on Fig. 8 against the vector truth (see fig8.validation). Recommended "
                            "+-0.005 for all points except those with a larger u listed above."
                            % ([450] + [dd for dd, _ in alt], float(np.median(uv)), float(np.max(uv))),
        "omitted_hidden_points": hidden,
        "overlay": "thin lines = fitted polylines (red FE, green vopt), dots 14 px above = traced run centres, "
                   "orange circles = CSV points, orange boxes = masked label boxes",
    }
    return R, alt


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--paper-dir", default=DEFAULT_PAPER_DIR)
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    ap.add_argument("--no-vector", action="store_true", help="raster only (no pdftocairo / PDF)")
    ap.add_argument("--no-resolution-check", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    pdf = os.path.join(a.paper_dir, "paper.pdf")
    have_pdf = os.path.exists(pdf) and shutil.which("pdftocairo") is not None and not a.no_vector
    hi09 = ensure_render(a.paper_dir, 9, os.path.join(a.paper_dir, "hi-09.png"))
    hi13 = ensure_render(a.paper_dir, 13, os.path.join(a.paper_dir, "hi-13.png"))
    meta = {"paper": "Ceron, Cecilio, Linn & Maghous, IJNAMG 2025, doi:10.1002/nag.3993",
            "script": os.path.abspath(__file__), "renders_dpi": DPI}
    with tempfile.TemporaryDirectory() as work:
        run_fig5(pdf if have_pdf else None, hi09, a.out_dir, work, meta)
        run_fig8(pdf if have_pdf else None, hi13, a.out_dir, work, meta)
        extra = []
        if not a.no_resolution_check and os.path.exists(pdf) and shutil.which("pdftoppm"):
            for dpi in (300, 600):
                root = os.path.join(work, "p13_%d" % dpi)
                subprocess.run(["pdftoppm", "-r", str(dpi), "-f", "13", "-l", "13", "-png", pdf, root], check=True)
                f = [x for x in os.listdir(work) if x.startswith("p13_%d" % dpi) and x.endswith(".png")][0]
                extra.append((os.path.join(work, f), dpi))
        run_fig9(hi13, a.out_dir, meta, extra)
    with open(os.path.join(a.out_dir, "digitize_meta.json"), "w") as f:
        json.dump(meta, f, indent=1, default=float)
    print("wrote", ", ".join(sorted(x for x in os.listdir(a.out_dir) if x.startswith(("paper_fig", "digitize_")))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
