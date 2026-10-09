#!/usr/bin/env python3
"""Extract the Fig. 5 curves of Ceron et al. (2025) from the vector page (SVG of PDF page 9).

The grey band of each panel is a matplotlib fill_between polygon: its lower boundary is the solid curve
-J*(v'_opt)/(k_h H^2 gamma_w^2) and its upper boundary the dashed curve J(u'_FE)/(k_h H^2 gamma_w^2).
Axis calibration uses the tick marks of the page (x: 15..90 deg, y: 0..ymax per panel).

usage: python3 extract_fig5_vector.py <page9.svg> [out.csv]
(the SVG was produced with `pdftocairo -svg -f 9 -l 9 paper.pdf p09.svg`)
"""
import os
import re
import sys

import numpy as np

# (alpha, x at 15 deg, x at 90 deg, y at 0, y at top, value at top) in SVG points, from the tick marks
CAL = [(1, 149.598, 296.648, 129.176, 35.365, 1.0), (2, 326.055, 473.105, 129.176, 35.365, 0.70),
       (4, 149.598, 296.648, 241.746, 147.938, 0.5), (10, 326.055, 473.105, 241.746, 147.938, 0.30)]


def extract(svg_path):
    s = open(svg_path).read()
    body = s[s.find("</defs>"):]
    polys = []
    for p in re.findall(r"<path[^>]*>", body):
        if 'fill="rgb(82.745361%' not in p:      # light-grey fill of the bands
            continue
        d = re.search(r' d="([^"]*)"', p).group(1)
        polys.append(np.array([float(t) for t in re.findall(r"-?\d+\.?\d*", d)]).reshape(-1, 2))
    polys = polys[0::2]                            # each band appears twice (fill + clipped copy)
    rows = []
    betas = np.arange(15, 91, 5.0)
    for P, (a, x15, x90, y0, yt, vt) in zip(polys, CAL):
        imax = int(np.argmax(P[:, 0]))
        low = P[1:imax + 1]
        up = P[imax:][::-1]
        up = up[np.argsort(up[:, 0], kind="stable")]
        xb = x15 + (betas - 15) / 75 * (x90 - x15)
        val = lambda yy: (y0 - yy) / (y0 - yt) * vt
        sol = np.interp(xb, low[:, 0], val(low[:, 1]))
        das = np.interp(xb, up[:, 0], val(up[:, 1]))
        rows += [(a, b, s1, d1) for b, s1, d1 in zip(betas, sol, das)]
    return np.array(rows)


if __name__ == "__main__":
    svg = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                                            "data", "fig5_vector_fill_polygons.csv")
    arr = extract(svg)
    np.savetxt(out, arr, delimiter=",", fmt="%.5f",
               header="Fig.5 of Ceron et al. 2025 extracted from the vector fill_between polygons (PDF page 9)\n"
                      "alpha,beta_deg,minusJstar_vopt_solid,J_uFE_dashed  (normalised by k_h H^2 gamma_w^2)")
    print(f"written {out} ({len(arr)} rows)")
