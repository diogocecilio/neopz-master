#!/usr/bin/env python3
"""Comparison tables of the FEM batch (SlopeSeepageForces fembatch, scripts/run_fem_batch.sh): for every row of
results/cpp/fem_fig9.csv and fem_fig8.csv, the FEM gravity-increase factor of the last refinement cycle, its order-1
extrapolation to h -> 0 (FEMStability.h ExtrapolateCycles; the Gamma_FEM / Hcrit_FEM_m column), the limit-analysis value of the same case (same seepage
field) and the paper's curve (data/paper_fig9.csv; data/paper_fig8.csv interpolated in h_w / H), with the relative
differences. Prints the tables and writes results/cpp/comparison_fem_fig9.csv / comparison_fem_fig8.csv. When
results/fem/outliers/flags.csv exists (scripts/run_fem_outliers.sh, fem_outliers_table.py: diagnosis of the rows where
Gamma_FEM is well below Gamma_LA), its flag and best FEM estimate are shown next to the rows concerned ('flag' column:
the batch value of that row should not be used as it stands; the note says why and what to use instead).

    python3 Projects/SlopeSeepageForces/scripts/fem_batch_table.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
CPP = os.path.join(ROOT, 'results', 'cpp')


def read(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [r for r in csv.DictReader(f)]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def rel(a, b):
    return a / b - 1. if (math.isfinite(a) and math.isfinite(b) and b != 0.) else math.nan


def pct(v):
    return '' if not math.isfinite(v) else f'{100 * v:+.2f} %'


def paper9():
    out = {}
    for r in read(os.path.join(ROOT, 'data', 'paper_fig9.csv')):
        out[(num(r['alpha']), num(r['beta_deg']), r['curve'])] = num(r['Gamma'])
    return out


def paper8():
    """(panel, beta, curve) -> sorted [(hw, Hcrit)] of data/paper_fig8.csv"""
    out = {}
    rows = read(os.path.join(ROOT, 'data', 'paper_fig8.csv'))
    if not rows:
        return out
    keys = list(rows[0].keys())
    for r in rows:
        try:
            k = (r[keys[0]], num(r[keys[1]]), r[keys[2]])
            out.setdefault(k, []).append((num(r[keys[3]]), num(r[keys[4]])))
        except (KeyError, IndexError):
            continue
    for k in out:
        out[k].sort()
    return out


def interp_log(pts, x):
    """log-linear interpolation of a digitized curve (H_crit spans decades), NaN outside its range"""
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1 and y0 > 0 and y1 > 0 and x1 > x0:
            t = (x - x0) / (x1 - x0)
            return math.exp((1 - t) * math.log(y0) + t * math.log(y1))
    for xx, yy in pts:
        if abs(xx - x) < 1e-9:
            return yy
    return math.nan


def flags():
    """results/fem/outliers/flags.csv: (fig, key) -> row (flag, Gamma_FEM_best, Gamma_FEM_best_unc, note);
    key = alpha|beta|curve (fig 9) or panel|beta|curve|hw (fig 8)"""
    out = {}
    for r in read(os.path.join(ROOT, 'results', 'fem', 'outliers', 'flags.csv')):
        out[(r['fig'], r['key'])] = r
    return out


def main():
    fl = flags()
    p9 = paper9()
    rows9 = read(os.path.join(CPP, 'fem_fig9.csv'))
    out9 = []
    if rows9:
        print('Fig. 9 (H = 5 m, h_w = H): FEM last cycle / order-1 extrapolation (h -> 0) / limit analysis / paper')
        print(f'{"alpha":>5} {"beta":>5} {"curve":>5} {"neq":>6} {"FEM last":>9} {"FEM h->0":>9} {"LA":>9} {"paper":>8}'
              f' {"h0/LA-1":>9} {"last/LA-1":>9} {"LA/paper-1":>10} {"time s":>7}  flag (results/fem/outliers/flags.csv)')
    for r in rows9:
        a, b, cv = num(r['alpha']), num(r['beta_deg']), r['curve']
        last, ext, la = num(r['Gamma_FEM_last']), num(r['Gamma_FEM']), num(r['Gamma_LA'])
        pap = p9.get((a, b, cv), math.nan)
        t = num(r['t_fem_s']) + num(r['t_field_s']) + num(r['t_la_s'])
        f = fl.get(('9', f'{a:g}|{b:g}|{cv}'), {})
        print(f'{a:5g} {b:5g} {cv:>5} {r["neq"]:>6} {last:9.4f} {ext:9.4f} {la:9.4f} {pap:8.4f} {pct(rel(ext, la)):>9}'
              f' {pct(rel(last, la)):>9} {pct(rel(la, pap)):>10} {t:7.0f}  {f.get("flag", "")}'
              + (f' (best FEM {f["Gamma_FEM_best"]} +- {f["Gamma_FEM_best_unc"]})' if f.get('Gamma_FEM_best') else ''))
        out9.append(dict(alpha=a, beta_deg=b, curve=cv, neq=r['neq'], Gamma_FEM_last=last, Gamma_FEM_h0=ext,
                         Gamma_LA=la, Gamma_paper=pap, h0_over_LA_minus_1=rel(ext, la),
                         last_over_LA_minus_1=rel(last, la), h0_over_paper_minus_1=rel(ext, pap), time_s=t,
                         flag=f.get('flag', ''), Gamma_FEM_best=f.get('Gamma_FEM_best', ''),
                         Gamma_FEM_best_unc=f.get('Gamma_FEM_best_unc', ''), note=f.get('note', '')))
    if out9:
        with open(os.path.join(CPP, 'comparison_fem_fig9.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(out9[0].keys()))
            w.writeheader()
            w.writerows(out9)
    p8 = paper8()
    rows8 = read(os.path.join(CPP, 'fem_fig8.csv'))
    out8 = []
    if rows8:
        print('\nFig. 8 (alpha = 1): H_crit (m), FEM last cycle / order-1 extrapolation / limit analysis / paper')
        print(f'{"panel":>8} {"beta":>5} {"curve":>5} {"hw/H":>5} {"neq":>6} {"FEM last":>9} {"FEM h->0":>9} {"LA":>9}'
              f' {"paper":>9} {"h0/LA-1":>9} {"last/LA-1":>9} {"time s":>7}  flag (results/fem/outliers/flags.csv)')
    for r in rows8:
        panel, b, cv, hw = r['soil'], num(r['beta_deg']), r['curve'], num(r['hw_over_H'])
        last, ext, la = num(r['Hcrit_FEM_last_m']), num(r['Hcrit_FEM_m']), num(r['Hcrit_LA_m'])
        pap = math.nan
        for k, pts in p8.items():
            if k[0].lower().startswith(panel.lower()) and abs(k[1] - b) < 1e-9 and k[2] == cv:
                pap = interp_log(pts, hw)
        t = num(r['t_fem_s']) + num(r['t_field_s']) + num(r['t_la_s'])
        f = fl.get(('8', f'{panel}|{b:g}|{cv}|{hw:g}'), {})
        print(f'{panel:>8} {b:5g} {cv:>5} {hw:5g} {r["neq"]:>6} {last:9.3f} {ext:9.3f} {la:9.3f} {pap:9.3f}'
              f' {pct(rel(ext, la)):>9} {pct(rel(last, la)):>9} {t:7.0f}  {f.get("flag", "")}'
              + (f' (best FEM {f["Gamma_FEM_best"]} +- {f["Gamma_FEM_best_unc"]} m)' if f.get('Gamma_FEM_best') else ''))
        out8.append(dict(panel=panel, beta_deg=b, curve=cv, hw_over_H=hw, neq=r['neq'], Hcrit_FEM_last=last,
                         Hcrit_FEM_h0=ext, Hcrit_LA=la, Hcrit_paper=pap, h0_over_LA_minus_1=rel(ext, la),
                         last_over_LA_minus_1=rel(last, la), h0_over_paper_minus_1=rel(ext, pap), time_s=t,
                         flag=f.get('flag', ''), Hcrit_FEM_best=f.get('Gamma_FEM_best', ''),
                         Hcrit_FEM_best_unc=f.get('Gamma_FEM_best_unc', ''), note=f.get('note', '')))
    if out8:
        with open(os.path.join(CPP, 'comparison_fem_fig8.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(out8[0].keys()))
            w.writeheader()
            w.writerows(out8)


if __name__ == '__main__':
    main()
