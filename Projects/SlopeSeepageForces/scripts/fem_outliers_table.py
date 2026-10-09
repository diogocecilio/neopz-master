#!/usr/bin/env python3
"""Summary of the FEM-outlier diagnosis runs (scripts/run_fem_outliers.sh, logs in results/fem/outliers/): for every
'fs' log, lambda of each refinement cycle, the successive differences and their ratios, the extrapolations to h -> 0
(order 1 from the last two cycles; Richardson with the observed order and the geometric limit from the last three),
the 1 % plastic zone of the last cycle against the stability box, and the comparison with the limit analysis of the
same case (the batch row of results/cpp/fem_fig9.csv / fem_fig8.csv for the production box; the 'la' logs of the
same directory for the other boxes). Also checks that the first cycles reproduce the batch row. Prints the table and
writes results/fem/outliers/summary.csv (one row per run) and results/fem/outliers/flags.csv (one row per batch row
concerned: the note that scripts/fem_batch_table.py shows next to it).

    python3 Projects/SlopeSeepageForces/scripts/fem_outliers_table.py
"""
import csv
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
OUT = os.path.join(ROOT, 'results', 'fem', 'outliers')
CPP = os.path.join(ROOT, 'results', 'cpp')

# run name -> batch row (csv, key) and, for a run with another hydraulic box, the 'la' log of the same box
RUNS = {
    'A_fs_a10_b30_n4': dict(csv='fem_fig9.csv', key=dict(alpha='10', beta_deg='30', curve='FE'), box='50/10/30 m'),
    'A_fs_a10_b30_box30': dict(csv='fem_fig9.csv', key=dict(alpha='10', beta_deg='30', curve='FE'), box='50/30/30 m',
                               la='A_la_a10_b30_box30'),
    'A_fs_a5_b30_n4': dict(csv='fem_fig9.csv', key=dict(alpha='5', beta_deg='30', curve='FE'), box='50/10/30 m'),
    'B_fs_isr35_hw0_n5': dict(csv='fem_fig8.csv', key=dict(soil='Israeli', beta_deg='35', curve='FE', hw_over_H='0'),
                              box='50/10/30 H'),
    'B_fs_lon30_hw0_n4': dict(csv='fem_fig8.csv', key=dict(soil='London', beta_deg='30', curve='FE', hw_over_H='0'),
                              box='50/10/30 H'),
    'B_fs_lon30_vopt02_n4': dict(csv='fem_fig8.csv', key=dict(soil='London', beta_deg='30', curve='vopt', hw_over_H='0.2'),
                                 box='50/10/30 H'),
}

CYCLE = re.compile(r'\[FS\] cycle (\d+): (\d+) equations, lambda GI ([-\d.e+]+) \(([\d.e+-]+) s\).*?'
                   r'\(1 %: x in \[([-\d.e+]+), ([-\d.e+]+)\], y in \[([-\d.e+]+), ([-\d.e+]+)\]\)')
BOX = re.compile(r'stability domain: H ([\d.e+-]+) m, beta ([\d.e+-]+) deg, h_w ([\d.e+-]+) m, box x in \[([-\d.e+]+), '
                 r'([-\d.e+]+)\], y in \[([-\d.e+]+), ([-\d.e+]+)\]')


def read_fs_log(path):
    r = dict(cycles=[], done=False)
    with open(path) as f:
        for line in f:
            m = CYCLE.search(line)
            if m:
                r['cycles'].append(dict(k=int(m.group(1)), neq=int(m.group(2)), lam=float(m.group(3)), t=float(m.group(4)),
                                        zone=[float(m.group(i)) for i in range(5, 9)]))
            m = BOX.search(line)
            if m:
                r['H'], r['beta'], r['hw'] = float(m.group(1)), float(m.group(2)), float(m.group(3))
                r['box'] = [float(m.group(i)) for i in range(4, 8)]
            if line.startswith('total time'):
                r['done'] = True
            if 'no collapse up to pmax' in line:
                r['cap'] = r.get('cap', 0) + 1
            if 'failed (100 it.)' in line:
                r['newton_cap'] = r.get('newton_cap', 0) + 1
    return r


def read_la_log(path):
    if not os.path.exists(path):
        return math.nan
    with open(path) as f:
        for line in f:
            m = re.match(r'Gamma = ([-\d.e+]+)', line)
            if m:
                return float(m.group(1))
    return math.nan


def batch_row(csvname, key):
    path = os.path.join(CPP, csvname)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for r in csv.DictReader(f):
            if all(r.get(k) == v for k, v in key.items()):
                return r
    return None


def extrapolations(lams):
    """order 1 (last two), Richardson with the observed order (last three), geometric limit (last three): NaN if n/a"""
    e = dict(h1=math.nan, rich=math.nan, order=math.nan, geom=math.nan, ratio=math.nan)
    n = len(lams)
    if n >= 2:
        e['h1'] = 2. * lams[-1] - lams[-2]
    if n >= 3:
        d1, d2 = lams[-3] - lams[-2], lams[-2] - lams[-1]
        if d1 > 0 and d2 > 0 and d2 < d1:
            e['order'] = math.log2(d1 / d2)
            e['rich'] = lams[-1] - d2 / (2 ** e['order'] - 1)
            e['ratio'] = d2 / d1
            e['geom'] = lams[-1] - d2 * e['ratio'] / (1 - e['ratio'])  # = rich (same formula, written as a series)
    return e


def fmt(v, p=4):
    return '' if (isinstance(v, float) and not math.isfinite(v)) else (f'{v:.{p}f}' if isinstance(v, float) else str(v))


def pct(v):
    return '' if not math.isfinite(v) else f'{100 * v:+.1f} %'


def main():
    rows = []
    flags = {}
    for name, info in RUNS.items():
        path = os.path.join(OUT, name + '.log')
        if not os.path.exists(path):
            continue
        r = read_fs_log(path)
        br = batch_row(info['csv'], info['key'])
        lams = [c['lam'] for c in r['cycles']]
        neqs = [c['neq'] for c in r['cycles']]
        # limit analysis of the same case and box: the batch row (production box) or the la log of the other box
        H = r.get('H', math.nan)
        if 'la' in info:
            la = read_la_log(os.path.join(OUT, info['la'] + '.log'))
        elif br is not None and 'Gamma_LA' in br:
            la = float(br['Gamma_LA'])
        elif br is not None:
            la = float(br['Hcrit_LA_m']) / float(br['H_ref'])  # lambda of the limit analysis at H_ref
        else:
            la = math.nan
        la_prod = math.nan
        if br is not None:
            la_prod = float(br['Gamma_LA']) if 'Gamma_LA' in br else float(br['Hcrit_LA_m']) / float(br['H_ref'])
        # reproduction of the batch row: the common cycles must agree (same mesh sequence, same driver)
        repro = ''
        if br is not None and 'la' not in info:
            bl = [float(x) for x in br['lambda_cycles'].split(';')]
            n = min(len(bl), len(lams))
            dev = max(abs(a - b) / b for a, b in zip(lams[:n], bl[:n])) if n else math.nan
            repro = f'{n} cycles, max dev {dev:.1e}'
        e = extrapolations(lams)
        diffs = [lams[i] - lams[i + 1] for i in range(len(lams) - 1)]
        ratios = [diffs[i + 1] / diffs[i] if diffs[i] else math.nan for i in range(len(diffs) - 1)]
        # 1 % plastic zone of the last cycle against the stability box (units of H; x from O, y from the crest)
        zone = r['cycles'][-1]['zone'] if r['cycles'] else [math.nan] * 4
        box = r.get('box', [math.nan] * 4)
        gap = dict(left=(zone[0] - box[0]) / H, right=(box[1] - zone[1]) / H, bottom=(zone[2] - box[2]) / H)
        rows.append(dict(run=name, box=info['box'], done=r['done'], ncycles=len(lams), neq_last=neqs[-1] if neqs else '',
                         lambda_cycles=';'.join(f'{v:.6g}' for v in lams), neq_cycles=';'.join(str(v) for v in neqs),
                         diff_ratios=';'.join(f'{v:.3f}' for v in ratios), lambda_last=lams[-1] if lams else math.nan,
                         h1=e['h1'], richardson=e['rich'], order=e['order'],
                         lambda_LA_same_box=la, lambda_LA_production=la_prod,
                         last_over_LA=lams[-1] / la - 1 if lams and math.isfinite(la) else math.nan,
                         h1_over_LA=e['h1'] / la - 1 if math.isfinite(la) else math.nan,
                         rich_over_LA=e['rich'] / la - 1 if math.isfinite(la) else math.nan,
                         zone1_gap_left_H=gap['left'], zone1_gap_right_H=gap['right'], zone1_gap_bottom_H=gap['bottom'],
                         continuation_cap_hits=r.get('cap', 0), newton_cap_hits=r.get('newton_cap', 0),
                         reproduces_batch=repro, time_s=sum(c['t'] for c in r['cycles'])))
    if not rows:
        print('no logs in', OUT)
        return
    print('FEM outlier runs (results/fem/outliers): lambda per cycle, extrapolations, limit analysis of the same box')
    for r in rows:
        print(f"\n{r['run']} (box {r['box']}, {'finished' if r['done'] else 'RUNNING'}, {r['time_s']:.0f} s)")
        print(f"  lambda: {r['lambda_cycles']}   neq: {r['neq_cycles']}   ratios of successive differences: {r['diff_ratios']}")
        print(f"  last {fmt(r['lambda_last'])}  order-1 {fmt(r['h1'])}  Richardson {fmt(r['richardson'])} (order {fmt(r['order'], 2)})"
              f"  LA same box {fmt(r['lambda_LA_same_box'])}  LA production box {fmt(r['lambda_LA_production'])}")
        print(f"  last/LA {pct(r['last_over_LA'])}  order-1/LA {pct(r['h1_over_LA'])}  Richardson/LA {pct(r['rich_over_LA'])}"
              f"  1 % zone gap to the box (H): left {fmt(r['zone1_gap_left_H'], 2)} right {fmt(r['zone1_gap_right_H'], 2)}"
              f" bottom {fmt(r['zone1_gap_bottom_H'], 2)}  cap hits {r['continuation_cap_hits']} / Newton {r['newton_cap_hits']}"
              f"  batch: {r['reproduces_batch']}")
    with open(os.path.join(OUT, 'summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print('\nwritten', os.path.join(OUT, 'summary.csv'))


if __name__ == '__main__':
    main()
