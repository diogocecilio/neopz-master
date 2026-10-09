#!/usr/bin/env python3
"""Summary of the FEM convergence study (scripts/run_fem_convergence.sh): one table row per refinement cycle of each
case log in results/fem/convergence/*.log (output of 'SlopeSeepageForces fs'), with the extrapolations of
FEMStability.h (ExtrapolateCycles) recomputed here for the last cycles, and a CSV of all the cycles
(results/fem/convergence/summary.csv).

    python3 Projects/SlopeSeepageForces/scripts/fem_convergence_table.py [name filter]
"""
import csv
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(HERE, '..', 'results', 'fem', 'convergence')

CYCLE = re.compile(r'^\[FS\] cycle (\d+): (\d+) equations, lambda GI ([0-9.eE+-]+) \(([0-9.eE+-]+) s\); GI plastic zone '
                   r'x in \[([^,]+), ([^\]]+)\], y in \[([^,]+), ([^\]]+)\] \(1 %: x in \[([^,]+), ([^\]]+)\], '
                   r'y in \[([^,]+), ([^\]]+)\]\)')


def parse(path):
    cycles, info = [], {}
    with open(path) as f:
        for line in f:
            m = CYCLE.match(line)
            if m:
                g = m.groups()
                cycles.append(dict(cycle=int(g[0]), neq=int(g[1]), lam=float(g[2]), t=float(g[3]),
                                   zone=[float(v) for v in g[4:8]], zone1=[float(v) for v in g[8:12]]))
            elif line.startswith('stability domain:'):
                info['stab'] = line.split(':', 1)[1].strip()
            elif line.startswith('stability mesh'):
                info['mesh'] = line.strip()
            elif line.startswith('total time'):
                info['total'] = float(line.split()[2])
    return cycles, info


def extrapolations(cyc):
    """lambda_inf by (a) linear in 1/sqrt(neq) through the last two cycles, (b) Richardson in h ~ 2^-k with the
    observed order of the last three cycles (when the differences decrease monotonically)"""
    out = {}
    if len(cyc) >= 2:
        x1, x2 = 1 / math.sqrt(cyc[-2]['neq']), 1 / math.sqrt(cyc[-1]['neq'])
        l1, l2 = cyc[-2]['lam'], cyc[-1]['lam']
        out['sqrtneq2'] = l2 - (l1 - l2) * x2 / (x1 - x2)
    if len(cyc) >= 3:
        d1, d2 = cyc[-3]['lam'] - cyc[-2]['lam'], cyc[-2]['lam'] - cyc[-1]['lam']
        if d1 > 0 and d2 > 0 and d2 < d1:
            p = math.log2(d1 / d2)
            out['order'] = p
            out['richardson'] = cyc[-1]['lam'] - d2 / (2 ** p - 1)
    return out


def main():
    filt = sys.argv[1] if len(sys.argv) > 1 else ''
    rows = []
    for path in sorted(glob.glob(os.path.join(DIR, '*.log'))):
        name = os.path.basename(path)[:-4]
        if filt not in name:
            continue
        cyc, info = parse(path)
        if not cyc:
            continue
        done = 'total' in info
        print(f"\n{name}{'' if done else '  (incomplete)'}: {info.get('stab', '')}; {info.get('mesh', '')}")
        print('  cycle  neq     lambda    d(lambda)   time(s)  zone 1%: x_min  x_max  y_min')
        prev = None
        for c in cyc:
            d = '' if prev is None else f'{c["lam"] / prev - 1:+.4f}'
            z = c['zone1']
            print(f'  {c["cycle"]:5d} {c["neq"]:6d} {c["lam"]:10.5f} {d:>10s} {c["t"]:9.1f}  {z[0]:8.2f} {z[1]:6.2f} {z[2]:6.2f}')
            prev = c['lam']
            rows.append(dict(case=name, cycle=c['cycle'], neq=c['neq'], lam=c['lam'], t=c['t'], x1min=z[0], x1max=z[1],
                             y1min=z[2], complete=int(done)))
        ex = extrapolations(cyc)
        if ex:
            print('  extrapolated: ' + ', '.join(f'{k} {v:.5f}' for k, v in ex.items()))
        if 'total' in info:
            print(f'  total {info["total"]:.0f} s')
    with open(os.path.join(DIR, 'summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['case', 'cycle', 'neq', 'lam', 't', 'x1min', 'x1max', 'y1min', 'complete'])
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
