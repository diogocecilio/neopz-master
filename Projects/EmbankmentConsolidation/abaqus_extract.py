"""Compact reference files of the Abaqus/Standard analysis of the embankment (reference/abaqus) from the CSV export
of its ODB (aterro_camclay_tabelas.zip, extracted: historicos.csv, frames.csv, nos.csv, elementos.csv,
metadados.json, campos/*.csv and arquivos_do_job/aterro_camclay.inp).

Usage:
    python3 abaqus_extract.py <directory of the extracted zip> [-o <output directory>]

The Abaqus model (aterro_camclay.inp) is the plane-strain adaptation of the FLAC3D example: 20 x 10 CPE8RP elements
of 1 m (quadratic displacements, bilinear pore pressure, reduced integration), units m, N, Pa, s, y vertical from the
base. Steps: GEOSTATIC (1 s), UNDRAINED_LOADING (50 kPa ramped in 1 s on 0 <= x <= 4 m, drained top),
CONSOLIDATION_EARLY (to 1e5 s) and CONSOLIDATION_LATE (to 1e8 s after the loading). Histories: U2 at the top nodes
x = 0, 2, 4, 6 m (sets MON_U0..MON_U6) and POR at the four corner nodes of the elements centred at (0.5, 9.5) m
(MON_P1) and (1.5, 7.5) m (MON_P2).

Output (kPa, m, s; t = time after the start of the loading, i.e. total_time - 1):
- abaqus_history.csv: t, s_x0, s_x2, s_x4, s_x6 (settlements, positive downwards, relative to the end of the
  geostatic step), pp1, pp2 (total pore pressures, means of the four corner nodes), one row per history frame from
  the start of the loading;
- abaqus_nodal_<tag>.csv (tags undrained, t1e5, t1e6, t1e8): x, y, ux, uy (m, relative to the geostatic step) at every
  node and, at the corner nodes, p (total pore pressure, kPa) and p_excess = p - 10 (10 - y); the mid-side nodes have
  p = p_excess = nan; the tag t1e6 is the frame closest to 1e6 s (9.85e5 s);
- abaqus_nodes.csv, abaqus_elements.csv: the mesh;
- abaqus_frames.csv: the frames of the four steps with t after the loading;
- aterro_camclay.inp and metadados.json: copied.
"""
import argparse
import csv
import json
import os
import shutil
import sys

import numpy as np

GAMMA_W = 10.0  # kN/m^3, as in the FLAC3D data file
SETTLEMENT_NODES = {'s_x0': 621, 's_x2': 625, 's_x4': 629, 's_x6': 633}
PP_NODES = {'pp1': (559, 561, 623, 621), 'pp2': (437, 439, 501, 499)}
TAGS = {'undrained': ('UNDRAINED_LOADING', None), 't1e5': ('CONSOLIDATION_EARLY', None),
        't1e6': ('CONSOLIDATION_LATE', 1.0e6), 't1e8': ('CONSOLIDATION_LATE', None)}


def read_rows(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def histories(src, out):
    rows = read_rows(os.path.join(src, 'historicos.csv'))
    series = {}
    for r in rows:
        if not r['region'].startswith('Node'):
            continue
        node = int(r['region'].split('.')[-1])
        series.setdefault((node, r['variable']), []).append((r['step'], float(r['total_time']), float(r['value'])))
    for k in series:
        series[k].sort(key=lambda e: e[1])
    # settlements relative to the end of the geostatic step; the loading starts at total_time = 1
    times = [t for st, t, v in series[(621, 'U2')] if st != 'GEOSTATIC']
    cols = {'t': np.array(times) - 1.0}
    for name, node in SETTLEMENT_NODES.items():
        s = series[(node, 'U2')]
        u_geo = [v for st, t, v in s if st == 'GEOSTATIC'][-1]
        cols[name] = np.array([-(v - u_geo) for st, t, v in s if st != 'GEOSTATIC'])
    for name, nodes in PP_NODES.items():
        vals = [np.array([v for st, t, v in series[(n, 'POR')] if st != 'GEOSTATIC']) for n in nodes]
        cols[name] = np.mean(vals, axis=0) / 1000.0
    names = ['t', 's_x0', 's_x2', 's_x4', 's_x6', 'pp1', 'pp2']
    with open(os.path.join(out, 'abaqus_history.csv'), 'w') as f:
        f.write(','.join(names) + '\n')
        for i in range(len(times)):
            f.write(','.join('%.10g' % cols[n][i] for n in names) + '\n')
    print('abaqus_history.csv: %d rows; end of loading s = %.4f %.4f %.4f %.4f m, pp1 = %.2f, pp2 = %.2f kPa; '
          'final s = %.4f %.4f %.4f %.4f m, pp1 = %.2f, pp2 = %.2f kPa'
          % ((len(times),) + tuple(cols[n][np.nonzero(cols['t'] <= 1.0)[0][-1]] for n in names[1:])
             + tuple(cols[n][-1] for n in names[1:])))
    return cols


def frames(src, out):
    rows = read_rows(os.path.join(src, 'frames.csv'))
    with open(os.path.join(out, 'abaqus_frames.csv'), 'w') as f:
        f.write('step,frame,t\n')
        for r in rows:
            f.write('%s,%s,%.10g\n' % (r['step'], r['frame'], float(r['total_time']) - 1.0))
    return rows


def field_table(src, step, var):
    index = json.load(open(os.path.join(src, 'metadados.json')))
    for s in index['steps']:
        if s['name'] == step:
            return os.path.join(src, s['variables'][var]['table'])
    sys.exit('step %s not found in metadados.json' % step)


def last_frame(src, step, target):
    """Frame number of the step: the last one, or the one whose time is closest to target (s after the loading)."""
    fr = [r for r in read_rows(os.path.join(src, 'frames.csv')) if r['step'] == step]
    if target is None:
        return int(fr[-1]['frame']), float(fr[-1]['total_time']) - 1.0
    best = min(fr, key=lambda r: abs(float(r['total_time']) - 1.0 - target))
    return int(best['frame']), float(best['total_time']) - 1.0


def nodal_fields(src, out):
    nodes = {int(r['node']): (float(r['x']), float(r['y'])) for r in read_rows(os.path.join(src, 'nos.csv'))}
    # geostatic displacements (frame 1 of GEOSTATIC) are subtracted
    u_geo = {}
    for r in read_rows(field_table(src, 'GEOSTATIC', 'U')):
        if int(r['frame']) == 1:
            u_geo[(int(r['node']), r['component'])] = float(r['value'])
    for tag, (step, target) in TAGS.items():
        frame, t = last_frame(src, step, target)
        u = {}
        for r in read_rows(field_table(src, step, 'U')):
            if int(r['frame']) == frame:
                n, c = int(r['node']), r['component']
                u[(n, c)] = float(r['value']) - u_geo.get((n, c), 0.0)
        p = {}
        for r in read_rows(field_table(src, step, 'POR')):
            if int(r['frame']) == frame:
                p[int(r['node'])] = float(r['value']) / 1000.0
        with open(os.path.join(out, 'abaqus_nodal_%s.csv' % tag), 'w') as f:
            f.write('x,y,ux,uy,p,p_excess\n')
            for n in sorted(nodes):
                x, y = nodes[n]
                if n in p:
                    pe = p[n] - GAMMA_W * (10.0 - y)
                    f.write('%.6g,%.6g,%.10g,%.10g,%.10g,%.10g\n' % (x, y, u[(n, 'U1')], u[(n, 'U2')], p[n], pe))
                else:
                    f.write('%.6g,%.6g,%.10g,%.10g,nan,nan\n' % (x, y, u[(n, 'U1')], u[(n, 'U2')]))
        top = [(nodes[n][0], -u[(n, 'U2')]) for n in nodes if nodes[n][1] == 10.0]
        print('abaqus_nodal_%s.csv: step %s frame %d, t = %.6g s; settlement at x = 0: %.4f m; largest excess pore '
              'pressure %.2f kPa' % (tag, step, frame, t, dict(top)[0.0], max(pe for n in p for pe in [p[n] - GAMMA_W * (10.0 - nodes[n][1])])))


def mesh(src, out):
    shutil.copy(os.path.join(src, 'nos.csv'), os.path.join(out, 'abaqus_nodes.csv'))
    shutil.copy(os.path.join(src, 'elementos.csv'), os.path.join(out, 'abaqus_elements.csv'))
    shutil.copy(os.path.join(src, 'metadados.json'), os.path.join(out, 'metadados.json'))
    shutil.copy(os.path.join(src, 'arquivos_do_job', 'aterro_camclay.inp'), os.path.join(out, 'aterro_camclay.inp'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('src', help='directory of the extracted zip')
    ap.add_argument('-o', '--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference', 'abaqus'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    histories(a.src, a.out)
    frames(a.src, a.out)
    nodal_fields(a.src, a.out)
    mesh(a.src, a.out)
    print('written to', a.out)
