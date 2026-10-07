"""Generates IntegrationSchemesReference.h (the results of the Python code used by IntegrationSchemes for its
comparison) from data_rivais.pkl, written by gen_data.py rivais() of the Python code of the article (v0.6), and a
flat CSV table of the same results.

usage: python3 make_reference.py <data_rivais.pkl> <out.csv> <IntegrationSchemesReference.h>
"""
import pickle
import sys

import numpy as np

pkl, outcsv, outh = sys.argv[1:4]
R = pickle.load(open(pkl, 'rb'))
rows = []  # (test, scheme, control, n, tol, converged, v0, v1, v2, work)


def add(test, scheme, control, n, tol, val, nvals):
    if val is None:
        rows.append((test, scheme, control, n, tol, 0, np.nan, np.nan, np.nan, -1))
        return
    v = [float(x) for x in val[:nvals]] + [np.nan] * (3 - nvals)
    rows.append((test, scheme, control, n, tol, 1, v[0], v[1], v[2], int(val[nvals])))


xb = R['xieB']
rows.append(('xieB', 'reference', '-', 0, 0.0, 1, float(xb['ref'][0]), float(xb['ref'][1]), np.nan, 0))
for (k, n), v in sorted(xb['be'].items()):
    add('xieB', k, '-', n, 0.0, v, 1)
for (k, tol), v in xb['rk'].items():
    add('xieB', k, 'sloan', 1, tol, v, 1)
for (k, tol), v in xb['rk_xie'].items():
    add('xieB', k, 'xie', 1, tol, v, 1)
for (kind, ocr), res in R['kl'].items():
    test = f'kl_{kind}_ocr{ocr}'
    ref = [float(x) for x in res['ref']] + [np.nan] * (3 - len(res['ref']))
    rows.append((test, 'reference', '-', 0, 0.0, 1, ref[0], ref[1], ref[2], 0))
    nv = 2 if kind == 'undrained' else 3
    for (k, n), v in sorted(res['be'].items()):
        add(test, k, '-', n, 0.0, v, nv)
    for (k, par), v in res['rk'].items():
        if kind == 'undrained':
            add(test, k, 'sloan', 10, par, v, nv)
        else:
            add(test, k, 'sloan', par, 1e-4, v, nv)

with open(outcsv, 'w') as f:
    f.write('test,scheme,control,n,tol,converged,v0,v1,v2,work\n')
    for r in rows:
        f.write(','.join([r[0], r[1], r[2], str(r[3]), repr(r[4]), str(r[5])] + [repr(x) for x in r[6:9]] +
                         [str(r[9])]) + '\n')


def cnum(x):
    return 'NAN' if not np.isfinite(x) else repr(float(x))


with open(outh, 'w') as f:
    f.write('''/**
 * @file IntegrationSchemesReference.h
 * @brief Results of the Python code of the article (v0.6, gen_data.py rivais(), file data_rivais.pkl) used by
 * IntegrationSchemes to check the C++ port of rivais.py. Generated from the pickle file; do not edit.
 */
#pragma once

#include "pzreal.h"

#include <cmath>

/**
 * @ingroup mccpaper
 * @brief One result of gen_data.py rivais() (v0.6 Python code), see IntegrationSchemes::TRecord for the meaning
 * of the values: (relative error, -, -) for test B of Xie et al., (p' - p'_exact, q - q_exact, -) for the
 * undrained tests and (q - q_exact, eps_v - eps_v,exact, largest q) for the drained tests; scheme "reference":
 * the closed-form values (p', q) or (q, eps_v, peak q). work = -1 and converged = 0: no solution (None).
 */
struct TIntegrationSchemesPythonRow {
    const char *fTest;    ///< xieB, kl_undrained_ocr1, kl_undrained_ocr10, kl_drained_ocr1 or kl_drained_ocr10
    const char *fScheme;  ///< exact_n, exact_secant, frozen_n, ME2(1), RKDP5(4) or reference
    const char *fControl; ///< step control of the RK schemes (sloan or xie), "-" for the BE schemes
    int fN;               ///< number of increments
    REAL fTol;            ///< tolerance of the RK schemes (0 for the BE schemes)
    int fConverged;       ///< 1 if the test has a solution
    REAL fV[3];           ///< values (see above)
    int fWork;            ///< local Newton iterations (BE) or evaluations of the elastoplastic operator (RK)
};

/** @brief The rows of data_rivais.pkl */
static const TIntegrationSchemesPythonRow gIntegrationSchemesPython[] = {
''')
    for r in rows:
        f.write('    {"%s", "%s", "%s", %d, %s, %d, {%s, %s, %s}, %d},\n'
                % (r[0], r[1], r[2], r[3], repr(r[4]), r[5], cnum(r[6]), cnum(r[7]), cnum(r[8]), r[9]))
    f.write('};\n\n/** @brief Number of rows of gIntegrationSchemesPython */\n'
            'static const int gIntegrationSchemesPythonSize = sizeof(gIntegrationSchemesPython) / '
            'sizeof(gIntegrationSchemesPython[0]);\n')
print(len(rows), 'rows')
