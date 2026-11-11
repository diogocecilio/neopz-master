#!/usr/bin/env python3
"""monta_ceron_dig.py: junta as curvas extraidas do PDF de Vargas Ceron et al. (2025) (uma pasta por grupo de figuras,
cada uma com manifest.json e CSV) num ceron_dig.json no formato que ceron_dig.carrega le:

  fig24  -> [x, y] da superficie critica, em coordenadas do modelo (crista em (10, 10), pe em (15, 5));
  fig25a -> {"hist": {"bordas", "alturas"}, "pdf": [x, y]}; fig25b -> [N, P_f]; fig26 -> {"k", "phi", "c"};
  pdfs   -> {familia: {valor: [x, y]}} das Figs. 27a-c, 31a, 32a e 33a.

Uso: python3 -I monta_ceron_dig.py <pasta digit/out> [saida.json, padrao <pasta>/ceron_dig.json]"""

import csv
import json
import math
import os
import re
import sys

FAMILIA = {("27", "a"): "covk", ("27", "b"): "covc", ("27", "c"): "covphi", ("31", "a"): "s", ("32", "a"): "alfa",
           ("33", "a"): "hw"}


def le(arq):
    with open(arq, newline="") as f:
        return list(csv.DictReader(f))


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def fig24_25_26(pasta):
    out = {}
    p = os.path.join(pasta, "fig24_failure_surface.csv")
    if os.path.isfile(p):
        L = le(p)
        out["fig24"] = [[num(r["x_m"]) for r in L], [num(r["y_m"]) for r in L]]
    p = os.path.join(pasta, "fig25a_histogram.csv")
    if os.path.isfile(p):
        L = le(p)
        bordas = [num(L[0]["Gamma_left"])] + [num(r["Gamma_right"]) for r in L]
        alturas = [num(r["density"]) if r.get("bar", "drawn") == "drawn" else 0.0 for r in L]
        out.setdefault("fig25a", {})["hist"] = {"bordas": bordas, "alturas": alturas}
    p = os.path.join(pasta, "fig25a_lognormal_pdf.csv")
    if os.path.isfile(p):
        L = [r for r in le(p) if num(r["Gamma"]) is not None and num(r["Gamma"]) >= 0]
        out.setdefault("fig25a", {})["pdf"] = [[num(r["Gamma"]) for r in L], [max(0.0, num(r["pdf"])) for r in L]]
    p = os.path.join(pasta, "fig25b_pf_convergence.csv")
    if os.path.isfile(p):
        L = le(p)
        out["fig25b"] = [[num(r["N"]) for r in L], [num(r["Pf"]) for r in L]]
    p = os.path.join(pasta, "fig26_spearman.csv")
    if os.path.isfile(p):
        nomes = {"k_v": "k", "kv": "k", "k": "k", "phi": "phi", "c": "c"}
        out["fig26"] = {nomes[r["variable"].strip()]: num(r["rho"]) for r in le(p) if r["variable"].strip() in nomes}
    return out


def _valor(entrada):
    for k in ("value", "valor", "parameter_value"):
        if num(entrada.get(k)) is not None:
            return num(entrada[k])
    txt = " ".join(str(entrada.get(k, "")) for k in ("legend_label", "legend", "label"))
    m = re.search(r"(-?\d+(?:[.,]\d+)?)\s*(%?)", txt)
    if not m:
        return None
    v = float(m.group(1).replace(",", "."))
    return v / 100. if m.group(2) == "%" else v


def curvas_pdf(pasta):
    """as densidades dos paineis (a) dos manifests das Figs. 27, 31, 32 e 33"""
    arq = os.path.join(pasta, "manifest.json")
    if not os.path.isfile(arq):
        return {}
    m = json.load(open(arq))
    entradas = (m.get("curves") or m.get("files") or []) if isinstance(m, dict) else m
    out = {}
    for e in entradas:
        fig, painel = str(e.get("figure", "")).strip(), str(e.get("panel", "")).strip().strip("()")
        fam = FAMILIA.get((fig, painel))
        if not fam:
            continue
        v = _valor(e)
        f = e.get("file") or e.get("csv")
        if v is None or not f:
            continue
        f = f if os.path.isabs(f) else os.path.join(pasta, f)
        L = le(f)
        if not L:
            continue
        cols = [c for c in L[0].keys() if num(L[0][c]) is not None]
        cx = next((c for c in cols if c.lower() in ("gamma", "x")), cols[0])
        cy = next((c for c in cols if c.lower() in ("pdf", "density", "y")), cols[1])
        pts = [(num(r[cx]), num(r[cy])) for r in L if num(r[cx]) is not None and num(r[cy]) is not None]
        pts.sort()
        out.setdefault(fam, {})[repr(v)] = [[p[0] for p in pts], [max(0.0, p[1]) for p in pts]]
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    base = sys.argv[1]
    saida = sys.argv[2] if len(sys.argv) > 2 else os.path.join(base, "ceron_dig.json")
    d = {"pdfs": {}}
    if os.path.isdir(os.path.join(base, "fig24_25_26")):
        d.update(fig24_25_26(os.path.join(base, "fig24_25_26")))
    for g in ("fig27", "fig31", "fig32", "fig33"):
        if os.path.isdir(os.path.join(base, g)):
            for fam, cur in curvas_pdf(os.path.join(base, g)).items():
                d["pdfs"].setdefault(fam, {}).update(cur)
    with open(saida, "w") as f:
        json.dump(d, f)
    print(saida)
    print("fig24: %s pontos; fig25a: %s; fig25b: %s pontos; fig26: %s" % (
        len(d["fig24"][0]) if d.get("fig24") else 0, sorted(d.get("fig25a", {})),
        len(d["fig25b"][0]) if d.get("fig25b") else 0, d.get("fig26")))
    for fam, cur in sorted(d["pdfs"].items()):
        print("densidades %s: %s" % (fam, sorted(float(v) for v in cur)))


if __name__ == "__main__":
    main()
