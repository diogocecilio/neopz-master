#!/usr/bin/env python3
"""resumo_form.py: uma linha por FORM do form_kl.py (Projects/Cho2010), pronto (form.json) ou em andamento
(form_parcial.json): beta, beta corrigido, P_f = Phi(-beta), iteracoes, avaliacoes do FS e os maiores alfa^2.

Uso: python3 -I resumo_form.py <pasta form_M...> [<pasta> ...]"""
import json
import math
import os
import sys


def phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.))


def maiores(a2, M, k=4):
    top = sorted(range(len(a2)), key=lambda i: -a2[i])[:k]
    return ", ".join(f"{'c' if i < M else 'phi'}{i % M} {a2[i]:.2f}" for i in top)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for pasta in sys.argv[1:]:
        nome = os.path.basename(os.path.normpath(pasta))
        pronto, parcial = os.path.join(pasta, "form.json"), os.path.join(pasta, "form_parcial.json")
        if os.path.exists(pronto):
            d = json.load(open(pronto))
            M = d["parametros"]["M"]
            bc = d.get("beta_corrigido")
            print(f"{nome}: beta = {d['beta']:.4f}" + (f" (corrigido {bc:.4f}, P_f {100 * phi(-bc):.3f}%)" if bc else "") +
                  f", P_f = {100 * d['P_f']:.3f}%; {d.get('iteracoes', len(d['historico']))} iteracoes, "
                  f"{d['avaliacoes']} avaliacoes; alfa^2 por campo c' {d['alfa2_c']:.2f}, phi' {d['alfa2_phi']:.2f}; "
                  f"maiores: {maiores(d['alfa2'], M)}")
        elif os.path.exists(parcial):
            d = json.load(open(parcial))
            M = d["parametros"]["M"]
            h = d["historico"][-1]
            n2 = sum(x * x for x in h["grad"])
            a2 = [x * x / n2 for x in h["grad"]] if n2 > 0 else [0.] * len(h["grad"])
            print(f"{nome}: em andamento, iteracao {h['it']}: g = {h['g']:+.5f}, |u| = {h['beta']:.4f}, "
                  f"{d['avaliacoes']} avaliacoes; maiores alfa^2: {maiores(a2, M)}")
        else:
            print(f"{nome}: sem form.json nem form_parcial.json")


if __name__ == "__main__":
    main()
