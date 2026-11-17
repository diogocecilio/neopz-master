#!/usr/bin/env python3
"""valida_gradiente.py: o gradiente modal (Cho2010 gradiente=M, uma analise) contra o das diferencas centrais (4M + 1
analises) num ponto xi, no SRM ou no GIM.

O fator e g = FS - 1 (SRM) ou g = Gamma - 1 (GIM). O ponto e xi = 0 (os campos nas medianas, com M termos) ou um arquivo
com os 2M numeros (os M de c' e os M de phi', como o u de um form.json). Os dois gradientes saem com o mesmo fator, a mesma
malha e as mesmas opcoes do form_kl.py (a classe Problema dele).

Uso: python3 -I valida_gradiente.py <exe Cho2010> <pasta> [ponto=<arquivo>|0] [M=10] [metodo=gim] [h=0.5] [paralelos=12]
     [fator=tri2] [nref=4] [rcphi=-0.5]
Saida: <pasta>/valida.json (fator, os dois gradientes, o angulo entre eles, |grad| e o termo dominante) e a tabela na tela."""
import json
import math
import os
import runpy
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
FORM = runpy.run_path(os.path.join(AQUI, "form_kl.py"), run_name="form_kl")


def angulo(a, b):
    na, nb = math.sqrt(sum(x * x for x in a)), math.sqrt(sum(x * x for x in b))
    c = sum(x * y for x, y in zip(a, b)) / (na * nb)
    return math.degrees(math.acos(max(-1., min(1., c))))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    exe, pasta = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    extra = {"ponto": "0"}
    resto = []
    for a in sys.argv[3:]:
        k, v = a.split("=", 1)
        (extra.__setitem__(k, v) if k in extra else resto.append(a))
    p = FORM["parametros"](["metodo=gim", "paralelos=12"] + resto)
    M = p["M"]
    u = [0.] * (2 * M) if extra["ponto"] == "0" else [float(x) for x in open(extra["ponto"]).read().split()]
    if len(u) != 2 * M:
        sys.exit(f"o ponto tem {len(u)} numeros; M = {M} pede {2 * M}")
    prob = FORM["Problema"](exe, pasta, p)
    g_m, grad_m, _ = prob.g_modal(u, "modal")
    p["gradiente"] = "dif"
    g_d, grad_d, fs = prob.g_e_gradiente(u, "dif")
    nome = lambda i: f"{'c' if i < M else 'phi'}{i % M}"
    top = sorted(range(2 * M), key=lambda i: -abs(grad_d[i]))[:6]
    res = {"parametros": p, "ponto": extra["ponto"], "beta_ponto": math.sqrt(sum(x * x for x in u)),
           "fator_modal": g_m + 1., "fator_centro_dif": g_d + 1., "grad_modal": grad_m, "grad_dif": grad_d,
           "angulo_graus": angulo(grad_m, grad_d), "norma_modal": math.sqrt(sum(x * x for x in grad_m)),
           "norma_dif": math.sqrt(sum(x * x for x in grad_d)),
           "dominante": {nome(i): [grad_m[i], grad_d[i]] for i in top}, "avaliacoes": prob.n}
    with open(os.path.join(pasta, "valida.json"), "w") as f:
        json.dump(res, f, indent=1)
    rot = "Gamma" if p["metodo"] == "gim" else "FS"
    print(f"{p['metodo'].upper()}, M = {M}, ponto {extra['ponto']} (|xi| = {res['beta_ponto']:.3f}): {rot} = "
          f"{res['fator_modal']:.4f}; angulo modal x diferencas centrais = {res['angulo_graus']:.1f} graus; "
          f"|grad| {res['norma_modal']:.4f} (modal) e {res['norma_dif']:.4f} (dif.)")
    for k, (a, b) in res["dominante"].items():
        print(f"  {k:6s} modal {a:+.4f}   dif. {b:+.4f}   razao {a / b if b else float('nan'):.3f}")


if __name__ == "__main__":
    main()
