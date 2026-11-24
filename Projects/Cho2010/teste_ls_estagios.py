#!/usr/bin/env python3
"""teste_ls_estagios.py: teste T0 do LS em estagios (opcao C, nota viabilidade_ls_adaptativo.md, secao 10) com o
executavel falso coluna_falsa.py (estado limite linear em 100 dimensoes, P_f exata = Phi(-beta)).

    python3 teste_ls_estagios.py <pasta> [reps=10] [n1=50] [n2=100] [ang=60] [beta=2.2] [modos=5] [paralelos=4]
                                 [baixos=0.1] [semente=1]

A direcao verdadeira a fica nos modos 0-4 (como o sinal da coluna do Cho). A direcao inicial ruim fica a ang graus
dela; o desvio tem a fracao 'baixos' da norma^2 nos modos 0-4 e o resto nos modos 5-99 (o ruido de um piloto pequeno).
Cada repeticao r: estagio 1 = line_sampling.py na direcao inicial, n1 linhas, semente 1000+2r; direcao_ls.py nas
linhas do estagio 1 (truncada nos modos 0..modos-1); estagio 2 = n2 linhas, semente 1001+2r (Z novas); combina_ls.py
com os pesos n1:n2 (os planejados) e so o estagio 2 (w=0,1); direcao final = direcao_ls.py nas linhas dos 2 estagios e
so nas do estagio 2 (a de um 3o estagio).
Criterios (nota, T0 (b)): a media das P_f combinadas a < 2 EP (da media) da exata; |P - P_exata| < 2 EP na maioria
das repeticoes (~95%); a razao entre a variancia estimada (media dos EP^2) e a empirica entre as repeticoes perto de 1
(so indicativo com poucas repeticoes: o EP relativo da variancia empirica e ~raiz(2/(reps-1))); angulo final < 10 graus.
Tambem o CoV por linha previsto pelo modelo linear (quadratura) para o angulo de cada estagio, contra o medido.
Saidas: <pasta>/r<r>_* (as rodadas) e <pasta>/resumo.json.
"""
import json
import math
import os
import random
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
LS, DL, CB, FALSO = (os.path.join(AQUI, x) for x in ("line_sampling.py", "direcao_ls.py", "combina_ls.py",
                                                      "coluna_falsa.py"))


def cauda(c):
    return 0.5 * math.erfc(c / math.sqrt(2.))


def unit(v):
    n = math.sqrt(sum(x * x for x in v))
    return [x / n for x in v]


def angulo(a, b):
    a, b = unit(a), unit(b)
    return math.degrees(math.acos(max(-1., min(1., sum(x * y for x, y in zip(a, b))))))


def cov_linear(ang, beta):
    """CoV por linha do LS num estado limite linear (indice beta) com a direcao a ang graus da verdadeira: c_i =
    (beta - sen Z)/cos, Z ~ N(0, 1)"""
    s, c = math.sin(math.radians(ang)), math.cos(math.radians(ang))
    if s < 1e-12:
        return 0.
    h, m1, m2 = 0.002, 0., 0.
    for k in range(-5000, 5001):
        z = k * h
        p, w = cauda((beta - s * z) / c), math.exp(-0.5 * z * z) / math.sqrt(2. * math.pi) * h
        m1 += p * w
        m2 += p * p * w
    return math.sqrt(max(m2 - m1 * m1, 0.)) / m1


def roda(cmd, log):
    with open(log, "w") as f:
        subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    pasta = os.path.abspath(sys.argv[1])
    op = {"reps": 10, "n1": 50, "n2": 100, "ang": 60., "beta": 2.2, "modos": 5, "paralelos": 4, "baixos": 0.1,
          "semente": 1}
    for a in sys.argv[2:]:
        k, v = a.split("=", 1)
        op[k] = type(op[k])(v)
    if not os.access(FALSO, os.X_OK):
        sys.exit(f"{FALSO} sem permissao de execucao (chmod +x)")
    os.makedirs(pasta, exist_ok=True)
    rng, K, beta = random.Random(op["semente"]), 100, op["beta"]
    a = unit([0.25, -0.85, 0.10, -0.15, 0.30] + [0.] * (K - 5))
    wl = [rng.gauss(0., 1.) for _ in range(5)]
    d = sum(x * y for x, y in zip(wl, a))
    wl = unit([x - d * y for x, y in zip(wl, a[:5])])
    wh = unit([rng.gauss(0., 1.) for _ in range(K - 5)])
    w = [math.sqrt(op["baixos"]) * x for x in wl] + [math.sqrt(1. - op["baixos"]) * x for x in wh]
    t = math.radians(op["ang"])
    e0 = [math.cos(t) * x + math.sin(t) * y for x, y in zip(a, w)]
    alvo = os.path.join(pasta, "alvo.txt")
    open(alvo, "w").write("\n".join(repr(x) for x in a) + "\n")
    form1 = os.path.join(pasta, "form_inicial.json")
    json.dump({"parametros": {"prog": "coluna", "falha": "fs", "T": 18., "M": K, "threads": 0,
                              "exe": [f"alvo={alvo}", f"beta={beta!r}", "escala=0.1"]},
               "beta": beta, "u": [beta * x for x in e0],
               "origem": f"teste_ls_estagios.py: direcao inicial a {op['ang']:g} graus da verdadeira"},
              open(form1, "w"), indent=1)
    exato = cauda(beta)
    e0t = e0[:op["modos"]] + [0.] * (K - op["modos"])
    print(f"P_f exata = {exato:.6e}; direcao inicial: {angulo(e0, a):.1f} graus (truncada nos modos "
          f"0-{op['modos'] - 1}: {angulo(e0t, a):.1f}); CoV por linha linear: {cov_linear(angulo(e0, a), beta):.2f}",
          flush=True)
    py, comum = sys.executable, [f"paralelos={op['paralelos']}", f"K={K}", "cov=0"]
    res = []
    for r in range(op["reps"]):
        p1, p2 = os.path.join(pasta, f"r{r}_e1"), os.path.join(pasta, f"r{r}_e2")
        f2, f3 = os.path.join(pasta, f"r{r}_form2.json"), os.path.join(pasta, f"r{r}_form3.json")
        roda([py, "-I", LS, FALSO, form1, p1, f"n={op['n1']}", f"semente={1000 + 2 * r}"] + comum, p1 + ".log")
        roda([py, "-I", DL, p1, f"saida={f2}", f"modos={op['modos']}"], f2 + ".log")
        roda([py, "-I", LS, FALSO, f2, p2, f"n={op['n2']}", f"semente={1001 + 2 * r}"] + comum, p2 + ".log")
        cj, cj2 = os.path.join(pasta, f"r{r}_comb.json"), os.path.join(pasta, f"r{r}_comb_so2.json")
        roda([py, "-I", CB, p1, p2, f"json={cj}", "boot=0"], cj + ".log")
        roda([py, "-I", CB, p1, p2, "w=0,1", f"json={cj2}", "boot=0"], cj2 + ".log")
        roda([py, "-I", DL, p1, p2, f"saida={f3}", f"modos={op['modos']}"], f3 + ".log")
        f3b = os.path.join(pasta, f"r{r}_form3_so2.json")
        roda([py, "-I", DL, p2, f"saida={f3b}", f"modos={op['modos']}"], f3b + ".log")
        c, c2 = json.load(open(cj)), json.load(open(cj2))
        d2, d3 = json.load(open(f2)), json.load(open(f3))
        l1, l2 = json.load(open(os.path.join(p1, "ls.json"))), json.load(open(os.path.join(p2, "ls.json")))
        x = {"P1": l1["P_f"], "cov1": l1["CoV_linha"], "an1": l1["analises_linha"], "ess1": d2["direcao_ls"]["ESS"],
             "P2": l2["P_f"], "cov2": l2["CoV_linha"], "an2": l2["analises_linha"],
             "Pc": c["combinada"]["P_f"], "EPc": c["combinada"]["EP"],
             "Ps2": c2["combinada"]["P_f"], "EPs2": c2["combinada"]["EP"],
             "ang2": angulo(d2["u"], a), "ang3": angulo(d3["u"], a),
             "ang3b": angulo(json.load(open(f3b))["u"], a)}
        x["cov2_lin"] = cov_linear(x["ang2"], beta)
        res.append(x)
        print(f"r{r:2d}: est.1 P {x['P1']:.3e} CoV/linha {x['cov1']:5.2f} ESS {x['ess1']:5.1f} | dir. 2 a "
              f"{x['ang2']:5.1f} graus | est.2 P {x['P2']:.3e} CoV/linha {x['cov2']:5.2f} (linear {x['cov2_lin']:4.2f}) "
              f"an/linha {x['an2']:.2f} | comb. {x['Pc']:.3e} +- {x['EPc']:.1e} (z {(x['Pc'] - exato) / x['EPc']:+.2f})"
              f" | so 2: {x['Ps2']:.3e} +- {x['EPs2']:.1e} (z {(x['Ps2'] - exato) / x['EPs2']:+.2f}) | final a "
              f"{x['ang3']:5.1f} graus (so com o est. 2: {x['ang3b']:5.1f})", flush=True)
    R = len(res)
    resumo = {"op": op, "P_exata": exato, "repeticoes": res}
    for nome, kp, ke in (("combinada (pesos n1:n2)", "Pc", "EPc"), ("so o estagio 2", "Ps2", "EPs2")):
        P = [x[kp] for x in res]
        m = sum(P) / R
        sd = math.sqrt(sum((p - m) ** 2 for p in P) / (R - 1)) if R > 1 else float("nan")
        v_est = sum(x[ke] ** 2 for x in res) / R
        dentro = sum(abs(x[kp] - exato) < 2. * x[ke] for x in res)
        z = (m - exato) / (sd / math.sqrt(R)) if R > 1 and sd > 0 else float("nan")
        razao = v_est / sd ** 2 if R > 1 and sd > 0 else float("nan")
        print(f"{nome}: media {m:.4e} (exata {exato:.4e}), z da media {z:+.2f}; |P - exata| < 2 EP em {dentro}/{R}; "
              f"variancia estimada / empirica = {razao:.2f} (EP relativo da empirica ~{math.sqrt(2. / max(R - 1, 1)):.0%})",
              flush=True)
        resumo[kp] = {"media": m, "sd": sd, "z_media": z, "dentro_2EP": dentro, "razao_var": razao}
    for k, nome in (("ang2", "direcao do estagio 2"), ("ang3", "direcao final (estagios 1+2)"),
                    ("ang3b", "direcao final (so o estagio 2)")):
        v = sorted(x[k] for x in res)
        print(f"{nome}: angulo mediano {v[R // 2]:.1f}, maximo {v[-1]:.1f} graus; < 10 graus em "
              f"{sum(t < 10. for t in v)}/{R}", flush=True)
    print(f"CoV por linha do estagio 2: medido mediano {sorted(x['cov2'] for x in res)[R // 2]:.2f}, linear mediano "
          f"{sorted(x['cov2_lin'] for x in res)[R // 2]:.2f}; estagio 1: {sorted(x['cov1'] for x in res)[R // 2]:.2f}"
          f" (linear {cov_linear(op['ang'], beta):.2f})", flush=True)
    json.dump(resumo, open(os.path.join(pasta, "resumo.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
