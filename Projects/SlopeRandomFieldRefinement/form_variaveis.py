#!/usr/bin/env python3
"""form_variaveis.py: FORM (HLRF) da chuva com c', phi' e k_s como variaveis aleatorias lognormais (solo homogeneo em
cada avaliacao: o limite s -> infinito dos campos aleatorios), com o FS do SlopeRainMonteCarlo (opcao valores=).

Variaveis no espaco normal padrao: u_i = (ln x_i - mu_ln_i) / sigma_ln_i, com sigma_ln = sqrt(ln(1 + CoV^2)) e
mu_ln = ln(media) - sigma_ln^2 / 2 (as mesmas lognormais dos campos). Funcao de estado limite: g(u) = FS(24 h) - 1.
Gradiente por diferencas centrais com passo h em u (o FS tem ruido de ~0,001: o patamar do arc-length e a malha
adaptativa), as 2n + 1 avaliacoes de cada iteracao em paralelo. Passo HLRF:
  u_{k+1} = (grad . u_k - g) / |grad|^2 grad
com amortecimento (metade do passo, ate 3 vezes) quando |g| piora muito. Para quando |u_{k+1} - u_k| < tol e
|g| < 0,005.

Uso: python3 form_variaveis.py <exe SlopeRainMonteCarlo> <pasta de trabalho> [I=31.05] [covc=0.3] [covphi=0.1]
     [covk=0.6] [solo=f50] [h=0.2] [paralelos=7] [threads=2]
Saida: historico na tela e <pasta>/form.json (beta, P_f = Phi(-beta), ponto de projeto, cossenos diretores alfa).
"""
import concurrent.futures as cf
import csv
import json
import math
import os
import subprocess
import sys

MEDIAS = {"c": 10.0, "phi": 26.0, "k": 1.0}  # c' (kPa), phi' (graus), k / k_s


def parametros(argv):
    p = {"I": 31.05, "covc": 0.3, "covphi": 0.1, "covk": 0.6, "solo": "f50", "h": 0.2, "paralelos": 7, "threads": 2}
    for a in argv:
        k, v = a.split("=", 1)
        if k not in p:
            sys.exit(f"opcao desconhecida: {a}")
        p[k] = type(p[k])(v)
    return p


class Problema:
    def __init__(self, exe, pasta, p):
        self.exe, self.pasta, self.p = exe, pasta, p
        self.n = 0
        self.ln = {}
        for nome, cov in (("c", p["covc"]), ("phi", p["covphi"]), ("k", p["covk"])):
            s = math.sqrt(math.log(1. + cov * cov))
            self.ln[nome] = (math.log(MEDIAS[nome]) - s * s / 2., s)
        os.makedirs(pasta, exist_ok=True)

    def x(self, u):
        return [math.exp(self.ln[nome][0] + self.ln[nome][1] * ui) for nome, ui in zip(("c", "phi", "k"), u)]

    def fs(self, u, rotulo):
        """FS(24 h) do solo homogeneo com os valores de u (uma rodada do SlopeRainMonteCarlo)"""
        c, phi, k = self.x(u)
        saida = os.path.join(self.pasta, rotulo)
        cmd = [self.exe, "n=0", "det=1", f"I={self.p['I']}", f"solo={self.p['solo']}", "tempos=24",
               f"valores={c:.6g},{phi:.6g},{k:.6g}", f"saida={saida}"]
        env = dict(os.environ, PZ_THREADS=str(self.p["threads"]))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, cwd=self.pasta)
        with open(saida + "_0_-1.csv", newline="") as f:
            linhas = [r for r in csv.DictReader(f) if r["t_h"] in ("24", "24.0")]
        return float(linhas[0]["FS_min"])

    def g_e_gradiente(self, u, it):
        h = self.p["h"]
        pontos = [("centro", list(u))]
        for i in range(len(u)):
            for sinal in (+1, -1):
                v = list(u)
                v[i] += sinal * h
                pontos.append((f"d{i}{'+' if sinal > 0 else '-'}", v))
        with cf.ThreadPoolExecutor(max_workers=self.p["paralelos"]) as ex:
            fs = list(ex.map(lambda pt: self.fs(pt[1], f"it{it}_{pt[0]}"), pontos))
        self.n += len(pontos)
        g = fs[0] - 1.
        grad = [(fs[1 + 2 * i] - fs[2 + 2 * i]) / (2. * h) for i in range(len(u))]
        return g, grad, fs


def phi_normal(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    exe, pasta = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    p = parametros(sys.argv[3:])
    prob = Problema(exe, pasta, p)
    u = [0., 0., 0.]
    hist = []
    print(f"FORM com variaveis aleatorias: solo {p['solo']}, I = {p['I']} mm/h, CoV c' {p['covc']}, phi' {p['covphi']}, "
          f"k {p['covk']}; passo h = {p['h']} em u")
    print("  it   g = FS - 1    beta    u (c', phi', k)               x (c' kPa, phi' graus, k/k_s)")
    gAnterior = None
    for it in range(15):
        g, grad, fs = prob.g_e_gradiente(u, it)
        beta = math.sqrt(sum(ui * ui for ui in u))
        x = prob.x(u)
        print(f"  {it:2d}   {g:+.5f}    {beta:.4f}  ({u[0]:+.3f}, {u[1]:+.3f}, {u[2]:+.3f})   "
              f"({x[0]:.3f}, {x[1]:.3f}, {x[2]:.3f})   grad = ({grad[0]:+.4f}, {grad[1]:+.4f}, {grad[2]:+.4f})")
        hist.append({"it": it, "u": u, "x": x, "g": g, "grad": grad, "beta": beta, "fs_pontos": fs})
        n2 = sum(gi * gi for gi in grad)
        gu = sum(gi * ui for gi, ui in zip(grad, u))
        novo = [(gu - g) / n2 * gi for gi in grad]
        passo = math.sqrt(sum((a - b) ** 2 for a, b in zip(novo, u)))
        if it > 0 and abs(g) < 5.e-3 and passo < 2.e-2:
            break
        if gAnterior is not None and abs(g) > 2. * abs(gAnterior) and it > 1:
            novo = [ui + 0.5 * (ni - ui) for ui, ni in zip(u, novo)]  # amortece: metade do passo
        gAnterior = g
        u = novo
    beta = math.sqrt(sum(ui * ui for ui in u))
    n = math.sqrt(sum(gi * gi for gi in grad))
    alfa = [-gi / n for gi in grad]  # cossenos diretores (do lado da falha)
    res = {"parametros": p, "beta": beta, "P_f": phi_normal(-beta), "u": u, "x": prob.x(u), "alfa": alfa,
           "alfa2": [a * a for a in alfa], "avaliacoes": prob.n, "historico": hist}
    with open(os.path.join(pasta, "form.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"\n  beta = {beta:.4f}  P_f = Phi(-beta) = {100 * res['P_f']:.2f}%  ({prob.n} avaliacoes do FS)")
    print(f"  ponto de projeto: c' = {res['x'][0]:.3f} kPa, phi' = {res['x'][1]:.3f} graus, k/k_s = {res['x'][2]:.3f}")
    print("  importancia alfa^2: c' %.3f, phi' %.3f, k %.3f" % tuple(res["alfa2"]))


if __name__ == "__main__":
    main()
