#!/usr/bin/env python3
"""form_kl.py: FORM (HL-RF) do Exemplo 2 de Cho (2010) com os campos de c' e phi' truncados em M termos da KL cada.

As variaveis sao os 2M coeficientes xi da KL: normais padrao independentes (a correlacao cruzada r(c, phi) e a
transformacao lognormal ficam dentro do Cho2010), entao nao ha transformacao de Nataf. A funcao de estado limite e
    g(u) = FS(u) - 1,
com o FS do Cho2010 (SRM com refinamento adaptativo) nos xi dados (opcao xi=).
  - gradiente por diferencas centrais com passo h em u: 2n + 1 avaliacoes por iteracao (n = 2M), em paralelo. O FS tem
    ruido de uns 0,001 (malha adaptativa e patamar do arc-length): o passo e grande de proposito;
  - passo HL-RF u_{k+1} = (grad . u_k - g) / |grad|^2 grad, com amortecimento (metade do passo) quando |g| mais que dobra;
  - para quando |g| < tol_g (0,001) e o beta do ponto novo muda menos que tol_beta (0,01), ou no limite de iteracoes. O
    passo |u_{k+1} - u_k| nao serve de criterio: o ruido do FS nos termos fracos o mantem em ~0,1 (mais com M grande);
  - o beta corrigido, a distancia da origem ao plano tangente no ultimo ponto (|u| + g / |grad| com u alinhado ao
    gradiente: o ponto ainda nao esta em g = 0), vai junto no form.json.
Saida: <pasta>/form.json (beta, P_f = Phi(-beta), u*, alfa e alfa^2 por modo e por campo, o historico) e o historico na
tela.

Uso: python3 -I form_kl.py <exe Cho2010> <pasta de trabalho> [M=20] [rcphi=-0.5] [h=0.5] [paralelos=2] [iteracoes=10]
     [fator=tri2] [nref=4] [srm_f0=0.8] [tol_g=0.001] [tol_beta=0.01]"""
import concurrent.futures as cf
import csv
import json
import math
import os
import subprocess
import sys


def parametros(argv):
    p = {"M": 20, "rcphi": -0.5, "h": 0.5, "paralelos": 2, "iteracoes": 10, "fator": "tri2", "nref": 4, "srm_f0": 0.8,
         "tol_g": 0.001, "tol_beta": 0.01}
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
        os.makedirs(pasta, exist_ok=True)

    def fs(self, u, rotulo):
        """FS do Cho2010 com os xi u (os M de c' e os M de phi')"""
        arq = os.path.join(self.pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(x) for x in u) + "\n")
        cmd = [self.exe, "n=1", "primeira=0", f"xi={arq}", f"M={self.p['M']}", f"rcphi={self.p['rcphi']}",
               f"fator={self.p['fator']}", f"nref={self.p['nref']}", f"srm_f0={self.p['srm_f0']}", f"saida={rotulo}"]
        env = dict(os.environ, PZ_THREADS="0", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, cwd=self.pasta)
        with open(os.path.join(self.pasta, rotulo + "_0_0.csv"), newline="") as f:
            linha = next(csv.DictReader(f))
        if linha["situacao"] != "ok":
            raise RuntimeError(f"{rotulo}: situacao {linha['situacao']}")
        return float(linha["FS_min"])

    def g_e_gradiente(self, u, it):
        h = self.p["h"]
        pontos = [("centro", list(u))]
        for i in range(len(u)):
            for sinal in (+1, -1):
                v = list(u)
                v[i] += sinal * h
                pontos.append((f"d{i}{'p' if sinal > 0 else 'm'}", v))
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
    M = p["M"]
    u = [0.] * (2 * M)
    hist = []
    print(f"FORM do Cho (2010), Exemplo 2: KL com M = {M} termos por campo ({2 * M} variaveis), r(c, phi) = {p['rcphi']}, "
          f"passo h = {p['h']} em u, {2 * 2 * M + 1} avaliacoes por iteracao, {p['paralelos']} em paralelo", flush=True)
    print("  it   g = FS - 1    beta     |grad|   maiores alfa^2 (campo, termo)", flush=True)
    gAnterior = None
    convergiu = False
    for it in range(p["iteracoes"]):
        g, grad, fs = prob.g_e_gradiente(u, it)
        beta = math.sqrt(sum(x * x for x in u))
        n2 = sum(x * x for x in grad)
        alfa2 = [x * x / n2 for x in grad] if n2 > 0 else [0.] * len(grad)
        top = sorted(range(len(grad)), key=lambda i: -alfa2[i])[:4]
        print(f"  {it:2d}   {g:+.5f}    {beta:.4f}   {math.sqrt(n2):.4f}   " +
              ", ".join(f"{'c' if i < M else 'phi'}{i % M}: {alfa2[i]:.2f}" for i in top), flush=True)
        hist.append({"it": it, "u": u, "g": g, "grad": grad, "beta": beta, "fs_pontos": fs})
        with open(os.path.join(pasta, "form_parcial.json"), "w") as f:
            json.dump({"parametros": p, "historico": hist, "avaliacoes": prob.n}, f, indent=1)
        gu = sum(gi * ui for gi, ui in zip(grad, u))
        novo = [(gu - g) / n2 * gi for gi in grad]
        beta_novo = math.sqrt(sum(x * x for x in novo))
        if it > 0 and abs(g) < p["tol_g"] and abs(beta_novo - beta) < p["tol_beta"]:
            convergiu = True
            break
        if gAnterior is not None and abs(g) > 2. * abs(gAnterior) and it > 1:
            novo = [ui + 0.5 * (ni - ui) for ui, ni in zip(u, novo)]  # amortece: metade do passo
        gAnterior = g
        u = novo
    beta = math.sqrt(sum(x * x for x in u))
    n = math.sqrt(sum(x * x for x in grad))
    alfa = [-x / n for x in grad]
    a2 = [x * x for x in alfa]
    beta_corr = beta_novo  # distancia ao plano tangente no ultimo ponto avaliado (o proximo passo do HL-RF)
    res = {"parametros": p, "beta": beta, "P_f": phi_normal(-beta), "beta_corrigido": beta_corr,
           "P_f_corrigida": phi_normal(-beta_corr), "convergiu": convergiu, "iteracoes": len(hist), "u": u, "alfa": alfa, "alfa2": a2,
           "alfa2_c": sum(a2[:M]), "alfa2_phi": sum(a2[M:]), "avaliacoes": prob.n, "historico": hist}
    with open(os.path.join(pasta, "form.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"\n  beta = {beta:.4f}  P_f = Phi(-beta) = {100 * res['P_f']:.3f}%  ({prob.n} avaliacoes do FS)")
    print(f"  beta corrigido = {beta_corr:.4f}  P_f = {100 * res['P_f_corrigida']:.3f}%")
    print("  convergiu" if convergiu else f"  NAO convergiu em {p['iteracoes']} iteracoes")
    print(f"  importancia por campo: c' {res['alfa2_c']:.3f}, phi' {res['alfa2_phi']:.3f}")


if __name__ == "__main__":
    main()
