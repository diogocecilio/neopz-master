#!/usr/bin/env python3
"""form_kl.py: FORM (HL-RF) do Exemplo 2 de Cho (2010) com os campos de c' e phi' truncados em M termos da KL cada.

As variaveis sao os 2M coeficientes xi da KL: normais padrao independentes (a correlacao cruzada r(c, phi) e a
transformacao lognormal ficam dentro do Cho2010), entao nao ha transformacao de Nataf. A funcao de estado limite e
    g(u) = FS(u) - 1   (metodo=srm, o padrao)   ou   g(u) = Gamma(u) - 1   (metodo=gim),
com o fator do Cho2010 (SRM ou GIM, com refinamento adaptativo) nos xi dados (opcao xi=). Os dois tem a mesma regiao de
falha (FS < 1 se e so se Gamma < 1), entao o mesmo ponto de projeto; o g e o gradiente e que mudam.
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
     [fator=tri2] [nref=4] [srm_f0=0.8] [tol_g=0.001] [tol_beta=0.01] [gradiente=dif|modal] [nystrom=<base>]
     [metodo=srm|gim]
prog=rebaixamento: o rebaixamento do Ceron (SlopeDrawdownMonteCarlo, Projects/SlopeRandomFieldRefinement), com 3 campos
(c', phi' e k: 3M variaveis) e g = Gamma - 1 (GIM); as opcoes rcphi, srm_f0, metodo e nystrom nao valem, e transiente=
escolhe a malha do transiente. O gradiente modal dele da 0 nos termos de k (a analise limite congela as pressoes).
As opcoes do rebaixamento contorno= (reservatorio|ceron), agua= (unsat|sat), trebaix= (h), tfim= (h; -1: = trebaix),
dt= (h) e hw= (h_w/H) e as dos casos das Secoes 6.2 a 6.4 do artigo, covc=, covphi=, covk= (CoV de c', phi' e k_v),
s= (escala dos comprimentos de autocorrelacao: 20 s e 2 s m) e alfa= (k_h/k_v), tem os padroes do
SlopeDrawdownMonteCarlo (os de antes) e so vao para o executavel quando diferentes do padrao; vao no form.json, e o
line_sampling.py usa as mesmas. Ex.: a agua do artigo (percolacao permanente depois do rebaixamento instantaneo):
contorno=ceron agua=sat trebaix=0.0001 tfim=20 dt=5.
gradiente=modal: uma analise por iteracao, com o gradiente pela projecao modal da sensibilidade de limite superior
(Cho2010 gradiente=M; no colapso do SRM ou do GIM, conforme o metodo), em vez das 4M + 1 das diferencas centrais."""
import concurrent.futures as cf
import csv
import json
import math
import os
import subprocess
import sys

# as opcoes do rebaixamento passadas ao SlopeDrawdownMonteCarlo, com os padroes dele (so as diferentes vao na linha)
REBAIX_PADRAO = {"contorno": "reservatorio", "agua": "unsat", "trebaix": 0.1, "tfim": -1.0, "dt": 0.1, "hw": 1.0,
                 "covc": 0.3, "covphi": 0.1, "covk": 0.6, "s": 1.0, "alfa": 1.0}


def parametros(argv):
    p = {"M": 20, "rcphi": -0.5, "h": 0.5, "paralelos": 2, "iteracoes": 10, "fator": "tri2", "nref": 4, "srm_f0": 0.8,
         "tol_g": 0.001, "tol_beta": 0.01, "gradiente": "dif", "nystrom": "", "metodo": "srm", "prog": "cho",
         "transiente": "tri3", **REBAIX_PADRAO}
    for a in argv:
        k, v = a.split("=", 1)
        if k not in p:
            sys.exit(f"opcao desconhecida: {a}")
        p[k] = type(p[k])(v)
    if p["prog"] not in ("cho", "rebaixamento"):
        sys.exit("prog=cho ou prog=rebaixamento")
    if p["prog"] == "rebaixamento":
        p["metodo"] = "gim"  # o fator do rebaixamento e o Gamma do GIM
    return p


def campos(p):
    """os nomes dos campos (a ordem dos xi: os M de cada um)"""
    return ["c", "phi", "k"] if p["prog"] == "rebaixamento" else ["c", "phi"]


def comando(exe, p, arq, rotulo):
    """a linha de comando de uma analise nos xi do arquivo"""
    if p["prog"] == "rebaixamento":
        extra = [f"{k}={p[k]}" for k, v in REBAIX_PADRAO.items() if k in p and p[k] != v]
        return [exe, "n=1", "primeira=0", f"xi={arq}", f"transiente={p['transiente']}", f"fator={p['fator']}",
                f"nref={p['nref']}"] + extra + [f"saida={rotulo}"]
    return [exe, "n=1", "primeira=0", f"xi={arq}", f"M={p['M']}", f"rcphi={p['rcphi']}", f"fator={p['fator']}",
            f"nref={p['nref']}", f"srm_f0={p['srm_f0']}", f"metodo={p['metodo']}", f"saida={rotulo}"] + \
        ([f"nystrom={p['nystrom']}"] if p["nystrom"] else [])


class Problema:
    def __init__(self, exe, pasta, p):
        self.exe, self.pasta, self.p = exe, pasta, p
        self.n = 0
        os.makedirs(pasta, exist_ok=True)

    def fs(self, u, rotulo):
        """o fator do Cho2010 (FS ou Gamma) com os xi u (os M de c' e os M de phi')"""
        arq = os.path.join(self.pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(x) for x in u) + "\n")
        cmd = comando(self.exe, self.p, arq, rotulo)
        env = dict(os.environ, PZ_THREADS="0", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, cwd=self.pasta)
        with open(os.path.join(self.pasta, rotulo + "_0_0.csv"), newline="") as f:
            linha = next(csv.DictReader(f))
        if linha["situacao"] != "ok":
            raise RuntimeError(f"{rotulo}: situacao {linha['situacao']}")
        return float(linha["FS_min"])

    def g_modal(self, u, it):
        """uma analise: o fator e o gradiente pela projecao modal (Cho2010 gradiente=M: limite superior no colapso)"""
        rotulo = f"it{it}_centro"
        arq = os.path.join(self.pasta, rotulo + ".xi")
        with open(arq, "w") as f:
            f.write("\n".join(repr(x) for x in u) + "\n")
        cmd = comando(self.exe, self.p, arq, rotulo) + [f"gradiente={self.p['M']}"]
        env = dict(os.environ, PZ_THREADS="0", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, cwd=self.pasta)
        with open(os.path.join(self.pasta, rotulo + "_0_0.csv"), newline="") as f:
            linha = next(csv.DictReader(f))
        if linha["situacao"] != "ok":  # sem o _grad_0.txt: a mensagem clara em vez de FileNotFoundError
            raise RuntimeError(f"{rotulo}: situacao {linha['situacao']} (FS_min {linha['FS_min']}), sem gradiente")
        v = [float(x) for x in open(os.path.join(self.pasta, rotulo + "_grad_0.txt")).read().split("\n")[1:] if x.strip()]
        self.n += 1
        return v[0] - 1., v[1:], [v[0]]

    def g_e_gradiente(self, u, it):
        if self.p["gradiente"] == "modal":
            return self.g_modal(u, it)
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
    nomes = campos(p)
    nc = len(nomes)
    u = [0.] * (nc * M)
    hist = []
    titulo = ("rebaixamento do Ceron (2025), " + (" ".join(comando("", p, "", "")[4:-1]) or "padrao")) \
        if p["prog"] == "rebaixamento" else \
        f"Cho (2010), Exemplo 2, r(c, phi) = {p['rcphi']}"
    tipo = "modal (1 avaliacao por iteracao)" if p["gradiente"] == "modal" else \
        f"por diferencas centrais (h = {p['h']}, {2 * nc * M + 1} avaliacoes por iteracao)"
    print(f"FORM do {titulo}, {p['metodo'].upper()}: KL com M = {M} termos por campo ({nc * M} variaveis), gradiente "
          f"{tipo}, {p['paralelos']} em paralelo", flush=True)
    print(f"  it   g = {'Gamma' if p['metodo'] == 'gim' else 'FS'} - 1    beta     |grad|   maiores alfa^2 (campo, termo)",
          flush=True)
    gAnterior = None
    convergiu = False
    for it in range(p["iteracoes"]):
        g, grad, fs = prob.g_e_gradiente(u, it)
        beta = math.sqrt(sum(x * x for x in u))
        n2 = sum(x * x for x in grad)
        alfa2 = [x * x / n2 for x in grad] if n2 > 0 else [0.] * len(grad)
        top = sorted(range(len(grad)), key=lambda i: -alfa2[i])[:4]
        print(f"  {it:2d}   {g:+.5f}    {beta:.4f}   {math.sqrt(n2):.4f}   " +
              ", ".join(f"{nomes[i // M]}{i % M}: {alfa2[i]:.2f}" for i in top), flush=True)
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
           "alfa2_c": sum(a2[:M]), "alfa2_phi": sum(a2[M:2 * M]), "alfa2_k": sum(a2[2 * M:]), "avaliacoes": prob.n,
           "historico": hist}
    with open(os.path.join(pasta, "form.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"\n  beta = {beta:.4f}  P_f = Phi(-beta) = {100 * res['P_f']:.3f}%  ({prob.n} avaliacoes do FS)")
    print(f"  beta corrigido = {beta_corr:.4f}  P_f = {100 * res['P_f_corrigida']:.3f}%")
    print("  convergiu" if convergiu else f"  NAO convergiu em {p['iteracoes']} iteracoes")
    print(f"  importancia por campo: c' {res['alfa2_c']:.3f}, phi' {res['alfa2_phi']:.3f}"
          + (f", k {res['alfa2_k']:.3f}" if nc == 3 else ""))


if __name__ == "__main__":
    main()
