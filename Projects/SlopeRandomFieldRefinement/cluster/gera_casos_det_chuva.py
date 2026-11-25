#!/usr/bin/env python3
"""gera_casos_det_chuva.py: a lista de casos da validacao deterministica da chuva (casos_det_chuva.txt), para o
det_chuva.sbatch (jobs detChuva). O plano: Projects/PlanoDeTrabalho/chuva2/plano_deterministico.md.

Todos os casos: o deterministico (as medias, det=1 n=0, sem campos: covc=covphi=covk=0) do SlopeRainMonteCarlo nas tres
variantes do transiente (iii desacoplado e ar passivo, ii acoplado e ar passivo, i acoplado e ar ativo), com o FS (SRM,
adaptativo do zero em cada instante) de hora em hora, o diagnostico do transiente em cada instante (diag=diag.csv:
passos, Newton, nos na face, p_w e p_a extremos, |u| max, agua infiltrada) e as pressoes nos nos (grava_pressoes=pw.txt,
com p_a na variante i). Uma linha por caso: <id> <opcoes>. Na ordem de prioridade:
  T1  as chuvas da Serie B do Rahardjo (24 h) nos tres solos
  T2  o mesmo volume (864 mm) em 3, 6, 12, 24 e 72 h nos tres solos (o A6 do relatorio de outubro), mais as duracoes
      criticas t* = V / k_s fora dessa faixa (f10 em 1 h, f100 em 240 h)
  T3  a sensibilidade em torno da chuva base de cada solo (f10 90, f50 36, f100 9 mm/h; 24 h): c', phi' e k_s em 3
      niveis (baixo, base, alto), a anisotropia de k (k_h/k_v com o tensor diagonal), a rotacao do tensor e as duas
      juntas
  T4  as verificacoes numericas e de convencao nas chuvas base: malha do transiente (tri3, tri4, camadas de 0,25 m ate
      5 m), passo de tempo (0,05 h), malha do fator (nref 5) e a chuva na projecao horizontal (chuva_face=vertical)
  T5  as series do relatorio de outubro e do artigo (figuras e contraexemplos): Serie A e A5 (angulo), Serie C
      (lencol), Serie D (altura), A2 (f100 a 63,4 graus e com H_s = 20 m), A3 (tensor girado em todas as chuvas) e A4
      (teto da succao inicial)
Casos repetidos (as mesmas opcoes) ficam so uma vez; a linha '#= <id> <id do caso que roda>' guarda o apelido.
Uso: python3 gera_casos_det_chuva.py > casos_det_chuva.txt
"""
import math

SOLOS = {"f10": 1e-4, "f50": 1e-5, "f100": 1e-6}  # k_s (m/s)
VARIANTES = ["iii", "ii", "i"]
BASE = {"f10": 90, "f50": 36, "f100": 9}  # as chuvas base da sensibilidade (24 h), perto da ruptura na Serie B
SERIE_B = {"f10": [9, 36, 54, 72, 90, 120, 180, 360, 900],
           "f50": [3.6, 9, 18, 27, 36, 45, 54, 90, 180],
           "f100": [0.9, 3.6, 9, 18, 36, 80]}
COMUM = "det=1 n=0 covc=0 covphi=0 covk=0 saida=det diag=diag.csv grava_pressoes=pw.txt"


def num(x):
    """o numero sem zeros inuteis (para os ids e as opcoes)"""
    s = f"{x:.6g}"
    return s


def tempos(D):
    """os instantes do FS (h): de hora em hora ate max(48, D + 24); D <= 6: de meia em meia hora ate D + 6"""
    fim = max(48, D + 24)
    ts = set()
    if D <= 6:
        ts |= {k * 0.5 for k in range(int(2 * (D + 6)) + 1)}
    passo = 2 if fim > 120 else 1
    ts |= set(range(0, int(fim) + 1, passo))
    return ",".join(num(t) for t in sorted(ts))


def geometria(alfa_graus, H):
    """o talude do Rahardjo (Fig. 2): base em 0, pe em 2H, crista em 3H, platos de 3H; face com inclinacao alfa"""
    b = H / math.tan(math.radians(alfa_graus))
    x2 = 3 * H + b
    x1 = x2 + 3 * H
    pts = [0, 0, x1, 0, x1, 2 * H, x2, 2 * H, 3 * H, 3 * H, 0, 3 * H]
    return "geometria=" + ",".join(num(round(v, 4)) for v in pts)


casos = []  # (id, opcoes)
vistos = {}  # opcoes -> id
apelidos = []


def caso(cid, opc):
    chave = " ".join(sorted(opc.split()))
    if chave in vistos:
        apelidos.append((cid, vistos[chave]))
        return
    vistos[chave] = cid
    casos.append((cid, opc))


def base_opc(solo, I, D=24, extra=""):
    return f"solo={solo} I={num(I)} chuva={num(D)} tempos={tempos(D)} {extra}".strip()


# T1: a Serie B
for solo, Is in SERIE_B.items():
    for I in Is:
        for v in VARIANTES:
            caso(f"T1_{solo}_I{num(I)}_{v}", f"{base_opc(solo, I)} variante={v}")

# T2: o mesmo volume em duracoes diferentes
V = 864.0  # mm (36 mm/h x 24 h, o A6)
for solo in SOLOS:
    Ds = [3, 6, 12, 24, 72] + ({"f10": [1], "f100": [240]}.get(solo, []))
    for D in Ds:
        for v in VARIANTES:
            caso(f"T2_{solo}_D{num(D)}_{v}", f"{base_opc(solo, V / D, D)} variante={v}")

# T3: a sensibilidade em torno da chuva base
for solo, ks in SOLOS.items():
    I = BASE[solo]
    niveis = [("c5", "c=5"), ("c15", "c=15"), ("phi21", "phi=21"), ("phi31", "phi=31"),
              ("ksd3", f"ks={num(ks / 3)}"), ("ksx3", f"ks={num(ks * 3)}"),
              ("a2", "alfa=2"), ("a5", "alfa=5"), ("a10", "alfa=10"),
              ("a5t45", "alfa=5 k_theta=45"), ("a5t90", "alfa=5 k_theta=90"), ("a5t135", "alfa=5 k_theta=135"),
              ("a2t135", "alfa=2 k_theta=135"), ("a10t135", "alfa=10 k_theta=135")]
    for nome, extra in niveis:
        for v in VARIANTES:
            caso(f"T3_{solo}_{nome}_{v}", f"{base_opc(solo, I, 24, extra)} variante={v}")

# T4: verificacoes numericas e de convencao nas chuvas base
for solo in SOLOS:
    I = BASE[solo]
    for nome, extra in [("tri3", "transiente=tri3"), ("tri4", "transiente=tri4"), ("cam", "transiente=cam0.25p5"),
                        ("dt005", "dt=0.05"), ("nref5", "nref=5"), ("vert", "chuva_face=vertical")]:
        for v in VARIANTES:
            caso(f"T4_{solo}_{nome}_{v}", f"{base_opc(solo, I, 24, extra)} variante={v}")

# T5: as series do relatorio e do artigo
PARES_A = {"f10": [9, 80, 360], "f50": [9, 36, 80], "f100": [3.6, 9, 80]}
for solo, Is in PARES_A.items():  # Serie A
    for a in [26.565, 33.690, 45, 63.435]:
        for I in Is:
            for v in VARIANTES:
                caso(f"T5A_{solo}_a{num(round(a, 1))}_I{num(I)}_{v}",
                     f"{base_opc(solo, I, 24, geometria(a, 10) if a != 45 else '')} variante={v}")
for solo, I in [("f10", 80), ("f50", 36), ("f100", 9)]:  # A5: o angulo de 20 a 70 graus
    for a in range(20, 71, 5):
        for v in VARIANTES:
            caso(f"T5A5_{solo}_a{a}_I{num(I)}_{v}",
                 f"{base_opc(solo, I, 24, geometria(a, 10) if a != 45 else '')} variante={v}")
for hw in [2.5, 5, 7.5, 10, 15]:  # Serie C (f50)
    for I in [9, 36, 80]:
        for v in VARIANTES:
            caso(f"T5C_f50_hw{num(hw)}_I{num(I)}_{v}", f"{base_opc('f50', I, 24, f'hw={num(hw)}')} variante={v}")
for H in [5, 10, 20, 40]:  # Serie D (f50; hw = 5 m)
    for I in [9, 36, 80]:
        for v in VARIANTES:
            caso(f"T5D_f50_H{num(H)}_I{num(I)}_{v}",
                 f"{base_opc('f50', I, 24, geometria(45, H) if H != 10 else '')} variante={v}")
for I in [9, 18, 36, 80]:  # A2: o f100 a 63,4 graus (G1) e com H_s = 20 m (G2)
    for v in VARIANTES:
        caso(f"T5A2G1_f100_I{num(I)}_{v}", f"{base_opc('f100', I, 24, geometria(63.435, 10))} variante={v}")
        caso(f"T5A2G2_f100_I{num(I)}_{v}", f"{base_opc('f100', I, 24, geometria(45, 20))} variante={v}")
for t in [0, 45, 135]:  # A3: k_1 = 10 k_s girado, f50, todas as chuvas
    for I in [9, 18, 36, 54, 90, 180]:
        for v in VARIANTES:
            extra = "alfa=10" + (f" k_theta={t}" if t else "")
            caso(f"T5A3_f50_t{t}_I{num(I)}_{v}", f"{base_opc('f50', I, 24, extra)} variante={v}")
for smax in [75, 150, 300]:  # A4: o teto da succao inicial (f50)
    for hw in [2.5, 5, 7.5, 10, 15]:
        for I in [9, 36]:
            for v in VARIANTES:
                extra = f"hw={num(hw)}" + (f" succao_max={smax}" if smax != 75 else "")
                caso(f"T5A4_f50_s{smax}_hw{num(hw)}_I{num(I)}_{v}", f"{base_opc('f50', I, 24, extra)} variante={v}")

print(f"# casos_det_chuva.txt: gerado por gera_casos_det_chuva.py ({len(casos)} casos, {len(apelidos)} apelidos). "
      "Opcoes comuns (o det_chuva.sbatch poe):")
print(f"#   {COMUM}")
for cid, opc in casos:
    print(f"{cid} {opc}")
for a, b in apelidos:
    print(f"#= {a} {b}")
