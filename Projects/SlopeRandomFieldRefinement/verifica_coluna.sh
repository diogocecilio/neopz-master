#!/bin/bash
# verifica_coluna.sh - as verificacoes do SlopeRainColumn (coluna.cpp), da secao C da spec dos passos 1 a 3:
#   v01  hidrostatico parado (I=0 nivel0=1): |p_w(t) - p_w(0)| < 1e-8 kPa (UnsatWaterPhase 2.1: 1,8e-15)
#   v03  simetria com campos (n=5, tri e quad): max |p_L - p_R| <= 1e-8 kPa e nos_face em {0, 2}
#   v04  balanco de agua (I=5.4 e 10.8): |dW - infiltrado| <= 1e-4 infiltrado
#   v05  parada do transiente: a realizacao 9 do caso base (lv 0,5) em tri, ny 80, com o Newton antigo 1e-6,40,4 para
#        em ~20,83 h; o ultimo sub-passo aceito entra na travessia: t_passagem em [0, t_parou] (~20,80 h) e FS_passagem < 1
#   v12  FS em forma fechada no estado inicial: succao 20 kPa (1,5065/z + 0,7447) e coluna cheia (0,5204/z + 0,1067),
#        FS(2) = 1,4979 e 0,3669 +- 1e-3
#   v13  o FS recalculado pelo 'coluna_ref.py fs' no perfil gravado: 1e-10
#   v21  o deterministico do Cho (Fig. 3) em ny 40/80/200/400 x dt 0.1/0.05/0.025 x tri/quad x busca 4/0, contra o
#        'coluna_ref.py det N=400 dt=0.025' nos mesmos pontos; na malha mais fina, FS +-0,003 (6 e 24 h) e +-0,01 (12 e
#        18 h), frente +-1 elemento, t_passagem +-0,05 h; as ordens observadas
#   v31  a Fig. 3 na malha de producao (NY, dt 0.1, EL): 1,50/1,47/1,06/0,37 +-0,02 em 0/6/18/24 h, 12 h em [1,22; 1,27],
#        frentes 0,60-0,65 / 1,15-1,25 / 1,75-1,85 m (usa a rodada do v21 se houver; senao, roda)
#   v23  custo: n=20 com o campo do caso base em ny 40/80/200 x pardiso/skyline x tri/quad (s por realizacao)
# Uso: verifica_coluna.sh <pasta> [testes...]   (padrao: v01 v03 v04 v05 v12 v13; 'todos': os nove)
# Variaveis: EXE (o SlopeRainColumn; padrao, o do build release), PY (python3), P (processos em paralelo; 8),
#            NY e EL (a malha de producao do v31; 80 e quad), PC (processos do v23; 4: o custo e medido com metade dos nucleos)
# Cada teste grava em <pasta>/<teste>/ e escreve 'Vxx passou|falhou: ...' em <pasta>/resumo.txt
set -u
AQUI="$(cd "$(dirname "$0")" && pwd)"
EXE="${EXE:-$HOME/projects/neopz-master-build-release/Projects/SlopeRandomFieldRefinement/SlopeRainColumn}"
PY="${PY:-python3}"
P="${P:-8}"
PC="${PC:-4}"
NY="${NY:-80}"
EL="${EL:-quad}"
REF="$AQUI/referencias_1d/coluna_ref.py"
NYS="$AQUI/../Cho2010/nystrom_base.py"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
[ $# -ge 1 ] || { sed -n '2,20p' "$0"; exit 1; }
PASTA="$(mkdir -p "$1" && cd "$1" && pwd)"
shift
TESTES="${*:-v01 v03 v04 v05 v12 v13}"
[ "$TESTES" = todos ] && TESTES="v01 v03 v04 v05 v12 v13 v21 v31 v23"
[ -x "$EXE" ] || { echo "sem o executavel: $EXE"; exit 1; }
RESUMO="$PASTA/resumo.txt"
echo "# $(date '+%F %T') $EXE ($(md5sum < "$EXE" | cut -c1-8)), testes: $TESTES" >> "$RESUMO"
anota() { echo "$1" | tee -a "$RESUMO"; }

# a base de Nystrom 1D do caso base (lv = 0,5 m, M = 100), uma vez
base() {
    local b="$PASTA/bases/base_lv$1_M100.txt"
    if [ ! -s "$b" ]; then
        mkdir -p "$PASTA/bases"
        "$PY" -I "$NYS" "$b" coluna=2 lv="$1" dy=0.0025 M=100 > "${b%.txt}.log" 2>&1 || { echo "falhou: $b"; exit 1; }
    fi
    echo "$b"
}
# roda comandos (um por linha de um arquivo) com P processos
paralelo() { xargs -d '\n' -P "$1" -I{} bash -c '{}' < "$2"; }

for T in $TESTES; do
    D="$PASTA/$T"
    mkdir -p "$D" && cd "$D" || exit 1
    case $T in
    v01)
        "$EXE" det=1 n=0 I=0 nivel0=1 tempos=0,1,10 saida=v01 solver=skyline > v01.log 2>&1
        anota "$("$PY" -I - <<'EOF'
import csv, collections
d = collections.defaultdict(dict)
for r in csv.DictReader(open('v01_perfil_-1.csv')):
    d[float(r['t_h'])][float(r['z'])] = float(r['pw'])
e = max(abs(d[t][z] - d[0.][z]) for t in (1., 10.) for z in d[0.])
print(f"V0.1 {'passou' if e < 1e-8 else 'falhou'}: hidrostatico (nivel em y = 1 m, I = 0) em 1 e 10 h: "
      f"max |p_w - p_w0| = {e:.2e} kPa (< 1e-8; UnsatWaterPhase 2.1: 1,8e-15)")
EOF
)" ;;
    v03)
        B="$(base 0.5)"
        for e in tri quad; do
            echo "\"$EXE\" n=5 nystrom=$B elementos=$e saida=v03_$e solver=skyline > v03_$e.log 2>&1"
        done > cmds.txt
        paralelo 2 cmds.txt
        anota "$("$PY" -I - <<'EOF'
import csv, re
ok, txt = True, []
for e in ('tri', 'quad'):
    s = [float(x) for x in re.findall(r'max \|p_L - p_R\| = (\S+) kPa', open(f'v03_{e}.log').read())]
    faces = {int(r['nos_face']) for r in csv.DictReader(open(f'v03_{e}_0_4.csv'))}
    ok = ok and len(s) == 5 and max(s) <= 1e-8 and faces <= {0, 2}
    txt.append(f"{e}: max |p_L - p_R| = {max(s):.1e} kPa, nos_face {sorted(faces)}")
print(f"V0.3 {'passou' if ok else 'falhou'}: simetria com 5 campos (lv 0,5, CoV 0,5): " + "; ".join(txt))
EOF
)" ;;
    v04)
        for I in 5.4 10.8; do
            "$EXE" det=1 n=0 I=$I saida=v04_$I solver=skyline > v04_$I.log 2>&1
        done
        anota "$("$PY" -I - <<'EOF'
import re
ok, txt = True, []
for I in ('5.4', '10.8'):
    m = re.search(r'\|dW - infiltrado\| = (\S+) m3/m \((\S+) do infiltrado (\S+) m3/m\)', open(f'v04_{I}.log').read())
    ok = ok and m is not None and float(m.group(2)) <= 1e-4
    txt.append(f"I = {I}: |dW - infiltrado| = {m.group(1)} m3/m ({m.group(2)} de {m.group(3)})")
print(f"V0.4 {'passou' if ok else 'falhou'}: balanco de agua em 24 h: " + "; ".join(txt))
EOF
)" ;;
    v05)
        B="$(base 0.5)"
        "$EXE" primeira=9 n=1 nystrom=$B elementos=tri ny=80 newton=1e-6,40,4 saida=v05 solver=skyline > v05.log 2>&1
        anota "$("$PY" -I - <<'EOF'
import csv
r = list(csv.DictReader(open('v05_9_9.csv')))[-1]
tp, tq, fp = float(r['t_parou']), float(r['t_passagem']), float(r['FS_passagem'])
ok = r['situacao'] == 'transiente' and 0. <= tq <= tp and fp < 1.
print(f"V0.5 {'passou' if ok else 'falhou'}: parada do transiente (realizacao 9, tri ny 80, Newton 1e-6,40,4): "
      f"t_parou {tp:.5f} h, t_passagem {tq:.4f} h (em [0, t_parou]), FS_passagem {fp:.4f} (< 1)")
EOF
)" ;;
    v12)
        "$EXE" det=1 n=0 I=0 tempos=0 saida=v12a solver=skyline > v12a.log 2>&1
        "$EXE" det=1 n=0 I=0 nivel0=2 tempos=0 saida=v12b solver=skyline > v12b.log 2>&1
        anota "$("$PY" -I - <<'EOF'
import csv, math
gw, gd, c, phi, beta = 9.81, 16., 5., math.radians(32), math.radians(40)
ts, tr, a, n = 0.358, 3.58e-4, 0.149, 1.289
Se = (1 + (a * 20) ** n) ** (-(1 - 1 / n))
gt, gs, sc = gd + (tr + (ts - tr) * Se) * gw, gd + ts * gw, math.sin(beta) * math.cos(beta)
casos = (('v12a', (c + Se * 20 * math.tan(phi)) / (gt * sc), math.tan(phi) / math.tan(beta), 1.4979),
         ('v12b', c / (gs * sc), (gs * math.cos(beta) ** 2 - gw) * math.tan(phi) / (gs * sc), 0.3669))
ok, txt = True, []
for nome, A, B, alvo in casos:
    e, f2 = 0., None
    for r in csv.DictReader(open(f'{nome}_perfil_-1.csv')):
        z = float(r['z'])
        if z > 0:
            e = max(e, abs(float(r['FS']) - (A / z + B)))
        if abs(z - 2) < 1e-12:
            f2 = float(r['FS'])
    ok = ok and e < 1e-10 and abs(f2 - alvo) <= 1e-3
    txt.append(f"{A:.4f}/z + {B:.4f}: FS(2) = {f2:.4f}, max |FS - formula| = {e:.1e}")
print(f"V1.2 {'passou' if ok else 'falhou'}: FS fechado, succao 20 kPa e coluna cheia: " + "; ".join(txt))
EOF
)" ;;
    v13)
        [ -s ../v04/v04_10.8_perfil_-1.csv ] || "$EXE" det=1 n=0 saida=../v04/v04_10.8 solver=skyline > /dev/null 2>&1
        "$PY" -I "$REF" fs ../v04/v04_10.8_perfil_-1.csv > v13_det.log 2>&1
        "$PY" -I "$REF" fs ../v12/v12b_perfil_-1.csv > v13_cheia.log 2>&1
        anota "$("$PY" -I - <<'EOF'
import re
v = [float(x) for f in ('v13_det.log', 'v13_cheia.log') for x in re.findall(r'FS do arquivo\| = (\S+)', open(f).read())]
print(f"V1.3 {'passou' if v and max(v) <= 1e-10 else 'falhou'}: o FS do perfil recalculado pelo coluna_ref.py "
      f"(deterministico do Cho em 6/12/18/24 h e coluna cheia): max |dFS| = {max(v):.1e} ({len(v)} instantes)")
EOF
)" ;;
    v21)
        # o C++: 48 deterministicos; a referencia: o coluna_ref.py det N=400 dt=0.025 nos pontos de cada malha
        : > cmds.txt
        for k in 40 80 160 200 400 800; do # pontos z = j H/k (tri: k = 2 ny; quad: k = ny)
            "$PY" -I -c "print('\n'.join(repr(2. * j / $k) for j in range($k + 1)))" > nos_$k.txt
            [ -s ref_$k\_0_-1.csv ] || echo "\"$PY\" -I \"$REF\" det N=400 dt=0.025 tempos=0,6,12,18,24 nos=nos_$k.txt saida=ref_$k > ref_$k.log 2>&1" >> cmds.txt
        done
        # as rodadas longas primeiro
        for ny in 400 200 80 40; do for dt in 0.025 0.05 0.1; do for e in tri quad; do for b in 4 0; do
            s="${e}_b${b}_ny${ny}_dt${dt}"
            [ -s "${s}_0_-1.csv" ] || echo "\"$EXE\" det=1 n=0 ny=$ny dt=$dt elementos=$e newton=1e-6,40,$b tempos=0,6,12,18,24 saida=$s solver=skyline > $s.log 2>&1" >> cmds.txt
        done; done; done; done
        paralelo "$P" cmds.txt
        anota "$("$PY" -I - <<'EOF'
import csv, math, itertools
H = 2.
def le(arq):
    d = {}
    for r in csv.DictReader(open(arq)):
        d[float(r['t_h'])] = r
    return d
T = (6., 12., 18., 24.)
NYS, DTS = (40, 80, 200, 400), (0.1, 0.05, 0.025)
res, ref = {}, {}
for e, b, ny, dt in itertools.product(('tri', 'quad'), (4, 0), NYS, DTS):
    try:
        res[e, b, ny, dt] = le(f"{e}_b{b}_ny{ny}_dt{dt}_0_-1.csv")
    except FileNotFoundError:
        pass
for k in (40, 80, 160, 200, 400, 800):
    ref[k] = le(f"ref_{k}_0_-1.csv")
f = lambda r, t, c: float(r[t][c])
out = open('v21_tabela.txt', 'w')
def p(s=''):
    print(s, file=out)
p("elementos busca ny dt | FS 6 12 18 24 h | frente 6 12 18 h (m) | t_passagem (h) | passos newton | s | "
  "dFS contra o ref nos mesmos pontos (6 12 18 24) | dfrente/dy | dt_pass")
for (e, b, ny, dt), r in sorted(res.items()):
    k = 2 * ny if e == 'tri' else ny
    R = ref[k]
    parou = any(r[t]['situacao'] != 'ok' for t in r)
    dfs = [f(r, t, 'FS_min') - f(R, t, 'FS_min') for t in T]
    dfr = [(f(r, t, 'frente_m') - f(R, t, 'frente_m')) / (H / ny) for t in T[:3]]
    tp, tpr = f(r, 24., 't_passagem'), f(R, 24., 't_passagem')
    p(f"{e:4s} {b} {ny:3d} {dt:5.3f} | " + " ".join(f"{f(r, t, 'FS_min'):.4f}" for t in T) + " | " +
      " ".join(f"{f(r, t, 'frente_m'):.4f}" for t in T[:3]) + f" | {tp:.4f} | {r[24.]['passos']} {r[24.]['newton']} | "
      f"{f(r, 24., 's_transiente'):.1f} | " + " ".join(f"{x:+.4f}" for x in dfs) + " | " +
      " ".join(f"{x:+.1f}" for x in dfr) + f" | {tp - tpr:+.4f}" + (" | PAROU" if parou else ""))
p()
p("referencia (coluna_ref.py det N=400 dt=0.025) nos pontos z = j H/k:")
for k, R in ref.items():
    p(f"  k = {k}: FS " + " ".join(f"{f(R, t, 'FS_min'):.4f}" for t in T) + "; frentes " +
      " ".join(f"{f(R, t, 'frente_m'):.4f}" for t in T[:3]) + f"; t_passagem {f(R, 24., 't_passagem'):.4f}")
# ordens observadas: em ny (dt = 0,025; ny 80, 200, 400, razoes 2,5 e 2) e em dt (ny = 400; razao 2)
def ordem_h(q1, q2, q3, h1, h2, h3):
    d1, d2 = q1 - q2, q2 - q3
    if d1 * d2 <= 0 or d2 == 0:
        return float('nan')
    g = lambda s: (h1 ** s - h2 ** s) / (h2 ** s - h3 ** s) - d1 / d2
    lo, hi = 0.05, 8.
    if g(lo) * g(hi) > 0:
        return float('nan')
    for _ in range(100):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if g(lo) * g(m) > 0 else (lo, m)
    return 0.5 * (lo + hi)
p()
p("ordens observadas (nan: nao monotono): em ny com dt = 0,025 (ny 80/200/400) e em dt com ny = 400 (0,1/0,05/0,025)")
for e, b in itertools.product(('tri', 'quad'), (4, 0)):
    lin = []
    for nome, q in [(f"FS {t:g} h", lambda r, t=t: f(r, t, 'FS_min')) for t in T] + [("t_passagem", lambda r: f(r, 24., 't_passagem'))]:
        try:
            oh = ordem_h(*(q(res[e, b, ny, 0.025]) for ny in (80, 200, 400)), *(H / ny for ny in (80, 200, 400)))
            ot = ordem_h(*(q(res[e, b, 400, dt]) for dt in DTS), *DTS)
        except KeyError:
            continue
        lin.append(f"{nome}: {oh:.2f} / {ot:.2f}")
    p(f"  {e} busca {b}: " + "; ".join(lin))
out.close()
# o criterio, na malha mais fina (ny = 400, dt = 0,025), busca 4, tri e quad
ok, txt = True, []
for e in ('tri', 'quad'):
    r = res.get((e, 4, 400, 0.025))
    if r is None:
        ok = False
        txt.append(f"{e}: falta a rodada mais fina")
        continue
    R = ref[800 if e == 'tri' else 400]
    dfs = [abs(f(r, t, 'FS_min') - f(R, t, 'FS_min')) for t in T]
    dfr = max(abs(f(r, t, 'frente_m') - f(R, t, 'frente_m')) for t in T[:3]) / (H / 400)
    dtp = abs(f(r, 24., 't_passagem') - f(R, 24., 't_passagem'))
    ok_e = dfs[0] <= 0.003 and dfs[3] <= 0.003 and dfs[1] <= 0.01 and dfs[2] <= 0.01 and dfr <= 1. + 1e-9 and dtp <= 0.05
    ok = ok and ok_e
    txt.append(f"{e}: |dFS| " + "/".join(f"{x:.4f}" for x in dfs) + f", frente {dfr:.1f} dy, t_passagem {dtp:.4f} h")
print(f"V2.1 {'passou' if ok else 'falhou'}: ny 400, dt 0,025 contra o coluna_ref.py N=400 dt=0.025 nos mesmos pontos: " +
      "; ".join(txt) + " (tabela e ordens em v21/v21_tabela.txt)")
EOF
)" ;;
    v31)
        r="../v21/${EL}_b4_ny${NY}_dt0.1_0_-1.csv" # a rodada do v21; sem ela, roda aqui
        if [ ! -s "$r" ]; then
            r="${EL}_ny${NY}_dt0.1_0_-1.csv"
            "$EXE" det=1 n=0 ny=$NY elementos=$EL dt=0.1 tempos=0,6,12,18,24 saida="${r%_0_-1.csv}" solver=skyline > "${r%_0_-1.csv}.log" 2>&1
        fi
        anota "$("$PY" -I - "$r" "$NY" "$EL" <<'EOF'
import csv, sys
d = {float(r['t_h']): r for r in csv.DictReader(open(sys.argv[1]))}
fs = {t: float(d[t]['FS_min']) for t in d}
fr = {t: float(d[t]['frente_m']) for t in d}
alvo = {0.: 1.50, 6.: 1.47, 18.: 1.06, 24.: 0.37}
ok = all(abs(fs[t] - v) <= 0.02 for t, v in alvo.items()) and 1.22 <= fs[12.] <= 1.27
faixas = {6.: (0.60, 0.65), 12.: (1.15, 1.25), 18.: (1.75, 1.85)}
ok = ok and all(a - 1e-9 <= fr[t] <= b + 1e-9 for t, (a, b) in faixas.items())
print(f"V3.1 {'passou' if ok else 'falhou'}: Fig. 3 com ny = {sys.argv[2]}, dt = 0,1, {sys.argv[3]}: FS " +
      "/".join(f"{fs[t]:.4f}" for t in (0., 6., 12., 18., 24.)) + " em 0/6/12/18/24 h; frentes " +
      "/".join(f"{fr[t]:.4f}" for t in (6., 12., 18.)) + f" m; t_passagem {float(d[24.]['t_passagem']):.4f} h")
EOF
)" ;;
    v23)
        B="$(base 0.5)"
        : > cmds.txt
        for e in tri quad; do for s in pardiso skyline; do for ny in 40 80 200; do
            c="${e}_${s}_ny$ny"
            [ -s "$c.log" ] && grep -q "tempo medio" "$c.log" || echo "\"$EXE\" n=20 nystrom=$B elementos=$e solver=$s ny=$ny saida=$c > $c.log 2>&1" >> cmds.txt
        done; done; done
        paralelo "$PC" cmds.txt
        anota "$("$PY" -I - "$PC" <<'EOF'
import re, sys, itertools
lin = []
for e, s, ny in itertools.product(('tri', 'quad'), ('pardiso', 'skyline'), (40, 80, 200)):
    log = open(f'{e}_{s}_ny{ny}.log').read()
    t = re.findall(r'tempo medio por realizacao: (\S+) s', log)
    # as realizacoes cujo transiente parou (o Newton nao convergiu num sub-passo): numero e lista (real@t_parou)
    par = re.findall(r'realizacao (\d+):.*o transiente parou em (\S+) h', log)
    lp = f", {len(par)} paradas" + (" (" + " ".join(f"{r}@{tp}" for r, tp in par) + ")" if par else "")
    lin.append((f"{e}/{s}/{ny}: {float(t[0]):.1f}" if t else f"{e}/{s}/{ny}: ?") + lp)
print(f"V2.3 (custo, {sys.argv[1]} processos juntos): s por realizacao e paradas (n = 20, CoV 0,5, lv 0,5): "
      + "; ".join(lin))
EOF
)" ;;
    *) echo "teste desconhecido: $T" ;;
    esac
done
