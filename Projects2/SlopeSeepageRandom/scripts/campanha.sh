#!/bin/bash
# Campanha completa do artigo: Monte Carlo de referência e de sensibilidade (Tabelas 3-6), exemplos de Cho (2010),
# Cam-Clay e rebaixamento acoplado. Gera <dir>/jobs.txt e o executa com P processos (fila.sh); no fim escreve
# <dir>/tabelas.md (as colunas determinísticas vêm de deterministico.sh, rodado em <dir>/det).
# Uso: campanha.sh <executável> <diretório> [N amostras por caso = 300] [P processos = 3] [N referência = 1000]
set -u
E=$(readlink -f "$1"); OUT=$2; N=${3:-300}; P=${4:-3}; NREF=${5:-1000}
S=$(dirname "$(readlink -f "$0")")
mkdir -p "$OUT"; cd "$OUT"
B="caso=percolacao modelo=mc h=1 adapt=2 hkl=1 covc=0.3 covphi=0.1 covk=0.6 Lx=20 Ly=2"
mc() { local nome=$1; shift; echo "$S/montecarlo.sh $E mc/$nome $nome $N 1 $*"; }
rb() {
    local nome=$1; shift
    echo "mkdir -p rebaix && cd rebaix && $E rebaixamento caso=percolacao h=1 adapt=2 gamma=1 $* saida=$nome.csv > $nome.log 2>&1"
}
{
    echo "$S/montecarlo.sh $E mc_ref ref $NREF $P $B"
    # Tabela 3: coeficientes de variação
    mc covk0 $B covk=0;     mc covk100 $B covk=1.0
    mc covc10 $B covc=0.1;  mc covc50 $B covc=0.5;  mc covc70 $B covc=0.7
    mc covphi5 $B covphi=0.05; mc covphi20 $B covphi=0.2
    # Tabela 4: escala s das distâncias de autocorrelação (Lx = 20 s, Ly = 2 s)
    mc s2 $B Lx=40 Ly=4; mc s5 $B Lx=100 Ly=10; mc s20 $B Lx=400 Ly=40
    # Tabela 5: anisotropia; Tabela 6: rebaixamento h_w/H
    mc alfa2 $B alpha=2; mc alfa3 $B alpha=3; mc alfa5 $B alpha=5
    mc hw0.5 $B hw=2.5;  mc hw0.7 $B hw=3.5;  mc hw0.9 $B hw=4.5
    # seção 5.3: Cho (2010), sem percolação
    mc cho_coesivo caso=cho_coesivo modelo=mc h=2 adapt=2 hkl=1 covc=0.3 covphi=0 covk=0 Lx=20 Ly=2
    mc cho_cphi caso=cho_cphi modelo=mc h=2 adapt=2 hkl=1 covc=0.3 covphi=0.2 covk=0 Lx=20 Ly=2
    # Cam-Clay modificado na configuração de referência
    mc mcc_ref caso=percolacao modelo=mcc h=1 adapt=2 hkl=1 covc=0.3 covphi=0.1 covk=0.6 Lx=20 Ly=2
    # rebaixamento acoplado (u-p): Mohr-Coulomb e Cam-Clay com OCR = 1 e 2
    for td in 0.001 0.01 0.1 1 10; do rb mc_Td$td modelo=mc Td=$td; done
    for ocr in 1 2; do for td in 0.1 1 10; do rb mcc_ocr${ocr}_Td$td modelo=mcc OCR=$ocr Td=$td; done; done
} > jobs.txt
"$S/fila.sh" jobs.txt "$P"
python3 "$S/tabelas.py" . > tabelas.md
echo "tabelas em $OUT/tabelas.md"
