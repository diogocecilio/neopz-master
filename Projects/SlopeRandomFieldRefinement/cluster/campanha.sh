# campanha.sh: as campanhas desta pasta, escolhidas pela variavel CAMPANHA (lido por submeter.sh, job.sbatch, status.sh,
# baixar.sh, compilar.sh e teste.sbatch):
#   chuva (padrao)  Monte Carlo da chuva no talude de Rahardjo (SlopeRainMonteCarlo, casos.txt)
#   ceron           rebaixamento rapido do Ceron com o transiente, Monte Carlo direto e amostragem por importancia
#                   (SlopeDrawdownMonteCarlo, casos_ceron.txt)
#   chuva_ai        chuva com amostragem por importancia: piloto direto e AI, com o FS de 6 em 6 h ate 48 h
#                   (SlopeRainMonteCarlo, casos_chuva_ai.txt; blocos de 1 realizacao)
#   solos           os tres solos do Rahardjo com I = 1 k_s (a Fig. 5), Monte Carlo direto com o FS de 6 em 6 h ate
#                   48 h (SlopeRainMonteCarlo, casos_solos.txt; blocos de 5 realizacoes)
#   etapa1          etapa 1 do plano de trabalho (Projects/PlanoDeTrabalho): o SRM novo contra o antigo e a malha de teste
#                   do Cho2010 (Cho2010 e, com prog=, SlopeRainMonteCarlo; casos_etapa1.txt; blocos de 10)
# Define NOME (jobs no SLURM), PROG (executavel em bin/), CASOS_TXT, PASTA (em $DATA), IREF_PADRAO, BLOCO_PADRAO e DESTINO
# (pasta no PC, em ~/projects/runs)
CAMPANHA="${CAMPANHA:-chuva}"
case "$CAMPANHA" in
chuva) NOME=chuvaMC PROG=SlopeRainMonteCarlo CASOS_TXT=casos.txt PASTA=chuva_mc IREF_PADRAO=31.05 BLOCO_PADRAO=20
    DESTINO=fermi_chuva ;;
ceron) NOME=ceronMC PROG=SlopeDrawdownMonteCarlo CASOS_TXT=casos_ceron.txt PASTA=ceron_ai IREF_PADRAO=0 BLOCO_PADRAO=20
    DESTINO=fermi_ceron_ai ;;
chuva_ai) NOME=chuvaAI PROG=SlopeRainMonteCarlo CASOS_TXT=casos_chuva_ai.txt PASTA=chuva_ai IREF_PADRAO=31.05
    BLOCO_PADRAO=1 DESTINO=fermi_chuva_ai ;;
solos) NOME=solosMC PROG=SlopeRainMonteCarlo CASOS_TXT=casos_solos.txt PASTA=chuva_solos IREF_PADRAO=0
    BLOCO_PADRAO=5 DESTINO=fermi_chuva_solos ;;
etapa1) NOME=etapa1 PROG=Cho2010 CASOS_TXT=casos_etapa1.txt PASTA=etapa1 IREF_PADRAO=0 BLOCO_PADRAO=10
    DESTINO=fermi_etapa1 ;;
*) echo "CAMPANHA=$CAMPANHA desconhecida: chuva, ceron, chuva_ai, solos ou etapa1" >&2; exit 2 ;;
esac
export CAMPANHA
# todos os nomes de job (o compilar.sh espera todos acabarem: a libpz.so e a mesma)
NOMES_TODOS=chuvaMC,ceronMC,chuvaAI,solosMC,etapa1
