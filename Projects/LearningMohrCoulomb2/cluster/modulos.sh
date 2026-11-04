# Ajuste para o seu cluster (veja `module avail`). Lido pela compilacao e pelos jobs.
# fermi (CESUP/UFRGS): gcc do sistema + MKL do oneAPI 2025.1 (o modulo define MKLROOT, que o FindMKL usa)
module load mkl/2025.1
# export USAR_MKL=OFF      # descomente se o cluster nao tiver MKL
# export PARTICAO=normal   # particao do SLURM; vazio usa a padrao
