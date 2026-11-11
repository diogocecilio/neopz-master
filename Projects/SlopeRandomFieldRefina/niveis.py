#!/usr/bin/env python3
"""niveis.py: o estudo dos niveis de refinamento em volta da banda de ruptura de cada realizacao (refina_banda.py), para
escolher o tamanho h dos elementos em volta da banda, o numero de niveis (um passo ou dois) e a malha inicial (a da
etapa 1), que pode ser mais grossa que a foco375 da comparacao (compara.py).

Modelo e casos: o exemplo saturado do Ceron exatamente como em compara.py (MODELO: GIM com arc-length em 0,01 h), com
SOLVER=skyline e 1 thread, no deterministico (det, CAMPOS=0) e nas realizacoes 0 a 7 (REAL=n), com a SEMENTE padrao.

Malhas iniciais (etapa 1), em ../LearningMohrCoulomb2/cluster/malhas, todas com o mesmo contorno, os mesmos 6 cantos e
os mesmos grupos fisicos: o retangulo do Karhunen-Loeve (0 a 25 m por 0 a 10 m) e o campo aleatorio de cada realizacao
sao os mesmos em todas as malhas, iniciais e refinadas.
  foco375  ceron_h5_a45_foco375 (720 triangulos): 0,375 m na regiao das bandas, 0,9 m na superficie e 3 m no fundo;
  g09      ceron_h5_a45_g09 (235 triangulos): a foco375 sem o foco, hSup = 0,9 m na superficie e hMax = 3 m;
  g15      ceron_h5_a45_g15 (126 triangulos): a mesma com hSup = 1,5 m.
Niveis: as configuracoes <inicial>/<tipo>, com folga 0,375 m e limiar 0,1 (os padroes de srf_refinado.sh). Toda malha
refinada e a de refina_banda.py com o .geo da malha INICIAL da cadeia como base (os campos de tamanho dela, mais o da
banda) e a banda da rodada pai:
  <inicial>/inicial              a rodada na malha inicial (a etapa 1);
  <inicial>/um_passo_h<h>        h = 0,25, 0,1875 e 0,125 m, da banda da rodada <inicial>/inicial;
  g09|g15/dois_niveis_h0.125     h = 0,125 m, da banda da rodada <inicial>/um_passo_h0.25;
  foco375/referencia_h0.0625     h = 0,0625 m, da banda da rodada foco375/um_passo_h0.125: a referencia do erro;
  foco375/extra_h0.046875        so o det, por ultimo na fila (--sem-extra tira): h = 0,046875 m, da banda da
                                 referencia, para estimar o erro de discretizacao da propria referencia.
O custo de um metodo num caso e o da cadeia: a etapa 1 + o refina_banda.py (com o gmsh) + a rodada refinada [+ o
segundo refina_banda.py + a segunda rodada refinada].

niveis.py rodar <saida> [--casos det,0-7] [--paralelo 2] [--sem-extra] [--configs <c1,c2,...>] [--listar]
    Roda (ou retoma) a campanha, cada tarefa em <saida>/<inicial>/<tipo>/<caso>/ (caso: det ou r0000n).
    - Uma tarefa e uma rodada do programa (na malha inicial ou numa refinada) e, nas refinadas, antes dela, o
      refina_banda.py na pasta da tarefa. Uma tarefa refinada espera a pai terminar. Se o refina_banda.py sai com
      codigo 3 (a pai nao tem banda), a tarefa fica com sem_banda.txt e as que dependem dela sao puladas.
    - O tempo, medido como em compara.py: uma fila so, PARALELO processos ao mesmo tempo (padrao 2), e as tarefas
      prontas tomadas na ordem (caso, configuracao), a da lista acima. Assim toda configuracao roda ao lado do mesmo
      numero de processos da fila (menos no fim da fila). tempo.txt de cada tarefa: parede_s (relogio, de logo antes de
      criar o processo do programa ate o wait4 que o recolhe), cpu_s (usuario + sistema do processo e dos filhos que ele
      esperou, do wait4), inicio e fim (s desde 1970), paralelo, carga_inicio e carga_fim (a carga media de 1 min,
      /proc/loadavg) e codigo (a saida do programa); nas refinadas, tambem o refina_banda.py com o gmsh, medido a parte
      do mesmo jeito: seg_malha (relogio), cpu_malha_s, inicio_malha e carga_malha.
    - Executavel congelado: na primeira vez, o SlopeRandomField (SRF_BIN, padrao ../SlopeRandomField/bin), o seu
      versao.txt e as bibliotecas de fora do sistema que ele usa (a libpz.so, pelo ldd) sao copiados para
      <saida>/bin_congelado (com manifesto.txt: sha256 e origem), e toda rodada usa essa copia, com LD_LIBRARY_PATH; o
      ldd confere, a cada inicio, que as bibliotecas sao as da copia. As malhas iniciais tambem sao copiadas, para
      <saida>/malhas_iniciais, e sao a base de toda malha refinada.
    - Arquivos de cada tarefa: os do programa (resultado.txt, resumo.csv, fs.csv, fs_tentativas.csv, passos.csv,
      realizacao.txt, triangulos.csv e mecanismo.csv so da banda: MEC_CSV=1 e MEC_LIMIAR=0,1, como na etapa 1 de
      srf_refinado.sh, sem pw.csv, campos.csv e VTK), run.log e run.err, malha.geo e malha.msh (a malha da rodada),
      malha.txt (triangulos e, nas refinadas, a linha de refina_banda.py: celulas_banda, fator_grossa (o da pai), t_h,
      h, folga e limiar, mais a pai e a base) com malha.log e malha.err, e tempo.txt. O checkpoint.bin e apagado
      depois do sucesso.
    - Retomavel: uma tarefa com resultado.txt e tempo.txt (e com a pai feita) e pulada; uma pasta incompleta e apagada e
      a tarefa roda do zero (o tempo e o de uma execucao inteira), e as pastas das tarefas que dependem dela tambem.
      Ctrl-C ou SIGTERM no niveis.py param a fila e matam os filhos; rodar de novo continua. Uma falha (sem
      resultado.txt) fica em falhou.txt (a etapa, o codigo e o fim do .err), as tarefas independentes continuam e as
      que dependem dela sao puladas; numa retomada, ela roda de novo. Um so niveis.py por pasta (trava em
      <saida>/.trava).
    - Progresso: uma linha por tarefa terminada em <saida>/progresso.log (data, tarefa, fator, parede_s, cpu_s,
      terminadas/total); a carga da maquina a cada 30 s em <saida>/carga.log; cada inicio (com a CPU da maquina) em
      <saida>/sessoes.log; <saida>/CONCLUIDO quando a fila acaba, tambem com falhas (listadas nele).
    - --configs: so essas configuracoes (e as suas ancestrais), por exemplo --configs g15/um_passo_h0.25 para um teste;
      --listar: so mostra as tarefas (ordem, pai e estado) e sai.
niveis.py tabela <saida> [--csv <arquivo>] [--comparacao <pasta>]
    Uma linha por rodada (configuracao e caso): caso, inicial, tipo, h, estado, fator, e_ref (a diferenca relativa
    para a referencia do caso), triangulos, neq, celulas_banda, seg_malha, parede_s, cpu_s, cadeia_parede_s e
    cadeia_cpu_s (o custo do metodo: a soma das etapas da cadeia), carga_media (a media da carga no inicio e no fim),
    simultaneos (processos da fila rodando junto, media no tempo, contando ele mesmo) e, nas refinadas, a contencao da
    banda: dist_max, a maior distancia (m) dos pontos da banda da rodada (sqrt(J2) >= 0,1 do maximo no instante do
    fator, de mecanismo.csv) aos pontos da banda com que a malha dela foi feita (os Point que refina_banda.py pos no
    malha.geo), e fora, a fracao dos pontos da banda a mais de 'folga' deles. Mais as malhas fixas da comparacao
    (--comparacao, padrao <saida>/../refina/comparacao: foco375, foco25, foco1875 e foco125; fator, tempos e
    triangulos). --csv grava as linhas. Imprime tabelas em Markdown (decimais com virgula): o fator e o erro de cada
    metodo por caso, o resumo por metodo (|e| medio e maximo, e medio com sinal, tempos medios da cadeia e a razao para
    os da foco125 da comparacao), as etapas e as conferencias com a comparacao.
Ambiente: SRF_BIN (o executavel a congelar, so na primeira vez), SOLVER (padrao skyline), PYTHON (o python de
refina_banda.py, com o modulo gmsh; padrao o que roda este)."""

import argparse
import csv
import fcntl
import hashlib
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
import traceback

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import compara  # noqa: E402  o MODELO e as funcoes de leitura e de formato da comparacao
import refina_banda  # noqa: E402  a contagem de triangulos de um .msh e o codigo de saida "sem banda"

MALHAS = compara.MALHAS
REFINA = os.path.join(AQUI, "refina_banda.py")
SRF_BIN = os.environ.get("SRF_BIN",
                         os.path.normpath(os.path.join(AQUI, "..", "SlopeRandomField", "bin", "SlopeRandomField")))
PYTHON = os.environ.get("PYTHON", sys.executable)
SOLVER = os.environ.get("SOLVER", "skyline")
FOLGA, LIMIAR = 0.375, 0.1  # os padroes de srf_refinado.sh
INICIAIS = ("foco375", "g09", "g15")
H_UM = (0.25, 0.1875, 0.125)
H_DOIS, H_REF, H_EXTRA = 0.125, 0.0625, 0.046875
REF = "foco375/referencia_h%g" % H_REF
# as malhas fixas da comparacao (compara.py) e o hFoco de cada uma
FIXAS = (("foco375", 0.375), ("foco25", 0.25), ("foco1875", 0.1875), ("foco125", 0.125))
# as saidas de toda rodada: as da etapa 1 de srf_refinado.sh (mecanismo.csv so da banda; sem pw.csv, campos.csv e VTK)
SAIDAS = dict(MEC_CSV="1", MEC_LIMIAR="%g" % LIMIAR, PW_CSV="0", CAMPOS_CSV="0", CAMPOS_VTK="0", VTK_H="")
# o ambiente das rodadas: so estas variaveis vem de quem chama (nenhuma variavel do modelo vaza para as rodadas)
HERDADAS = ("PATH", "HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TZ", "TMPDIR")
# um filho que morre por um destes sinais foi parado (Ctrl-C no terminal, kill): a fila para, sem contar falha
SINAIS_PARADA = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
TERMINAIS = ("feito", "sem_banda", "falhou", "pulado")


# ---------------------------------------------------------------------------------------------------------------------
# as tarefas

def nome_caso(caso):
    return "det" if caso == "det" else "r%05d" % int(caso)


def nome_config(inicial, tipo, h=None):
    return "%s/%s" % (inicial, tipo if h is None else "%s_h%g" % (tipo, h))


def configuracoes(extra=True):
    """[(configuracao, inicial, tipo, h, configuracao pai)] na ordem da fila (a pai sempre antes)"""
    cs = [(nome_config(i, "inicial"), i, "inicial", None, None) for i in INICIAIS]
    for i in INICIAIS:
        cs += [(nome_config(i, "um_passo", h), i, "um_passo", h, nome_config(i, "inicial")) for h in H_UM]
    cs += [(nome_config(i, "dois_niveis", H_DOIS), i, "dois_niveis", H_DOIS, nome_config(i, "um_passo", 0.25))
           for i in ("g09", "g15")]
    cs.append((REF, "foco375", "referencia", H_REF, nome_config("foco375", "um_passo", 0.125)))
    if extra:  # so o det, por ultimo na fila
        cs.append((nome_config("foco375", "extra", H_EXTRA), "foco375", "extra", H_EXTRA, REF))
    return cs


class Tarefa:
    def __init__(self, base, cfg, caso, ordem):
        self.config, self.inicial, self.tipo, self.h, self.config_pai = cfg
        self.caso, self.ordem = caso, ordem
        self.pasta = os.path.join(base, self.config, nome_caso(caso))
        self.pai = None
        self.estado = "pendente"
        self.motivo = ""

    def arq(self, nome):
        return os.path.join(self.pasta, nome)


def tarefas(base, lista_casos, extra=True, so=None):
    """as tarefas, na ordem da fila: por caso e, em cada caso, por configuracao; a extra (so o det) depois de todas.
    so: so essas configuracoes e as suas ancestrais"""
    cfgs = configuracoes(extra)
    if so:
        pais = {c[0]: c[4] for c in cfgs}
        erradas = [c for c in so if c not in pais]
        if erradas:
            sys.exit("niveis: configuracao desconhecida: %s (as configuracoes: %s)" % (
                ", ".join(erradas), ", ".join(pais)))
        quero = set()
        for c in so:
            while c:
                quero.add(c)
                c = pais[c]
        cfgs = [c for c in cfgs if c[0] in quero]
    ts, por = [], {}
    for k, caso in enumerate(lista_casos):
        for j, cfg in enumerate(cfgs):
            if cfg[2] == "extra" and caso != "det":
                continue
            t = Tarefa(base, cfg, caso, (len(lista_casos), j) if cfg[2] == "extra" else (k, j))
            ts.append(t)
            por[t.config, caso] = t
    for t in ts:
        if t.config_pai:
            t.pai = por[t.config_pai, t.caso]
    ts.sort(key=lambda t: t.ordem)
    return ts


def completa(pasta):
    return os.path.isfile(os.path.join(pasta, "resultado.txt")) and os.path.isfile(os.path.join(pasta, "tempo.txt"))


def estado_inicial(ts):
    """o estado de cada tarefa pelo que esta no disco; uma tarefa so conta como feita se a pai tambem conta"""
    for t in ts:  # a pai vem antes da filha
        if t.pai is not None and t.pai.estado != "feito":
            t.estado = "pendente"
        elif completa(t.pasta):
            t.estado = "feito"
        elif t.pai is not None and os.path.isfile(t.arq("sem_banda.txt")):
            t.estado = "sem_banda"
            t.motivo = "a pai %s nao tem banda" % t.config_pai
        else:
            t.estado = "pendente"


def descende(t, a):
    """t depende (direta ou indiretamente) de a"""
    p = t.pai
    while p is not None:
        if p is a:
            return True
        p = p.pai
    return False


# ---------------------------------------------------------------------------------------------------------------------
# arquivos, executavel congelado e malhas iniciais

def escreve(arq, txt):
    """grava com troca atomica"""
    with open(arq + ".tmp", "w") as f:
        f.write(txt)
    os.replace(arq + ".tmp", arq)


def fim_de(arq, n=2000):
    try:
        with open(arq, errors="replace") as f:
            return f.read()[-n:]
    except OSError:
        return ""


def sha256(arq):
    h = hashlib.sha256()
    with open(arq, "rb") as f:
        for bloco in iter(lambda: f.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


def bibliotecas(exe, env):
    """{nome: caminho} das bibliotecas que o ldd resolve para o executavel fora do sistema (/lib*, /usr)"""
    r = subprocess.run(["ldd", exe], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    libs = {}
    for linha in r.stdout.splitlines():
        p = linha.split()
        if len(p) >= 3 and p[1] == "=>":
            if p[2] == "not":
                sys.exit("niveis: o ldd de %s nao encontra %s" % (exe, p[0]))
            if not p[2].startswith(("/lib/", "/lib64/", "/usr/")):
                libs[p[0]] = p[2]
    return libs


def congela(base):
    """(executavel, versao, LD_LIBRARY_PATH) da copia congelada em <base>/bin_congelado (copiada na primeira vez)"""
    d = os.path.join(base, "bin_congelado")
    exe = os.path.join(d, "SlopeRandomField")
    if not os.path.isfile(exe):
        if not os.access(SRF_BIN, os.X_OK):
            sys.exit("niveis: executavel %s nao encontrado (SRF_BIN)" % SRF_BIN)
        libs = bibliotecas(SRF_BIN, {k: v for k, v in os.environ.items() if k != "LD_LIBRARY_PATH"})
        if not any(n.startswith("libpz") for n in libs):
            sys.exit("niveis: o ldd de %s nao mostra a libpz" % SRF_BIN)
        tmp = d + ".tmp"
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        origens = [SRF_BIN, os.path.join(os.path.dirname(SRF_BIN), "versao.txt")] + sorted(libs.values())
        man = []
        for o in origens:
            shutil.copy2(o, tmp)
            man.append("%s %s %s %s" % (sha256(o), os.path.basename(o), o,
                                        time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime(os.path.getmtime(o)))))
        escreve(os.path.join(tmp, "manifesto.txt"), "# sha256 arquivo origem data_da_origem (copiado em %s)\n%s\n" % (
            time.strftime("%Y-%m-%d %H:%M:%S"), "\n".join(man)))
        os.rename(tmp, d)
    # o ldd da copia, com LD_LIBRARY_PATH, tem de resolver cada biblioteca de fora do sistema para a copia
    ld = d + (":" + os.environ["LD_LIBRARY_PATH"] if os.environ.get("LD_LIBRARY_PATH") else "")
    libs = bibliotecas(exe, dict(os.environ, LD_LIBRARY_PATH=ld))
    erradas = ["%s => %s" % (n, c) for n, c in libs.items()
               if os.path.dirname(os.path.realpath(c)) != os.path.realpath(d)]
    if erradas or not any(n.startswith("libpz") for n in libs):
        sys.exit("niveis: a copia congelada %s nao usa as bibliotecas copiadas: %s" % (exe, "; ".join(erradas) or
                                                                                      "sem libpz"))
    with open(os.path.join(d, "versao.txt")) as f:
        return exe, f.read().strip(), ld


def copia_malhas(base):
    """as malhas iniciais em <base>/malhas_iniciais (copiadas na primeira vez): a base de toda malha refinada"""
    d = os.path.join(base, "malhas_iniciais")
    os.makedirs(d, exist_ok=True)
    for i in INICIAIS:
        for ext in (".geo", ".msh"):
            dst = os.path.join(d, "ceron_h5_a45_%s%s" % (i, ext))
            if not os.path.isfile(dst):
                shutil.copy2(os.path.join(MALHAS, "ceron_h5_a45_%s%s" % (i, ext)), dst + ".tmp")
                os.replace(dst + ".tmp", dst)


def mata_orfao(pasta):
    """um processo de uma vez anterior (de um niveis.py morto com SIGKILL) que ainda rode nesta pasta: SIGTERM e, depois
    de 10 s, SIGKILL (a trava garante que nenhum outro niveis.py roda nesta campanha)"""
    try:
        pid = int(compara.chaves(os.path.join(pasta, "pid.txt")).get("pid", ""))
        if os.path.realpath(os.readlink("/proc/%d/cwd" % pid)) != os.path.realpath(pasta):
            return
    except (ValueError, OSError):
        return
    for sinal in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.kill(pid, sinal)
        except ProcessLookupError:
            return
        limite = time.time() + 10
        while time.time() < limite:
            try:
                os.readlink("/proc/%d/cwd" % pid)
            except OSError:
                return
            time.sleep(0.2)


def limpa(pasta):
    if os.path.lexists(pasta):
        mata_orfao(pasta)
        shutil.rmtree(pasta)


def cpu_modelo():
    try:
        with open("/proc/cpuinfo") as f:
            for linha in f:
                if linha.startswith("model name"):
                    return linha.split(":", 1)[1].strip()
    except OSError:
        pass
    return "-"


# ---------------------------------------------------------------------------------------------------------------------
# rodar

class Campanha:
    def __init__(self, base, ts, paralelo, exe, versao, ld):
        self.base, self.ts, self.paralelo = base, ts, paralelo
        self.exe, self.versao, self.ld = exe, versao, ld
        self.sessao = int(time.time())
        self.cond = threading.Condition()  # o estado das tarefas
        self.parar = threading.Event()
        self.trava_filhos = threading.RLock()  # RLock: o tratador de sinal pode entrar com a thread principal dentro
        self.filhos = {}  # pid -> Popen dos processos em andamento (ate serem recolhidos)
        env = {k: os.environ[k] for k in HERDADAS if k in os.environ}
        env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", LC_ALL="C")
        self.env_base = env

    # ---- processos
    def processo(self, cmd, pasta, env, saida, erro):
        """roda cmd na pasta e espera: (codigo, inicio, fim, cpu_s); (None, ...) se a parada ja foi pedida"""
        with open(os.path.join(pasta, saida), "w") as out, open(os.path.join(pasta, erro), "w") as err:
            with self.trava_filhos:
                if self.parar.is_set():
                    return None, None, None, None
                t0 = time.time()
                p = subprocess.Popen(cmd, cwd=pasta, env=env, stdout=out, stderr=err, stdin=subprocess.DEVNULL)
                self.filhos[p.pid] = p
            escreve(os.path.join(pasta, "pid.txt"), "pid=%d\n" % p.pid)
            # espera sem recolher (WNOWAIT): ate sair de self.filhos o pid nao pode ser reusado, e um kill nele e seguro
            os.waitid(os.P_PID, p.pid, os.WEXITED | os.WNOWAIT)
            t1 = time.time()
            with self.trava_filhos:
                del self.filhos[p.pid]
                _, status, uso = os.wait4(p.pid, 0)
        os.remove(os.path.join(pasta, "pid.txt"))
        p.returncode = os.waitstatus_to_exitcode(status)
        return p.returncode, t0, t1, uso.ru_utime + uso.ru_stime

    def mata_filhos(self, sinal):
        with self.trava_filhos:
            for pid in list(self.filhos):
                try:
                    os.kill(pid, sinal)
                except ProcessLookupError:
                    pass

    def interrompido(self, rc):
        """a parada foi pedida (sinal no niveis.py) ou o filho morreu por um sinal de parada (Ctrl-C no terminal)"""
        if rc is not None and (-rc in SINAIS_PARADA or rc - 128 in SINAIS_PARADA):
            self.parar.set()
        return self.parar.is_set()

    # ---- uma tarefa
    def falhou(self, t, etapa, rc, arq_err):
        escreve(t.arq("falhou.txt"), "etapa=%s codigo=%s\n%s" % (etapa, rc, fim_de(t.arq(arq_err))))
        t.motivo = "%s, codigo %s" % (etapa, rc)
        return "falhou"

    def executa(self, t):
        for u in self.ts:  # as que dependem desta usam a banda dela: refeitas depois dela
            if descende(u, t):
                limpa(u.pasta)
        limpa(t.pasta)
        os.makedirs(t.pasta)
        base_geo = os.path.join(self.base, "malhas_iniciais", "ceron_h5_a45_%s.geo" % t.inicial)
        extra = ""
        if t.pai is None:
            for ext in (".geo", ".msh"):
                shutil.copy2(base_geo[:-4] + ext, t.arq("malha" + ext))
            try:
                tri = "%d" % refina_banda.triangulos(t.arq("malha.msh"))
            except (OSError, ValueError, IndexError):  # um .msh que nao se le: o programa vai falhar e dizer por que
                tri = "-"
            escreve(t.arq("malha.txt"), "malha=%s triangulos=%s\n" % (os.path.basename(base_geo[:-4]) + ".msh", tri))
        else:
            # o refina_banda.py (com o gmsh): a banda da pai, a base e a malha inicial da cadeia
            c0 = compara.carga()
            cmd = [PYTHON, REFINA, t.pai.pasta, base_geo, t.arq("malha"), "%g" % t.h, "%g" % FOLGA, "%g" % LIMIAR]
            rc, m0, m1, mcpu = self.processo(cmd, t.pasta, self.env_base, "malha.log", "malha.err")
            if rc is None or self.interrompido(rc):
                return "interrompido"
            if rc == refina_banda.SEM_BANDA:
                escreve(t.arq("sem_banda.txt"), "pai=%s codigo=%d\n%s" % (t.config_pai, rc, fim_de(t.arq("malha.err"))))
                t.motivo = "a pai %s nao tem banda" % t.config_pai
                return "sem_banda"
            if rc != 0 or not os.path.isfile(t.arq("malha.msh")):
                return self.falhou(t, "refina_banda.py", rc, "malha.err")
            with open(t.arq("malha.log")) as f:
                chaves = [p for p in f.read().split() if "=" in p and not p.startswith("malha=")]
            escreve(t.arq("malha.txt"), "malha=malha.msh %s pai=%s base=%s\n" % (
                " ".join(chaves), t.config_pai, os.path.basename(base_geo)))
            extra = " seg_malha=%.3f cpu_malha_s=%.3f inicio_malha=%.3f carga_malha=%s" % (m1 - m0, mcpu, m0, c0)
        # o programa, na malha da tarefa
        env = dict(self.env_base)
        env.update(compara.MODELO)
        env.update(SAIDAS)
        env.update(SOLVER=SOLVER, CAMPOS="0" if t.caso == "det" else "1", REAL="0" if t.caso == "det" else t.caso,
                   SRF_VERSAO=self.versao, MESH=t.arq("malha.msh"), LD_LIBRARY_PATH=self.ld)
        c0 = compara.carga()
        rc, t0, t1, cpu = self.processo([self.exe], t.pasta, env, "run.log", "run.err")
        if rc is None:
            return "interrompido"
        c1 = compara.carga()
        if not (os.path.isfile(t.arq("resultado.txt")) and 0 <= rc < 128):
            if self.interrompido(rc):
                return "interrompido"
            return self.falhou(t, "programa", rc, "run.err")
        for f in ("checkpoint.bin", "checkpoint.bin.tmp"):
            if os.path.exists(t.arq(f)):
                os.remove(t.arq(f))
        # tempo.txt por ultimo: com resultado.txt, marca a tarefa como feita
        escreve(t.arq("tempo.txt"), "parede_s=%.3f cpu_s=%.3f inicio=%.3f fim=%.3f paralelo=%d carga_inicio=%s "
                                    "carga_fim=%s codigo=%d%s sessao=%d\n" % (
                                        t1 - t0, cpu, t0, t1, self.paralelo, c0, c1, rc, extra, self.sessao))
        return "feito"

    # ---- a fila
    def registra(self, t):
        """uma linha no progresso.log e na saida (com self.cond)"""
        if t.estado == "feito":
            r, tm = compara.chaves(t.arq("resultado.txt")), compara.chaves(t.arq("tempo.txt"))
            txt = "fator=%s parede_s=%s cpu_s=%s" % (r.get("FS_min", "-"), tm.get("parede_s", "-"), tm.get("cpu_s", "-"))
            if "seg_malha" in tm:
                txt += " seg_malha=%s" % tm["seg_malha"]
        else:
            txt = {"sem_banda": "SEM BANDA", "falhou": "FALHOU", "pulado": "PULADA"}[t.estado] + ": " + t.motivo
        n = sum(u.estado in TERMINAIS for u in self.ts)
        linha = "%s %s %s: %s (%d/%d)" % (time.strftime("%Y-%m-%d %H:%M:%S"), t.config, t.caso, txt, n, len(self.ts))
        with open(os.path.join(self.base, "progresso.log"), "a") as f:
            f.write(linha + "\n")
        print(linha, flush=True)

    def proxima(self):
        """a primeira tarefa pronta na ordem da fila (com self.cond); as filhas de uma tarefa sem banda, que falhou ou
        foi pulada sao puladas"""
        for t in self.ts:
            if t.estado != "pendente":
                continue
            if t.pai is None or t.pai.estado == "feito":
                return t
            if t.pai.estado in ("sem_banda", "falhou", "pulado"):
                t.estado = "pulado"
                t.motivo = "a pai %s: %s" % (t.config_pai, t.pai.estado)
                self.registra(t)
        return None

    def trabalhador(self):
        while True:
            with self.cond:
                while True:
                    if self.parar.is_set():
                        return
                    t = self.proxima()
                    if t is not None:
                        t.estado = "rodando"
                        break
                    if not any(u.estado in ("pendente", "rodando") for u in self.ts):
                        self.cond.notify_all()
                        return
                    self.cond.wait(1.)
            try:
                estado = self.executa(t)
            except Exception:  # um erro do proprio niveis.py nesta tarefa: conta como falha, as outras continuam
                estado = "falhou"
                t.motivo = "erro no niveis.py"
                try:
                    escreve(t.arq("falhou.txt"), "etapa=niveis.py codigo=-\n" + traceback.format_exc())
                except OSError:
                    pass
            with self.cond:
                t.estado = "pendente" if estado == "interrompido" else estado
                if estado != "interrompido":
                    self.registra(t)
                self.cond.notify_all()

    def amostra_carga(self, fim):
        """a carga da maquina (/proc/loadavg) a cada 30 s, em carga.log"""
        while not fim.wait(30.):
            try:
                with open("/proc/loadavg") as f, open(os.path.join(self.base, "carga.log"), "a") as g:
                    g.write("%d %s %s" % (time.time(), time.strftime("%H:%M:%S"), f.read()))
            except OSError:
                pass

    def roda(self):
        sinal = []

        def ao_sinal(signum, _frame):  # na thread principal; sem print (reentrancia) e sem esperar trava da fila
            sinal.append(signum)
            self.parar.set()
            self.mata_filhos(signal.SIGTERM)

        signal.signal(signal.SIGINT, ao_sinal)
        signal.signal(signal.SIGTERM, ao_sinal)
        fim = threading.Event()
        threading.Thread(target=self.amostra_carga, args=(fim,), daemon=True).start()
        ts = [threading.Thread(target=self.trabalhador, daemon=True) for _ in range(self.paralelo)]
        for th in ts:
            th.start()
        t_parada = None
        while any(th.is_alive() for th in ts):
            time.sleep(0.5)
            if self.parar.is_set():
                if t_parada is None:
                    t_parada = time.time()
                    print("niveis: parando%s; as tarefas em andamento ficam incompletas e rodam do zero na proxima vez"
                          % (" (sinal %d)" % sinal[0] if sinal else " (um filho foi parado por um sinal)"), flush=True)
                elif time.time() - t_parada > 15:
                    self.mata_filhos(signal.SIGKILL)
        fim.set()
        if self.parar.is_set():
            return 128 + (sinal[0] if sinal else signal.SIGINT)
        conta = {e: sum(t.estado == e for t in self.ts) for e in TERMINAIS}
        txt = "concluido=%s tarefas=%d feitas=%d sem_banda=%d falharam=%d puladas=%d\n" % (
            time.strftime("%Y-%m-%d_%H:%M:%S"), len(self.ts), conta["feito"], conta["sem_banda"], conta["falhou"],
            conta["pulado"])
        txt += "".join("%s %s: %s (%s)\n" % (t.config, t.caso, t.estado, t.motivo) for t in self.ts
                       if t.estado != "feito")
        escreve(os.path.join(self.base, "CONCLUIDO"), txt)
        print("niveis: " + txt, end="", flush=True)
        return 1 if conta["falhou"] else 0


def listar(ts):
    print("%-5s %-28s %-6s %-28s %s" % ("ordem", "configuracao", "caso", "pai", "estado"))
    for k, t in enumerate(ts):
        print("%-5d %-28s %-6s %-28s %s" % (k, t.config, t.caso, t.config_pai or "-", t.estado))
    conta = {}
    for t in ts:
        conta[t.estado] = conta.get(t.estado, 0) + 1
    print("%d tarefas: %s" % (len(ts), ", ".join("%d %s" % (n, e) for e, n in sorted(conta.items()))))


def rodar(a):
    base = os.path.abspath(a.saida)
    ts = tarefas(base, compara.casos(a.casos), extra=not a.sem_extra, so=a.configs.split(",") if a.configs else None)
    if a.listar:
        estado_inicial(ts)
        listar(ts)
        return 0
    if a.paralelo < 1:
        sys.exit("niveis: --paralelo deve ser >= 1")
    os.makedirs(base, exist_ok=True)
    trava = open(os.path.join(base, ".trava"), "w")
    try:
        fcntl.flock(trava, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit("niveis: outro niveis.py ja roda em %s" % base)
    exe, versao, ld = congela(base)
    copia_malhas(base)
    estado_inicial(ts)
    if os.path.exists(os.path.join(base, "CONCLUIDO")):
        os.remove(os.path.join(base, "CONCLUIDO"))
    c = Campanha(base, ts, a.paralelo, exe, versao, ld)
    pend = sum(t.estado == "pendente" for t in ts)
    with open(os.path.join(base, "sessoes.log"), "a") as f:
        f.write("sessao=%d inicio=%s pid=%d paralelo=%d casos=%s configs=%s extra=%d a_rodar=%d tarefas=%d "
                "versao=%s cpu=%s nucleos=%d\n" % (
                    c.sessao, time.strftime("%Y-%m-%d_%H:%M:%S"), os.getpid(), a.paralelo, a.casos, a.configs or "todas",
                    0 if a.sem_extra else 1, pend, len(ts), versao.replace(" ", "_"), cpu_modelo().replace(" ", "_"),
                    os.cpu_count() or 0))
    print("%s: %d de %d tarefas a rodar, %d ao mesmo tempo, programa %s (%s), SOLVER=%s" % (
        base, pend, len(ts), a.paralelo, exe, versao, SOLVER), flush=True)
    return c.roda()


# ---------------------------------------------------------------------------------------------------------------------
# tabela

COLUNAS = ["caso", "inicial", "tipo", "h", "config", "pai", "origem", "estado", "fator", "e_ref", "triangulos", "neq",
           "celulas_banda", "fator_pai", "seg_malha", "cpu_malha_s", "parede_s", "cpu_s", "seg", "cadeia_parede_s",
           "cadeia_cpu_s", "carga_media", "simultaneos", "newton", "newton_fs", "memoria_MB", "situacao", "dist_max",
           "fora", "inicio", "fim", "inicio_malha"]


def num(dic, k):
    try:
        return float(dic[k])
    except (KeyError, TypeError, ValueError):
        return None


def contencao(pasta):
    """(dist_max, fora) de uma rodada refinada: dos pontos da banda dela aos pontos da banda com que a malha dela foi
    feita (os Point de refina_banda.py no malha.geo)"""
    import numpy as np
    cel, folga = compara.celulas_refino(os.path.join(pasta, "malha.geo"))
    b = np.array(compara.banda(pasta, LIMIAR))
    c = np.array(cel)
    d = np.empty(len(b))
    for i in range(0, len(b), 2000):  # em blocos, para nao montar uma matriz grande demais
        bb = b[i:i + 2000]
        d[i:i + 2000] = np.sqrt(((bb[:, None, :] - c[None, :, :]) ** 2).sum(-1)).min(1)
    return float(d.max()), float((d > folga + 1e-9).mean())


def le_rodada(p):
    """os dados de uma pasta de rodada (desta campanha ou da comparacao)"""
    r, tm = compara.chaves(os.path.join(p, "resultado.txt")), compara.chaves(os.path.join(p, "tempo.txt"))
    info, mt = compara.chaves(os.path.join(p, "realizacao.txt")), compara.chaves(os.path.join(p, "malha.txt"))
    rs = compara.linha_csv(os.path.join(p, "resumo.csv"))
    cargas = [x for x in (num(tm, "carga_inicio"), num(tm, "carga_fim")) if x is not None]
    inteiro = lambda x: None if x is None else int(x)
    return dict(fator=num(r, "FS_min"), seg=num(r, "seg"), newton=inteiro(num(r, "newton")), situacao=r.get("situacao", ""),
                neq=inteiro(num(info, "neq")), triangulos=inteiro(num(mt, "triangulos")),
                celulas_banda=inteiro(num(mt, "celulas_banda")),
                fator_pai=num(mt, "fator_grossa"), seg_malha=num(tm, "seg_malha"), cpu_malha_s=num(tm, "cpu_malha_s"),
                parede_s=num(tm, "parede_s"), cpu_s=num(tm, "cpu_s"), inicio=num(tm, "inicio"), fim=num(tm, "fim"),
                inicio_malha=num(tm, "inicio_malha"), carga_media=sum(cargas) / len(cargas) if cargas else None,
                newton_fs=inteiro(num(rs, "newton_fs")), memoria_MB=num(rs, "memoria_MB"))


def casos_presentes(base):
    cs = set()
    for c in configuracoes(True):
        d = os.path.join(base, c[0])
        for n in (os.listdir(d) if os.path.isdir(d) else []):
            if n == "det":
                cs.add("det")
            elif n[:1] == "r" and n[1:].isdigit():
                cs.add(str(int(n[1:])))
    return sorted(cs, key=lambda c: -1 if c == "det" else int(c))


def simultaneos(linhas):
    """processos da mesma fila rodando junto com cada rodada (media no tempo do programa, contando ela mesma); o
    refina_banda.py ocupa o lugar na fila desde inicio_malha"""
    ivs = [(v["inicio_malha"] or v["inicio"], v["fim"]) for v in linhas if v["inicio"] is not None]
    for v in linhas:
        if v["inicio"] is not None and v["fim"] > v["inicio"]:
            sob = sum(max(0., min(v["fim"], b) - max(v["inicio"], a)) for a, b in ivs)
            v["simultaneos"] = sob / (v["fim"] - v["inicio"])


def tabela(a):
    base = os.path.abspath(a.saida)
    comp = os.path.abspath(a.comparacao) if a.comparacao else os.path.join(os.path.dirname(base), "refina",
                                                                           "comparacao")
    lista = casos_presentes(base)
    if not lista:
        sys.exit("niveis: nenhuma rodada em %s" % base)
    ts = tarefas(base, lista, extra=True)
    estado_inicial(ts)
    linhas, por = [], {}
    for t in ts:
        v = dict.fromkeys(COLUNAS)
        v.update(caso=t.caso, inicial=t.inicial, tipo=t.tipo, h=t.h, config=t.config, pai=t.config_pai or "",
                 origem="niveis", estado=t.estado)
        if t.estado == "pendente" and os.path.isfile(t.arq("falhou.txt")):
            v["estado"] = "falhou"
        if t.estado == "feito":
            v.update(le_rodada(t.pasta))
            if t.pai is not None:
                try:
                    v["dist_max"], v["fora"] = contencao(t.pasta)
                except (OSError, ValueError, IndexError) as e:
                    print("niveis: contencao de %s %s: %s" % (t.config, t.caso, e), file=sys.stderr)
        elif t.estado == "sem_banda":
            v.update({k: x for k, x in le_rodada(t.pasta).items() if k in ("seg_malha", "cpu_malha_s")})
        linhas.append(v)
        por[t.config, t.caso] = v
    # o custo do metodo: a soma das etapas da cadeia (cada refinada soma o seu refina_banda.py)
    for t in ts:  # a pai vem antes
        v = por[t.config, t.caso]
        for k, km, kp in (("cadeia_parede_s", "seg_malha", "parede_s"), ("cadeia_cpu_s", "cpu_malha_s", "cpu_s")):
            if v[kp] is None or (t.pai is not None and v[km] is None):
                continue
            proprio = v[kp] + (v[km] or 0.)
            v[k] = proprio if t.pai is None else (por[t.config_pai, t.caso][k] + proprio
                                                   if por[t.config_pai, t.caso][k] is not None else None)
    simultaneos([v for v in linhas if v["estado"] == "feito"])
    # as malhas fixas da comparacao
    fixas, refinada_comp = [], {}
    if os.path.isdir(comp):
        for m, h in FIXAS:
            tri = refina_banda.triangulos(os.path.join(MALHAS, "ceron_h5_a45_%s.msh" % m))
            for c in lista:
                p = os.path.join(comp, m, nome_caso(c))
                d = le_rodada(p)
                if d["fator"] is None:
                    continue
                v = dict.fromkeys(COLUNAS)
                v.update(d)
                v.update(caso=c, inicial=m, tipo="fixa", h=h, config="%s/fixa" % m, pai="", origem="comparacao",
                         estado="feito", triangulos=tri, cadeia_parede_s=d["parede_s"], cadeia_cpu_s=d["cpu_s"])
                fixas.append(v)
        for c in lista:  # o metodo de duas etapas da comparacao, para conferir foco375/um_passo_h0.125
            p = os.path.join(comp, "refinada", nome_caso(c))
            d, rf = le_rodada(p), compara.linha_csv(os.path.join(p, "refino.csv"))
            if d["fator"] is not None and num(rf, "seg_fina") and num(rf, "triangulos"):
                refinada_comp[c] = dict(fator=d["fator"], seg_fina=num(rf, "seg_fina"),
                                        triangulos=int(num(rf, "triangulos")))
        comp_todas = []
        for m in os.listdir(comp):
            if os.path.isdir(os.path.join(comp, m)):
                for c in os.listdir(os.path.join(comp, m)):
                    d = le_rodada(os.path.join(comp, m, c))
                    d["simultaneos"] = None
                    comp_todas.append(d)
        simultaneos(comp_todas)
        sim = {(round(d["inicio"], 3)): d["simultaneos"] for d in comp_todas if d["inicio"] is not None}
        for v in fixas:
            v["simultaneos"] = sim.get(round(v["inicio"], 3)) if v["inicio"] is not None else None
    else:
        print("niveis: sem a comparacao em %s: sem as malhas fixas" % comp, file=sys.stderr)
    linhas = fixas + linhas
    ref = {c: por[REF, c]["fator"] for c in lista if (REF, c) in por and por[REF, c]["fator"]}
    for v in linhas:
        v["e_ref"] = (v["fator"] - ref[v["caso"]]) / ref[v["caso"]] if v["fator"] is not None and v["caso"] in ref \
            else None
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(COLUNAS)
            for v in linhas:
                w.writerow(["" if v[k] is None else ("%g" % v[k] if k == "h" else v[k]) for k in COLUNAS])
    markdown(linhas, lista, ref, refinada_comp)


def markdown(linhas, lista, ref, refinada_comp):
    fmt, media = compara.fmt, compara.media
    metodos = []
    for v in linhas:
        if v["config"] not in metodos:
            metodos.append(v["config"])
    dados = {(v["config"], v["caso"]): v for v in linhas}
    ok = lambda m, c: (m, c) in dados and dados[m, c]["fator"] is not None
    print("Referência: %s (Γ do caso). e = (Γ − Γ_ref) / Γ_ref. Casos: %s; com referência: %d.\n" % (
        REF, ", ".join(lista), len(ref)))
    print("Fator (Γ) de cada caso:\n")
    print("| método | " + " | ".join(lista) + " |")
    print("|---|" + "---|" * len(lista))
    for m in metodos:
        print("| %s | %s |" % (m, " | ".join(fmt(dados[m, c]["fator"]) if ok(m, c) else "-" for c in lista)))
    print("\nDiferença para a referência, e (%), caso a caso:\n")
    print("| método | " + " | ".join(lista) + " |")
    print("|---|" + "---|" * len(lista))
    for m in metodos:
        cel = []
        for c in lista:
            e = dados[m, c]["e_ref"] if (m, c) in dados else None
            cel.append("-" if e is None else fmt(100 * e, 2, True))
        print("| %s | %s |" % (m, " | ".join(cel)))
    # resumo por metodo: os casos com referencia em que o metodo rodou
    tf = {c: dados["foco125/fixa", c] for c in lista if ok("foco125/fixa", c)}
    print("\nResumo por método (médias por caso; tempos da cadeia: etapa 1 + refina_banda.py + rodadas refinadas; a "
          "razão é para a foco125 fixa da comparação, nos mesmos casos):\n")
    print("| método | casos | triângulos | células de banda | \\|e\\| médio | \\|e\\| máximo | e médio | "
          "parede da cadeia (s) | CPU da cadeia (s) | parede / foco125 | CPU / foco125 | processos juntos | carga | "
          "banda fora da região refinada |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for m in metodos:
        vs = [dados[m, c] for c in lista if ok(m, c) and dados[m, c]["e_ref"] is not None]
        if not vs:
            continue
        es = [v["e_ref"] for v in vs]
        tp, tc = media([v["cadeia_parede_s"] for v in vs]), media([v["cadeia_cpu_s"] for v in vs])
        # razao das medias = razao das somas, nos casos com a foco125 e com os dois tempos da cadeia
        comf = [v for v in vs if v["caso"] in tf and v["cadeia_parede_s"] is not None and v["cadeia_cpu_s"] is not None]
        rp = sum(v["cadeia_parede_s"] for v in comf) / sum(tf[v["caso"]]["parede_s"] for v in comf) if comf else None
        rc = sum(v["cadeia_cpu_s"] for v in comf) / sum(tf[v["caso"]]["cpu_s"] for v in comf) if comf else None
        fr = [v for v in vs if v["fora"] is not None]
        banda = "%d de %d (no pior, %s%% dos pontos; distância máxima %s m)" % (
            sum(v["fora"] > 0 for v in fr), len(fr), fmt(100 * max(v["fora"] for v in fr), 1),
            fmt(max(v["dist_max"] for v in fr), 2)) if fr else "-"
        print("| %s | %d | %s | %s | %s%% | %s%% | %s%% | %s | %s | %s | %s | %s | %s | %s |" % (
            m, len(vs), fmt(media([v["triangulos"] for v in vs]), 0), fmt(media([v["celulas_banda"] for v in vs]), 0),
            fmt(100 * media([abs(e) for e in es]), 2), fmt(100 * max(abs(e) for e in es), 2), fmt(100 * media(es), 2, True),
            fmt(tp, 1), fmt(tc, 1), fmt(rp, 3), fmt(rc, 3), fmt(media([v["simultaneos"] for v in vs]), 2),
            fmt(media([v["carga_media"] for v in vs]), 2), banda))
    # as etapas de cada cadeia (tempo de parede, medias por caso)
    print("\nEtapas da cadeia (tempo de parede, média por caso, s): a rodada na malha inicial e, em cada nível, o "
          "refina_banda.py (com o gmsh) e a rodada refinada:\n")
    print("| método | casos | etapa 1 | malha 1 | refinada 1 | malha 2 | refinada 2 | malha 3 | refinada 3 | total |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    pais = {c[0]: c[4] for c in configuracoes(True)}
    for m in metodos:
        if m not in pais:
            continue
        cadeia = [m]
        while pais[cadeia[0]]:
            cadeia.insert(0, pais[cadeia[0]])
        cs = [c for c in lista if all(ok(x, c) for x in cadeia)]
        if not cs:
            continue
        cel = [fmt(media([dados[cadeia[0], c]["parede_s"] for c in cs]), 1)]
        for x in cadeia[1:]:
            cel += [fmt(media([dados[x, c]["seg_malha"] for c in cs]), 2),
                    fmt(media([dados[x, c]["parede_s"] for c in cs]), 1)]
        cel += ["-"] * (7 - len(cel))
        print("| %s | %d | %s | %s |" % (m, len(cs), " | ".join(cel),
                                         fmt(media([dados[m, c]["cadeia_parede_s"] for c in cs]), 1)))
    # conferencias com a comparacao (outro executavel e talvez outra maquina) e a extra
    print()
    soma = lambda m, k, cs: sum(dados[m, c][k] for c in cs)
    e1 = lambda x: ("%.1e" % x).replace(".", ",")
    m, mc = "foco375/inicial", "foco375/fixa"
    cs = [c for c in lista if ok(m, c) and ok(mc, c)]
    if cs:
        print("- %s e a %s da comparação, %d casos: maior |ΔΓ|/Γ = %s; aqui / comparação: seg do programa %s, parede "
              "%s, CPU %s." % (m, mc, len(cs), e1(max(abs(dados[m, c]["fator"] / dados[mc, c]["fator"] - 1) for c in cs)),
                               fmt(soma(m, "seg", cs) / soma(mc, "seg", cs), 3),
                               fmt(soma(m, "parede_s", cs) / soma(mc, "parede_s", cs), 3),
                               fmt(soma(m, "cpu_s", cs) / soma(mc, "cpu_s", cs), 3)))
    m = "foco375/um_passo_h0.125"
    cs = [c for c in lista if ok(m, c) and c in refinada_comp]
    if cs:
        print("- %s e a `refinada` da comparação, %d casos: maior |ΔΓ|/Γ = %s; o mesmo número de triângulos em %d; "
              "seg da rodada refinada, aqui / comparação: %s." % (
                  m, len(cs), e1(max(abs(dados[m, c]["fator"] / refinada_comp[c]["fator"] - 1) for c in cs)),
                  sum(dados[m, c]["triangulos"] == refinada_comp[c]["triangulos"] for c in cs),
                  fmt(soma(m, "seg", cs) / sum(refinada_comp[c]["seg_fina"] for c in cs), 3)))
    mx = "foco375/extra_h%g" % H_EXTRA
    if ok(mx, "det") and ok(REF, "det"):
        print("- det: Γ da referência %s, Γ da extra (h = %s m) %s: a referência fica %s%% acima da extra." % (
            fmt(dados[REF, "det"]["fator"]), fmt(H_EXTRA, 6), fmt(dados[mx, "det"]["fator"]),
            fmt(100 * (dados[REF, "det"]["fator"] / dados[mx, "det"]["fator"] - 1), 3, True)))
    nao = [v for v in linhas if v["origem"] == "niveis" and v["estado"] != "feito"]
    if nao:
        print("- Rodadas não feitas: " + "; ".join("%s %s (%s)" % (v["config"], v["caso"], v["estado"]) for v in nao) + ".")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("rodar")
    r.add_argument("saida")
    r.add_argument("--casos", default="det,0-7")
    r.add_argument("--paralelo", type=int, default=2)
    r.add_argument("--sem-extra", action="store_true")
    r.add_argument("--configs", help="so estas configuracoes (e as ancestrais), separadas por virgula")
    r.add_argument("--listar", action="store_true")
    t = sub.add_parser("tabela")
    t.add_argument("saida")
    t.add_argument("--csv")
    t.add_argument("--comparacao")
    a = ap.parse_args()
    sys.exit(rodar(a) if a.cmd == "rodar" else tabela(a))


if __name__ == "__main__":
    main()
