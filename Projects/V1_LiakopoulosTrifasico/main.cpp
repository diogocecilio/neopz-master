// V1_LiakopoulosTrifasico: executavel do caso completo (ver README.md). Coluna de Liakopoulos drenando ate 600 min com
// ar passivo (bifasico) e ar ativo (trifasico, trifasico_laloui), como na Secao 5 de Laloui, Klubertanz e Vulliet
// (2003). Grava na pasta corrente:
//   historico.csv  uma linha por passo: vazao de agua na base, vazao de ar no topo e na base, recalque do topo, p_w no
//                  meio e no topo, minimo de p_a e a sua altura;
//   perfis.csv     p_w, p_a, succao, S_w e u_y na vertical x = 0,05 m nos instantes dos perfis;
//   resumo.csv     uma linha por caso com os valores comparados com as figuras do artigo;
//   resultado.txt  "chave=valor" (situacao=ok, seg=...), gravado por ultimo.
// Variaveis de ambiente:
//   CASO      bifasico | trifasico | trifasico_laloui | todos (padrao; tambem uma lista separada por virgulas)
//   T_FIM     fim do transiente em s (padrao 36000 = 600 min)
//   NY        0 = malha do Gmsh 4 x 40 (padrao); n > 0 = malha estruturada 4 x n (estudo de malha)
//   DT_FATOR  fator dos passos de tempo (padrao 1; 0.5 = metade dos passos)
//   VTK       1 grava VTK nos instantes dos perfis (padrao 0)
#include "problema.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string Env(const char *nome, const std::string &padrao)
{
    const char *v = std::getenv(nome);
    return (v && *v) ? std::string(v) : padrao;
}
} // namespace

int main()
{
    using namespace v1;
    const std::string casoEnv = Env("CASO", "todos");
    Opcoes op;
    op.tFim = std::stod(Env("T_FIM", "36000"));
    op.ny = std::stoi(Env("NY", "0"));
    op.passos = PassosReferencia(std::stod(Env("DT_FATOR", "1")));
    op.vtk = Env("VTK", "0") == "1";

    std::vector<std::string> nomes;
    if (casoEnv == "todos") {
        for (const Caso &c : Casos()) nomes.push_back(c.nome);
    } else {
        std::stringstream ss(casoEnv);
        std::string n;
        while (std::getline(ss, n, ','))
            if (!n.empty()) nomes.push_back(n);
    }
    std::cout << "V1_LiakopoulosTrifasico: casos " << casoEnv << ", t_fim = " << op.tFim << " s, malha "
              << (op.ny > 0 ? "estruturada 4 x " + std::to_string(op.ny) : std::string("do Gmsh 4 x 40"))
              << ", fator dos passos " << op.passos.fator << std::endl;

    std::ofstream hist("historico.csv"), prof("perfis.csv"), resumo("resumo.csv");
    CabecalhoHistorico(hist);
    CabecalhoPerfis(prof);
    // resumo: os valores comparados com as figuras (NaN quando o instante nao foi simulado)
    resumo << "caso,ok,passos,seg,neq,newton_por_passo,cortes,pontos_plasticos,qb_10min_1e6,qb_30min_1e6,qb_60min_1e6,"
              "pw_topo_5min,pw_topo_10min,pw_topo_30min,pw_topo_600min,w_topo_600min_mm,pa_min_5min,y_pa_min_5min,"
              "pa_min_30min,y_pa_min_30min,pa_min_600min,y_pa_min_600min\n";
    bool ok = true;
    double segTotal = 0.;
    std::string paMin30 = "nan";
    for (const std::string &nome : nomes) {
        const Caso caso = CasoPorNome(nome);
        const Resultado r = RodaCaso(caso, op);
        GravaHistorico(hist, r);
        GravaPerfis(prof, r);
        hist.flush();
        prof.flush();
        ok = ok && r.ok;
        segTotal += r.seg;
        auto em = [&](REAL t, auto f) -> std::string {
            const Perfil *p = r.PerfilEm(t);
            if (!p) return "nan";
            std::ostringstream o;
            o.precision(6);
            o << f(*p);
            return o.str();
        };
        auto qb = [&](REAL t) -> std::string {
            if (r.hist.empty() || r.hist.back().t < t - 1.e-6) return "nan";
            std::ostringstream o;
            o.precision(6);
            o << r.HistEm(&Linha::qb, t) * 1.e6;
            return o.str();
        };
        auto topo = [](const Perfil &p) { return p.pw.back(); };
        auto pamin = [](const Perfil &p) { return p.Minimo(p.pa).first; };
        auto ypamin = [](const Perfil &p) { return p.Minimo(p.pa).second; };
        resumo << r.caso << ',' << (r.ok ? 1 : 0) << ',' << r.passos << ',' << r.seg << ',' << r.neq << ','
               << REAL(r.newtonTotal) / std::max(r.passos, 1) << ',' << r.cortes << ',' << r.pontosPlasticos << ','
               << qb(600.) << ',' << qb(1800.) << ',' << qb(3600.) << ',' << em(300., topo) << ',' << em(600., topo)
               << ',' << em(1800., topo) << ',' << em(36000., topo) << ','
               << em(36000., [](const Perfil &p) { return p.uy.back() * 1.e3; }) << ',' << em(300., pamin) << ','
               << em(300., ypamin) << ',' << em(1800., pamin) << ',' << em(1800., ypamin) << ',' << em(36000., pamin)
               << ',' << em(36000., ypamin) << '\n';
        resumo.flush();
        if (caso.arAtivo && paMin30 == "nan") paMin30 = em(1800., pamin);
    }
    // resultado.txt por ultimo (um .tmp renomeado)
    {
        std::ofstream o("resultado.txt.tmp");
        o << "casos=" << casoEnv << " situacao=" << (ok ? "ok" : "falhou") << " pa_min_30min=" << paMin30
          << " seg=" << segTotal << '\n';
    }
    std::rename("resultado.txt.tmp", "resultado.txt");
    return ok ? 0 : 1;
}
