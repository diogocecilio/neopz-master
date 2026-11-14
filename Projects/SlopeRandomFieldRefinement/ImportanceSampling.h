// Amostragem por importancia (AI) nos campos aleatorios do RandomField: a direcao da banda de ruptura de um
// deterministico, o piloto (um Monte Carlo direto ja feito) e o deslocamento da media dos xi pela entropia cruzada
// (Rubinstein e Kroese, 2004). O peso exato de cada realizacao e RandomField::LogWeight(), e o estimador
// P_f = soma(peso * falha) / N nao tem vies para qualquer deslocamento fixado antes da amostra
// (posproc/amostragem_importancia.py; relatorio_chuva.tex, Secao 7).
// Usado pelo montecarlo.cpp (chuva). O rebaixamento_mc.cpp tem a mesma logica no proprio arquivo (a campanha do Ceron
// no fermi usa aquela versao; juntar as duas depois que ela acabar).
// Com varios instantes por realizacao (chuva), a falha e FS < 1 em algum instante: no piloto vale o menor FS_min dos
// instantes, e a realizacao so conta se todos os instantes tiverem fator.
#pragma once

#include "RandomField.h"

#include <dirent.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ai {

// os campos de uma linha de CSV
inline std::vector<std::string> Divide(const std::string &linha)
{
    std::vector<std::string> v;
    std::stringstream sl(linha);
    for (std::string c; std::getline(sl, c, ',');) v.push_back(c);
    return v;
}

// pontos da banda do mecanismo.csv de um deterministico (colunas x, y, epsp_raizJ2 e, se houver, t_h): os do instante t
// (t < 0: todos) com sqrt(J2) >= 0,1 do maior, com peso 1. Vazio se o arquivo nao abre ou nao tem o instante
inline std::vector<std::array<double, 3>> LeBanda(const std::string &arq, double t = -1.)
{
    std::ifstream f(arq);
    std::string linha;
    if (!std::getline(f, linha)) return {};
    const std::vector<std::string> cab = Divide(linha);
    auto col = [&cab](const char *nome) { return int(std::find(cab.begin(), cab.end(), nome) - cab.begin()); };
    const int ix = col("x"), iy = col("y"), ij = col("epsp_raizJ2"), it = col("t_h");
    const bool temT = it < int(cab.size());
    std::vector<std::array<double, 3>> p; // x, y, sqrt(J2)
    while (std::getline(f, linha)) {
        const std::vector<std::string> c = Divide(linha);
        if (int(c.size()) <= std::max({ix, iy, ij}) || (temT && int(c.size()) <= it)) continue;
        if (t >= 0. && temT && std::abs(std::stod(c[it]) - t) > 1e-6) continue;
        p.push_back({std::stod(c[ix]), std::stod(c[iy]), std::stod(c[ij])});
    }
    double jmax = 0.;
    for (const auto &q : p) jmax = std::max(jmax, q[2]);
    std::vector<std::array<double, 3>> banda;
    for (const auto &q : p)
        if (q[2] >= 0.1 * jmax) banda.push_back({q[0], q[1], 1.});
    return banda;
}

// O piloto: as realizacoes 0..np-1 do Monte Carlo direto da pasta (os <prefixo>_<a>_<b>.csv da fila do cluster, com as
// colunas real, FS_min e situacao; com varios instantes, o menor FS_min deles) e o conjunto elite da entropia cruzada: as
// que falharam (FS < 1) ou, se forem menos que a fracao rho das validas, as rho com menor FS (um nivel intermediario:
// FS <= nivel). Ordenado por (FS, realizacao): o mesmo em todos os blocos que leem o mesmo piloto completo
struct Piloto {
    int presentes = 0, validas = 0, falhas = 0;
    double nivel = 1.;
    std::vector<int64_t> elite;
};

inline Piloto LePiloto(const std::string &pasta, int np, double rho)
{
    std::map<int64_t, double> g; // realizacao -> menor FS_min dos instantes (-1: algum instante sem fator)
    DIR *d = opendir(pasta.c_str());
    if (!d) {
        std::fprintf(stderr, "piloto: nao abri a pasta %s\n", pasta.c_str());
        return {};
    }
    for (dirent *e; (e = readdir(d));) {
        const std::string nome = e->d_name;
        // <prefixo>_<a>_<b>.csv com o prefixo sem '_' (nao pega mc_mapa_*, mc_proj_*)
        if (nome.size() < 5 || nome.compare(nome.size() - 4, 4, ".csv")) continue;
        const std::string base = nome.substr(0, nome.size() - 4);
        const size_t u = base.find('_');
        if (u == std::string::npos || std::count(base.begin(), base.end(), '_') != 2 ||
            base.find_first_not_of("0123456789_", u) != std::string::npos)
            continue;
        std::ifstream f(pasta + "/" + nome);
        std::string linha;
        if (!std::getline(f, linha)) continue;
        const std::vector<std::string> cab = Divide(linha);
        auto col = [&cab](const char *c) { return int(std::find(cab.begin(), cab.end(), c) - cab.begin()); };
        const int ir = col("real"), ig = col("FS_min"), is = col("situacao");
        if (std::max({ir, ig, is}) >= int(cab.size())) continue;
        while (std::getline(f, linha)) {
            const std::vector<std::string> v = Divide(linha);
            if (int(v.size()) <= std::max({ir, ig, is})) continue; // linha incompleta
            try {
                const int64_t r = std::stoll(v[ir]);
                if (r < 0 || r >= np) continue;
                const double fs = std::stod(v[ig]);
                const bool ok = v[is] == "ok" && fs > 0.;
                const auto it = g.find(r);
                if (it == g.end()) g[r] = ok ? fs : -1.;
                else if (it->second > 0.) it->second = ok ? std::min(it->second, fs) : -1.;
            } catch (const std::exception &) {
            }
        }
    }
    closedir(d);
    Piloto p;
    p.presentes = g.size();
    std::vector<std::pair<double, int64_t>> v;
    for (const auto &x : g)
        if (x.second > 0.) v.push_back({x.second, x.first});
    std::sort(v.begin(), v.end());
    p.validas = v.size();
    p.falhas = std::count_if(v.begin(), v.end(), [](const std::pair<double, int64_t> &x) { return x.first < 1.; });
    const int m = std::min(p.validas, std::max(1, int(std::ceil(rho * p.validas))));
    const int ne = p.falhas >= m ? p.falhas : m;
    for (int i = 0; i < ne; i++) p.elite.push_back(v[i].second);
    if (p.falhas < m && m > 0) p.nivel = v[m - 1].first;
    return p;
}

// opcoes da AI pelo piloto: a pasta do Monte Carlo direto, quantas realizacoes dele, a semente dele, o tipo de
// deslocamento (banda: s_c d por campo, com d a direcao da banda; modos: os 'modos' primeiros modos de cada campo) e a
// fracao elite
struct Opcoes {
    std::string piloto, tipo = "banda";
    int n = 1000, modos = 30;
    uint64_t semente = 0;
    double elite = 0.1;
};

// O deslocamento da entropia cruzada pelo piloto: a media dos xi do conjunto elite, sorteados de novo com a semente do
// piloto (sem deslocamento: chame antes de qualquer SetShift), aplicada em campos (SetShift). Grava o registro em
// 'registro'. Devolve false se o piloto tiver menos de 90% das realizacoes 0..n-1 ou nenhuma valida
inline bool DeslocaPeloPiloto(RandomField &campos, const std::vector<double> &dir, const Opcoes &o,
                              const std::string &registro)
{
    const Piloto pl = LePiloto(o.piloto, o.n, o.elite);
    if (pl.presentes < 0.9 * o.n || pl.elite.empty()) {
        std::fprintf(stderr, "piloto incompleto em %s: %d das realizacoes 0..%d, %d validas\n", o.piloto.c_str(),
                     pl.presentes, o.n - 1, pl.validas);
        return false;
    }
    const int nc = campos.NCampos();
    const double ne = pl.elite.size();
    double mu2 = 0.;
    std::ostringstream reg; // o registro do deslocamento
    reg << std::setprecision(10);
    if (o.tipo == "banda") { // s_c = media de xi_c . d
        std::vector<double> s(nc, 0.);
        for (int64_t r : pl.elite) {
            campos.Sample(o.semente, r);
            const std::vector<double> pr = campos.Projection(dir);
            for (int c = 0; c < nc; c++) s[c] += pr[c] / ne;
        }
        campos.SetShift(dir, s);
        reg << "s=";
        for (int c = 0; c < nc; c++) reg << (c ? "," : "") << s[c], mu2 += s[c] * s[c];
        reg << "\nmu2=" << mu2 << "\n";
    } else { // a media dos K primeiros modos de cada campo
        const int K = std::min(o.modos, campos.M());
        std::vector<std::vector<double>> mu(nc, std::vector<double>(K, 0.));
        for (int64_t r : pl.elite) {
            campos.Sample(o.semente, r);
            for (int c = 0; c < nc; c++)
                for (int k = 0; k < K; k++) mu[c][k] += campos.Xi()[c][k] / ne;
        }
        campos.SetShift(mu);
        for (const auto &m : mu)
            for (double x : m) mu2 += x * x;
        reg << "modos=" << K << "\nmu2=" << mu2 << "\ntermo";
        for (int c = 0; c < nc; c++) reg << ",campo" << c;
        reg << "\n";
        for (int k = 0; k < K; k++) {
            reg << k;
            for (int c = 0; c < nc; c++) reg << "," << mu[c][k];
            reg << "\n";
        }
    }
    std::ofstream df(registro);
    df << "# deslocamento da amostragem por importancia: a media dos xi do conjunto elite do piloto\nai=" << o.tipo
       << "\npiloto=" << o.piloto << "\npiloto_n=" << o.n << "\npiloto_semente=" << o.semente
       << "\npresentes=" << pl.presentes << "\nvalidas=" << pl.validas << "\nfalhas=" << pl.falhas
       << "\nelite=" << pl.elite.size() << "\nnivel=" << pl.nivel << "\n"
       << reg.str();
    std::printf("  amostragem por importancia (%s) pelo piloto %s: %d validas de %d, %d falhas, elite de %zu "
                "(FS < %.4f); |mu|^2 = %.4f\n",
                o.tipo.c_str(), o.piloto.c_str(), pl.validas, o.n, pl.falhas, pl.elite.size(), pl.nivel, mu2);
    return true;
}

} // namespace ai
