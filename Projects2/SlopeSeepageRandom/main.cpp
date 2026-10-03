// main.cpp — SlopeSeepageRandom
//
// Reprodução, por elementos finitos elastoplásticos no NeoPZ, de
//   M. Vargas Ceron, D. L. Cecílio, R. V. Linn, S. Maghous, "Stability Analysis of Slope Subjected to Seepage
//   Forces Considering Spatial Variability of Soil Properties", Int J Numer Anal Methods Geomech 49 (2025)
//   2459-2491,
// com Mohr-Coulomb (TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>) e Cam-Clay modificado (TPZModifiedCamClay), campos
// aleatórios de c, φ e kv pela expansão de Karhunen-Loève (TPZMatKLKernel + pzdoublestrmatriz + LAPACK) e
// forças de percolação do problema de Darcy desacoplado (TPZDarcyFlow).
//
// Uso: SlopeSeepageRandom <comando> [opções]   (ver Usage())
//
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <vector>

#include <unistd.h>

#include "CoupledDrawdown.h"
#include "FileIO.h"
#include "KLRandomField.h"
#include "SeepageProblem.h"
#include "SlopeGeometry.h"
#include "SlopeStability.h"
#include "TPZBFileStream.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeoel.h"

namespace {

using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

/// Opções de linha de comando no formato chave=valor
struct TArgs {
    std::map<std::string, std::string> kv;
    TArgs(int argc, char **argv, int first) {
        for (int i = first; i < argc; i++) {
            std::string a = argv[i];
            auto p = a.find('=');
            if (p == std::string::npos) kv[a] = "1";
            else kv[a.substr(0, p)] = a.substr(p + 1);
        }
    }
    REAL Get(const std::string &k, REAL def) const {
        auto it = kv.find(k);
        return it == kv.end() ? def : std::atof(it->second.c_str());
    }
    int GetI(const std::string &k, int def) const {
        auto it = kv.find(k);
        return it == kv.end() ? def : std::atoi(it->second.c_str());
    }
    std::string GetS(const std::string &k, const std::string &def) const {
        auto it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
};

/// Gradiente do excesso de poropressão do problema de Darcy nos pontos de integração da malha mecânica
template <class T>
void TransferSeepage(TSeepageProblem &seep, TSlopeFEM<T> &fem) {
    const auto &pts = fem.Points();
    std::vector<REAL> u(pts.size(), 0.);
    std::vector<TPZManVector<REAL, 2>> g(pts.size(), TPZManVector<REAL, 2>(2, 0.));
    for (size_t i = 0; i < pts.size(); i++) {
        if (pts[i].gel < 0) continue;
        seep.Evaluate(pts[i].gel, pts[i].qsi, u[i], g[i]);
    }
    fem.SetSeepage(u, g);
}

struct TCase {
    std::string name;
    TSlopeGeometry geo;
    TSoil soil;
    bool seepage = false;
    REAL hw = 5.;
    REAL alpha = 1.;
};

TCase MakeCase(const std::string &name, REAL h) {
    TCase c;
    c.name = name;
    if (name == "cho_coesivo") {  // Cho (2010) / seção 5.3.1: cu = 23 kPa, φu = 0, γ = 20, 2:1, H = 5 m
        c.geo = TSlopeGeometry::Cho2H1V(h);
        c.soil.c = 23.;
        c.soil.phiDeg = 0.;
        c.soil.gamma = 20.;
        c.soil.buoyant = false;
    } else if (name == "cho_cphi") {
        // seção 5.3.2 / Cho (2010): c = 10 kPa, φ = 30°, γ = 20, 1:1 (seco). O FS = 1.204 de Cho (e Γ = 1.777 do
        // artigo) corresponde a H = 10 m: com H = 5 m o próprio Bishop simplificado dá FS = 1.61 (ver README).
        c.geo = TSlopeGeometry::Cho1H1V(h);
        c.geo.H = 10.;
        c.soil.c = 10.;
        c.soil.phiDeg = 30.;
        c.soil.gamma = 20.;
        c.soil.buoyant = false;
    } else {  // "percolacao": seção 6 (Tabela 2): c = 10, φ = 30°, γ = 20, β = 45°, H = hw = 5 m, α = 1
        c.geo = TSlopeGeometry::Cho1H1V(h);
        c.soil.c = 10.;
        c.soil.phiDeg = 30.;
        c.soil.gamma = 20.;
        c.soil.buoyant = true;
        c.seepage = true;
    }
    return c;
}

// ============================================================ gravação e retomada (read/write)

/// Erro com código de saída (1: uso, 2: arquivo inválido, 3: parâmetros diferentes dos gravados, 4: CSV em uso)
struct TExitError : std::runtime_error {
    int code;
    TExitError(const std::string &msg, int c) : std::runtime_error(msg), code(c) {}
};

/// Real com 15 algarismos significativos (texto dos parâmetros, comparado entre execuções)
std::string Fmt(REAL v) {
    std::ostringstream s;
    s << std::setprecision(15) << v;
    return s.str();
}

std::string Now() {
    const std::time_t t = std::time(nullptr);
    char b[32];
    std::strftime(b, sizeof(b), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    return b;
}

std::vector<std::string> Split(const std::string &s, char sep) {
    std::vector<std::string> v;
    size_t start = 0;
    while (true) {
        const size_t p = s.find(sep, start);
        v.push_back(s.substr(start, p == std::string::npos ? std::string::npos : p - start));
        if (p == std::string::npos) return v;
        start = p + 1;
    }
}

bool ParseInteger(const std::string &s, int64_t &v) {
    if (s.empty() || std::isspace((unsigned char)s[0])) return false;
    char *end = nullptr;
    errno = 0;
    const long long x = std::strtoll(s.c_str(), &end, 10);
    if (errno != 0 || *end != '\0') return false;
    v = (int64_t)x;
    return true;
}

bool ParseReal(const std::string &s, REAL &v) {
    if (s.empty() || std::isspace((unsigned char)s[0])) return false;
    char *end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return *end == '\0';
}

/// Lista de amostras "3,17,40-45" (ordenada, sem repetições)
std::vector<int64_t> ParseSampleList(const std::string &str) {
    std::set<int64_t> ids;
    for (const std::string &item : Split(str, ',')) {
        if (item.empty()) continue;
        const size_t d = item.find('-', 1);
        int64_t a = -1, b = -1;
        const bool ok = (d == std::string::npos) ? (ParseInteger(item, a) && (b = a) >= 0)
                                                 : (ParseInteger(item.substr(0, d), a) &&
                                                    ParseInteger(item.substr(d + 1), b));
        if (!ok || a < 0 || b < a || b - a > 10000000)
            throw TExitError("amostras=: item inválido \"" + item + "\" (ex.: amostras=3,17,40-45)", 1);
        for (int64_t s = a; s <= b; s++) ids.insert(s);
    }
    return std::vector<int64_t>(ids.begin(), ids.end());
}

// ----- CSV do Monte Carlo: uma linha por amostra, gravada inteira (uma escrita + fdatasync)
const char *const kCsvHeader =
    "amostra,fator,limite_superior,status,passos,iteracoes,cortes,tempo_s,c_medio,phi_medio,kv_medio";
const int kCsvTimeColumn = 7;  ///< tempo_s: a única coluna que muda quando a amostra é recalculada

struct TMCRow {
    int64_t sample = -1;
    REAL factor = 0., upper = 0., time = 0.;
    std::string status;
    std::vector<std::string> fields;  ///< os 11 campos, como gravados
    std::string Line() const {
        std::string l;
        for (size_t i = 0; i < fields.size(); i++) l += (i ? "," : "") + fields[i];
        return l;
    }
};

/// Linha completa e válida: 11 campos, amostra inteira >= 0, status [A-Za-z0-9_]+, inteiros (passos, iterações,
/// cortes) e reais (aceita nan/inf) nos demais
bool ParseRow(const std::string &line, TMCRow &r) {
    std::vector<std::string> f = Split(line, ',');
    if (f.size() != 11) return false;
    int64_t iv;
    REAL rv;
    if (!ParseInteger(f[0], r.sample) || r.sample < 0 || f[3].empty()) return false;
    for (char ch : f[3])
        if (!std::isalnum((unsigned char)ch) && ch != '_') return false;
    for (int k : {4, 5, 6})
        if (!ParseInteger(f[k], iv)) return false;
    for (int k : {1, 2, 7, 8, 9, 10})
        if (!ParseReal(f[k], rv)) return false;
    ParseReal(f[1], r.factor);
    ParseReal(f[2], r.upper);
    ParseReal(f[kCsvTimeColumn], r.time);
    r.status = f[3];
    r.fields = std::move(f);
    return true;
}

/// Linha do CSV (mesmo texto que as versões anteriores: precisão 8)
TMCRow MakeRow(int64_t s, const TFactorResult &r, double dt, REAL cm, REAL pm, REAL km) {
    std::ostringstream o;
    o << std::setprecision(8) << s << "," << r.factor << "," << r.upper << "," << (r.status.empty() ? "ok" : r.status)
      << "," << r.steps << "," << r.iterations << "," << r.cuts << "," << dt << "," << cm << "," << pm << "," << km;
    TMCRow row;
    if (!ParseRow(o.str(), row)) throw std::runtime_error("linha do CSV inválida: " + o.str());
    return row;
}

struct TCsvContent {
    bool exists = false, headerOk = false;
    std::vector<TMCRow> rows;                                  ///< válidas, sem repetição, na ordem do arquivo
    std::vector<std::pair<std::string, std::string>> dropped;  ///< (linha, motivo)
    bool NeedsRewrite() const { return !exists || !headerOk || !dropped.empty(); }
};

/// Lê o CSV: mantém as linhas completas e válidas (a primeira de cada amostra); as demais (linha final sem '\n'
/// de um processo interrompido, lixo, repetições) vão para "dropped"
TCsvContent ReadCsv(const std::string &csv) {
    TCsvContent c;
    if (!fileio::Exists(csv)) return c;
    std::string text;
    if (!fileio::ReadAll(csv, text)) throw TExitError("não foi possível ler " + csv, 2);
    c.exists = true;
    std::set<int64_t> seen;
    size_t pos = 0;
    bool first = true;
    while (pos < text.size()) {
        const size_t nl = text.find('\n', pos);
        const bool complete = nl != std::string::npos;
        std::string line = text.substr(pos, complete ? nl - pos : std::string::npos);
        pos = complete ? nl + 1 : text.size();
        if (complete && !line.empty() && line.back() == '\r') line.pop_back();
        if (first) {
            first = false;
            const std::string h = kCsvHeader;
            if (complete && line == h) {
                c.headerOk = true;
                continue;
            }
            if (!complete && h.compare(0, line.size(), line) == 0) continue;  // cabeçalho interrompido
            throw TExitError("o arquivo " + csv + " não começa pelo cabeçalho do comando mc (" + h +
                                 "); nada foi alterado",
                             2);
        }
        TMCRow r;
        if (!complete) c.dropped.emplace_back(line, "linha final incompleta (processo interrompido)");
        else if (!ParseRow(line, r)) c.dropped.emplace_back(line, "linha inválida");
        else if (!seen.insert(r.sample).second) c.dropped.emplace_back(line, "amostra repetida (mantida a primeira)");
        else c.rows.push_back(std::move(r));
    }
    return c;
}

/// Regrava o CSV (cabeçalho + linhas válidas) atomicamente; as linhas retiradas são guardadas em <csv>.descartadas
void RewriteCsv(const std::string &csv, const TCsvContent &c) {
    if (!c.dropped.empty()) {
        std::ostringstream d;
        d << "# " << Now() << ": " << c.dropped.size() << " linha(s) retirada(s) de " << csv << "\n";
        for (auto &[line, why] : c.dropped) d << "# " << why << ":\n" << line << "\n";
        fileio::TAppendFile(csv + ".descartadas").Append(d.str());
        std::cout << "[MC] " << csv << ": " << c.dropped.size()
                  << " linha(s) incompleta(s), inválida(s) ou repetida(s) retirada(s) (guardadas em " << csv
                  << ".descartadas)\n";
    }
    std::string s = std::string(kCsvHeader) + "\n";
    for (const TMCRow &r : c.rows) s += r.Line() + "\n";
    fileio::WriteAtomic(csv, s);
}

// ----- parâmetros da execução (<csv>.param): tudo o que muda as amostras, "chave=valor"
using TParamList = std::vector<std::pair<std::string, std::string>>;

bool ReadParamFile(const std::string &file, std::map<std::string, std::string> &kv, std::string &raw) {
    if (!fileio::Exists(file)) return false;
    if (!fileio::ReadAll(file, raw)) throw TExitError("não foi possível ler " + file, 2);
    std::istringstream in(raw);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const size_t p = line.find('=');
        if (p != std::string::npos) kv[line.substr(0, p)] = line.substr(p + 1);
    }
    return true;
}

std::string ParamText(const std::string &csv, const TParamList &p) {
    std::ostringstream s;
    s << "# SlopeSeepageRandom mc: parâmetros das amostras de " << csv << " (gravado em " << Now() << ")\n"
      << "# conferidos a cada retomada; uma diferença interrompe a execução (forcar=1 aceita e anota aqui)\n";
    for (auto &[k, v] : p) s << k << "=" << v << "\n";
    return s.str();
}

/// Diferenças entre os parâmetros gravados e os atuais. requireAll: chave atual ausente do arquivo também é
/// diferença; extraSaved: chaves do arquivo que podem não estar em cur (as demais, ausentes de cur, são diferenças)
std::vector<std::string> DiffParams(const std::map<std::string, std::string> &saved, const TParamList &cur,
                                    bool requireAll, const std::set<std::string> *extraSaved) {
    std::vector<std::string> d;
    std::set<std::string> keys;
    for (auto &[k, v] : cur) {
        keys.insert(k);
        auto it = saved.find(k);
        if (it == saved.end()) {
            if (requireAll) d.push_back(k + ": ausente do arquivo, agora " + v);
        } else if (it->second != v) {
            d.push_back(k + ": " + it->second + " no arquivo, agora " + v);
        }
    }
    if (extraSaved)
        for (auto &[k, v] : saved)
            if (!keys.count(k) && !extraSaved->count(k)) d.push_back(k + ": " + v + " no arquivo, agora não se aplica");
    return d;
}

void CheckParams(const std::string &parFile, const std::vector<std::string> &diffs, bool force,
                 std::vector<std::string> &forced) {
    if (diffs.empty()) return;
    std::ostringstream m;
    m << "parâmetros diferentes dos gravados em " << parFile << ":\n";
    for (auto &d : diffs) m << "    " << d << "\n";
    if (!force)
        throw TExitError(m.str() + "  as novas amostras não seriam do mesmo experimento: use outro saida= ou, "
                                   "conscientemente, forcar=1",
                         3);
    std::cout << "[MC] AVISO (forcar=1): " << m.str();
    forced.insert(forced.end(), diffs.begin(), diffs.end());
}

// ----- campos realizados das amostras (<csv>.campos, campos=1)
//
// Binário nativo (o mesmo de TPZBFileStream), little-endian no x86:
//   cabeçalho: char[8] "SSR-CMP\n" | int versão (1) | int64 np | int64 ponto[np] (índice de memória) |
//              double x[np] | double y[np] | int64 ne | int64 elemento[ne] | double xc[ne] | double yc[ne] |
//              uint64 FNV-1a do cabeçalho
//   registro:  char[8] "SSR-AMO\n" | int64 amostra | int64 nc | double c[nc] (kPa) | int64 nφ | double φ[nφ]
//              (graus) | int64 nk | double kv[nk] (nk = 0 se kv não é aleatório) | uint64 FNV-1a do registro
// O cabeçalho é gravado por TPZBFileStream (tmp + rename); cada registro é acrescentado com uma escrita +
// fdatasync. Um registro final truncado (SIGKILL) é descartado na retomada. Leitores: comando "campos"
// (ReadFields) e scripts/le_campos.py.
const char kFieldsMagic[8] = {'S', 'S', 'R', '-', 'C', 'M', 'P', '\n'};
const char kRecordMagic[8] = {'S', 'S', 'R', '-', 'A', 'M', 'O', '\n'};
const int kFieldsVersion = 1;

struct TFieldsHeader {
    std::vector<int64_t> point;  ///< índice de memória dos pontos de integração (c, φ)
    std::vector<REAL> xp, yp;    ///< coordenadas dos pontos
    std::vector<int64_t> gel;    ///< elementos geométricos (kv)
    std::vector<REAL> xe, ye;    ///< centróides
    bool operator==(const TFieldsHeader &o) const {
        return point == o.point && xp == o.xp && yp == o.yp && gel == o.gel && xe == o.xe && ye == o.ye;
    }
    uint64_t Hash() const {
        const int64_t np = point.size(), ne = gel.size();
        uint64_t h = fileio::HashWords(&np, 1, fileio::Fnv1a(kFieldsMagic, 8));
        h = fileio::HashWords(point.data(), point.size(), h);
        h = fileio::HashWords(xp.data(), xp.size(), h);
        h = fileio::HashWords(yp.data(), yp.size(), h);
        h = fileio::HashWords(&ne, 1, h);
        h = fileio::HashWords(gel.data(), gel.size(), h);
        h = fileio::HashWords(xe.data(), xe.size(), h);
        return fileio::HashWords(ye.data(), ye.size(), h);
    }
};

struct TFieldsRecord {
    int64_t sample = -1;
    std::vector<REAL> c, phi, kv;
    uint64_t Hash() const {
        const int64_t n[3] = {(int64_t)c.size(), (int64_t)phi.size(), (int64_t)kv.size()};
        uint64_t h = fileio::HashWords(&sample, 1, fileio::Fnv1a(kRecordMagic, 8));
        h = fileio::HashWords(n, 3, h);
        h = fileio::HashWords(c.data(), c.size(), h);
        h = fileio::HashWords(phi.data(), phi.size(), h);
        return fileio::HashWords(kv.data(), kv.size(), h);
    }
};

/// Lê um arquivo de campos (TPZBFileStream) e chama f para cada registro completo e íntegro, parando no primeiro
/// truncado ou corrompido (msg diz o motivo). Devolve o número de bytes válidos (0: cabeçalho inválido).
uint64_t ReadFields(const std::string &file, TFieldsHeader &h, const std::function<void(const TFieldsRecord &)> &f,
                    std::string &msg) {
    uint64_t fsize = 0, pos = 0;
    msg.clear();
    if (!fileio::FileSize(file, fsize)) {
        msg = "inexistente";
        return 0;
    }
    TPZBFileStream in;
    in.OpenRead(file);
    // TPZBFileStream não acusa leitura além do fim: cada leitura é conferida antes com o tamanho do arquivo
    auto have = [&](uint64_t n) { return pos + n <= fsize; };
    auto readI = [&](int64_t &v) {
        if (!have(8)) return false;
        in.Read(&v, 1);
        pos += 8;
        return true;
    };
    auto readVec = [&](auto &v, int64_t n) {
        using V = typename std::remove_reference_t<decltype(v)>::value_type;
        if (n < 0 || !have((uint64_t)n * sizeof(V))) return false;
        v.resize(n);
        for (int64_t i = 0; i < n; i += (1 << 24)) in.Read(v.data() + i, (int)std::min<int64_t>(n - i, 1 << 24));
        pos += (uint64_t)n * sizeof(V);
        return true;
    };
    uint64_t magic = 0;
    int version = 0;
    if (!have(12)) {
        msg = "cabeçalho truncado";
        return 0;
    }
    in.Read(&magic, 1);
    in.Read(&version, 1);
    pos += 12;
    if (magic != fileio::Word8(kFieldsMagic) || version != kFieldsVersion) {
        msg = "não é um arquivo de campos (versão " + std::to_string(kFieldsVersion) + ")";
        return 0;
    }
    int64_t np = -1, ne = -1;
    uint64_t hs = 0;
    if (!readI(np) || !readVec(h.point, np) || !readVec(h.xp, np) || !readVec(h.yp, np) || !readI(ne) ||
        !readVec(h.gel, ne) || !readVec(h.xe, ne) || !readVec(h.ye, ne) || !have(8)) {
        msg = "cabeçalho truncado";
        return 0;
    }
    in.Read(&hs, 1);
    pos += 8;
    if (hs != h.Hash()) {
        msg = "cabeçalho corrompido";
        return 0;
    }
    uint64_t valid = pos;
    while (pos < fsize) {
        TFieldsRecord r;
        uint64_t rm = 0;
        int64_t n = -1;
        uint64_t sum = 0;
        if (!have(16)) {
            msg = "registro final truncado";
            break;
        }
        in.Read(&rm, 1);
        in.Read(&r.sample, 1);
        pos += 16;
        if (rm != fileio::Word8(kRecordMagic)) {
            msg = "registro corrompido";
            break;
        }
        if (!readI(n) || n != np || !readVec(r.c, n) || !readI(n) || n != np || !readVec(r.phi, n) || !readI(n) ||
            (n != 0 && n != ne) || !readVec(r.kv, n) || !have(8)) {
            msg = "registro final truncado ou corrompido";
            break;
        }
        in.Read(&sum, 1);
        pos += 8;
        if (sum != r.Hash()) {
            msg = "registro corrompido (soma de verificação)";
            break;
        }
        valid = pos;
        if (f) f(r);
    }
    return valid;
}

/// Acrescenta os campos das amostras a <csv>.campos (cria o cabeçalho; na retomada confere a malha, descarta um
/// registro final incompleto e não repete amostras já gravadas)
class TFieldsWriter {
public:
    TFieldsWriter(const std::string &file, const TFieldsHeader &h) : fFile(file) {
        if (fileio::Exists(file)) {
            TFieldsHeader old;
            std::string msg;
            uint64_t fsize = 0;
            fileio::FileSize(file, fsize);
            const uint64_t valid =
                ReadFields(file, old, [this](const TFieldsRecord &r) { fIds.insert(r.sample); }, msg);
            if (valid == 0)
                throw TExitError("arquivo de campos " + file + " inválido (" + msg + "); renomeie-o ou apague-o", 2);
            if (!(old == h))
                throw TExitError("o arquivo de campos " + file + " é de outra malha (pontos ou elementos diferentes)", 3);
            if (valid < fsize) {
                std::cout << "[MC] " << file << ": " << msg << "; " << fsize - valid << " bytes finais descartados\n";
                if (::truncate(file.c_str(), (off_t)valid) != 0) throw TExitError("não foi possível truncar " + file, 2);
            }
            std::cout << "[MC] campos: " << fIds.size() << " amostra(s) já gravada(s) em " << file << "\n";
        } else {
            WriteHeader(h);
        }
        fOut = std::make_unique<fileio::TAppendFile>(file);
    }
    bool Has(int64_t s) const { return fIds.count(s) > 0; }
    void Append(const TFieldsRecord &r) {
        std::string b;
        auto put = [&b](const void *p, size_t n) { b.append(static_cast<const char *>(p), n); };
        const int64_t nc = r.c.size(), nphi = r.phi.size(), nk = r.kv.size();
        const uint64_t sum = r.Hash();
        put(kRecordMagic, 8);
        put(&r.sample, 8);
        put(&nc, 8);
        put(r.c.data(), 8 * nc);
        put(&nphi, 8);
        put(r.phi.data(), 8 * nphi);
        put(&nk, 8);
        put(r.kv.data(), 8 * nk);
        put(&sum, 8);
        fOut->Append(b);
        fIds.insert(r.sample);
    }

private:
    void WriteHeader(const TFieldsHeader &h) {
        const std::string tmp = fileio::TmpName(fFile);
        {
            TPZBFileStream out;
            out.OpenWrite(tmp);
            if (!out.AmIOpenForWrite()) throw TExitError("não foi possível criar " + tmp, 2);
            const int version = kFieldsVersion;
            const int64_t np = h.point.size(), ne = h.gel.size();
            const uint64_t hs = h.Hash();
            const uint64_t magic = fileio::Word8(kFieldsMagic);
            out.Write(&magic, 1);
            out.Write(&version, 1);
            out.Write(&np, 1);
            out.Write(h.point.data(), (int)np);
            out.Write(h.xp.data(), (int)np);
            out.Write(h.yp.data(), (int)np);
            out.Write(&ne, 1);
            out.Write(h.gel.data(), (int)ne);
            out.Write(h.xe.data(), (int)ne);
            out.Write(h.ye.data(), (int)ne);
            out.Write(&hs, 1);
            out.CloseWrite();
        }
        fileio::Commit(tmp, fFile);
    }
    std::string fFile;
    std::set<int64_t> fIds;
    std::unique_ptr<fileio::TAppendFile> fOut;
};

/// Comando "campos": lê <csv>.campos; sem amostra=, uma linha "amostra,c_medio,phi_medio,kv_medio" por registro
/// (médias como no CSV do mc); com amostra=k, os campos dessa amostra (x,y,c,phi por ponto e xc,yc,kv por
/// elemento) em saida= (ou na tela)
int RunFields(const TArgs &args) {
    const std::string file = args.GetS("arquivo", "");
    if (file.empty()) throw TExitError("campos: informe arquivo=<csv>.campos", 1);
    const bool one = args.kv.count("amostra") > 0;
    const int64_t want = args.GetI("amostra", -1);
    TFieldsHeader h;
    std::string msg;
    int64_t nrec = 0, found = 0;
    std::ostringstream dump;
    dump << std::setprecision(17);
    std::cout << std::setprecision(8);
    auto onRecord = [&](const TFieldsRecord &r) {
        nrec++;
        if (one) {
            if (r.sample != want || found++) return;
            dump << "x,y,c,phi\n";
            for (size_t i = 0; i < r.c.size(); i++) dump << h.xp[i] << "," << h.yp[i] << "," << r.c[i] << "," << r.phi[i] << "\n";
            dump << "xc,yc,kv\n";
            for (size_t e = 0; e < r.kv.size(); e++) dump << h.xe[e] << "," << h.ye[e] << "," << r.kv[e] << "\n";
            return;
        }
        if (nrec == 1) std::cout << "amostra,c_medio,phi_medio,kv_medio\n";
        REAL cm = 0., pm = 0., km = r.kv.empty() ? 1. : 0.;
        for (REAL v : r.c) cm += v;
        for (REAL v : r.phi) pm += v;
        for (REAL v : r.kv) km += v / r.kv.size();
        cm /= r.c.size();
        pm /= r.phi.size();
        std::cout << r.sample << "," << cm << "," << pm << "," << km << "\n";
    };
    uint64_t fsize = 0;
    fileio::FileSize(file, fsize);
    const uint64_t valid = ReadFields(file, h, onRecord, msg);
    if (valid == 0) throw TExitError("campos: " + file + ": " + msg, 2);
    std::cerr << "[campos] " << file << ": " << h.point.size() << " pontos de integração, " << h.gel.size()
              << " elementos, " << nrec << " registro(s)" << (valid < fsize ? " (" + msg + ": final ignorado)" : "")
              << "\n";
    if (one) {
        if (!found) throw TExitError("campos: amostra " + std::to_string(want) + " não está em " + file, 2);
        const std::string out = args.GetS("saida", "");
        if (out.empty()) std::cout << dump.str();
        else fileio::WriteAtomic(out, dump.str());
    }
    return 0;
}

// ----- malha adaptada (malha=arquivo): elementos divididos em cada nível, refeitos sem as análises
struct TAdaptLevel {
    int64_t neq = 0, nmarked = 0, nel = 0, nnod = 0;  ///< equações e marcados do nível; elementos e nós depois
    REAL factor = 0.;                                 ///< Γ (ou FS) do problema médio no nível
    uint64_t sig = 0;                                 ///< assinatura da malha depois do nível
    std::vector<int64_t> divided;                     ///< elementos divididos (marcados + camadas)
};

/// Parâmetros de que a adaptação depende (problema médio com Mohr-Coulomb, ver AdaptedMesh); o número de níveis
/// não entra: os primeiros níveis de uma adaptação mais profunda são os mesmos
std::string AdaptKey(const TCase &cs, const TArgs &args) {
    const bool fs = args.GetS("medida", "gamma") == "fs";
    std::ostringstream k;
    k << "versao=1 H=" << Fmt(cs.geo.H) << " beta=" << Fmt(cs.geo.betaDeg) << " Lc=" << Fmt(cs.geo.Lc)
      << " Lt=" << Fmt(cs.geo.Lt) << " Hb=" << Fmt(cs.geo.Hb) << " h=" << Fmt(cs.geo.h) << " tri=" << cs.geo.triangles
      << " c=" << Fmt(cs.soil.c) << " phi=" << Fmt(cs.soil.phiDeg) << " gam=" << Fmt(cs.soil.gamma)
      << " gamw=" << Fmt(cs.soil.gammaW) << " flutuante=" << cs.soil.buoyant << " E=" << Fmt(cs.soil.E)
      << " nu=" << Fmt(cs.soil.nu) << " percolacao=" << cs.seepage;
    if (cs.seepage) k << " hw=" << Fmt(cs.hw) << " alpha=" << Fmt(cs.alpha);
    k << " medida=" << (fs ? "fs" : "gamma");
    if (fs) k << " F0=" << Fmt(args.Get("F0", 0.5));
    k << " p=" << args.GetI("p", 2) << " reltol=" << Fmt(args.Get("reltol", 5.e-3)) << " maxit=" << args.GetI("maxit", 20)
      << " estagnacao=" << (args.GetI("estagnacao", 1) != 0) << " frac=" << Fmt(args.Get("frac", 0.02))
      << " camadas=" << args.GetI("camadas", 1) << " incremento=" << (args.GetI("incremento", 1) != 0);
    return k.str();
}

/// malha=0: sem arquivo; malha=<arquivo>; padrão malha_<caso>_h<h>_<FNV-1a dos parâmetros>.malha
std::string AdaptFile(const TCase &cs, const TArgs &args, const std::string &key) {
    const std::string f = args.GetS("malha", "");
    if (f == "0") return "";
    if (!f.empty() && f != "1") return f;
    return "malha_" + cs.name + "_h" + Fmt(cs.geo.h) + "_" +
           fileio::Hex(fileio::Fnv1a(key.data(), key.size())).substr(0, 8) + ".malha";
}

void WriteAdaptFile(const std::string &file, const std::string &key, int64_t nel0, int64_t nnod0, uint64_t sig0,
                    const std::vector<TAdaptLevel> &lv) {
    std::ostringstream s;
    s << "# SlopeSeepageRandom: malha adaptada = elementos divididos em cada nível (AdaptedMesh em main.cpp)\n"
      << "formato 1\n"
      << "chave " << key << "\n"
      << "base " << nel0 << " " << nnod0 << " " << fileio::Hex(sig0) << "\n"
      << "niveis " << lv.size() << "\n"
      << std::setprecision(17);
    for (size_t l = 0; l < lv.size(); l++) {
        const TAdaptLevel &L = lv[l];
        s << "nivel " << l << " " << L.neq << " " << L.factor << " " << L.nmarked << " " << L.divided.size() << " "
          << L.nel << " " << L.nnod << " " << fileio::Hex(L.sig) << "\n";
        for (size_t i = 0; i < L.divided.size(); i++)
            s << L.divided[i] << ((i + 1) % 20 == 0 || i + 1 == L.divided.size() ? "\n" : " ");
    }
    s << "fim\n";
    fileio::WriteAtomic(file, s.str());
}

bool ReadAdaptFile(const std::string &file, const std::string &key, int64_t nel0, int64_t nnod0, uint64_t sig0,
                   std::vector<TAdaptLevel> &lv, std::string &why) {
    std::string content;
    if (!fileio::ReadAll(file, content)) {
        why = "ilegível";
        return false;
    }
    std::istringstream in(content);
    std::string line, tok, hx;
    std::getline(in, line);
    if (!std::getline(in, line) || line != "formato 1") {
        why = "de formato desconhecido";
        return false;
    }
    if (!std::getline(in, line) || line != "chave " + key) {
        why = "é de outro caso (parâmetros da adaptação diferentes)";
        return false;
    }
    auto hex = [](const std::string &s, uint64_t &v) {
        char *end = nullptr;
        v = std::strtoull(s.c_str(), &end, 16);
        return s.size() == 16 && *end == '\0';
    };
    int64_t a = -1, b = -1, nl = -1;
    uint64_t sig = 0;
    if (!(in >> tok >> a >> b >> hx) || tok != "base" || !hex(hx, sig)) {
        why = "corrompido";
        return false;
    }
    if (a != nel0 || b != nnod0 || sig != sig0) {
        why = "parte de outra malha inicial";
        return false;
    }
    if (!(in >> tok >> nl) || tok != "niveis" || nl < 0 || nl > 1000) {
        why = "corrompido";
        return false;
    }
    lv.assign(nl, TAdaptLevel());
    for (int64_t l = 0; l < nl; l++) {
        TAdaptLevel &L = lv[l];
        int64_t idx = -1, ndiv = -1;
        if (!(in >> tok >> idx >> L.neq >> L.factor >> L.nmarked >> ndiv >> L.nel >> L.nnod >> hx) || tok != "nivel" ||
            idx != l || ndiv < 0 || ndiv > 100000000 || !hex(hx, L.sig)) {
            why = "corrompido (nível " + std::to_string(l) + ")";
            return false;
        }
        L.divided.resize(ndiv);
        for (int64_t &g : L.divided)
            if (!(in >> g) || g < 0) {
                why = "truncado ou corrompido (nível " + std::to_string(l) + ")";
                return false;
            }
    }
    if (!(in >> tok) || tok != "fim") {
        why = "truncado";
        return false;
    }
    return true;
}

/// Malha do caso, adaptada (ver a definição)
std::unique_ptr<TPZGeoMesh> AdaptedMesh(const TCase &cs, const TArgs &args);

/// Tensão efetiva geostática: análise elástica (Mohr-Coulomb com c muito alta) com a força de corpo (γ') e sem
/// percolação, na mesma malha e ordem; devolve σ'0 nos pontos dados (que devem coincidir com os da análise)
template <class TPoints>
std::vector<TPZTensor<REAL>> GeostaticStress(const TCase &cs, TPZGeoMesh *gmesh, const TSolverOptions &opt,
                                             const TPoints &target) {
    TSoil elastic = cs.soil;
    elastic.c = 1.e6;
    TSlopeFEM<TMohrCoulomb> fe(gmesh, cs.geo, elastic, opt);
    fe.ResetState();
    int its = 0;
    if (!fe.Solve(1., 0., 1., its)) throw std::runtime_error("análise geostática não convergiu");
    std::vector<TPZTensor<REAL>> s;
    fe.Stresses(s);
    const auto &pa = fe.Points();
    if (pa.size() != target.size()) throw std::runtime_error("pontos de integração diferentes (geostática)");
    for (size_t i = 0; i < pa.size(); i++)
        if (target[i].gel >= 0 && std::hypot(pa[i].x[0] - target[i].x[0], pa[i].x[1] - target[i].x[1]) > 1.e-9)
            throw std::runtime_error("pontos de integração diferentes (geostática)");
    return s;
}

/// Tensão inicial do Cam-Clay (Mohr-Coulomb: nada a fazer)
template <class T>
void GeostaticStress(const TCase &cs, TPZGeoMesh *gmesh, const TSolverOptions &opt, TSlopeFEM<T> &target) {
    if constexpr (std::is_same_v<T, TPZModifiedCamClay>) {
        target.SetInitialStress(GeostaticStress(cs, gmesh, opt, target.Points()));
    } else {
        (void)cs;
        (void)gmesh;
        (void)opt;
        (void)target;
    }
}

/// Estatísticas de um campo lognormal (média, CoV; CoV = 0: determinístico)
struct TFieldSpec {
    REAL mean = 1., cov = 0.;
};

/// Parâmetros que definem as amostras do Monte Carlo (valores efetivos; vão para <csv>.param). Não entram inicio,
/// n, amostras, saida, vtk, campos, malha, klcache, verbose e forcar, que não mudam o resultado de uma amostra.
TParamList MCInputParams(const TCase &cs, const TArgs &args, const std::string &model) {
    TParamList p;
    auto add = [&p](const std::string &k, const std::string &v) { p.emplace_back(k, v); };
    const std::string medida = args.GetS("medida", "gamma");
    const REAL Ly = args.Get("Ly", 2.);
    const int adapt = args.GetI("adapt", 0);
    add("versao_param", "1");
    add("comando", "mc");
    add("caso", cs.name);
    add("modelo", model);
    add("medida", medida);
    add("seed", std::to_string((uint64_t)args.GetI("seed", 2025)));
    add("covc", Fmt(args.Get("covc", 0.3)));
    add("covphi", Fmt(args.Get("covphi", 0.1)));
    if (cs.seepage) add("covk", Fmt(args.Get("covk", 0.6)));
    add("Lx", Fmt(args.Get("Lx", 20.)));
    add("Ly", Fmt(Ly));
    add("hkl", Fmt(args.Get("hkl", std::min(Ly / 2., 1.))));
    add("M", std::to_string(args.GetI("M", -1)));
    add("epsM", Fmt(args.Get("epsM", -1.)));
    add("normvar", std::to_string(args.GetI("normvar", 1) != 0));
    add("H", Fmt(cs.geo.H));
    add("beta", Fmt(cs.geo.betaDeg));
    add("Lc", Fmt(cs.geo.Lc));
    add("Lt", Fmt(cs.geo.Lt));
    add("Hb", Fmt(cs.geo.Hb));
    add("h", Fmt(cs.geo.h));
    add("tri", std::to_string(cs.geo.triangles));
    add("c", Fmt(cs.soil.c));
    add("phi", Fmt(cs.soil.phiDeg));
    add("gam", Fmt(cs.soil.gamma));
    add("gamw", Fmt(cs.soil.gammaW));
    add("E", Fmt(cs.soil.E));
    add("nu", Fmt(cs.soil.nu));
    if (cs.seepage) {
        add("hw", Fmt(cs.hw));
        add("alpha", Fmt(cs.alpha));
    }
    add("p", std::to_string(args.GetI("p", 2)));
    add("reltol", Fmt(args.Get("reltol", 5.e-3)));
    add("maxit", std::to_string(args.GetI("maxit", 20)));
    add("estagnacao", std::to_string(args.GetI("estagnacao", 1) != 0));
    if (medida == "fs") add("F0", Fmt(args.Get("F0", 0.5)));
    add("adapt", std::to_string(adapt));
    if (adapt > 0) {
        add("frac", Fmt(args.Get("frac", 0.02)));
        add("camadas", std::to_string(args.GetI("camadas", 1)));
        add("incremento", std::to_string(args.GetI("incremento", 1) != 0));
    }
    if (model == "mcc") {
        add("lambda", Fmt(cs.soil.lambda));
        add("kappa", Fmt(cs.soil.kappa));
        add("v0", Fmt(cs.soil.v0));
        add("OCR", Fmt(cs.soil.OCR));
        add("mapeamento", args.GetS("mapeamento", "deformacao_plana"));
    }
    return p;
}

/// Estatísticas de todas as linhas do CSV (não só das desta execução), na tela e em <csv>.resumo (atômico)
void Summarize(const std::string &csv, const std::vector<TMCRow> &rows) {
    struct TStat {
        int64_t n = 0;
        REAL mean = 0., sd = 0., cov = 0., pf = 0., covpf = 0., min = 0., max = 0.;
    };
    auto stat = [](const std::vector<REAL> &v) {
        TStat s;
        s.n = (int64_t)v.size();
        if (!s.n) return s;
        s.min = s.max = v[0];
        int64_t nf = 0;
        for (REAL x : v) {
            s.mean += x;
            s.min = std::min(s.min, x);
            s.max = std::max(s.max, x);
            if (x < 1.) nf++;
        }
        s.mean /= s.n;
        REAL ss = 0.;
        for (REAL x : v) ss += (x - s.mean) * (x - s.mean);
        s.sd = s.n > 1 ? std::sqrt(ss / (s.n - 1)) : 0.;
        s.cov = s.mean != 0. ? s.sd / s.mean : 0.;
        s.pf = REAL(nf) / s.n;
        s.covpf = s.pf > 0. ? std::sqrt((1. - s.pf) / (s.n * s.pf)) : INFINITY;
        return s;
    };
    std::vector<REAL> lo, mid;
    std::map<std::string, int64_t> status;
    int64_t nonfinite = 0, smin = -1, smax = -1;
    REAL ttot = 0.;
    for (const TMCRow &r : rows) {
        status[r.status]++;
        if (std::isfinite(r.time)) ttot += r.time;
        smin = (smin < 0) ? r.sample : std::min(smin, r.sample);
        smax = std::max(smax, r.sample);
        if (!std::isfinite(r.factor)) {
            nonfinite++;
            continue;
        }
        lo.push_back(r.factor);
        mid.push_back(r.status == "ok" && r.upper > r.factor ? 0.5 * (r.factor + r.upper) : r.factor);
    }
    const TStat a = stat(lo), b = stat(mid);
    std::ostringstream o;
    o << std::setprecision(8);
    o << "# SlopeSeepageRandom mc: resumo de TODAS as linhas de " << csv << " (" << Now() << ")\n"
      << "# fator = último valor convergido (coluna fator); ponto_medio = (fator + limite_superior)/2 se status=ok\n"
      << "# (convenção de analisa_mc.py); desvio com n-1; Pf = P(fator < 1); cov_Pf = sqrt((1 - Pf)/(N Pf))\n"
      << "arquivo=" << csv << "\n"
      << "N=" << rows.size() << "\n"
      << "amostra_min=" << smin << "\n"
      << "amostra_max=" << smax << "\n"
      << "faltando_entre_min_e_max=" << (rows.empty() ? 0 : smax - smin + 1 - (int64_t)rows.size()) << "\n";
    for (auto &[st, cnt] : status) o << "status_" << st << "=" << cnt << "\n";
    o << "nao_finitos=" << nonfinite << "\n"
      << "media=" << a.mean << "\ndesvio=" << a.sd << "\ncov=" << a.cov << "\nminimo=" << a.min << "\nmaximo=" << a.max
      << "\nPf=" << a.pf << "\ncov_Pf=" << a.covpf << "\n"
      << "media_ponto_medio=" << b.mean << "\ndesvio_ponto_medio=" << b.sd << "\nPf_ponto_medio=" << b.pf
      << "\ncov_Pf_ponto_medio=" << b.covpf << "\n"
      << "tempo_total_s=" << ttot << "\n";
    fileio::WriteAtomic(csv + ".resumo", o.str());
    std::cout << "[MC] resumo de " << csv << " (todas as " << rows.size() << " linhas):";
    if (a.n) {
        std::cout << " média " << a.mean << ", desvio " << a.sd << " (CoV " << a.cov << "), Pf = P(fator < 1) = " << a.pf
                  << " (CoV(Pf) = " << a.covpf << "); ponto médio do intervalo de colapso: média " << b.mean
                  << ", desvio " << b.sd << ", Pf = " << b.pf << "; status:";
        for (auto &[st, cnt] : status) std::cout << " " << st << " " << cnt;
        if (nonfinite) std::cout << " (" << nonfinite << " fatores não finitos)";
    } else {
        std::cout << " sem amostras";
    }
    std::cout << " -> " << csv << ".resumo\n";
}

/// Monte Carlo. Retomada nativa: o CSV (saida=) é lido e reparado (linhas incompletas, inválidas ou repetidas são
/// retiradas, regravação atômica), as amostras de [inicio, inicio+n) já gravadas são puladas e cada nova linha é
/// gravada inteira com fdatasync; os parâmetros ficam em <csv>.param e são conferidos na retomada; as
/// estatísticas de todo o CSV vão para <csv>.resumo. amostras=lista recalcula amostras e confere com o CSV;
/// campos=1 grava os campos realizados em <csv>.campos. Devolve o código de saída.
template <class T>
int RunMonteCarlo(const TCase &cs, const TArgs &args, const std::string &modelName) {
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.verbose = args.GetI("verbose", 0);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    const std::string medida = args.GetS("medida", "gamma");  // gamma (fator de carga) ou fs (redução)
    const int64_t n = std::max(args.GetI("n", 100), 0), first = args.GetI("inicio", 0);
    const uint64_t seed = (uint64_t)args.GetI("seed", 2025);
    TFieldSpec fc{cs.soil.c, args.Get("covc", 0.3)}, fphi{cs.soil.phiDeg, args.Get("covphi", 0.1)},
        fk{1., args.Get("covk", cs.seepage ? 0.6 : 0.)};
    const REAL Lx = args.Get("Lx", 20.), Ly = args.Get("Ly", 2.);
    const std::string csv = args.GetS("saida", "mc_" + cs.name + "_" + modelName + ".csv");
    const bool force = args.GetI("forcar", 0) != 0, replay = args.kv.count("amostras") > 0,
               withFields = args.GetI("campos", 0) != 0;
    auto t0 = Clock::now();

    // 1) trava, CSV existente e parâmetros gravados (antes de qualquer cálculo)
    fileio::TFileLock lock;
    if (!lock.Acquire(csv + ".lock"))
        throw TExitError("outro processo está gravando " + csv + " (trava em " + csv + ".lock)", 4);
    TCsvContent content = ReadCsv(csv);
    const std::string parFile = csv + ".param";
    const TParamList inputs = MCInputParams(cs, args, modelName);
    const std::set<std::string> derivedKeys = {"equacoes", "pontos_integracao", "assinatura_malha", "equacoes_KL",
                                               "modos_KL"};
    std::map<std::string, std::string> saved;
    std::string savedRaw;
    const bool hasPar = ReadParamFile(parFile, saved, savedRaw);
    const bool checkPar = hasPar && !content.rows.empty();
    std::vector<std::string> forced;
    if (checkPar) CheckParams(parFile, DiffParams(saved, inputs, true, &derivedKeys), force, forced);
    else if (!content.rows.empty())
        std::cout << "[MC] AVISO: " << csv << " tem " << content.rows.size() << " amostras e não tem " << parFile
                  << " (CSV de versão anterior?): parâmetros não conferidos; os atuais serão gravados\n";
    if (content.NeedsRewrite()) RewriteCsv(csv, content);
    std::vector<TMCRow> rows = std::move(content.rows);
    std::map<int64_t, size_t> rowOf;
    int64_t nfAll = 0;
    for (size_t i = 0; i < rows.size(); i++) {
        rowOf[rows[i].sample] = i;
        if (rows[i].factor < 1.) nfAll++;
    }

    // 2) amostras a calcular
    std::vector<int64_t> todo;
    if (replay) {
        todo = ParseSampleList(args.GetS("amostras", ""));
        int64_t present = 0;
        for (int64_t s : todo) present += rowOf.count(s);
        std::cout << "[MC] amostras=" << args.GetS("amostras", "") << ": " << todo.size() << " amostra(s), " << present
                  << " já no CSV (recalculadas e conferidas, sem nova linha), " << todo.size() - present
                  << " nova(s)\n";
    } else {
        int64_t skipped = 0;
        for (int64_t s = first; s < first + n; s++) {
            if (rowOf.count(s)) skipped++;
            else todo.push_back(s);
        }
        if (n > 0)
            std::cout << "[MC] amostras " << first << ".." << first + n - 1 << ": " << skipped << " já em " << csv
                      << " (puladas), " << todo.size() << " a calcular\n";
    }
    if (todo.empty() && (replay || n > 0)) {
        Summarize(csv, rows);
        return 0;
    }

    std::unique_ptr<TPZGeoMesh> gmesh = AdaptedMesh(cs, args);
    TSlopeFEM<T> fem(gmesh.get(), cs.geo, cs.soil, opt);
    GeostaticStress(cs, gmesh.get(), opt, fem);
    std::cout << "[MC] malha mecânica: " << fem.NEquations() << " equações, " << fem.NPoints()
              << " pontos de integração\n";

    // campo aleatório: malha KL própria (quadriláteros de 9 nós, h_KL)
    TSlopeGeometry geoKL = cs.geo;
    geoKL.h = args.Get("hkl", std::min(Ly / 2., 1.));
    geoKL.triangles = false;
    TPZKLRandomField::TOptions klopt;
    klopt.Lx = Lx;
    klopt.Ly = Ly;
    klopt.porder = 2;
    klopt.nModes = args.GetI("M", -1);
    klopt.targetVarianceError = args.Get("epsM", -1.);
    klopt.normalizeVariance = args.GetI("normvar", 1) != 0;
    {
        std::ostringstream f;
        f << "kl_H" << cs.geo.H << "_b" << std::lround(cs.geo.betaDeg * 100) << "_Lx" << Lx << "_Ly" << Ly << "_h"
          << geoKL.h << "_v3.bin";
        klopt.cacheFile = args.GetS("klcache", f.str());
    }
    auto tkl = Clock::now();
    TPZKLRandomField kl(geoKL.CreateGeoMesh(), klopt);
    kl.Compute();
    std::cout << "[KL] " << kl.NEquations() << " equações, M = " << kl.NModes() << ", ε_M = "
              << kl.VarianceError(kl.NModes()) << " (" << Seconds(tkl) << " s)\n";

    // pontos alvo: pontos de integração da malha mecânica (c, φ) e centróides dos elementos (kv)
    std::vector<TPZManVector<REAL, 3>> ipts;
    for (auto &p : fem.Points()) ipts.push_back(p.gel >= 0 ? p.x : TPZManVector<REAL, 3>(3, 0.));
    const int setIP = kl.AddTargetSet(ipts);
    std::vector<int64_t> elems;
    std::vector<TPZManVector<REAL, 3>> cpts;
    for (TPZGeoEl *gel : gmesh->ElementVec()) {
        if (!gel || gel->Dimension() != 2 || gel->HasSubElement()) continue;
        TPZManVector<REAL, 3> qc(2, 0.), xc(3, 0.);
        gel->CenterPoint(gel->NSides() - 1, qc);
        gel->X(qc, xc);
        elems.push_back(gel->Index());
        cpts.push_back(xc);
    }
    const int setEl = kl.AddTargetSet(cpts);
    {
        REAL vmin = 1., vmean = 0.;
        for (REAL v : kl.PointVariance(setIP)) {
            vmin = std::min(vmin, v);
            vmean += v / kl.PointVariance(setIP).size();
        }
        std::cout << "[KL] variância truncada nos pontos de integração: média " << vmean << ", mínima " << vmin
                  << (klopt.normalizeVariance ? " (compensada)" : " (não compensada)") << "\n";
    }

    std::unique_ptr<TSeepageProblem> seep;
    if (cs.seepage) {
        TSeepageProblem::TParams sp;
        sp.hw = cs.hw;
        sp.alpha = cs.alpha;
        sp.gammaW = cs.soil.gammaW;
        seep = std::make_unique<TSeepageProblem>(gmesh.get(), cs.geo, sp);
        if (fk.cov == 0.) {
            seep->Solve();
            TransferSeepage(*seep, fem);
        }
    }
    const bool randomKv = seep && fk.cov > 0.;

    // 3) parâmetros derivados (malha e KL: mudam se o código ou a biblioteca mudarem a discretização) e .param
    const TParamList derived = {{"equacoes", std::to_string(fem.NEquations())},
                                {"pontos_integracao", std::to_string(fem.NPoints())},
                                {"assinatura_malha", fileio::Hex(TSlopeGeometry::Signature(*gmesh))},
                                {"equacoes_KL", std::to_string(kl.NEquations())},
                                {"modos_KL", std::to_string(kl.NModes())}};
    if (checkPar) CheckParams(parFile, DiffParams(saved, derived, false, nullptr), force, forced);
    if (!checkPar) {
        TParamList all = inputs;
        all.insert(all.end(), derived.begin(), derived.end());
        fileio::WriteAtomic(parFile, ParamText(csv, all));
    } else if (!forced.empty()) {
        std::string note = savedRaw;
        if (!note.empty() && note.back() != '\n') note += "\n";
        for (auto &d : forced) note += "# " + Now() + " forcar=1, amostras seguintes com " + d + "\n";
        fileio::WriteAtomic(parFile, note);
    }

    // campos realizados (opcional)
    std::unique_ptr<TFieldsWriter> fieldsOut;
    if (withFields) {
        TFieldsHeader fh;
        for (int64_t i = 0; i < fem.NPoints(); i++) {
            if (fem.Points()[i].gel < 0) continue;
            fh.point.push_back(i);
            fh.xp.push_back(fem.Points()[i].x[0]);
            fh.yp.push_back(fem.Points()[i].x[1]);
        }
        if (randomKv)
            for (size_t e = 0; e < elems.size(); e++) {
                fh.gel.push_back(elems[e]);
                fh.xe.push_back(cpts[e][0]);
                fh.ye.push_back(cpts[e][1]);
            }
        fieldsOut = std::make_unique<TFieldsWriter>(csv + ".campos", fh);
    }

    fileio::TAppendFile out(csv);
    const int vtkEvery = args.GetI("vtk", 0);
    std::cout << "[MC] " << cs.name << " (" << modelName << "), medida = " << medida << ", "
              << (replay ? "amostras " + args.GetS("amostras", "")
                         : "amostras " + std::to_string(first) + ".." + std::to_string(first + n - 1))
              << ", CoV(c) = " << fc.cov << ", CoV(phi) = " << fphi.cov << ", CoV(kv) = " << fk.cov << ", Lx = " << Lx
              << ", Ly = " << Ly << " -> " << csv << "\n";
    std::vector<REAL> c(fem.NPoints()), phi(fem.NPoints());
    std::vector<REAL> kvEl(gmesh->NElements(), 1.);
    int64_t added = 0, same = 0, differ = 0;
    for (size_t k = 0; k < todo.size(); k++) {
        const int64_t s = todo[k];
        auto ts = Clock::now();
        std::vector<std::vector<REAL>> xi(3);
        for (int f = 0; f < 3; f++) kl.Xi(seed, s, f, xi[f]);
        std::vector<std::vector<std::vector<REAL>>> H;
        kl.Evaluate(xi, {setIP, setEl}, H);
        TFieldsRecord rec;
        rec.sample = s;
        REAL cm = 0., pm = 0., km = 0.;
        int64_t np = 0;
        for (int64_t i = 0; i < fem.NPoints(); i++) {
            if (fem.Points()[i].gel < 0) continue;
            c[i] = TPZKLRandomField::Lognormal(H[0][setIP][i], fc.mean, fc.cov);
            const REAL phideg = TPZKLRandomField::Lognormal(H[1][setIP][i], fphi.mean, fphi.cov);
            phi[i] = phideg * M_PI / 180.;
            cm += c[i];
            pm += phideg;
            np++;
            if (fieldsOut) {
                rec.c.push_back(c[i]);
                rec.phi.push_back(phideg);
            }
        }
        cm /= np;
        pm /= np;
        fem.SetStrength(c, phi);
        if (randomKv) {
            for (size_t e = 0; e < elems.size(); e++) {
                kvEl[elems[e]] = TPZKLRandomField::Lognormal(H[2][setEl][e], fk.mean, fk.cov);
                km += kvEl[elems[e]] / elems.size();
                if (fieldsOut) rec.kv.push_back(kvEl[elems[e]]);
            }
            seep->SetElementPermeability(kvEl);
            seep->Solve();
            TransferSeepage(*seep, fem);
        } else {
            km = 1.;
        }
        fem.ResetState();
        TFactorResult r = (medida == "fs") ? fem.StrengthReduction(args.Get("F0", 0.5)) : fem.LoadFactor(0.);
        const double dt = Seconds(ts);
        const TMCRow row = MakeRow(s, r, dt, cm, pm, km);
        // campos antes da linha do CSV: uma amostra no CSV tem os campos gravados (com campos=1)
        if (fieldsOut && !fieldsOut->Has(s)) fieldsOut->Append(rec);
        auto it = rowOf.find(s);
        if (it != rowOf.end()) {  // amostras=: confere com a linha gravada (todas as colunas menos tempo_s)
            const TMCRow &old = rows[it->second];
            bool eq = old.fields.size() == row.fields.size();
            for (size_t j = 0; eq && j < row.fields.size(); j++)
                if ((int)j != kCsvTimeColumn && old.fields[j] != row.fields[j]) eq = false;
            (eq ? same : differ)++;
            std::cout << "  amostra " << s << (eq ? ": reproduz a linha do CSV" : ": DIFERE da linha do CSV") << "\n";
            if (!eq) std::cout << "    CSV:   " << old.Line() << "\n    agora: " << row.Line() << "\n";
        } else {
            out.Append(row.Line() + "\n");
            rowOf[s] = rows.size();
            rows.push_back(row);
            added++;
            if (r.factor < 1.) nfAll++;
        }
        if (vtkEvery > 0 && (replay ? (int64_t)k : s - first) % vtkEvery == 0) {
            const std::string base = cs.name + "_" + modelName + "_amostra" + std::to_string(s);
            fem.DefineVTK(base + ".vtk");
            fem.WriteVTK(0);
            if (seep) {
                seep->DefineVTK(base + "_darcy.vtk");
                seep->WriteVTK(0);
            }
        }
        if (opt.verbose || k % 10 == 0)
            std::cout << "  amostra " << s << ": fator = " << r.factor << " (" << dt << " s), Pf parcial do CSV = "
                      << REAL(nfAll) / std::max<size_t>(rows.size(), 1) << "\n";
    }
    std::cout << "[MC] " << added << " amostra(s) nova(s) gravada(s) nesta execução";
    if (replay) std::cout << "; recalculadas do CSV: " << same << " reproduzem, " << differ << " diferem";
    std::cout << " (tempo " << Seconds(t0) << " s)\n";
    Summarize(csv, rows);
    return differ ? 5 : 0;
}

/// Malha geométrica do caso adaptada ao mecanismo de colapso do problema com propriedades médias (Mohr-Coulomb):
/// em cada um dos adapt=N níveis resolve Γ (ou FS), marca os elementos com ||Δε^p|| > frac max e "camadas" de
/// vizinhos e os divide. Os elementos divididos em cada nível vão para um arquivo (malha=<arquivo>; padrão
/// malha_<caso>_h<h>_<FNV-1a dos parâmetros da adaptação>.malha; malha=0 desliga), gravado atomicamente; numa
/// nova execução os níveis são refeitos por TSlopeGeometry::Refine a partir de CreateGeoMesh, sem as análises,
/// e conferidos nível a nível (número de elementos e de nós, assinatura exata). Se algo não confere (arquivo
/// truncado, corrompido ou de outro caso), a adaptação é recalculada e o arquivo regravado. Um arquivo com
/// menos níveis é completado; um com mais níveis serve aos primeiros.
std::unique_ptr<TPZGeoMesh> AdaptedMesh(const TCase &cs, const TArgs &args) {
    std::unique_ptr<TPZGeoMesh> gmesh(cs.geo.CreateGeoMesh());
    const int levels = args.GetI("adapt", 0);
    if (levels <= 0) return gmesh;
    const REAL frac = args.Get("frac", 0.02);
    const int layers = args.GetI("camadas", 1);
    const bool fs = args.GetS("medida", "gamma") == "fs";
    const char *label = fs ? "FS" : "Gamma";
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    const std::string key = AdaptKey(cs, args), file = AdaptFile(cs, args, key);
    const int64_t nel0 = gmesh->NElements(), nnod0 = gmesh->NNodes();
    const uint64_t sig0 = TSlopeGeometry::Signature(*gmesh);
    std::vector<TAdaptLevel> hist;
    int done = 0;
    if (!file.empty() && fileio::Exists(file)) {
        std::string why;
        bool ok = ReadAdaptFile(file, key, nel0, nnod0, sig0, hist, why);
        const int nrep = ok ? std::min<int>(levels, (int)hist.size()) : 0;
        for (int l = 0; ok && l < nrep; l++) {
            const TAdaptLevel &L = hist[l];
            TSlopeGeometry::Refine(gmesh.get(), L.divided);
            if (gmesh->NElements() != L.nel || gmesh->NNodes() != L.nnod || TSlopeGeometry::Signature(*gmesh) != L.sig) {
                ok = false;
                why = "não confere no nível " + std::to_string(l) + " (elementos, nós ou assinatura)";
                break;
            }
            std::cout << "[adapt] nível " << l << ": " << L.neq << " equações, " << label << " = " << L.factor
                      << "; dividindo " << L.nmarked << " + " << (int64_t)L.divided.size() - L.nmarked
                      << " elementos (lido de " << file << ")\n";
        }
        if (ok) {
            done = nrep;
        } else {
            std::cout << "[adapt] arquivo de malha " << file << " " << why << "; recalculando a adaptação\n";
            hist.clear();
            gmesh.reset(cs.geo.CreateGeoMesh());
        }
    }
    for (int l = done; l < levels; l++) {
        auto t0 = Clock::now();
        TAdaptLevel L;
        {
            TSlopeFEM<TMohrCoulomb> fem(gmesh.get(), cs.geo, cs.soil, opt);
            std::unique_ptr<TSeepageProblem> seep;
            if (cs.seepage) {
                TSeepageProblem::TParams sp;
                sp.hw = cs.hw;
                sp.alpha = cs.alpha;
                sp.gammaW = cs.soil.gammaW;
                seep = std::make_unique<TSeepageProblem>(gmesh.get(), cs.geo, sp);
                seep->Solve();
                TransferSeepage(*seep, fem);
            }
            fem.ResetState();
            TFactorResult r = fs ? fem.StrengthReduction(args.Get("F0", 0.5)) : fem.LoadFactor(0.);
            std::vector<REAL> ind;
            fem.PlasticIndicator(ind, args.GetI("incremento", 1) != 0);
            REAL vmax = 0.;
            for (REAL v : ind) vmax = std::max(vmax, v);
            std::vector<int64_t> marked;
            for (TPZGeoEl *gel : gmesh->ElementVec())
                if (gel && gel->Dimension() == 2 && !gel->HasSubElement() && ind[gel->Index()] > frac * vmax)
                    marked.push_back(gel->Index());
            const size_t nmarked = marked.size();
            TSlopeGeometry::Grow(gmesh.get(), marked, layers);
            std::cout << "[adapt] nível " << l << ": " << fem.NEquations() << " equações, " << label << " = "
                      << r.factor << "; dividindo " << nmarked << " + " << marked.size() - nmarked << " elementos ("
                      << Seconds(t0) << " s)\n";
            TSlopeGeometry::Refine(gmesh.get(), marked);
            L.neq = fem.NEquations();
            L.factor = r.factor;
            L.nmarked = (int64_t)nmarked;
            L.divided = marked;
        }
        L.nel = gmesh->NElements();
        L.nnod = gmesh->NNodes();
        L.sig = TSlopeGeometry::Signature(*gmesh);
        hist.resize(l);
        hist.push_back(L);
    }
    if (done < levels && !file.empty()) {
        try {
            WriteAdaptFile(file, key, nel0, nnod0, sig0, hist);
            std::cout << "[adapt] malha adaptada (" << hist.size() << " níveis) gravada em " << file << "\n";
        } catch (std::exception &e) {
            std::cout << "[adapt] malha adaptada não gravada: " << e.what() << "\n";
        }
    }
    return gmesh;
}

void Usage(const char *prog) {
    std::cout << "uso: " << prog << " <comando> [chave=valor ...]\n"
              << "  det caso=cho_coesivo|cho_cphi|percolacao modelo=mc|mcc h=0.5 p=2 hw=5 alpha=1 beta=45 H= c= phi= gam=\n"
              << "      adapt=0 frac=0.02 camadas=1 gamma=1 fs=1 vtk=0 (Cam-Clay: lambda= kappa= v0= OCR=)\n"
              << "      Lc= Lt= Hb= (dimensões do domínio: crista, pé, base abaixo do pé)\n"
              << "      malha=<arquivo>|0: malha adaptada gravada/relida (padrão malha_<caso>_h<h>_<hash>.malha)\n"
              << "  mc  caso=... modelo=mc|mcc medida=gamma|fs n=100 inicio=0 seed=2025 covc=0.3 covphi=0.1 covk=0.6\n"
              << "      Lx=20 Ly=2 hkl=1 M=-1 epsM=-1 normvar=1 saida=arquivo.csv vtk=0 klcache=<arquivo> malha=\n"
              << "      retomada nativa: as amostras de [inicio, inicio+n) que já estão no CSV são puladas (rodar de\n"
              << "      novo o mesmo comando completa o CSV, sem repetições); linhas incompletas são retiradas.\n"
              << "      <csv>.param: parâmetros, conferidos na retomada (forcar=1 aceita diferenças);\n"
              << "      <csv>.resumo: média, desvio, Pf e CoV(Pf) de todo o CSV; <csv>.lock: trava do processo.\n"
              << "      amostras=3,17,40-45: recalcula essas amostras (com vtk=1, figuras) e confere com o CSV\n"
              << "      campos=1: grava c, phi (pontos de integração) e kv (elementos) de cada amostra em <csv>.campos\n"
              << "  campos arquivo=<csv>.campos [amostra=k saida=arquivo.csv]: lê o arquivo de campos (médias por\n"
              << "      amostra, como no CSV, ou os campos da amostra k)\n"
              << "  rebaixamento caso=percolacao modelo=mc|mcc Td=0.1 k=1e-5 Se=0 hw=5 alpha=1 tempos=0.01,0.1,1,3\n"
              << "      fs=1 gamma=0 saida=rebaixamento.csv vtk=0 malha= (u-p acoplado; T = c_v t / H^2)\n"
              << "  códigos de saída: 0 ok, 1 uso, 2 erro/arquivo inválido, 3 parâmetros diferentes dos do CSV,\n"
              << "      4 CSV em uso por outro processo, 5 amostra recalculada diferente da do CSV\n";
}

/// Lista de números separados por vírgula
std::vector<REAL> ParseList(const std::string &str) {
    std::vector<REAL> v;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty()) v.push_back(std::atof(item.c_str()));
    return v;
}

/// Rebaixamento rápido/lento com o problema u-p acoplado (TPZMatPoroElastoPlastic3DMem em deformação plana):
/// em tempos escolhidos a poropressão é congelada e transferida à análise de estabilidade (Γ e/ou FS), como no
/// artigo, mas com p(x, t) do adensamento em vez do fluxo estacionário desacoplado.
template <class T>
void RunDrawdown(const TCase &cs, const TArgs &args, const std::string &modelName) {
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.verbose = args.GetI("verbose", 0);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    const bool doFS = args.GetI("fs", 1) != 0, doGamma = args.GetI("gamma", 0) != 0;
    std::unique_ptr<TPZGeoMesh> gmesh = AdaptedMesh(cs, args);

    typename TCoupledDrawdown<T>::TParams par;
    par.hw = cs.hw;
    par.alpha = cs.alpha;
    par.kv = args.Get("k", 1.e-5);  // artigo: k_v/γw = 1e-6 m⁴/(kN s)
    par.Se = args.Get("Se", 0.);
    par.Td = args.Get("Td", 0.1);
    par.porderU = opt.porder;
    par.nDrawdown = args.GetI("nrebaixamento", 10);
    par.growth = args.Get("crescimento", 1.5);
    par.tol = args.Get("tolup", 1.e-7);
    par.verbose = args.GetI("verbose", 0);

    TSlopeFEM<T> fem(gmesh.get(), cs.geo, cs.soil, opt);
    GeostaticStress(cs, gmesh.get(), opt, fem);
    TCoupledDrawdown<T> cd(gmesh.get(), cs.geo, cs.soil, par);
    if constexpr (std::is_same_v<T, TPZModifiedCamClay>) cd.SetInitialStress(GeostaticStress(cs, gmesh.get(), opt, cd.Points()));
    std::cout << "\n=== rebaixamento acoplado (" << modelName << "): " << cs.geo.Describe() << "\n"
              << "    h_w = " << cs.hw << " m em T_d = " << par.Td << " (t_d = " << par.Td * cd.TimeScale() / 3600.
              << " h), k_v = " << par.kv << " m/s, k_h/k_v = " << par.alpha << ", c_v = " << cd.Cv()
              << " m2/s, H2/c_v = " << cd.TimeScale() / 3600. << " h\n"
              << "    u-p: " << cd.NEquations() << " equacoes; estabilidade: " << fem.NEquations() << " equacoes\n";

    const std::string out = args.GetS("saida", "rebaixamento_" + modelName + ".csv");
    std::ofstream csv(out);
    csv << "Td,T,t_h,zw,u_A,umax_desloc,pontos_plasticos,FS,FS_sup,FS_status,Gamma,Gamma_sup,Gamma_status,tempo_s\n";
    const bool vtk = args.GetI("vtk", 0) != 0;
    if (vtk) cd.DefineVTK(cs.name + "_" + modelName + "_Td" + args.GetS("Td", "0.1"));
    // ponto A: abaixo da borda da crista, a meia altura do talude
    TPZManVector<REAL, 3> xA = {cs.geo.Lc, cs.geo.Hb + 0.5 * cs.geo.H, 0.};
    int vtkStep = 0;
    auto evaluate = [&](const std::string &label) {
        auto t1 = Clock::now();
        cd.TransferSeepage(fem);
        TFactorResult rf, rg;
        if (doFS) {
            fem.ResetState();
            rf = fem.StrengthReduction(args.Get("F0", 0.5));
        }
        if (doGamma) {
            fem.ResetState();
            rg = fem.LoadFactor(0.);
        }
        const REAL Tnow = cd.Time();
        const REAL uA = cd.PorePressureAt(xA) - cs.soil.gammaW * (cs.geo.D() - xA[1]);
        std::cout << "    " << std::setw(10) << label << " T = " << std::setw(9) << Tnow << "  z_w = " << std::setw(6)
                  << cd.WaterLevel() << "  u_A = " << std::setw(8) << uA << "  |du|max = " << std::setw(10)
                  << cd.MaxDisplacementIncrement() << "  plast. " << std::setw(5) << cd.NPlasticPoints();
        if (doFS) std::cout << "  FS = " << rf.factor << " (" << rf.status << ")";
        if (doGamma) std::cout << "  Gamma = " << rg.factor << " (" << rg.status << ")";
        std::cout << "  [" << Seconds(t1) << " s]" << std::endl;
        csv << par.Td << "," << Tnow << "," << Tnow * cd.TimeScale() / 3600. << "," << cd.WaterLevel() << "," << uA << ","
            << cd.MaxDisplacementIncrement() << "," << cd.NPlasticPoints() << "," << rf.factor << "," << rf.upper
            << "," << rf.status << "," << rg.factor << "," << rg.upper << "," << rg.status << "," << Seconds(t1)
            << std::endl;
        if (vtk) cd.WriteVTK(vtkStep++);
    };

    auto t0 = Clock::now();
    // referência do artigo: fluxo estacionário desacoplado (Darcy) com as mesmas condições de contorno finais
    {
        TSeepageProblem::TParams sp;
        sp.hw = cs.hw;
        sp.alpha = cs.alpha;
        sp.gammaW = cs.soil.gammaW;
        TSeepageProblem seep(gmesh.get(), cs.geo, sp);
        seep.Solve();
        TransferSeepage(seep, fem);
        TFactorResult rf, rg;
        if (doFS) {
            fem.ResetState();
            rf = fem.StrengthReduction(args.Get("F0", 0.5));
        }
        if (doGamma) {
            fem.ResetState();
            rg = fem.LoadFactor(0.);
        }
        std::cout << "    estacionario desacoplado (artigo):";
        if (doFS) std::cout << " FS = " << rf.factor << " (" << rf.status << ")";
        if (doGamma) std::cout << " Gamma = " << rg.factor << " (" << rg.status << ")";
        std::cout << std::endl;
        csv << par.Td << ",inf,inf," << cs.geo.D() - cs.hw << ",nan,nan,0," << rf.factor << "," << rf.upper << ","
            << rf.status << "," << rg.factor << "," << rg.upper << "," << rg.status << ",0" << std::endl;
    }
    if (!cd.Initialize()) {
        std::cout << "    equilibrio inicial nao convergiu\n";
        return;
    }
    evaluate("inicial");
    std::vector<REAL> times;
    for (REAL f : {0.25, 0.5, 0.75, 1.}) times.push_back(f * par.Td);
    for (REAL t : ParseList(args.GetS("tempos", "0.01,0.03,0.1,0.3,1,3"))) times.push_back(par.Td + t);
    for (REAL Tout : times) {
        if (!cd.AdvanceTo(Tout, args.Get("dtmax", 0.5))) {
            std::cout << "    colapso no adensamento acoplado em T = " << cd.Time() << " (z_w = " << cd.WaterLevel()
                      << ")\n";
            csv << par.Td << "," << cd.Time() << "," << cd.Time() * cd.TimeScale() / 3600. << "," << cd.WaterLevel()
                << ",nan,nan,0,nan,nan,colapso_acoplado,nan,nan,colapso_acoplado,0" << std::endl;
            break;
        }
        evaluate(Tout <= par.Td * (1. + 1.e-9) ? "rebaixando" : "dissipacao");
    }
    std::cout << "    tempo total " << Seconds(t0) << " s\n";
}

template <class T>
void RunDeterministic(const TCase &cs, const TArgs &args, const std::string &modelName) {
    TSolverOptions opt;
    opt.porder = args.GetI("p", 2);
    opt.verbose = args.GetI("verbose", 0);
    opt.relTol = args.Get("reltol", 5.e-3);
    opt.maxIter = args.GetI("maxit", 20);
    opt.stagnation = args.GetI("estagnacao", 1) != 0;
    std::unique_ptr<TPZGeoMesh> gmesh = AdaptedMesh(cs, args);
    std::cout << "\n=== " << cs.name << " (" << modelName << "): " << cs.geo.Describe() << "\n";
    std::cout << "    c = " << cs.soil.c << " kPa, phi = " << cs.soil.phiDeg << " graus, gamma = " << cs.soil.gamma
              << (cs.soil.buoyant ? " (forca de corpo gamma')" : "") << (cs.seepage ? ", percolacao hw = " : "")
              << (cs.seepage ? std::to_string(cs.hw) : std::string()) << "\n";
    auto t0 = Clock::now();
    TSlopeFEM<T> fem(gmesh.get(), cs.geo, cs.soil, opt);
    GeostaticStress(cs, gmesh.get(), opt, fem);
    std::cout << "    " << fem.NEquations() << " equacoes, " << fem.NPoints() << " pontos de integracao\n";
    if (args.GetI("vtk", 0)) {  // malha (elementos computacionais = folhas da malha adaptada)
        std::ofstream f(cs.name + "_malha.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(fem.Mesh(), f, true);
    }
    std::unique_ptr<TSeepageProblem> seep;
    if (cs.seepage) {
        TSeepageProblem::TParams sp;
        sp.hw = cs.hw;
        sp.alpha = cs.alpha;
        sp.gammaW = cs.soil.gammaW;
        seep = std::make_unique<TSeepageProblem>(gmesh.get(), cs.geo, sp);
        seep->Solve();
        TransferSeepage(*seep, fem);
        if (args.GetI("vtk", 0)) {
            seep->DefineVTK(cs.name + "_darcy.vtk");
            seep->WriteVTK(0);
        }
    }
    if (args.GetI("gamma", 1)) {
        fem.ResetState();
        auto t1 = Clock::now();
        TFactorResult r = fem.LoadFactor(0.);
        std::cout << "    Gamma (fator de carga) = " << r.factor << " (colapso em " << r.upper << "; " << r.steps
                  << " passos, " << r.iterations << " it., " << r.cuts << " cortes, " << Seconds(t1) << " s) "
                  << r.status << "\n";
        if (args.GetI("vtk", 0)) {
            fem.DefineVTK(cs.name + "_" + modelName + "_gamma.vtk");
            fem.WriteVTK(0);
        }
    }
    if (args.GetI("fs", 1)) {
        fem.ResetState();
        auto t1 = Clock::now();
        TFactorResult r = fem.StrengthReduction(args.Get("F0", 0.5));
        std::cout << "    FS (reducao de resistencia) = " << r.factor << " (colapso em " << r.upper << "; "
                  << r.steps << " passos, " << r.iterations << " it., " << r.cuts << " cortes, " << Seconds(t1)
                  << " s) " << r.status << "\n";
        if (args.GetI("vtk", 0)) {
            fem.DefineVTK(cs.name + "_" + modelName + "_fs.vtk");
            fem.WriteVTK(0);
        }
    }
    std::cout << "    tempo total " << Seconds(t0) << " s\n";
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        Usage(argv[0]);
        return 1;
    }
    try {
        const std::string cmd = argv[1];
        TArgs args(argc, argv, 2);
        if (cmd == "campos") return RunFields(args);
        if (cmd == "det" || cmd == "mc" || cmd == "rebaixamento") {
            TCase cs = MakeCase(args.GetS("caso", "percolacao"), args.Get("h", 0.5));
            cs.hw = args.Get("hw", cs.hw);
            cs.alpha = args.Get("alpha", cs.alpha);
            if (args.kv.count("beta")) cs.geo.betaDeg = args.Get("beta", 45.);
            if (args.kv.count("H")) {  // escala a geometria com H (crista, pé e base proporcionais)
                const REAL s = args.Get("H", cs.geo.H) / cs.geo.H;
                cs.geo.H *= s;
                cs.geo.Lc *= s;
                cs.geo.Lt *= s;
                cs.geo.Hb *= s;
                cs.geo.h *= s;
                if (!args.kv.count("hw")) cs.hw *= s;
            }
            // dimensões do domínio (crista, pé e base abaixo do pé), em m
            if (args.kv.count("Lc")) cs.geo.Lc = args.Get("Lc", cs.geo.Lc);
            if (args.kv.count("Lt")) cs.geo.Lt = args.Get("Lt", cs.geo.Lt);
            if (args.kv.count("Hb")) cs.geo.Hb = args.Get("Hb", cs.geo.Hb);
            if (args.kv.count("c")) cs.soil.c = args.Get("c", cs.soil.c);
            if (args.kv.count("phi")) cs.soil.phiDeg = args.Get("phi", cs.soil.phiDeg);
            if (args.kv.count("nu")) cs.soil.nu = args.Get("nu", cs.soil.nu);
            cs.soil.gamma = args.Get("gam", cs.soil.gamma);
            cs.soil.lambda = args.Get("lambda", cs.soil.lambda);
            cs.soil.kappa = args.Get("kappa", cs.soil.kappa);
            cs.soil.v0 = args.Get("v0", cs.soil.v0);
            cs.soil.OCR = args.Get("OCR", cs.soil.OCR);
            if (args.GetS("mapeamento", "deformacao_plana") == "triaxial")
                cs.soil.mapping = TPZModifiedCamClay::ETriaxialCompression;
            if (args.kv.count("E")) cs.soil.E = args.Get("E", cs.soil.E);
            cs.geo.triangles = args.GetI("tri", 0) != 0;
            const std::string model = args.GetS("modelo", "mc");
            if (cmd == "det") {
                if (model == "mcc") RunDeterministic<TPZModifiedCamClay>(cs, args, model);
                else RunDeterministic<TMohrCoulomb>(cs, args, model);
            } else if (cmd == "rebaixamento") {
                if (model == "mcc") RunDrawdown<TPZModifiedCamClay>(cs, args, model);
                else RunDrawdown<TMohrCoulomb>(cs, args, model);
            } else {
                if (model == "mcc") return RunMonteCarlo<TPZModifiedCamClay>(cs, args, model);
                return RunMonteCarlo<TMohrCoulomb>(cs, args, model);
            }
            return 0;
        }
    } catch (TExitError &e) {
        std::cout.flush();
        std::cerr << "erro: " << e.what() << std::endl;
        return e.code;
    } catch (std::exception &e) {
        std::cout.flush();
        std::cerr << "erro: " << e.what() << std::endl;
        return 2;
    }
    Usage(argv[0]);
    return 1;
}
