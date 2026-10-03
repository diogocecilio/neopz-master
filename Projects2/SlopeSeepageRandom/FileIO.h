// FileIO.h — gravação robusta dos arquivos do SlopeSeepageRandom (CSV do Monte Carlo, parâmetros, malha
// adaptada, cache da KL, campos das amostras):
//
//   * escrita atômica: grava <arquivo>.tmp<pid>, fsync e std::rename sobre o destino (um processo que lê o
//     arquivo, ou que é interrompido, nunca vê um arquivo pela metade);
//   * acréscimo durável: cada linha/registro é escrito com uma única chamada write() em O_APPEND, seguido de
//     fdatasync (um SIGKILL deixa no máximo a última linha incompleta, que a retomada descarta);
//   * trava exclusiva (flock) em <arquivo>.lock: dois processos nunca acrescentam ao mesmo CSV (a trava é
//     liberada pelo sistema quando o processo termina, inclusive por SIGKILL);
//   * FNV-1a de 64 bits para assinaturas (malhas, conteúdo dos arquivos binários).
//
#ifndef SSR_FILEIO_H
#define SSR_FILEIO_H

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>

namespace fileio {

/// FNV-1a de 64 bits (bytes)
inline uint64_t Fnv1a(const void *data, size_t n, uint64_t h = 1469598103934665603ULL) {
    const unsigned char *p = static_cast<const unsigned char *>(data);
    for (size_t i = 0; i < n; i++) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

/// FNV-1a sobre palavras de 64 bits (padrão de bits dos valores: assinatura exata, rápida para vetores grandes)
template <class T>
inline uint64_t HashWords(const T *v, size_t n, uint64_t h) {
    static_assert(sizeof(T) == 8, "HashWords: tipo de 8 bytes");
    for (size_t i = 0; i < n; i++) {
        uint64_t w;
        std::memcpy(&w, v + i, 8);
        h ^= w;
        h *= 1099511628211ULL;
    }
    return h;
}

/// 8 caracteres como uma palavra de 64 bits (marcas de formato: o TPZBFileStream só instancia leitura/escrita de
/// int, int64, uint64 e double; os bytes gravados são os mesmos 8 caracteres)
inline uint64_t Word8(const char *m) {
    uint64_t w;
    std::memcpy(&w, m, 8);
    return w;
}

inline std::string Hex(uint64_t h) {
    char s[17];
    std::snprintf(s, sizeof(s), "%016llx", (unsigned long long)h);
    return s;
}

inline bool Exists(const std::string &path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

inline bool FileSize(const std::string &path, uint64_t &size) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return false;
    size = (uint64_t)st.st_size;
    return true;
}

/// Nome temporário exclusivo do processo, no mesmo diretório do destino (rename atômico)
inline std::string TmpName(const std::string &path) { return path + ".tmp" + std::to_string((long)::getpid()); }

inline bool FsyncPath(const std::string &path) {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) return false;
    const bool ok = ::fsync(fd) == 0;
    ::close(fd);
    return ok;
}

inline std::string DirName(const std::string &path) {
    const size_t p = path.find_last_of('/');
    if (p == std::string::npos) return ".";
    return p == 0 ? "/" : path.substr(0, p);
}

/// fsync do arquivo temporário, rename sobre o destino e fsync do diretório
inline void Commit(const std::string &tmp, const std::string &dst) {
    if (!FsyncPath(tmp)) {
        ::unlink(tmp.c_str());
        throw std::runtime_error("não foi possível gravar " + tmp + ": " + std::strerror(errno));
    }
    if (std::rename(tmp.c_str(), dst.c_str()) != 0) {
        const int e = errno;
        ::unlink(tmp.c_str());
        throw std::runtime_error("rename " + tmp + " -> " + dst + ": " + std::strerror(e));
    }
    FsyncPath(DirName(dst));
}

/// Substitui o conteúdo do arquivo de forma atômica
inline void WriteAtomic(const std::string &path, const std::string &content) {
    const std::string tmp = TmpName(path);
    {
        std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
        f << content;
        f.flush();
        if (!f) {
            ::unlink(tmp.c_str());
            throw std::runtime_error("não foi possível gravar " + tmp);
        }
    }
    Commit(tmp, path);
}

inline bool ReadAll(const std::string &path, std::string &content) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    std::ostringstream s;
    s << f.rdbuf();
    content = s.str();
    return !f.bad();
}

/// Arquivo aberto em O_APPEND: cada Append é uma única chamada write() (repetida só se parcial) + fdatasync
class TAppendFile {
public:
    explicit TAppendFile(const std::string &path) : fPath(path) {
        fFd = ::open(path.c_str(), O_WRONLY | O_APPEND | O_CREAT, 0644);
        if (fFd < 0) throw std::runtime_error("não foi possível abrir " + path + ": " + std::strerror(errno));
    }
    ~TAppendFile() {
        if (fFd >= 0) ::close(fFd);
    }
    TAppendFile(const TAppendFile &) = delete;
    TAppendFile &operator=(const TAppendFile &) = delete;
    void Append(const void *data, size_t n, bool sync = true) {
        const char *p = static_cast<const char *>(data);
        while (n > 0) {
            const ssize_t w = ::write(fFd, p, n);
            if (w < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error("erro ao gravar " + fPath + ": " + std::strerror(errno));
            }
            p += w;
            n -= (size_t)w;
        }
        if (sync) ::fdatasync(fFd);
    }
    void Append(const std::string &s, bool sync = true) { Append(s.data(), s.size(), sync); }

private:
    std::string fPath;
    int fFd = -1;
};

/// Trava exclusiva (flock, não bloqueante) em <path>; liberada no destrutor ou no fim do processo
class TFileLock {
public:
    /// devolve false (sem lançar) se outro processo detém a trava
    bool Acquire(const std::string &path) {
        fFd = ::open(path.c_str(), O_RDWR | O_CREAT, 0644);
        if (fFd < 0) throw std::runtime_error("não foi possível criar " + path + ": " + std::strerror(errno));
        if (::flock(fFd, LOCK_EX | LOCK_NB) != 0) {
            ::close(fFd);
            fFd = -1;
            return false;
        }
        // pid do dono (só informativo)
        const std::string pid = std::to_string((long)::getpid()) + "\n";
        if (::ftruncate(fFd, 0) == 0) {
            const ssize_t w = ::write(fFd, pid.data(), pid.size());
            (void)w;
        }
        return true;
    }
    ~TFileLock() {
        if (fFd >= 0) ::close(fFd);
    }

private:
    int fFd = -1;
};

} // namespace fileio

#endif
