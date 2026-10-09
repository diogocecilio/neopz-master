// Output of the long production runs (fig5 out=, fig8, fig9, fembatch): resumable CSV files and progress logs, so that
// a run killed at any point (container restarts) is simply restarted and skips the cases already done.
#ifndef SLOPESEEPAGE_RESUMABLECSV_H
#define SLOPESEEPAGE_RESUMABLECSV_H

#include "pzerror.h"
#include "pzreal.h"

#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace slope {

/// number of the CSVs: prec significant digits, "inf" / "-inf" / "nan" when not finite
inline std::string CsvNum(REAL v, int prec = 10) {
    if (!std::isfinite(v)) return v > 0. ? "inf" : (v < 0. ? "-inf" : "nan");
    std::ostringstream s;
    s << std::setprecision(prec) << v;
    return s.str();
}

/// comma-separated list without its empty items
inline std::vector<std::string> SplitList(const std::string &s) {
    std::vector<std::string> v;
    std::string item;
    std::istringstream ss(s);
    while (std::getline(ss, item, ','))
        if (!item.empty()) v.push_back(item);
    return v;
}

/// Resumable CSV: the header is written when the file is new or empty (an existing file must have the same header);
/// a row is identified by its key columns (indices into the header, the settings string being one of them). Rows are
/// appended and flushed one at a time, so a run killed at any point loses at most the case in progress; a row cut by
/// the kill (wrong number of columns) is ignored, and a missing final newline is restored before the next append.
class ResumableCSV {
    std::string fFile;
    std::vector<size_t> fKey;
    size_t fNCols = 0;
    std::set<std::string> fDone;
    bool fNeedNewline = false;
    int fIgnored = 0;

public:
    /// the columns of a line (empty columns kept)
    static std::vector<std::string> Split(const std::string &line) {
        std::vector<std::string> v;
        std::string item;
        std::istringstream ss(line);
        while (std::getline(ss, item, ',')) v.push_back(item);
        if (!line.empty() && line.back() == ',') v.push_back("");
        return v;
    }
    std::string Key(const std::vector<std::string> &cols) const {
        std::string k;
        for (size_t i : fKey) k += (i < cols.size() ? cols[i] : std::string("?")) + "|";
        return k;
    }
    ResumableCSV(const std::string &file, const std::string &header, const std::vector<size_t> &key)
        : fFile(file), fKey(key), fNCols(Split(header).size()) {
        const std::filesystem::path dir = std::filesystem::path(file).parent_path();
        if (!dir.empty()) std::filesystem::create_directories(dir);
        bool haveHeader = false;
        {
            std::ifstream in(file);
            std::string line;
            while (std::getline(in, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;
                if (!haveHeader) {
                    if (line != header) {
                        std::cerr << file << ": the header differs from this version's\n  " << header << "\n(use another out=)\n";
                        exit(1);
                    }
                    haveHeader = true;
                    continue;
                }
                const std::vector<std::string> cols = Split(line);
                if (cols.size() == fNCols) fDone.insert(Key(cols));
                else fIgnored++;
            }
        }
        if (!haveHeader) std::ofstream(file, std::ios::trunc) << header << "\n";
        else { // the last row may have been cut by a kill: the next one starts on a new line
            std::ifstream f(file, std::ios::binary | std::ios::ate);
            const std::streamoff size = f.tellg();
            if (size > 0) {
                f.seekg(size - 1);
                char c = '\n';
                f.get(c);
                fNeedNewline = c != '\n';
            }
        }
    }
    const std::string &File() const { return fFile; }
    int Ignored() const { return fIgnored; }
    size_t NDone() const { return fDone.size(); }
    bool Done(const std::vector<std::string> &cols) const { return fDone.count(Key(cols)) > 0; }
    void Append(const std::vector<std::string> &cols) {
        if (cols.size() != fNCols) {
            std::cerr << "ResumableCSV: " << cols.size() << " columns for a header of " << fNCols << "\n";
            DebugStop();
        }
        std::string line = fNeedNewline ? "\n" : "";
        for (size_t i = 0; i < cols.size(); i++) line += (i ? "," : "") + cols[i];
        std::ofstream f(fFile, std::ios::app);
        f << line << "\n";
        f.flush();
        fNeedNewline = false;
        fDone.insert(Key(cols));
    }
};

/// progress log: each message to stdout and, with a UTC time stamp, appended to a file
struct ProgressLog {
    std::string file;
    void operator()(const std::string &msg) const {
        const std::time_t t = std::time(nullptr);
        char stamp[32];
        std::strftime(stamp, sizeof stamp, "%Y-%m-%d %H:%M:%S", std::gmtime(&t));
        std::cout << msg << std::endl;
        if (!file.empty()) std::ofstream(file, std::ios::app) << stamp << "  " << msg << "\n";
    }
};

/// the .log file next to a .csv
inline std::string LogFileOf(const std::string &csv) {
    const size_t n = csv.size();
    return (n > 4 && csv.compare(n - 4, 4, ".csv") == 0 ? csv.substr(0, n - 4) : csv) + ".log";
}

} // namespace slope

#endif
