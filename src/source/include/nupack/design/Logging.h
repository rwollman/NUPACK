#pragma once
#include "../reflect/Reflection.h"
#include "../reflect/Serialize.h"

namespace nupack { namespace newdesign {

template <class T, NUPACK_IF(!can_call<T>)>
T result(T &&t) {return static_cast<T&&>(t);}

template <class T, NUPACK_IF(can_call<T>)>
auto result(T &&t) {return t();}


struct Logger {
    string filename;
    std::shared_ptr<std::ostream> out;
    bool active;

    Logger() : Logger("") {}
    Logger(string filename) : filename(filename), out(), active(false) {
        if (filename == "") {
            return;
        } else if (filename == "stdout") {
            out = std::shared_ptr<std::ostream>(&std::cout, [](auto){});
        } else if (filename == "stderr") {
            out = std::shared_ptr<std::ostream>(&std::cerr, [](auto){});
        } else {
            out = std::make_shared<std::ofstream>(filename);
        }
        active = true;
    }

    template <class ...Ts>
    void log(Ts &&...ts) {
        if (active) print_os<io::character<';'>, io::endl>(*out, result(fw<Ts>(ts))...);
    }

    NUPACK_REFLECT(Logger, filename, out, active);

    auto save_repr() const { return filename; }

    void load_repr(string filename) { *this = Logger(filename); }
};


struct Logs {
    std::map<string, Logger> loggers;

    Logs() = default;
    Logs(std::map<string, string> active_logs) {
        for (auto const & [k, v] : active_logs)
            loggers.emplace(k, Logger(v));
    }

    template <class ...Ts>
    void log(string particular, Ts &&...ts) {
        if (loggers.count(particular))
            at(loggers, particular).log(fw<Ts>(ts)...);
    }

    auto save_repr() const { return loggers; }

    void load_repr(std::map<string, Logger> loggers) { this->loggers = loggers; }
};

static_assert(nupack::is_nupack<Logs> && nupack::has_save_repr<Logs> && !nupack::has_repr_names<Logs>);
// static_assert(std::is_convertible_v<Logs, json>);
// static_assert(!std::is_convertible_v<Logs, json>);
static_assert(std::is_constructible_v<json, Logs>);
// static_assert(!std::is_constructible_v<json, Logs>);


struct EngineObserver {
    uint slowdown {0};
    Logs *logs;

    template <class ...Ts>
    void log(Ts &&...ts) { if (logs) logs->log(fw<Ts>(ts)...); }
};

extern EngineObserver NullEngineObserver;

}}
