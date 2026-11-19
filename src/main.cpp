// netevolve command-line entry point.
#include <iostream>
#include <string_view>

namespace {

void print_usage(std::ostream& os) {
    os << "usage:\n"
          "  netevolve evolve   --config FILE [--seed N] [--out DIR]\n"
          "  netevolve simulate --model FILE.ant [--t-end T] [--points N] [--csv OUT.csv]\n"
          "  netevolve fitness  --model FILE.ant [--config FILE]\n"
          "  netevolve prune    --model FILE.ant [--config FILE] [--out FILE.ant]\n"
          "  netevolve --version\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(std::cerr);
        return 2;
    }
    const std::string_view cmd = argv[1];
    if (cmd == "--version") {
        std::cout << "netevolve " << NETEVOLVE_VERSION << "\n";
        return 0;
    }
    if (cmd == "--help" || cmd == "-h") {
        print_usage(std::cout);
        return 0;
    }
    std::cerr << "error: unknown command '" << cmd << "'\n";
    print_usage(std::cerr);
    return 2;
}
