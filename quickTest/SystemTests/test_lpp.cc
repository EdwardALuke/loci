#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

namespace {

  std::string read_file(const std::string &path) {
    std::ifstream in(path.c_str());
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
  }

  void write_file(const std::string &path, const std::string &content) {
    std::ofstream out(path.c_str());
    out << content;
  }

  struct LppRun {
    std::string input;
    std::string output;
    std::string log;
    int rc;

    ~LppRun() {
      std::remove(input.c_str());
      std::remove(output.c_str());
      std::remove(log.c_str());
    }
  };

  LppRun run_lpp(const std::string &name,
                 const std::string &source,
                 const std::string &extra_args = "") {
    const std::string base =
      std::string("test_lpp_") + std::to_string(getpid()) + "_" + name;
    LppRun run;
    run.input = base + ".loci";
    run.output = base + ".cc";
    run.log = base + ".log";

    write_file(run.input, source);

    const char *loci_base = std::getenv("LOCI_BASE");
    std::string command =
      std::string(loci_base ? loci_base : "../..") + "/bin/lpp -x ";
    if(!extra_args.empty())
      command += extra_args + " ";
    command += run.input + " -o " + run.output + " > " + run.log + " 2>&1";
    run.rc = std::system(command.c_str());
    return run;
  }

  const char *basic_rule_source =
    "#include <Loci.h>\n"
    "namespace Loci {\n"
    "$type lpp_src store<int>;\n"
    "$type lpp_tgt store<int>;\n"
    "$rule pointwise(lpp_tgt<-lpp_src) {\n"
    "  $lpp_tgt = $lpp_src + 1;\n"
    "}\n"
    "}\n";

} // namespace

//----------------------------------------------------------------------------
// Regression Tests For Known Bugs
//----------------------------------------------------------------------------

// Current regression:
// The `-DNAME=value` option records the system variable name as `NAME=value`
// instead of `NAME`, so `$ifdef NAME` does not see it.
//
// Suggested fix in `src/lpp/main.cc`:
// Stop appending characters to `var` after the `=` is encountered, or parse the
// option into separate name/value substrings before calling `setSystemVar()`.
TEST_CASE("lpp -DNAME=value defines NAME for conditional parsing"
          * doctest::may_fail()) {
  LppRun run = run_lpp(
    "dash_d",
    "#include <Loci.h>\n"
    "namespace Loci {\n"
    "$ifdef LPP_ENABLE_TEST\n"
    "$type lpp_dash_d_enabled store<int>;\n"
    "$endif\n"
    "}\n",
    "-DLPP_ENABLE_TEST=1");

  REQUIRE(run.rc == 0);
  CHECK(read_file(run.output).find("lpp_dash_d_enabled") != std::string::npos);
}

//----------------------------------------------------------------------------
// Additional behavior/invariant tests
//----------------------------------------------------------------------------

TEST_CASE("lpp translates a pointwise rule into registered C++ rule code") {
  LppRun run = run_lpp("pointwise", basic_rule_source);

  REQUIRE(run.rc == 0);
  const std::string generated = read_file(run.output);

  CHECK(generated.find("public Loci::pointwise_rule") != std::string::npos);
  CHECK(generated.find("name_store(\"lpp_src\"") != std::string::npos);
  CHECK(generated.find("name_store(\"lpp_tgt\"") != std::string::npos);
  CHECK(generated.find("input(\"lpp_src\")") != std::string::npos);
  CHECK(generated.find("output(\"lpp_tgt\")") != std::string::npos);
  CHECK(generated.find("register_rule") != std::string::npos);
  CHECK(generated.find("L_lpp_tgt_[_e_]= L_lpp_src_[_e_]+ 1") != std::string::npos);
}

TEST_CASE("lpp in-file define controls conditional blocks") {
  LppRun run = run_lpp(
    "define",
    "#include <Loci.h>\n"
    "namespace Loci {\n"
    "$define LPP_LOCAL_TEST \"1\"\n"
    "$ifdef LPP_LOCAL_TEST\n"
    "$type lpp_local_enabled store<int>;\n"
    "$else\n"
    "$type lpp_local_disabled store<int>;\n"
    "$endif\n"
    "}\n");

  REQUIRE(run.rc == 0);
  const std::string generated = read_file(run.output);

  CHECK(generated.find("lpp_local_enabled") != std::string::npos);
  CHECK(generated.find("lpp_local_disabled") == std::string::npos);
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
