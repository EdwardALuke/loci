#include <Loci.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

namespace {

  /*
    These tests exercise the lpp executable as a small integration boundary:
    write tiny .loci inputs, run bin/lpp -x, then inspect stable generated C++
    markers without compiling the generated files.
  */

  std::string read_file(const std::string &path) {
    CAPTURE(path);
    std::ifstream in(path.c_str());
    REQUIRE(in.good());
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
  }

  void write_file(const std::string &path, const std::string &content) {
    CAPTURE(path);
    std::ofstream out(path.c_str());
    REQUIRE(out.good());
    out << content;
    REQUIRE(out.good());
  }

  std::string run_base(const std::string &name) {
    return std::string("test_lpp_") + std::to_string(getpid()) + "_" + name;
  }

  std::string shell_quote(const std::string &value) {
    std::string quoted = "'";
    for(size_t i = 0; i < value.size(); ++i) {
      if(value[i] == '\'')
        quoted += "'\\''";
      else
        quoted += value[i];
    }
    quoted += "'";
    return quoted;
  }

  std::string compact_whitespace(const std::string &text) {
    std::string compact;
    for(size_t i = 0; i < text.size(); ++i) {
      const unsigned char ch = static_cast<unsigned char>(text[i]);
      if(!std::isspace(ch))
        compact += text[i];
    }
    return compact;
  }

  void check_contains(const std::string &text, const std::string &needle) {
    CAPTURE(needle);
    CHECK(text.find(needle) != std::string::npos);
  }

  void check_not_contains(const std::string &text, const std::string &needle) {
    CAPTURE(needle);
    CHECK(text.find(needle) == std::string::npos);
  }

  struct LppRun {
    std::string input;
    std::string output;
    std::string log;
    std::vector<std::string> extra_files;
    int rc;

    LppRun() : rc(-1) {}

    LppRun(const LppRun&) = delete;
    LppRun& operator=(const LppRun&) = delete;

    LppRun(LppRun &&other) {
      move_from(other);
    }

    LppRun& operator=(LppRun &&other) {
      if(this != &other) {
        cleanup();
        move_from(other);
      }
      return *this;
    }

    ~LppRun() {
      cleanup();
    }

  private:
    void cleanup() {
      for(size_t i = 0; i < extra_files.size(); ++i) {
        if(!extra_files[i].empty())
          std::remove(extra_files[i].c_str());
      }
      if(!input.empty())
        std::remove(input.c_str());
      if(!output.empty())
        std::remove(output.c_str());
      if(!log.empty())
        std::remove(log.c_str());
    }

    void move_from(LppRun &other) {
      input = std::move(other.input);
      output = std::move(other.output);
      log = std::move(other.log);
      extra_files = std::move(other.extra_files);
      rc = other.rc;

      other.input.clear();
      other.output.clear();
      other.log.clear();
      other.extra_files.clear();
    }
  };

  LppRun run_lpp(const std::string &name,
                 const std::string &source,
                 const std::string &extra_args = "") {
    const std::string base = run_base(name);
    LppRun run;
    run.input = base + ".loci";
    run.output = base + ".cc";
    run.log = base + ".log";

    write_file(run.input, source);

    const char *loci_base = std::getenv("LOCI_BASE");
    std::string command =
      shell_quote(std::string(loci_base ? loci_base : "../..") + "/bin/lpp")
      + " -x ";
    if(!extra_args.empty())
      command += extra_args + " ";
    command += shell_quote(run.input) + " -o " + shell_quote(run.output)
      + " > " + shell_quote(run.log) + " 2>&1";
    run.rc = std::system(command.c_str());
    return run;
  }

  void require_lpp_success(const LppRun &run) {
    CAPTURE(run.input);
    INFO("lpp log:\n" << read_file(run.log));
    REQUIRE(run.rc == 0);
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

// Current expected-failing regression: `-DNAME=value` records `NAME=value` as
// the variable name, so `$ifdef NAME` does not see the command-line definition.
TEST_CASE("lpp documents current -DNAME=value conditional bug"
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

  require_lpp_success(run);
  check_contains(read_file(run.output), "lpp_dash_d_enabled");
}

//----------------------------------------------------------------------------
// Additional behavior/invariant tests
//----------------------------------------------------------------------------

TEST_CASE("lpp translates a pointwise rule into registered C++ rule code") {
  LppRun run = run_lpp("pointwise", basic_rule_source);

  require_lpp_success(run);
  const std::string generated = read_file(run.output);
  const std::string compact = compact_whitespace(generated);

  check_contains(generated, "public Loci::pointwise_rule");
  check_contains(generated, "name_store(\"lpp_src\"");
  check_contains(generated, "name_store(\"lpp_tgt\"");
  check_contains(generated, "input(\"lpp_src\")");
  check_contains(generated, "output(\"lpp_tgt\")");
  check_contains(generated, "register_rule");
  check_contains(compact, "L_lpp_tgt_[_e_]=L_lpp_src_[_e_]+1;");
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

  require_lpp_success(run);
  const std::string generated = read_file(run.output);

  check_contains(generated, "lpp_local_enabled");
  check_not_contains(generated, "lpp_local_disabled");
}

TEST_CASE("lpp -DNAME without value enables conditional parsing") {
  LppRun run = run_lpp(
    "dash_d_flag",
    "#include <Loci.h>\n"
    "namespace Loci {\n"
    "$ifdef LPP_COMMAND_LINE_TEST\n"
    "$type lpp_cli_enabled store<int>;\n"
    "$else\n"
    "$type lpp_cli_disabled store<int>;\n"
    "$endif\n"
    "$ifndef LPP_COMMAND_LINE_TEST\n"
    "$type lpp_cli_unexpected store<int>;\n"
    "$endif\n"
    "}\n",
    "-DLPP_COMMAND_LINE_TEST");

  require_lpp_success(run);
  const std::string generated = read_file(run.output);

  check_contains(generated, "lpp_cli_enabled");
  check_not_contains(generated, "lpp_cli_disabled");
  check_not_contains(generated, "lpp_cli_unexpected");
}

TEST_CASE("lpp preserves rule constraints conditionals and comments") {
  LppRun run = run_lpp(
    "rule_metadata",
    "#include <Loci.h>\n"
    "namespace Loci {\n"
    "$type lpp_src store<int>;\n"
    "$type lpp_tgt store<int>;\n"
    "$type lpp_guard store<int>;\n"
    "$type lpp_cond param<bool>;\n"
    "$rule pointwise(lpp_tgt<-lpp_src), constraint(lpp_guard), "
    "conditional(lpp_cond), comments(\"metadata comment\") {\n"
    "  $lpp_tgt = $lpp_src;\n"
    "}\n"
    "}\n");

  require_lpp_success(run);
  const std::string generated = read_file(run.output);

  check_contains(generated, "constraint(\"lpp_guard\")");
  check_contains(generated, "conditional(\"lpp_cond\")");
  check_contains(generated, "comments(\"metadata comment\")");
}

TEST_CASE("lpp includes type declarations only once") {
  const std::string include_file = run_base("include_once") + ".lh";
  write_file(include_file,
             "$type lpp_included_src store<int>;\n"
             "$type lpp_included_tgt store<int>;\n");

  LppRun run = run_lpp(
    "include_once",
    "#include <Loci.h>\n"
    "namespace Loci {\n"
    "$include \"" + include_file + "\"\n"
    "$include \"" + include_file + "\"\n"
    "$rule pointwise(lpp_included_tgt<-lpp_included_src) {\n"
    "  $lpp_included_tgt = $lpp_included_src;\n"
    "}\n"
    "}\n");
  run.extra_files.push_back(include_file);

  require_lpp_success(run);
  const std::string generated = read_file(run.output);

  check_contains(generated, "input(\"lpp_included_src\")");
  check_contains(generated, "output(\"lpp_included_tgt\")");
  check_contains(generated, "//$include \"" + include_file + "\"");
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  Loci::Finalize();
  return result;
}
