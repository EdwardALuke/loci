#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/parse.h>

#ifndef LOCI_ENABLE_PARSE_KNOWN_BUG_TESTS
#define LOCI_ENABLE_PARSE_KNOWN_BUG_TESTS 0
#endif

namespace {
std::string read_remaining(std::istream &in) {
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}
} // namespace

//----------------------------------------------------------------------------
// Potential Bugs
//----------------------------------------------------------------------------
//
// These are cases that appear to demonstrate behavior that is not intended.
//
// Enable them with LOCI_ENABLE_PARSE_KNOWN_BUG_TESTS=1 when we want to
// revisit these behaviors.

// The parser already treats "12ex" as the valid prefix "12" with "ex" left in
// the stream. Malformed exponent tails such as "e+" and "e-" should behave the
// same way instead of being silently consumed.
TEST_CASE("known bug: get_real should not consume incomplete exponent suffixes [known-bug]" *
          doctest::skip(LOCI_ENABLE_PARSE_KNOWN_BUG_TESTS == 0)) {
  struct RealCase {
    const char *text;
    const char *rest;
  };

  const std::vector<RealCase> cases = {
      {"1e+", "e+"},
      {"1e-", "e-"},
  };

  for(const auto &tc : cases) {
    CAPTURE(tc.text);
    std::istringstream in(tc.text);

    CHECK(Loci::parse::is_real(in));
    CHECK(Loci::parse::get_real(in) == doctest::Approx(1.0));
    CHECK(read_remaining(in) == tc.rest);
  }
}

// is_int() and is_real() currently treat a leading + or - as sufficient to
// begin a number, even when no digits follow. That leaks upward into callers
// like options_list and UNIT_type, which can accept "+" or "-" as zero.
TEST_CASE("known bug: bare signs should not be classified as numeric literals [known-bug]" *
          doctest::skip(LOCI_ENABLE_PARSE_KNOWN_BUG_TESTS == 0)) {
  const std::vector<std::string> cases = {
      "+",
      "-",
      "+foo",
      "-foo",
  };

  for(const auto &tc : cases) {
    CAPTURE(tc);
    std::istringstream int_in(tc);
    std::istringstream real_in(tc);

    CHECK_FALSE(Loci::parse::is_int(int_in));
    CHECK_FALSE(Loci::parse::is_real(real_in));
  }
}

TEST_CASE("kill_white_space skips leading whitespace and line comments") {
  std::istringstream in(" \t\n// first comment\n  // second comment\r\nname");

  Loci::parse::kill_white_space(in);

  CHECK(in.peek() == 'n');
  CHECK(Loci::parse::get_name(in) == "name");
  CHECK(in.peek() == EOF);
}

TEST_CASE("kill_white_space preserves a slash that does not start a line comment") {
  std::istringstream in("   /value");

  Loci::parse::kill_white_space(in);

  CHECK(in.peek() == '/');
  CHECK_FALSE(Loci::parse::is_token(in, "//"));
  CHECK(read_remaining(in) == "/value");
}

TEST_CASE("is_name and get_name parse identifiers and leave trailing delimiters") {
  std::istringstream in("  // comment\n_alpha42+rest");

  CHECK(Loci::parse::is_name(in));
  CHECK(Loci::parse::get_name(in) == "_alpha42");
  CHECK(in.peek() == '+');
}

TEST_CASE("get_name returns empty string without consuming non-name input") {
  std::istringstream in(" 9alpha");

  CHECK_FALSE(Loci::parse::is_name(in));
  CHECK(Loci::parse::get_name(in).empty());
  CHECK(in.peek() == '9');
}

TEST_CASE("is_int and get_int parse signed integers") {
  std::istringstream in(" \t+42,rest");

  CHECK(Loci::parse::is_int(in));
  CHECK(Loci::parse::get_int(in) == 42);
  CHECK(in.peek() == ',');
}

TEST_CASE("get_int returns zero without consuming invalid input") {
  std::istringstream in("alpha");

  CHECK_FALSE(Loci::parse::is_int(in));
  CHECK(Loci::parse::get_int(in) == 0);
  CHECK(in.peek() == 'a');
}

TEST_CASE("is_real and get_real parse a broad set of numeric forms") {
  struct RealCase {
    const char *text;
    double value;
    const char *rest;
  };

  const std::vector<RealCase> cases = {
      {"0", 0.0, ""},
      {"-1.25e2", -125.0, ""},
      {"+.5", 0.5, ""},
      {"7.", 7.0, ""},
      {".75", 0.75, ""},
      {"3E-2", 0.03, ""},
      {".", 0.0, ""},
      {"12.5x", 12.5, "x"},
  };

  for(const auto &tc : cases) {
    CAPTURE(tc.text);
    std::istringstream in(tc.text);

    CHECK(Loci::parse::is_real(in));
    CHECK(Loci::parse::get_real(in) == doctest::Approx(tc.value));
    CHECK(read_remaining(in) == tc.rest);
  }
}

TEST_CASE("get_real ignores invalid exponents and leaves the exponent marker in the stream") {
  std::istringstream in("12ex");

  CHECK(Loci::parse::get_real(in) == doctest::Approx(12.0));
  CHECK(in.peek() == 'e');
  CHECK(read_remaining(in) == "ex");
}

TEST_CASE("get_real returns zero without consuming invalid non-numeric input") {
  std::istringstream in("name");

  CHECK_FALSE(Loci::parse::is_real(in));
  CHECK(Loci::parse::get_real(in) == doctest::Approx(0.0));
  CHECK(in.peek() == 'n');
}

TEST_CASE("is_string and get_string parse quoted text including slashes and spaces") {
  std::istringstream in(" // comment\n\"a // b c\" tail");

  CHECK(Loci::parse::is_string(in));
  CHECK(Loci::parse::get_string(in) == "a // b c");

  Loci::parse::kill_white_space(in);
  CHECK(Loci::parse::get_name(in) == "tail");
}

TEST_CASE("get_string supports empty quoted strings and leaves following text") {
  std::istringstream in("\"\"next");

  CHECK(Loci::parse::get_string(in).empty());
  CHECK(Loci::parse::get_name(in) == "next");
}

TEST_CASE("get_string returns empty string without consuming non-string input") {
  std::istringstream in("name");

  CHECK_FALSE(Loci::parse::is_string(in));
  CHECK(Loci::parse::get_string(in).empty());
  CHECK(in.peek() == 'n');
}

TEST_CASE("get_string returns collected text when the closing quote is missing") {
  std::istringstream in("\"unterminated");

  CHECK(Loci::parse::get_string(in) == "unterminated");
  CHECK(in.eof());
}

TEST_CASE("is_token checks for a token without consuming it") {
  std::istringstream in("  &&rest");

  CHECK(Loci::parse::is_token(in, "&&"));
  CHECK(in.peek() == '&');
  CHECK(Loci::parse::get_token(in, "&&"));
  CHECK(Loci::parse::get_name(in) == "rest");
}

TEST_CASE("get_token skips leading whitespace and comments before matching") {
  std::istringstream in(" // comment\n::name");

  CHECK(Loci::parse::get_token(in, "::"));
  CHECK(Loci::parse::get_name(in) == "name");
}

TEST_CASE("is_token and get_token restore the stream after a partial mismatch") {
  {
    std::istringstream in("=>");
    CHECK_FALSE(Loci::parse::is_token(in, "=="));
    CHECK(read_remaining(in) == "=>");
  }

  {
    std::istringstream in("=>");
    CHECK_FALSE(Loci::parse::get_token(in, "=="));
    CHECK(read_remaining(in) == "=>");
  }
}
