#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/parse.h>

namespace {
/// @brief Returns every unread character left in `in`.
std::string read_remaining(std::istream &in) {
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}
} // namespace

//----------------------------------------------------------------------------
// Potential Bugs
//----------------------------------------------------------------------------
//
// These stay in the suite as visible regressions without breaking the overall
// test run.

// `get_real()` currently consumes incomplete exponent suffixes such as `e+`
// instead of leaving them unread for the caller.
TEST_CASE("get_real leaves incomplete exponent suffixes unread [known-bug]" *
          doctest::may_fail()) {
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

    // The numeric probe should accept the leading real value.
    CHECK(Loci::parse::is_real(in));
    CHECK(Loci::parse::get_real(in) == doctest::Approx(1.0));

    // The incomplete exponent suffix should remain unread for the caller.
    CHECK(read_remaining(in) == tc.rest);
  }
}

// `is_int()` and `is_real()` currently accept bare signs and sign-prefixed
// identifiers even when no digits follow.
TEST_CASE("numeric probes reject bare signs without digits [known-bug]" *
          doctest::may_fail()) {
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

    // Bare signs and sign-prefixed identifiers are not valid integers or reals.
    CHECK_FALSE(Loci::parse::is_int(int_in));
    CHECK_FALSE(Loci::parse::is_real(real_in));
  }
}

TEST_CASE("kill_white_space skips leading whitespace and line comments") {
  std::istringstream in(" \t\n// first comment\n  // second comment\r\nname");

  // Leading whitespace and whole-line comments should be skipped together.
  Loci::parse::kill_white_space(in);

  // The next token in the stream should now be the identifier itself.
  CHECK(in.peek() == 'n');
  CHECK(Loci::parse::get_name(in) == "name");
  CHECK(in.peek() == EOF);
}

TEST_CASE("kill_white_space preserves a slash that does not start a line comment") {
  std::istringstream in("   /value");

  // A slash that is not the start of `//` should be left alone.
  Loci::parse::kill_white_space(in);

  CHECK(in.peek() == '/');
  CHECK_FALSE(Loci::parse::is_token(in, "//"));
  CHECK(read_remaining(in) == "/value");
}

TEST_CASE("is_name and get_name parse identifiers and leave trailing delimiters") {
  std::istringstream in("  // comment\n_alpha42+rest");

  // Identifier parsing should skip comments and stop cleanly at the delimiter.
  CHECK(Loci::parse::is_name(in));
  CHECK(Loci::parse::get_name(in) == "_alpha42");
  CHECK(in.peek() == '+');
}

TEST_CASE("get_name returns empty string without consuming non-name input") {
  std::istringstream in(" 9alpha");

  // A non-identifier head should be rejected without advancing the stream.
  CHECK_FALSE(Loci::parse::is_name(in));
  CHECK(Loci::parse::get_name(in).empty());
  CHECK(in.peek() == '9');
}

TEST_CASE("is_int and get_int parse signed integers") {
  std::istringstream in(" \t+42,rest");

  // Signed integers should parse after leading whitespace.
  CHECK(Loci::parse::is_int(in));
  CHECK(Loci::parse::get_int(in) == 42);

  // Parsing should stop at the first non-integer delimiter.
  CHECK(in.peek() == ',');
}

TEST_CASE("get_int returns zero without consuming invalid input") {
  std::istringstream in("alpha");

  // Invalid integer input should leave the stream untouched.
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

    // The parser should accept each supported numeric spelling.
    CHECK(Loci::parse::is_real(in));
    CHECK(Loci::parse::get_real(in) == doctest::Approx(tc.value));

    // Any trailing non-numeric suffix should remain available to the caller.
    CHECK(read_remaining(in) == tc.rest);
  }
}

TEST_CASE("get_real ignores invalid exponents and leaves the exponent marker in the stream") {
  std::istringstream in("12ex");

  // Parsing should stop before the invalid exponent marker.
  CHECK(Loci::parse::get_real(in) == doctest::Approx(12.0));
  CHECK(in.peek() == 'e');
  CHECK(read_remaining(in) == "ex");
}

TEST_CASE("get_real returns zero without consuming invalid non-numeric input") {
  std::istringstream in("name");

  // Non-numeric input should be rejected without consuming the leading text.
  CHECK_FALSE(Loci::parse::is_real(in));
  CHECK(Loci::parse::get_real(in) == doctest::Approx(0.0));
  CHECK(in.peek() == 'n');
}

TEST_CASE("is_string and get_string parse quoted text including slashes and spaces") {
  std::istringstream in(" // comment\n\"a // b c\" tail");

  // Quoted strings should preserve internal spaces and slash sequences.
  CHECK(Loci::parse::is_string(in));
  CHECK(Loci::parse::get_string(in) == "a // b c");

  // After consuming the string, the trailing identifier should still parse.
  Loci::parse::kill_white_space(in);
  CHECK(Loci::parse::get_name(in) == "tail");
}

TEST_CASE("get_string supports empty quoted strings and leaves following text") {
  std::istringstream in("\"\"next");

  // Empty quoted strings should parse successfully and leave following text.
  CHECK(Loci::parse::get_string(in).empty());
  CHECK(Loci::parse::get_name(in) == "next");
}

TEST_CASE("get_string returns empty string without consuming non-string input") {
  std::istringstream in("name");

  // Non-string input should remain untouched when no opening quote is present.
  CHECK_FALSE(Loci::parse::is_string(in));
  CHECK(Loci::parse::get_string(in).empty());
  CHECK(in.peek() == 'n');
}

TEST_CASE("get_string returns collected text when the closing quote is missing") {
  std::istringstream in("\"unterminated");

  // Unterminated strings currently return the collected payload and consume the stream.
  CHECK(Loci::parse::get_string(in) == "unterminated");
  CHECK(in.eof());
}

TEST_CASE("is_token checks for a token without consuming it") {
  std::istringstream in("  &&rest");

  // Token probing should not consume the token on success.
  CHECK(Loci::parse::is_token(in, "&&"));
  CHECK(in.peek() == '&');

  // A subsequent `get_token()` should then consume it normally.
  CHECK(Loci::parse::get_token(in, "&&"));
  CHECK(Loci::parse::get_name(in) == "rest");
}

TEST_CASE("get_token skips leading whitespace and comments before matching") {
  std::istringstream in(" // comment\n::name");

  // Token matching should honor the same whitespace/comment skipping rules.
  CHECK(Loci::parse::get_token(in, "::"));
  CHECK(Loci::parse::get_name(in) == "name");
}

TEST_CASE("is_token and get_token restore the stream after a partial mismatch") {
  {
    std::istringstream in("=>");

    // A failed probe should restore the stream to its original position.
    CHECK_FALSE(Loci::parse::is_token(in, "=="));
    CHECK(read_remaining(in) == "=>");
  }

  {
    std::istringstream in("=>");

    // A failed consuming match should likewise leave the input untouched.
    CHECK_FALSE(Loci::parse::get_token(in, "=="));
    CHECK(read_remaining(in) == "=>");
  }
}
