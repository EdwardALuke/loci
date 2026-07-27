#include <sstream>
#include <string>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/parse.h>

/// A typical input line should skip leading comments and read each value and
/// separator in order.
TEST_CASE("names, numbers, strings, and punctuation parse in order") {
  std::istringstream input(
    " // settings\n"
    "field = -1.25e2, count = +42, label = \"wall face\"");

  CHECK(Loci::parse::get_name(input) == "field");
  CHECK(Loci::parse::get_token(input, "="));
  CHECK(Loci::parse::get_real(input) == doctest::Approx(-125.0));
  CHECK(Loci::parse::get_token(input, ","));

  CHECK(Loci::parse::get_name(input) == "count");
  CHECK(Loci::parse::get_token(input, "="));
  CHECK(Loci::parse::get_int(input) == 42);
  CHECK(Loci::parse::get_token(input, ","));

  CHECK(Loci::parse::get_name(input) == "label");
  CHECK(Loci::parse::get_token(input, "="));
  CHECK(Loci::parse::get_string(input) == "wall face");
  CHECK(input.peek() == std::char_traits<char>::eof());
}

/// Checking an operator must leave it available, and rejecting a longer
/// operator must preserve the shorter valid one.
TEST_CASE("operator checks leave the next valid operator available") {
  std::istringstream checked_parenthesis(" // comment\n(name");

  CHECK(Loci::parse::is_token(checked_parenthesis, "("));
  CHECK(Loci::parse::get_token(checked_parenthesis, "("));
  CHECK(Loci::parse::get_name(checked_parenthesis) == "name");

  std::istringstream checked_operator(">value");
  CHECK_FALSE(Loci::parse::is_token(checked_operator, ">="));
  CHECK(Loci::parse::get_token(checked_operator, ">"));
  CHECK(Loci::parse::get_name(checked_operator) == "value");

  std::istringstream rejected_operator(">value");
  CHECK_FALSE(Loci::parse::get_token(rejected_operator, ">="));
  CHECK(Loci::parse::get_token(rejected_operator, ">"));
  CHECK(Loci::parse::get_name(rejected_operator) == "value");
}
