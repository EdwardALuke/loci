#include <map>
#include <string>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/expr.h>

TEST_CASE("expression evaluate handles unary operators with double variables") {
  std::map<std::string, double> vars;
  vars["y"] = 0.00075;

  Loci::exprP expr = Loci::expression::create("-y-0.00025");
  CHECK(expr->evaluate(vars) == doctest::Approx(-0.001));

  Loci::exprP plus_expr = Loci::expression::create("+y");
  CHECK(plus_expr->evaluate(vars) == doctest::Approx(0.00075));
}

TEST_CASE("expression evaluate handles unary operators with int variables") {
  std::map<std::string, int> vars;
  vars["n"] = 7;

  Loci::exprP expr = Loci::expression::create("-n+2");
  CHECK(expr->evaluate(vars) == -5);

  Loci::exprP plus_expr = Loci::expression::create("+n");
  CHECK(plus_expr->evaluate(vars) == 7);
}
