#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/variable.h>

#ifndef LOCI_ENABLE_VARIABLE_KNOWN_BUG_TESTS
#define LOCI_ENABLE_VARIABLE_KNOWN_BUG_TESTS 0
#endif

namespace {
Loci::variable parse_variable(const std::string &text) {
  return Loci::variable(text);
}

std::string time_text(const Loci::time_ident &t) {
  std::ostringstream out;
  out << t;
  return out.str();
}

bool contains(const std::string &text, const std::string &token) {
  return text.find(token) != std::string::npos;
}

template <typename Fn>
std::string capture_cerr(Fn fn) {
  std::ostringstream err;
  std::streambuf *old = std::cerr.rdbuf(err.rdbuf());
  try {
    fn();
  } catch(...) {
    std::cerr.rdbuf(old);
    throw;
  }
  std::cerr.rdbuf(old);
  return err.str();
}
} // namespace

//----------------------------------------------------------------------------
// Potential Bugs
//----------------------------------------------------------------------------
//
// These are cases that appear to demonstrate behavior that is not intended.
//
// Enable them with LOCI_ENABLE_VARIABLE_KNOWN_BUG_TESTS=1 when we want to
// revisit these behaviors.

// The time-list parser mutates partial state before it has validated the full
// list, so malformed input such as rho{n,1,m} can retain a broken partial time
// label instead of falling back to an empty time.
TEST_CASE("known bug: malformed time list should not retain partial parse state [known-bug]" *
          doctest::skip(LOCI_ENABLE_VARIABLE_KNOWN_BUG_TESTS == 0)) {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho{n,1,m}"); });
  CHECK(contains(err, "syntax error interpreting time list"));
  CHECK(v.time() == Loci::time_ident());
}

// drop_priority() assumes a priority scope exists and unconditionally pops from
// the vector. On an unscoped variable this should be a no-op, but the current
// implementation can crash instead.
TEST_CASE("known bug: drop_priority on an unscoped variable should be a no-op [known-bug]" *
          doctest::skip(LOCI_ENABLE_VARIABLE_KNOWN_BUG_TESTS == 0)) {
  Loci::variable v = parse_variable("rho");
  CHECK(v.get_info().priority.empty());
  Loci::variable dropped = v.drop_priority();
  CHECK(dropped == v);
}

// variable::info::operator= copies most fields but currently omits time_id, so
// assigning metadata from a timed variable silently loses its time label.
TEST_CASE("known bug: variable::info assignment should copy time_id [known-bug]" *
          doctest::skip(LOCI_ENABLE_VARIABLE_KNOWN_BUG_TESTS == 0)) {
  Loci::variable::info src = parse_variable("rho{n,it}").get_info();
  Loci::variable::info dst;
  dst = src;
  CHECK(dst.time_id == src.time_id);
}

TEST_CASE("time_ident hierarchy supports parent/child traversal and before ordering") {
  Loci::time_ident t_n = parse_variable("rho{n}").time();
  const Loci::time_ident t_nit = parse_variable("rho{n,it}").time();
  const Loci::time_ident t_nrk = parse_variable("rho{n,rk}").time();

  CHECK(time_text(t_n) == "n");
  CHECK(time_text(t_nit) == "n,it");
  CHECK(time_text(t_nrk) == "n,rk");
  CHECK(t_n.before(t_nit));
  CHECK_FALSE(t_nit.before(t_n));
  CHECK(t_nit.parent() == t_n);

  std::vector<Loci::time_ident> children = t_n.children();
  bool saw_it = false;
  bool saw_rk = false;
  for(size_t i = 0; i < children.size(); ++i) {
    if(children[i].level_name() == "it") {
      saw_it = true;
    } else if(children[i].level_name() == "rk") {
      saw_rk = true;
    }
  }
  CHECK(saw_it);
  CHECK(saw_rk);
}

TEST_CASE("time_ident append/prepend constructors and prepend_time helper compose paths") {
  const Loci::time_ident stationary;
  const Loci::time_ident t_n("n", stationary);
  const Loci::time_ident t_nit(t_n, "it");
  const Loci::time_ident t_root_nit("root", t_nit);
  const Loci::time_ident t_p = parse_variable("q{p}").time();
  const Loci::time_ident t_pnit = Loci::prepend_time(t_p, t_nit);

  CHECK(time_text(t_n) == "n");
  CHECK(time_text(t_nit) == "n,it");
  CHECK(time_text(t_root_nit) == "root,n,it");
  CHECK(time_text(t_pnit) == "p,n,it");
}

TEST_CASE("variable parser handles plain names, priority scope, namespace and $ variables") {
  CHECK(parse_variable("rho").str() == "rho");
  CHECK(parse_variable("i::j::rho").str() == "i::j::rho");
  CHECK(parse_variable("ns1@ns2@rho").str() == "ns1@ns2@rho");

  const Loci::variable tv = parse_variable("$n");
  CHECK(tv.is_time_variable());
  CHECK(tv.str() == "$n");
}

TEST_CASE("variable parser handles time offsets and assignment forms in braces") {
  const Loci::variable plus = parse_variable("rho{n+2}");
  const Loci::variable minus = parse_variable("rho{n-2}");
  const Loci::variable assign = parse_variable("rho{n=3}");
  const Loci::variable assign_neg = parse_variable("rho{n=-2}");

  CHECK(time_text(plus.time()) == "n");
  CHECK(plus.get_info().offset == 2);
  CHECK_FALSE(plus.get_info().assign);
  CHECK(plus.str() == "rho{n+2}");

  CHECK(minus.get_info().offset == -2);
  CHECK_FALSE(minus.get_info().assign);
  CHECK(minus.str() == "rho{n-2}");

  CHECK(assign.get_info().offset == 3);
  CHECK(assign.get_info().assign);
  CHECK(assign.str() == "rho{n=3}");

  CHECK(assign_neg.get_info().offset == -2);
  CHECK(assign_neg.get_info().assign);
  CHECK(assign_neg.str() == "rho{n=-2}");
}

TEST_CASE("function and function-brace variable forms parse arguments and time labels") {
  Loci::variable vf = parse_variable("f(x,y)");
  REQUIRE(vf.get_arg_list().size() == 2);
  CHECK(Loci::variable(vf.get_arg_list()[0]).str() == "x");
  CHECK(Loci::variable(vf.get_arg_list()[1]).str() == "y");
  CHECK(vf.str() == "f(x,y)");

  Loci::variable vfb = parse_variable("f(x){n+1}");
  REQUIRE(vfb.get_arg_list().size() == 1);
  CHECK(Loci::variable(vfb.get_arg_list()[0]).str() == "x");
  CHECK(time_text(vfb.time()) == "n");
  CHECK(vfb.get_info().offset == 1);
  CHECK_FALSE(vfb.get_info().assign);
  CHECK(vfb.str() == "f(x){n+1}");
}

TEST_CASE("variable constructors rewrite argument list and time labels") {
  const Loci::variable vf = parse_variable("f(x,y)");
  const Loci::variable x = parse_variable("x");
  const Loci::variable y = parse_variable("y");
  const Loci::variable swapped(vf, std::vector<int>{y.ident(), x.ident()});
  CHECK(swapped.str() == "f(y,x)");

  const Loci::time_ident t_it = parse_variable("q{it}").time();
  const Loci::variable rho_n2 = parse_variable("rho{n+2}");
  const Loci::variable rho_it(rho_n2, t_it);
  CHECK(rho_it.str() == "rho{it}");

  const Loci::variable rho_nit(parse_variable("q{n}").time(), parse_variable("rho{it}"));
  CHECK(rho_nit.str() == "rho{n,it}");

  const Loci::variable tvar_from_time(t_it);
  CHECK(tvar_from_time.is_time_variable());
  CHECK(tvar_from_time.time() == t_it);
}

TEST_CASE("variable transform helpers modify only the requested metadata") {
  const Loci::variable timed = parse_variable("rho{n,it}");
  CHECK(timed.parent().str() == "rho{n}");

  const Loci::variable assigned = parse_variable("rho{n=3}");
  CHECK(assigned.drop_assign().str() == "rho{n+3}");
  CHECK(assigned.new_offset(-2).str() == "rho{n=-2}");

  const Loci::variable scoped = parse_variable("a::b::rho");
  CHECK(scoped.drop_priority().str() == "b::rho");
  CHECK(scoped.drop_priority().drop_priority().str() == "rho");

  const Loci::variable namespaced = parse_variable("x@y@rho");
  CHECK(namespaced.drop_namespace().str() == "y@rho");
  CHECK(namespaced.drop_namespace().drop_namespace().str() == "rho");

  const Loci::variable f = parse_variable("f(x,g(y))");
  const Loci::variable f_ns = f.add_namespace("chem");
  CHECK(f_ns.str() == "chem@f(chem@x,chem@g(chem@y))");
}

TEST_CASE("variableSet constructors and print are deterministic and lexicographically sorted") {
  {
    const Loci::variableSet empty;
    std::ostringstream out;
    out << empty;
    CHECK(out.str() == "()");
  }

  {
    const Loci::variableSet one(Loci::expression::create("a"));
    std::ostringstream out;
    out << one;
    CHECK(out.str() == "a");
  }

  {
    const Loci::variableSet many(Loci::expression::create("c,a,b"));
    CHECK(many.inSet(parse_variable("a")));
    CHECK(many.inSet(parse_variable("b")));
    CHECK(many.inSet(parse_variable("c")));

    const std::vector<std::string> sorted = many.lexico_sort();
    REQUIRE(sorted.size() == 3);
    CHECK(sorted[0] == "a");
    CHECK(sorted[1] == "b");
    CHECK(sorted[2] == "c");

    std::ostringstream out;
    out << many;
    CHECK(out.str() == "(a,b,c)");
  }
}

TEST_CASE("vmap_info parses mapping chain, output vars and assignment pairs") {
  const Loci::vmap_info vm(Loci::expression::create("a->b->(x,y=z)"));

  REQUIRE(vm.mapping.size() == 2);
  CHECK(vm.mapping[0].inSet(parse_variable("a")));
  CHECK(vm.mapping[1].inSet(parse_variable("b")));

  CHECK(vm.var.inSet(parse_variable("x")));
  CHECK(vm.var.inSet(parse_variable("y")));

  REQUIRE(vm.assign.size() == 1);
  CHECK(vm.assign[0].first == parse_variable("y"));
  CHECK(vm.assign[0].second == parse_variable("z"));

  std::ostringstream out;
  out << vm;
  CHECK(out.str() == "a->b->(x,y)");
}

TEST_CASE("failure mode: non-integer time offset logs error and keeps zero offset") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho{n+foo}"); });

  CHECK(contains(err, "time level should be of form NAME <op> int"));
  CHECK(time_text(v.time()) == "n");
  CHECK(v.get_info().offset == 0);
  CHECK_FALSE(v.get_info().assign);
  CHECK(v.str() == "rho{n}");
}

TEST_CASE("failure mode: malformed n+1+2 time offset logs errors and drops time label") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho{n+1+2}"); });

  CHECK(contains(err, "syntax error in time label"));
  CHECK(v.time() == Loci::time_ident());
  CHECK(v.get_info().offset == 0);
  CHECK_FALSE(v.get_info().assign);
  CHECK(v.str() == "rho");
}

TEST_CASE("failure mode: non-name time level head logs errors and drops time label") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho{1+2}"); });

  CHECK(contains(err, "time level should be of form NAME <op> int"));
  CHECK(v.time() == Loci::time_ident());
  CHECK(v.str() == "rho");
}

TEST_CASE("failure mode: invalid namespace head logs error and falls back to parsed tail name") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("1@rho"); });

  CHECK(contains(err, "unable to interpret namespace list in expression"));
  CHECK(v.str() == "rho");
}

TEST_CASE("failure mode: non-variable expressions log error and map to *NONAME* sentinel") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho+1"); });

  CHECK(contains(err, "unable to interpret expression"));
  CHECK(v.get_info().name == "*NONAME*");
  CHECK(v.str() == "*NONAME*");
}

TEST_CASE("failure mode: unterminated name-brace expression is tolerated as plain name") {
  CHECK_NOTHROW(parse_variable("rho{"));
  const Loci::variable v = parse_variable("rho{");
  CHECK(v.time() == Loci::time_ident());
  CHECK(v.str() == "rho");
}
