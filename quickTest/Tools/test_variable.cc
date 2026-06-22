#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/variable.h>

namespace {
/// @brief Parses a variable expression through the public `variable` API.
Loci::variable parse_variable(const std::string &text) {
  return Loci::variable(text);
}

/// @brief Renders a `time_ident` into the user-facing textual form.
std::string time_text(const Loci::time_ident &t) {
  std::ostringstream out;
  out << t;
  return out.str();
}

/// @brief Returns true when `token` appears anywhere inside `text`.
bool contains(const std::string &text, const std::string &token) {
  return text.find(token) != std::string::npos;
}

template <typename Fn>
/// @brief Captures `std::cerr` while `fn` runs so tests can assert on
///   diagnostics.
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
// Variable Expression Behavior
//----------------------------------------------------------------------------
// These tests exercise the public variable APIs that rules and scheduler
// construction depend on: time labels, priority/namespace identity, parametric
// facts, mapped signatures, and deterministic variable-set behavior.

TEST_CASE("time_ident hierarchy supports parent traversal and before ordering") {
  Loci::time_ident t_n = parse_variable("rho{n}").time();
  const Loci::time_ident t_nit = parse_variable("rho{n,it}").time();
  const Loci::time_ident t_nrk = parse_variable("rho{n,rk}").time();

  // Time labels define loop identity.  `{n,it}` should sort after and report
  // `{n}` as its parent so promoted and loop-local facts remain distinct.
  CHECK(time_text(t_n) == "n");
  CHECK(time_text(t_nit) == "n,it");
  CHECK(time_text(t_nrk) == "n,rk");
  CHECK(t_n.before(t_nit));
  CHECK_FALSE(t_nit.before(t_n));
  CHECK(t_nit.parent() == t_n);
}

TEST_CASE("time_ident append/prepend constructors and prepend_time helper compose paths") {
  const Loci::time_ident stationary;
  const Loci::time_ident t_n("n", stationary);
  const Loci::time_ident t_nit(t_n, "it");
  const Loci::time_ident t_root_nit("root", t_nit);
  const Loci::time_ident t_p = parse_variable("q{p}").time();
  const Loci::time_ident t_pnit = Loci::prepend_time(t_p, t_nit);

  // Appending and prepending should build the paths used when scheduler-created
  // rules carry facts into nested iteration levels.
  CHECK(time_text(t_n) == "n");
  CHECK(time_text(t_nit) == "n,it");
  CHECK(time_text(t_root_nit) == "root,n,it");
  CHECK(time_text(t_pnit) == "p,n,it");
}

TEST_CASE("variable parser handles plain names, priority scope, namespace and $ variables") {
  // Plain facts, priority-qualified override facts, and namespaced module-local
  // facts should round trip through the public parser and printer.
  CHECK(parse_variable("rho").str() == "rho");
  CHECK(parse_variable("i::j::rho").str() == "i::j::rho");
  CHECK(parse_variable("ns1@ns2@rho").str() == "ns1@ns2@rho");

  // `$` variables carry loop indices rather than ordinary fact names.
  const Loci::variable tv = parse_variable("$n");
  CHECK(tv.is_time_variable());
  CHECK(tv.str() == "$n");
}

TEST_CASE("variable parser handles time offsets and assignment forms in braces") {
  const Loci::variable plus = parse_variable("rho{n+2}");
  const Loci::variable minus = parse_variable("rho{n-2}");
  const Loci::variable assign = parse_variable("rho{n=3}");
  const Loci::variable assign_neg = parse_variable("rho{n=-2}");

  // Relative and fixed time references should preserve the exact textual form
  // that rule signatures use for history and seed facts.
  CHECK(time_text(plus.time()) == "n");
  CHECK(plus.str() == "rho{n+2}");

  CHECK(time_text(minus.time()) == "n");
  CHECK(minus.str() == "rho{n-2}");

  CHECK(time_text(assign.time()) == "n");
  CHECK(assign.str() == "rho{n=3}");

  CHECK(time_text(assign_neg.time()) == "n");
  CHECK(assign_neg.str() == "rho{n=-2}");
}

TEST_CASE("function and function-brace variable forms parse arguments and time labels") {
  // Function-style variables represent parametric fact families.  The concrete
  // arguments should be available to callers and preserved in printed form.
  Loci::variable vf = parse_variable("f(x,y)");
  REQUIRE(vf.get_arg_list().size() == 2);
  CHECK(Loci::variable(vf.get_arg_list()[0]).str() == "x");
  CHECK(Loci::variable(vf.get_arg_list()[1]).str() == "y");
  CHECK(vf.str() == "f(x,y)");

  // Time-qualified parametric facts should preserve both the concrete argument
  // and loop identity.
  Loci::variable vfb = parse_variable("f(x){n+1}");
  REQUIRE(vfb.get_arg_list().size() == 1);
  CHECK(Loci::variable(vfb.get_arg_list()[0]).str() == "x");
  CHECK(time_text(vfb.time()) == "n");
  CHECK(vfb.str() == "f(x){n+1}");
}

TEST_CASE("variable constructors rewrite argument list and time labels") {
  const Loci::variable vf = parse_variable("f(x,y)");
  const Loci::variable x = parse_variable("x");
  const Loci::variable y = parse_variable("y");

  // Rebuilding a parametric fact with a new argument order should not change
  // the fact family name.
  const Loci::variable swapped(vf, std::vector<int>{y.ident(), x.ident()});
  CHECK(swapped.str() == "f(y,x)");

  const Loci::time_ident t_it = parse_variable("q{it}").time();
  const Loci::variable rho_n2 = parse_variable("rho{n+2}");

  // Replacing or prepending time labels should rewrite only loop identity, not
  // the base fact name.
  const Loci::variable rho_it(rho_n2, t_it);
  CHECK(rho_it.str() == "rho{it}");

  const Loci::variable rho_nit(parse_variable("q{n}").time(), parse_variable("rho{it}"));
  CHECK(rho_nit.str() == "rho{n,it}");

  // Constructing from a `time_ident` alone yields the `$` variable used to
  // represent loop indices.
  const Loci::variable tvar_from_time(t_it);
  CHECK(tvar_from_time.is_time_variable());
  CHECK(tvar_from_time.time() == t_it);
}

TEST_CASE("variable transform helpers modify only the requested metadata") {
  const Loci::variable timed = parse_variable("rho{n,it}");

  // Hierarchy navigation should peel off only the last time level, preserving
  // the outer loop identity.
  CHECK(timed.parent().str() == "rho{n}");

  const Loci::variable assigned = parse_variable("rho{n=3}");

  // Assignment and offset helpers should preserve the same fact while changing
  // only the time-reference mode.
  CHECK(assigned.drop_assign().str() == "rho{n+3}");
  CHECK(assigned.new_offset(-2).str() == "rho{n=-2}");

  const Loci::variable scoped = parse_variable("a::b::rho");

  // Priority and namespace helpers should strip one override/module layer at a
  // time, leaving the underlying fact identity intact.
  CHECK(scoped.drop_priority().str() == "b::rho");
  CHECK(scoped.drop_priority().drop_priority().str() == "rho");

  const Loci::variable namespaced = parse_variable("x@y@rho");
  CHECK(namespaced.drop_namespace().str() == "y@rho");
  CHECK(namespaced.drop_namespace().drop_namespace().str() == "rho");

  const Loci::variable f = parse_variable("f(x,g(y))");
  const Loci::variable f_ns = f.add_namespace("chem");

  // Adding a namespace should propagate through the full nested parametric
  // expression so every referenced fact moves into the module-local scope.
  CHECK(f_ns.str() == "chem@f(chem@x,chem@g(chem@y))");
}

TEST_CASE("variableSet constructors and print are deterministic and lexicographically sorted") {
  {
    // The empty set should print as an empty tuple-like form.
    const Loci::variableSet empty;
    std::ostringstream out;
    out << empty;
    CHECK(out.str() == "()");
  }

  {
    // A singleton set should print without tuple decoration.
    const Loci::variableSet one(Loci::expression::create("a"));
    std::ostringstream out;
    out << one;
    CHECK(out.str() == "a");
  }

  {
    // Multi-variable sets should sort deterministically regardless of input
    // order so diagnostics and generated rule descriptions are stable.
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

  // Mapping links should be recorded in traversal order; `a->b->x` means the
  // rule crosses those maps before reading or writing the payload facts.
  REQUIRE(vm.mapping.size() == 2);
  CHECK(vm.mapping[0].inSet(parse_variable("a")));
  CHECK(vm.mapping[1].inSet(parse_variable("b")));

  // The output variable set should include both plain and assigned outputs so a
  // rule can advertise every fact it derives.
  CHECK(vm.var.inSet(parse_variable("x")));
  CHECK(vm.var.inSet(parse_variable("y")));

  // Assignment pairs should retain the source and destination variables used by
  // aliasing target forms such as `y=z`.
  REQUIRE(vm.assign.size() == 1);
  CHECK(vm.assign[0].first == parse_variable("y"));
  CHECK(vm.assign[0].second == parse_variable("z"));

  // Printing should expose the scheduler-visible mapping chain and derived
  // outputs without leaking the internal alias source.
  std::ostringstream out;
  out << vm;
  CHECK(out.str() == "a->b->(x,y)");
}

TEST_CASE("failure mode: non-integer time offset logs error and keeps zero offset") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho{n+foo}"); });

  // A malformed relative offset should not invent a different loop level; the
  // recoverable base time label remains visible to callers.
  CHECK(contains(err, "time level should be of form NAME <op> int"));
  CHECK(time_text(v.time()) == "n");
  CHECK(v.str() == "rho{n}");
}

TEST_CASE("failure mode: malformed n+1+2 time offset logs errors and drops time label") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho{n+1+2}"); });

  // A malformed chained offset should not leave a phantom partial time label
  // attached to the fact.
  CHECK(contains(err, "syntax error in time label"));
  CHECK(v.time() == Loci::time_ident());
  CHECK(v.str() == "rho");
}

TEST_CASE("failure mode: non-name time level head logs errors and drops time label") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho{1+2}"); });

  // A non-name time head should be rejected and leave a plain variable behind.
  CHECK(contains(err, "time level should be of form NAME <op> int"));
  CHECK(v.time() == Loci::time_ident());
  CHECK(v.str() == "rho");
}

TEST_CASE("failure mode: invalid namespace head logs error and falls back to parsed tail name") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("1@rho"); });

  // Namespace parse failure should still preserve the recoverable tail name.
  CHECK(contains(err, "unable to interpret namespace list in expression"));
  CHECK(v.str() == "rho");
}

TEST_CASE("failure mode: non-variable expressions log error and map to *NONAME* sentinel") {
  Loci::variable v;
  const std::string err = capture_cerr([&]() { v = parse_variable("rho+1"); });

  // Non-variable expressions should report an error and fall back to the public
  // sentinel variable name.
  CHECK(contains(err, "unable to interpret expression"));
  CHECK(v.str() == "*NONAME*");
}

TEST_CASE("failure mode: unterminated name-brace expression is tolerated as plain name") {
  // This malformed input is currently tolerated by falling back to a plain
  // variable name rather than throwing.
  CHECK_NOTHROW(parse_variable("rho{"));
  const Loci::variable v = parse_variable("rho{");
  CHECK(v.time() == Loci::time_ident());
  CHECK(v.str() == "rho");
}
