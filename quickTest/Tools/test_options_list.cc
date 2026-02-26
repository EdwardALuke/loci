#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/options_list.h>

namespace {
Loci::options_list parse_options(const std::string &text) {
  Loci::options_list opts;
  std::istringstream in(text);
  in >> opts;
  return opts;
}

Loci::option_values parse_value(const std::string &text) {
  Loci::option_values v;
  std::istringstream in(text);
  in >> v;
  return v;
}
} // namespace

TEST_CASE("options_list parses mixed scalar, bool, name and string values") {
  Loci::options_list opts =
      parse_options("<flag,enabled=$false,num=3.5,name=alpha,msg=\"hello world\">");

  CHECK(opts.getOptionValueType("flag") == Loci::BOOLEAN);
  CHECK(opts.getOptionValueType("enabled") == Loci::BOOLEAN);
  CHECK(opts.getOptionValueType("num") == Loci::REAL);
  CHECK(opts.getOptionValueType("name") == Loci::NAME);
  CHECK(opts.getOptionValueType("msg") == Loci::STRING);

  bool flag = false;
  bool enabled = true;
  double num = 0.0;
  std::string name;
  std::string msg;

  opts.getOption("flag", flag);
  opts.getOption("enabled", enabled);
  opts.getOption("num", num);
  opts.getOption("name", name);
  opts.getOption("msg", msg);

  CHECK(flag == true);
  CHECK(enabled == false);
  CHECK(num == doctest::Approx(3.5));
  CHECK(name == "alpha");
  CHECK(msg == "hello world");
}

TEST_CASE("options_list parses function and list values") {
  Loci::options_list opts = parse_options(
      "<f=foo(1,bar,$true),vals=[1,2,3],cfg=[mode=fast,tag=\"x\",n=3]>");

  std::string fname;
  Loci::options_list::arg_list fargs;
  opts.getOption("f", fname, fargs);
  CHECK(fname == "foo");
  REQUIRE(fargs.size() == 3);
  CHECK(fargs[0].type_of() == Loci::REAL);
  CHECK(fargs[1].type_of() == Loci::NAME);
  CHECK(fargs[2].type_of() == Loci::BOOLEAN);

  Loci::options_list::arg_list vals;
  opts.getOption("vals", vals);
  REQUIRE(vals.size() == 3);
  CHECK(vals[0].type_of() == Loci::REAL);
  CHECK(vals[1].type_of() == Loci::REAL);
  CHECK(vals[2].type_of() == Loci::REAL);

  Loci::options_list::arg_list cfg;
  opts.getOption("cfg", cfg);
  REQUIRE(cfg.size() == 3);
  CHECK(cfg[0].type_of() == Loci::NAME_ASSIGN);
  CHECK(cfg[1].type_of() == Loci::NAME_ASSIGN);
  CHECK(cfg[2].type_of() == Loci::NAME_ASSIGN);
}

TEST_CASE("options_list stream parser rejects duplicate keys and malformed header") {
  {
    Loci::options_list opts;
    std::istringstream in("<a=1,a=2>");
    CHECK_THROWS(in >> opts);
  }
  {
    Loci::options_list opts;
    std::istringstream in("a=1>");
    CHECK_THROWS(in >> opts);
  }
}

TEST_CASE("options_list stream parser accepts implicit booleans and rejects bad separators") {
  {
    Loci::options_list opts = parse_options("<enabled,rate=2>");
    bool enabled = false;
    double rate = 0.0;
    opts.getOption("enabled", enabled);
    opts.getOption("rate", rate);
    CHECK(enabled == true);
    CHECK(rate == doctest::Approx(2.0));
  }
  {
    Loci::options_list opts;
    std::istringstream in("<a=1 b=2>");
    CHECK_THROWS(in >> opts);
  }
}

TEST_CASE("options_list enforces restricted option set for stream Input") {
  {
    Loci::options_list opts("a:b");
    std::istringstream in("<a=1,b=2>");
    CHECK_NOTHROW(in >> opts);
  }
  {
    Loci::options_list opts("a:b");
    std::istringstream in("<a=1,c=2>");
    CHECK_THROWS(in >> opts);
  }
  {
    Loci::options_list opts;
    std::istringstream in("<newopt=1>");
    CHECK_NOTHROW(in >> opts);
    CHECK(opts.optionExists("newopt"));
  }
}

TEST_CASE("options_list setOption/getOption and checkOption basics") {
  Loci::options_list opts;
  opts.setOption("flag", true);
  opts.setOption("num", 2.5);
  opts.setOption("mode", std::string("steady"));
  opts.setOption("vel", Loci::UNIT_type(Loci::UNIT_type::MKS, "general", 10.0, "m/s"));

  bool flag = false;
  double num = 0.0;
  std::string mode;
  double vel = 0.0;

  opts.getOption("flag", flag);
  opts.getOption("num", num);
  opts.getOption("mode", mode);
  opts.getOptionUnits("vel", "m/s", vel);

  CHECK(flag == true);
  CHECK(num == doctest::Approx(2.5));
  CHECK(mode == "steady");
  CHECK(vel == doctest::Approx(10.0));

  CHECK(opts.checkOption("mode", "steady"));
  CHECK_FALSE(opts.checkOption("mode", "transient"));
  CHECK_FALSE(opts.optionExists("missing"));
}

TEST_CASE("options_list Print/Input roundtrip keeps parseable semantics") {
  Loci::options_list orig = parse_options("<a=1,b=$false,c=name,d=[1,2],e=f(3),u=2 m/s>");

  std::ostringstream out;
  out << orig;

  Loci::options_list copy;
  std::istringstream in(out.str());
  in >> copy;

  double a = 0.0;
  bool b = true;
  std::string c;
  double u = 0.0;
  std::string fname;
  Loci::options_list::arg_list fargs;
  Loci::options_list::arg_list dvals;

  copy.getOption("a", a);
  copy.getOption("b", b);
  copy.getOption("c", c);
  copy.getOptionUnits("u", "m/s", u);
  copy.getOption("e", fname, fargs);
  copy.getOption("d", dvals);

  CHECK(a == doctest::Approx(1.0));
  CHECK(b == false);
  CHECK(c == "name");
  CHECK(u == doctest::Approx(2.0));
  CHECK(fname == "f");
  REQUIRE(fargs.size() == 1);
  CHECK(fargs[0].type_of() == Loci::REAL);
  REQUIRE(dvals.size() == 2);
  CHECK(dvals[0].type_of() == Loci::REAL);
  CHECK(dvals[1].type_of() == Loci::REAL);
}

TEST_CASE("option_values parser handles derivative syntax with and without lists") {
  {
    Loci::option_values v = parse_value("3^2^^4 m/s");
    REQUIRE(v.type_of() == Loci::UNIT_VALUE);

    Loci::UNIT_type u;
    v.get_value(u);
    double value = 0.0;
    double grad = 0.0;
    double grad2 = 0.0;
    std::vector<double> grads;
    u.get_values_in("m/s", value, grad, grad2, grads);

    CHECK(value == doctest::Approx(3.0));
    CHECK(grad == doctest::Approx(2.0));
    CHECK(grad2 == doctest::Approx(4.0));
    CHECK(grads.empty());
  }
  {
    Loci::option_values v = parse_value("5^[2,3,4] m");
    REQUIRE(v.type_of() == Loci::UNIT_VALUE);

    Loci::UNIT_type u;
    v.get_value(u);
    double value = 0.0;
    double grad = 0.0;
    double grad2 = 0.0;
    std::vector<double> grads;
    u.get_values_in("m", value, grad, grad2, grads);

    CHECK(value == doctest::Approx(5.0));
    CHECK(grad == doctest::Approx(2.0));
    CHECK(grad2 == doctest::Approx(0.0));
    REQUIRE(grads.size() == 2);
    CHECK(grads[0] == doctest::Approx(3.0));
    CHECK(grads[1] == doctest::Approx(4.0));
  }
}

TEST_CASE("option_values parser accepts broad valid grammar corpus") {
  struct ValidCase {
    const char *text;
    Loci::option_value_type type;
  };

  const std::vector<ValidCase> cases = {
      {"0", Loci::REAL},
      {"-1.25e2", Loci::REAL},
      {"+.5", Loci::REAL},
      {"2 m/s", Loci::UNIT_VALUE},
      {"5 one/s", Loci::UNIT_VALUE},
      {"1 one/(Pa*s)", Loci::UNIT_VALUE},
      {"name_1", Loci::NAME},
      {"\"quoted text\"", Loci::STRING},
      {"$true", Loci::BOOLEAN},
      {"$false", Loci::BOOLEAN},
      {"func(1,name,$false)", Loci::FUNCTION},
      {"id=7", Loci::NAME_ASSIGN},
      {"[1,name,$true,k=4,\"s\"=5]", Loci::LIST},
  };

  for(const auto &tc : cases) {
    CAPTURE(tc.text);
    Loci::option_values v;
    std::istringstream in(tc.text);
    CHECK_NOTHROW(in >> v);
    CHECK(v.type_of() == tc.type);
  }
}

TEST_CASE("option_values parser rejects malformed grammar corpus") {
  const std::vector<std::string> invalid_cases = {
      "$maybe",
      "[1=2]",
      "1^[x] m",
      "1^[2,,3] m",
      "1 m/",
      "1 one/",
      "1 (m/s",
      "foo(>",
      "name=",
      "[a=]",
  };

  for(const auto &text : invalid_cases) {
    CAPTURE(text);
    Loci::option_values v;
    std::istringstream in(text);
    CHECK_THROWS(in >> v);
  }
}

TEST_CASE("option_values parser rejects unterminated quoted strings") {
  Loci::option_values v;
  std::istringstream in("\"datafile.dat");
  CHECK_THROWS(in >> v);
}

TEST_CASE("options_list parses quoted filename values and roundtrips them") {
  Loci::options_list opts =
      parse_options("<myVar=\"datafile.dat\",path=\"inputs/run_01/input-file.dat\">");

  std::string my_var;
  std::string path;
  opts.getOption("myVar", my_var);
  opts.getOption("path", path);
  CHECK(my_var == "datafile.dat");
  CHECK(path == "inputs/run_01/input-file.dat");

  std::ostringstream out;
  out << opts;
  Loci::options_list copy;
  std::istringstream in(out.str());
  in >> copy;

  std::string my_var_copy;
  std::string path_copy;
  copy.getOption("myVar", my_var_copy);
  copy.getOption("path", path_copy);
  CHECK(my_var_copy == my_var);
  CHECK(path_copy == path);
}

TEST_CASE("nested list NAME_ASSIGN roundtrips as nested structure") {
  Loci::option_values v = parse_value("[a=[1,2]]");
  REQUIRE(v.type_of() == Loci::LIST);

  std::ostringstream out;
  CHECK_NOTHROW(out << v);
  CHECK(out.str() == "[a=[1,2]]");

  Loci::option_values reparsed = parse_value(out.str());
  std::ostringstream out2;
  out2 << reparsed;
  CHECK(out2.str() == out.str());
}

TEST_CASE("options_list Input(arg_list) handles NAME and NAME_ASSIGN") {
  Loci::option_values list_value = parse_value("[a=1,b,c=3]");
  REQUIRE(list_value.type_of() == Loci::LIST);

  Loci::options_list::arg_list l;
  list_value.get_value(l);

  Loci::options_list opts;
  opts.Input(l);

  double a = 0.0, c = 0.0;
  bool b = false;
  opts.getOption("a", a);
  opts.getOption("c", c);
  opts.getOption("b", b);

  CHECK(a == doctest::Approx(1.0));
  CHECK(c == doctest::Approx(3.0));
  CHECK(b == true);
}

TEST_CASE("options_list Input(arg_list) rejects duplicates and ignores unknown when restricted") {
  {
    Loci::option_values list_value = parse_value("[a=1,a=2]");
    Loci::options_list::arg_list l;
    list_value.get_value(l);
    Loci::options_list opts;
    CHECK_THROWS(opts.Input(l));
  }
  {
    Loci::option_values list_value = parse_value("[a=1,b=2]");
    Loci::options_list::arg_list l;
    list_value.get_value(l);
    Loci::options_list opts("a");
    CHECK_NOTHROW(opts.Input(l));
    CHECK(opts.optionExists("a"));
    CHECK_FALSE(opts.optionExists("b"));
  }
}

TEST_CASE("options_list Input(arg_list) rejects non-name entries") {
  Loci::option_values list_value = parse_value("[1]");
  REQUIRE(list_value.type_of() == Loci::LIST);

  Loci::options_list::arg_list l;
  list_value.get_value(l);
  Loci::options_list opts;
  CHECK_THROWS(opts.Input(l));
}

TEST_CASE("options_list stream Input failure is non-atomic and may partially update values") {
  Loci::options_list opts = parse_options("<a=1,b=2>");

  std::istringstream bad("<a=3 b=4>");
  CHECK_THROWS(bad >> opts);

  double b = 0.0;
  opts.getOption("b", b);

  CHECK(opts.getOptionValueType("a") == Loci::UNIT_VALUE);
  double a_barns = 0.0;
  opts.getOptionUnits("a", "b", a_barns);
  CHECK(a_barns == doctest::Approx(3.0));
  CHECK(b == doctest::Approx(2.0));
}

TEST_CASE("restricted option constructor handles delimiter edge cases") {
  {
    Loci::options_list opts("a::b:");
    std::istringstream in("<a=1,b=2>");
    CHECK_NOTHROW(in >> opts);
  }
  {
    Loci::options_list opts(":");
    std::istringstream in("<a=1>");
    CHECK_THROWS(in >> opts);
  }
  {
    Loci::options_list opts(" a :b ");
    std::istringstream in("<a=1>");
    CHECK_THROWS(in >> opts);
  }
  {
    Loci::options_list opts("");
    std::istringstream in("<a=1>");
    CHECK_NOTHROW(in >> opts);
  }
}

TEST_CASE("options_list parser accepts broad valid format corpus") {
  struct ValidCase {
    const char *text;
    size_t expected_count;
  };

  const std::vector<ValidCase> cases = {
      {"<a=1>", 1},
      {"<a,b=$false,c=name>", 3},
      {"<a=1,// comment\nb=2>", 2},
      {"<a=[1,2,3],f=foo(1,$true),s=\"x y\">", 3},
      {"<u=5 one/s,v=1 one/(Pa*s)>", 2},
      {"<spaced = 1 , tabs = [ 1 , 2 ] >", 2},
  };

  for(const auto &tc : cases) {
    CAPTURE(tc.text);
    Loci::options_list opts;
    std::istringstream in(tc.text);
    CHECK_NOTHROW(in >> opts);
    CHECK(opts.getOptionNameList().size() == tc.expected_count);
  }
}

TEST_CASE("options_list parser rejects malformed format corpus") {
  const std::vector<std::string> invalid_cases = {
      "<",
      "<=1>",
      "<1=2>",
      "<a==1>",
      "<a=,b=2>",
      "<a b=2>",
      "<a=1;;b=2>",
      "<a=1,,b=2>",
      "<a=foo(1,2>",
      "<a=[1,2>",
  };

  for(const auto &text : invalid_cases) {
    CAPTURE(text);
    Loci::options_list opts;
    std::istringstream in(text);
    CHECK_THROWS(in >> opts);
  }
}

TEST_CASE("options_list getOption type mismatch behavior is non-throw and leaves outputs defaulted") {
  Loci::options_list opts = parse_options("<flag=$true,num=2,name=alpha,lst=[1,2],fn=f(3),u=5 m>");

  bool as_bool = false;
  CHECK_NOTHROW(opts.getOption("num", as_bool));
  CHECK(as_bool == false);

  double as_real = -99.0;
  CHECK_NOTHROW(opts.getOption("flag", as_real));
  CHECK(as_real == doctest::Approx(-99.0));

  Loci::UNIT_type as_units(Loci::UNIT_type::MKS, "general", 9.0, "m");
  CHECK_NOTHROW(opts.getOption("num", as_units));
  CHECK(as_units.get_value_in("m") == doctest::Approx(9.0));

  std::string as_name = "seed";
  CHECK_NOTHROW(opts.getOption("num", as_name));
  CHECK(as_name.empty());

  Loci::options_list::arg_list as_list;
  as_list.push_back(parse_value("9"));
  CHECK_NOTHROW(opts.getOption("num", as_list));
  CHECK(as_list.empty());

  std::string func_name = "seed";
  Loci::options_list::arg_list func_args;
  func_args.push_back(parse_value("9"));
  CHECK_NOTHROW(opts.getOption("num", func_name, func_args));
  CHECK(func_name.empty());
  CHECK(func_args.empty());
}

TEST_CASE("options_list stream Input updates existing option values on later parse") {
  Loci::options_list opts;

  {
    std::istringstream in("<a=1,b=$false>");
    CHECK_NOTHROW(in >> opts);
  }
  {
    std::istringstream in("<a=4,b>");
    CHECK_NOTHROW(in >> opts);
  }

  double a = 0.0;
  bool b = false;
  opts.getOption("a", a);
  opts.getOption("b", b);
  CHECK(a == doctest::Approx(4.0));
  CHECK(b == true);
}

TEST_CASE("minimal forms parse and roundtrip: empty options, empty list and empty function") {
  {
    Loci::options_list empty_opts = parse_options("<>");
    CHECK(empty_opts.getOptionNameList().empty());

    std::ostringstream out;
    out << empty_opts;
    CHECK(out.str() == "<>");
  }
  {
    Loci::option_values empty_list = parse_value("[]");
    REQUIRE(empty_list.type_of() == Loci::LIST);
    Loci::options_list::arg_list values;
    empty_list.get_value(values);
    CHECK(values.empty());
  }
  {
    Loci::option_values empty_fn = parse_value("foo()");
    REQUIRE(empty_fn.type_of() == Loci::FUNCTION);
    std::string fn_name;
    Loci::options_list::arg_list fn_args;
    empty_fn.get_value(fn_name);
    empty_fn.get_value(fn_args);
    CHECK(fn_name == "foo");
    CHECK(fn_args.empty());
  }
}

TEST_CASE("options_list vector getOptionUnits supports scalar, list and polar forms") {
  {
    Loci::options_list opts = parse_options("<v=2>");
    Loci::vector3d<double> vec;
    opts.getOptionUnits("v", "m/s", vec, 3.0);
    CHECK(vec.x == doctest::Approx(6.0));
    CHECK(vec.y == doctest::Approx(0.0));
    CHECK(vec.z == doctest::Approx(0.0));
  }
  {
    Loci::options_list opts = parse_options("<v=[1 m,2 m,3 m]>");
    Loci::vector3d<double> vec;
    opts.getOptionUnits("v", "cm", vec);
    CHECK(vec.x == doctest::Approx(100.0));
    CHECK(vec.y == doctest::Approx(200.0));
    CHECK(vec.z == doctest::Approx(300.0));
  }
  {
    Loci::options_list opts = parse_options("<v=polar(2 m,90 deg,0 deg)>");
    Loci::vector3d<double> vec;
    opts.getOptionUnits("v", "m", vec);
    CHECK(vec.x == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(vec.y == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(vec.z == doctest::Approx(0.0).epsilon(1e-12));
  }
  {
    Loci::options_list opts = parse_options("<v=2 m>");
    Loci::vector3d<double> vec;
    opts.getOptionUnits("v", "cm", vec);
    CHECK(vec.x == doctest::Approx(200.0));
    CHECK(vec.y == doctest::Approx(0.0));
    CHECK(vec.z == doctest::Approx(0.0));
  }
  {
    Loci::options_list opts = parse_options("<v=[1,2 m,3]>");
    Loci::vector3d<double> vec;
    opts.getOptionUnits("v", "m", vec, 10.0);
    CHECK(vec.x == doctest::Approx(10.0));
    CHECK(vec.y == doctest::Approx(2.0));
    CHECK(vec.z == doctest::Approx(30.0));
  }
}

TEST_CASE("options_list AD getOptionUnits overloads work for scalar and vector paths") {
  {
    Loci::options_list opts = parse_options("<q=10^[2,3,4] m/s>");

    Loci::FAD2d qd;
    opts.getOptionUnits("q", "km/h", qd);
    CHECK(qd.value == doctest::Approx(36.0));
    CHECK(qd.grad == doctest::Approx(7.2));
    CHECK(qd.grad2 == doctest::Approx(0.0));

    Loci::VFAD qv;
    opts.getOptionUnits("q", "km/h", qv, 0);
    CHECK(qv.data.value == doctest::Approx(36.0));
    CHECK(qv.data.grad[0] == doctest::Approx(7.2));
    CHECK(qv.data.grad[1] == doctest::Approx(10.8));
    CHECK(qv.data.grad[2] == doctest::Approx(14.4));
  }
  {
    Loci::options_list opts = parse_options("<vec=[1^[2],2^[3],3^[4]]>");

    Loci::vector3d<Loci::FAD2d> vec_d;
    opts.getOptionUnits("vec", "m", vec_d, Loci::FAD2d(10.0, 5.0, 0.0));
    CHECK(vec_d.x.value == doctest::Approx(10.0));
    CHECK(vec_d.y.value == doctest::Approx(20.0));
    CHECK(vec_d.z.value == doctest::Approx(30.0));
    CHECK(vec_d.x.grad == doctest::Approx(25.0));
    CHECK(vec_d.y.grad == doctest::Approx(40.0));
    CHECK(vec_d.z.grad == doctest::Approx(55.0));

    Loci::vector3d<Loci::FAD2d> vec_d_default;
    opts.getOptionUnits("vec", "m", vec_d_default);
    CHECK(vec_d_default.x.value == doctest::Approx(1.0));
    CHECK(vec_d_default.y.value == doctest::Approx(2.0));
    CHECK(vec_d_default.z.value == doctest::Approx(3.0));
    CHECK(vec_d_default.x.grad == doctest::Approx(2.0));
    CHECK(vec_d_default.y.grad == doctest::Approx(3.0));
    CHECK(vec_d_default.z.grad == doctest::Approx(4.0));

    Loci::vector3d<Loci::VFAD> vec_v;
    opts.getOptionUnits("vec", "m", vec_v, Loci::VFAD(10.0), 0);
    CHECK(vec_v.x.data.value == doctest::Approx(10.0));
    CHECK(vec_v.y.data.value == doctest::Approx(20.0));
    CHECK(vec_v.z.data.value == doctest::Approx(30.0));
    CHECK(vec_v.x.data.grad[0] == doctest::Approx(20.0));
    CHECK(vec_v.y.data.grad[0] == doctest::Approx(30.0));
    CHECK(vec_v.z.data.grad[0] == doctest::Approx(40.0));

    Loci::vector3d<Loci::VFAD> vec_v_default;
    opts.getOptionUnits("vec", "m", vec_v_default, 0);
    CHECK(vec_v_default.x.data.value == doctest::Approx(1.0));
    CHECK(vec_v_default.y.data.value == doctest::Approx(2.0));
    CHECK(vec_v_default.z.data.value == doctest::Approx(3.0));
    CHECK(vec_v_default.x.data.grad[0] == doctest::Approx(2.0));
    CHECK(vec_v_default.y.data.grad[0] == doctest::Approx(3.0));
    CHECK(vec_v_default.z.data.grad[0] == doctest::Approx(4.0));
  }
}

TEST_CASE("options_list vector getOptionUnits rejects malformed vector specifications") {
  {
    Loci::options_list opts = parse_options("<v=[1,2]>");
    Loci::vector3d<double> vec;
    CHECK_THROWS(opts.getOptionUnits("v", "m", vec));
  }
  {
    Loci::options_list opts = parse_options("<v=[1,$true,3]>");
    Loci::vector3d<double> vec;
    CHECK_THROWS(opts.getOptionUnits("v", "m", vec));
  }
  {
    Loci::options_list opts = parse_options("<v=cyl(1,2,3)>");
    Loci::vector3d<double> vec;
    CHECK_THROWS(opts.getOptionUnits("v", "m", vec));
  }
  {
    Loci::options_list opts = parse_options("<v=[1 s,2 s,3 s]>");
    Loci::vector3d<double> vec;
    CHECK_THROWS(opts.getOptionUnits("v", "m", vec));
  }
  {
    Loci::options_list opts = parse_options("<v=polar(1,2)>");
    Loci::vector3d<double> vec;
    CHECK_THROWS(opts.getOptionUnits("v", "m", vec));
  }
  {
    Loci::options_list opts = parse_options("<v=polar(1,$true,0)>");
    Loci::vector3d<double> vec;
    CHECK_THROWS(opts.getOptionUnits("v", "m", vec));
  }
  {
    Loci::options_list opts = parse_options("<v=name>");
    Loci::vector3d<double> vec;
    CHECK_THROWS(opts.getOptionUnits("v", "m", vec));
  }
}

TEST_CASE("option_values parser rejects invalid boolean and invalid list assignment lhs") {
  {
    Loci::option_values v;
    std::istringstream in("$maybe");
    CHECK_THROWS(in >> v);
  }
  {
    Loci::option_values v;
    std::istringstream in("[1=2]");
    CHECK_THROWS(in >> v);
  }
}

TEST_CASE("options_list getOption throws for missing option") {
  Loci::options_list opts;
  bool value = false;
  CHECK_THROWS(opts.getOption("missing", value));
}

TEST_CASE("options_list getOptionUnits throws for missing option and non-unit types") {
  {
    Loci::options_list opts;
    double value = 0.0;
    CHECK_THROWS(opts.getOptionUnits("missing", "m", value));
  }
  {
    Loci::options_list opts = parse_options("<x=name>");
    double value = 0.0;
    CHECK_THROWS(opts.getOptionUnits("x", "m", value));
  }
}

TEST_CASE("parse-print-parse idempotence holds for randomized printable options corpus") {
  std::mt19937 rng(424242);

  auto make_value = [](std::mt19937 &gen, int idx) {
    switch(static_cast<int>(gen() % 8)) {
    case 0:
      return std::to_string(static_cast<int>(gen() % 200) - 100);
    case 1:
      return std::string((gen() % 2) ? "$true" : "$false");
    case 2:
      return std::string((gen() % 2) ? "alpha" : "beta_2");
    case 3:
      return std::string("\"file_") + std::to_string(idx) + ".dat\"";
    case 4:
      return std::string("[") + std::to_string(idx % 5) + "," + std::to_string((idx + 1) % 7) + "]";
    case 5:
      return std::string("f(") + std::to_string((idx % 3) + 1) + ",$true)";
    case 6:
      return std::to_string((idx % 9) + 1) + " m/s";
    default:
      return std::to_string((idx % 9) + 1) + " one/s";
    }
  };

  for(int trial = 0; trial < 120; ++trial) {
    const int nopts = 1 + static_cast<int>(rng() % 6);
    std::string text = "<";
    for(int i = 0; i < nopts; ++i) {
      if(i > 0)
        text += ",";
      text += "o" + std::to_string(i) + "=" + make_value(rng, trial * 10 + i);
    }
    text += ">";

    CAPTURE(text);
    Loci::options_list first;
    std::istringstream in1(text);
    CHECK_NOTHROW(in1 >> first);

    std::ostringstream out1;
    out1 << first;
    const std::string p1 = out1.str();

    Loci::options_list second;
    std::istringstream in2(p1);
    CHECK_NOTHROW(in2 >> second);

    std::ostringstream out2;
    out2 << second;
    const std::string p2 = out2.str();

    CHECK(p1 == p2);
  }
}

TEST_CASE("randomized malformed options corpus is rejected") {
  std::mt19937 rng(1337);

  for(int trial = 0; trial < 60; ++trial) {
    std::string base = "<a=" + std::to_string(trial + 1) + ",b=\"file_" + std::to_string(trial) + ".dat\">";
    std::string bad = base;

    switch(static_cast<int>(rng() % 4)) {
    case 0:
      bad.erase(0, 1);
      break;
    case 1:
      bad.pop_back();
      break;
    case 2: {
      size_t pos = bad.find('=');
      if(pos != std::string::npos)
        bad.insert(pos, "=");
      break;
    }
    default: {
      size_t pos = bad.find(',');
      if(pos != std::string::npos)
        bad[pos] = ';';
      break;
    }
    }

    CAPTURE(base);
    CAPTURE(bad);
    Loci::options_list opts;
    std::istringstream in(bad);
    CHECK_THROWS(in >> opts);
  }
}
