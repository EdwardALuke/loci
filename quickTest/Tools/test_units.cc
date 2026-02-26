#include <cmath>
#include <string_view>
#include <sstream>
#include <string>
#include <vector>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Tools/options_list.h>
#include <Tools/unit_type.h>

namespace {
Loci::options_list parse_options(const std::string &text) {
  Loci::options_list opts;
  std::istringstream in(text);
  in >> opts;
  return opts;
}

bool is_temperature_alias(std::string_view name) {
  return name == "Fahrenheit" || name == "F" ||
         name == "Rankine" || name == "R" ||
         name == "Celsius" || name == "C";
}

template <typename TableT>
int find_table_index(const TableT *table, std::string_view name) {
  for(int i = 0; table[i].name != 0; ++i) {
    if(name == table[i].name) {
      return i;
    }
  }
  return -1;
}
} // namespace

TEST_CASE("options_list parses pressure units and converts to Pa") {
  Loci::options_list opts = parse_options("<p=101.325 kPa>");

  CHECK(opts.getOptionValueType("p") == Loci::UNIT_VALUE);
  double p_pa = 0.0;
  opts.getOptionUnits("p", "Pa", p_pa);
  CHECK(p_pa == doctest::Approx(101325.0));
}

TEST_CASE("options_list parses reciprocal rates with explicit one/s") {
  Loci::options_list opts = parse_options("<sigma=5 one/s>");

  double sigma = 0.0;
  opts.getOptionUnits("sigma", "one/s", sigma);
  CHECK(sigma == doctest::Approx(5.0));

  double sigma_hz = 0.0;
  opts.getOptionUnits("sigma", "Hz", sigma_hz);
  CHECK(sigma_hz == doctest::Approx(5.0));
}

TEST_CASE("options_list rejects reciprocal rates using 5/s shorthand") {
  CHECK_THROWS(parse_options("<sigma=5/s>"));
}

TEST_CASE("options_list rejects reciprocal rates using spaced 5 / s shorthand") {
  CHECK_THROWS(parse_options("<sigma=5 / s>"));
}

TEST_CASE("options_list rejects explicit 1/s unit syntax in favor of one/s") {
  CHECK_THROWS(parse_options("<sigma=5 1/s>"));
}

TEST_CASE("options_list parses kWh, kW*h and perm reference aliases") {
  Loci::options_list opts = parse_options("<energy=2 kW*h,energy2=3 kWh,permv=1 perm>");

  double energy_j = 0.0;
  opts.getOptionUnits("energy", "joule", energy_j);
  CHECK(energy_j == doctest::Approx(7.2e6));

  double energy2_j = 0.0;
  opts.getOptionUnits("energy2", "joule", energy2_j);
  CHECK(energy2_j == doctest::Approx(1.08e7));

  double perm_si = 0.0;
  opts.getOptionUnits("permv", "kg/Pa/s/m/m", perm_si);
  CHECK(perm_si == doctest::Approx(5.72135e-11));
}

TEST_CASE("options_list parses parenthesized reciprocal unit expressions") {
  Loci::options_list opts = parse_options("<muinv=2 one/(Pa*s)>");

  double muinv = 0.0;
  opts.getOptionUnits("muinv", "one/(Pa*s)", muinv);
  CHECK(muinv == doctest::Approx(2.0));

  double muinv_rhe = 0.0;
  opts.getOptionUnits("muinv", "rhe", muinv_rhe);
  CHECK(muinv_rhe == doctest::Approx(0.2));
}

TEST_CASE("options_list parses angle units in degrees") {
  Loci::options_list opts = parse_options("<angle=180 deg>");

  double angle_rad = 0.0;
  opts.getOptionUnits("angle", "radians", angle_rad);
  CHECK(angle_rad == doctest::Approx(std::acos(-1.0)));
}

TEST_CASE("options_list preserves derivative lists with unit conversion") {
  Loci::options_list opts = parse_options("<vel=10^[2,3,4] m/s>");

  Loci::UNIT_type vel;
  opts.getOption("vel", vel);

  double value = 0.0;
  double grad = 0.0;
  double grad2 = 0.0;
  std::vector<double> grad_list;
  vel.get_values_in("km/h", value, grad, grad2, grad_list);

  CHECK(value == doctest::Approx(36.0));
  CHECK(grad == doctest::Approx(7.2));
  CHECK(grad2 == doctest::Approx(0.0));
  REQUIRE(grad_list.size() == 2);
  CHECK(grad_list[0] == doctest::Approx(10.8));
  CHECK(grad_list[1] == doctest::Approx(14.4));
}

TEST_CASE("options_list parses absolute Fahrenheit values") {
  Loci::options_list opts = parse_options("<T=32 F>");

  double t_kelvin = 0.0;
  opts.getOptionUnits("T", "kelvin", t_kelvin);
  CHECK(t_kelvin == doctest::Approx(273.15));
}

TEST_CASE("options_list rejects numeric coefficients inside unit expressions") {
  Loci::options_list opts;
  std::istringstream in("<mu=1 1.5/(Pa*s)>");
  CHECK_THROWS(in >> opts);
}

TEST_CASE("options_list rejects malformed trailing unit operators") {
  Loci::options_list opts;
  std::istringstream in("<bad=3 m/>");
  CHECK_THROWS(in >> opts);
}

TEST_CASE("options_list rejects unknown unit symbol") {
  Loci::options_list opts;
  std::istringstream in("<x=1 blarg>");
  CHECK_THROWS(in >> opts);
}

TEST_CASE("options_list rejects unterminated parenthesized unit expression") {
  Loci::options_list opts;
  std::istringstream in("<x=1 one/(Pa*s>");
  CHECK_THROWS(in >> opts);
}

TEST_CASE("options_list rejects incompatible target units in getOptionUnits") {
  Loci::options_list opts = parse_options("<len=1 m>");
  double value = 0.0;
  CHECK_THROWS(opts.getOptionUnits("len", "s", value));
}

TEST_CASE("UNIT_type parser supports derivative blocks and parenthesized units") {
  Loci::UNIT_type unit;
  std::istringstream in("1^[2 3] one/(Pa*s)");
  in >> unit;

  CHECK(unit.is_compatible("one/(Pa*s)"));

  double value = 0.0;
  double grad = 0.0;
  double grad2 = 0.0;
  std::vector<double> grad_list;
  unit.get_values_in("one/(Pa*s)", value, grad, grad2, grad_list);

  CHECK(value == doctest::Approx(1.0));
  CHECK(grad == doctest::Approx(2.0));
  CHECK(grad2 == doctest::Approx(0.0));
  REQUIRE(grad_list.size() == 1);
  CHECK(grad_list[0] == doctest::Approx(3.0));
}

TEST_CASE("UNIT_type parser accepts reciprocal units written as one/s") {
  Loci::UNIT_type unit;
  std::istringstream in("5 one/s");
  in >> unit;

  CHECK(unit.is_compatible("one/s"));
  CHECK(unit.get_value_in("one/s") == doctest::Approx(5.0));
  CHECK(unit.get_value_in("Hz") == doctest::Approx(5.0));
}

TEST_CASE("UNIT_type parser rejects numeric coefficients in unit expression") {
  Loci::UNIT_type unit;
  std::istringstream in("1 1.5/s");
  CHECK_THROWS(in >> unit);
}

TEST_CASE("UNIT_type parser rejects numeric reciprocal shorthand syntaxes") {
  {
    Loci::UNIT_type unit;
    std::istringstream in("5/s");
    CHECK_THROWS(in >> unit);
  }
  {
    Loci::UNIT_type unit;
    std::istringstream in("5 / s");
    CHECK_THROWS(in >> unit);
  }
  {
    Loci::UNIT_type unit;
    std::istringstream in("1 1/s");
    CHECK_THROWS(in >> unit);
  }
}

TEST_CASE("UNIT_type parser rejects unknown unit symbols") {
  Loci::UNIT_type unit;
  std::istringstream in("1 blarg");
  CHECK_THROWS(in >> unit);
}

TEST_CASE("UNIT_type parser rejects unterminated parenthesized units") {
  Loci::UNIT_type unit;
  std::istringstream in("1 one/(Pa*s");
  CHECK_THROWS(in >> unit);
}

TEST_CASE("UNIT_type one alias behaves as nondimensional in MKS and CGS") {
  {
    Loci::UNIT_type u(Loci::UNIT_type::MKS, "general", 7.5, "one");
    CHECK(u.is_compatible("meter/meter"));
    CHECK(u.is_compatible("second/second"));
    CHECK(u.get_value_in("one") == doctest::Approx(7.5));
    CHECK(u.get_value_in("meter/meter") == doctest::Approx(7.5));
  }
  {
    Loci::UNIT_type u(Loci::UNIT_type::CGS, "general", 2.25, "one");
    CHECK(u.is_compatible("centimeter/centimeter"));
    CHECK(u.get_value_in("one") == doctest::Approx(2.25));
    CHECK(u.get_value_in("centimeter/centimeter") == doctest::Approx(2.25));
  }
}

TEST_CASE("UNIT_type constructor rejects empty unit strings") {
  CHECK_THROWS(Loci::UNIT_type(Loci::UNIT_type::MKS, "length", 3.0, ""));
  CHECK_THROWS(Loci::UNIT_type(Loci::UNIT_type::MKS, "general", 1.0, ""));
}

TEST_CASE("UNIT_type temperature reverse conversion scales derivative outputs") {
  Loci::UNIT_type temp(Loci::UNIT_type::MKS, "general", "K", 300.0, 2.0, 3.0, {4.0, 5.0});

  double value = 0.0;
  double grad = 0.0;
  double grad2 = 0.0;
  std::vector<double> grad_list;
  temp.get_values_in("F", value, grad, grad2, grad_list);

  CHECK(value == doctest::Approx(80.33).epsilon(1e-4));
  CHECK(grad == doctest::Approx(3.6));
  CHECK(grad2 == doctest::Approx(5.4));
  REQUIRE(grad_list.size() == 2);
  CHECK(grad_list[0] == doctest::Approx(7.2));
  CHECK(grad_list[1] == doctest::Approx(9.0));
}

TEST_CASE("UNIT_type rhe alias remains consistent with one/(Pa*s)") {
  Loci::UNIT_type rhe(Loci::UNIT_type::MKS, "general", 1.0, "rhe");
  CHECK(rhe.is_compatible("one/(Pa*s)"));
  CHECK(rhe.get_value_in("one/(Pa*s)") == doctest::Approx(10.0));

  Loci::UNIT_type inv_visc(Loci::UNIT_type::MKS, "general", 1.0, "one/(Pa*s)");
  CHECK(inv_visc.get_value_in("rhe") == doctest::Approx(0.1));
}

TEST_CASE("UNIT_type conversion round-trip is stable for representative compatible pairs") {
  struct Case {
    const char *src;
    const char *dst;
    double value;
  };
  const std::vector<Case> cases = {
      {"one/s", "Hz", 2.75},
      {"kHz", "one/s", 0.5},
      {"kWh", "J", 1.25},
      {"perm", "kg/Pa/s/m/m", 3.0},
      {"m/s", "km/h", 6.5},
  };

  for(const auto &tc : cases) {
    CAPTURE(tc.src);
    CAPTURE(tc.dst);
    CAPTURE(tc.value);
    Loci::UNIT_type u1(Loci::UNIT_type::MKS, "general", tc.value, tc.src);
    const double in_dst = u1.get_value_in(tc.dst);
    Loci::UNIT_type u2(Loci::UNIT_type::MKS, "general", in_dst, tc.dst);
    CHECK(u2.get_value_in(tc.src) == doctest::Approx(tc.value).epsilon(1e-9));
  }
}

TEST_CASE("UNIT_type stream roundtrip preserves one-based reciprocal units and derivatives") {
  Loci::UNIT_type original;
  std::istringstream in("5^[2 3 4] one/s");
  in >> original;

  std::ostringstream out;
  out << original;

  Loci::UNIT_type parsed;
  std::istringstream in2(out.str());
  in2 >> parsed;

  double value = 0.0;
  double grad = 0.0;
  double grad2 = 0.0;
  std::vector<double> grad_list;
  parsed.get_values_in("Hz", value, grad, grad2, grad_list);

  CHECK(value == doctest::Approx(5.0));
  CHECK(grad == doctest::Approx(2.0));
  CHECK(grad2 == doctest::Approx(0.0));
  REQUIRE(grad_list.size() == 2);
  CHECK(grad_list[0] == doctest::Approx(3.0));
  CHECK(grad_list[1] == doctest::Approx(4.0));
}

TEST_CASE("new one/kWh/perm aliases are present with expected canonical definitions") {
  const int mks_one = find_table_index(Loci::UNIT_type::composite_unit_table, "one");
  REQUIRE(mks_one >= 0);
  CHECK(std::string(Loci::UNIT_type::composite_unit_table[mks_one].derived_unit) == "meter/meter");
  CHECK(Loci::UNIT_type::composite_unit_table[mks_one].convert_factor == doctest::Approx(1.0));

  const int cgs_one = find_table_index(Loci::UNIT_type::cgs_composite_unit_table, "one");
  REQUIRE(cgs_one >= 0);
  CHECK(std::string(Loci::UNIT_type::cgs_composite_unit_table[cgs_one].derived_unit) == "centimeter/centimeter");
  CHECK(Loci::UNIT_type::cgs_composite_unit_table[cgs_one].convert_factor == doctest::Approx(1.0));

  const int kwh = find_table_index(Loci::UNIT_type::reference_unit_table, "kWh");
  REQUIRE(kwh >= 0);
  CHECK(std::string(Loci::UNIT_type::reference_unit_table[kwh].refer_unit) == "joule");
  CHECK(Loci::UNIT_type::reference_unit_table[kwh].convert_factor == doctest::Approx(3.6e6));

  const int kwh_star = find_table_index(Loci::UNIT_type::reference_unit_table, "kW*h");
  REQUIRE(kwh_star >= 0);
  CHECK(std::string(Loci::UNIT_type::reference_unit_table[kwh_star].refer_unit) == "joule");
  CHECK(Loci::UNIT_type::reference_unit_table[kwh_star].convert_factor == doctest::Approx(3.6e6));

  const int perm = find_table_index(Loci::UNIT_type::reference_unit_table, "perm");
  REQUIRE(perm >= 0);
  CHECK(std::string(Loci::UNIT_type::reference_unit_table[perm].refer_unit) == "kg/Pa/s/m/m");
  CHECK(Loci::UNIT_type::reference_unit_table[perm].convert_factor == doctest::Approx(5.72135e-11));
}

TEST_CASE("UNIT_type constructor rejects numeric-leading unit expressions") {
  CHECK_THROWS(Loci::UNIT_type(Loci::UNIT_type::MKS, "general", 1.0, "1/s"));
  CHECK_THROWS(Loci::UNIT_type(Loci::UNIT_type::MKS, "general", 1.0, "5/s"));
}

TEST_CASE("UNIT_type parser accepts one-based spacing variants") {
  const std::vector<std::string> accepted = {
      "5 one / s",
      "5 one/ s",
      "5 one /s",
      "2 one /( Pa * s )",
  };

  for(const auto &text : accepted) {
    CAPTURE(text);
    Loci::UNIT_type unit;
    std::istringstream in(text);
    CHECK_NOTHROW(in >> unit);
  }

  {
    Loci::UNIT_type unit;
    std::istringstream in("5 one / s");
    in >> unit;
    CHECK(unit.get_value_in("Hz") == doctest::Approx(5.0));
  }
  {
    Loci::UNIT_type unit;
    std::istringstream in("2 one /( Pa * s )");
    in >> unit;
    CHECK(unit.get_value_in("rhe") == doctest::Approx(0.2));
  }
}

TEST_CASE("options_list treats one/s, Hz and hertz as equivalent reciprocal-rate units") {
  Loci::options_list opts = parse_options("<a=5 one/s,b=5 Hz,c=5 hertz>");

  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
  opts.getOptionUnits("a", "one/s", a);
  opts.getOptionUnits("b", "one/s", b);
  opts.getOptionUnits("c", "one/s", c);
  CHECK(a == doctest::Approx(5.0));
  CHECK(b == doctest::Approx(5.0));
  CHECK(c == doctest::Approx(5.0));
}

TEST_CASE("UNIT_type stream parser handles every registered alias used in conversion tables") {
  int tested_composite = 0;
  int tested_reference = 0;

  for(int i = 0; Loci::UNIT_type::composite_unit_table[i].name != 0; ++i) {
    const std::string name = Loci::UNIT_type::composite_unit_table[i].name;
    if(is_temperature_alias(name)) {
      continue;
    }
    const std::string derived = Loci::UNIT_type::composite_unit_table[i].derived_unit;
    const double factor = Loci::UNIT_type::composite_unit_table[i].convert_factor;
    CAPTURE(name);
    Loci::UNIT_type u;
    std::istringstream in(std::string("1 ") + name);
    CHECK_NOTHROW(in >> u);
    CHECK(u.get_value_in(derived) == doctest::Approx(factor).epsilon(1e-9));
    ++tested_composite;
  }

  for(int i = 0; Loci::UNIT_type::reference_unit_table[i].name != 0; ++i) {
    const std::string name = Loci::UNIT_type::reference_unit_table[i].name;
    const std::string refer = Loci::UNIT_type::reference_unit_table[i].refer_unit;
    const double factor = Loci::UNIT_type::reference_unit_table[i].convert_factor;
    CAPTURE(name);
    Loci::UNIT_type u;
    std::istringstream in(std::string("1 ") + name);
    CHECK_NOTHROW(in >> u);
    CHECK(u.get_value_in(refer) == doctest::Approx(factor).epsilon(1e-9));
    ++tested_reference;
  }

  CHECK(tested_composite > 0);
  CHECK(tested_reference > 0);
}

TEST_CASE("MKS composite unit table entries are conversion-consistent") {
  int tested = 0;
  int skipped_temperature = 0;
  for(int i=0; Loci::UNIT_type::composite_unit_table[i].name != 0; ++i) {
    const std::string name = Loci::UNIT_type::composite_unit_table[i].name;
    if(is_temperature_alias(name)) {
      ++skipped_temperature;
      continue;
    }
    const std::string derived = Loci::UNIT_type::composite_unit_table[i].derived_unit;
    const double factor = Loci::UNIT_type::composite_unit_table[i].convert_factor;
    CAPTURE(name);
    CAPTURE(derived);
    CAPTURE(factor);

    Loci::UNIT_type u(Loci::UNIT_type::MKS, "general", 1.0, name);
    CHECK(u.get_value_in(derived) == doctest::Approx(factor).epsilon(1e-9));
    ++tested;
  }
  CHECK(tested > 0);
  CHECK(skipped_temperature == 6);
}

TEST_CASE("reference unit table entries are conversion-consistent") {
  int tested = 0;
  for(int i=0; Loci::UNIT_type::reference_unit_table[i].name != 0; ++i) {
    const std::string name = Loci::UNIT_type::reference_unit_table[i].name;
    const std::string refer = Loci::UNIT_type::reference_unit_table[i].refer_unit;
    const double factor = Loci::UNIT_type::reference_unit_table[i].convert_factor;
    CAPTURE(name);
    CAPTURE(refer);
    CAPTURE(factor);

    Loci::UNIT_type u(Loci::UNIT_type::MKS, "general", 1.0, name);
    CHECK(u.get_value_in(refer) == doctest::Approx(factor).epsilon(1e-9));
    ++tested;
  }
  CHECK(tested > 0);
}

TEST_CASE("CGS composite unit table entries are conversion-consistent") {
  int tested = 0;
  int skipped_temperature = 0;
  for(int i=0; Loci::UNIT_type::cgs_composite_unit_table[i].name != 0; ++i) {
    const std::string name = Loci::UNIT_type::cgs_composite_unit_table[i].name;
    if(is_temperature_alias(name)) {
      ++skipped_temperature;
      continue;
    }
    const std::string derived = Loci::UNIT_type::cgs_composite_unit_table[i].derived_unit;
    const double factor = Loci::UNIT_type::cgs_composite_unit_table[i].convert_factor;
    CAPTURE(name);
    CAPTURE(derived);
    CAPTURE(factor);

    Loci::UNIT_type u(Loci::UNIT_type::CGS, "general", 1.0, name);
    CHECK(u.get_value_in(derived) == doctest::Approx(factor).epsilon(1e-9));
    ++tested;
  }
  CHECK(tested > 0);
  CHECK(skipped_temperature == 6);
}
