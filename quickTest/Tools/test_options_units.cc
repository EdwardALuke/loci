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
