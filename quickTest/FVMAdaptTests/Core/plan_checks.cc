//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#include "fvmadapt_core.h"

namespace {

/// One split-plan case and its expected canonical form.
struct PlanCase {
  PlanCase(const std::string& n, const std::vector<char>& p,
           const std::vector<char>& c, bool example):
           name(n), plan(p), canonical(c), write_example(example) {}

  std::string name ;
  std::vector<char> plan ;
  std::vector<char> canonical ;
  bool write_example ;
} ;

/// Build a one-entry plan for a root split code.
std::vector<char> root_plan(int split_code) {
  return std::vector<char>(1, char(split_code)) ;
}

/// Return the canonical plan for one root split code.
std::vector<char> canonical_root_plan(int split_code) {
  if(split_code == 0) {
    return std::vector<char>() ;
  }
  return root_plan(split_code) ;
}

/// Build a filesystem-safe token for generated case artifacts.
std::string plan_file_token(const std::string& shape, const std::string& name) {
  return shape + "_" + name ;
}

/// Return the first split code, or zero for an empty plan.
int root_split_code(const std::vector<char>& plan) {
  return plan.empty() ? 0 : int(plan[0]) ;
}

/// Replay one HexCell plan and check count and canonical-plan invariants.
CaseResult run_hex_case(const PlanCase& plan_case,
                        bool write_example_files,
                        const std::string& output_dir) {
  HexFixture fixture ;
  build_unit_hex(fixture) ;

  std::vector<HexCell*> resplit_leaves ;
  fixture.root->resplit(plan_case.plan, fixture.nodes, fixture.edges,
                        fixture.faces, resplit_leaves) ;

  std::list<HexCell*> sorted_leaves ;
  fixture.root->sort_leaves(sorted_leaves) ;

  HexCell counter ;
  const int counted = counter.num_fine_cells(plan_case.plan) ;

  if(int(resplit_leaves.size()) != counted) {
    std::ostringstream msg ;
    msg << "hex plan " << plan_case.name << " produced "
        << resplit_leaves.size() << " resplit leaves, counted " << counted ;
    throw std::runtime_error(msg.str()) ;
  }

  if(int(sorted_leaves.size()) != counted) {
    std::ostringstream msg ;
    msg << "hex plan " << plan_case.name << " produced "
        << sorted_leaves.size() << " sorted leaves, counted " << counted ;
    throw std::runtime_error(msg.str()) ;
  }

  const std::vector<char> round_trip = fixture.root->make_cellplan() ;
  if(round_trip != plan_case.canonical) {
    std::ostringstream msg ;
    msg << "hex plan " << plan_case.name << " round-tripped as "
        << plan_string(round_trip) << ", expected "
        << plan_string(plan_case.canonical) ;
    throw std::runtime_error(msg.str()) ;
  }

  CaseResult result ;
  result.shape = "hex" ;
  result.case_name = plan_case.name ;
  result.input_plan = plan_case.plan ;
  result.canonical_plan = round_trip ;
  result.counted_leaves = counted ;
  result.resplit_leaves = int(resplit_leaves.size()) ;
  result.sorted_leaves = int(sorted_leaves.size()) ;
  result.point_count = 0 ;

  if(write_example_files && plan_case.write_example) {
    std::ostringstream name ;
    name << plan_file_token(result.shape, plan_case.name) << ".vtk" ;
    result.vtk_file = join_path(output_dir, name.str()) ;
    std::ostringstream title ;
    title << "FVMAdapt HexCell " << plan_case.name
          << " plan " << plan_string(plan_case.plan) ;
    result.point_count = write_hex_vtk(result.vtk_file, title.str(),
                                       root_split_code(plan_case.plan),
                                       resplit_leaves) ;
  }

  return result ;
}

/// Run one Prism split plan through the real Prism library paths. This checks
/// that resplit(), sort_leaves(), empty_resplit(), and make_cellplan() agree on
/// the leaf count and canonical plan; example mode also writes a VTK wireframe.
CaseResult run_prism_case(const PlanCase& plan_case, bool write_example_files,
                          const std::string& output_dir) {
  PrismFixture fixture ;
  build_unit_prism(fixture) ;

  std::vector<Prism*> resplit_leaves ;
  fixture.root->resplit(plan_case.plan, fixture.nodes, fixture.edges,
                        fixture.qfaces, fixture.gfaces, resplit_leaves) ;

  std::list<Prism*> sorted_leaves ;
  fixture.root->sort_leaves(sorted_leaves) ;

  PrismFixture counter ;
  build_unit_prism(counter) ;
  const int counted = counter.root->empty_resplit(plan_case.plan) ;

  if(int(resplit_leaves.size()) != counted) {
    std::ostringstream msg ;
    msg << "prism plan " << plan_case.name << " produced "
        << resplit_leaves.size() << " resplit leaves, counted " << counted ;
    throw std::runtime_error(msg.str()) ;
  }

  if(int(sorted_leaves.size()) != counted) {
    std::ostringstream msg ;
    msg << "prism plan " << plan_case.name << " produced "
        << sorted_leaves.size() << " sorted leaves, counted " << counted ;
    throw std::runtime_error(msg.str()) ;
  }

  const std::vector<char> round_trip = fixture.root->make_cellplan() ;
  if(round_trip != plan_case.canonical) {
    std::ostringstream msg ;
    msg << "prism plan " << plan_case.name << " round-tripped as "
        << plan_string(round_trip) << ", expected "
        << plan_string(plan_case.canonical) ;
    throw std::runtime_error(msg.str()) ;
  }

  CaseResult result ;
  result.shape = "prism" ;
  result.case_name = plan_case.name ;
  result.input_plan = plan_case.plan ;
  result.canonical_plan = round_trip ;
  result.counted_leaves = counted ;
  result.resplit_leaves = int(resplit_leaves.size()) ;
  result.sorted_leaves = int(sorted_leaves.size()) ;
  result.point_count = 0 ;

  if(write_example_files && plan_case.write_example) {
    std::ostringstream name ;
    name << plan_file_token(result.shape, plan_case.name) << "_wire.vtk" ;
    result.vtk_file = join_path(output_dir, name.str()) ;
    std::ostringstream title ;
    title << "FVMAdapt Prism " << plan_case.name
          << " plan " << plan_string(plan_case.plan) << " wireframe" ;
    result.point_count = write_wireframe_vtk(result.vtk_file, title.str(),
                                             fixture.nodes, fixture.edges) ;
  }

  return result ;
}

/// Replay one general Cell plan and check count and canonical-plan invariants.
CaseResult run_general_case(const std::string& shape,
                            const std::string& title_shape,
                            GeneralBuilder build_cell,
                            const PlanCase& plan_case,
                            bool write_example_files,
                            const std::string& output_dir) {
  GeneralFixture fixture ;
  build_cell(fixture) ;

  std::vector<DiamondCell*> resplit_leaves ;
  if(!plan_case.plan.empty()) {
    fixture.root->resplit(plan_case.plan,
                          fixture.nodes,
                          fixture.edges,
                          fixture.faces,
                          resplit_leaves) ;
  }

  std::list<DiamondCell*> sorted_leaves ;
  if(!plan_case.plan.empty()) {
    fixture.root->sort_leaves(sorted_leaves) ;
  }

  GeneralFixture counter ;
  build_cell(counter) ;
  const int counted = counter.root->empty_resplit(plan_case.plan) ;
  const int resplit_count = plan_case.plan.empty() ? 1 : int(resplit_leaves.size()) ;
  const int sorted_count = plan_case.plan.empty() ? 1 : int(sorted_leaves.size()) ;

  if(resplit_count != counted) {
    std::ostringstream msg ;
    msg << title_shape << " plan " << plan_case.name << " produced "
        << resplit_count << " resplit leaves, counted " << counted ;
    throw std::runtime_error(msg.str()) ;
  }

  if(sorted_count != counted) {
    std::ostringstream msg ;
    msg << title_shape << " plan " << plan_case.name << " produced "
        << sorted_count << " sorted leaves, counted " << counted ;
    throw std::runtime_error(msg.str()) ;
  }

  const std::vector<char> round_trip = fixture.root->make_cellplan() ;
  if(round_trip != plan_case.canonical) {
    std::ostringstream msg ;
    msg << title_shape << " plan " << plan_case.name << " round-tripped as "
        << plan_string(round_trip) << ", expected "
        << plan_string(plan_case.canonical) ;
    throw std::runtime_error(msg.str()) ;
  }

  CaseResult result ;
  result.shape = shape ;
  result.case_name = plan_case.name ;
  result.input_plan = plan_case.plan ;
  result.canonical_plan = round_trip ;
  result.counted_leaves = counted ;
  result.resplit_leaves = resplit_count ;
  result.sorted_leaves = sorted_count ;
  result.point_count = 0 ;

  if(write_example_files && plan_case.write_example) {
    std::ostringstream name ;
    name << plan_file_token(shape, plan_case.name) << "_wire.vtk" ;
    result.vtk_file = join_path(output_dir, name.str()) ;
    std::ostringstream title ;
    title << "FVMAdapt " << title_shape << " " << plan_case.name
          << " plan " << plan_string(plan_case.plan)
          << " wireframe" ;
    result.point_count = write_wireframe_vtk(result.vtk_file, title.str(),
                                             fixture.nodes, fixture.edges) ;
  }

  return result ;
}

/// Return representative HexCell split-plan cases.
std::vector<PlanCase> hex_plan_cases() {
  std::vector<PlanCase> cases ;
  for(int split_code = 0; split_code <= 7; ++split_code) {
    std::ostringstream name ;
    name << "split_" << split_code ;
    cases.push_back(PlanCase(name.str(),
                             root_plan(split_code),
                             canonical_root_plan(split_code),
                             true)) ;
  }

  cases.push_back(PlanCase("split_7_trailing_zeros",
                           {char(7), char(0), char(0), char(0), char(0),
                            char(0), char(0), char(0), char(0)},
                           {char(7)},
                           false)) ;
  cases.push_back(PlanCase("split_7_child0_z",
                           {char(7), char(1)},
                           {char(7), char(1)},
                           true)) ;
  cases.push_back(PlanCase("split_7_child7_z",
                           {char(7), char(0), char(0), char(0), char(0),
                            char(0), char(0), char(0), char(1)},
                           {char(7), char(0), char(0), char(0), char(0),
                            char(0), char(0), char(0), char(1)},
                           true)) ;
  cases.push_back(PlanCase("split_4_child0_yz",
                           {char(4), char(3)},
                           {char(4), char(3)},
                           true)) ;
  return cases ;
}

/// Return representative Prism split-plan cases.
std::vector<PlanCase> prism_plan_cases() {
  std::vector<PlanCase> cases ;
  for(int split_code = 0; split_code <= 3; ++split_code) {
    std::ostringstream name ;
    name << "split_" << split_code ;
    cases.push_back(PlanCase(name.str(),
                             root_plan(split_code),
                             canonical_root_plan(split_code),
                             true)) ;
  }

  cases.push_back(PlanCase("split_3_trailing_zeros",
                           {char(3), char(0), char(0), char(0), char(0),
                            char(0), char(0)},
                           {char(3)},
                           false)) ;
  cases.push_back(PlanCase("split_3_child0_z",
                           {char(3), char(1)},
                           {char(3), char(1)},
                           true)) ;
  cases.push_back(PlanCase("split_2_child0_z",
                           {char(2), char(1)},
                           {char(2), char(1)},
                           true)) ;
  return cases ;
}

/// Return representative general Cell split-plan cases.
std::vector<PlanCase> general_plan_cases() {
  std::vector<PlanCase> cases ;
  cases.push_back(PlanCase("split_0", std::vector<char>(), std::vector<char>(), true)) ;
  cases.push_back(PlanCase("split_1", {char(1)}, {char(1)}, true)) ;
  cases.push_back(PlanCase("split_1_trailing_zeros",
                           {char(1), char(0), char(0), char(0), char(0)},
                           {char(1)},
                           false)) ;
  cases.push_back(PlanCase("split_1_child0",
                           {char(1), char(1)},
                           {char(1), char(1)},
                           true)) ;
  return cases ;
}

} // namespace

/// Check that general-face and quad-face plan conversions preserve splits.
void check_plan_transfer_round_trips() {
  const std::vector<std::vector<char> > general_face_plans = {
    std::vector<char>(),
    {char(1)},
    {char(1), char(1)},
    {char(1), char(0), char(1)}
  } ;

  for(size_t i = 0; i < general_face_plans.size(); ++i) {
    std::vector<char> general_plan = general_face_plans[i] ;
    Face canonical_face(4) ;
    canonical_face.empty_resplit(general_plan) ;
    const std::vector<char> canonical = canonical_face.make_faceplan() ;
    const std::vector<char> quad_plan = transfer_plan_g2q(general_plan) ;
    const std::vector<char> round_trip = transfer_plan_q2g(quad_plan) ;
    if(round_trip != canonical) {
      std::ostringstream msg ;
      msg << "general face plan " << plan_string(general_face_plans[i])
          << " transferred through quad as " << plan_string(round_trip)
          << ", expected " << plan_string(canonical) ;
      throw std::runtime_error(msg.str()) ;
    }
  }

  const std::vector<std::vector<char> > quad_face_plans = {
    std::vector<char>(),
    {char(3)},
    {char(3), char(3)},
    {char(3), char(0), char(0), char(0), char(3)}
  } ;

  for(size_t i = 0; i < quad_face_plans.size(); ++i) {
    std::vector<char> general_plan = transfer_plan_q2g(quad_face_plans[i]) ;
    const std::vector<char> round_trip = transfer_plan_g2q(general_plan) ;
    if(round_trip != quad_face_plans[i]) {
      std::ostringstream msg ;
      msg << "quad face plan " << plan_string(quad_face_plans[i])
          << " transferred through general as " << plan_string(round_trip)
          << ", expected " << plan_string(quad_face_plans[i]) ;
      throw std::runtime_error(msg.str()) ;
    }
  }
}

/// Run all HexCell plan checks.
std::vector<CaseResult> run_hex_cell_plan_tests(bool write_example_files,
                                                const std::string& output_dir) {
  std::vector<CaseResult> results ;

  const std::vector<PlanCase> hex_cases = hex_plan_cases() ;
  for(size_t i = 0; i < hex_cases.size(); ++i) {
    results.push_back(run_hex_case(hex_cases[i],
                                   write_example_files,
                                   output_dir)) ;
  }

  return results ;
}

/// Run all Prism plan checks.
std::vector<CaseResult> run_prism_cell_plan_tests(bool write_example_files,
                                                  const std::string& output_dir) {
  std::vector<CaseResult> results ;

  const std::vector<PlanCase> prism_cases = prism_plan_cases() ;
  for(size_t i = 0; i < prism_cases.size(); ++i) {
    results.push_back(run_prism_case(prism_cases[i],
                                     write_example_files,
                                     output_dir)) ;
  }

  return results ;
}

/// Run all general tetrahedron plan checks.
std::vector<CaseResult> run_tetra_cell_plan_tests(bool write_example_files,
                                                  const std::string& output_dir) {
  std::vector<CaseResult> results ;

  const std::vector<PlanCase> general_cases = general_plan_cases() ;
  for(size_t i = 0; i < general_cases.size(); ++i) {
    results.push_back(run_general_case("general_tet",
                                       "general tetra",
                                       build_tetra_cell,
                                       general_cases[i],
                                       write_example_files,
                                       output_dir)) ;
  }

  return results ;
}

/// Run all general pyramid plan checks.
std::vector<CaseResult> run_pyramid_cell_plan_tests(bool write_example_files,
                                                    const std::string& output_dir) {
  std::vector<CaseResult> results ;

  const std::vector<PlanCase> general_cases = general_plan_cases() ;
  for(size_t i = 0; i < general_cases.size(); ++i) {
    results.push_back(run_general_case("general_pyramid",
                                       "general pyramid",
                                       build_pyramid_cell,
                                       general_cases[i],
                                       write_example_files,
                                       output_dir)) ;
  }

  return results ;
}

/// Run every cell-family plan check.
std::vector<CaseResult> run_cell_plan_tests(bool write_example_files,
                                            const std::string& output_dir) {
  std::vector<CaseResult> results ;

  std::vector<CaseResult> group =
    run_hex_cell_plan_tests(write_example_files, output_dir) ;
  results.insert(results.end(), group.begin(), group.end()) ;

  group = run_prism_cell_plan_tests(write_example_files, output_dir) ;
  results.insert(results.end(), group.begin(), group.end()) ;

  group = run_tetra_cell_plan_tests(write_example_files, output_dir) ;
  results.insert(results.end(), group.begin(), group.end()) ;

  group = run_pyramid_cell_plan_tests(write_example_files, output_dir) ;
  results.insert(results.end(), group.begin(), group.end()) ;

  return results ;
}
