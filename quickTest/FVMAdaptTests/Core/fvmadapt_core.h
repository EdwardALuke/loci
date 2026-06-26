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

#ifndef FVMADAPT_CORE_H
#define FVMADAPT_CORE_H

#include <FVMAdapt/diamondcell.h>
#include <FVMAdapt/hexcell.h>
#include <FVMAdapt/prism.h>

#include <fstream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/// FVMAdapt library functions used directly by these core tests.
void extract_quad_edge(const std::vector<char>& facePlan,
                       std::vector<char>& edgePlan,
                       unsigned int dd) ;
std::vector<char> extract_prism_face(const std::vector<char>& cellPlan, int dd) ;
std::vector<char> transfer_plan_g2q(std::vector<char>& facePlan) ;
std::vector<char> transfer_plan_q2g(const std::vector<char>& facePlan) ;

/// Owns the unit hexahedron used by the HexCell tests.
struct HexFixture {
  HexFixture() ;
  ~HexFixture() ;

  std::list<Node*> nodes ;
  std::list<Edge*> edges ;
  std::list<QuadFace*> faces ;
  HexCell* root ;
} ;

/// Owns the unit triangular prism used by the Prism tests.
struct PrismFixture {
  PrismFixture() ;
  ~PrismFixture() ;

  std::list<Node*> nodes ;
  std::list<Edge*> edges ;
  std::list<QuadFace*> qfaces ;
  std::list<Face*> gfaces ;
  Prism* root ;
} ;

/// Owns one general Cell shape used by the Cell/DiamondCell tests.
struct GeneralFixture {
  GeneralFixture() ;
  ~GeneralFixture() ;

  std::list<Node*> nodes ;
  std::list<Edge*> edges ;
  std::list<Face*> faces ;
  Cell* root ;
} ;

typedef void (*GeneralBuilder)(GeneralFixture&) ;

/// Records the observable result of replaying one cell plan.
struct CaseResult {
  std::string shape ;
  std::string case_name ;
  std::vector<char> input_plan ;
  std::vector<char> canonical_plan ;
  int counted_leaves ;
  int resplit_leaves ;
  int sorted_leaves ;
  int point_count ;
  std::string vtk_file ;
} ;

/// Build a unit hexahedron test topology.
void build_unit_hex(HexFixture& fixture) ;

/// Build a unit triangular-prism test topology.
void build_unit_prism(PrismFixture& fixture) ;

/// Build a tetrahedron as a general Cell test topology.
void build_tetra_cell(GeneralFixture& fixture) ;

/// Build a pyramid as a general Cell test topology.
void build_pyramid_cell(GeneralFixture& fixture) ;

/// Format a split plan as a compact integer list.
std::string plan_string(const std::vector<char>& plan) ;

/// Create one output directory if it does not already exist.
void make_directory(const std::string& path) ;

/// Join a directory and file name without normalizing the path.
std::string join_path(const std::string& dir, const std::string& file) ;

/// Check fixed general-face and quad-face transfer round trips.
void check_plan_transfer_round_trips() ;

/// Replay the HexCell plan cases and return their observed counts.
std::vector<CaseResult> run_hex_cell_plan_tests(bool write_example_files,
                                                const std::string& output_dir) ;

/// Replay the Prism plan cases and return their observed counts.
std::vector<CaseResult> run_prism_cell_plan_tests(bool write_example_files,
                                                  const std::string& output_dir) ;

/// Replay the general tetrahedron plan cases and return their observed counts.
std::vector<CaseResult> run_tetra_cell_plan_tests(bool write_example_files,
                                                  const std::string& output_dir) ;

/// Replay the general pyramid plan cases and return their observed counts.
std::vector<CaseResult> run_pyramid_cell_plan_tests(bool write_example_files,
                                                    const std::string& output_dir) ;

/// Replay all cell plan cases and return their observed counts.
std::vector<CaseResult> run_cell_plan_tests(bool write_example_files,
                                            const std::string& output_dir) ;

/// Check replay-trace leaf counts against the plan check results.
void check_plan_replay_counts(const std::vector<CaseResult>& results) ;

/// Write a simple VTK file for HexCell leaves.
int write_hex_vtk(const std::string& path,
                  const std::string& title,
                  int root_split_code,
                  const std::vector<HexCell*>& leaves) ;

/// Write a VTK wireframe for the owned fixture topology.
int write_wireframe_vtk(const std::string& path,
                        const std::string& title,
                        std::list<Node*>& root_nodes,
                        std::list<Edge*>& root_edges) ;

/// Write the breadth-first replay trace artifact.
void write_plan_replay(const std::string& output_dir,
                       const std::vector<CaseResult>& results) ;

/// Write all optional human-inspection artifacts for the plan cases.
void write_example_artifacts(const std::string& output_dir,
                             const std::vector<CaseResult>& results) ;

#endif
