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

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {

/// Format cell-neighbor index lists as readable integers.
std::string int_vector_string(const std::vector<int32>& values) {
  std::ostringstream out ;
  out << "[" ;
  for(size_t i = 0; i < values.size(); ++i) {
    if(i != 0) {
      out << "," ;
    }
    out << values[i] ;
  }
  out << "]" ;
  return out.str() ;
}

/// Return x, y, or z from a node coordinate.
double coord(const Node* node, int axis) {
  if(axis == 0) {
    return node->p.x ;
  }
  if(axis == 1) {
    return node->p.y ;
  }
  return node->p.z ;
}

/// Collect the eight distinct corner nodes owned by one HexCell leaf.
std::vector<Node*> unique_leaf_nodes(HexCell* leaf) {
  std::map<Node*, bool> seen ;
  std::vector<Node*> nodes ;
  std::vector<Edge*> edges = leaf->get_edges() ;

  for(size_t i = 0; i < edges.size(); ++i) {
    if(edges[i] == 0 || edges[i]->head == 0 || edges[i]->tail == 0) {
      throw std::runtime_error("leaf has an incomplete edge") ;
    }
    if(seen.find(edges[i]->head) == seen.end()) {
      seen[edges[i]->head] = true ;
      nodes.push_back(edges[i]->head) ;
    }
    if(seen.find(edges[i]->tail) == seen.end()) {
      seen[edges[i]->tail] = true ;
      nodes.push_back(edges[i]->tail) ;
    }
  }

  if(nodes.size() != 8) {
    std::ostringstream msg ;
    msg << "expected 8 nodes for a hex leaf, found " << nodes.size() ;
    throw std::runtime_error(msg.str()) ;
  }

  return nodes ;
}

/// Compare coordinates with a small tolerance for corner matching.
bool close_to(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <= 1.0e-10 ;
}

/// Find the node at one min/max corner of a hex leaf.
Node* find_corner(const std::vector<Node*>& nodes,
                  const double min_coord[3],
                  const double max_coord[3],
                  const int corner[3]) {
  for(size_t i = 0; i < nodes.size(); ++i) {
    bool match = true ;
    for(int axis = 0; axis < 3; ++axis) {
      const double target = corner[axis] ? max_coord[axis] : min_coord[axis] ;
      match = match && close_to(coord(nodes[i], axis), target) ;
    }
    if(match) {
      return nodes[i] ;
    }
  }

  throw std::runtime_error("could not identify a VTK hex corner") ;
}

/// Order a HexCell leaf's corner nodes in VTK hexahedron order.
std::vector<Node*> vtk_hex_nodes(HexCell* leaf) {
  std::vector<Node*> nodes = unique_leaf_nodes(leaf) ;
  double min_coord[3] = {coord(nodes[0], 0), coord(nodes[0], 1), coord(nodes[0], 2)} ;
  double max_coord[3] = {min_coord[0], min_coord[1], min_coord[2]} ;

  for(size_t i = 1; i < nodes.size(); ++i) {
    for(int axis = 0; axis < 3; ++axis) {
      min_coord[axis] = std::min(min_coord[axis], coord(nodes[i], axis)) ;
      max_coord[axis] = std::max(max_coord[axis], coord(nodes[i], axis)) ;
    }
  }

  const int corner_bits[8][3] = {
    {0, 0, 0},
    {1, 0, 0},
    {1, 1, 0},
    {0, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {1, 1, 1},
    {0, 1, 1}
  } ;

  std::vector<Node*> ordered(8) ;
  for(int i = 0; i < 8; ++i) {
    ordered[i] = find_corner(nodes, min_coord, max_coord, corner_bits[i]) ;
  }
  return ordered ;
}

} // namespace

/// Write HexCell leaves as a VTK unstructured-grid file.
int write_hex_vtk(const std::string& path,
                  const std::string& title,
                  int root_split_code,
                  const std::vector<HexCell*>& leaves) {
  std::vector<std::vector<Node*> > cell_nodes(leaves.size()) ;
  std::vector<Node*> points ;
  std::map<Node*, int> point_id ;

  for(size_t cell = 0; cell < leaves.size(); ++cell) {
    cell_nodes[cell] = vtk_hex_nodes(leaves[cell]) ;
    for(size_t n = 0; n < cell_nodes[cell].size(); ++n) {
      Node* node = cell_nodes[cell][n] ;
      if(point_id.find(node) == point_id.end()) {
        const int id = int(points.size()) ;
        point_id[node] = id ;
        points.push_back(node) ;
      }
    }
  }

  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open VTK output file") ;
  }

  out << "# vtk DataFile Version 3.0\n" ;
  out << title << "\n" ;
  out << "ASCII\n" ;
  out << "DATASET UNSTRUCTURED_GRID\n" ;
  out << "POINTS " << points.size() << " float\n" ;
  out << std::setprecision(17) ;
  for(size_t i = 0; i < points.size(); ++i) {
    out << points[i]->p.x << " " << points[i]->p.y << " " << points[i]->p.z << "\n" ;
  }

  out << "CELLS " << cell_nodes.size() << " " << cell_nodes.size() * 9 << "\n" ;
  for(size_t cell = 0; cell < cell_nodes.size(); ++cell) {
    out << "8" ;
    for(size_t n = 0; n < cell_nodes[cell].size(); ++n) {
      out << " " << point_id[cell_nodes[cell][n]] ;
    }
    out << "\n" ;
  }

  out << "CELL_TYPES " << cell_nodes.size() << "\n" ;
  for(size_t cell = 0; cell < cell_nodes.size(); ++cell) {
    out << "12\n" ;
  }

  out << "CELL_DATA " << cell_nodes.size() << "\n" ;
  out << "SCALARS root_split_code int 1\n" ;
  out << "LOOKUP_TABLE default\n" ;
  for(size_t cell = 0; cell < cell_nodes.size(); ++cell) {
    out << root_split_code << "\n" ;
  }
  out << "SCALARS fvmadapt_cell_index int 1\n" ;
  out << "LOOKUP_TABLE default\n" ;
  for(size_t cell = 0; cell < leaves.size(); ++cell) {
    out << leaves[cell]->getCellIndex() << "\n" ;
  }

  return int(points.size()) ;
}

/// Collect unique leaf edges from a fixture's root edge list.
std::vector<Edge*> edge_leaves(std::list<Edge*>& root_edges) {
  std::set<Edge*> seen ;
  std::vector<Edge*> leaves ;

  for(std::list<Edge*>::iterator edge = root_edges.begin() ;
      edge != root_edges.end(); ++edge) {
    std::list<Edge*> local_leaves ;
    (*edge)->sort_leaves(local_leaves) ;
    for(std::list<Edge*>::iterator leaf = local_leaves.begin() ;
        leaf != local_leaves.end(); ++leaf) {
      if(seen.insert(*leaf).second) {
        leaves.push_back(*leaf) ;
      }
    }
  }

  return leaves ;
}

/// Return an existing point id or assign a new one.
int point_id(Node* node, std::map<Node*, int>& ids, std::vector<Node*>& points) {
  std::map<Node*, int>::const_iterator existing = ids.find(node) ;
  if(existing != ids.end()) {
    return existing->second ;
  }

  const int id = int(points.size()) ;
  ids[node] = id ;
  points.push_back(node) ;
  return id ;
}

/// Write fixture edge leaves as a VTK wireframe file.
int write_wireframe_vtk(const std::string& path,
                        const std::string& title,
                        std::list<Node*>& root_nodes,
                        std::list<Edge*>& root_edges) {
  std::vector<Node*> points ;
  std::map<Node*, int> ids ;
  for(std::list<Node*>::iterator node = root_nodes.begin() ;
      node != root_nodes.end(); ++node) {
    point_id(*node, ids, points) ;
  }

  const std::vector<Edge*> leaves = edge_leaves(root_edges) ;
  for(size_t i = 0; i < leaves.size(); ++i) {
    point_id(leaves[i]->head, ids, points) ;
    point_id(leaves[i]->tail, ids, points) ;
  }

  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open wireframe VTK output file") ;
  }

  out << "# vtk DataFile Version 3.0\n" ;
  out << title << "\n" ;
  out << "ASCII\n" ;
  out << "DATASET POLYDATA\n" ;
  out << "POINTS " << points.size() << " float\n" ;
  out << std::setprecision(17) ;
  for(size_t i = 0; i < points.size(); ++i) {
    out << points[i]->p.x << " " << points[i]->p.y << " " << points[i]->p.z << "\n" ;
  }

  out << "LINES " << leaves.size() << " " << leaves.size() * 3 << "\n" ;
  for(size_t i = 0; i < leaves.size(); ++i) {
    out << "2 " << ids[leaves[i]->head] << " " << ids[leaves[i]->tail] << "\n" ;
  }

  return int(points.size()) ;
}

namespace {

/// Count leaves in a QuadFace plan without writing geometry.
int count_quad_face_leaves(const std::vector<char>& face_plan) {
  QuadFace face(4) ;
  std::vector<QuadFace*> leaves ;
  face.empty_resplit(face_plan, 0, leaves) ;
  return int(leaves.size()) ;
}

/// Count leaves in an isotropic Face plan without writing geometry.
int count_general_face_leaves(const std::vector<char>& face_plan, int num_edges) {
  if(face_plan.empty()) {
    return 1 ;
  }

  std::vector<char> work ;
  work.push_back(char(0)) ;
  size_t index = 0 ;
  int leaves = 0 ;
  for(size_t cursor = 0; cursor < work.size(); ++cursor) {
    const char code = index < face_plan.size() ? face_plan[index++] : char(0) ;
    if(code == 0) {
      ++leaves ;
    } else if(code == 1) {
      for(int child = 0; child < num_edges; ++child) {
        work.push_back(char(0)) ;
      }
    } else if(code == 8) {
      work.push_back(char(0)) ;
    } else {
      std::ostringstream msg ;
      msg << "unsupported general face plan code " << int(code)
          << " in " << plan_string(face_plan) ;
      throw std::runtime_error(msg.str()) ;
    }
  }
  return leaves ;
}

/// Count leaves in an Edge plan without keeping the generated nodes.
int count_edge_leaves(const std::vector<char>& edge_plan) {
  Node head(vect3d(0.0, 0.0, 0.0)) ;
  Node tail(vect3d(1.0, 0.0, 0.0)) ;
  Edge edge(&head, &tail) ;
  std::list<Node*> nodes ;
  edge.resplit(edge_plan, nodes) ;
  std::list<Edge*> leaves ;
  edge.sort_leaves(leaves) ;
  const int count = int(leaves.size()) ;
  for(std::list<Node*>::iterator node = nodes.begin(); node != nodes.end(); ++node) {
    delete *node ;
  }
  nodes.clear() ;
  return count ;
}

/// Write the HexCell section of core_plan_audit.dat for one plan case.
/// This is a human-inspection artifact: it shows the split plan induced on
/// each parent face, the fine-cell-side c1 index for each face leaf, and the
/// edge plans induced by each face plan.
void write_hex_plan_audit(std::ofstream& out, const CaseResult& result) {
  out << "case " << result.shape << " " << result.case_name << "\n" ;
  out << "  input_cell_plan " << plan_string(result.input_plan) << "\n" ;
  out << "  canonical_cell_plan " << plan_string(result.canonical_plan) << "\n" ;
  out << "  counted_leaves " << result.counted_leaves << "\n" ;
  if(!result.vtk_file.empty()) {
    out << "  vtk_file " << result.vtk_file << "\n" ;
  }
  for(int face = 0; face < 6; ++face) {
    const std::vector<char> face_plan =
      extract_hex_face(result.input_plan, DIRECTION(face)) ;
    const std::vector<int32> c1 =
      get_c1_hex(result.input_plan, face_plan, 0, char(face)) ;
    out << "  face " << face
        << " plan " << plan_string(face_plan)
        << " leaves " << count_quad_face_leaves(face_plan)
        << " c1 " << int_vector_string(c1) << "\n" ;
    for(int edge = 0; edge < 4; ++edge) {
      std::vector<char> edge_plan ;
      extract_quad_edge(face_plan, edge_plan, unsigned(edge)) ;
      out << "    edge " << edge
          << " plan " << plan_string(edge_plan)
          << " leaves " << count_edge_leaves(edge_plan) << "\n" ;
    }
  }
}

/// Write face and c1 projections for one Prism plan case.
void write_prism_plan_audit(std::ofstream& out, const CaseResult& result) {
  out << "case " << result.shape << " " << result.case_name << "\n" ;
  out << "  input_cell_plan " << plan_string(result.input_plan) << "\n" ;
  out << "  canonical_cell_plan " << plan_string(result.canonical_plan) << "\n" ;
  out << "  counted_leaves " << result.counted_leaves << "\n" ;
  if(!result.vtk_file.empty()) {
    out << "  vtk_file " << result.vtk_file << "\n" ;
  }
  for(int face = 0; face < 5; ++face) {
    const std::vector<char> face_plan = extract_prism_face(result.input_plan, face) ;
    const std::vector<int32> c1 =
      get_c1_prism(result.input_plan, face_plan, 0, face) ;
    out << "  face " << face
        << " plan " << plan_string(face_plan)
        << " leaves "
        << (face < 2 ? count_general_face_leaves(face_plan, 3)
                     : count_quad_face_leaves(face_plan))
        << " c1 " << int_vector_string(c1) << "\n" ;
  }
}

/// Write the plan-level audit rows for one general Cell case.
void write_general_plan_audit(std::ofstream& out, const CaseResult& result) {
  out << "case " << result.shape << " " << result.case_name << "\n" ;
  out << "  input_cell_plan " << plan_string(result.input_plan) << "\n" ;
  out << "  canonical_cell_plan " << plan_string(result.canonical_plan) << "\n" ;
  out << "  counted_leaves " << result.counted_leaves << "\n" ;
  if(!result.vtk_file.empty()) {
    out << "  vtk_file " << result.vtk_file << "\n" ;
  }
}

/// Write one summary row per plan case.
void write_summary(const std::string& output_dir,
                   const std::vector<CaseResult>& results) {
  const std::string path = join_path(output_dir, "core_split_summary.dat") ;
  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open summary output file") ;
  }

  out << "# shape case input_plan canonical_plan counted_leaves "
      << "resplit_leaves sorted_leaves vtk_points vtk_file\n" ;
  for(size_t i = 0; i < results.size(); ++i) {
    out << results[i].shape << " "
        << results[i].case_name << " "
        << plan_string(results[i].input_plan) << " "
        << plan_string(results[i].canonical_plan) << " "
        << results[i].counted_leaves << " "
        << results[i].resplit_leaves << " "
        << results[i].sorted_leaves << " "
        << results[i].point_count << " "
        << results[i].vtk_file << "\n" ;
  }
}

/// Write the detailed plan projection audit.
void write_plan_audit(const std::string& output_dir,
                      const std::vector<CaseResult>& results) {
  const std::string path = join_path(output_dir, "core_plan_audit.dat") ;
  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open plan audit output file") ;
  }

  out << "# FVMAdapt core plan audit\n" ;
  out << "# Plans are char vectors printed as integer codes in breadth-first order.\n" ;
  for(size_t i = 0; i < results.size(); ++i) {
    if(results[i].shape == "hex") {
      write_hex_plan_audit(out, results[i]) ;
    } else if(results[i].shape == "prism") {
      write_prism_plan_audit(out, results[i]) ;
    } else {
      write_general_plan_audit(out, results[i]) ;
    }
    out << "\n" ;
  }
}

} // namespace

/// Write all optional human-inspection artifacts for the core tests.
void write_example_artifacts(const std::string& output_dir,
                             const std::vector<CaseResult>& results) {
  write_summary(output_dir, results) ;
  write_plan_audit(output_dir, results) ;
  write_plan_replay(output_dir, results) ;
  for(size_t i = 0; i < results.size(); ++i) {
    if(!results[i].vtk_file.empty()) {
      std::cout << "wrote " << results[i].vtk_file << "\n" ;
    }
  }
  std::cout << "wrote " << join_path(output_dir, "core_split_summary.dat") << "\n" ;
  std::cout << "wrote " << join_path(output_dir, "core_plan_audit.dat") << "\n" ;
  std::cout << "wrote " << join_path(output_dir, "core_plan_replay.dat") << "\n" ;
}
