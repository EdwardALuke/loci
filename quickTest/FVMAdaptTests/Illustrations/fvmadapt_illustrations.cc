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

#include "../Core/fvmadapt_core.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <sys/stat.h>

namespace {

struct Point {
  Point(): x(0.0), y(0.0), z(0.0) {}
  Point(double x_in, double y_in, double z_in): x(x_in), y(y_in), z(z_in) {}

  double x ;
  double y ;
  double z ;
} ;

struct HexBlock {
  HexBlock(): root(0) {}
  ~HexBlock() {
    delete root ;
    cleanup_list(nodes, edges, faces) ;
  }

  std::list<Node*> nodes ;
  std::list<Edge*> edges ;
  std::list<QuadFace*> faces ;
  HexCell* root ;
} ;

struct PrismBlock {
  PrismBlock(): root(0) {}
  ~PrismBlock() {
    delete root ;
    cleanup_list(nodes, edges) ;
    cleanup_list(qfaces) ;
    cleanup_list(gfaces) ;
  }

  std::list<Node*> nodes ;
  std::list<Edge*> edges ;
  std::list<QuadFace*> qfaces ;
  std::list<Face*> gfaces ;
  Prism* root ;
} ;

struct HexView {
  HexView(HexCell* c, int source, int code): cell(c), source_cell(source),
                                             root_split_code(code) {}

  HexCell* cell ;
  int source_cell ;
  int root_split_code ;
} ;

struct QuadPatch {
  QuadPatch(): u0(0.0), u1(1.0), v0(0.0), v1(1.0) {}
  QuadPatch(double u0_in, double u1_in, double v0_in, double v1_in):
            u0(u0_in), u1(u1_in), v0(v0_in), v1(v1_in) {}

  double u0 ;
  double u1 ;
  double v0 ;
  double v1 ;
} ;

std::vector<char> plan(std::initializer_list<int> codes) {
  std::vector<char> result ;
  for(std::initializer_list<int>::const_iterator code = codes.begin() ;
      code != codes.end(); ++code) {
    result.push_back(char(*code)) ;
  }
  return result ;
}

std::string plan_string(const std::vector<char>& codes) {
  std::ostringstream out ;
  out << "[" ;
  for(size_t i = 0; i < codes.size(); ++i) {
    if(i != 0) {
      out << "," ;
    }
    out << int(codes[i]) ;
  }
  out << "]" ;
  return out.str() ;
}

void make_directory(const std::string& path) {
  if(path.empty()) {
    return ;
  }
  if(mkdir(path.c_str(), 0775) != 0 && errno != EEXIST) {
    std::ostringstream msg ;
    msg << "could not create directory '" << path << "': "
        << std::strerror(errno) ;
    throw std::runtime_error(msg.str()) ;
  }
}

void make_directories(const std::string& path) {
  if(path.empty()) {
    return ;
  }

  std::string current ;
  size_t start = 0 ;
  if(path[0] == '/') {
    current = "/" ;
    start = 1 ;
  }

  while(start <= path.size()) {
    const size_t slash = path.find('/', start) ;
    const std::string part =
      path.substr(start,
                  slash == std::string::npos ?
                  std::string::npos : slash - start) ;
    if(!part.empty() && part != ".") {
      if(!current.empty() && current[current.size() - 1] != '/') {
        current += "/" ;
      }
      current += part ;
      make_directory(current) ;
    }
    if(slash == std::string::npos) {
      break ;
    }
    start = slash + 1 ;
  }
}

std::string join_path(const std::string& dir, const std::string& file) {
  if(dir.empty()) {
    return file ;
  }
  if(dir[dir.size() - 1] == '/') {
    return dir + file ;
  }
  return dir + "/" + file ;
}

std::string output_group_dir(const std::string& output_dir,
                             const std::string& group) {
  const std::string group_dir = join_path(output_dir, group) ;
  make_directories(group_dir) ;
  return group_dir ;
}

Point mix(const Point& p00, const Point& p10,
          const Point& p11, const Point& p01,
          double u, double v) {
  const double w00 = (1.0 - u) * (1.0 - v) ;
  const double w10 = u * (1.0 - v) ;
  const double w11 = u * v ;
  const double w01 = (1.0 - u) * v ;
  return Point(w00 * p00.x + w10 * p10.x + w11 * p11.x + w01 * p01.x,
               w00 * p00.y + w10 * p10.y + w11 * p11.y + w01 * p01.y,
               w00 * p00.z + w10 * p10.z + w11 * p11.z + w01 * p01.z) ;
}

Node* add_node(std::list<Node*>& nodes, double x, double y, double z) {
  vect3d p(x, y, z) ;
  Node* node = new Node(p, int32(nodes.size() + 1)) ;
  nodes.push_back(node) ;
  return node ;
}

Edge* add_edge(std::list<Edge*>& edges, Node* head, Node* tail) {
  Edge* edge = new Edge(head, tail) ;
  edges.push_back(edge) ;
  return edge ;
}

void build_hex_box(HexBlock& block,
                   double x0, double x1,
                   double y0, double y1,
                   double z0, double z1) {
  Node* node[8] ;
  node[0] = add_node(block.nodes, x0, y0, z0) ;
  node[1] = add_node(block.nodes, x0, y0, z1) ;
  node[2] = add_node(block.nodes, x0, y1, z0) ;
  node[3] = add_node(block.nodes, x0, y1, z1) ;
  node[4] = add_node(block.nodes, x1, y0, z0) ;
  node[5] = add_node(block.nodes, x1, y0, z1) ;
  node[6] = add_node(block.nodes, x1, y1, z0) ;
  node[7] = add_node(block.nodes, x1, y1, z1) ;

  const int edge_head[12] = {0, 1, 2, 3, 0, 1, 4, 5, 0, 2, 4, 6} ;
  const int edge_tail[12] = {4, 5, 6, 7, 2, 3, 6, 7, 1, 3, 5, 7} ;

  Edge* edge[12] ;
  for(int i = 0; i < 12; ++i) {
    edge[i] = add_edge(block.edges, node[edge_head[i]], node[edge_tail[i]]) ;
  }

  const int face_to_edge[6][4] = {
    {6, 11, 7, 10},
    {4, 9, 5, 8},
    {2, 11, 3, 9},
    {0, 10, 1, 8},
    {1, 7, 3, 5},
    {0, 6, 2, 4}
  } ;

  QuadFace** face = new QuadFace*[6] ;
  for(int i = 0; i < 6; ++i) {
    face[i] = new QuadFace(4) ;
    block.faces.push_back(face[i]) ;
    for(int j = 0; j < 4; ++j) {
      face[i]->edge[j] = edge[face_to_edge[i][j]] ;
    }
  }

  block.root = new HexCell(face) ;
}

void build_prism(PrismBlock& block, const Point base[3], double height) {
  Node* node[6] ;
  for(int i = 0; i < 3; ++i) {
    node[i] = add_node(block.nodes, base[i].x, base[i].y, base[i].z) ;
    node[i + 3] = add_node(block.nodes, base[i].x, base[i].y,
                           base[i].z + height) ;
  }

  const int edge_head[9] = {0, 1, 2, 3, 4, 5, 0, 1, 2} ;
  const int edge_tail[9] = {1, 2, 0, 4, 5, 3, 3, 4, 5} ;

  Edge* edge[9] ;
  for(int i = 0; i < 9; ++i) {
    edge[i] = add_edge(block.edges, node[edge_head[i]], node[edge_tail[i]]) ;
  }

  const int gface_to_edge[2][3] = {{0, 1, 2}, {3, 4, 5}} ;
  const int qface_to_edge[3][4] = {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}} ;

  Prism* prism = new Prism(3) ;
  for(int i = 0; i < 2; ++i) {
    Face* face = new Face(3) ;
    block.gfaces.push_back(face) ;
    for(int j = 0; j < 3; ++j) {
      face->edge[j] = edge[gface_to_edge[i][j]] ;
      face->needReverse[j] = false ;
    }
    prism->setFace(i, face) ;
  }

  for(int i = 0; i < 3; ++i) {
    QuadFace* face = new QuadFace(4) ;
    block.qfaces.push_back(face) ;
    for(int j = 0; j < 4; ++j) {
      face->edge[j] = edge[qface_to_edge[i][j]] ;
    }
    prism->setFace(i, face) ;
  }

  block.root = prism ;
}

bool close_to(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <= 1.0e-10 ;
}

double coord(const Node* node, int axis) {
  if(axis == 0) {
    return node->p.x ;
  }
  if(axis == 1) {
    return node->p.y ;
  }
  return node->p.z ;
}

std::vector<Node*> unique_hex_nodes(HexCell* leaf) {
  std::map<Node*, bool> seen ;
  std::vector<Node*> nodes ;
  std::vector<Edge*> edges = leaf->get_edges() ;

  for(size_t i = 0; i < edges.size(); ++i) {
    if(edges[i]->head == 0 || edges[i]->tail == 0) {
      throw std::runtime_error("hex leaf has an incomplete edge") ;
    }
    if(!seen[edges[i]->head]) {
      seen[edges[i]->head] = true ;
      nodes.push_back(edges[i]->head) ;
    }
    if(!seen[edges[i]->tail]) {
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

Node* find_corner(const std::vector<Node*>& nodes,
                  const double min_coord[3],
                  const double max_coord[3],
                  const bool corner[3]) {
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

std::vector<Node*> vtk_hex_nodes(HexCell* leaf) {
  std::vector<Node*> nodes = unique_hex_nodes(leaf) ;
  double min_coord[3] = {coord(nodes[0], 0), coord(nodes[0], 1),
                         coord(nodes[0], 2)} ;
  double max_coord[3] = {min_coord[0], min_coord[1], min_coord[2]} ;

  for(size_t i = 1; i < nodes.size(); ++i) {
    for(int axis = 0; axis < 3; ++axis) {
      min_coord[axis] = std::min(min_coord[axis], coord(nodes[i], axis)) ;
      max_coord[axis] = std::max(max_coord[axis], coord(nodes[i], axis)) ;
    }
  }

  const bool corner_bits[8][3] = {
    {false, false, false},
    {true, false, false},
    {true, true, false},
    {false, true, false},
    {false, false, true},
    {true, false, true},
    {true, true, true},
    {false, true, true}
  } ;

  std::vector<Node*> ordered(8) ;
  for(int i = 0; i < 8; ++i) {
    ordered[i] = find_corner(nodes, min_coord, max_coord, corner_bits[i]) ;
  }
  return ordered ;
}

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

void write_hex_vtk(const std::string& path,
                   const std::string& title,
                   const std::vector<HexView>& cells) {
  std::vector<std::vector<Node*> > cell_nodes(cells.size()) ;
  std::vector<Node*> points ;
  std::map<Node*, int> ids ;

  for(size_t cell = 0; cell < cells.size(); ++cell) {
    cell_nodes[cell] = vtk_hex_nodes(cells[cell].cell) ;
    for(size_t n = 0; n < cell_nodes[cell].size(); ++n) {
      point_id(cell_nodes[cell][n], ids, points) ;
    }
  }

  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open hex VTK output file") ;
  }

  out << "# vtk DataFile Version 3.0\n" ;
  out << title << "\n" ;
  out << "ASCII\n" ;
  out << "DATASET UNSTRUCTURED_GRID\n" ;
  out << "POINTS " << points.size() << " float\n" ;
  out << std::setprecision(17) ;
  for(size_t i = 0; i < points.size(); ++i) {
    out << points[i]->p.x << " " << points[i]->p.y << " "
        << points[i]->p.z << "\n" ;
  }

  out << "CELLS " << cell_nodes.size() << " " << cell_nodes.size() * 9 << "\n" ;
  for(size_t cell = 0; cell < cell_nodes.size(); ++cell) {
    out << "8" ;
    for(size_t n = 0; n < cell_nodes[cell].size(); ++n) {
      out << " " << ids[cell_nodes[cell][n]] ;
    }
    out << "\n" ;
  }

  out << "CELL_TYPES " << cell_nodes.size() << "\n" ;
  for(size_t cell = 0; cell < cell_nodes.size(); ++cell) {
    out << "12\n" ;
  }

  out << "CELL_DATA " << cell_nodes.size() << "\n" ;
  out << "SCALARS source_cell int 1\n" ;
  out << "LOOKUP_TABLE default\n" ;
  for(size_t cell = 0; cell < cells.size(); ++cell) {
    out << cells[cell].source_cell << "\n" ;
  }
  out << "SCALARS root_split_code int 1\n" ;
  out << "LOOKUP_TABLE default\n" ;
  for(size_t cell = 0; cell < cells.size(); ++cell) {
    out << cells[cell].root_split_code << "\n" ;
  }
  out << "SCALARS fvmadapt_cell_index int 1\n" ;
  out << "LOOKUP_TABLE default\n" ;
  for(size_t cell = 0; cell < cells.size(); ++cell) {
    out << cells[cell].cell->getCellIndex() << "\n" ;
  }
}

void collect_edge_leaves(const std::list<Edge*>& root_edges,
                         std::vector<Edge*>& leaves) {
  std::set<Edge*> seen ;
  for(std::list<Edge*>::const_iterator edge = root_edges.begin();
      edge != root_edges.end(); ++edge) {
    std::list<Edge*> local_leaves ;
    (*edge)->sort_leaves(local_leaves) ;
    for(std::list<Edge*>::iterator leaf = local_leaves.begin();
        leaf != local_leaves.end(); ++leaf) {
      if(seen.insert(*leaf).second) {
        leaves.push_back(*leaf) ;
      }
    }
  }
}

void add_wire_inputs(HexBlock& block,
                     std::vector<std::list<Node*>*>& node_lists,
                     std::vector<std::list<Edge*>*>& edge_lists) {
  node_lists.push_back(&block.nodes) ;
  edge_lists.push_back(&block.edges) ;
}

void collect_nodes(const std::list<Node*>& roots, std::vector<Node*>& points,
                   std::map<Node*, int>& ids) {
  for(std::list<Node*>::const_iterator node = roots.begin();
      node != roots.end(); ++node) {
    point_id(*node, ids, points) ;
  }
}

void write_wireframe_vtk(const std::string& path,
                         const std::string& title,
                         const std::vector<std::list<Node*>*>& node_lists,
                         const std::vector<std::list<Edge*>*>& edge_lists) {
  std::vector<Node*> points ;
  std::map<Node*, int> ids ;
  std::vector<Edge*> leaves ;

  for(size_t i = 0; i < node_lists.size(); ++i) {
    collect_nodes(*node_lists[i], points, ids) ;
  }
  for(size_t i = 0; i < edge_lists.size(); ++i) {
    collect_edge_leaves(*edge_lists[i], leaves) ;
  }
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
    out << points[i]->p.x << " " << points[i]->p.y << " "
        << points[i]->p.z << "\n" ;
  }

  out << "LINES " << leaves.size() << " " << leaves.size() * 3 << "\n" ;
  for(size_t i = 0; i < leaves.size(); ++i) {
    out << "2 " << ids[leaves[i]->head] << " "
        << ids[leaves[i]->tail] << "\n" ;
  }
}

void replay_quad_face_plan(const std::vector<char>& face_plan,
                           std::vector<QuadPatch>& leaves) {
  std::vector<QuadPatch> queue ;
  queue.push_back(QuadPatch()) ;

  size_t plan_index = 0 ;
  for(size_t cursor = 0; cursor < queue.size(); ++cursor) {
    const char code = plan_index < face_plan.size() ?
      face_plan[plan_index++] : char(0) ;
    const QuadPatch patch = queue[cursor] ;
    const double um = 0.5 * (patch.u0 + patch.u1) ;
    const double vm = 0.5 * (patch.v0 + patch.v1) ;

    if(code == 0) {
      leaves.push_back(patch) ;
    } else if(code == 1) {
      queue.push_back(QuadPatch(patch.u0, um, patch.v0, patch.v1)) ;
      queue.push_back(QuadPatch(um, patch.u1, patch.v0, patch.v1)) ;
    } else if(code == 2) {
      queue.push_back(QuadPatch(patch.u0, patch.u1, patch.v0, vm)) ;
      queue.push_back(QuadPatch(patch.u0, patch.u1, vm, patch.v1)) ;
    } else if(code == 3) {
      queue.push_back(QuadPatch(patch.u0, um, patch.v0, vm)) ;
      queue.push_back(QuadPatch(um, patch.u1, patch.v0, vm)) ;
      queue.push_back(QuadPatch(patch.u0, um, vm, patch.v1)) ;
      queue.push_back(QuadPatch(um, patch.u1, vm, patch.v1)) ;
    } else {
      std::ostringstream msg ;
      msg << "unsupported quad face split code " << int(code)
          << " in " << plan_string(face_plan) ;
      throw std::runtime_error(msg.str()) ;
    }
  }
}

void write_quad_face_vtk(const std::string& path,
                         const std::string& title,
                         const Point& p00,
                         const Point& p10,
                         const Point& p11,
                         const Point& p01,
                         const std::vector<char>& face_plan) {
  std::vector<QuadPatch> leaves ;
  replay_quad_face_plan(face_plan, leaves) ;
  std::vector<Point> points ;
  std::map<std::pair<int, int>, int> ids ;
  std::vector<std::vector<int> > quads ;

  for(size_t i = 0; i < leaves.size(); ++i) {
    const double us[2] = {leaves[i].u0, leaves[i].u1} ;
    const double vs[2] = {leaves[i].v0, leaves[i].v1} ;
    std::vector<int> quad ;
    const int corner[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}} ;
    for(int c = 0; c < 4; ++c) {
      const int ukey = int(std::floor(us[corner[c][0]] * 1048576.0 + 0.5)) ;
      const int vkey = int(std::floor(vs[corner[c][1]] * 1048576.0 + 0.5)) ;
      const std::pair<int, int> key = std::make_pair(ukey, vkey) ;
      std::map<std::pair<int, int>, int>::const_iterator existing =
        ids.find(key) ;
      if(existing == ids.end()) {
        const int id = int(points.size()) ;
        ids[key] = id ;
        points.push_back(mix(p00, p10, p11, p01,
                             us[corner[c][0]], vs[corner[c][1]])) ;
        quad.push_back(id) ;
      } else {
        quad.push_back(existing->second) ;
      }
    }
    quads.push_back(quad) ;
  }

  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open shared-face VTK output file") ;
  }

  out << "# vtk DataFile Version 3.0\n" ;
  out << title << "\n" ;
  out << "ASCII\n" ;
  out << "DATASET UNSTRUCTURED_GRID\n" ;
  out << "POINTS " << points.size() << " float\n" ;
  out << std::setprecision(17) ;
  for(size_t i = 0; i < points.size(); ++i) {
    out << points[i].x << " " << points[i].y << " " << points[i].z << "\n" ;
  }

  out << "CELLS " << quads.size() << " " << quads.size() * 5 << "\n" ;
  for(size_t i = 0; i < quads.size(); ++i) {
    out << "4 " << quads[i][0] << " " << quads[i][1] << " "
        << quads[i][2] << " " << quads[i][3] << "\n" ;
  }

  out << "CELL_TYPES " << quads.size() << "\n" ;
  for(size_t i = 0; i < quads.size(); ++i) {
    out << "9\n" ;
  }
}

void write_plan_note(const std::string& path,
                     const std::string& case_name,
                     const std::vector<char>& left_cell_plan,
                     const std::vector<char>& shared_face_plan,
                     const std::string& note) {
  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open plan note output file") ;
  }
  out << "case " << case_name << "\n" ;
  out << "left_cell_plan " << plan_string(left_cell_plan) << "\n" ;
  out << "shared_face_plan " << plan_string(shared_face_plan) << "\n" ;
  out << "note " << note << "\n" ;
}

void add_manifest_case(std::ofstream& manifest,
                       const std::string& group,
                       const std::string& name,
                       const std::string& description,
                       const std::vector<std::string>& files) {
  manifest << "case " << group << "/" << name << "\n" ;
  manifest << "  " << description << "\n" ;
  for(size_t i = 0; i < files.size(); ++i) {
    manifest << "  file " << files[i] << "\n" ;
  }
  manifest << "\n" ;
}

void write_hex_pair_case(const std::string& output_dir,
                         std::ofstream& manifest,
                         const std::string& group,
                         const std::string& name,
                         const std::vector<char>& left_plan,
                         int split_code,
                         double z_extent) {
  const std::string group_dir = output_group_dir(output_dir, group) ;

  const std::string before = join_path(group_dir, name + "_before.vtk") ;
  const std::string after = join_path(group_dir, name + "_after.vtk") ;
  const std::string before_wire = join_path(group_dir, name + "_before_wire.vtk") ;
  const std::string after_wire = join_path(group_dir, name + "_after_wire.vtk") ;
  const std::string face = join_path(group_dir, name + "_shared_face.vtk") ;
  const std::string note = join_path(group_dir, name + "_plans.dat") ;

  {
    HexBlock left ;
    HexBlock right ;
    build_hex_box(left, 0.0, 1.0, 0.0, 1.0, 0.0, z_extent) ;
    build_hex_box(right, 1.0, 2.0, 0.0, 1.0, 0.0, z_extent) ;
    left.root->empty_resplit(std::vector<char>()) ;
    right.root->empty_resplit(std::vector<char>()) ;
    std::vector<HexView> cells ;
    cells.push_back(HexView(left.root, 0, 0)) ;
    cells.push_back(HexView(right.root, 1, 0)) ;
    write_hex_vtk(before, "before " + name, cells) ;

    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    add_wire_inputs(left, nodes, edges) ;
    add_wire_inputs(right, nodes, edges) ;
    write_wireframe_vtk(before_wire, "before wireframe " + name, nodes, edges) ;
  }

  {
    HexBlock left ;
    HexBlock right ;
    build_hex_box(left, 0.0, 1.0, 0.0, 1.0, 0.0, z_extent) ;
    build_hex_box(right, 1.0, 2.0, 0.0, 1.0, 0.0, z_extent) ;
    std::vector<HexCell*> left_leaves ;
    left.root->resplit(left_plan, left.nodes, left.edges, left.faces,
                       left_leaves) ;
    right.root->empty_resplit(std::vector<char>()) ;

    std::vector<HexView> cells ;
    for(size_t i = 0; i < left_leaves.size(); ++i) {
      cells.push_back(HexView(left_leaves[i], 0, split_code)) ;
    }
    cells.push_back(HexView(right.root, 1, 0)) ;
    write_hex_vtk(after, "after " + name, cells) ;

    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    add_wire_inputs(left, nodes, edges) ;
    add_wire_inputs(right, nodes, edges) ;
    write_wireframe_vtk(after_wire, "after wireframe " + name, nodes, edges) ;
  }

  std::vector<char> shared_face_plan = extract_hex_face(left_plan, RIGHT) ;
  write_quad_face_vtk(face, "shared face " + name,
                      Point(1.0, 0.0, 0.0),
                      Point(1.0, 1.0, 0.0),
                      Point(1.0, 1.0, z_extent),
                      Point(1.0, 0.0, z_extent),
                      shared_face_plan) ;
  write_plan_note(note, name, left_plan, shared_face_plan,
                  "The right cell remains coarse; this file shows the split "
                  "pattern induced on the original shared face.") ;

  std::vector<std::string> files ;
  files.push_back(join_path(group, name + "_before.vtk")) ;
  files.push_back(join_path(group, name + "_after.vtk")) ;
  files.push_back(join_path(group, name + "_before_wire.vtk")) ;
  files.push_back(join_path(group, name + "_after_wire.vtk")) ;
  files.push_back(join_path(group, name + "_shared_face.vtk")) ;
  files.push_back(join_path(group, name + "_plans.dat")) ;
  add_manifest_case(manifest, group, name,
                    "Two neighboring hex cells with the left cell refined.",
                    files) ;
}

void build_hex_grid(std::vector<HexBlock*>& blocks,
                    int nx,
                    int ny,
                    double z_extent) {
  for(int j = 0; j < ny; ++j) {
    for(int i = 0; i < nx; ++i) {
      HexBlock* block = new HexBlock ;
      build_hex_box(*block,
                    double(i),
                    double(i + 1),
                    double(j),
                    double(j + 1),
                    0.0,
                    z_extent) ;
      blocks.push_back(block) ;
    }
  }
}

void delete_hex_grid(std::vector<HexBlock*>& blocks) {
  for(size_t i = 0; i < blocks.size(); ++i) {
    delete blocks[i] ;
  }
  blocks.clear() ;
}

void add_hex_grid_wire_inputs(std::vector<HexBlock*>& blocks,
                              std::vector<std::list<Node*>*>& nodes,
                              std::vector<std::list<Edge*>*>& edges) {
  for(size_t i = 0; i < blocks.size(); ++i) {
    add_wire_inputs(*blocks[i], nodes, edges) ;
  }
}

void write_hex_grid_plan_note(const std::string& path,
                              const std::string& case_name,
                              const std::vector<char>& center_plan,
                              const std::vector<char>& right_plan,
                              const std::vector<char>& left_plan,
                              const std::vector<char>& front_plan,
                              const std::vector<char>& back_plan,
                              const std::string& note) {
  std::ofstream out(path.c_str()) ;
  if(!out) {
    throw std::runtime_error("could not open 5x5 plan note output file") ;
  }
  out << "case " << case_name << "\n" ;
  out << "center_cell ij 2 2\n" ;
  out << "center_cell_plan " << plan_string(center_plan) << "\n" ;
  out << "right_face_plan " << plan_string(right_plan) << "\n" ;
  out << "left_face_plan " << plan_string(left_plan) << "\n" ;
  out << "front_face_plan " << plan_string(front_plan) << "\n" ;
  out << "back_face_plan " << plan_string(back_plan) << "\n" ;
  out << "note " << note << "\n" ;
}

void write_hex_grid_focus_case(const std::string& output_dir,
                               std::ofstream& manifest) {
  const std::string group = "2d/hex" ;
  const std::string name = "grid_5x5_refine_center_twice" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;

  const int nx = 5 ;
  const int ny = 5 ;
  const int center = 2 + nx * 2 ;
  const double z_extent = 0.08 ;
  const std::vector<char> center_plan = plan({6, 6, 6, 6, 6}) ;

  const std::string before = join_path(group_dir, name + "_before.vtk") ;
  const std::string after = join_path(group_dir, name + "_after.vtk") ;
  const std::string before_wire = join_path(group_dir, name + "_before_wire.vtk") ;
  const std::string after_wire = join_path(group_dir, name + "_after_wire.vtk") ;
  const std::string right_face = join_path(group_dir, name + "_center_right_face.vtk") ;
  const std::string left_face = join_path(group_dir, name + "_center_left_face.vtk") ;
  const std::string front_face = join_path(group_dir, name + "_center_front_face.vtk") ;
  const std::string back_face = join_path(group_dir, name + "_center_back_face.vtk") ;
  const std::string note = join_path(group_dir, name + "_plans.dat") ;

  {
    std::vector<HexBlock*> blocks ;
    build_hex_grid(blocks, nx, ny, z_extent) ;

    std::vector<HexView> cells ;
    for(size_t i = 0; i < blocks.size(); ++i) {
      blocks[i]->root->empty_resplit(std::vector<char>()) ;
      cells.push_back(HexView(blocks[i]->root, int(i), 0)) ;
    }
    write_hex_vtk(before, "before " + name, cells) ;

    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    add_hex_grid_wire_inputs(blocks, nodes, edges) ;
    write_wireframe_vtk(before_wire, "before wireframe " + name, nodes, edges) ;
    delete_hex_grid(blocks) ;
  }

  {
    std::vector<HexBlock*> blocks ;
    build_hex_grid(blocks, nx, ny, z_extent) ;

    std::vector<HexView> cells ;
    for(size_t i = 0; i < blocks.size(); ++i) {
      if(int(i) == center) {
        std::vector<HexCell*> leaves ;
        blocks[i]->root->resplit(center_plan,
                                 blocks[i]->nodes,
                                 blocks[i]->edges,
                                 blocks[i]->faces,
                                 leaves) ;
        for(size_t leaf = 0; leaf < leaves.size(); ++leaf) {
          cells.push_back(HexView(leaves[leaf], int(i), 6)) ;
        }
      } else {
        blocks[i]->root->empty_resplit(std::vector<char>()) ;
        cells.push_back(HexView(blocks[i]->root, int(i), 0)) ;
      }
    }
    write_hex_vtk(after, "after " + name, cells) ;

    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    add_hex_grid_wire_inputs(blocks, nodes, edges) ;
    write_wireframe_vtk(after_wire, "after wireframe " + name, nodes, edges) ;
    delete_hex_grid(blocks) ;
  }

  const std::vector<char> right_plan = extract_hex_face(center_plan, RIGHT) ;
  const std::vector<char> left_plan = extract_hex_face(center_plan, LEFT) ;
  const std::vector<char> front_plan = extract_hex_face(center_plan, FRONT) ;
  const std::vector<char> back_plan = extract_hex_face(center_plan, BACK) ;
  write_quad_face_vtk(right_face, "center right face " + name,
                      Point(3.0, 2.0, 0.0),
                      Point(3.0, 3.0, 0.0),
                      Point(3.0, 3.0, z_extent),
                      Point(3.0, 2.0, z_extent),
                      right_plan) ;
  write_quad_face_vtk(left_face, "center left face " + name,
                      Point(2.0, 2.0, 0.0),
                      Point(2.0, 3.0, 0.0),
                      Point(2.0, 3.0, z_extent),
                      Point(2.0, 2.0, z_extent),
                      left_plan) ;
  write_quad_face_vtk(front_face, "center front face " + name,
                      Point(2.0, 3.0, 0.0),
                      Point(3.0, 3.0, 0.0),
                      Point(3.0, 3.0, z_extent),
                      Point(2.0, 3.0, z_extent),
                      front_plan) ;
  write_quad_face_vtk(back_face, "center back face " + name,
                      Point(2.0, 2.0, 0.0),
                      Point(3.0, 2.0, 0.0),
                      Point(3.0, 2.0, z_extent),
                      Point(2.0, 2.0, z_extent),
                      back_plan) ;
  write_hex_grid_plan_note(note,
                           name,
                           center_plan,
                           right_plan,
                           left_plan,
                           front_plan,
                           back_plan,
                           "This is a direct library illustration. It shows "
                           "the refined center-cell tree and the face plans "
                           "induced on its four neighbors; it does not run "
                           "the full marker/refmesh scheduling path.") ;

  std::vector<std::string> files ;
  files.push_back(join_path(group, name + "_before.vtk")) ;
  files.push_back(join_path(group, name + "_after.vtk")) ;
  files.push_back(join_path(group, name + "_before_wire.vtk")) ;
  files.push_back(join_path(group, name + "_after_wire.vtk")) ;
  files.push_back(join_path(group, name + "_center_right_face.vtk")) ;
  files.push_back(join_path(group, name + "_center_left_face.vtk")) ;
  files.push_back(join_path(group, name + "_center_front_face.vtk")) ;
  files.push_back(join_path(group, name + "_center_back_face.vtk")) ;
  files.push_back(join_path(group, name + "_plans.dat")) ;
  add_manifest_case(manifest, group, name,
                    "A 5x5 thin hex grid with the center cell refined twice.",
                    files) ;
}

void write_hex_anisotropic_grid_case(const std::string& output_dir,
                                     std::ofstream& manifest,
                                     const std::string& name,
                                     const std::vector<char>& center_plan,
                                     int split_code,
                                     const std::string& description) {
  const std::string group = "anisotropic/hex" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;

  const int nx = 5 ;
  const int ny = 5 ;
  const int center = 2 + nx * 2 ;
  const double z_extent = 0.08 ;

  const std::string before = join_path(group_dir, name + "_before.vtk") ;
  const std::string after = join_path(group_dir, name + "_after.vtk") ;
  const std::string before_wire = join_path(group_dir, name + "_before_wire.vtk") ;
  const std::string after_wire = join_path(group_dir, name + "_after_wire.vtk") ;
  const std::string right_face = join_path(group_dir, name + "_center_right_face.vtk") ;
  const std::string left_face = join_path(group_dir, name + "_center_left_face.vtk") ;
  const std::string front_face = join_path(group_dir, name + "_center_front_face.vtk") ;
  const std::string back_face = join_path(group_dir, name + "_center_back_face.vtk") ;
  const std::string note = join_path(group_dir, name + "_plans.dat") ;

  {
    std::vector<HexBlock*> blocks ;
    build_hex_grid(blocks, nx, ny, z_extent) ;

    std::vector<HexView> cells ;
    for(size_t i = 0; i < blocks.size(); ++i) {
      blocks[i]->root->empty_resplit(std::vector<char>()) ;
      cells.push_back(HexView(blocks[i]->root, int(i), 0)) ;
    }
    write_hex_vtk(before, "before " + name, cells) ;

    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    add_hex_grid_wire_inputs(blocks, nodes, edges) ;
    write_wireframe_vtk(before_wire, "before wireframe " + name, nodes, edges) ;
    delete_hex_grid(blocks) ;
  }

  {
    std::vector<HexBlock*> blocks ;
    build_hex_grid(blocks, nx, ny, z_extent) ;

    std::vector<HexView> cells ;
    for(size_t i = 0; i < blocks.size(); ++i) {
      if(int(i) == center) {
        std::vector<HexCell*> leaves ;
        blocks[i]->root->resplit(center_plan,
                                 blocks[i]->nodes,
                                 blocks[i]->edges,
                                 blocks[i]->faces,
                                 leaves) ;
        for(size_t leaf = 0; leaf < leaves.size(); ++leaf) {
          cells.push_back(HexView(leaves[leaf], int(i), split_code)) ;
        }
      } else {
        blocks[i]->root->empty_resplit(std::vector<char>()) ;
        cells.push_back(HexView(blocks[i]->root, int(i), 0)) ;
      }
    }
    write_hex_vtk(after, "after " + name, cells) ;

    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    add_hex_grid_wire_inputs(blocks, nodes, edges) ;
    write_wireframe_vtk(after_wire, "after wireframe " + name, nodes, edges) ;
    delete_hex_grid(blocks) ;
  }

  const std::vector<char> right_plan = extract_hex_face(center_plan, RIGHT) ;
  const std::vector<char> left_plan = extract_hex_face(center_plan, LEFT) ;
  const std::vector<char> front_plan = extract_hex_face(center_plan, FRONT) ;
  const std::vector<char> back_plan = extract_hex_face(center_plan, BACK) ;
  write_quad_face_vtk(right_face, "center right face " + name,
                      Point(3.0, 2.0, 0.0),
                      Point(3.0, 3.0, 0.0),
                      Point(3.0, 3.0, z_extent),
                      Point(3.0, 2.0, z_extent),
                      right_plan) ;
  write_quad_face_vtk(left_face, "center left face " + name,
                      Point(2.0, 2.0, 0.0),
                      Point(2.0, 3.0, 0.0),
                      Point(2.0, 3.0, z_extent),
                      Point(2.0, 2.0, z_extent),
                      left_plan) ;
  write_quad_face_vtk(front_face, "center front face " + name,
                      Point(2.0, 3.0, 0.0),
                      Point(3.0, 3.0, 0.0),
                      Point(3.0, 3.0, z_extent),
                      Point(2.0, 3.0, z_extent),
                      front_plan) ;
  write_quad_face_vtk(back_face, "center back face " + name,
                      Point(2.0, 2.0, 0.0),
                      Point(3.0, 2.0, 0.0),
                      Point(3.0, 2.0, z_extent),
                      Point(2.0, 2.0, z_extent),
                      back_plan) ;
  write_hex_grid_plan_note(note,
                           name,
                           center_plan,
                           right_plan,
                           left_plan,
                           front_plan,
                           back_plan,
                           description + " This is a direct library "
                           "illustration of one anisotropically refined "
                           "center cell embedded in coarse neighbors.") ;

  std::vector<std::string> files ;
  files.push_back(join_path(group, name + "_before.vtk")) ;
  files.push_back(join_path(group, name + "_after.vtk")) ;
  files.push_back(join_path(group, name + "_before_wire.vtk")) ;
  files.push_back(join_path(group, name + "_after_wire.vtk")) ;
  files.push_back(join_path(group, name + "_center_right_face.vtk")) ;
  files.push_back(join_path(group, name + "_center_left_face.vtk")) ;
  files.push_back(join_path(group, name + "_center_front_face.vtk")) ;
  files.push_back(join_path(group, name + "_center_back_face.vtk")) ;
  files.push_back(join_path(group, name + "_plans.dat")) ;
  add_manifest_case(manifest, group, name, description, files) ;
}

void write_prism_pair_case(const std::string& output_dir,
                           std::ofstream& manifest) {
  const std::string group = "2d/prism" ;
  const std::string name = "square_refine_lower" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;

  const std::string before = join_path(group_dir, name + "_before_wire.vtk") ;
  const std::string after = join_path(group_dir, name + "_after_wire.vtk") ;
  const std::string face = join_path(group_dir, name + "_shared_face.vtk") ;
  const std::string note = join_path(group_dir, name + "_plans.dat") ;

  const Point lower[3] = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0)
  } ;
  const Point upper[3] = {
    Point(1.0, 1.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(1.0, 0.0, 0.0)
  } ;
  const double height = 0.08 ;
  const std::vector<char> lower_plan = plan({3}) ;

  {
    PrismBlock lower_prism ;
    PrismBlock upper_prism ;
    build_prism(lower_prism, lower, height) ;
    build_prism(upper_prism, upper, height) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&lower_prism.nodes) ;
    nodes.push_back(&upper_prism.nodes) ;
    edges.push_back(&lower_prism.edges) ;
    edges.push_back(&upper_prism.edges) ;
    write_wireframe_vtk(before, "before " + name, nodes, edges) ;
  }

  {
    PrismBlock lower_prism ;
    PrismBlock upper_prism ;
    build_prism(lower_prism, lower, height) ;
    build_prism(upper_prism, upper, height) ;
    std::vector<Prism*> leaves ;
    lower_prism.root->resplit(lower_plan,
                              lower_prism.nodes,
                              lower_prism.edges,
                              lower_prism.qfaces,
                              lower_prism.gfaces,
                              leaves) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&lower_prism.nodes) ;
    nodes.push_back(&upper_prism.nodes) ;
    edges.push_back(&lower_prism.edges) ;
    edges.push_back(&upper_prism.edges) ;
    write_wireframe_vtk(after, "after " + name, nodes, edges) ;
  }

  const std::vector<char> shared_face_plan = extract_prism_face(lower_plan, 3) ;
  write_quad_face_vtk(face, "shared prism face " + name,
                      Point(1.0, 0.0, 0.0),
                      Point(0.0, 1.0, 0.0),
                      Point(0.0, 1.0, height),
                      Point(1.0, 0.0, height),
                      shared_face_plan) ;
  write_plan_note(note, name, lower_plan, shared_face_plan,
                  "The triangular prism split induces this plan on the "
                  "diagonal quad face shared with the neighboring prism.") ;

  std::vector<std::string> files ;
  files.push_back(join_path(group, name + "_before_wire.vtk")) ;
  files.push_back(join_path(group, name + "_after_wire.vtk")) ;
  files.push_back(join_path(group, name + "_shared_face.vtk")) ;
  files.push_back(join_path(group, name + "_plans.dat")) ;
  add_manifest_case(manifest, group, name,
                    "Two extruded triangular prisms forming one square.",
                    files) ;
}

void write_prism_anisotropic_pair_case(const std::string& output_dir,
                                       std::ofstream& manifest,
                                       const std::string& name,
                                       const std::vector<char>& lower_plan,
                                       const std::string& description) {
  const std::string group = "anisotropic/prism" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;

  const std::string before = join_path(group_dir, name + "_before_wire.vtk") ;
  const std::string after = join_path(group_dir, name + "_after_wire.vtk") ;
  const std::string face = join_path(group_dir, name + "_shared_face.vtk") ;
  const std::string note = join_path(group_dir, name + "_plans.dat") ;

  const Point lower[3] = {
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0)
  } ;
  const Point upper[3] = {
    Point(1.0, 1.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(1.0, 0.0, 0.0)
  } ;
  const double height = 0.08 ;

  {
    PrismBlock lower_prism ;
    PrismBlock upper_prism ;
    build_prism(lower_prism, lower, height) ;
    build_prism(upper_prism, upper, height) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&lower_prism.nodes) ;
    nodes.push_back(&upper_prism.nodes) ;
    edges.push_back(&lower_prism.edges) ;
    edges.push_back(&upper_prism.edges) ;
    write_wireframe_vtk(before, "before " + name, nodes, edges) ;
  }

  {
    PrismBlock lower_prism ;
    PrismBlock upper_prism ;
    build_prism(lower_prism, lower, height) ;
    build_prism(upper_prism, upper, height) ;
    std::vector<Prism*> leaves ;
    lower_prism.root->resplit(lower_plan,
                              lower_prism.nodes,
                              lower_prism.edges,
                              lower_prism.qfaces,
                              lower_prism.gfaces,
                              leaves) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&lower_prism.nodes) ;
    nodes.push_back(&upper_prism.nodes) ;
    edges.push_back(&lower_prism.edges) ;
    edges.push_back(&upper_prism.edges) ;
    write_wireframe_vtk(after, "after " + name, nodes, edges) ;
  }

  const std::vector<char> shared_face_plan = extract_prism_face(lower_plan, 3) ;
  write_quad_face_vtk(face, "shared prism face " + name,
                      Point(1.0, 0.0, 0.0),
                      Point(0.0, 1.0, 0.0),
                      Point(0.0, 1.0, height),
                      Point(1.0, 0.0, height),
                      shared_face_plan) ;
  write_plan_note(note, name, lower_plan, shared_face_plan,
                  description + " The upper prism remains coarse so the "
                  "shared-face artifact shows the induced interface plan.") ;

  std::vector<std::string> files ;
  files.push_back(join_path(group, name + "_before_wire.vtk")) ;
  files.push_back(join_path(group, name + "_after_wire.vtk")) ;
  files.push_back(join_path(group, name + "_shared_face.vtk")) ;
  files.push_back(join_path(group, name + "_plans.dat")) ;
  add_manifest_case(manifest, group, name, description, files) ;
}

void write_single_prism_case(const std::string& output_dir,
                             std::ofstream& manifest) {
  const std::string group = "3d/prism" ;
  const std::string name = "full_split" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;

  const std::string before = join_path(group_dir, name + "_before_wire.vtk") ;
  const std::string after = join_path(group_dir, name + "_after_wire.vtk") ;
  const std::vector<char> cell_plan = plan({3}) ;

  {
    PrismFixture fixture ;
    build_unit_prism(fixture) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&fixture.nodes) ;
    edges.push_back(&fixture.edges) ;
    write_wireframe_vtk(before, "before " + name, nodes, edges) ;
  }

  {
    PrismFixture fixture ;
    build_unit_prism(fixture) ;
    std::vector<Prism*> leaves ;
    fixture.root->resplit(cell_plan,
                          fixture.nodes,
                          fixture.edges,
                          fixture.qfaces,
                          fixture.gfaces,
                          leaves) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&fixture.nodes) ;
    edges.push_back(&fixture.edges) ;
    write_wireframe_vtk(after, "after " + name, nodes, edges) ;
  }

  std::vector<std::string> files ;
  files.push_back(join_path(group, name + "_before_wire.vtk")) ;
  files.push_back(join_path(group, name + "_after_wire.vtk")) ;
  add_manifest_case(manifest, group, name,
                    "A single 3D prism before and after a full prism split.",
                    files) ;
}

void write_general_cell_case(const std::string& output_dir,
                             std::ofstream& manifest,
                             const std::string& group,
                             const std::string& name,
                             GeneralBuilder builder,
                             const std::string& description) {
  const std::string group_dir = output_group_dir(output_dir, group) ;
  const std::string before = join_path(group_dir, name + "_before_wire.vtk") ;
  const std::string after = join_path(group_dir, name + "_after_wire.vtk") ;
  const std::vector<char> cell_plan = plan({1}) ;

  {
    GeneralFixture fixture ;
    builder(fixture) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&fixture.nodes) ;
    edges.push_back(&fixture.edges) ;
    write_wireframe_vtk(before, "before " + name, nodes, edges) ;
  }

  {
    GeneralFixture fixture ;
    builder(fixture) ;
    std::vector<DiamondCell*> leaves ;
    fixture.root->resplit(cell_plan,
                          fixture.nodes,
                          fixture.edges,
                          fixture.faces,
                          leaves) ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&fixture.nodes) ;
    edges.push_back(&fixture.edges) ;
    write_wireframe_vtk(after, "after " + name, nodes, edges) ;
  }

  std::vector<std::string> files ;
  files.push_back(join_path(group, name + "_before_wire.vtk")) ;
  files.push_back(join_path(group, name + "_after_wire.vtk")) ;
  add_manifest_case(manifest, group, name, description, files) ;
}

void write_hex_split_catalog(const std::string& output_dir,
                             std::ofstream& manifest) {
  const std::string group = "split_catalog/hex" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;
  std::vector<std::string> files ;

  for(int code = 0; code <= 7; ++code) {
    HexFixture fixture ;
    build_unit_hex(fixture) ;
    std::vector<HexCell*> leaves ;
    fixture.root->resplit(plan({code}),
                          fixture.nodes,
                          fixture.edges,
                          fixture.faces,
                          leaves) ;
    std::vector<HexView> cells ;
    for(size_t i = 0; i < leaves.size(); ++i) {
      cells.push_back(HexView(leaves[i], 0, code)) ;
    }

    std::ostringstream cell_name ;
    cell_name << "split_" << code << "_cells.vtk" ;
    write_hex_vtk(join_path(group_dir, cell_name.str()),
                  "HexCell root split " + plan_string(plan({code})),
                  cells) ;
    files.push_back(join_path(group, cell_name.str())) ;

    std::ostringstream wire_name ;
    wire_name << "split_" << code << "_wire.vtk" ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&fixture.nodes) ;
    edges.push_back(&fixture.edges) ;
    write_wireframe_vtk(join_path(group_dir, wire_name.str()),
                        "HexCell root split wireframe " + plan_string(plan({code})),
                        nodes,
                        edges) ;
    files.push_back(join_path(group, wire_name.str())) ;
  }

  add_manifest_case(manifest,
                    group,
                    "root_split_strategies",
                    "Root split visuals for HexCell split codes 0 through 7.",
                    files) ;
}

void write_prism_split_catalog(const std::string& output_dir,
                               std::ofstream& manifest) {
  const std::string group = "split_catalog/prism" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;
  std::vector<std::string> files ;

  for(int code = 0; code <= 3; ++code) {
    PrismFixture fixture ;
    build_unit_prism(fixture) ;
    std::vector<Prism*> leaves ;
    fixture.root->resplit(plan({code}),
                          fixture.nodes,
                          fixture.edges,
                          fixture.qfaces,
                          fixture.gfaces,
                          leaves) ;

    std::ostringstream name ;
    name << "split_" << code << "_wire.vtk" ;
    std::vector<std::list<Node*>*> nodes ;
    std::vector<std::list<Edge*>*> edges ;
    nodes.push_back(&fixture.nodes) ;
    edges.push_back(&fixture.edges) ;
    write_wireframe_vtk(join_path(group_dir, name.str()),
                        "Prism root split " + plan_string(plan({code})),
                        nodes,
                        edges) ;
    files.push_back(join_path(group, name.str())) ;
  }

  add_manifest_case(manifest,
                    group,
                    "root_split_strategies",
                    "Root split visuals for Prism split codes 0 through 3.",
                    files) ;
}

void write_general_split_catalog_case(const std::string& output_dir,
                                      const std::string& shape,
                                      GeneralBuilder builder,
                                      const std::vector<char>& cell_plan,
                                      std::vector<std::string>& files) {
  const std::string group = "split_catalog/general" ;
  const std::string group_dir = output_group_dir(output_dir, group) ;

  GeneralFixture fixture ;
  builder(fixture) ;
  if(!cell_plan.empty()) {
    std::vector<DiamondCell*> leaves ;
    fixture.root->resplit(cell_plan,
                          fixture.nodes,
                          fixture.edges,
                          fixture.faces,
                          leaves) ;
  }

  std::ostringstream name ;
  name << shape << "_split_" << (cell_plan.empty() ? 0 : int(cell_plan[0]))
       << "_wire.vtk" ;
  std::vector<std::list<Node*>*> nodes ;
  std::vector<std::list<Edge*>*> edges ;
  nodes.push_back(&fixture.nodes) ;
  edges.push_back(&fixture.edges) ;
  write_wireframe_vtk(join_path(group_dir, name.str()),
                      shape + " split " + plan_string(cell_plan),
                      nodes,
                      edges) ;
  files.push_back(join_path(group, name.str())) ;
}

void write_split_catalog(const std::string& output_dir,
                         std::ofstream& manifest) {
  write_hex_split_catalog(output_dir, manifest) ;
  write_prism_split_catalog(output_dir, manifest) ;

  std::vector<std::string> files ;
  write_general_split_catalog_case(output_dir,
                                   "tetra",
                                   build_tetra_cell,
                                   std::vector<char>(),
                                   files) ;
  write_general_split_catalog_case(output_dir,
                                   "tetra",
                                   build_tetra_cell,
                                   plan({1}),
                                   files) ;
  write_general_split_catalog_case(output_dir,
                                   "pyramid",
                                   build_pyramid_cell,
                                   std::vector<char>(),
                                   files) ;
  write_general_split_catalog_case(output_dir,
                                   "pyramid",
                                   build_pyramid_cell,
                                   plan({1}),
                                   files) ;

  add_manifest_case(manifest,
                    "split_catalog/general",
                    "root_split_strategies",
                    "Root split visuals for tetrahedral and pyramidal "
                    "general Cell shapes.",
                    files) ;
}

void write_anisotropic_cases(const std::string& output_dir,
                             std::ofstream& manifest) {
  write_hex_anisotropic_grid_case(output_dir,
                                  manifest,
                                  "center_x_split",
                                  plan({4}),
                                  4,
                                  "A 5x5 thin hex grid with the center cell "
                                  "split in the x direction.") ;
  write_hex_anisotropic_grid_case(output_dir,
                                  manifest,
                                  "center_y_split",
                                  plan({2}),
                                  2,
                                  "A 5x5 thin hex grid with the center cell "
                                  "split in the y direction.") ;
  write_hex_anisotropic_grid_case(output_dir,
                                  manifest,
                                  "center_z_split",
                                  plan({1}),
                                  1,
                                  "A 5x5 thin hex grid with the center cell "
                                  "split through the thin z direction.") ;
  write_hex_anisotropic_grid_case(output_dir,
                                  manifest,
                                  "center_xy_split",
                                  plan({6}),
                                  6,
                                  "A 5x5 thin hex grid with the center cell "
                                  "split in both x and y.") ;

  write_prism_anisotropic_pair_case(output_dir,
                                    manifest,
                                    "pair_code1_two_child_split",
                                    plan({1}),
                                    "Two neighboring triangular prisms with "
                                    "the lower prism using split code 1, "
                                    "which creates two children.") ;
  write_prism_anisotropic_pair_case(output_dir,
                                    manifest,
                                    "pair_code2_three_child_split",
                                    plan({2}),
                                    "Two neighboring triangular prisms with "
                                    "the lower prism using split code 2, "
                                    "which creates three children.") ;
  write_prism_anisotropic_pair_case(output_dir,
                                    manifest,
                                    "pair_code3_six_child_split",
                                    plan({3}),
                                    "Two neighboring triangular prisms with "
                                    "the lower prism using split code 3, "
                                    "which creates six children.") ;
}

void write_all(const std::string& output_dir) {
  make_directories(output_dir) ;

  const std::string manifest_path = join_path(output_dir, "manifest.dat") ;
  std::ofstream manifest(manifest_path.c_str()) ;
  if(!manifest) {
    throw std::runtime_error("could not open illustration manifest") ;
  }
  manifest << "# FVMAdapt illustration artifacts\n" ;
  manifest << "# Run from quickTest/FVMAdaptTests/Illustrations with make.\n\n" ;

  write_hex_pair_case(output_dir, manifest,
                      "2d/hex",
                      "strip_refine_left_xy",
                      plan({6}),
                      6,
                      0.08) ;
  write_hex_grid_focus_case(output_dir, manifest) ;
  write_prism_pair_case(output_dir, manifest) ;
  write_hex_pair_case(output_dir, manifest,
                      "3d/hex",
                      "pair_refine_left_xyz",
                      plan({7}),
                      7,
                      1.0) ;
  write_single_prism_case(output_dir, manifest) ;
  write_general_cell_case(output_dir, manifest,
                          "3d/general/tetra",
                          "split",
                          build_tetra_cell,
                          "A tetrahedral general Cell before and after "
                          "isotropic split.") ;
  write_general_cell_case(output_dir, manifest,
                          "3d/general/pyramid",
                          "split",
                          build_pyramid_cell,
                          "A pyramidal general Cell before and after "
                          "isotropic split.") ;
  write_anisotropic_cases(output_dir, manifest) ;
  write_split_catalog(output_dir, manifest) ;

  std::cout << "wrote " << manifest_path << "\n" ;
  std::cout << "wrote illustration VTK files under "
            << join_path(output_dir, "2d") << ", "
            << join_path(output_dir, "3d") << ", "
            << join_path(output_dir, "anisotropic") << ", and "
            << join_path(output_dir, "split_catalog") << "\n" ;
}

} // namespace

int main(int argc, char** argv) {
  std::string output_dir = "output" ;

  for(int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]) ;
    if(arg == "--write-dir") {
      if(i + 1 >= argc) {
        std::cerr << "--write-dir requires a path\n" ;
        return 1 ;
      }
      output_dir = argv[++i] ;
    } else {
      std::cerr << "unknown argument: " << arg << "\n" ;
      return 1 ;
    }
  }

  try {
    write_all(output_dir) ;
    return 0 ;
  } catch(const std::exception& error) {
    std::cerr << error.what() << "\n" ;
    return 1 ;
  }
}
