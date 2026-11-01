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

#include <FVMAdapt/parametric/native_cell.h>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {

    typedef std::vector<EntityId> Loop ;
    typedef std::pair<EntityId,EntityId> Edge ;

    Loop sorted(Loop vertices) {
      std::sort(vertices.begin(),vertices.end()) ;
      return vertices ;
    }

    Shape cellShape(const std::vector<Loop> &faces) {
      unsigned triangles = 0, quads = 0 ;
      for(const Loop &face : faces) {
        triangles += face.size() == 3 ;
        quads += face.size() == 4 ;
      }
      if(faces.size() == 6 && quads == 6) return Shape::hex ;
      if(faces.size() == 5 && triangles == 2 && quads == 3) return Shape::prism ;
      if(faces.size() == 4 && triangles == 4) return Shape::tet ;
      if(faces.size() == 5 && triangles == 4 && quads == 1) return Shape::pyramid ;
      throw std::invalid_argument("native cell requires tet, pyramid, prism or hex faces") ;
    }

    bool matchFaces(NativeCell &cell, const std::vector<Loop> &faces,
                    const std::map<Loop,unsigned> &faceIndex) {
      const ReferenceElement &element = referenceElement(cell.shape) ;
      cell.faces.clear() ;
      cell.faceOrientations.clear() ;
      for(unsigned f = 0; f < element.boundaryCount; ++f) {
        const Boundary &boundary = element.boundaries[f] ;
        Loop loop ;
        for(unsigned v = 0; v < boundary.vertexCount; ++v)
          loop.push_back(cell.vertices[boundary.vertices[v]]) ;
        const auto found = faceIndex.find(sorted(loop)) ;
        if(found == faceIndex.end()) return false ;
        const Loop &input = faces[found->second] ;
        unsigned orientation = 0 ;
        while(orientation < input.size() && input[orientation] != loop[0]) ++orientation ;
        for(unsigned v = 0; v < loop.size(); ++v)
          if(input[(v+orientation)%input.size()] != loop[v]) return false ;
        cell.faces.push_back(found->second) ;
        cell.faceOrientations.push_back(orientation) ;
      }
      return true ;
    }

  }

  NativeCell recognizeNativeCell(const std::vector<std::vector<EntityId> > &outwardFaces) {
    const Shape shape = cellShape(outwardFaces) ;
    const ReferenceElement &element = referenceElement(shape) ;
    std::map<Loop,unsigned> faceIndex ;
    std::map<EntityId,std::set<EntityId> > neighbors ;
    std::map<Edge,std::pair<unsigned,int> > edgeUses ;
    for(unsigned f = 0; f < outwardFaces.size(); ++f) {
      const Loop &face = outwardFaces[f] ;
      const Loop key = sorted(face) ;
      if(std::adjacent_find(key.begin(),key.end()) != key.end())
        throw std::invalid_argument("native cell face repeats a vertex") ;
      if(!faceIndex.emplace(key,f).second)
        throw std::invalid_argument("native cell repeats a face") ;
      for(unsigned v = 0; v < face.size(); ++v) {
        const EntityId a = face[v], b = face[(v+1)%face.size()] ;
        neighbors[a].insert(b) ;
        neighbors[b].insert(a) ;
        auto &use = edgeUses[std::minmax(a,b)] ;
        ++use.first ;
        use.second += a < b ? 1 : -1 ;
      }
    }
    if(neighbors.size() != element.vertexCount || edgeUses.size() != element.edgeCount)
      throw std::invalid_argument("native cell vertex or edge count is invalid") ;
    unsigned degreeFour = 0 ;
    for(const auto &entry : neighbors) {
      degreeFour += entry.second.size() == 4 ;
      if(entry.second.size() != 3 && !(shape == Shape::pyramid && entry.second.size() == 4))
        throw std::invalid_argument("native cell vertex has invalid incident edges") ;
    }
    if(degreeFour != (shape == Shape::pyramid ? 1u : 0u))
      throw std::invalid_argument("native pyramid requires one four-edge apex") ;
    for(const auto &entry : edgeUses)
      if(entry.second.first != 2 || entry.second.second != 0)
        throw std::invalid_argument("native cell edges require two opposite face uses") ;

    NativeCell best{shape,{},{},{}} ;
    const Boundary &base = element.boundaries[0] ;
    for(const Loop &face : outwardFaces) {
      if(face.size() != base.vertexCount) continue ;
      const std::set<EntityId> baseVertices(face.begin(),face.end()) ;
      for(unsigned rotation = 0; rotation < face.size(); ++rotation) {
        NativeCell candidate{shape,Loop(element.vertexCount),{},{}} ;
        for(unsigned v = 0; v < base.vertexCount; ++v)
          candidate.vertices[base.vertices[v]] = face[(v+rotation)%face.size()] ;
        bool valid = true ;
        if(shape == Shape::tet || shape == Shape::pyramid) {
          for(const auto &vertex : neighbors)
            if(!baseVertices.count(vertex.first)) candidate.vertices[base.vertexCount] = vertex.first ;
        } else {
          // Hexes and prisms connect each base vertex to one opposite-base vertex.
          for(unsigned v = 0; v < base.vertexCount; ++v) {
            unsigned count = 0 ;
            for(EntityId neighbor : neighbors.at(candidate.vertices[v]))
              if(!baseVertices.count(neighbor)) {
                candidate.vertices[v+base.vertexCount] = neighbor ;
                ++count ;
              }
            valid = valid && count == 1 ;
          }
        }
        const Loop vertexIds = sorted(candidate.vertices) ;
        if(!valid || std::adjacent_find(vertexIds.begin(),vertexIds.end()) != vertexIds.end()) continue ;
        if(matchFaces(candidate,outwardFaces,faceIndex) &&
           (best.vertices.empty() || candidate.vertices < best.vertices))
          best = std::move(candidate) ;
      }
    }
    if(best.vertices.empty())
      throw std::invalid_argument("native cell faces do not match their reference topology") ;
    return best ;
  }

}
}
