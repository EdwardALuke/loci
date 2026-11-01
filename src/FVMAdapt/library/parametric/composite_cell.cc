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

#include <FVMAdapt/parametric/composite_cell.h>

#include "exact_geometry.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {

    using detail::Exact ;
    using detail::ExactPoint ;
    using Loop = std::vector<unsigned> ;
    using Triangle = std::array<unsigned,3> ;
    using Edge = std::pair<unsigned,unsigned> ;

    ExactPoint exactPoint(const Point &point) {
      for(double value : point)
        if(!std::isfinite(value))
          throw std::invalid_argument("composite geometry coordinates must be finite") ;
      return {{Exact(point[0]),Exact(point[1]),Exact(point[2])}} ;
    }

    ExactPoint subtract(const ExactPoint &a, const ExactPoint &b) {
      return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}} ;
    }

    ExactPoint cross(const ExactPoint &a, const ExactPoint &b) {
      return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
    }

    Exact dot(const ExactPoint &a, const ExactPoint &b) {
      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] ;
    }

    bool nonzero(const ExactPoint &point) {
      return point[0] != 0 || point[1] != 0 || point[2] != 0 ;
    }

    Exact determinant(const ExactPoint &a, const ExactPoint &b, const ExactPoint &c) {
      return dot(a,cross(b,c)) ;
    }

    Loop sorted(Loop loop) {
      std::sort(loop.begin(),loop.end()) ;
      return loop ;
    }

    Loop canonicalLoop(const std::vector<EntityId> &ids) {
      const unsigned count = ids.size() ;
      const unsigned first = std::min_element(ids.begin(),ids.end())-ids.begin() ;
      const bool forward = ids[(first+1)%count] < ids[(first+count-1)%count] ;
      Loop loop ;
      for(unsigned i = 0; i < count; ++i)
        loop.push_back((first+(forward ? i : count-i))%count) ;
      return loop ;
    }

    void validateShell(const std::vector<Point> &points,
                       const std::vector<EntityId> &vertexIds,
                       const std::vector<CompositeFace> &faces) {
      if(points.size() < 4 || points.size() != vertexIds.size() ||
         points.size() >= noCompositeIndex || faces.size() >= noCompositeIndex || faces.size() < 4)
        throw std::invalid_argument("composite cell requires matching point IDs and a closed shell") ;
      if(std::set<EntityId>(vertexIds.begin(),vertexIds.end()).size() != vertexIds.size())
        throw std::invalid_argument("composite cell repeats a vertex ID") ;
      std::set<Point> uniquePoints ;
      for(const Point &point : points) {
        exactPoint(point) ;
        if(!uniquePoints.insert(point).second)
          throw std::invalid_argument("composite cell repeats a vertex position") ;
      }
      std::set<EntityId> faceIds ;
      std::set<Loop> faceLoops ;
      std::map<Edge,std::vector<std::pair<unsigned,bool> > > edges ;
      std::vector<std::map<unsigned,std::vector<unsigned> > > links(points.size()) ;
      for(unsigned f = 0; f < faces.size(); ++f) {
        const CompositeFace &face = faces[f] ;
        if(!faceIds.insert(face.id).second)
          throw std::invalid_argument("composite cell repeats a face ID") ;
        if(face.vertices.size() < 3 || face.vertices.size() >= noCompositeIndex)
          throw std::invalid_argument("composite face requires at least three vertices") ;
        const Loop key = sorted(face.vertices) ;
        if(key.back() >= points.size() || std::adjacent_find(key.begin(),key.end()) != key.end())
          throw std::invalid_argument("composite face has invalid or repeated vertex indices") ;
        if(!faceLoops.insert(key).second)
          throw std::invalid_argument("composite cell repeats a face") ;
        for(unsigned v = 0; v < face.vertices.size(); ++v) {
          const unsigned a = face.vertices[v], b = face.vertices[(v+1)%face.vertices.size()] ;
          const unsigned previous = face.vertices[(v+face.vertices.size()-1)%face.vertices.size()] ;
          edges[std::minmax(a,b)].push_back({f,a < b}) ;
          links[a][previous].push_back(b) ;
          links[a][b].push_back(previous) ;
        }
      }
      std::vector<Loop> neighbors(faces.size()) ;
      for(const auto &edge : edges) {
        const auto &uses = edge.second ;
        if(uses.size() != 2 || uses[0].second == uses[1].second)
          throw std::invalid_argument("composite shell edges require two opposite face uses") ;
        neighbors[uses[0].first].push_back(uses[1].first) ;
        neighbors[uses[1].first].push_back(uses[0].first) ;
      }
      Loop pending{0} ;
      std::set<unsigned> visited ;
      while(!pending.empty()) {
        const unsigned face = pending.back() ;
        pending.pop_back() ;
        if(!visited.insert(face).second) continue ;
        pending.insert(pending.end(),neighbors[face].begin(),neighbors[face].end()) ;
      }
      if(visited.size() != faces.size() || points.size()+faces.size() != edges.size()+2)
        throw std::invalid_argument("composite shell must be one connected sphere") ;
      for(const auto &link : links) {
        if(link.empty()) throw std::invalid_argument("composite cell has an unused vertex") ;
        for(const auto &vertex : link)
          if(vertex.second.size() != 2)
            throw std::invalid_argument("composite vertex link must be one cycle") ;
        visited.clear() ;
        pending = {link.begin()->first} ;
        while(!pending.empty()) {
          const unsigned vertex = pending.back() ;
          pending.pop_back() ;
          if(!visited.insert(vertex).second) continue ;
          const auto &next = link.at(vertex) ;
          pending.insert(pending.end(),next.begin(),next.end()) ;
        }
        if(visited.size() != link.size())
          throw std::invalid_argument("composite vertex link must be one cycle") ;
      }
    }

    struct RadialTriangle {
      std::array<ExactPoint,3> vertices ;
      std::array<ExactPoint,3> normals ;
    } ;

    RadialTriangle radialTriangle(const Triangle &triangle,
                                  const std::vector<ExactPoint> &radial) {
      const ExactPoint &a = radial[triangle[0]], &b = radial[triangle[1]], &c = radial[triangle[2]] ;
      if(determinant(a,b,c) <= 0)
        throw std::invalid_argument("composite interior point must strictly see every face patch") ;
      return {{{a,b,c}},{{cross(a,b),cross(b,c),cross(c,a)}}} ;
    }

    bool separated(const RadialTriangle &first, const RadialTriangle &second) {
      for(const ExactPoint &normal : first.normals) {
        bool outside = true ;
        for(const ExactPoint &vertex : second.vertices)
          if(dot(normal,vertex) > 0) { outside = false ; break ; }
        if(outside) return true ;
      }
      return false ;
    }

    bool overlap(const RadialTriangle &first, const RadialTriangle &second) {
      if(separated(first,second) || separated(second,first)) return false ;
      std::array<ExactPoint,6> normals ;
      std::copy(first.normals.begin(),first.normals.end(),normals.begin()) ;
      std::copy(second.normals.begin(),second.normals.end(),normals.begin()+3) ;
      ExactPoint firstRay{}, plane{} ;
      unsigned rank = 0 ;
      // Intersect the two pointed cones. Three independent feasible extreme rays
      // certify a volume overlap; shared edge and vertex sectors have lower rank.
      for(unsigned i = 0; i < normals.size(); ++i) {
        for(unsigned j = i+1; j < normals.size(); ++j) {
          ExactPoint ray = cross(normals[i],normals[j]) ;
          if(!nonzero(ray)) continue ;
          for(unsigned direction = 0; direction < 2; ++direction) {
            bool feasible = true ;
            for(const ExactPoint &normal : normals)
              if(dot(normal,ray) < 0) { feasible = false ; break ; }
            if(feasible) {
              if(rank == 0) { firstRay = ray ; rank = 1 ; }
              else if(rank == 1) {
                plane = cross(firstRay,ray) ;
                if(nonzero(plane)) rank = 2 ;
              } else if(dot(plane,ray) != 0) return true ;
            }
            for(Exact &coordinate : ray) coordinate = -coordinate ;
          }
        }
      }
      return false ;
    }

    bool opposite(const Loop &first, const Loop &second) {
      if(first.size() != second.size()) return false ;
      const auto start = std::find(second.begin(),second.end(),first[0]) ;
      if(start == second.end()) return false ;
      const unsigned offset = start-second.begin(), count = first.size() ;
      for(unsigned i = 0; i < count; ++i)
        if(first[i] != second[(offset+count-i)%count]) return false ;
      return true ;
    }

  }

  std::vector<std::array<unsigned,3> > triangulateFace(
      const std::vector<EntityId> &vertexIds, const std::vector<Point> &points) {
    if(points.size() < 3 || points.size() != vertexIds.size() || points.size() >= noCompositeIndex ||
       std::set<EntityId>(vertexIds.begin(),vertexIds.end()).size() != points.size())
      throw std::invalid_argument("face triangulation requires distinct IDs for at least three points") ;
    std::vector<ExactPoint> exactPoints ;
    std::set<Point> unique ;
    for(const Point &point : points) {
      exactPoints.push_back(exactPoint(point)) ;
      if(!unique.insert(point).second)
        throw std::invalid_argument("face triangulation repeats a point") ;
    }
    Loop loop = canonicalLoop(vertexIds) ;
    const bool reversed = loop[1] != (loop[0]+1)%loop.size() ;
    ExactPoint normal{} ;
    const auto &origin = exactPoints[loop[0]] ;
    for(unsigned i = 1; i+1 < loop.size(); ++i) {
      const auto area = cross(subtract(exactPoints[loop[i]],origin),
                              subtract(exactPoints[loop[i+1]],origin)) ;
      for(unsigned d = 0; d < 3; ++d) normal[d] += area[d] ;
    }
    if(!nonzero(normal)) throw std::invalid_argument("face triangulation requires nonzero area") ;
    // Ear selection uses the area-normal projection; returned triangles retain the 3D vertices.
    for(unsigned i = 0; i < loop.size(); ++i) {
      const ExactPoint &a = exactPoints[loop[i]], &b = exactPoints[loop[(i+1)%loop.size()]] ;
      for(const ExactPoint &point : exactPoints)
        if(dot(normal,cross(subtract(b,a),subtract(point,a))) < 0)
          throw std::invalid_argument("face triangulation requires a convex polygon loop") ;
    }
    std::vector<Triangle> triangles ;
    while(loop.size() > 3) {
      unsigned ear = noCompositeIndex ;
      for(unsigned i = 0; i < loop.size(); ++i) {
        const unsigned previous = loop[(i+loop.size()-1)%loop.size()] ;
        const unsigned middle = loop[i], next = loop[(i+1)%loop.size()] ;
        const ExactPoint &a = exactPoints[previous], &b = exactPoints[middle], &c = exactPoints[next] ;
        if(dot(normal,cross(subtract(b,a),subtract(c,a))) <= 0) continue ;
        bool empty = true ;
        for(unsigned other : loop) {
          if(other == previous || other == middle || other == next) continue ;
          const ExactPoint &point = exactPoints[other] ;
          if(dot(normal,cross(subtract(b,a),subtract(point,a))) >= 0 &&
             dot(normal,cross(subtract(c,b),subtract(point,b))) >= 0 &&
             dot(normal,cross(subtract(a,c),subtract(point,c))) >= 0) {
            empty = false ; break ;
          }
        }
        if(empty && (ear == noCompositeIndex || vertexIds[middle] < vertexIds[loop[ear]])) ear = i ;
      }
      if(ear == noCompositeIndex)
        throw std::invalid_argument("face triangulation cannot retain all boundary vertices") ;
      triangles.push_back({{loop[(ear+loop.size()-1)%loop.size()],loop[ear],loop[(ear+1)%loop.size()]}}) ;
      loop.erase(loop.begin()+ear) ;
    }
    triangles.push_back({{loop[0],loop[1],loop[2]}}) ;
    for(Triangle &triangle : triangles)
      if(reversed) std::swap(triangle[1],triangle[2]) ;
    std::sort(triangles.begin(),triangles.end(),[&](const Triangle &first, const Triangle &second) {
      std::array<EntityId,3> a{{vertexIds[first[0]],vertexIds[first[1]],vertexIds[first[2]]}} ;
      std::array<EntityId,3> b{{vertexIds[second[0]],vertexIds[second[1]],vertexIds[second[2]]}} ;
      std::sort(a.begin(),a.end()) ; std::sort(b.begin(),b.end()) ;
      return a < b ;
    }) ;
    return triangles ;
  }

  CompositeCell makeCompositeCell(const std::vector<Point> &points,
                                 const std::vector<EntityId> &vertexIds,
                                 const std::vector<CompositeFace> &faces,
                                 std::optional<Point> interiorPoint) {
    validateShell(points,vertexIds,faces) ;
    CompositeCell result ;
    result.vertices = points ;
    result.facePatches.resize(faces.size()) ;
    for(unsigned f = 0; f < faces.size(); ++f) {
      const Loop &loop = faces[f].vertices ;
      if(loop.size() <= 4) result.facePatches[f].push_back(loop) ;
      else {
        std::vector<Point> positions ;
        std::vector<EntityId> ids ;
        for(unsigned vertex : loop) { positions.push_back(points[vertex]) ; ids.push_back(vertexIds[vertex]) ; }
        for(const Triangle &triangle : triangulateFace(ids,positions))
          result.facePatches[f].push_back({loop[triangle[0]],loop[triangle[1]],loop[triangle[2]]}) ;
      }
    }
    std::size_t pieceCount = 0 ;
    for(const auto &patches : result.facePatches) {
      if(patches.size() >= noCompositeIndex-pieceCount)
        throw std::length_error("composite cell has too many geometry pieces") ;
      pieceCount += patches.size() ;
    }
    result.pieces.reserve(pieceCount) ;
    Point center{} ;
    if(interiorPoint) center = *interiorPoint ;
    else {
      ExactPoint sum{} ;
      for(const Point &point : points)
        for(unsigned d = 0; d < 3; ++d) sum[d] += Exact(point[d]) ;
      for(unsigned d = 0; d < 3; ++d) center[d] = (sum[d]/points.size()).convert_to<double>() ;
    }
    const ExactPoint origin = exactPoint(center) ;
    std::vector<ExactPoint> radial ;
    for(const Point &point : points) radial.push_back(subtract(exactPoint(point),origin)) ;
    std::vector<RadialTriangle> sectors ;
    for(unsigned f = 0; f < faces.size(); ++f) {
      for(const Loop &patch : result.facePatches[f]) {
        if(patch.size() == 3) sectors.push_back(radialTriangle({{patch[0],patch[1],patch[2]}},radial)) ;
        else {
          // The pyramid determinant is bilinear on its collapsed square. Positive
          // corner values certify its sign throughout, without flattening the quad.
          for(unsigned v = 0; v < 4; ++v)
            if(determinant(radial[patch[v]],radial[patch[(v+1)%4]],radial[patch[(v+3)%4]]) <= 0)
              throw std::invalid_argument("composite quad cone requires positive corner Jacobians") ;
          unsigned anchor = 0 ;
          for(unsigned v = 1; v < 4; ++v)
            if(vertexIds[patch[v]] < vertexIds[patch[anchor]]) anchor = v ;
          sectors.push_back(radialTriangle({{patch[anchor],patch[(anchor+1)%4],patch[(anchor+2)%4]}},radial)) ;
          sectors.push_back(radialTriangle({{patch[anchor],patch[(anchor+2)%4],patch[(anchor+3)%4]}},radial)) ;
        }
      }
    }
    for(unsigned i = 0; i < sectors.size(); ++i)
      for(unsigned j = i+1; j < sectors.size(); ++j)
        if(overlap(sectors[i],sectors[j]))
          throw std::invalid_argument("composite face sectors overlap from the interior point") ;
    const unsigned centerIndex = result.vertices.size() ;
    result.vertices.push_back(center) ;
    for(unsigned f = 0; f < faces.size(); ++f) {
      for(unsigned p = 0; p < result.facePatches[f].size(); ++p) {
        const Loop &patch = result.facePatches[f][p] ;
        GeometryPiece piece ;
        piece.shape = patch.size() == 3 ? Shape::tet : Shape::pyramid ;
        piece.vertices.push_back(patch[0]) ;
        for(unsigned v = patch.size()-1; v > 0; --v) piece.vertices.push_back(patch[v]) ;
        piece.vertices.push_back(centerIndex) ;
        std::vector<Point> corners ;
        for(unsigned vertex : piece.vertices) corners.push_back(result.vertices[vertex]) ;
        piece.geometry = std::make_shared<CornerMap>(piece.shape,std::move(corners)) ;
        piece.boundaries.resize(referenceElement(piece.shape).boundaryCount) ;
        piece.boundaries[0].sourceFace = f ;
        piece.boundaries[0].sourcePatch = p ;
        result.pieces.push_back(std::move(piece)) ;
      }
    }
    struct Use { unsigned piece, boundary ; Loop loop ; } ;
    std::map<Loop,std::vector<Use> > interior ;
    for(unsigned p = 0; p < result.pieces.size(); ++p) {
      const GeometryPiece &piece = result.pieces[p] ;
      const ReferenceElement &element = referenceElement(piece.shape) ;
      for(unsigned f = 1; f < element.boundaryCount; ++f) {
        Loop loop ;
        const Boundary &face = element.boundaries[f] ;
        for(unsigned v = 0; v < face.vertexCount; ++v) loop.push_back(piece.vertices[face.vertices[v]]) ;
        interior[sorted(loop)].push_back({p,f,loop}) ;
      }
    }
    for(const auto &face : interior) {
      const auto &uses = face.second ;
      if(uses.size() != 2 || !opposite(uses[0].loop,uses[1].loop))
        throw std::invalid_argument("composite piece interfaces require two opposite uses") ;
      for(unsigned i = 0; i < 2; ++i) {
        PieceBoundary &boundary = result.pieces[uses[i].piece].boundaries[uses[i].boundary] ;
        boundary.neighborPiece = uses[1-i].piece ;
        boundary.neighborBoundary = uses[1-i].boundary ;
      }
    }
    return result ;
  }

}
}
