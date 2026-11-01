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

#include "subdivision_internal.h"
#include "exact_geometry.h"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {
namespace detail {
  namespace {

    struct UnpairedFace {
      FaceUse face ;
      Shape shape ;
      std::vector<ReferencePoint> vertices ;
    } ;

    unsigned matchingOrientation(Shape shape, const std::vector<ReferencePoint> &first,
                                 const std::vector<ReferencePoint> &second) {
      for(unsigned orientation = 0; orientation < orientationCount(shape); ++orientation) {
        bool matches = true ;
        for(unsigned v = 0; v < first.size(); ++v)
          if(first[v] != second[orientedVertex(shape,orientation,v)]) matches = false ;
        if(matches) return orientation ;
      }
      throw std::logic_error("subdivision faces have incompatible orientations") ;
    }

    bool edgeCoordinate(const ReferencePoint &point, const ReferencePoint &first,
                        const ReferencePoint &second, Fraction &coordinate) {
      unsigned axis = 0 ;
      while(axis < 3 && first[axis] == second[axis]) ++axis ;
      if(axis == 3) throw std::logic_error("reference edge has coincident endpoints") ;
      const Exact t = (exact(point[axis])-exact(first[axis]))/
                      (exact(second[axis])-exact(first[axis])) ;
      if(t < 0 || t > 1) return false ;
      for(unsigned d = 0; d < 3; ++d)
        if(exact(point[d]) != exact(first[d])+t*(exact(second[d])-exact(first[d]))) return false ;
      coordinate = Fraction(t) ;
      return true ;
    }

  }

  void completeSubdivision(SplitResult &result) {
    const ReferenceElement &parent = referenceElement(result.shape) ;
    result.internalFaces.clear() ; result.boundaryPieces.clear() ; result.edgePieces.clear() ;
    std::map<std::vector<ReferencePoint>,UnpairedFace> pending ;
    for(unsigned child = 0; child < result.children.size(); ++child) {
      const Child &current = result.children[child] ;
      const ReferenceElement &element = referenceElement(current.shape) ;
      std::vector<ReferencePoint> vertices ; vertices.reserve(element.vertexCount) ;
      for(unsigned v = 0; v < element.vertexCount; ++v)
        vertices.push_back(current.placement.apply(referenceVertex(current.shape,v))) ;
      for(unsigned b = 0; b < element.boundaryCount; ++b) {
        const Boundary &face = element.boundaries[b] ;
        std::vector<ReferencePoint> points ; points.reserve(face.vertexCount) ;
        for(unsigned v = 0; v < face.vertexCount; ++v) points.push_back(vertices[face.vertices[v]]) ;
        unsigned boundary = 0 ;
        for(; boundary < parent.boundaryCount; ++boundary) {
          bool belongs = true ;
          for(const ReferencePoint &point : points)
            if(!onBoundary(result.shape,boundary,point)) { belongs = false ; break ; }
          if(belongs) break ;
        }
        if(boundary < parent.boundaryCount) {
          BoundaryPiece piece ; piece.face = {child,b} ; piece.boundary = boundary ;
          for(const ReferencePoint &point : points)
            piece.vertices.push_back(boundaryCoordinates(result.shape,boundary,point)) ;
          result.boundaryPieces.push_back(std::move(piece)) ;
        } else {
          std::vector<ReferencePoint> key = points ; std::sort(key.begin(),key.end()) ;
          auto match = pending.find(key) ;
          if(match == pending.end()) pending.emplace(std::move(key),UnpairedFace{{child,b},face.shape,std::move(points)}) ;
          else {
            if(match->second.shape != face.shape)
              throw std::logic_error("subdivision faces have different shapes") ;
            const unsigned orientation = matchingOrientation(face.shape,match->second.vertices,points) ;
            result.internalFaces.push_back({match->second.face,{child,b},orientation}) ;
            pending.erase(match) ;
          }
        }
      }
      for(unsigned e = 0; e < element.edgeCount; ++e) {
        const ReferencePoint &first = vertices[element.edges[e][0]], &second = vertices[element.edges[e][1]] ;
        for(unsigned parentEdge = 0; parentEdge < parent.edgeCount; ++parentEdge) {
          const ReferencePoint a = referenceVertex(result.shape,parent.edges[parentEdge][0]) ;
          const ReferencePoint b = referenceVertex(result.shape,parent.edges[parentEdge][1]) ;
          Fraction t, u ;
          if(edgeCoordinate(first,a,b,t) && edgeCoordinate(second,a,b,u)) {
            result.edgePieces.push_back({child,e,parentEdge,{{t,u}}}) ;
            break ;
          }
        }
      }
    }
    if(!pending.empty()) throw std::logic_error("subdivision has an unpaired interior face") ;
  }

}
}
}
