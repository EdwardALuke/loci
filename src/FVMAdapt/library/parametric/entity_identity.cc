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

#include <FVMAdapt/parametric/entity_identity.h>

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <map>
#include <set>

namespace Loci {
namespace FVMAdapt {
  bool faceBalanceRequired(const std::vector<std::vector<unsigned>> &faces,
      const std::vector<unsigned> &divided, FaceBalance policy) {
    if(policy==FaceBalance::none) return false ;
    if(policy!=FaceBalance::majority && policy!=FaceBalance::opposite) throw std::invalid_argument("invalid face balance policy") ;
    auto selected=divided ; std::sort(selected.begin(),selected.end()) ;
    selected.erase(std::unique(selected.begin(),selected.end()),selected.end()) ;
    for(auto f:selected) if(f>=faces.size()) throw std::invalid_argument("divided face is outside the cell") ;
    if(selected.size()>faces.size()/2) return true ;
    if(policy==FaceBalance::opposite) for(unsigned i=0;i<selected.size();++i) for(unsigned j=0;j<i;++j) {
      const auto &a=faces[selected[i]], &b=faces[selected[j]] ; bool touching=false ;
      for(auto v:a) touching=touching || std::find(b.begin(),b.end(),v)!=b.end() ;
      if(!touching) return true ;
    }
    return false ;
  }
  namespace {

    void checkIds(const std::vector<EntityId> &ids, unsigned count) {
      if(ids.size() != count)
        throw std::invalid_argument("entity ID count does not match reference topology") ;
      std::vector<EntityId> ordered = ids ;
      std::sort(ordered.begin(),ordered.end()) ;
      if(std::adjacent_find(ordered.begin(),ordered.end()) != ordered.end())
        throw std::invalid_argument("original entity IDs must be unique within each kind") ;
    }

    unsigned canonicalOrientation(unsigned vertices, const std::vector<EntityId> &ids) {
      unsigned chosen = 0 ;
      std::vector<EntityId> least = ids, candidate(ids.size()) ;
      for(unsigned orientation = 1; orientation < 2*vertices; ++orientation) {
        for(unsigned vertex = 0; vertex < ids.size(); ++vertex)
          candidate[orientedVertex(vertices,orientation,vertex)] = ids[vertex] ;
        if(candidate < least) { least = candidate ; chosen = orientation ; }
      }
      return chosen ;
    }

    bool edgeCoordinate(Shape shape, unsigned edge, const ReferencePoint &point, Fraction &coordinate) {
      const auto &vertices = referenceElement(shape).edges[edge] ;
      const ReferencePoint start = referenceVertex(shape,vertices[0]), end = referenceVertex(shape,vertices[1]) ;
      for(unsigned d = 0; d < 3; ++d)
        if(start[d] != end[d]) {
          coordinate = (point[d]-start[d])/(end[d]-start[d]) ;
          break ;
        }
      if(coordinate < 0 || coordinate > 1) return false ;
      for(unsigned d = 0; d < 3; ++d)
        if(start[d]+coordinate*(end[d]-start[d]) != point[d]) return false ;
      return true ;
    }

  }

  bool operator==(const VertexKey &a, const VertexKey &b) {
    return a.kind == b.kind && a.entity == b.entity && a.point == b.point ;
  }

  bool operator!=(const VertexKey &a, const VertexKey &b) { return !(a == b) ; }

  bool operator<(const VertexKey &a, const VertexKey &b) {
    if(a.kind != b.kind) return a.kind < b.kind ;
    if(a.entity != b.entity) return a.entity < b.entity ;
    return a.point < b.point ;
  }

  CellTopology::CellTopology(Shape shape, EntityId cell, std::vector<EntityId> vertices,
                             std::vector<EntityId> edges, std::vector<EntityId> faces)
    : shape_(shape), id_(cell), vertices_(std::move(vertices)),
      edges_(std::move(edges)), faces_(std::move(faces)) {
    const ReferenceElement &element = referenceElement(shape) ;
    checkIds(vertices_,element.vertexCount) ;
    checkIds(edges_,element.dimension >= 2 ? element.edgeCount : 0) ;
    checkIds(faces_,element.dimension == 3 ? element.boundaryCount : 0) ;
    for(unsigned edge = 0; edge < edges_.size(); ++edge)
      edgeOrientations_.push_back(vertices_[element.edges[edge][0]] < vertices_[element.edges[edge][1]] ? 0 : 1) ;
    for(unsigned face = 0; face < faces_.size(); ++face) {
      const Boundary &boundary = element.boundaries[face] ;
      std::vector<EntityId> ids ;
      for(unsigned vertex = 0; vertex < boundary.vertexCount; ++vertex)
        ids.push_back(vertices_[boundary.vertices[vertex]]) ;
      faceOrientations_.push_back(canonicalOrientation(boundary.vertexCount,ids)) ;
    }
  }

  struct CellTopology::Incidence {
    std::vector<std::array<unsigned,2>> edges ;
    std::vector<std::vector<unsigned>> faces ;
  } ;

  CellTopology::CellTopology(EntityId cell, std::vector<EntityId> vertices, std::vector<EntityId> edges,
      std::vector<EntityId> faces, std::vector<std::array<unsigned,2>> edgeVertices,
      std::vector<std::vector<unsigned>> faceVertices)
    : id_(cell),vertices_(std::move(vertices)),edges_(std::move(edges)),faces_(std::move(faces)) {
    checkIds(vertices_,vertices_.size()) ; checkIds(edges_,edgeVertices.size()) ; checkIds(faces_,faceVertices.size()) ;
    if(vertices_.size()<4 || faces_.size()<4 || vertices_.size()+faces_.size()!=edges_.size()+2)
      throw std::invalid_argument("cell boundary must be a closed spherical shell") ;
    std::map<std::pair<unsigned,unsigned>,unsigned> byEnds ;
    for(unsigned e=0;e<edgeVertices.size();++e) {
      const auto a=edgeVertices[e][0],b=edgeVertices[e][1] ;
      if(a>=vertices_.size() || b>=vertices_.size() || a==b || !byEnds.emplace(std::minmax(a,b),e).second)
        throw std::invalid_argument("invalid or repeated cell edge") ;
      edgeOrientations_.push_back(vertices_[a]<vertices_[b] ? 0:1) ;
    }
    std::set<std::vector<unsigned>> uniqueFaces ; std::set<unsigned> usedVertices ;
    std::vector<std::vector<std::pair<unsigned,int>>> uses(edges_.size()) ;
    for(unsigned f=0;f<faceVertices.size();++f) {
      const auto &loop=faceVertices[f] ; auto sorted=loop ; std::sort(sorted.begin(),sorted.end()) ;
      if(loop.size()<3 || sorted.back()>=vertices_.size() ||
         std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end() || !uniqueFaces.insert(sorted).second)
        throw std::invalid_argument("cell requires distinct polygon faces") ;
      std::vector<EntityId> ids ;
      for(unsigned i=0;i<loop.size();++i) {
        const auto a=loop[i],b=loop[(i+1)%loop.size()] ; const auto edge=byEnds.find(std::minmax(a,b)) ;
        if(edge==byEnds.end()) throw std::invalid_argument("cell face edge is missing") ;
        uses[edge->second].push_back({f,a<b ? 1:-1}) ; ids.push_back(vertices_[a]) ; usedVertices.insert(a) ;
      }
      faceOrientations_.push_back(canonicalOrientation(loop.size(),ids)) ;
    }
    std::vector<std::vector<unsigned>> adjacent(faces_.size()) ;
    for(const auto &edge:uses) {
      if(edge.size()!=2 || edge[0].second+edge[1].second!=0)
        throw std::invalid_argument("cell edges require two opposite face uses") ;
      adjacent[edge[0].first].push_back(edge[1].first) ; adjacent[edge[1].first].push_back(edge[0].first) ;
    }
    std::set<unsigned> reached{0} ; std::vector<unsigned> pending{0} ;
    while(!pending.empty()) {
      const auto face=pending.back() ; pending.pop_back() ;
      for(auto other:adjacent[face]) if(reached.insert(other).second) pending.push_back(other) ;
    }
    if(reached.size()!=faces_.size() || usedVertices.size()!=vertices_.size())
      throw std::invalid_argument("cell boundary is disconnected or has unused vertices") ;
    incidence_=std::make_shared<const Incidence>(Incidence{std::move(edgeVertices),std::move(faceVertices)}) ;
  }

  bool CellTopology::operator==(const CellTopology &other) const {
    return shape_==other.shape_ && id_==other.id_ && vertices_==other.vertices_ && edges_==other.edges_ && faces_==other.faces_ &&
      (incidence_==other.incidence_ || (incidence_ && other.incidence_ &&
        incidence_->edges==other.incidence_->edges && incidence_->faces==other.incidence_->faces)) ;
  }

  Shape CellTopology::shape() const {
    if(!shape_) throw std::invalid_argument("cell has no standard reference element") ;
    return *shape_ ;
  }
  std::array<unsigned,2> CellTopology::edgeVertices(unsigned edge) const {
    edges_.at(edge) ;
    if(incidence_) return incidence_->edges.at(edge) ;
    const auto &ends=referenceElement(shape()).edges[edge] ; return {{unsigned(ends[0]),unsigned(ends[1])}} ;
  }
  std::vector<unsigned> CellTopology::faceVertices(unsigned face) const {
    faces_.at(face) ;
    if(incidence_) return incidence_->faces.at(face) ;
    const auto &loop=referenceElement(shape()).boundaries[face] ;
    return {loop.vertices.begin(),loop.vertices.begin()+loop.vertexCount} ;
  }
  Shape CellTopology::faceShape(unsigned face) const {
    faces_.at(face) ;
    return incidence_ ? (incidence_->faces.at(face).size()==3 ? Shape::triangle:
      incidence_->faces.at(face).size()==4 ? Shape::quad:Shape::polygon):referenceElement(shape()).boundaries[face].shape ;
  }

  unsigned CellTopology::edgeOrientation(unsigned edge) const {
    return edgeOrientations_.at(edge) ;
  }

  unsigned CellTopology::faceOrientation(unsigned face) const {
    return faceOrientations_.at(face) ;
  }

  VertexKey CellTopology::vertexKey(const ReferencePoint &point) const {
    if(!containsExact(shape(),point))
      throw std::invalid_argument("vertex is outside the original reference cell") ;
    const ReferenceElement &element = referenceElement(shape()) ;
    for(unsigned vertex = 0; vertex < element.vertexCount; ++vertex)
      if(point == referenceVertex(shape(),vertex))
        return {EntityKind::vertex,vertices_[vertex],{{0,0,0}}} ;
    for(unsigned edge = 0; edge < edges_.size(); ++edge) {
      Fraction coordinate ;
      if(edgeCoordinate(shape(),edge,point,coordinate)) {
        if(edgeOrientations_[edge]) coordinate = 1-coordinate ;
        return {EntityKind::edge,edges_[edge],{{coordinate,0,0}}} ;
      }
    }
    for(unsigned face = 0; face < faces_.size(); ++face)
      if(onBoundary(shape(),face,point)) {
        const ReferencePoint local = boundaryCoordinates(shape(),face,point) ;
        return {EntityKind::face,faces_[face],
                orientReferencePoint(element.boundaries[face].shape,faceOrientations_[face],local)} ;
      }
    return {EntityKind::cell,id_,point} ;
  }

  bool CellTopology::locateVertex(const VertexKey &key, ReferencePoint &point) const {
    const ReferenceElement &element = referenceElement(shape()) ;
    ReferencePoint located ;
    switch(key.kind) {
    case EntityKind::vertex: {
      const auto found = std::find(vertices_.begin(),vertices_.end(),key.entity) ;
      if(found == vertices_.end()) return false ;
      located = referenceVertex(shape(),unsigned(found-vertices_.begin())) ;
      break ;
    }
    case EntityKind::edge: {
      const auto found = std::find(edges_.begin(),edges_.end(),key.entity) ;
      if(found == edges_.end()) return false ;
      const unsigned edge = unsigned(found-edges_.begin()) ;
      const Fraction coordinate = orientReferencePoint(Shape::segment,edgeOrientations_[edge],key.point)[0] ;
      const ReferencePoint first = referenceVertex(shape(),element.edges[edge][0]) ;
      const ReferencePoint second = referenceVertex(shape(),element.edges[edge][1]) ;
      for(unsigned d = 0; d < 3; ++d)
        located[d] = first[d]+coordinate*(second[d]-first[d]) ;
      break ;
    }
    case EntityKind::face: {
      const auto found = std::find(faces_.begin(),faces_.end(),key.entity) ;
      if(found == faces_.end()) return false ;
      const unsigned face = unsigned(found-faces_.begin()) ;
      const Shape shape = element.boundaries[face].shape ;
      const ReferencePoint local = orientReferencePoint(shape,inverseOrientation(shape,faceOrientations_[face]),key.point) ;
      located = boundaryReferencePoint(this->shape(),face,local) ;
      break ;
    }
    case EntityKind::cell:
      if(key.entity != id_) return false ;
      located = key.point ;
      break ;
    default:
      throw std::invalid_argument("unknown vertex entity kind") ;
    }
    if(vertexKey(located) != key)
      throw std::invalid_argument("vertex key is not canonical for its entity") ;
    point = located ;
    return true ;
  }

}
}
