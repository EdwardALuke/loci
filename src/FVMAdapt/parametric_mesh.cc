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

#include <FVMAdapt/parametric_mesh.h>
#include <FVMAdapt/parametric/native_cell.h>
#include "library/parametric/exact_geometry.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {
    const BoundaryLimits meshLimits = {1000000,8000000} ;
    unsigned faceIndex(const CellTopology &topology, EntityId face) {
      for(unsigned f = 0; f < referenceElement(topology.shape()).boundaryCount; ++f)
        if(topology.faceId(f) == face) return f ;
      throw std::invalid_argument("face is absent from parametric root") ;
    }
    std::map<TreeNodeId,int> leafNumbers(const AffineSubdivision &history) {
      std::map<TreeNodeId,int> result ;
      const auto leaves = history.history().leafIds() ;
      if(leaves.size() > std::size_t(std::numeric_limits<int>::max())) throw std::length_error("too many parametric leaves") ;
      for(std::size_t i = 0; i < leaves.size(); ++i) result.emplace(leaves[i],int(i)) ;
      return result ;
    }
    std::vector<int> leafDepths(const AffineSubdivision &history) {
      std::vector<int> result ;
      for(auto id: history.history().leafIds()) result.push_back(history.history().path(id).size()) ;
      return result ;
    }
    std::vector<EntityId> faceVertices(const CellTopology &topology, unsigned face) {
      const auto &boundary = referenceElement(topology.shape()).boundaries[face] ;
      std::vector<EntityId> result(boundary.vertexCount) ;
      for(unsigned v = 0; v < boundary.vertexCount; ++v)
        result[orientedVertex(boundary.shape,topology.faceOrientation(face),v)] = topology.vertexId(boundary.vertices[v]) ;
      return result ;
    }
    std::vector<BoundaryUse> faceCover(const CellTopology &topology, const AffineSubdivision &history, unsigned face) {
      const Shape frame = referenceElement(topology.shape()).boundaries[face].shape ;
      std::vector<BoundaryUse> result ; unsigned index = 0 ;
      for(auto id: history.history().leafIds()) {
        const auto &child = history.node(id).cell ; const auto &element = referenceElement(child.shape) ;
        for(unsigned b = 0; b < element.boundaryCount; ++b) {
          const auto &side = element.boundaries[b] ; std::vector<ReferencePoint> points ;
          for(unsigned v = 0; v < side.vertexCount; ++v)
            points.push_back(child.placement.apply(referenceVertex(child.shape,side.vertices[v]))) ;
          if(!std::all_of(points.begin(),points.end(),[&](const ReferencePoint &q) { return onBoundary(topology.shape(),face,q) ; })) continue ;
          for(auto &q: points) q = orientReferencePoint(frame,topology.faceOrientation(face),boundaryCoordinates(topology.shape(),face,q)) ;
          Placement placement ; placement.origin = points.front() ;
          for(unsigned d = 0; d < 3; ++d) {
            placement.columns[0][d] = points[1][d]-points[0][d] ;
            placement.columns[1][d] = points.back()[d]-points[0][d] ; placement.columns[2][d] = 0 ;
          }
          result.push_back({{index,b},side.shape,placement}) ;
        }
        ++index ;
      }
      return result ;
    }
    VertexKey faceKey(const CellTopology &topology, unsigned face, const ReferencePoint &canonical) {
      const Shape shape = referenceElement(topology.shape()).boundaries[face].shape ;
      const auto local = orientReferencePoint(shape,inverseOrientation(shape,topology.faceOrientation(face)),canonical) ;
      return topology.vertexKey(boundaryReferencePoint(topology.shape(),face,local)) ;
    }

    /// A leaf touching an original face at a vertex or along an edge.
    struct LeafContact { int leaf ; std::vector<VertexKey> nodes ; } ;

    std::vector<int> gradeFaces(const std::vector<ParametricFace> &faces, const std::vector<int> &depths,
                                unsigned difference, NeighborGrading grading, const std::vector<LeafContact> &contacts = {}) {
      std::set<int> required ;
      const auto grade = [&](const std::set<int> &cells) {
        int maximum = 0 ; for(int cell: cells) maximum = std::max(maximum,depths.at(cell)) ;
        for(int cell: cells) if(unsigned(maximum-depths.at(cell)) > difference) required.insert(cell) ;
      } ;
      std::map<VertexKey,std::set<int>> vertices ;
      std::map<std::pair<VertexKey,VertexKey>,std::set<int>> edges ;
      for(const auto &face: faces) {
        // Parallel node records require signed-int numerators and denominators.
        for(const auto &vertex: face.vertices) for(const auto &coordinate: vertex.point) {
          checkedExactInteger<int>(coordinate.numerator()) ;
          checkedExactInteger<int>(coordinate.denominator()) ;
        }
        std::set<int> cells{face.left} ; if(face.right >= 0) cells.insert(face.right) ; grade(cells) ;
        for(std::size_t v = 0; v < face.vertices.size(); ++v) {
          if(grading == NeighborGrading::vertices) vertices[face.vertices[v]].insert(cells.begin(),cells.end()) ;
          if(grading == NeighborGrading::edges)
            edges[std::minmax(face.vertices[v],face.vertices[(v+1)%face.vertices.size()])].insert(cells.begin(),cells.end()) ;
        }
      }
      for(const auto &contact: contacts) {
        if(grading == NeighborGrading::vertices)
          for(const auto &key: contact.nodes) vertices[key].insert(contact.leaf) ;
        if(grading == NeighborGrading::edges && contact.nodes.size() > 1)
          for(std::size_t v = 0; v < contact.nodes.size(); ++v)
            edges[std::minmax(contact.nodes[v],contact.nodes[(v+1)%contact.nodes.size()])].insert(contact.leaf) ;
      }
      for(const auto &entry: vertices) grade(entry.second) ;
      for(const auto &entry: edges) grade(entry.second) ;
      return {required.begin(),required.end()} ;
    }
    std::map<ReferencePoint,VertexKey> referenceNodes(const CellTopology &topology, const std::vector<VertexKey> &nodes) {
      std::map<ReferencePoint,VertexKey> points ;
      for(const auto &key: nodes) {
        ReferencePoint point ;
        if(!topology.locateVertex(key,point)) throw std::invalid_argument("supplied boundary node is absent from parametric root") ;
        points.emplace(point,key) ;
      }
      return points ;
    }
    // A remote refinement point belongs on every emitted polygon edge that contains it.
    void insertNodes(std::vector<VertexKey> &loop, const CellTopology &topology,
                     const std::map<ReferencePoint,VertexKey> &points) {
      if(points.empty()) return ;
      std::vector<VertexKey> result ;
      for(std::size_t e = 0; e < loop.size(); ++e) {
        ReferencePoint first, second ;
        if(!topology.locateVertex(loop[e],first) ||
           !topology.locateVertex(loop[(e+1)%loop.size()],second))
          throw std::logic_error("emitted vertex is absent from parametric root") ;
        unsigned axis = 0 ; while(axis < 3 && first[axis] == second[axis]) ++axis ;
        if(axis == 3) throw std::logic_error("zero-length parametric polygon edge") ;
        std::map<detail::Exact,VertexKey> interior ;
        for(const auto &entry : points) {
          const auto &point = entry.first ;
          bool outside = false ;
          for(unsigned d = 0; d < 3; ++d)
            if(point[d] < std::min(first[d],second[d]) || point[d] > std::max(first[d],second[d])) { outside = true ; break ; }
          if(outside) continue ;
          const detail::Exact t = (detail::exact(point[axis])-detail::exact(first[axis])) /
                                  (detail::exact(second[axis])-detail::exact(first[axis])) ;
          if(t <= 0 || t >= 1) continue ;
          bool on = true ;
          for(unsigned d = 0; d < 3; ++d)
            if(detail::exact(point[d]) != detail::exact(first[d])+t*(detail::exact(second[d])-detail::exact(first[d]))) {
              on = false ; break ;
            }
          if(on) interior.emplace(t,entry.second) ;
        }
        result.push_back(loop[e]) ;
        for(const auto &entry : interior) result.push_back(entry.second) ;
        if(result.size() > meshLimits.vertices) throw std::length_error("parametric polygon vertex limit exceeded") ;
      }
      loop.swap(result) ;
    }
    std::vector<LeafContact> faceContacts(EntityId face, const CellTopology &topology,
        const AffineSubdivision &history, int offset) {
      std::vector<LeafContact> contacts ; const auto local = faceIndex(topology,face) ;
      for(auto id: history.history().leafIds()) {
        const auto &child = history.node(id).cell ; const auto &element = referenceElement(child.shape) ;
        std::vector<VertexKey> keys ; std::vector<bool> onFace ;
        for(unsigned v = 0; v < element.vertexCount; ++v) {
          const auto point = child.placement.apply(referenceVertex(child.shape,v)) ;
          keys.push_back(topology.vertexKey(point)) ; onFace.push_back(onBoundary(topology.shape(),local,point)) ;
          if(onFace.back()) contacts.push_back({offset,{keys.back()}}) ;
        }
        for(unsigned e = 0; e < element.edgeCount; ++e) {
          const auto a = element.edges[e][0], b = element.edges[e][1] ;
          if(onFace[a] && onFace[b]) contacts.push_back({offset,{keys[a],keys[b]}}) ;
        }
        ++offset ;
      }
      return contacts ;
    }

    // Express a root-coordinate face piece in the leaf face's own coordinates.
    std::vector<ReferencePoint> leafFaceCoordinates(const Child &cell,unsigned face,const std::vector<ReferencePoint> &points) {
      const auto &side=referenceElement(cell.shape).boundaries[face] ;
      const auto a=cell.placement.apply(referenceVertex(cell.shape,side.vertices[0])),
                 b=cell.placement.apply(referenceVertex(cell.shape,side.vertices[1])),
                 c=cell.placement.apply(referenceVertex(cell.shape,side.vertices[side.vertexCount-1])) ;
      for(unsigned x=0;x<3;++x) for(unsigned y=x+1;y<3;++y)
        if((b[x]-a[x])*(c[y]-a[y])!=(b[y]-a[y])*(c[x]-a[x])) {
          Placement placement ; placement.origin={a[x],a[y],0} ;
          placement.columns={{{b[x]-a[x],b[y]-a[y],0},{c[x]-a[x],c[y]-a[y],0},{0,0,0}}} ;
          const BoundaryUse use{{0,face},side.shape,placement} ; std::vector<ReferencePoint> local ;
          for(const auto &p:points) local.push_back(use.coordinates({p[x],p[y],0})) ;
          return local ;
        }
      throw std::logic_error("degenerate leaf face in directional balancing") ;
    }
    unsigned faceDirections(Shape shape,unsigned face,const std::vector<ReferencePoint> &points) {
      if(shape!=Shape::hex && shape!=Shape::prism) return 0 ;
      const auto &side=referenceElement(shape).boundaries[face] ;
      auto low=referenceVertex(shape,side.vertices[0]), high=low ;
      for(unsigned v=1;v<side.vertexCount;++v) for(unsigned d=0;d<3;++d) {
        const auto q=referenceVertex(shape,side.vertices[v])[d] ; low[d]=std::min(low[d],q) ; high[d]=std::max(high[d],q) ;
      }
      auto partLow=high, partHigh=low ;
      for(const auto &p:points) {
        const auto q=boundaryReferencePoint(shape,face,p) ;
        for(unsigned d=0;d<3;++d) { partLow[d]=std::min(partLow[d],q[d]) ; partHigh[d]=std::max(partHigh[d],q[d]) ; }
      }
      unsigned directions=0 ;
      for(unsigned d=0;d<3;++d) if(partLow[d]>low[d] || partHigh[d]<high[d]) directions|=1u<<d ;
      if(shape==Shape::prism && (directions&3)) directions|=3 ;
      return directions ;
    }
    struct FacePieces { unsigned count=0, directions=0 ; } ;
    using DividedFaces=std::map<std::pair<int,int>,FacePieces> ;
    std::vector<FaceBalanceRequest> dividedFaceRequests(EntityId root,const DividedFaces &pieces,bool directions) {
      std::vector<FaceBalanceRequest> out ;
      for(const auto &p:pieces) if(p.second.count>1) {
        // A division without a restricted coordinate range needs the full split.
        const auto axes=directions && !p.second.directions ? 7:p.second.directions ;
        out.push_back({int(root),p.first.first,p.first.second,axes}) ;
      }
      return out ;
    }
  }

  SplitOperation parametricGradingSplit(const GeometryMap &root, int axis) {
    if(axis == -1) return isotropicSplit(root.shape()) ;
    if(axis < 0 || axis > 2) throw std::invalid_argument("parametricExtrusionAxis must be -1, 0, 1 or 2") ;
    if(root.shape() != Shape::hex && root.shape() != Shape::prism)
      throw std::invalid_argument("parametric extrusion requires hexes or prisms") ;
    if(!dynamic_cast<const CornerMap *>(&root))
      throw std::invalid_argument("parametric extrusion requires a native corner map") ;
    const auto &element = referenceElement(root.shape()) ;
    std::vector<Point> corners ; double scale = 0 ;
    for(unsigned v = 0; v < element.vertexCount; ++v) {
      corners.push_back(root.evaluate(element.vertices[v]).position) ;
      for(unsigned d = 0; d < 3; ++d) scale = std::max(scale,std::abs(corners.back()[d]-corners.front()[d])) ;
    }
    const double tolerance = 64*std::numeric_limits<double>::epsilon()*scale ;
    for(unsigned direction = root.shape() == Shape::prism ? 2 : 0; direction < 3; ++direction) {
      bool matches = true ; double thickness = 0, plane = 0 ; bool first = true ;
      for(unsigned v = 0; v < element.vertexCount; ++v) if(element.vertices[v][direction] == 0) {
        unsigned other = 0 ; Point opposite = element.vertices[v] ; opposite[direction] = 1 ;
        while(other < element.vertexCount && element.vertices[other] != opposite) ++other ;
        if(other == element.vertexCount) { matches = false ; break ; }
        const double length = corners[other][axis]-corners[v][axis] ;
        if(first) { thickness = length ; plane = corners[v][axis] ; first = false ; }
        if(thickness == 0 || !std::isfinite(thickness) || std::abs(length-thickness) > tolerance ||
           std::abs(corners[v][axis]-plane) > tolerance) matches = false ;
        for(unsigned d = 0; d < 3; ++d)
          if(int(d) != axis && std::abs(corners[v][d]-corners[other][d]) > tolerance) matches = false ;
      }
      if(!matches) continue ;
      if(root.shape() == Shape::prism) { PrismSplit split ; split.base = TriangleSplit::uniform4() ; return split ; }
      TensorSplit split ;
      for(unsigned d = 0; d < 3; ++d) if(d != direction) split.cuts[d] = {Fraction(1,2)} ;
      return split ;
    }
    throw std::invalid_argument("cell is not a straight extrusion along parametricExtrusionAxis") ;
  }

  SplitOperation parametricGradingSplit(const AffineSubdivision &history, int axis) {
    const auto &saved = history.history() ;
    const auto root = saved.geometry().map(saved.node(0).map) ;
    auto request = parametricGradingSplit(*root,axis) ;
    if(axis == -1) return request ;
    unsigned direction = 2 ;
    if(const auto *tensor = std::get_if<TensorSplit>(&request))
      for(unsigned d = 0; d < 3; ++d) if(tensor->cuts[d].empty()) direction = d ;
    for(auto id: saved.leafIds()) {
      const auto &cell = history.node(id).cell ; bool lower = false, upper = false ;
      for(unsigned v = 0; v < referenceElement(cell.shape).vertexCount; ++v) {
        const auto q = cell.placement.apply(referenceVertex(cell.shape,v))[direction] ;
        if(q == Fraction(0)) lower = true ;
        else if(q == Fraction(1)) upper = true ;
        else throw std::invalid_argument("parametric split divides the protected extrusion direction") ;
      }
      if(!lower || !upper) throw std::invalid_argument("parametric child does not span the extrusion thickness") ;
    }
    return request ;
  }


  NeighborGrading parametricGrading(const std::string &name) {
    if(name == "faces") return NeighborGrading::faces ;
    if(name == "edges") return NeighborGrading::edges ;
    if(name == "vertices") return NeighborGrading::vertices ;
    throw std::invalid_argument("parametricGrading must be faces, edges or vertices") ;
  }

  CellTopology parametricCellTopology(EntityId cell, const std::vector<EntityId> &faces,
      const std::vector<std::vector<EntityId>> &outwardFaces,
      const std::vector<EntityId> &edges, const std::vector<std::array<EntityId,2>> &edgeEnds) {
    if(faces.size() != outwardFaces.size() || edges.size() != edgeEnds.size())
      throw std::invalid_argument("parametric topology ID and connectivity counts differ") ;
    std::optional<NativeCell> recognized ;
    try { recognized=recognizeNativeCell(outwardFaces) ; }
    catch(const std::invalid_argument &) { /* A general shell is checked below. */ }
    std::map<std::array<EntityId,2>,EntityId> byEnds ; std::map<EntityId,std::array<EntityId,2>> byId ;
    for(std::size_t i = 0; i < edges.size(); ++i) {
      auto ends = edgeEnds[i] ; std::sort(ends.begin(),ends.end()) ;
      const auto a = byEnds.emplace(ends,edges[i]) ;
      const auto b = byId.emplace(edges[i],ends) ;
      if(ends[0] == ends[1] || (!a.second && a.first->second != edges[i]) || (!b.second && b.first->second != ends))
        throw std::invalid_argument("parametric edge IDs and endpoints disagree") ;
    }
    if(!recognized) {
      std::set<EntityId> used ; for(const auto &face:outwardFaces) used.insert(face.begin(),face.end()) ;
      std::vector<EntityId> vertices(used.begin(),used.end()),sortedEdges,sortedFaces ;
      std::map<EntityId,unsigned> index ; for(unsigned v=0;v<vertices.size();++v) index.emplace(vertices[v],v) ;
      std::vector<std::array<unsigned,2>> localEdges ; std::vector<std::vector<unsigned>> localFaces ;
      for(const auto &edge:byId) {
        sortedEdges.push_back(edge.first) ; localEdges.push_back({{index.at(edge.second[0]),index.at(edge.second[1])}}) ;
      }
      std::map<EntityId,std::vector<EntityId>> ordered ;
      for(unsigned f=0;f<faces.size();++f) if(!ordered.emplace(faces[f],outwardFaces[f]).second)
        throw std::invalid_argument("original cell repeats a face ID") ;
      for(const auto &entry:ordered) {
        sortedFaces.push_back(entry.first) ; std::vector<unsigned> loop ;
        for(auto v:entry.second) loop.push_back(index.at(v)) ;
        if(!loop.empty()) std::rotate(loop.begin(),std::min_element(loop.begin(),loop.end()),loop.end()) ;
        localFaces.push_back(std::move(loop)) ;
      }
      return CellTopology(cell,std::move(vertices),std::move(sortedEdges),std::move(sortedFaces),std::move(localEdges),std::move(localFaces)) ;
    }
    const auto &native=*recognized ;
    const auto &element = referenceElement(native.shape) ; std::vector<EntityId> nativeEdges, nativeFaces ;
    for(unsigned e = 0; e < element.edgeCount; ++e) {
      std::array<EntityId,2> ends{{native.vertices[element.edges[e][0]],native.vertices[element.edges[e][1]]}} ;
      std::sort(ends.begin(),ends.end()) ; const auto found = byEnds.find(ends) ;
      if(found == byEnds.end()) throw std::invalid_argument("parametric root edge is missing") ;
      nativeEdges.push_back(found->second) ;
    }
    for(unsigned f: native.faces) nativeFaces.push_back(faces[f]) ;
    return CellTopology(native.shape,cell,native.vertices,nativeEdges,nativeFaces) ;
  }

  std::vector<int> encodeParametricTopology(const CellTopology &topology) {
    std::vector<int> result=topology.isStandard() ? std::vector<int>{1,int(topology.shape())}:std::vector<int>{2} ;
    const auto append = [&](EntityId id) {
      if(id > EntityId(std::numeric_limits<int>::max())) throw std::overflow_error("parametric topology ID exceeds int range") ;
      result.push_back(int(id)) ;
    } ;
    append(topology.id()) ;
    if(!topology.isStandard()) { append(topology.vertexCount()) ; append(topology.edgeCount()) ; append(topology.faceCount()) ; }
    for(unsigned v=0;v<topology.vertexCount();++v) append(topology.vertexId(v)) ;
    for(unsigned e=0;e<topology.edgeCount();++e) {
      append(topology.edgeId(e)) ; if(!topology.isStandard()) for(auto v:topology.edgeVertices(e)) append(v) ;
    }
    for(unsigned f=0;f<topology.faceCount();++f) {
      append(topology.faceId(f)) ;
      if(!topology.isStandard()) { const auto loop=topology.faceVertices(f) ; append(loop.size()) ; for(auto v:loop) append(v) ; }
    }
    return result ;
  }
  CellTopology decodeParametricTopology(const std::vector<int> &record) {
    if(record.size() < 3 || (record[0] != 1 && record[0] != 2) || std::any_of(record.begin(),record.end(),[](int id) { return id < 0 ; }))
      throw std::invalid_argument("invalid parametric topology record") ;
    if(record[0]==2) {
      std::size_t cursor=1 ;
      const auto take=[&]() {
        if(cursor==record.size()) throw std::invalid_argument("truncated cell topology") ;
        return unsigned(record[cursor++]) ;
      } ;
      const auto id=take(),nv=take(),ne=take(),nf=take() ;
      if(std::size_t(nv)+3*std::size_t(ne)+5*std::size_t(nf)>record.size()-cursor)
        throw std::invalid_argument("cell topology counts exceed record length") ;
      std::vector<EntityId> vertices(nv),edges(ne),faces(nf) ;
      std::vector<std::array<unsigned,2>> ends(ne) ; std::vector<std::vector<unsigned>> loops(nf) ;
      for(auto &v:vertices) v=take() ;
      for(unsigned e=0;e<ne;++e) { edges[e]=take() ; for(auto &v:ends[e]) v=take() ; }
      for(unsigned f=0;f<nf;++f) {
        faces[f]=take() ; const auto count=take() ;
        if(count<3 || count>record.size()-cursor) throw std::invalid_argument("invalid cell topology face size") ;
        loops[f].resize(count) ; for(auto &v:loops[f]) v=take() ;
      }
      if(cursor!=record.size()) throw std::invalid_argument("trailing cell topology data") ;
      return CellTopology(id,std::move(vertices),std::move(edges),std::move(faces),std::move(ends),std::move(loops)) ;
    }
    const auto shape = Shape(record[1]) ;
    if(shape != Shape::tet && shape != Shape::pyramid && shape != Shape::prism && shape != Shape::hex)
      throw std::invalid_argument("parametric topology requires a native volume cell") ;
    const auto &element = referenceElement(shape) ;
    if(record.size() != 3+element.vertexCount+element.edgeCount+element.boundaryCount)
      throw std::invalid_argument("invalid parametric topology record length") ;
    auto v = record.begin()+3, e = v+element.vertexCount, f = e+element.edgeCount ;
    return CellTopology(shape,EntityId(record[2]),{v,e},{e,f},{f,record.end()}) ;
  }

  std::vector<int> parametricLocalRequests(const CellTopology &topology, const AffineSubdivision &history,
      unsigned difference, NeighborGrading grading) {
    const auto boundary = assembleBoundaries({{topology,history}},meshLimits) ;
    const auto numbers = leafNumbers(history) ; std::vector<ParametricFace> faces ;
    for(const auto &face: boundary.faces)
      faces.push_back({face.vertices,numbers.at(face.sides[0].cell.leaf),
        face.sides.size() == 2 ? numbers.at(face.sides[1].cell.leaf) : -1}) ;
    return gradeFaces(faces,leafDepths(history),difference,grading) ;
  }

  std::vector<FaceBalanceRequest> parametricDividedFaces(const CellTopology &topology, const AffineSubdivision &history,bool directions) {
    const auto boundary=assembleBoundaries({{topology,history}},meshLimits) ; const auto numbers=leafNumbers(history) ;
    DividedFaces pieces ;
    for(const auto &face:boundary.faces) for(const auto &side:face.sides)
      ++pieces[{numbers.at(side.cell.leaf),int(side.cell.boundary)}].count ;
    if(directions) for(const auto &face:boundary.faces) for(const auto &side:face.sides) {
      auto &part=pieces.at({numbers.at(side.cell.leaf),int(side.cell.boundary)}) ;
      const auto &cell=history.node(side.cell.leaf).cell ;
      if(part.count>1 && (cell.shape==Shape::hex || cell.shape==Shape::prism))
        part.directions|=faceDirections(cell.shape,side.cell.boundary,leafFaceCoordinates(cell,side.cell.boundary,side.points)) ;
    }
    return dividedFaceRequests(topology.id(),pieces,directions) ;
  }
  std::vector<FaceBalanceRequest> parametricDividedFaces(EntityId face, const CellTopology &left, const AffineSubdivision &lh,
      const CellTopology &right, const AffineSubdivision &rh,bool directions) {
    const auto lf=faceIndex(left,face), rf=faceIndex(right,face) ;
    const auto matches=reconcileBoundary(referenceElement(left.shape()).boundaries[lf].shape,
      faceCover(left,lh,lf),faceCover(right,rh,rf),meshLimits.faces) ;
    DividedFaces a,b ;
    for(const auto &m:matches) { ++a[{m.first.face.child,m.first.face.boundary}].count ; ++b[{m.second.face.child,m.second.face.boundary}].count ; }
    if(directions) {
      const auto ll=lh.leafIds(), rl=rh.leafIds() ;
      const auto collect=[&](DividedFaces &parts,const BoundaryUse &side,Shape shape,const std::vector<ReferencePoint> &points) {
        auto &part=parts.at({side.face.child,side.face.boundary}) ;
        if(part.count<=1 || (shape!=Shape::hex && shape!=Shape::prism)) return ;
        std::vector<ReferencePoint> local ; for(const auto &p:points) local.push_back(side.coordinates(p)) ;
        part.directions|=faceDirections(shape,side.face.boundary,local) ;
      } ;
      for(const auto &m:matches) {
        collect(a,m.first,lh.node(ll.at(m.first.face.child)).cell.shape,m.vertices) ;
        collect(b,m.second,rh.node(rl.at(m.second.face.child)).cell.shape,m.vertices) ;
      }
    }
    auto out=dividedFaceRequests(left.id(),a,directions) ; const auto other=dividedFaceRequests(right.id(),b,directions) ;
    out.insert(out.end(),other.begin(),other.end()) ;
    return out ;
  }

  ParametricTrace parametricFaceRequests(EntityId face, const CellTopology &left, const AffineSubdivision &lh,
      const CellTopology &right, const AffineSubdivision &rh, unsigned difference, NeighborGrading grading) {
    auto faces = parametricFaceFaces(face,left,lh,{},&right,&rh) ; auto depths = leafDepths(lh) ;
    const auto other = leafDepths(rh) ; const int offset = depths.size() ; depths.insert(depths.end(),other.begin(),other.end()) ;
    for(auto &piece: faces) piece.right += offset ;
    std::vector<LeafContact> contacts ;
    if(grading != NeighborGrading::faces) {
      contacts = faceContacts(face,left,lh,0) ; const auto otherContacts = faceContacts(face,right,rh,offset) ;
      contacts.insert(contacts.end(),otherContacts.begin(),otherContacts.end()) ;
      std::vector<VertexKey> points ;
      for(const auto &piece: faces) points.insert(points.end(),piece.vertices.begin(),piece.vertices.end()) ;
      for(const auto &contact: contacts) points.insert(points.end(),contact.nodes.begin(),contact.nodes.end()) ;
      const auto coordinates = referenceNodes(left,points) ;
      for(auto &piece: faces) insertNodes(piece.vertices,left,coordinates) ;
      for(auto &contact: contacts) if(contact.nodes.size() == 2) insertNodes(contact.nodes,left,coordinates) ;
    }
    ParametricTrace result ;
    for(int leaf: gradeFaces(faces,depths,difference,grading,contacts))
      result.push_back({int(leaf < offset ? left.id() : right.id()),leaf < offset ? leaf : leaf-offset}) ;
    return result ;
  }

  ParametricTrace parametricEdgeTrace(const CellTopology &topology, const AffineSubdivision &history) {
    const auto &element = referenceElement(topology.shape()) ; ParametricTrace result ; int index = 0 ;
    for(auto id: history.history().leafIds()) {
      const auto &child = history.node(id).cell ; std::vector<VertexKey> keys ;
      for(unsigned v = 0; v < referenceElement(child.shape).vertexCount; ++v)
        keys.push_back(topology.vertexKey(child.placement.apply(referenceVertex(child.shape,v)))) ;
      for(unsigned edge = 0; edge < element.edgeCount; ++edge) {
        const auto a = topology.vertexId(element.edges[edge][0]), b = topology.vertexId(element.edges[edge][1]) ;
        const auto ends = std::minmax(a,b) ;
        std::vector<Fraction> points ;
        for(const auto &key: keys) {
          if(key.kind == EntityKind::edge && key.entity == topology.edgeId(edge)) points.push_back(key.point[0]) ;
          if(key.kind == EntityKind::vertex && (key.entity == ends.first || key.entity == ends.second)) points.push_back(key.entity == ends.first ? 0 : 1) ;
        }
        if(points.empty()) continue ;
        const auto range = std::minmax_element(points.begin(),points.end()) ;
        result.push_back({int(topology.edgeId(edge)),int(topology.id()),index,int(history.history().path(id).size()),
          checkedExactInteger<int>(range.first->numerator()),checkedExactInteger<int>(range.first->denominator()),
          checkedExactInteger<int>(range.second->numerator()),checkedExactInteger<int>(range.second->denominator())}) ;
      }
      ++index ;
    }
    return result ;
  }

  std::vector<int> parametricEdgeRequests(EntityId root, const ParametricTrace &trace, unsigned difference, NeighborGrading grading) {
    if(grading == NeighborGrading::faces) return {} ;
    std::set<int> required ;
    for(const auto &own: trace) if(own.at(1) == int(root)) for(const auto &other: trace) {
      if(own[0] != other[0] || other[3] <= own[3] || unsigned(other[3]-own[3]) <= difference) continue ;
      const auto lower = std::max(Fraction(own[4],own[5]),Fraction(other[4],other[5])) ;
      const auto upper = std::min(Fraction(own[6],own[7]),Fraction(other[6],other[7])) ;
      if(lower < upper || (grading == NeighborGrading::vertices && lower == upper)) required.insert(own[2]) ;
    }
    return {required.begin(),required.end()} ;
  }

  std::vector<int> parametricVertexDepths(const CellTopology &topology, const AffineSubdivision &history) {
    std::vector<int> result(referenceElement(topology.shape()).vertexCount,0) ;
    for(auto id: history.history().leafIds()) {
      const auto &child = history.node(id).cell ; const int depth = history.history().path(id).size() ;
      for(unsigned v = 0; v < referenceElement(child.shape).vertexCount; ++v) {
        const auto key = topology.vertexKey(child.placement.apply(referenceVertex(child.shape,v))) ;
        if(key.kind == EntityKind::vertex) for(unsigned i = 0; i < result.size(); ++i)
          if(topology.vertexId(i) == key.entity) result[i] = std::max(result[i],depth) ;
      }
    }
    return result ;
  }
  std::vector<int> parametricVertexRequests(const CellTopology &topology, const AffineSubdivision &history,
      const std::map<EntityId,int> &depths, unsigned difference) {
    std::set<int> result ; int index = 0 ;
    for(auto id: history.history().leafIds()) {
      const auto &child = history.node(id).cell ; const int depth = history.history().path(id).size() ;
      for(unsigned v = 0; v < referenceElement(child.shape).vertexCount; ++v) {
        const auto key = topology.vertexKey(child.placement.apply(referenceVertex(child.shape,v))) ;
        if(key.kind == EntityKind::vertex && depths.at(key.entity) > depth && unsigned(depths.at(key.entity)-depth) > difference)
          result.insert(index) ;
      }
      ++index ;
    }
    return {result.begin(),result.end()} ;
  }

  std::vector<VertexKey> parametricCellNodes(const CellTopology &topology, const AffineSubdivision &history) {
    std::set<VertexKey> result ;
    for(const auto &face: assembleBoundaries({{topology,history}},meshLimits).faces)
      for(const auto &key: face.vertices) if(key.kind == EntityKind::cell) result.insert(key) ;
    return {result.begin(),result.end()} ;
  }
  std::vector<VertexKey> parametricFaceNodes(EntityId face, const CellTopology &left,
      const AffineSubdivision &lh, const CellTopology *right, const AffineSubdivision *rh) {
    std::set<VertexKey> result ;
    for(const auto &piece: parametricFaceFaces(face,left,lh,{},right,rh))
      for(const auto &key: piece.vertices) if(key.kind == EntityKind::face) result.insert(key) ;
    return {result.begin(),result.end()} ;
  }
  std::vector<ParametricFace> parametricCellFaces(const CellTopology &topology,
      const AffineSubdivision &history, const std::vector<VertexKey> &boundaryNodes) {
    const auto assembled = assembleBoundaries({{topology,history}},meshLimits) ; const auto numbers = leafNumbers(history) ;
    std::vector<ParametricFace> result ; const auto points = referenceNodes(topology,boundaryNodes) ;
    for(const auto &face: assembled.faces) {
      if(face.sides.size() != 2) continue ;
      ParametricFace output{face.vertices,numbers.at(face.sides[0].cell.leaf),numbers.at(face.sides[1].cell.leaf)} ;
      if(face.sides[0].reversed) std::reverse(output.vertices.begin(),output.vertices.end()) ;
      insertNodes(output.vertices,topology,points) ; result.push_back(std::move(output)) ;
    }
    return result ;
  }
  std::vector<ParametricFace> parametricFaceFaces(EntityId face, const CellTopology &left,
      const AffineSubdivision &lh, const std::vector<VertexKey> &edgeNodes,
      const CellTopology *right, const AffineSubdivision *rh) {
    if((right == nullptr) != (rh == nullptr)) throw std::invalid_argument("parametric neighbor needs topology and history") ;
    const unsigned original = faceIndex(left,face) ; const auto &boundary = referenceElement(left.shape()).boundaries[original] ;
    const auto first = faceCover(left,lh,original) ; auto second = first ;
    if(right) {
      const unsigned other = faceIndex(*right,face) ;
      if(faceVertices(left,original) != faceVertices(*right,other) ||
         (left.faceOrientation(original) < boundary.vertexCount) == (right->faceOrientation(other) < boundary.vertexCount))
        throw std::invalid_argument("shared parametric face connectivity or orientation differs") ;
      second = faceCover(*right,*rh,other) ;
    }
    const auto matches = reconcileBoundary(boundary.shape,first,second,meshLimits.faces) ;
    std::vector<ParametricFace> result ; auto points = edgeNodes ;
    for(const auto &match: matches) {
      ParametricFace output{{},int(match.first.face.child),right ? int(match.second.face.child) : -1} ;
      for(const auto &point: match.vertices) output.vertices.push_back(faceKey(left,original,point)) ;
      if(left.faceOrientation(original) >= boundary.vertexCount) std::reverse(output.vertices.begin(),output.vertices.end()) ;
      points.insert(points.end(),output.vertices.begin(),output.vertices.end()) ; result.push_back(std::move(output)) ;
    }
    std::sort(points.begin(),points.end()) ; points.erase(std::unique(points.begin(),points.end()),points.end()) ;
    const auto coordinates = referenceNodes(left,points) ;
    for(auto &piece: result) insertNodes(piece.vertices,left,coordinates) ;
    std::sort(result.begin(),result.end(),[](const ParametricFace &a, const ParametricFace &b) { return a.vertices < b.vertices ; }) ;
    return result ;
  }
}
}
