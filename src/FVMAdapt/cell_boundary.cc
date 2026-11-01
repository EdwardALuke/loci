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

#include <FVMAdapt/cell_boundary.h>
#include "library/parametric/exact_geometry.h"

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {

    using detail::Exact ;
    using detail::exact ;
    typedef std::array<Exact,4> Plane ;

    // Interior supports are canonical lines or planes in the original reference cell.
    struct Support {
      EntityKind kind ;
      EntityId entity ;
      Plane plane ;
      bool operator<(const Support &other) const {
        return std::tie(kind,entity,plane) < std::tie(other.kind,other.entity,other.plane) ;
      }
    } ;

    struct Cover {
      unsigned rootBoundary ;
      bool reversed ;
      std::vector<BoundaryUse> uses ;
    } ;

    struct Interface {
      Shape shape ;
      std::map<EntityId,Cover> covers ;
    } ;

    Plane supportingPlane(const std::vector<ReferencePoint> &vertices, unsigned dimension) {
      detail::ExactPoint u, v ; Plane plane = {{0,0,0,0}} ;
      for(unsigned d = 0; d < 3; ++d) {
        u[d] = exact(vertices[1][d])-exact(vertices[0][d]) ;
        v[d] = exact(vertices.back()[d])-exact(vertices[0][d]) ;
      }
      if(dimension == 2) { plane[0] = u[1] ; plane[1] = -u[0] ; }
      else
        for(unsigned d = 0; d < 3; ++d)
          plane[d] = u[(d+1)%3]*v[(d+2)%3]-u[(d+2)%3]*v[(d+1)%3] ;
      unsigned first = 0 ; while(first < dimension && plane[first] == 0) ++first ;
      if(first == dimension) throw std::logic_error("degenerate subdivision boundary") ;
      for(unsigned d = 0; d < dimension; ++d) plane[3] -= plane[d]*exact(vertices[0][d]) ;
      const Exact scale = plane[first] ;
      for(Exact &component : plane) component /= scale ;
      return plane ;
    }

    unsigned projectionAxis(const Plane &plane, unsigned dimension) {
      unsigned axis = 0 ;
      for(unsigned d = 1; d < dimension; ++d)
        if(abs(plane[d]) > abs(plane[axis])) axis = d ;
      return axis ;
    }

    ReferencePoint project(const ReferencePoint &point, unsigned axis, unsigned dimension) {
      ReferencePoint result = {{0,0,0}} ; unsigned i = 0 ;
      for(unsigned d = 0; d < dimension; ++d)
        if(d != axis) result[i++] = point[d] ;
      return result ;
    }

    ReferencePoint lift(const ReferencePoint &point, const Plane &plane, unsigned dimension) {
      const unsigned axis = projectionAxis(plane,dimension) ;
      ReferencePoint result = {{0,0,0}} ; unsigned i = 0 ;
      Exact value = -plane[3] ;
      for(unsigned d = 0; d < dimension; ++d)
        if(d != axis) {
          result[d] = point[i++] ; value -= plane[d]*exact(result[d]) ;
        }
      result[axis] = Fraction(value/plane[axis]) ;
      return result ;
    }

    typedef std::pair<EntityKind,EntityId> Entity ;

    void recordEntity(std::map<Entity,std::vector<EntityId>> &entities,
                      std::map<std::pair<EntityKind,std::vector<EntityId>>,EntityId> &ids,
                      Entity key, const std::vector<EntityId> &vertices) {
      const auto found = entities.emplace(key,vertices) ;
      if(!found.second && found.first->second != vertices)
        throw std::invalid_argument("shared entity IDs have inconsistent connectivity") ;
      const auto inverse = ids.emplace(std::make_pair(key.first,vertices),key.second) ;
      if(!inverse.second && inverse.first->second != key.second)
        throw std::invalid_argument("shared connectivity has inconsistent entity IDs") ;
    }

    void checkTopology(const std::vector<CellPartition> &cells) {
      std::set<EntityId> roots ;
      std::map<Entity,std::vector<EntityId>> entities ;
      std::map<std::pair<EntityKind,std::vector<EntityId>>,EntityId> ids ;
      unsigned dimension = 0 ;
      for(const CellPartition &cell : cells) {
        const CellTopology &topology = cell.topology ;
        const Shape shape = topology.shape() ;
        if(shape != Shape::triangle && shape != Shape::quad && shape != Shape::tet &&
           shape != Shape::prism && shape != Shape::pyramid && shape != Shape::hex)
          throw std::invalid_argument("boundary assembly requires a native area or volume cell") ;
        if(cell.subdivision.shape() != shape)
          throw std::invalid_argument("subdivision shape differs from its original topology") ;
        const ReferenceElement &element = referenceElement(shape) ;
        if(dimension && dimension != element.dimension)
          throw std::invalid_argument("boundary assembly requires one cell dimension") ;
        dimension = element.dimension ;
        if(!roots.insert(topology.id()).second)
          throw std::invalid_argument("original cell IDs must be unique") ;
        std::vector<EntityId> vertices ;
        for(unsigned v = 0; v < element.vertexCount; ++v)
          vertices.push_back(topology.vertexId(v)) ;
        for(unsigned e = 0; e < element.edgeCount; ++e) {
          std::vector<EntityId> ends = {vertices[element.edges[e][0]],vertices[element.edges[e][1]]} ;
          std::sort(ends.begin(),ends.end()) ;
          recordEntity(entities,ids,{EntityKind::edge,topology.edgeId(e)},ends) ;
        }
        if(dimension == 3)
          for(unsigned b = 0; b < element.boundaryCount; ++b) {
            const Boundary &face = element.boundaries[b] ;
            std::vector<EntityId> loop(face.vertexCount) ;
            const unsigned orientation = topology.faceOrientation(b) ;
            for(unsigned v = 0; v < face.vertexCount; ++v)
              loop[orientedVertex(face.shape,orientation,v)] = vertices[face.vertices[v]] ;
            recordEntity(entities,ids,{EntityKind::face,topology.faceId(b)},loop) ;
          }
      }
    }

    Placement facePlacement(const std::vector<ReferencePoint> &vertices) {
      Placement result ; result.origin = vertices[0] ;
      result.columns = {{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }} ;
      for(unsigned d = 0; d < 3; ++d) {
        result.columns[0][d] = Fraction(exact(vertices[1][d])-exact(vertices[0][d])) ;
        if(vertices.size() > 2) result.columns[1][d] = Fraction(exact(vertices.back()[d])-exact(vertices[0][d])) ;
      }
      return result ;
    }

    bool reverses(Shape shape, const Placement &placement) {
      if(shape == Shape::segment) return placement.columns[0][0] < 0 ;
      const Exact determinant = exact(placement.columns[0][0])*exact(placement.columns[1][1])
                              -exact(placement.columns[0][1])*exact(placement.columns[1][0]) ;
      if(determinant == 0) throw std::logic_error("degenerate boundary placement") ;
      return determinant < 0 ;
    }

    unsigned originalBoundary(Shape rootShape, const std::vector<ReferencePoint> &vertices) {
      const unsigned count = referenceElement(rootShape).boundaryCount ;
      for(unsigned b = 0; b < count; ++b) {
        bool contained = true ;
        for(const ReferencePoint &point : vertices)
          if(!onBoundary(rootShape,b,point)) { contained = false ; break ; }
        if(contained) return b ;
      }
      return count ;
    }

    void collect(const std::vector<CellPartition> &cells, std::vector<LeafUse> &leaves,
                 std::map<Support,Interface> &interfaces) {
      for(unsigned root = 0; root < cells.size(); ++root) {
        const CellPartition &cell = cells[root] ; const Shape rootShape = cell.topology.shape() ;
        const ReferenceElement &rootElement = referenceElement(rootShape) ;
        const unsigned dimension = rootElement.dimension ;
        for(TreeNodeId id : cell.subdivision.leafIds()) {
          const Child &child = cell.subdivision.node(id).cell ;
          const unsigned index = unsigned(leaves.size()) ; leaves.push_back({root,id,0}) ;
          const ReferenceElement &element = referenceElement(child.shape) ;
          for(unsigned b = 0; b < element.boundaryCount; ++b) {
            const Boundary &face = element.boundaries[b] ;
            std::vector<ReferencePoint> rootVertices ;
            for(unsigned v = 0; v < face.vertexCount; ++v)
              rootVertices.push_back(child.placement.apply(referenceVertex(child.shape,face.vertices[v]))) ;
            const unsigned parent = originalBoundary(rootShape,rootVertices) ;
            const bool original = parent < rootElement.boundaryCount ;
            Support support = {EntityKind::cell,cell.topology.id(),{{0,0,0,0}}} ;
            Shape frame = dimension == 2 ? Shape::segment : Shape::quad ;
            unsigned orientation = 0, axis = 0 ;
            if(original) {
              frame = rootElement.boundaries[parent].shape ;
              if(dimension == 3) {
                support.kind = EntityKind::face ; support.entity = cell.topology.faceId(parent) ;
                orientation = cell.topology.faceOrientation(parent) ;
              } else {
                support.kind = EntityKind::edge ; support.entity = cell.topology.edgeId(parent) ;
                orientation = cell.topology.edgeOrientation(parent) ;
              }
            } else {
              support.plane = supportingPlane(rootVertices,dimension) ;
              axis = projectionAxis(support.plane,dimension) ;
            }
            std::vector<ReferencePoint> vertices ;
            for(const ReferencePoint &point : rootVertices)
              vertices.push_back(original ? orientReferencePoint(frame,orientation,boundaryCoordinates(rootShape,parent,point))
                                          : project(point,axis,dimension)) ;
            const Placement placement = facePlacement(vertices) ;
            const bool reversed = reverses(face.shape,placement) ;
            auto found = interfaces.emplace(support,Interface{frame,{}}).first ;
            if(found->second.shape != frame) throw std::invalid_argument("shared boundary shapes disagree") ;
            auto &covers = found->second.covers ;
            const EntityId side = original ? cell.topology.id() : EntityId(reversed) ;
            auto cover = covers.emplace(side,Cover{parent,reversed,{}}).first ;
            if(cover->second.reversed != reversed || cover->second.rootBoundary != parent)
              throw std::invalid_argument("boundary cover has inconsistent incidence") ;
            cover->second.uses.push_back({{index,b},face.shape,placement}) ;
          }
        }
      }
    }

    FaceSide makeSide(const std::vector<CellPartition> &cells, const std::vector<LeafUse> &leaves,
                      FaceUse use, const Support &support, const Cover &cover,
                      const std::vector<ReferencePoint> &points) {
      LeafUse owner = leaves[use.child] ; owner.boundary = use.boundary ;
      const CellTopology &topology = cells[owner.root].topology ;
      const Shape shape = topology.shape() ;
      const unsigned dimension = referenceElement(shape).dimension ;
      FaceSide result ; result.cell = owner ; result.reversed = cover.reversed ;
      // Shared coordinates lift directly into the root, without child normalization.
      for(const ReferencePoint &point : points) {
        if(support.kind == EntityKind::cell) {
          result.points.push_back(lift(point,support.plane,dimension)) ;
        } else {
          const unsigned boundary = cover.rootBoundary ;
          const Shape face = referenceElement(shape).boundaries[boundary].shape ;
          const unsigned orientation = dimension == 3 ? topology.faceOrientation(boundary) : topology.edgeOrientation(boundary) ;
          result.points.push_back(boundaryReferencePoint(shape,boundary,
              orientReferencePoint(face,inverseOrientation(face,orientation),point))) ;
        }
      }
      return result ;
    }

    bool onSegment(const ReferencePoint &point, const ReferencePoint &first, const ReferencePoint &second, unsigned axis) {
      for(unsigned d = 0; d < 3; ++d)
        if(point[d] < std::min(first[d],second[d]) || point[d] > std::max(first[d],second[d])) return false ;
      const Exact coordinate = (exact(point[axis])-exact(first[axis]))/(exact(second[axis])-exact(first[axis])) ;
      for(unsigned d = 0; d < 3; ++d)
        if(exact(point[d])-exact(first[d]) != coordinate*(exact(second[d])-exact(first[d]))) return false ;
      return true ;
    }

    void canonicalLoop(CellFace &face, unsigned dimension) {
      const std::size_t count = face.vertices.size() ;
      std::size_t start = 0 ; bool reverse = false ;
      if(dimension == 2) {
        reverse = face.vertices.back() < face.vertices.front() ;
        if(reverse) start = count-1 ;
      } else {
        start = std::size_t(std::min_element(face.vertices.begin(),face.vertices.end())-face.vertices.begin()) ;
        reverse = face.vertices[(start+count-1)%count] < face.vertices[(start+1)%count] ;
      }
      std::vector<VertexKey> vertices ;
      for(std::size_t i = 0; i < count; ++i)
        vertices.push_back(face.vertices[reverse ? (start+count-i)%count : (start+i)%count]) ;
      face.vertices.swap(vertices) ;
      for(FaceSide &side : face.sides) {
        std::vector<ReferencePoint> points ;
        for(std::size_t i = 0; i < count; ++i)
          points.push_back(side.points[reverse ? (start+count-i)%count : (start+i)%count]) ;
        side.points.swap(points) ; side.reversed = side.reversed != reverse ;
      }
    }

    struct PointIndex {
      std::vector<ReferencePoint> points ;
      std::array<std::vector<std::size_t>,3> order ;

      void sort() {
        std::sort(points.begin(),points.end()) ;
        points.erase(std::unique(points.begin(),points.end()),points.end()) ;
        for(unsigned d = 0; d < 3; ++d) {
          for(std::size_t p = 0; p < points.size(); ++p) order[d].push_back(p) ;
          std::sort(order[d].begin(),order[d].end(),[&](std::size_t a, std::size_t b) {
            return points[a][d] < points[b][d] ;
          }) ;
        }
      }

      std::vector<ReferencePoint> onEdge(const ReferencePoint &first, const ReferencePoint &second, unsigned axis) const {
        typedef std::vector<std::size_t>::const_iterator Iterator ;
        Iterator begin = order[0].begin(), end = order[0].end() ; std::size_t count = points.size() ;
        // Query the narrowest coordinate range before testing exact collinearity.
        for(unsigned d = 0; d < 3; ++d) {
          const auto lower = std::lower_bound(order[d].begin(),order[d].end(),std::min(first[d],second[d]),
              [&](std::size_t p, const Fraction &value) { return points[p][d] < value ; }) ;
          const auto upper = std::upper_bound(lower,order[d].end(),std::max(first[d],second[d]),
              [&](const Fraction &value, std::size_t p) { return value < points[p][d] ; }) ;
          if(std::size_t(upper-lower) < count) {
            begin = lower ; end = upper ; count = std::size_t(upper-lower) ;
          }
        }
        std::vector<ReferencePoint> result ;
        for(auto p = begin; p != end; ++p)
          if(points[*p] != second && onSegment(points[*p],first,second,axis)) result.push_back(points[*p]) ;
        std::sort(result.begin(),result.end(),[&](const ReferencePoint &a, const ReferencePoint &b) {
          return first[axis] < second[axis] ? a[axis] < b[axis] : b[axis] < a[axis] ;
        }) ;
        return result ;
      }
    } ;

    std::vector<PointIndex> indexPoints(const std::vector<CellPartition> &cells, const std::vector<CellFace> &faces) {
      std::map<Entity,std::vector<unsigned>> incident ;
      for(unsigned root = 0; root < cells.size(); ++root) {
        const CellTopology &topology = cells[root].topology ;
        const ReferenceElement &element = referenceElement(topology.shape()) ;
        incident[{EntityKind::cell,topology.id()}].push_back(root) ;
        for(unsigned v = 0; v < element.vertexCount; ++v)
          incident[{EntityKind::vertex,topology.vertexId(v)}].push_back(root) ;
        for(unsigned e = 0; e < element.edgeCount; ++e)
          incident[{EntityKind::edge,topology.edgeId(e)}].push_back(root) ;
        if(element.dimension == 3)
          for(unsigned f = 0; f < element.boundaryCount; ++f)
            incident[{EntityKind::face,topology.faceId(f)}].push_back(root) ;
      }
      std::set<VertexKey> vertices ;
      for(const CellFace &face : faces) vertices.insert(face.vertices.begin(),face.vertices.end()) ;
      std::vector<PointIndex> result(cells.size()) ;
      // The key's smallest original entity identifies every root that can contain it.
      for(const VertexKey &key : vertices)
        for(unsigned root : incident.at({key.kind,key.entity})) {
          ReferencePoint point ;
          if(!cells[root].topology.locateVertex(key,point))
            throw std::logic_error("vertex entity is absent from an incident root") ;
          result[root].points.push_back(point) ;
        }
      for(PointIndex &index : result) index.sort() ;
      return result ;
    }

    void closeEdges(const std::vector<CellPartition> &cells, std::vector<CellFace> &faces, std::size_t limit) {
      const std::vector<PointIndex> points = indexPoints(cells,faces) ;
      std::size_t count = 0 ;
      for(CellFace &face : faces) {
        const unsigned root = face.sides[0].cell.root ;
        const CellTopology &topology = cells[root].topology ;
        const unsigned dimension = referenceElement(topology.shape()).dimension ;
        const std::vector<ReferencePoint> &old = face.sides[0].points ;
        std::vector<VertexKey> loop ;
        const std::size_t edges = dimension == 2 ? old.size()-1 : old.size() ;
        for(std::size_t e = 0; e < edges; ++e) {
          const ReferencePoint &first = old[e], &second = old[(e+1)%old.size()] ;
          unsigned axis = 0 ; while(axis < 3 && first[axis] == second[axis]) ++axis ;
          if(axis == 3) throw std::logic_error("zero-length assembled face edge") ;
          for(const ReferencePoint &point : points[root].onEdge(first,second,axis)) {
            if(count == limit) throw std::length_error("assembled face vertex limit exceeded") ;
            loop.push_back(topology.vertexKey(point)) ; ++count ;
          }
        }
        if(dimension == 2) {
          if(count == limit) throw std::length_error("assembled face vertex limit exceeded") ;
          loop.push_back(topology.vertexKey(old.back())) ; ++count ;
        }
        face.vertices.swap(loop) ;
        for(FaceSide &side : face.sides) {
          side.points.clear() ;
          for(const VertexKey &key : face.vertices) {
            ReferencePoint point ;
            if(!cells[side.cell.root].topology.locateVertex(key,point))
              throw std::logic_error("shared face vertex is absent from an incident root") ;
            side.points.push_back(point) ;
          }
        }
        canonicalLoop(face,dimension) ;
      }
    }

  }

  CellBoundaries assembleBoundaries(const std::vector<CellPartition> &cells, BoundaryLimits limits) {
    checkTopology(cells) ;
    if(cells.size() > std::numeric_limits<unsigned>::max())
      throw std::length_error("too many original cells") ;
    std::size_t leafCount = 0 ;
    for(const CellPartition &cell : cells) {
      const std::size_t count = cell.subdivision.leafIds().size() ;
      if(count > std::numeric_limits<unsigned>::max()-leafCount || count > limits.faces-leafCount)
        throw std::length_error("too many leaves for boundary assembly") ;
      leafCount += count ;
    }
    std::vector<LeafUse> leaves ;
    std::map<Support,Interface> interfaces ;
    collect(cells,leaves,interfaces) ;
    CellBoundaries result ; std::size_t vertexCount = 0 ;
    for(const auto &entry : interfaces) {
      const Interface &interface = entry.second ;
      if(interface.covers.size() > 2 || (entry.first.kind == EntityKind::cell && interface.covers.size() != 2))
        throw std::invalid_argument("interface has an invalid number of incident roots or sides") ;
      const Cover &first = interface.covers.begin()->second ;
      const bool paired = interface.covers.size() == 2 ;
      const Cover &second = paired ? std::next(interface.covers.begin())->second : first ;
      if(paired && first.reversed == second.reversed)
        throw std::invalid_argument("shared boundaries require opposite outward orientations") ;
      const auto matches = reconcileBoundary(interface.shape,first.uses,second.uses,limits.faces-result.faces.size()) ;
      for(const BoundaryMatch &match : matches) {
        if(match.vertices.size() > limits.vertices-vertexCount)
          throw std::length_error("assembled face vertex limit exceeded") ;
        vertexCount += match.vertices.size() ;
        CellFace face ;
        face.sides.push_back(makeSide(cells,leaves,match.first.face,entry.first,first,match.vertices)) ;
        if(paired) face.sides.push_back(makeSide(cells,leaves,match.second.face,entry.first,second,match.vertices)) ;
        const CellTopology &topology = cells[face.sides[0].cell.root].topology ;
        for(const ReferencePoint &point : face.sides[0].points) face.vertices.push_back(topology.vertexKey(point)) ;
        if(paired) {
          const CellTopology &other = cells[face.sides[1].cell.root].topology ;
          for(std::size_t v = 0; v < face.vertices.size(); ++v)
            if(other.vertexKey(face.sides[1].points[v]) != face.vertices[v])
              throw std::invalid_argument("shared boundary vertices disagree") ;
        }
        result.faces.push_back(std::move(face)) ;
      }
    }
    closeEdges(cells,result.faces,limits.vertices) ;
    std::sort(result.faces.begin(),result.faces.end(),[](const CellFace &a, const CellFace &b) { return a.vertices < b.vertices ; }) ;
    std::map<std::pair<unsigned,TreeNodeId>,std::size_t> indices ;
    for(const LeafUse &leaf : leaves) {
      indices[{leaf.root,leaf.leaf}] = result.leaves.size() ;
      result.leaves.push_back({leaf.root,leaf.leaf,{}}) ;
    }
    for(std::size_t f = 0; f < result.faces.size(); ++f)
      for(const FaceSide &side : result.faces[f].sides)
        result.leaves[indices.at({side.cell.root,side.cell.leaf})].faces.push_back({f,side.reversed}) ;
    return result ;
  }

}
}
