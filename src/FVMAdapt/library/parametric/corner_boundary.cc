// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/corner_boundary.h>
#include <type_traits>
namespace Loci {
namespace FVMAdapt {
  namespace {
    using VertexId = BoundaryId ;
    using EdgeId = BoundaryId ;
    using SurfaceId = BoundaryId ;
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    void checkMethod(CenterMethod method) {
      require(method == CenterMethod::reference || method == CenterMethod::weighted || method == CenterMethod::integrated,
              "unknown vertex volume center policy") ;
    }
    unsigned orientation(Shape shape, const std::vector<VertexId> &local, const std::vector<VertexId> &canonical) {
      for(unsigned o = 0; o < orientationCount(shape); ++o) {
        bool match = true ;
        for(unsigned v = 0; v < local.size(); ++v) match = match && local[v] == canonical[orientedVertex(shape,o,v)] ;
        if(match) return o ;
      }
      throw std::invalid_argument("inconsistent vertex volume face incidence") ;
    }
    std::shared_ptr<const GeometryMap> boundaryMap(const std::shared_ptr<const GeometryMap> &map, unsigned boundary) {
      const auto &face = referenceElement(map->shape()).boundaries[boundary] ;
      const auto &element = referenceElement(face.shape) ;
      std::vector<Point> corners ;
      for(unsigned v = 0; v < face.vertexCount; ++v) corners.push_back(boundaryPoint(map->shape(),boundary,element.vertices[v])) ;
      return std::make_shared<RestrictedMap>(map,std::make_shared<CornerMap>(face.shape,corners)) ;
    }
  }
  BoundaryRecordLimits CornerBoundary::boundaryLimits() const {
    BoundaryRecordLimits limits ; limits.maxEntities = options.maxEntities ;
    limits.maxSplits = options.maxEntities ; limits.maxBytes = options.maxSnapshotBytes ; return limits ;
  }

  RecordLimits CornerBoundary::recordLimits() const {
    RecordLimits limits ; limits.maxNodes = options.maxNodes ; limits.maxSplits = options.maxNodes ;
    limits.maxBytes = options.maxSnapshotBytes ; limits.jacobian = options.jacobianLimits ;
    limits.maxDepth = HistoryLimits().maxDepth ; return limits ;
  }

  void CornerBoundary::checkEntity() const {
    const auto &b = boundaries ; const auto limit = options.maxEntities ;
    if(b.points.size() > limit || b.edges.size() > limit-b.points.size() || b.faces.size() >= limit-b.points.size()-b.edges.size())
      throw std::length_error("vertex volume shared entity limit exceeded") ;
  }

  BoundaryId CornerBoundary::addPoint(const Point &point) {
    checkEntity() ; for(double x: point) require(std::isfinite(x),"nonfinite vertex volume point") ;
    boundaries.points.push_back(point) ; return boundaries.points.size()-1 ;
  }

  BoundaryId CornerBoundary::addEdge(VertexId a, VertexId b, std::shared_ptr<const GeometryMap> geometry, EdgeId parent) {
    if(b < a) {
      std::swap(a,b) ;
      geometry = std::make_shared<RestrictedMap>(geometry,
        std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{{{1,0,0}},{{0,0,0}}})) ;
    }
    checkEntity() ;
    const auto map = boundaries.geometry.capture(geometry) ;
    boundaries.edges.push_back({{{a,b}},map,parent,{},{}}) ; return boundaries.edges.size()-1 ;
  }

  BoundaryId CornerBoundary::addSurface(std::vector<VertexId> vertices, std::vector<EdgeId> edges,
                         std::shared_ptr<const GeometryMap> geometry, int boundary,
                       BoundaryId parent, std::size_t parentSplit) {
    checkEntity() ; const auto map = boundaries.geometry.capture(geometry) ;
    boundaries.faces.push_back({parent,parentSplit,map,std::move(vertices),std::move(edges),{},boundary}) ;
    return boundaries.faces.size()-1 ;
  }

  BoundaryId CornerBoundary::splitEdge(EdgeId id, Transaction &transaction) {
    if(!boundaries.edges.at(id).cuts.empty()) return boundaries.edges[id].cuts.front().vertex ;
    const auto edge = boundaries.edges[id] ; const auto geometry = boundaries.geometry.map(edge.geometry) ;
    const auto middle = addPoint(geometry->evaluate({{.5,0,0}}).position) ;
    const auto segment = [&](double a, double b) {
      return std::make_shared<RestrictedMap>(geometry,
        std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{{{a,0,0}},{{b,0,0}}})) ;
    } ;
    const auto first = addEdge(edge.vertices[0],middle,segment(0,.5),id) ;
    const auto second = addEdge(middle,edge.vertices[1],segment(.5,1),id) ;
    transaction.saveEdge(id) ;
    boundaries.edges[id].cuts = {{Fraction(1,2),middle}} ;
    boundaries.edges[id].children = {{Fraction(0),Fraction(1,2),first},{Fraction(1,2),Fraction(1),second}} ;
    return middle ;
  }

  BoundaryId CornerBoundary::edgeSegment(EdgeId parent, VertexId endpoint) const {
    const auto &edge = boundaries.edges.at(parent) ;
    return edge.children.at(edge.vertices[0] == endpoint ? 0 : 1).edge ;
  }

  std::size_t CornerBoundary::splitSurface(SurfaceId id, CenterMethod method, const std::vector<std::size_t> &active,
      Transaction &transaction, const FaceSplitRecord *accepted) {
    const auto parent = boundaries.faces.at(id) ;
    if(id < active.size() && active[id] != absent) {
      require(active[id]<parent.splits.size(),"selected face construction does not exist") ;
      return active[id] ;
    }
    for(std::size_t i = 0; i < parent.splits.size(); ++i) if(!accepted && parent.splits[i].centerMethod == method) return i ;
    const auto geometry = boundaries.geometry.map(parent.geometry) ;
    const auto center = accepted ? accepted->referenceCenter : faceCenterCoordinates(*geometry,method,options.assessment.measure) ;
    for(std::size_t i = 0; i < parent.splits.size(); ++i) if(parent.splits[i].referenceCenter == center) return i ;
    const auto centerVertex = addPoint(geometry->evaluate(center).position) ;
    const auto coordinates=geometry->referenceCorners() ; const unsigned n=parent.vertices.size() ;
    const auto midpoint=[&](unsigned e) {
      Point q ; for(unsigned d=0;d<3;++d) q[d]=(coordinates[e][d]+coordinates[(e+1)%n][d])/2 ; return q ;
    } ;
    std::vector<VertexId> midpoints ; std::vector<EdgeId> radial ; std::vector<SurfaceId> children ;
    for(auto edge: parent.edges) midpoints.push_back(splitEdge(edge,transaction)) ;
    for(unsigned e = 0; e < n; ++e) {
      const auto mid = midpoint(e) ;
      radial.push_back(addEdge(midpoints[e],centerVertex,std::make_shared<RestrictedMap>(geometry,
        std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{mid,center})))) ;
    }
    for(unsigned v = 0; v < n; ++v) {
      const unsigned prev = (v+n-1)%n ;
      auto placement = std::make_shared<CornerMap>(Shape::quad,std::vector<Point>{coordinates[v],midpoint(v),center,midpoint(prev)}) ;
      children.push_back(addSurface({parent.vertices[v],midpoints[v],centerVertex,midpoints[prev]},
        {edgeSegment(parent.edges[v],parent.vertices[v]),radial[v],radial[prev],edgeSegment(parent.edges[prev],parent.vertices[v])},
        std::make_shared<RestrictedMap>(geometry,placement),parent.boundary,id,parent.splits.size())) ;
    }
    transaction.saveSurface(id) ;
    boundaries.faces[id].splits.push_back({method,center,centerVertex,std::move(children)}) ;
    return parent.splits.size() ;
  }

  void CornerBoundary::edgeVertices(EdgeId id, VertexId start, const std::vector<bool> &split, std::vector<VertexId> &vertices) const {
    const auto &edge = boundaries.edges.at(id) ;
    if(!split[id]) { vertices.push_back(start) ; return ; }
    const bool forward = edge.vertices[0] == start ;
    if(forward) {
      for(std::size_t i = 0; i < edge.children.size(); ++i)
        edgeVertices(edge.children[i].edge,i ? edge.cuts[i-1].vertex : start,split,vertices) ;
    } else {
      for(std::size_t i = edge.children.size(); i-- > 0;)
        edgeVertices(edge.children[i].edge,i+1 == edge.children.size() ? start : edge.cuts[i].vertex,split,vertices) ;
    }
  }

  std::map<std::pair<VertexId,VertexId>,EdgeId> CornerBoundary::cellEdges(const CornerCell &cell) const {
    std::map<std::pair<VertexId,VertexId>,EdgeId> result ;
    for(const auto use: cell.faces) for(auto edge: boundaries.faces.at(use.surface).edges) {
      const auto &ends = boundaries.edges.at(edge).vertices ;
      const auto inserted = result.emplace(endpoints(ends[0],ends[1]),edge) ;
      require(inserted.second || inserted.first->second == edge,"cell faces disagree on a shared edge") ;
    }
    return result ;
  }


  void CornerBoundary::validateCells(const MappedVolumeSubdivision &subdivision, const std::vector<CornerCell> &cells) const {
    const auto close = [](const Point &a, const Point &b) {
      double error = 0, scale = 1 ;
      for(unsigned d = 0; d < 3; ++d) { error = std::max(error,std::abs(a[d]-b[d])) ; scale = std::max({scale,std::abs(a[d]),std::abs(b[d])}) ; }
      require(error <= 2e-11*scale,"saved shared geometry disagrees with cell geometry") ;
    } ;
    require(cells.size() == subdivision.nodes().size(),"saved cell connectivity count does not match history") ;
    for(TreeNodeId id = 0; id < cells.size(); ++id) {
      const auto &cell = cells[id] ; const auto geometry = subdivision.geometry(id) ; const auto &element = referenceElement(geometry->shape()) ;
      require(cell.vertices.size() == element.vertexCount && cell.faces.size() == element.boundaryCount,"invalid saved cell connectivity") ;
      for(unsigned v = 0; v < element.vertexCount; ++v)
        close(boundaries.points.at(cell.vertices[v]),geometry->evaluate(element.vertices[v]).position) ;
      for(unsigned f = 0; f < element.boundaryCount; ++f) {
        const auto &boundary = element.boundaries[f] ; const auto use = cell.faces[f] ; const auto &face = boundaries.faces.at(use.surface) ;
        require(use.orientation < orientationCount(boundary.shape) && face.vertices.size() == boundary.vertexCount,"invalid saved face orientation") ;
        for(unsigned v = 0; v < boundary.vertexCount; ++v)
          require(cell.vertices[boundary.vertices[v]] == face.vertices[orientedVertex(boundary.shape,use.orientation,v)],"saved cell face has inconsistent vertices") ;
        const auto map = boundaries.geometry.map(face.geometry) ;
        require(map->shape() == boundary.shape,"saved cell face has incorrect shape") ;
        std::vector<Point> probes ; const auto &reference = referenceElement(boundary.shape) ;
        for(unsigned v = 0; v < reference.vertexCount; ++v) probes.push_back(reference.vertices[v]) ;
        probes.push_back(boundary.shape == Shape::triangle ? Point{{1./3,1./3,0}} : Point{{.5,.5,0}}) ;
        for(const auto &point: probes)
          close(geometry->evaluate(boundaryPoint(geometry->shape(),f,point)).position,
                map->evaluate(orientPoint(boundary.shape,use.orientation,point)).position) ;
      }
      cellEdges(cell) ;
    }
  }

  std::vector<TreeNodeId> CornerBoundary::splitCell(MappedVolumeSubdivision &subdivision, std::vector<CornerCell> &cells,
      TreeNodeId leaf, CenterMethod method, const std::vector<std::size_t> &active, const SplitQualityLimits &limits, SplitQuality *report) {
    return constructCell(subdivision,cells,leaf,method,active,limits,report,nullptr) ;
  }
  void CornerBoundary::restoreHexCell(MappedVolumeSubdivision &subdivision, std::vector<CornerCell> &cells,
      HexPlan::Cell leaf, const HexPlan &plan) {
    const auto split = plan.construction(leaf) ;
    constructCell(subdivision,cells,leaf,CenterMethod::reference,{}, {},nullptr,&split) ;
  }

  std::vector<TreeNodeId> CornerBoundary::constructCell(MappedVolumeSubdivision &subdivision, std::vector<CornerCell> &cells,
      TreeNodeId leafId, CenterMethod method, const std::vector<std::size_t> &active,
      const SplitQualityLimits &limits, SplitQuality *report, const VertexVolumeSplit *accepted) {
    if(report) *report = SplitQuality{} ;
    checkMethod(method) ;
    const auto node = subdivision.node(leafId) ;
    require(subdivision.history().active(leafId) && node.isLeaf(),"vertex volume refinement requires an active leaf") ;
    const auto geometry = subdivision.geometry(leafId) ; const auto &element = referenceElement(geometry->shape()) ;
    const CornerCell cell = cells.at(leafId) ; Transaction transaction(*this) ;
    std::vector<Point> centers ; std::vector<std::size_t> choices ;
    for(unsigned f = 0; f < element.boundaryCount; ++f) {
      const auto use = cell.faces[f] ;
      if(accepted) require(boundaries.faces[use.surface].splits.size() == 1,"accepted hex face choice is missing") ;
      const auto choice = accepted ? 0 : splitSurface(use.surface,method,active,transaction) ; choices.push_back(choice) ;
      centers.push_back(orientPoint(element.boundaries[f].shape,inverseOrientation(element.boundaries[f].shape,use.orientation),
        boundaries.faces[use.surface].splits[choice].referenceCenter)) ;
    }
    const auto center = accepted ? accepted->points.back() : cellCenterCoordinates(*geometry,centers,method,options.assessment.measure) ;
    const auto split = std::make_shared<VertexVolumeSplit>(accepted ? *accepted :
      splitVertexVolume(geometry->shape(),centers,center,options.jacobianLimits)) ;
    MappedVolumeSubdivision::AppendGuard historyEdit(subdivision,leafId) ;
    auto ids = subdivision.append(leafId,split) ;
    std::vector<std::shared_ptr<const GeometryMap>> maps ;
    for(auto id: ids) maps.push_back(subdivision.geometry(id)) ;
    if(!accepted) {
      const auto assessment = splitQuality(*geometry,maps,options.assessment) ;
      if(report) *report = assessment ;
      requireSplitQuality(assessment,limits) ;
    }
    if(ids.front() < cells.size()) {
      historyEdit.committed = transaction.committed = true ; return ids ;
    }
    const auto parentEdges = cellEdges(cell) ; std::vector<VertexId> vertices = cell.vertices ;
    for(unsigned e = 0; e < element.edgeCount; ++e)
      vertices.push_back(splitEdge(parentEdges.at(endpoints(cell.vertices[element.edges[e][0]],cell.vertices[element.edges[e][1]])),transaction)) ;
    for(unsigned f = 0; f < cell.faces.size(); ++f) vertices.push_back(boundaries.faces[cell.faces[f].surface].splits[choices[f]].centerVertex) ;
    vertices.push_back(addPoint(geometry->evaluate(center).position)) ;
    std::vector<CornerCell> children(ids.size()) ;
    for(unsigned c = 0; c < children.size(); ++c) {
      for(unsigned v: split->children[c].vertices) children[c].vertices.push_back(vertices[v]) ;
      children[c].faces.resize(6,{absent,0}) ;
    }
    const auto childFaceVertices = [&](const FaceUse &face) {
      std::vector<VertexId> result ; const auto &boundary = referenceElement(Shape::hex).boundaries[face.boundary] ;
      for(unsigned v = 0; v < boundary.vertexCount; ++v) result.push_back(children[face.child].vertices[boundary.vertices[v]]) ;
      return result ;
    } ;
    const auto assign = [&](const FaceUse &face, SurfaceId surface) {
      children[face.child].faces[face.boundary] = {surface,orientation(Shape::quad,childFaceVertices(face),boundaries.faces[surface].vertices)} ;
    } ;
    std::map<std::pair<VertexId,VertexId>,EdgeId> constructionEdges ;
    for(const auto &piece: split->boundaryPieces) {
      const auto use = cell.faces[piece.boundary] ;
      const auto corner = orientedVertex(element.boundaries[piece.boundary].shape,use.orientation,piece.corner) ;
      const auto face = boundaries.faces[use.surface].splits[choices[piece.boundary]].children[corner] ; assign(piece.face,face) ;
      for(auto edge: boundaries.faces[face].edges) {
        const auto &ends = boundaries.edges[edge].vertices ; constructionEdges.emplace(endpoints(ends[0],ends[1]),edge) ;
      }
    }
    // Endpoint lookup is local to this accepted cell construction, never global geometry identity.
    for(const auto &face: split->internalFaces) {
      auto faceVertices = childFaceVertices(face.first) ; const auto map = boundaryMap(maps[face.first.child],face.first.boundary) ;
      std::vector<EdgeId> faceEdges ;
      for(unsigned e = 0; e < faceVertices.size(); ++e) {
        const auto a = faceVertices[e], b = faceVertices[(e+1)%faceVertices.size()] ;
        auto found = constructionEdges.find(endpoints(a,b)) ;
        if(found == constructionEdges.end()) found = constructionEdges.emplace(endpoints(a,b),addEdge(a,b,boundaryMap(map,e))).first ;
        faceEdges.push_back(found->second) ;
      }
      const auto surface = addSurface(std::move(faceVertices),std::move(faceEdges),map) ; assign(face.first,surface) ; assign(face.second,surface) ;
    }
    require(ids.front() == cells.size(),"noncontiguous new vertex volume children") ;
    static_assert(std::is_nothrow_move_constructible<CornerCell>::value,"cell publication requires nonthrowing moves") ;
    const auto needed = cells.size()+children.size() ;
    if(cells.capacity() < needed) {
      const auto maximum = std::min(subdivision.maxNodes(),cells.max_size()) ;
      const auto grown = cells.capacity() > maximum/2 ? maximum : 2*cells.capacity() ;
      cells.reserve(std::max(needed,grown)) ;
    }
    for(auto &child: children) cells.push_back(std::move(child)) ;
    historyEdit.committed = transaction.committed = true ; return ids ;
  }

}
}
