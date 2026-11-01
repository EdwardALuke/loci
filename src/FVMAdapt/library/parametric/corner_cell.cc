// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/corner_boundary.h>
#include <set>

namespace Loci {
namespace FVMAdapt {
  namespace {
    using Id=BoundaryId ;
    void require(bool ok,const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    std::vector<Id> outwardVertices(const BoundaryRecords &b,CornerFaceUse use) {
      const auto &face=b.faces.at(use.surface) ;
      require(use.orientation<2*face.vertices.size(),"invalid cell face orientation") ;
      std::vector<Id> result ;
      for(unsigned v=0; v<face.vertices.size(); ++v) result.push_back(face.vertices[orientedVertex(unsigned(face.vertices.size()),use.orientation,v)]) ;
      return result ;
    }
    Point average(const BoundaryRecords &b,const std::vector<Id> &vertices) {
      require(!vertices.empty(),"cell requires vertices") ; const auto origin=b.points.at(vertices[0]) ;
      std::array<long double,3> sum{} ;
      for(auto v:vertices) for(unsigned d=0; d<3; ++d) sum[d]+=static_cast<long double>(b.points.at(v)[d])-origin[d] ;
      Point result ; for(unsigned d=0; d<3; ++d) result[d]=origin[d]+sum[d]/vertices.size() ; return result ;
    }
    std::shared_ptr<const GeometryMap> directedEdge(const BoundaryRecords &b,Id edge,Id start) {
      const auto &record=b.edges.at(edge) ; auto map=b.geometry.map(record.geometry) ;
      require(start==record.vertices[0] || start==record.vertices[1],"edge does not contain its start vertex") ;
      if(record.vertices[0]!=start) map=std::make_shared<RestrictedMap>(map,
        std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{{{1,0,0}},{{0,0,0}}})) ;
      return map ;
    }
  }

  void CornerBoundary::validateCell(const CornerCell &cell) const {
    require(cell.vertices.size()>=4 && cell.faces.size()>=4,"cell boundary is incomplete") ;
    const std::set<Id> vertices(cell.vertices.begin(),cell.vertices.end()) ;
    require(vertices.size()==cell.vertices.size(),"cell repeats a vertex") ;
    for(auto v:vertices) for(double x:boundaries.points.at(v)) require(std::isfinite(x),"cell has a nonfinite vertex") ;
    std::set<Id> surfaces,used ;
    std::map<std::pair<Id,Id>,std::vector<std::pair<unsigned,bool>>> incidence ;
    std::map<Id,std::vector<std::pair<Id,Id>>> links ;
    for(unsigned f=0; f<cell.faces.size(); ++f) {
      const auto use=cell.faces[f] ; const auto &face=boundaries.faces.at(use.surface) ;
      require(surfaces.insert(use.surface).second,"cell repeats a face") ;
      const auto map=boundaries.geometry.map(face.geometry) ;
      require(map->dimension()==2 && face.vertices.size()==map->referenceCorners().size(),"cell boundary requires polygon surfaces") ;
      require(face.edges.size()==face.vertices.size(),"cell face has incorrect edge count") ;
      const auto loop=outwardVertices(boundaries,use) ;
      require(std::set<Id>(loop.begin(),loop.end()).size()==loop.size(),"cell face repeats a vertex") ;
      for(unsigned i=0; i<loop.size(); ++i) {
        const auto a=loop[i],b=loop[(i+1)%loop.size()] ; used.insert(a) ;
        incidence[endpoints(a,b)].push_back({f,a<b}) ;
        links[a].push_back({loop[(i+loop.size()-1)%loop.size()],b}) ;
        const auto &edge=boundaries.edges.at(face.edges[i]) ;
        require(endpoints(edge.vertices[0],edge.vertices[1])==endpoints(face.vertices[i],face.vertices[(i+1)%loop.size()]),
                "cell face edges do not follow its vertices") ;
      }
    }
    require(used==vertices,"cell vertex list differs from its faces") ; cellEdges(cell) ;
    std::vector<std::vector<unsigned>> adjacent(cell.faces.size()) ;
    for(const auto &item:incidence) {
      const auto &uses=item.second ;
      require(uses.size()==2 && uses[0].second!=uses[1].second,"cell edge requires two opposite face uses") ;
      adjacent[uses[0].first].push_back(uses[1].first) ; adjacent[uses[1].first].push_back(uses[0].first) ;
    }
    std::set<unsigned> reached{0} ; std::vector<unsigned> pending{0} ;
    for(std::size_t i=0; i<pending.size(); ++i) for(auto next:adjacent[pending[i]])
      if(reached.insert(next).second) pending.push_back(next) ;
    require(reached.size()==cell.faces.size() && vertices.size()+cell.faces.size()==incidence.size()+2,
            "cell boundary must be a connected sphere") ;
    for(const auto &entry:links) {
      std::map<Id,Id> next ;
      for(auto arc:entry.second) require(next.emplace(arc.first,arc.second).second,"cell vertex has a branched face fan") ;
      std::set<Id> visited ; auto v=next.begin()->first ;
      while(visited.insert(v).second) {
        const auto found=next.find(v) ; require(found!=next.end(),"cell vertex has an open face fan") ; v=found->second ;
      }
      require(v==next.begin()->first && visited.size()==next.size(),"cell vertex has disconnected face fans") ;
    }
  }

  std::vector<std::shared_ptr<const GeometryMap>> CornerBoundary::cellPieces(const CornerCell &cell,const Point &interior) const {
    validateCell(cell) ; std::vector<std::shared_ptr<const GeometryMap>> result ;
    for(auto use:cell.faces) {
      const auto &face=boundaries.faces[use.surface] ; const auto map=boundaries.geometry.map(face.geometry) ;
      const auto reference=map->referenceCorners() ; const unsigned n=reference.size() ; std::vector<Point> coordinates ;
      // Reverse the outward face so its normal points toward the cone apex.
      for(unsigned v=0;v<n;++v) coordinates.push_back(reference[orientedVertex(n,use.orientation,(n-v)%n)]) ;
      if(map->shape()==Shape::polygon) {
        for(unsigned v=0;v<n;++v) {
          auto base=std::make_shared<RestrictedMap>(map,std::make_shared<CornerMap>(Shape::triangle,
            std::vector<Point>{{{.5,.5,0}},coordinates[v],coordinates[(v+1)%n]})) ;
          result.push_back(std::make_shared<ConeMap>(base,interior)) ;
        }
      } else {
        auto base=std::make_shared<RestrictedMap>(map,std::make_shared<CornerMap>(map->shape(),coordinates)) ;
        result.push_back(std::make_shared<ConeMap>(base,interior)) ;
      }
    }
    return result ;
  }

  GeometryMeasure CornerBoundary::measureCell(const CornerCell &cell,const Point &interior) const {
    GeometryMeasure result ; result.converged=true ; result.estimatedMeasureError=result.estimatedCentroidError=0 ;
    long double volume=0, centroidError=0 ; std::array<long double,3> moment{} ;
    for(const auto &piece:cellPieces(cell,interior)) {
      auto samples=options.assessment.quality ;
      if(samples.samples.empty()) {
        // Cone coordinates have no unique derivative at the apex.
        for(double z:{0.,.5,.9375}) for(unsigned j=0; j<=2; ++j) for(unsigned i=0; i<=2; ++i)
          if(piece->shape()==Shape::pyramid || i+j<=2) samples.samples.push_back({{i*.5*(1-z),j*.5*(1-z),z}}) ;
      }
      require(geometryQuality(*piece,samples).status==GeometryQualityStatus::sampled,"cell cone has invalid or unresolved sampled geometry") ;
      const auto value=geometryMeasure(*piece,options.assessment.measure) ;
      if(!value.converged) throw std::length_error("cell cone integration did not converge") ;
      volume+=value.measure ;
      for(unsigned d=0; d<3; ++d) moment[d]+=value.measure*(static_cast<long double>(value.centroid[d])-interior[d]) ;
      result.order=std::max(result.order,value.order) ; result.evaluations+=value.evaluations ;
      result.estimatedMeasureError+=value.estimatedMeasureError ;
      centroidError+=value.measure*value.estimatedCentroidError+value.estimatedMeasureError*
        std::hypot(value.centroid[0]-interior[0],value.centroid[1]-interior[1],value.centroid[2]-interior[2]) ;
    }
    result.measure=volume ;
    require(std::isfinite(result.measure) && result.measure>0,"cell volume is not positive and finite") ;
    for(unsigned d=0; d<3; ++d) {
      result.centroid[d]=interior[d]+moment[d]/volume ; require(std::isfinite(result.centroid[d]),"cell centroid is not finite") ;
    }
    result.estimatedCentroidError=(centroidError+result.estimatedMeasureError*
      std::hypot(moment[0]/volume,moment[1]/volume,moment[2]/volume))/volume ;
    return result ;
  }

  CornerSplit CornerBoundary::splitCell(const CornerCell &cell,Point interior,CenterMethod method,
                                        const std::vector<std::size_t> &active) {
    Transaction transaction(*this) ;
    auto result=constructCell(cell,interior,method,active,transaction) ;
    transaction.committed=true ; return result ;
  }

  CornerSplit CornerBoundary::constructCell(const CornerCell &cell,Point interior,CenterMethod method,
                                            const std::vector<std::size_t> &active,Transaction &transaction) {
    require(method==CenterMethod::reference || method==CenterMethod::weighted || method==CenterMethod::integrated,
            "unknown cell center method") ;
    require(cell.vertices.size()<=options.assessment.maxChildren,"corner child count exceeds its limit") ;
    const auto parentMeasure=measureCell(cell,interior) ; CornerSplit result ; result.centerMethod=method ;
    std::vector<Id> centers ;
    for(auto use:cell.faces) {
      const auto choice=splitSurface(use.surface,method,active,transaction) ; result.faceSplits.push_back(choice) ;
      centers.push_back(boundaries.faces[use.surface].splits[choice].centerVertex) ;
    }
    Point center=interior ;
    if(method==CenterMethod::integrated) center=parentMeasure.centroid ;
    if(method==CenterMethod::weighted) {
      std::array<long double,3> sum{} ; long double total=0 ;
      for(unsigned f=0; f<cell.faces.size(); ++f) {
        const auto loop=outwardVertices(boundaries,cell.faces[f]) ; const auto origin=boundaries.points[loop[0]] ;
        std::array<long double,3> area{} ;
        for(unsigned i=0; i<loop.size(); ++i) {
          std::array<long double,3> a,b ;
          for(unsigned d=0; d<3; ++d) { a[d]=static_cast<long double>(boundaries.points[loop[i]][d])-origin[d] ;
            b[d]=static_cast<long double>(boundaries.points[loop[(i+1)%loop.size()]][d])-origin[d] ; }
          for(unsigned d=0; d<3; ++d) area[d]+=a[(d+1)%3]*b[(d+2)%3]-a[(d+2)%3]*b[(d+1)%3] ;
        }
        const auto weight=std::hypot(area[0],area[1],area[2]) ; total+=weight ;
        for(unsigned d=0; d<3; ++d) sum[d]+=weight*(static_cast<long double>(boundaries.points[centers[f]][d])-interior[d]) ;
      }
      require(total>0 && std::isfinite(total),"cell center requires finite nonzero face areas") ;
      for(unsigned d=0; d<3; ++d) center[d]=interior[d]+sum[d]/total ;
    }
    measureCell(cell,center) ; result.centerVertex=addPoint(center) ; result.children.resize(cell.vertices.size()) ;
    std::map<Id,unsigned> childIds ; for(unsigned v=0; v<cell.vertices.size(); ++v) childIds.emplace(cell.vertices[v],v) ;
    std::map<std::pair<Id,Id>,Id> edges ;
    const auto remember=[&](Id face) {
      for(auto edge:boundaries.faces[face].edges) { const auto ends=boundaries.edges[edge].vertices ; edges.emplace(endpoints(ends[0],ends[1]),edge) ; }
    } ;
    // Each parent face contributes one outward patch to each incident corner child.
    for(unsigned f=0; f<cell.faces.size(); ++f) {
      const auto use=cell.faces[f] ; const auto parent=boundaries.faces[use.surface] ; const bool reversed=use.orientation>=parent.vertices.size() ;
      const auto children=parent.splits[result.faceSplits[f]].children ;
      for(unsigned v=0; v<parent.vertices.size(); ++v) {
        result.children[childIds.at(parent.vertices[v])].faces.push_back({children[v],reversed ? 4u:0u}) ; remember(children[v]) ;
      }
    }
    for(auto faceCenter:centers) edges.emplace(endpoints(faceCenter,result.centerVertex),addEdge(faceCenter,result.centerVertex,
      std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{boundaries.points[faceCenter],center}))) ;
    // An original edge separates two corner children and meets two parent faces.
    const auto parentEdges=cellEdges(cell) ;
    for(const auto &entry:parentEdges) {
      const auto edge=entry.second ; const auto ends=boundaries.edges[edge].vertices ; std::vector<unsigned> incident ;
      for(unsigned f=0; f<cell.faces.size(); ++f) {
        const auto &list=boundaries.faces[cell.faces[f].surface].edges ;
        if(std::find(list.begin(),list.end(),edge)!=list.end()) incident.push_back(f) ;
      }
      require(incident.size()==2,"cell edge requires two faces") ;
      const auto midpoint=splitEdge(edge,transaction) ;
      std::vector<Id> vertices{midpoint,centers[incident[0]],result.centerVertex,centers[incident[1]]}, faceEdges ;
      std::array<std::shared_ptr<const GeometryMap>,4> curves ;
      for(unsigned e=0; e<4; ++e) {
        const auto id=edges.at(endpoints(vertices[e],vertices[(e+1)%4])) ; faceEdges.push_back(id) ;
        curves[e]=directedEdge(boundaries,id,vertices[e]) ;
      }
      const auto surface=addSurface(vertices,faceEdges,std::make_shared<BoundarySurfaceMap>(std::move(curves))) ;
      // Along [midpoint, first face center], the new face opposes the exterior patch.
      const auto loop=outwardVertices(boundaries,cell.faces[incident[0]]) ; bool forward=false ;
      for(unsigned v=0; v<loop.size(); ++v) if(loop[v]==ends[0] && loop[(v+1)%loop.size()]==ends[1]) forward=true ;
      result.children[childIds.at(ends[0])].faces.push_back({surface,forward ? 4u:0u}) ;
      result.children[childIds.at(ends[1])].faces.push_back({surface,forward ? 0u:4u}) ;
    }
    for(auto &child:result.children) {
      std::set<Id> vertices ;
      for(auto use:child.faces) for(auto v:boundaries.faces[use.surface].vertices) vertices.insert(v) ;
      child.vertices.assign(vertices.begin(),vertices.end()) ; const auto seed=average(boundaries,child.vertices) ;
      measureCell(child,seed) ; result.interiors.push_back(seed) ;
    }
    return result ;
  }
}
}
