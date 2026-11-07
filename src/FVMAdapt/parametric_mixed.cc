// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "parametric_mixed.h"
#include "library/parametric/face_plan_io.h"
#include <FVMAdapt/parametric/geometry_measure.h>
#include <algorithm>
#include <functional>
#include <set>

namespace Loci { namespace FVMAdapt {
  namespace {
    void require(bool ok,const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    bool segment(const ReferencePoint &p,const ReferencePoint &a,const ReferencePoint &b,Fraction &t) {
      unsigned axis=0 ; while(axis<3 && a[axis]==b[axis]) ++axis ; if(axis==3) return false ;
      t=(p[axis]-a[axis])/(b[axis]-a[axis]) ; if(t<0 || t>1) return false ;
      for(unsigned d=0;d<3;++d) if(p[d]!=a[d]+t*(b[d]-a[d])) return false ;
      return true ;
    }
    // Original-face geometry and node identities, independent of the incident cell type.
    struct SourceFace {
      EntityId id ; Shape shape ; bool reversed ;
      std::vector<EntityId> vertices,edges ; std::vector<Point> corners ;
      SourceFace(EntityId face,const CellTopology &topology,const ParametricCell &cell) : id(face) {
        unsigned f=0 ; while(f<topology.faceCount() && topology.faceId(f)!=face) ++f ;
        require(f<topology.faceCount(),"mixed face is absent from its source cell") ;
        const auto loop=topology.faceVertices(f) ; const unsigned n=loop.size(),orientation=topology.faceOrientation(f) ;
        require(n==3 || n==4,"mixed faces require three or four source vertices") ;
        shape=n==3 ? Shape::triangle:Shape::quad ; reversed=orientation>=n ; vertices.resize(n) ; corners.resize(n) ;
        const auto geometry=cell.native() ? cell.history().geometry().map(cell.history().node(0).map):nullptr ;
        for(unsigned i=0;i<n;++i) {
          const auto j=orientedVertex(n,orientation,i),v=loop[i] ; vertices[j]=topology.vertexId(v) ;
          corners[j]=cell.source() ? cell.source()->corners[v]:geometry ? geometry->evaluate(referenceElement(topology.shape()).vertices[v]).position:
            cell.cornerGeometry().boundaries.points.at(cell.cornerCell(0).vertices[v]) ;
        }
        for(unsigned i=0;i<n;++i) {
          bool found=false ;
          for(unsigned e=0;e<topology.edgeCount();++e) {
            const auto ends=topology.edgeVertices(e) ;
            if(std::minmax(topology.vertexId(ends[0]),topology.vertexId(ends[1]))==std::minmax(vertices[i],vertices[(i+1)%n])) {
              edges.push_back(topology.edgeId(e)) ; found=true ; break ;
            }
          }
          require(found,"mixed source face edge is missing") ;
        }
      }
      ParametricNodeKey key(const ReferencePoint &p) const {
        for(unsigned v=0;v<vertices.size();++v) if(p==referenceVertex(shape,v)) return {{EntityKind::vertex,vertices[v],{}},-1} ;
        for(unsigned e=0;e<edges.size();++e) {
          Fraction t ; if(segment(p,referenceVertex(shape,e),referenceVertex(shape,(e+1)%edges.size()),t)) {
            if(vertices[e]>vertices[(e+1)%vertices.size()]) t=Fraction(1)-t ;
            return {{EntityKind::edge,edges[e],{{t,0,0}}},-1} ;
          }
        }
        return {{EntityKind::face,id,p},-1} ;
      }
      bool point(const ParametricNodeKey &key,ReferencePoint &p) const {
        if(key.local>=0) return false ;
        const auto &v=key.reference ;
        if(v.kind==EntityKind::face && v.entity==id) { p=v.point ; return true ; }
        for(unsigned i=0;i<vertices.size();++i) if(v.kind==EntityKind::vertex && v.entity==vertices[i]) { p=referenceVertex(shape,i) ; return true ; }
        for(unsigned i=0;i<edges.size();++i) if(v.kind==EntityKind::edge && v.entity==edges[i]) {
          auto t=v.point[0] ; if(vertices[i]>vertices[(i+1)%vertices.size()]) t=Fraction(1)-t ;
          const auto a=referenceVertex(shape,i),b=referenceVertex(shape,(i+1)%edges.size()) ;
          for(unsigned d=0;d<3;++d) p[d]=a[d]+t*(b[d]-a[d]) ;
          return true ;
        }
        return false ;
      }
      Point position(const ReferencePoint &p) const { return CornerMap(shape,corners).evaluate(approximate(p)).position ; }
    } ;
    std::vector<BoundaryPolygon> faceCover(const SourceFace &face,const CellTopology &topology,const ParametricCell &cell) {
      std::vector<BoundaryPolygon> cover ; const auto leaves=cell.leafIds() ;
      if(const auto *native=cell.native()) {
        for(unsigned i=0;i<leaves.size();++i) {
          const auto &part=native->node(leaves[i]).cell ; const auto &element=referenceElement(part.shape) ;
          for(unsigned f=0;f<element.boundaryCount;++f) {
            BoundaryPolygon polygon{{i,f},{}} ; const auto &side=element.boundaries[f] ;
            for(unsigned v=0;v<side.vertexCount;++v) {
              const auto key=topology.vertexKey(part.placement.apply(referenceVertex(part.shape,side.vertices[v]))) ; ReferencePoint p ;
              if(!face.point({key,-1},p)) { polygon.vertices.clear() ; break ; } polygon.vertices.push_back(p) ;
            }
            if(!polygon.vertices.empty()) cover.push_back(std::move(polygon)) ;
          }
        }
      } else {
        const auto coordinates=cornerFaceCoordinates(face.id,cell) ;
        std::map<int,unsigned> localFaces ; const auto &b=cell.cornerGeometry().boundaries ; const auto root=cell.cornerFaceId(face.id) ;
        for(unsigned i=0;i<leaves.size();++i) {
          const auto &faces=cell.cornerCell(leaves[i]).faces ;
          for(unsigned f=0;f<faces.size();++f) {
            auto id=faces[f].surface ; while(b.faces[id].parent!=noBoundary) id=b.faces[id].parent ;
            if(id==root) { require(localFaces.emplace(i,f).second,"leaf uses a source face more than once") ; }
          }
        }
        for(const auto &part:cornerFaceFaces(face.id,cell,{})) {
          BoundaryPolygon polygon{{unsigned(part.left),localFaces.at(part.left)},{}} ;
          for(const auto &key:part.vertices) polygon.vertices.push_back(coordinates.at(key)) ;
          cover.push_back(std::move(polygon)) ;
        }
      }
      return cover ;
    }
    struct SharedFace {
      SourceFace face ; std::vector<BoundaryPolygon> left,right ; std::vector<PolygonMatch> pieces ;
      SharedFace(EntityId id,const CellTopology &lt,const ParametricCell &lc,const CellTopology &rt,const ParametricCell &rc)
        : face(id,lt,lc),left(faceCover(face,lt,lc)),right(faceCover(SourceFace(id,rt,rc),rt,rc)) {
        const SourceFace other(id,rt,rc) ;
        require(face.vertices==other.vertices && face.corners==other.corners && face.reversed!=other.reversed,
          "mixed source face geometry or orientation differs") ;
        pieces=reconcilePolygons(face.shape,left,right,1000000) ;
      }
    } ;
    void insertPoints(std::vector<ParametricNodeKey> &loop,const SourceFace &face,const std::vector<ParametricNodeKey> &points) {
      std::vector<ParametricNodeKey> out ;
      for(unsigned i=0;i<loop.size();++i) {
        out.push_back(loop[i]) ; ReferencePoint a,b ;
        if(!face.point(loop[i],a) || !face.point(loop[(i+1)%loop.size()],b)) continue ;
        std::map<Fraction,ParametricNodeKey> inside ;
        for(const auto &key:points) { ReferencePoint p ; Fraction t ;
          if(face.point(key,p) && segment(p,a,b,t) && t>0 && t<1) inside.emplace(t,key) ; }
        for(const auto &entry:inside) out.push_back(entry.second) ;
      }
      loop.swap(out) ;
    }
    FacePatch patch(const SourceFace &face,const std::vector<ReferencePoint> &polygon) {
      FacePatch out ; out.surface={2,static_cast<long long>(face.id)} ; out.orientation=face.reversed ? -1:1 ; out.mappedArea=true ;
      constexpr long long extent=1LL<<30 ; out.bounds={{extent,extent,0,0}} ;
      for(const auto &p:polygon) {
        out.polygon.push_back({{p[0].value(),p[1].value()}}) ;
        for(unsigned d=0;d<2;++d) { out.bounds[d]=std::min(out.bounds[d],static_cast<long long>(std::floor(p[d].value()*extent))) ;
          out.bounds[d+2]=std::max(out.bounds[d+2],static_cast<long long>(std::ceil(p[d].value()*extent))) ; }
      }
      const auto geometry=std::make_shared<CornerMap>(face.shape,face.corners) ;
      for(unsigned i=1;i+1<polygon.size();++i) {
        const auto placement=std::make_shared<CornerMap>(Shape::triangle,std::vector<Point>{approximate(polygon[0]),approximate(polygon[i]),approximate(polygon[i+1])}) ;
        const auto area=geometryMeasure(RestrictedMap(geometry,placement)) ;
        require(area.converged && area.measure>0,"mixed face area integration did not converge") ; out.area+=area.measure ;
      }
      return out ;
    }
    std::map<ParametricNodeKey,ParametricNodeKey> cornerKeys(const ParametricCell &cell,const std::vector<ParametricNodeKey> &required) {
      std::set<EntityId> faces ;
      for(const auto &key:required) if(key.reference.kind==EntityKind::face && key.local<0) faces.insert(key.reference.entity) ;
      std::map<ParametricNodeKey,ParametricNodeKey> keys ;
      for(auto id:faces) {
        const SourceFace face(id,cell.cornerTopology(),cell) ;
        for(const auto &entry:cornerFaceCoordinates(id,cell)) keys.emplace(entry.first,face.key(entry.second)) ;
      }
      return keys ;
    }
  }
  bool mixedFace(const ParametricCell &left,const ParametricCell *right) {
    return right && ((left.native() && right->hasCornerGeometry()) || (left.hasCornerGeometry() && right->native())) ;
  }
  void validateMixedFacePlan(const ParametricRecord &record) {
    if(record.empty()) return ;
    const BoundaryRecordLimits limits ; detail::FacePlanReader in(record,limits) ;
    require(in.vertices==3 || in.vertices==4,"mixed faces require three or four source vertices") ;
    std::function<void(unsigned)> edge=[&](unsigned depth) { if(in.edge(depth)) { edge(depth+1) ; edge(depth+1) ; } } ;
    std::function<void(unsigned,unsigned)> face=[&](unsigned n,unsigned depth) {
      const auto division=in.face(n,depth) ; for(unsigned e=0;e<n;++e) edge(depth) ;
      if(division.split) { require(division.center==detail::midpoint(n),"mixed faces currently require midpoint face centers") ;
        for(unsigned v=0;v<n;++v) face(4,depth+1) ; }
    } ;
    face(in.vertices,0) ; in.input.finish() ;
  }
  std::vector<ParametricOutputFace> mixedFaceFaces(EntityId id,const CellTopology &lt,const ParametricCell &lc,
      const CellTopology &rt,const ParametricCell &rc,const std::vector<ParametricNodeKey> &required,bool patches) {
    const SharedFace shared(id,lt,lc,rt,rc) ; std::vector<ParametricOutputFace> result ; auto points=required ;
    for(const auto &piece:shared.pieces) {
      ParametricOutputFace f{{},int(piece.first.child),int(piece.second.child)} ;
      for(const auto &p:piece.vertices) f.vertices.push_back(shared.face.key(p)) ;
      if(shared.face.reversed) std::reverse(f.vertices.begin(),f.vertices.end()) ;
      points.insert(points.end(),f.vertices.begin(),f.vertices.end()) ;
      if(patches) f.patch=patch(shared.face,piece.vertices) ;
      result.push_back(std::move(f)) ;
    }
    std::sort(points.begin(),points.end()) ; points.erase(std::unique(points.begin(),points.end()),points.end()) ;
    for(auto &f:result) insertPoints(f.vertices,shared.face,points) ;
    std::sort(result.begin(),result.end(),[](const ParametricOutputFace &a,const ParametricOutputFace &b) { return a.vertices<b.vertices ; }) ;
    return result ;
  }
  std::vector<ParametricNodeKey> mixedFaceNodes(EntityId id,const CellTopology &lt,const ParametricCell &lc,
      const CellTopology &rt,const ParametricCell &rc,std::vector<Point> *positions) {
    std::set<ParametricNodeKey> nodes ;
    for(const auto &face:mixedFaceFaces(id,lt,lc,rt,rc,{})) for(const auto &key:face.vertices)
      if(key.reference.kind==EntityKind::face) nodes.insert(key) ;
    if(positions) {
      const SourceFace face(id,lt,lc) ; positions->clear() ; positions->reserve(nodes.size()) ;
      for(const auto &key:nodes) positions->push_back(face.position(key.reference.point)) ;
    }
    return {nodes.begin(),nodes.end()} ;
  }
  void mixedCornerFaces(const ParametricCell &cell,const std::vector<ParametricNodeKey> &required,std::vector<ParametricOutputFace> &faces) {
    const auto keys=cornerKeys(cell,required) ; if(keys.empty()) return ; std::set<EntityId> originals ;
    for(const auto &entry:keys) if(entry.second.reference.kind==EntityKind::face) originals.insert(entry.second.reference.entity) ;
    for(auto &face:faces) for(auto &key:face.vertices) { const auto found=keys.find(key) ; if(found!=keys.end()) key=found->second ; }
    for(auto id:originals) { const SourceFace source(id,cell.cornerTopology(),cell) ; for(auto &face:faces) insertPoints(face.vertices,source,required) ; }
  }
  void mixedCornerVertices(const ParametricCell &cell,const std::vector<ParametricNodeKey> &required,std::vector<std::vector<ParametricNodeKey>> &vertices) {
    const auto keys=cornerKeys(cell,required) ;
    for(auto &leaf:vertices) for(auto &key:leaf) { const auto found=keys.find(key) ; if(found!=keys.end()) key=found->second ; }
  }
  std::vector<FaceBalanceRequest> mixedDividedFaces(EntityId id,const CellTopology &lt,const ParametricCell &lc,const CellTopology &rt,const ParametricCell &rc) {
    const SharedFace shared(id,lt,lc,rt,rc) ; using Index=std::pair<unsigned,unsigned> ; std::map<Index,unsigned> left,right ;
    for(const auto &p:shared.pieces) { ++left[{p.first.child,p.first.boundary}] ; ++right[{p.second.child,p.second.boundary}] ; }
    std::vector<FaceBalanceRequest> result ;
    for(const auto &entry:left) if(entry.second>1) result.push_back({int(lt.id()),int(entry.first.first),int(entry.first.second)}) ;
    for(const auto &entry:right) if(entry.second>1) result.push_back({int(rt.id()),int(entry.first.first),int(entry.first.second)}) ;
    return result ;
  }
  ParametricTrace mixedFaceRequests(EntityId id,const CellTopology &lt,const ParametricCell &lc,const CellTopology &rt,const ParametricCell &rc,
      unsigned difference,NeighborGrading grading) {
    const auto faces=mixedFaceFaces(id,lt,lc,rt,rc,{}) ; const auto ll=lc.leafIds(),rl=rc.leafIds() ;
    std::vector<int> depths ; for(auto leaf:ll) depths.push_back(lc.path(leaf).size()) ; for(auto leaf:rl) depths.push_back(rc.path(leaf).size()) ;
    const int offset=ll.size() ; std::set<int> required,covered ; std::map<std::vector<ParametricNodeKey>,std::set<int>> contacts ;
    const auto grade=[&](const std::set<int> &cells) {
      int maximum=0 ; for(int cell:cells) maximum=std::max(maximum,depths.at(cell)) ;
      for(int cell:cells) if(unsigned(maximum-depths.at(cell))>difference) required.insert(cell) ;
    } ;
    for(const auto &face:faces) {
      grade({face.left,offset+face.right}) ; covered.insert(face.left) ; covered.insert(offset+face.right) ;
      if(grading==NeighborGrading::faces) continue ;
      for(unsigned i=0;i<face.vertices.size();++i) {
        std::vector<ParametricNodeKey> key{face.vertices[i]} ;
        if(grading==NeighborGrading::edges) { key.push_back(face.vertices[(i+1)%face.vertices.size()]) ; std::sort(key.begin(),key.end()) ; }
        contacts[key].insert(face.left) ; contacts[key].insert(offset+face.right) ;
      }
    }
    // Some prism children touch a source side only along an edge.
    if(grading!=NeighborGrading::faces) {
      const auto add=[&](const CellTopology &topology,const ParametricCell &cell,int shift) {
        const SourceFace source(id,topology,cell) ; const auto own=parametricLeafVertices(topology,cell) ;
        const auto mapped=cell.native() ? std::map<ParametricNodeKey,ReferencePoint>{}:cornerFaceCoordinates(id,cell) ;
        for(unsigned i=0;i<own.size();++i) if(!covered.count(shift+i)) {
          std::vector<ReferencePoint> points ;
          for(const auto &key:own[i]) { ReferencePoint p ; const auto found=mapped.find(key) ;
            if(found!=mapped.end()) points.push_back(found->second) ; else if(source.point(key,p)) points.push_back(p) ; }
          if(points.empty()) continue ;
          std::sort(points.begin(),points.end()) ;
          for(auto &entry:contacts) {
            bool match=true ; for(const auto &key:entry.first) { ReferencePoint p ; Fraction t ;
              match=match && source.point(key,p) && (points.front()==points.back() ? p==points.front():segment(p,points.front(),points.back(),t)) ; }
            if(match) entry.second.insert(shift+i) ;
          }
        }
      } ;
      add(lt,lc,0) ; add(rt,rc,offset) ; for(const auto &entry:contacts) grade(entry.second) ;
    }
    ParametricTrace result ; for(int i:required) result.push_back({int(i<offset ? lt.id():rt.id()),i<offset ? i:i-offset}) ; return result ;
  }
} }
