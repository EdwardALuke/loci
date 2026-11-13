// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/surface_map.h>
#include "exact_geometry.h"
#include "geometry_validity.h"
#include "geometry_integration.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {
    using Map=std::shared_ptr<const GeometryMap> ;
    using detail::Exact ;
    using detail::ExactPoint ;
    using XY=std::array<Exact,2> ;
    using Polygon=std::vector<XY> ;
    void require(bool ok,const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    void finite(const Point &p) { for(double x:p) require(std::isfinite(x),"surface coordinates must be finite") ; }
    Point referenceCenter(Shape s) { return s==Shape::triangle ? Point{{1./3,1./3,0}}:Point{{.5,.5,0}} ; }
    std::vector<Point> reference(Shape s,unsigned n) {
      if(s==Shape::polygon) return polygonReferenceCorners(n) ;
      const auto &e=referenceElement(s) ; return {e.vertices.begin(),e.vertices.begin()+e.vertexCount} ;
    }
    Shape surfaceShape(unsigned n) { require(n>=3,"surface requires at least three corners") ; return n==3 ? Shape::triangle:n==4 ? Shape::quad:Shape::polygon ; }
    Map coordinates(const std::vector<Point> &p) {
      if(p.size()>4) return std::make_shared<PolygonMap>(p) ;
      return std::make_shared<CornerMap>(surfaceShape(p.size()),p) ;
    }
    Exact turn(const XY &a,const XY &b,const XY &c) { return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]) ; }
    Polygon exactPolygon(const std::vector<Point> &p) {
      Polygon out ; for(const auto &q:p) { finite(q) ; require(q[2]==0,"surface chart must lie in the XY plane") ; out.push_back({{Exact(q[0]),Exact(q[1])}}) ; } return out ;
    }
    Exact twiceSignedArea(const Polygon &p) {
      Exact a=0 ;
      for(unsigned i=0;i<p.size();++i)
        a+=p[i][0]*p[(i+1)%p.size()][1]-p[i][1]*p[(i+1)%p.size()][0] ;
      return a ;
    }
    void requireConvex(const Polygon &p) {
      const auto a=twiceSignedArea(p) ;
      require(a!=0,"surface chart has zero area") ;
      for(unsigned i=0;i<p.size();++i) {
        require(p[i]!=p[(i+1)%p.size()],"surface chart has a zero edge") ;
        for(const auto &q:p)
          require(turn(p[i],p[(i+1)%p.size()],q)*a>=0,"surface chart must be convex") ;
      }
    }
    void validateFan(const std::vector<Point> &corners,const Point &middle) {
      using EPoint=std::array<Exact,3> ;
      std::vector<EPoint> p ;
      for(const auto &v:corners)
        p.push_back({{Exact(v[0])-Exact(middle[0]),Exact(v[1])-Exact(middle[1]),Exact(v[2])-Exact(middle[2])}}) ;
      const auto cross=[](const EPoint &a,const EPoint &b) -> EPoint {
        return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
      } ;
      const auto dot=[](const EPoint &a,const EPoint &b) -> Exact {
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] ;
      } ;

      // Project along the summed face normal; the center must see every projected edge.
      EPoint normal{} ;
      for(unsigned i=0;i<p.size();++i) {
        const auto n=cross(p[i],p[(i+1)%p.size()]) ;
        for(unsigned d=0;d<3;++d) normal[d]+=n[d] ;
      }
      require(dot(normal,normal)>0,"face fan has zero area") ;
      EPoint axis{{1,0,0}} ;
      if(normal[1]==0 && normal[2]==0) axis={{0,1,0}} ;
      const auto u=cross(normal,axis),v=cross(normal,u) ;
      Polygon projected ;
      for(const auto &q:p) projected.push_back({{dot(u,q),dot(v,q)}}) ;
      const XY zero{{0,0}} ;
      for(unsigned i=0;i<p.size();++i)
        require(turn(zero,projected[i],projected[(i+1)%p.size()])>0,"face fan center must see each projected edge") ;

      // Nonadjacent boundary edges must not cross or touch in the projection.
      const auto withinSegmentBounds=[](const XY &a,const XY &b,const XY &q) {
        return q[0]>=std::min(a[0],b[0]) && q[0]<=std::max(a[0],b[0]) &&
               q[1]>=std::min(a[1],b[1]) && q[1]<=std::max(a[1],b[1]) ;
      } ;
      for(unsigned i=0;i<p.size();++i)
        for(unsigned j=i+1;j<p.size();++j) {
          if(j==i+1 || (i==0 && j+1==p.size())) continue ;
          const auto &a=projected[i],&b=projected[(i+1)%p.size()] ;
          const auto &c=projected[j],&d=projected[(j+1)%p.size()] ;
          const Exact x=turn(a,b,c),y=turn(a,b,d),z=turn(c,d,a),w=turn(c,d,b) ;
          require(!((x*y<0 && z*w<0) ||
                    (x==0 && withinSegmentBounds(a,b,c)) || (y==0 && withinSegmentBounds(a,b,d)) ||
                    (z==0 && withinSegmentBounds(c,d,a)) || (w==0 && withinSegmentBounds(c,d,b))),
                  "face fan boundary intersects itself") ;
        }
    }
    Polygon clip(Polygon p,const Polygon &boundary) {
      const auto sign=twiceSignedArea(boundary)>0 ? 1:-1 ;
      for(unsigned e=0;e<boundary.size() && !p.empty();++e) {
        Polygon next ; const auto &a=boundary[e],&b=boundary[(e+1)%boundary.size()] ;
        for(unsigned i=0;i<p.size();++i) {
          const auto &u=p[i],&v=p[(i+1)%p.size()] ; const Exact x=sign*turn(a,b,u),y=sign*turn(a,b,v) ;
          if(x>=0) next.push_back(u) ;
          if((x<0 && y>0) || (x>0 && y<0)) {
            const Exact t=x/(x-y) ; next.push_back({{u[0]+t*(v[0]-u[0]),u[1]+t*(v[1]-u[1])}}) ;
          }
        }
        p.clear() ; for(const auto &q:next) if(p.empty() || p.back()!=q) p.push_back(q) ;
        if(p.size()>1 && p.front()==p.back()) p.pop_back() ;
      }
      return p ;
    }
    MapValue compose(const MapValue &a,const MapValue &b) {
      MapValue out{a.position,{},a.hasJacobian && b.hasJacobian} ;
      if(out.hasJacobian) for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k)
        out.jacobian[i][j]+=a.jacobian[k][j]*b.jacobian[i][k] ;
      return out ;
    }
    MapValue affineTriangle(const std::array<Point,3> &q,const std::array<Point,3> &p,const Point &x) {
      const long double a=static_cast<long double>(q[1][0])-q[0][0],b=static_cast<long double>(q[2][0])-q[0][0],
        c=static_cast<long double>(q[1][1])-q[0][1],d=static_cast<long double>(q[2][1])-q[0][1],det=a*d-b*c ;
      require(det!=0,"surface triangle has zero chart area") ;
      const long double dx=static_cast<long double>(x[0])-q[0][0],dy=static_cast<long double>(x[1])-q[0][1],u=(d*dx-b*dy)/det,v=(a*dy-c*dx)/det ;
      MapValue out{{},{},true} ;
      for(unsigned k=0;k<3;++k) {
        const long double e=static_cast<long double>(p[1][k])-p[0][k],f=static_cast<long double>(p[2][k])-p[0][k] ;
        out.position[k]=p[0][k]+u*e+v*f ; out.jacobian[0][k]=(d*e-c*f)/det ; out.jacobian[1][k]=(a*f-b*e)/det ;
      }
      return out ;
    }
    // Expose one smooth evaluator without recursively partitioning its known restriction.
    class SmoothMap final : public GeometryMap {
      Map map_ ;
    public:
      explicit SmoothMap(Map map):map_(std::move(map)) {}
      const GeometryMap &map() const { return *map_ ; }
      Shape shape() const override { return map_->shape() ; }
      ReferenceCorners referenceCorners() const override { return map_->referenceCorners() ; }
      bool containsReference(const Point &q,double t=0) const override { return map_->containsReference(q,t) ; }
      MapValue evaluate(const Point &q) const override { return map_->evaluate(q) ; }
    } ;
    // Differences are taken before rounding, so small clipped pieces retain their derivatives.
    class AffineTriangleMap final : public GeometryMap {
      std::array<ExactPoint,3> exact_ ;
      std::array<long double,3> origin_ ;
      std::array<std::array<long double,3>,2> derivative_ ;
    public:
      explicit AffineTriangleMap(const std::array<ExactPoint,3> &p) {
        for(unsigned d=0;d<3;++d) {
          exact_[0][d]=p[0][d] ;
          origin_[d]=p[0][d].convert_to<long double>() ;
          for(unsigned i=0;i<2;++i) {
            exact_[i+1][d]=p[i+1][d]-p[0][d] ; derivative_[i][d]=exact_[i+1][d].convert_to<long double>() ;
          }
        }
      }
      detail::ExactMapView algebra() const {
        detail::ExactMapView view ; view.kind=detail::ExactMapView::Kind::affine ; view.origin=exact_[0] ;
        view.columns[0]=exact_[1] ; view.columns[1]=exact_[2] ; return view ;
      }
      Shape shape() const override { return Shape::triangle ; }
      MapValue evaluate(const Point &q) const override {
        finite(q) ; MapValue out{{},{},true} ;
        for(unsigned d=0;d<3;++d) {
          out.position[d]=origin_[d]+q[0]*derivative_[0][d]+q[1]*derivative_[1][d] ;
          for(unsigned i=0;i<2;++i) out.jacobian[i][d]=derivative_[i][d] ;
        }
        return out ;
      }
    } ;
    ExactPoint fanPosition(const FanMap &fan,unsigned sector,const XY &x) {
      const auto refs=fan.referenceCorners() ;
      const auto q=exactPolygon({referenceCenter(fan.shape()),refs[sector],refs[(sector+1)%refs.size()]}) ;
      const Exact a=q[1][0]-q[0][0],b=q[2][0]-q[0][0],c=q[1][1]-q[0][1],d=q[2][1]-q[0][1],det=a*d-b*c ;
      const Exact dx=x[0]-q[0][0],dy=x[1]-q[0][1],u=(d*dx-b*dy)/det,v=(a*dy-c*dx)/det ;
      ExactPoint p ;
      for(unsigned k=0;k<3;++k) p[k]=Exact(fan.center()[k])+u*(Exact(fan.corners()[sector][k])-Exact(fan.center()[k]))+
        v*(Exact(fan.corners()[(sector+1)%refs.size()][k])-Exact(fan.center()[k])) ;
      return p ;
    }
    ExactPoint exactChartPosition(const GeometryMap &map,const ExactPoint &q) {
      if(const auto *r=dynamic_cast<const RestrictedMap *>(&map)) return exactChartPosition(*r->parent(),exactChartPosition(*r->referenceMap(),q)) ;
      const auto *c=dynamic_cast<const CornerMap *>(&map) ;
      if(c && (c->shape()==Shape::segment || c->shape()==Shape::triangle || c->shape()==Shape::quad)) {
        const auto &p=c->corners() ; ExactPoint out ;
        for(unsigned d=0;d<3;++d) {
          out[d]=Exact(p[0][d])+q[0]*(Exact(p[1][d])-Exact(p[0][d])) ;
          if(c->shape()==Shape::triangle) out[d]+=q[1]*(Exact(p[2][d])-Exact(p[0][d])) ;
          if(c->shape()==Shape::quad) out[d]+=q[1]*(Exact(p[3][d])-Exact(p[0][d]))+
            q[0]*q[1]*(Exact(p[2][d])-Exact(p[1][d])-Exact(p[3][d])+Exact(p[0][d])) ;
        }
        return out ;
      }
      if(const auto *polygon=dynamic_cast<const PolygonMap *>(&map)) {
        const auto refs=polygon->referenceCorners() ; const auto r=exactPolygon({refs.begin(),refs.end()}) ; const XY x{{q[0],q[1]}} ;
        for(unsigned i=0;i<r.size();++i) {
          const auto &a=r[i],&b=r[(i+1)%r.size()] ;
          if(turn(a,b,x)!=0) continue ;
          const unsigned axis=a[0]!=b[0] ? 0:1 ; const Exact t=(x[axis]-a[axis])/(b[axis]-a[axis]) ;
          ExactPoint out ; for(unsigned d=0;d<3;++d) out[d]=Exact(polygon->corners()[i][d])+t*
            (Exact(polygon->corners()[(i+1)%r.size()][d])-Exact(polygon->corners()[i][d])) ; return out ;
        }
        Exact total=0 ; ExactPoint out{} ;
        for(unsigned i=0;i<r.size();++i) {
          const auto &a=r[(i+r.size()-1)%r.size()],&b=r[i],&c=r[(i+1)%r.size()] ;
          const Exact weight=turn(a,b,c)/(turn(a,b,x)*turn(b,c,x)) ; total+=weight ;
          for(unsigned d=0;d<3;++d) out[d]+=weight*Exact(polygon->corners()[i][d]) ;
        }
        require(total!=0,"surface chart has singular coordinates") ; for(auto &x:out) x/=total ; return out ;
      }
      throw std::invalid_argument("piecewise geometry requires a supported chart placement") ;
    }
    unsigned fanSector(ReferenceCorners refs,Shape shape,const ExactPoint &point) {
      const auto c=referenceCenter(shape) ; const XY x{{point[0],point[1]}} ;
      unsigned selected=0 ; Exact best ;
      for(unsigned i=0;i<refs.size();++i) {
        const auto triangle=exactPolygon({c,refs[i],refs[(i+1)%refs.size()]}) ;
        bool inside=true ; Exact score ;
        for(unsigned e=0;e<3;++e) {
          const Exact value=turn(triangle[e],triangle[(e+1)%3],x) ;
          inside=inside && value>=0 ; if(e==0 || value<score) score=value ;
        }
        if(inside) return i ;
        if(i==0 || score>best) { best=score ; selected=i ; }
      }
      // Accepted chart coordinates can round just beyond an outer edge. Continue
      // the same nearest sector as evaluation, without moving those coordinates.
      const Exact tolerance=64*std::numeric_limits<double>::epsilon() ;
      require(point[2]*point[2]<=tolerance*tolerance,"surface piece lies outside its retained fan") ;
      for(unsigned e=0;e<refs.size();++e) {
        const auto edge=exactPolygon({refs[e],refs[(e+1)%refs.size()]}) ;
        const Exact distance=turn(edge[0],edge[1],x),dx=edge[1][0]-edge[0][0],dy=edge[1][1]-edge[0][1] ;
        require(distance>=0 || distance*distance<=tolerance*tolerance*(dx*dx+dy*dy),"surface piece lies outside its retained fan") ;
      }
      return selected ;
    }
    ExactPoint affineValue(const detail::ExactMapView &view,const ExactPoint &q) {
      auto p=view.origin ;
      for(unsigned d=0;d<3;++d) for(unsigned c=0;c<3;++c) p[d]+=view.columns[c][d]*q[c] ;
      return p ;
    }
    ExactPoint radialPosition(const Point &center,const std::array<ExactPoint,2> &edge,unsigned sector,const XY &q) {
      const auto &corners=referenceElement(Shape::quad).vertices ;
      const auto &a=corners[sector],&b=corners[(sector+1)%4] ;
      const Exact ax=Exact(a[0])-Exact(.5),ay=Exact(a[1])-Exact(.5),bx=Exact(b[0])-Exact(.5),by=Exact(b[1])-Exact(.5),det=ax*by-bx*ay ;
      const Exact x=q[0]-Exact(.5),y=q[1]-Exact(.5),u=(by*x-bx*y)/det,v=(ax*y-ay*x)/det ;
      ExactPoint p ; for(unsigned d=0;d<3;++d) p[d]=Exact(center[d])+u*(edge[0][d]-Exact(center[d]))+v*(edge[1][d]-Exact(center[d])) ;
      return p ;
    }
    MapValue boundaryFanValue(const GeometryMap &edge,const Point &center,unsigned sector,const Point &q,long double centerParameter=.5L) {
      finite(q) ; const auto &corners=referenceElement(Shape::quad).vertices ;
      const auto &a=corners[sector],&b=corners[(sector+1)%4] ;
      const long double ax=a[0]-.5L,ay=a[1]-.5L,bx=b[0]-.5L,by=b[1]-.5L,det=ax*by-bx*ay ;
      const long double x=static_cast<long double>(q[0])-.5L,y=static_cast<long double>(q[1])-.5L,
        u=(by*x-bx*y)/det,v=(ax*y-ay*x)/det,r=u+v,t=r==0 ? centerParameter:v/r ;
      const auto value=edge.evaluate({{double(t),0,0}}) ;
      MapValue out{{},{},value.hasJacobian} ;
      const std::array<long double,2> du{{by/det,-bx/det}},dv{{-ay/det,ax/det}} ;
      for(unsigned d=0;d<3;++d) {
        const long double delta=static_cast<long double>(value.position[d])-center[d] ;
        out.position[d]=r==1 ? value.position[d]:center[d]+r*delta ;
        if(out.hasJacobian) for(unsigned i=0;i<2;++i) {
          const auto dr=du[i]+dv[i] ;
          out.jacobian[i][d]=dr*delta+(dv[i]-t*dr)*value.jacobian[0][d] ;
        }
      }
      return out ;
    }
    // A selected edge branch extends over its whole radial sector.
    class BoundaryFanSectorMap final : public GeometryMap {
      Map edge_ ; Point center_ ; unsigned sector_ ; long double centerParameter_ ; detail::ExactMapView affine_ ;
      std::array<long double,3> origin_{} ;
      std::array<std::array<long double,3>,2> derivative_{} ;
    public:
      BoundaryFanSectorMap(Map edge,Point center,unsigned sector,const Exact &centerParameter)
        :edge_(std::move(edge)),center_(center),sector_(sector),centerParameter_(centerParameter.convert_to<long double>()) {
        std::array<ExactPoint,2> ends{{{{0,0,0}},{{1,0,0}}}} ; std::size_t remaining=4096 ;
        if(!detail::affineTrace(*edge_,ends,remaining)) return ;
        affine_.kind=detail::ExactMapView::Kind::affine ; affine_.origin=radialPosition(center_,ends,sector_,{{0,0}}) ;
        const std::array<XY,2> q{{{{1,0}},{{0,1}}}} ;
        for(unsigned i=0;i<2;++i) {
          const auto p=radialPosition(center_,ends,sector_,q[i]) ;
          for(unsigned d=0;d<3;++d) { affine_.columns[i][d]=p[d]-affine_.origin[d] ; derivative_[i][d]=affine_.columns[i][d].convert_to<long double>() ; }
        }
        for(unsigned d=0;d<3;++d) origin_[d]=affine_.origin[d].convert_to<long double>() ;
      }
      const detail::ExactMapView &algebra() const { return affine_ ; }
      const Map &edge() const { return edge_ ; }
      const Point &center() const { return center_ ; }
      unsigned sector() const { return sector_ ; }
      Shape shape() const override { return Shape::quad ; }
      MapValue evaluate(const Point &q) const override {
        if(affine_.kind!=detail::ExactMapView::Kind::affine) return boundaryFanValue(*edge_,center_,sector_,q,centerParameter_) ;
        finite(q) ; MapValue out{{},{},true} ;
        for(unsigned d=0;d<3;++d) {
          out.position[d]=origin_[d]+q[0]*derivative_[0][d]+q[1]*derivative_[1][d] ;
          for(unsigned i=0;i<2;++i) out.jacobian[i][d]=derivative_[i][d] ;
        }
        return out ;
      }
    } ;
    class FanSectorMap final : public GeometryMap {
      std::shared_ptr<const FanMap> fan_ ; unsigned sector_ ;
    public:
      FanSectorMap(std::shared_ptr<const FanMap> fan,unsigned sector):fan_(std::move(fan)),sector_(sector) {}
      detail::ExactMapView algebra() const {
        detail::ExactMapView view ; view.kind=detail::ExactMapView::Kind::affine ;
        view.origin=fanPosition(*fan_,sector_,{{0,0}}) ;
        const std::array<XY,2> q{{{{1,0}},{{0,1}}}} ;
        for(unsigned i=0;i<2;++i) { const auto p=fanPosition(*fan_,sector_,q[i]) ; for(unsigned d=0;d<3;++d) view.columns[i][d]=p[d]-view.origin[d] ; }
        return view ;
      }
      Shape shape() const override { return fan_->shape() ; }
      ReferenceCorners referenceCorners() const override { return fan_->referenceCorners() ; }
      bool containsReference(const Point &q,double t=0) const override { return fan_->containsReference(q,t) ; }
      MapValue evaluate(const Point &q) const override {
        finite(q) ; const auto r=fan_->referenceCorners() ; const unsigned i=sector_,j=(i+1)%r.size() ;
        return affineTriangle({{referenceCenter(shape()),r[i],r[j]}},{{fan_->center(),fan_->corners()[i],fan_->corners()[j]}},q) ;
      }
    } ;
    class CoonsPieceMap final : public GeometryMap {
      std::array<Map,4> edges_ ; std::array<Point,4> corners_ ;
    public:
      CoonsPieceMap(std::array<Map,4> edges,const std::array<Point,4> &corners):edges_(std::move(edges)),corners_(corners) {}
      const auto &edges() const { return edges_ ; }
      const auto &corners() const { return corners_ ; }
      Shape shape() const override { return Shape::quad ; }
      MapValue evaluate(const Point &q) const override { return detail::evaluateBoundarySurface(edges_,corners_,q) ; }
    } ;
    // Select the smooth branch containing q, retaining its reference coordinates.
    Map selectMapPiece(Map map,const ExactPoint &q,unsigned depth=0) {
      if(depth>128) throw std::length_error("surface map graph exceeds its depth limit") ;
      if(const auto *fan=dynamic_cast<const FanMap *>(map.get()))
        return std::make_shared<FanSectorMap>(std::dynamic_pointer_cast<const FanMap>(map),fanSector(fan->referenceCorners(),fan->shape(),q)) ;
      if(const auto *fan=dynamic_cast<const BoundaryFanMap *>(map.get())) {
        const unsigned sector=fanSector(fan->referenceCorners(),fan->shape(),q) ;
        const auto refs=fan->referenceCorners() ; const auto a=exactPolygon({refs[sector]})[0],b=exactPolygon({refs[(sector+1)%4]})[0] ;
        const Exact ax=a[0]-Exact(.5),ay=a[1]-Exact(.5),bx=b[0]-Exact(.5),by=b[1]-Exact(.5),det=ax*by-bx*ay ;
        const Exact x=q[0]-Exact(.5),y=q[1]-Exact(.5),u=(by*x-bx*y)/det,v=(ax*y-ay*x)/det,r=u+v,t=r==0 ? Exact(.5):Exact(v/r) ;
        return std::make_shared<BoundaryFanSectorMap>(selectMapPiece(fan->edges()[sector],{{t,0,0}},depth+1),fan->center(),sector,t) ;
      }
      if(const auto *patch=dynamic_cast<const SurfacePatchMap *>(map.get()))
        return std::make_shared<RestrictedMap>(selectMapPiece(patch->source(),exactChartPosition(*patch->placement(),q),depth+1),patch->placement()) ;
      if(const auto *r=dynamic_cast<const RestrictedMap *>(map.get())) {
        if(!hasSurfacePieces(*map)) return map ;
        return std::make_shared<RestrictedMap>(selectMapPiece(r->parent(),exactChartPosition(*r->referenceMap(),q),depth+1),r->referenceMap()) ;
      }
      if(const auto *surface=dynamic_cast<const BoundarySurfaceMap *>(map.get())) {
        std::array<Map,4> edges ; const std::array<Exact,4> t{{q[0],q[1],1-q[0],1-q[1]}} ;
        for(unsigned i=0;i<4;++i) edges[i]=selectMapPiece(surface->edges()[i],{{t[i],0,0}},depth+1) ;
        return std::make_shared<CoonsPieceMap>(std::move(edges),surface->corners()) ;
      }
      return map ;
    }
    bool affine(const GeometryMap &map) {
      if(const auto *r=dynamic_cast<const RestrictedMap *>(&map)) return affine(*r->parent()) && affine(*r->referenceMap()) ;
      if(const auto *p=dynamic_cast<const PolygonMap *>(&map)) {
        const auto refs=p->referenceCorners() ; return p->corners()==std::vector<Point>(refs.begin(),refs.end()) ;
      }
      const auto *c=dynamic_cast<const CornerMap *>(&map) ; if(!c) return false ;
      if(c->shape()==Shape::triangle || c->shape()==Shape::segment) return true ;
      if(c->shape()!=Shape::quad) return false ;
      const auto &p=c->corners() ;
      for(unsigned d=0;d<3;++d) if(Exact(p[0][d])+Exact(p[2][d])!=Exact(p[1][d])+Exact(p[3][d])) return false ;
      return true ;
    }
    bool fullReorientation(const GeometryMap &parent,const CornerMap &placement) {
      if(placement.shape()!=parent.shape()) return false ;
      const auto refs=parent.referenceCorners() ; auto a=placement.corners() ; std::vector<Point> b(refs.begin(),refs.end()) ;
      std::sort(a.begin(),a.end()) ; std::sort(b.begin(),b.end()) ; return a==b ;
    }
    struct Chart {
      Map source ; ///< Original retained surface.
      Map placement ; ///< Requested surface coordinates to source coordinates.
      Polygon domain ; ///< Patch boundary in source coordinates.
    } ;
    Chart chart(Map map,unsigned depth=0) {
      if(depth>128) throw std::length_error("surface map graph exceeds its depth limit") ;
      if(const auto *patch=dynamic_cast<const SurfacePatchMap *>(map.get()))
        return {patch->source(),patch->placement(),exactPolygon(patch->chartCorners())} ;
      if(const auto *r=dynamic_cast<const RestrictedMap *>(map.get())) {
        const auto *p=dynamic_cast<const CornerMap *>(r->referenceMap().get()) ;
        require(p && p->dimension()==2 && r->parent()->dimension()==2 && affine(*p),"piecewise surface restriction requires an affine placement") ;
        auto base=chart(r->parent(),depth+1) ;
        require(affine(*base.placement) || fullReorientation(*r->parent(),*p),"nonlinear surface restrictions require a source-chart patch") ;
        auto placement=std::make_shared<RestrictedMap>(base.placement,r->referenceMap()) ;
        std::vector<Point> corners ; for(const auto &q:p->referenceCorners()) corners.push_back(placement->position(q)) ;
        base.domain=exactPolygon(corners) ; requireConvex(base.domain) ; base.placement=std::move(placement) ; return base ;
      }
      require(map->dimension()==2,"surface pieces require a surface") ;
      const auto r=map->referenceCorners() ; std::vector<Point> p(r.begin(),r.end()) ;
      return {map,coordinates(p),exactPolygon(p)} ;
    }
    std::vector<double> curveBreaks(Map map,std::size_t limit,unsigned depth) ;
    // Partition the source reference domain along every retained surface seam.
    std::vector<Polygon> regions(Map source,std::size_t limit,unsigned depth) {
      if(depth>128) throw std::length_error("surface map graph exceeds its depth limit") ;
      std::vector<Polygon> out ;
      const auto add=[&](Polygon p) { if(out.size()>=limit) throw std::length_error("surface piece count exceeds its limit") ; out.push_back(std::move(p)) ; } ;
      if(const auto *fan=dynamic_cast<const FanMap *>(source.get())) {
        const auto r=fan->referenceCorners() ; const auto c=referenceCenter(fan->shape()) ;
        for(unsigned i=0;i<r.size();++i) add(exactPolygon({c,r[i],r[(i+1)%r.size()]})) ;
      } else if(const auto *fan=dynamic_cast<const BoundaryFanMap *>(source.get())) {
        const auto refs=fan->referenceCorners() ; const XY center{{Exact(.5),Exact(.5)}} ;
        for(unsigned e=0;e<4;++e) {
          const auto a=exactPolygon({refs[e]})[0],b=exactPolygon({refs[(e+1)%4]})[0] ;
          const auto cuts=curveBreaks(fan->edges()[e],limit,depth+1) ;
          for(unsigned i=1;i<cuts.size();++i) {
            Polygon p{center,{}, {}} ;
            for(unsigned j=0;j<2;++j) for(unsigned d=0;d<2;++d) p[j+1][d]=a[d]+Exact(cuts[i-1+j])*(b[d]-a[d]) ;
            add(std::move(p)) ;
          }
        }
      } else if(const auto *surface=dynamic_cast<const BoundarySurfaceMap *>(source.get())) {
        std::vector<double> u{0,1},v{0,1} ;
        for(unsigned e=0;e<4;++e) {
          const auto cuts=curveBreaks(surface->edges()[e],limit,depth+1) ; auto &axis=e%2 ? v:u ;
          for(double t:cuts) axis.push_back(e>=2 ? 1-t:t) ;
        }
        for(auto *axis:{&u,&v}) { std::sort(axis->begin(),axis->end()) ; axis->erase(std::unique(axis->begin(),axis->end()),axis->end()) ; }
        for(unsigned j=1;j<v.size();++j) for(unsigned i=1;i<u.size();++i)
          add(exactPolygon({{{u[i-1],v[j-1],0}},{{u[i],v[j-1],0}},{{u[i],v[j],0}},{{u[i-1],v[j],0}}})) ;
      } else {
        const auto r=source->referenceCorners() ; add(exactPolygon({r.begin(),r.end()})) ;
      }
      return out ;
    }
    // Return ordered curve parameters at intersections with source-surface seams.
    std::vector<double> curveBreaks(Map map,std::size_t limit,unsigned depth) {
      if(depth>128) throw std::length_error("surface map graph exceeds its depth limit") ;
      std::vector<double> result{0,1} ;
      if(!hasSurfacePieces(*map)) return result ;
      const auto *r=dynamic_cast<const RestrictedMap *>(map.get()) ; if(!r) return result ;
      const auto *p=dynamic_cast<const CornerMap *>(r->referenceMap().get()) ;
      require(p && p->shape()==Shape::segment,"piecewise edge requires an affine segment placement") ;
      const auto &ends=p->corners() ;
      if(r->parent()->dimension()==1) {
        const double a=ends[0][0],b=ends[1][0] ; require(a!=b,"piecewise edge has a zero interval") ;
        for(double t:curveBreaks(r->parent(),limit,depth+1)) { const double q=(t-a)/(b-a) ; if(q>0 && q<1) result.push_back(q) ; }
      } else if(hasSurfacePieces(*r->parent())) {
        auto parent=chart(r->parent(),depth+1) ;
        require(affine(*parent.placement),"piecewise edge must use retained source coordinates") ;
        const auto a=exactPolygon({parent.placement->position(ends[0])})[0],b=exactPolygon({parent.placement->position(ends[1])})[0] ;
        for(const auto &region:regions(parent.source,limit,depth+1)) for(unsigned e=0;e<region.size();++e) {
          const auto &u=region[e],&v=region[(e+1)%region.size()] ; const Exact x=turn(u,v,a),y=turn(u,v,b) ;
          if(x==y) continue ;
          const Exact t=x/(x-y) ;
          if(t>0 && t<1) result.push_back(t.convert_to<double>()) ;
        }
      }
      std::sort(result.begin(),result.end()) ; result.erase(std::unique(result.begin(),result.end()),result.end()) ;
      if(result.size()>limit+1) throw std::length_error("edge piece count exceeds its limit") ;
      return result ;
    }
  }

  detail::ExactMapView detail::exactMapView(const GeometryMap &map) {
    if(const auto *affine=dynamic_cast<const AffineTriangleMap *>(&map)) return affine->algebra() ;
    if(const auto *fan=dynamic_cast<const FanSectorMap *>(&map)) return fan->algebra() ;
    if(const auto *fan=dynamic_cast<const BoundaryFanSectorMap *>(&map)) return fan->algebra() ;
    ExactMapView view ;
    if(const auto *smooth=dynamic_cast<const SmoothMap *>(&map)) {
      view.kind=ExactMapView::Kind::delegate ; view.delegate=&smooth->map() ;
    } else if(const auto *coons=dynamic_cast<const CoonsPieceMap *>(&map)) {
      view.kind=ExactMapView::Kind::coons ; view.edges=&coons->edges() ; view.corners=&coons->corners() ;
    }
    return view ;
  }

  bool detail::affineTrace(const GeometryMap &geometry,std::array<ExactPoint,2> &line,std::size_t &remaining) {
    const GeometryMap *map=&geometry ; std::vector<const GeometryMap *> parents ;
    const auto inside=[](Shape shape,const ExactPoint &p) {
      if(p[0]<0 || p[1]<0 || p[2]!=0) return false ;
      if(shape==Shape::segment) return p[0]<=1 && p[1]==0 ;
      if(shape==Shape::triangle) return p[0]+p[1]<=1 ;
      return shape==Shape::quad && p[0]<=1 && p[1]<=1 ;
    } ;
    while(true) {
      if(remaining==0) return false ;
      --remaining ;
      if(const auto *child=dynamic_cast<const RestrictedMap *>(map)) {
        parents.push_back(child->parent().get()) ; map=child->referenceMap().get() ; continue ;
      }
      const auto view=exactMapView(*map) ;
      if(view.kind==ExactMapView::Kind::delegate) { map=view.delegate ; continue ; }
      const auto *corners=dynamic_cast<const CornerMap *>(map) ;
      if(view.kind==ExactMapView::Kind::affine) {
        for(auto &q:line) q=affineValue(view,q) ;
      } else if(!corners || !inside(corners->shape(),line[0]) || !inside(corners->shape(),line[1])) return false ;
      else {
        if(corners->shape()==Shape::quad) {
          const auto &p=corners->corners() ; const Exact product=(line[1][0]-line[0][0])*(line[1][1]-line[0][1]) ;
          for(unsigned d=0;d<3;++d) if(product*(Exact(p[2][d])-Exact(p[1][d])-Exact(p[3][d])+Exact(p[0][d]))!=0) return false ;
        }
        for(auto &q:line) q=exactChartPosition(*corners,q) ;
      }
      if(parents.empty()) return true ;
      map=parents.back() ; parents.pop_back() ;
    }
  }

  bool detail::radialConeSamples(const ConeMap &cone,std::vector<std::function<MeasureSample(const Point&)>> &samples) {
    samples.clear() ;
    const GeometryMap *map=cone.base().get() ; if(map->shape()!=Shape::triangle) return false ;
    std::array<ExactPoint,3> points{{{{0,0,0}},{{1,0,0}},{{0,1,0}}}} ;
    for(unsigned depth=0;depth<128;++depth) {
      const auto view=exactMapView(*map) ;
      if(view.kind==ExactMapView::Kind::delegate) { map=view.delegate ; continue ; }
      const auto *restriction=dynamic_cast<const RestrictedMap *>(map) ; if(!restriction) break ;
      auto placement=exactMapView(*restriction->referenceMap()) ;
      if(placement.kind!=ExactMapView::Kind::affine) {
        const auto *corners=dynamic_cast<const CornerMap *>(restriction->referenceMap().get()) ;
        if(!corners || corners->shape()!=Shape::triangle) return false ;
        placement.kind=ExactMapView::Kind::affine ;
        for(unsigned d=0;d<3;++d) {
          placement.origin[d]=Exact(corners->corners()[0][d]) ;
          for(unsigned i=0;i<2;++i) placement.columns[i][d]=Exact(corners->corners()[i+1][d])-placement.origin[d] ;
        }
      }
      for(auto &q:points) q=affineValue(placement,q) ;
      map=restriction->parent().get() ;
    }
    const auto *radial=dynamic_cast<const BoundaryFanSectorMap *>(map) ;
    if(!radial || radial->algebra().kind==ExactMapView::Kind::affine) return false ;
    const auto &corners=referenceElement(Shape::quad).vertices ;
    const auto &a=corners[radial->sector()],&b=corners[(radial->sector()+1)%4] ;
    const Exact ax=Exact(a[0])-Exact(.5),ay=Exact(a[1])-Exact(.5),bx=Exact(b[0])-Exact(.5),by=Exact(b[1])-Exact(.5),det=ax*by-bx*ay ;

    // A ray through a triangle changes its intersected edges only at a vertex direction.
    std::array<XY,3> uv ; std::array<Exact,3> radius,parameter ; std::vector<Exact> cuts ;
    for(unsigned i=0;i<3;++i) {
      if(points[i][2]!=0) return false ;
      const Exact x=points[i][0]-Exact(.5),y=points[i][1]-Exact(.5) ;
      uv[i]={{(by*x-bx*y)/det,(ax*y-ay*x)/det}} ; radius[i]=uv[i][0]+uv[i][1] ;
      if(radius[i]<0) return false ;
      if(radius[i]==0) { if(x!=0 || y!=0) return false ; }
      else {
        parameter[i]=uv[i][1]/radius[i] ; if(parameter[i]<0 || parameter[i]>1) return false ;
        cuts.push_back(parameter[i]) ;
      }
    }
    std::sort(cuts.begin(),cuts.end()) ; cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end()) ;
    if(cuts.size()<2) return false ;
    const Exact orientation=turn({{points[0][0],points[0][1]}},{{points[1][0],points[1][1]}},{{points[2][0],points[2][1]}}) ;
    if(orientation==0) return false ;
    struct Bound { Exact n,b,m ; Exact at(const Exact &t) const { return n/(b+m*t) ; } } ;
    std::vector<std::function<MeasureSample(const Point&)>> result ;
    for(unsigned interval=0;interval+1<cuts.size();++interval) {
      const Exact first=cuts[interval],last=cuts[interval+1],middle=(first+last)/2,span=last-first ;
      std::vector<Bound> hits ; bool touchesCenter=false ;
      for(unsigned i=0;i<3;++i) {
        const unsigned j=(i+1)%3 ;
        if(radius[i]==0 || radius[j]==0) { touchesCenter=true ; continue ; }
        if(!(std::min(parameter[i],parameter[j])<middle && middle<std::max(parameter[i],parameter[j]))) continue ;
        const auto &a=uv[i],&b=uv[j] ;
        Bound bound{a[0]*b[1]-a[1]*b[0],b[1]-a[1],a[0]+a[1]-b[0]-b[1]} ;
        if(bound.at(middle)<=0) return false ;
        hits.push_back(std::move(bound)) ;
      }
      if(touchesCenter) hits.push_back({0,1,0}) ;
      if(hits.size()!=2) return false ;
      if(hits[0].at(middle)>hits[1].at(middle)) std::swap(hits[0],hits[1]) ;
      std::array<long double,3> lo,hi ;
      for(unsigned k=0;k<2;++k) {
        const auto &x=hits[k] ; auto &v=k ? hi:lo ;
        v={{x.n.convert_to<long double>(),(x.b+x.m*first).convert_to<long double>(),(x.m*span).convert_to<long double>()}} ;
      }
      // Preserve the thin radial width before rounding the two individual bounds.
      const Exact n0=hits[1].n*(hits[0].b+hits[0].m*first)-hits[0].n*(hits[1].b+hits[1].m*first) ;
      const Exact n1=(hits[1].n*hits[0].m-hits[0].n*hits[1].m)*span ;
      const long double width0=n0.convert_to<long double>(),width1=n1.convert_to<long double>() ;
      const long double begin=first.convert_to<long double>(),length=span.convert_to<long double>(),sign=orientation>0 ? 1:-1 ;
      const auto edge=radial->edge() ; const Point c=radial->center(),apex=cone.apex() ;
      result.push_back([edge,c,apex,lo,hi,width0,width1,begin,length,sign](const Point &q) {
        const long double p=q[0],dl=lo[1]+lo[2]*p,dh=hi[1]+hi[2]*p ;
        const long double lower=lo[0]/dl,upper=hi[0]/dh,width=(width0+width1*p)/(dl*dh),sum=upper+lower ;
        if(!(width>0 && lower>=0 && sum>0)) throw std::invalid_argument("clipped radial interval is invalid") ;
        const auto value=edge->evaluate({{double(begin+length*p),0,0}}) ;
        if(!value.hasJacobian) throw std::invalid_argument("radial cone measure requires edge derivatives") ;
        std::array<long double,3> delta,height,tangent ; Point centroid ;
        const long double moment=.5L*(upper*upper+upper*lower+lower*lower)/sum ;
        for(unsigned d=0;d<3;++d) {
          delta[d]=static_cast<long double>(value.position[d])-c[d] ; height[d]=static_cast<long double>(apex[d])-c[d] ; tangent[d]=value.jacobian[0][d] ;
          centroid[d]=apex[d]-.75L*height[d]+moment*delta[d] ;
        }
        return MeasureSample{centroid,sign*length*width*sum*jacobianDeterminant(delta,tangent,height)/6} ;
      }) ;
    }
    samples=std::move(result) ; return true ;
  }

  FanMap::FanMap(Shape shape,std::vector<Point> corners,Point middle):shape_(shape),corners_(std::move(corners)),center_(middle) {
    require(shape==Shape::triangle || shape==Shape::quad || shape==Shape::polygon,"fan requires a surface shape") ;
    require(shape==Shape::polygon ? corners_.size()>=3:corners_.size()==referenceElement(shape).vertexCount,"fan corner count differs from its shape") ;
    finite(center_) ; for(const auto &p:corners_) finite(p) ; validateFan(corners_,center_) ; reference_=reference(shape,corners_.size()) ;
  }
  bool FanMap::containsReference(const Point &p,double tolerance) const {
    return shape_==Shape::polygon ? containsReferencePolygon(corners_.size(),p,tolerance):GeometryMap::containsReference(p,tolerance) ;
  }
  MapValue FanMap::evaluate(const Point &p) const {
    finite(p) ; const auto c=referenceCenter(shape_) ; unsigned selected=0 ; long double best=-std::numeric_limits<long double>::infinity() ;
    const XY x{{Exact(p[0]),Exact(p[1])}} ;
    for(unsigned i=0;i<reference_.size();++i) {
      const auto triangle=exactPolygon({c,reference_[i],reference_[(i+1)%reference_.size()]}) ;
      const auto sign=twiceSignedArea(triangle)>0 ? 1:-1 ; bool inside=true ; long double score=std::numeric_limits<long double>::infinity() ;
      for(unsigned e=0;e<3;++e) { const Exact t=sign*turn(triangle[e],triangle[(e+1)%3],x) ; inside=inside && t>=0 ; score=std::min(score,t.convert_to<long double>()) ; }
      if(inside) { selected=i ; break ; } if(score>best) { best=score ; selected=i ; }
    }
    return affineTriangle({{c,reference_[selected],reference_[(selected+1)%reference_.size()]}},
      {{center_,corners_[selected],corners_[(selected+1)%corners_.size()]}},p) ;
  }
  BoundaryFanMap::BoundaryFanMap(std::array<Map,4> edges,Point center):edges_(std::move(edges)),center_(center) {
    finite(center_) ; std::array<Point,4> corners ;
    for(unsigned e=0;e<4;++e) {
      require(edges_[e] && edges_[e]->shape()==Shape::segment,"boundary fan requires four segment maps") ;
      corners[e]=edges_[e]->position({{0,0,0}}) ; finite(corners[e]) ;
    }
    for(unsigned e=0;e<4;++e) {
      const auto end=edges_[e]->position({{1,0,0}}) ; finite(end) ;
      for(unsigned d=0;d<3;++d) require(std::abs(end[d]-corners[(e+1)%4][d])<=64*std::numeric_limits<double>::epsilon()*
        std::max({1.,std::abs(end[d]),std::abs(corners[(e+1)%4][d])}),"boundary fan edges do not meet") ;
    }
  }
  MapValue BoundaryFanMap::evaluate(const Point &q) const {
    finite(q) ;
    const long double x=static_cast<long double>(q[0])-.5L,y=static_cast<long double>(q[1])-.5L ;
    const unsigned sector=y<=-std::abs(x) ? 0:x>=std::abs(y) ? 1:y>=std::abs(x) ? 2:3 ;
    return boundaryFanValue(*edges_[sector],center_,sector,q) ;
  }
  SurfacePatchMap::SurfacePatchMap(Map source,std::vector<Point> corners):source_(std::move(source)),corners_(std::move(corners)) {
    require(source_ && source_->dimension()==2 && !dynamic_cast<const SurfacePatchMap *>(source_.get()),"surface patch requires an unpatched surface source") ;
    require(!hasSurfacePieces(*source_) || dynamic_cast<const FanMap *>(source_.get()) || dynamic_cast<const BoundaryFanMap *>(source_.get()) || dynamic_cast<const BoundarySurfaceMap *>(source_.get()),
            "surface patch requires a direct retained surface source") ;
    require(corners_.size()>=3,"surface patch requires at least three corners") ; requireConvex(exactPolygon(corners_)) ;
    for(const auto &p:corners_) require(source_->containsReference(p,64*std::numeric_limits<double>::epsilon()),"surface patch lies outside its source") ;
    placement_=coordinates(corners_) ;
  }
  bool SurfacePatchMap::containsReference(const Point &p,double tolerance) const { return placement_->containsReference(p,tolerance) ; }
  Point SurfacePatchMap::chartPosition(const Point &p) const { return placement_->position(p) ; }
  MapValue SurfacePatchMap::evaluate(const Point &p) const { const auto q=placement_->evaluate(p) ; return compose(source_->evaluate(q.position),q) ; }
  Point SurfacePatchMap::position(const Point &p) const { return source_->position(chartPosition(p)) ; }

  bool hasSurfacePieces(const GeometryMap &map) {
    if(dynamic_cast<const FanMap *>(&map) || dynamic_cast<const BoundaryFanMap *>(&map) || dynamic_cast<const SurfacePatchMap *>(&map)) return true ;
    if(const auto *r=dynamic_cast<const RestrictedMap *>(&map)) return hasSurfacePieces(*r->parent()) || hasSurfacePieces(*r->referenceMap()) ;
    if(const auto *s=dynamic_cast<const BoundarySurfaceMap *>(&map)) for(const auto &e:s->edges()) if(hasSurfacePieces(*e)) return true ;
    return false ;
  }

  SurfacePieces smoothSurfacePieces(Map map,std::size_t maxPieces) {
    require(maxPieces>0,"surface pieces require a positive work limit") ;
    auto source=chart(std::move(map)) ;
    SurfacePieces out{source.source,source.placement,{}} ;
    const auto partitions=regions(source.source,maxPieces,0) ;
    const auto *fan=dynamic_cast<const FanMap *>(source.source.get()) ;
    const bool reversed=twiceSignedArea(source.domain)<0 ;
    for(unsigned sector=0;sector<partitions.size();++sector) {
      // Clip the requested patch against each retained surface piece.
      auto polygon=clip(source.domain,partitions[sector]) ;
      if(polygon.size()<3 || twiceSignedArea(polygon)==0) continue ;
      if((twiceSignedArea(polygon)<0)!=reversed) std::reverse(polygon.begin(),polygon.end()) ;
      ExactPoint middle{} ;
      for(const auto &p:partitions[sector])
        for(unsigned d=0;d<2;++d) middle[d]+=p[d]/partitions[sector].size() ;
      const auto smooth=fan ? Map{}:selectMapPiece(source.source,middle) ;
      for(unsigned i=1;i+1<polygon.size();++i) {
        if(turn(polygon[0],polygon[i],polygon[i+1])==0) continue ;
        if(out.pieces.size()>=maxPieces) throw std::length_error("surface piece count exceeds its limit") ;

        // Order each triangle for integration without changing its orientation.
        std::array<XY,3> q{{polygon[0],polygon[i],polygon[i+1]}} ;
        unsigned shortest=0 ;
        Exact minimum ;
        for(unsigned e=0;e<3;++e) {
          const Exact x=q[(e+1)%3][0]-q[e][0],y=q[(e+1)%3][1]-q[e][1],length=x*x+y*y ;
          if(e==0 || length<minimum) { shortest=e ; minimum=length ; }
        }
        if(shortest) std::rotate(q.begin(),q.begin()+shortest,q.end()) ;
        // Nonlinear radial pieces put their center at the quadrature collapse vertex.
        // Affine pieces retain a short first edge to avoid nearly parallel derivatives.
        if(dynamic_cast<const BoundaryFanMap *>(source.source.get()) && detail::exactMapView(*smooth).kind!=detail::ExactMapView::Kind::affine)
          for(unsigned j=0;j<3;++j) if(q[j]==XY{{Exact(.5),Exact(.5)}}) {
            std::rotate(q.begin(),q.begin()+(j+2)%3,q.end()) ;
            break ;
          }

        // Retain the source geometry on this triangle, with exact coordinates for affine pieces.
        std::array<ExactPoint,3> chartPoints ;
        for(unsigned j=0;j<3;++j) chartPoints[j]={{q[j][0],q[j][1],0}} ;
        auto placement=std::make_shared<AffineTriangleMap>(chartPoints) ;
        Map physical ;
        if(fan) {
          std::array<ExactPoint,3> p ;
          for(unsigned j=0;j<3;++j) p[j]=fanPosition(*fan,sector,q[j]) ;
          physical=std::make_shared<AffineTriangleMap>(p) ;
        } else {
          auto active=smooth ;
          if(dynamic_cast<const BoundaryFanMap *>(source.source.get())) {
            ExactPoint representative{} ;
            for(const auto &p:chartPoints)
              for(unsigned d=0;d<2;++d) representative[d]+=p[d]/3 ;
            active=selectMapPiece(source.source,representative) ;
          }
          const auto algebra=detail::exactMapView(*active) ;
          if(algebra.kind==detail::ExactMapView::Kind::affine) {
            std::array<ExactPoint,3> p ;
            for(unsigned j=0;j<3;++j) p[j]=affineValue(algebra,chartPoints[j]) ;
            physical=std::make_shared<AffineTriangleMap>(p) ;
          } else physical=std::make_shared<SmoothMap>(std::make_shared<RestrictedMap>(active,placement)) ;
        }
        out.pieces.push_back({std::move(physical),std::move(placement)}) ;
      }
    }
    require(!out.pieces.empty(),"surface pieces have zero area") ;
    return out ;
  }

  Map reverseSurfacePiece(Map map) {
    require(map && map->dimension()==2,"surface piece reversal requires a surface") ;
    const auto refs=map->referenceCorners() ; std::vector<Point> reversed ;
    if(map->shape()==Shape::triangle && !dynamic_cast<const AffineTriangleMap *>(map.get()))
      reversed={refs[2],refs[1],refs[0]} ;
    else for(unsigned i=0;i<refs.size();++i) reversed.push_back(refs[(refs.size()-i)%refs.size()]) ;
    return std::make_shared<RestrictedMap>(std::move(map),coordinates(reversed)) ;
  }

  std::vector<Map> smoothCurvePieces(Map map,std::size_t maxPieces) {
    require(map && map->dimension()==1 && maxPieces>0,"curve pieces require a segment and a positive work limit") ;
    const auto cuts=curveBreaks(map,maxPieces,0) ; std::vector<Map> out ;
    for(unsigned i=1;i<cuts.size();++i) {
      const Exact middle=(Exact(cuts[i-1])+Exact(cuts[i]))/2 ;
      out.push_back(std::make_shared<SmoothMap>(std::make_shared<RestrictedMap>(selectMapPiece(map,{{middle,0,0}}),
        std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{{{cuts[i-1],0,0}},{{cuts[i],0,0}}})))) ;
    }
    return out ;
  }

  Point inverseSurfaceCoordinates(const GeometryMap &placement,const Point &target) {
    finite(target) ; Point q=referenceCenter(placement.shape()) ;
    for(unsigned iteration=0;iteration<64;++iteration) {
      const auto value=placement.evaluate(q) ; require(value.hasJacobian,"surface chart inverse requires derivatives") ;
      const long double a=value.jacobian[0][0],b=value.jacobian[1][0],c=value.jacobian[0][1],d=value.jacobian[1][1],det=a*d-b*c ;
      require(det!=0 && std::isfinite(det),"surface chart inverse is singular") ;
      const long double x=static_cast<long double>(value.position[0])-target[0],y=static_cast<long double>(value.position[1])-target[1],u=(d*x-b*y)/det,v=(a*y-c*x)/det ;
      if(std::hypot(x,y)<=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(target[0]),std::abs(target[1])})) {
        require(placement.containsReference(q,64*std::numeric_limits<double>::epsilon()),"surface inverse lies outside its chart") ; return q ;
      }
      bool accepted=false ;
      for(unsigned backtrack=0;backtrack<48;++backtrack) {
        const auto t=std::ldexp(1.L,-int(backtrack)) ; const Point trial{{double(q[0]-t*u),double(q[1]-t*v),0}} ;
        if(!placement.containsReference(trial,64*std::numeric_limits<double>::epsilon())) continue ;
        const auto p=placement.position(trial) ;
        if(std::hypot(static_cast<long double>(p[0])-target[0],static_cast<long double>(p[1])-target[1])<std::hypot(x,y)) { q=trial ; accepted=true ; break ; }
      }
      if(!accepted) break ;
    }
    throw std::invalid_argument("surface chart inverse did not converge") ;
  }
}
}
