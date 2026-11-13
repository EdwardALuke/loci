// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/mesh_quality.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace Loci {
namespace FVMAdapt {
  namespace {
    Point subtract(const Point &a,const Point &b) { return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}} ; }
    double dot(const Point &a,const Point &b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] ; }
    Point cross(const Point &a,const Point &b) { return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ; }
    double length(const Point &p) { return std::hypot(p[0],p[1],p[2]) ; }
    void require(bool value,const char *message) { if(!value) throw std::invalid_argument(message) ; }
    double angle(const Point &a,const Point &b) {
      const double na=length(a),nb=length(b) ;
      require(na>0 && nb>0 && std::isfinite(na) && std::isfinite(nb),"undefined mesh quality angle") ;
      double cosine=0 ; for(unsigned d=0;d<3;++d) cosine+=(a[d]/na)*(b[d]/nb) ;
      return std::acos(std::clamp(cosine,-1.,1.))*180/std::acos(-1.) ;
    }
    bool locationLess(const QualityValue &a,const QualityValue &b) {
      return std::tie(a.root,a.leaf,a.face)<std::tie(b.root,b.leaf,b.face) ;
    }
  }
  void MeshQualityOptions::validate() const {
    require(mode=="off" || mode=="report" || mode=="reject","quality mode must be off, report, or reject") ;
    require(centroid=="wireframe" || centroid=="exact","quality centroid must be wireframe or exact") ;
    for(double limit:{maxNonOrthogonality,maxFaceWarpage})
      require(std::isfinite(limit) && (limit==-1 || (limit>=0 && limit<=180)),"quality angle limits must be -1 or in [0,180] degrees") ;
    require(std::isfinite(maxChildVolumeRatio) && (maxChildVolumeRatio==-1 || maxChildVolumeRatio>=1),"child volume ratio limit must be -1 or at least one") ;
  }
  FaceMeasure measurePolygon(const std::vector<Point> &points,const std::string &centroid) {
    require(points.size()>=3,"quality face needs at least three vertices") ;
    require(centroid=="wireframe" || centroid=="exact","unknown quality centroid") ;
    const auto origin=points.front() ; std::vector<Point> p ; p.reserve(points.size()) ;
    for(const auto &point:points) {
      for(double x:point) require(std::isfinite(x),"nonfinite quality vertex") ;
      p.push_back(subtract(point,origin)) ;
    }
    FaceMeasure result ; Point wire{},center{} ; double perimeter=0 ;
    for(unsigned i=0;i<p.size();++i) {
      const auto &a=p[i],&b=p[(i+1)%p.size()] ; const double edge=length(subtract(b,a)) ;
      require(edge>0 && std::isfinite(edge),"zero or nonfinite quality edge") ; perimeter+=edge ;
      for(unsigned d=0;d<3;++d) wire[d]+=.5*edge*(a[d]+b[d]) ;
    }
    for(unsigned d=0;d<3;++d) wire[d]/=perimeter ;
    center=wire ;
    if(centroid=="exact") {
      // Match the solver's four updates of the area-weighted triangle-fan center.
      center=Point{} ; for(const auto &a:p) for(unsigned d=0;d<3;++d) center[d]+=a[d]/p.size() ;
      for(unsigned iteration=0;iteration<4;++iteration) {
        Point sum{} ; double area=0 ;
        for(unsigned i=0;i<p.size();++i) {
          const auto &a=p[i],&b=p[(i+1)%p.size()] ; const double weight=length(cross(subtract(a,center),subtract(b,center))) ;
          for(unsigned d=0;d<3;++d) sum[d]+=weight*(a[d]+b[d]+center[d])/3 ;
          area+=weight ;
        }
        require(area>0 && std::isfinite(area),"zero or nonfinite quality face area") ;
        for(unsigned d=0;d<3;++d) center[d]=sum[d]/area ;
      }
    }
    std::vector<Point> triangles ;
    for(unsigned i=0;i<p.size();++i) {
      auto triangle=cross(subtract(p[i],center),subtract(p[(i+1)%p.size()],center)) ;
      for(unsigned d=0;d<3;++d) { triangle[d]*=.5 ; result.area[d]+=triangle[d] ; }
      triangles.push_back(triangle) ;
      const auto a=subtract(p[i],wire),b=subtract(p[(i+1)%p.size()],wire),normal=cross(a,b) ;
      for(unsigned d=0;d<3;++d) for(unsigned j=0;j<3;++j) result.moment[3*d+j]+=(a[d]+b[d])*normal[j]/24 ;
    }
    require(length(result.area)>0 && std::isfinite(length(result.area)),"zero or nonfinite quality face area") ;
    for(const auto &triangle:triangles) if(length(triangle)>0) result.warpage=std::max(result.warpage,angle(result.area,triangle)) ;
    for(unsigned d=0;d<3;++d) { result.center[d]=origin[d]+center[d] ; result.wireCenter[d]=origin[d]+wire[d] ; }
    return result ;
  }
  CellMeasure measurePolyhedron(const std::vector<MeasuredFace> &faces,const std::string &centroid) {
    require(faces.size()>=4,"quality cell has an incomplete boundary") ;
    require(centroid=="wireframe" || centroid=="exact","unknown quality centroid") ;
    const auto origin=faces.front().face->center ; Point base{} ; double areaSum=0 ;
    for(const auto &use:faces) {
      require(use.sign==1 || use.sign==-1,"invalid quality face orientation") ;
      const auto &f=*use.face ; const double area=length(f.area) ; areaSum+=area ;
      for(unsigned d=0;d<3;++d) base[d]+=area*(f.center[d]-origin[d]) ;
    }
    require(areaSum>0 && std::isfinite(areaSum),"invalid quality cell area") ;
    for(unsigned d=0;d<3;++d) base[d]=origin[d]+base[d]/areaSum ;
    CellMeasure result ; result.center=base ; long double volume=0,exactVolume=0 ; long double moment[3]{} ;
    for(const auto &use:faces) {
      const auto &f=*use.face ; const auto w=subtract(f.wireCenter,base) ;
      volume+=use.sign*dot(subtract(f.center,base),f.area)/3 ;
      if(centroid=="exact") {
        const double v=use.sign*dot(w,f.area)/3 ; exactVolume+=v ;
        for(unsigned d=0;d<3;++d) {
          moment[d]+=.75*v*w[d] ;
          for(unsigned j=0;j<3;++j) moment[d]+=use.sign*f.moment[3*d+j]*w[j] ;
        }
      }
    }
    result.volume=volume ; require(result.volume>0 && std::isfinite(result.volume),"nonpositive or nonfinite emitted cell volume") ;
    if(centroid=="exact") {
      require(exactVolume>0 && std::isfinite(exactVolume),"nonpositive exact-center volume") ;
      for(unsigned d=0;d<3;++d) result.center[d]+=moment[d]/exactVolume ;
    }
    for(double x:result.center) require(std::isfinite(x),"nonfinite emitted cell center") ;
    return result ;
  }
  double faceNonOrthogonality(const FaceMeasure &face,const Point &left,const Point &right) {
    return angle(face.area,subtract(right,left)) ;
  }
  void includeQualityValue(QualityValue &out,double value,int root,int leaf,int face) {
    require(std::isfinite(value),"nonfinite mesh quality measurement") ; const QualityValue next{value,root,leaf,face} ;
    if(out.root<0 || value>out.value || (value==out.value && locationLess(next,out))) out=next ;
  }
  void MeshQualityReport::include(const MeshQualityReport &other) {
    for(auto entry:{std::make_pair(&nonOrthogonality,&other.nonOrthogonality),std::make_pair(&faceWarpage,&other.faceWarpage),
                    std::make_pair(&childVolumeRatio,&other.childVolumeRatio)})
      if(entry.second->root>=0) includeQualityValue(*entry.first,entry.second->value,entry.second->root,entry.second->leaf,entry.second->face) ;
    if(other.invalid() && (!invalid() || std::tie(other.invalidRoot,other.invalidLeaf,other.invalidFace)<std::tie(invalidRoot,invalidLeaf,invalidFace))) {
      invalidRoot=other.invalidRoot ; invalidLeaf=other.invalidLeaf ; invalidFace=other.invalidFace ;
    }
  }
  bool MeshQualityReport::exceeds(const MeshQualityOptions &o) const {
    return (o.maxNonOrthogonality>=0 && nonOrthogonality.value>o.maxNonOrthogonality) ||
      (o.maxFaceWarpage>=0 && faceWarpage.value>o.maxFaceWarpage) ||
      (o.maxChildVolumeRatio>=0 && childVolumeRatio.value>o.maxChildVolumeRatio) ;
  }
  std::string MeshQualityReport::description(const MeshQualityOptions &o) const {
    std::ostringstream out ;
    if(invalid()) { out<<"invalid or unresolved output geometry at root "<<invalidRoot<<", leaf "<<invalidLeaf<<", face "<<invalidFace ; return out.str() ; }
    const char *names[]={"nonOrthogonality","faceWarpage","childVolumeRatio"} ;
    const QualityValue values[]={nonOrthogonality,faceWarpage,childVolumeRatio} ;
    const double limits[]={o.maxNonOrthogonality,o.maxFaceWarpage,o.maxChildVolumeRatio} ;
    for(unsigned i=0;i<3;++i) if(values[i].root>=0) {
      if(out.tellp()>0) out<<"; " ;
      out<<names[i]<<'='<<values[i].value ;
      if(limits[i]>=0) out<<" (limit "<<limits[i]<<')' ;
      out<<" at root "<<values[i].root<<", leaf "<<values[i].leaf ;
      if(values[i].face>=0) out<<", face "<<values[i].face ;
    }
    return out.str() ;
  }
  QualityValue childVolumeRatio(const std::vector<std::vector<unsigned>> &paths,const std::vector<CellMeasure> &cells,int root) {
    require(paths.size()==cells.size(),"quality leaf paths differ from emitted cells") ;
    struct Volume { double sum=0 ; int leaf=-1 ; } ;
    std::map<std::vector<unsigned>,std::map<unsigned,Volume>> families ;
    for(unsigned i=0;i<paths.size();++i) {
      std::vector<unsigned> parent ;
      for(auto child:paths[i]) {
        auto &v=families[parent][child] ; v.sum+=cells[i].volume ; if(v.leaf<0) v.leaf=i ; parent.push_back(child) ;
      }
    }
    QualityValue result ;
    for(const auto &family:families) {
      require(family.second.size()>1,"incomplete quality split family") ;
      const Volume *small=nullptr,*large=nullptr ;
      for(const auto &entry:family.second) { if(!small || entry.second.sum<small->sum) small=&entry.second ; if(!large || entry.second.sum>large->sum) large=&entry.second ; }
      require(small->sum>0,"nonpositive child family volume") ;
      includeQualityValue(result,large->sum/small->sum,root,small->leaf) ;
    }
    return result ;
  }
}
}
