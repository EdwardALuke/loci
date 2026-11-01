// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/xml_region.h>
#include <cmath>
#include <cfenv>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>
#ifdef USE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

namespace Loci {
namespace FVMAdapt {
#ifdef USE_LIBXML2
  namespace {
    using Predicate=std::function<bool(const Point &)> ;
    using Transform=std::function<Point(Point)> ;
    void require(bool ok,const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    std::string name(xmlNode *node) { return reinterpret_cast<const char *>(node->name) ; }
    std::vector<xmlNode *> children(xmlNode *node) {
      std::vector<xmlNode *> result ;
      for(auto child=node->children;child;child=child->next) if(child->type==XML_ELEMENT_NODE) result.push_back(child) ;
      return result ;
    }
    std::string content(xmlNode *node) {
      const auto text=xmlNodeGetContent(node) ;
      const std::string value=text ? reinterpret_cast<const char *>(text):"" ; xmlFree(text) ;
      const auto first=value.find_first_not_of(" \t\r\n"),last=value.find_last_not_of(" \t\r\n") ;
      return first==std::string::npos ? "":value.substr(first,last-first+1) ;
    }
    std::map<std::string,double> parameters(xmlNode *node,std::map<std::string,double> values) {
      std::map<std::string,bool> seen ;
      for(auto child:children(node)) {
        const auto key=name(child) ;
        require(values.count(key) && !seen[key],"unknown or duplicate XML shape parameter") ; seen[key]=true ;
        require(children(child).empty(),"XML shape parameter must be a number") ;
        std::istringstream input(content(child)) ; double value ; std::string extra ;
        require(bool(input>>value) && !(input>>extra) && std::isfinite(value),"XML shape parameter must be finite") ; values[key]=value ;
      }
      return values ;
    }
    Predicate shape(xmlNode *node) {
      const auto elements=children(node) ; require(elements.size()==1,"XML shape requires one primitive") ;
      node=elements.front() ; const auto kind=name(node) ;
      if(kind=="box") {
        const auto p=parameters(node,{{"x1",0},{"y1",0},{"z1",0},{"x2",1},{"y2",1},{"z2",1}}) ;
        const Point lo{{p.at("x1"),p.at("y1"),p.at("z1")}},hi{{p.at("x2"),p.at("y2"),p.at("z2")}} ;
        for(unsigned d=0;d<3;++d) require(lo[d]<=hi[d],"XML box bounds are reversed") ;
        return [lo,hi](const Point &q) { for(unsigned d=0;d<3;++d) if(q[d]<lo[d] || q[d]>hi[d]) return false ; return true ; } ;
      }
      if(kind=="sphere" || kind=="cylinder" || kind=="cone") {
        std::map<std::string,double> defaults={{"x0",0},{"y0",0},{"r",1}} ;
        if(kind!="cylinder") defaults["z0"]=0 ;
        if(kind!="sphere") { defaults["z1"]=0 ; defaults["z2"]=1 ; }
        const auto p=parameters(node,defaults) ;
        const Point center{{p.at("x0"),p.at("y0"),kind=="cylinder" ? 0:p.at("z0")}} ; const double radius=p.at("r") ;
        require(radius>=0,"XML radius must be nonnegative") ;
        if(kind=="sphere") return [center,radius](const Point &q) {
          double square=0 ; for(unsigned d=0;d<3;++d) square+=(q[d]-center[d])*(q[d]-center[d]) ; return square<=radius*radius ;
        } ;
        const double lo=p.at("z1"),hi=p.at("z2") ; require(lo<=hi,"XML axial bounds are reversed") ;
        return [center,radius,lo,hi,kind](const Point &q) {
          const double x=q[0]-center[0],y=q[1]-center[1],z=q[2]-center[2],r=kind=="cone" ? radius*z:radius ;
          return q[2]>=lo && q[2]<=hi && x*x+y*y<=r*r ;
        } ;
      }
      for(unsigned d=0;d<3;++d) for(bool plus:{false,true}) {
        const auto axis=std::string(1,"xyz"[d]) ;
        if(kind==axis+(plus ? "_plus_plane":"_minus_plane")) {
          const double cut=parameters(node,{{axis+"1",0}}).at(axis+"1") ;
          return [d,plus,cut](const Point &q) { return plus ? q[d]>=cut:q[d]<=cut ; } ;
        }
      }
      throw std::invalid_argument("unknown XML region shape: "+kind) ;
    }
    std::vector<Transform> transforms(xmlNode *node) {
      std::vector<Transform> result ;
      for(auto child:children(node)) {
        const auto kind=name(child) ;
        if(kind=="translate" || kind=="scale") {
          const double initial=kind=="scale" ? 1:0 ;
          const auto p=parameters(child,{{"x0",initial},{"y0",initial},{"z0",initial}}) ;
          const Point vector{{p.at("x0"),p.at("y0"),p.at("z0")}} ;
          if(kind=="scale") for(auto value:vector) require(value!=0,"XML scale must be nonzero") ;
          result.push_back([vector,kind](Point q) { for(unsigned d=0;d<3;++d) q[d]=kind=="scale" ? q[d]/vector[d]:q[d]-vector[d] ; return q ; }) ;
        } else {
          require(kind=="rotateX" || kind=="rotateY" || kind=="rotateZ","unknown XML transform") ;
          const unsigned axis=kind.back()-'X',a=(axis+1)%3,b=(axis+2)%3 ;
          const double angle=parameters(child,{{"theta",0}}).at("theta")*std::acos(-1.)/180.,c=std::cos(angle),s=std::sin(angle) ;
          result.push_back([a,b,c,s](Point q) { const double x=q[a],y=q[b] ; q[a]=c*x+s*y ; q[b]=-s*x+c*y ; return q ; }) ;
        }
      }
      return result ;
    }
    Predicate object(xmlNode *node) {
      const auto elements=children(node) ;
      require(elements.size()==1 || elements.size()==2,"XML object requires an optional transform and one shape") ;
      std::vector<Transform> mapping ;
      if(elements.size()==2) { require(name(elements[0])=="transform","XML transform must precede shape") ; mapping=transforms(elements[0]) ; }
      require(name(elements.back())=="shape","XML object requires a shape") ; const auto primitive=shape(elements.back()) ;
      // XML transforms map physical points into the shape's coordinates, in listed order.
      return [mapping,primitive](Point q) { for(const auto &map:mapping) q=map(q) ; return primitive(q) ; } ;
    }
    Predicate region(xmlNode *node) {
      require(node && name(node)=="region","XML root must be region") ;
      const auto elements=children(node) ; std::size_t i=0 ;
      const auto operand=[&]() {
        bool invert=false ;
        while(i<elements.size() && name(elements[i])=="op" && content(elements[i])=="complement") { invert=!invert ; ++i ; }
        require(i<elements.size(),"XML region is missing an operand") ; const auto next=elements[i++] ;
        require(name(next)=="object" || name(next)=="region","XML region expects an object or nested region") ;
        const auto value=name(next)=="region" ? region(next):object(next) ;
        return Predicate([value,invert](const Point &q) { return value(q)!=invert ; }) ;
      } ;
      auto result=operand() ;
      while(i<elements.size()) {
        require(name(elements[i])=="op","XML region is missing an operator") ; const auto operation=content(elements[i++]) ;
        require(operation=="union" || operation=="intersection" || operation=="difference","unknown XML region operator") ;
        const auto left=result,right=operand() ;
        result=[left,right,operation](const Point &q) {
          if(operation=="union") return left(q) || right(q) ;
          if(operation=="intersection") return left(q) && right(q) ;
          return left(q) && !right(q) ;
        } ;
      }
      return result ;
    }
  }
  XMLRegion::XMLRegion(const std::string &filename) {
    // libxml2 initializes its NaN/Inf constants with floating-point traps masked.
    std::fenv_t environment ; std::feholdexcept(&environment) ;
    xmlInitParser() ; std::fesetenv(&environment) ;
    std::unique_ptr<xmlDoc,decltype(&xmlFreeDoc)> doc(xmlReadFile(filename.c_str(),nullptr,XML_PARSE_NONET),xmlFreeDoc) ;
    require(bool(doc),"unable to read XML region") ; contains_=region(xmlDocGetRootElement(doc.get())) ;
  }
  XMLRegion::XMLRegion(_xmlNode *root) : contains_(region(root)) {}
#else
  XMLRegion::XMLRegion(const std::string &) { throw std::runtime_error("XML region support requires libxml2") ; }
  XMLRegion::XMLRegion(_xmlNode *) { throw std::runtime_error("XML region support requires libxml2") ; }
#endif
}
}
