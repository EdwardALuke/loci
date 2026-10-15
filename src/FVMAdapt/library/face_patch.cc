// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#include "../face_patch.h"
#include "diamondcell.h"
#include "hexcell.h"
#include "prism.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace Loci {
  namespace {
    // Keys identify a base face, or a root cell and its split surface.
    enum SurfaceKind { BASE_FACE = 0, CELL_FACE = 1 } ;
    enum CellKind { GENERAL_CELL = 0, HEX_CELL = 1, PRISM_CELL = 2 } ;
    enum PrismSurface { PRISM_END = 0, PRISM_SIDE = 1 } ;
    const long long PATH_END = -1 ; // Separates a child path from a local surface index.
    const long long extent = 1LL << MAXLEVEL ;
    typedef std::array<long long,2> Point ;
    typedef std::array<Point,4> Corners ;

    FacePatch whole(const std::vector<long long>& surface) {
      FacePatch result ;
      result.surface = surface ;
      result.bounds = {{0,0,extent,extent}} ;
      return result ;
    }

    Corners corners(const FacePatch& patch) {
      const auto& b = patch.bounds ;
      return {{{{b[0],b[1]}},{{b[2],b[1]}},{{b[2],b[3]}},{{b[0],b[3]}}}} ;
    }

    Point midpoint(Point a, Point b) {
      return {{a[0]+(b[0]-a[0])/2,a[1]+(b[1]-a[1])/2}} ;
    }

    void require_split(const FacePatch& patch) {
      if(patch.bounds[2]-patch.bounds[0] < 2 ||
         patch.bounds[3]-patch.bounds[1] < 2)
        throw std::runtime_error("Face history exceeds refinement-coordinate depth") ;
    }

    template<class FaceType>
    bool insert_patch(std::map<const FaceType*,FacePatch>& patches,
                      const FaceType* face, const FacePatch& patch) {
      const auto result = patches.emplace(face,patch) ;
      if(!result.second) {
        const auto& old = result.first->second ;
        if(old.surface != patch.surface || old.bounds != patch.bounds || old.sector != patch.sector)
          throw std::runtime_error("Conflicting refinement patches on a shared face") ;
      }
      return result.second ;
    }
  }

  void FacePatchMap::add(const Face* face, int faceNumber) {
    add(face,{BASE_FACE,faceNumber}) ;
  }

  void FacePatchMap::add(const Face* face, const std::vector<long long>& surface) {
    FacePatch patch = whole(surface) ;
    if(face->numEdge != 4) patch.sector = -face->numEdge ;
    polygon(face,patch,corners(patch)) ;
  }

  void FacePatchMap::polygon(const Face* face, const FacePatch& patch,
                            const Corners& points) {
    if(!insert_patch(faces,face,patch) || face->child == 0) return ;
    require_split(patch) ;
    for(int i=0; i<face->numEdge; ++i) {
      FacePatch part = patch ;
      Corners next ;
      if(patch.sector < 0) {
        part.sector = i ;
        next = corners(part) ;
      } else {
        next = {{points[i],midpoint(points[i],points[(i+1)%4]),
                 midpoint(points[0],points[2]),midpoint(points[(i+3)%4],points[i])}} ;
        part.bounds = {{std::min(next[0][0],next[2][0]),std::min(next[0][1],next[2][1]),
                        std::max(next[0][0],next[2][0]),std::max(next[0][1],next[2][1])}} ;
      }
      polygon(face->child[i],part,next) ;
    }
  }

  void FacePatchMap::add(const QuadFace* face, int faceNumber) {
    quad(face,whole({BASE_FACE,faceNumber})) ;
  }

  void FacePatchMap::quad(const QuadFace* face, const FacePatch& patch) {
    if(!insert_patch(quads,face,patch) || face->code == 0) return ;
    for(int axis=0; axis<2; ++axis) {
      QuadFace** children = axis == 0 ? face->childx : face->childy ;
      if(children == 0) continue ;
      if(patch.bounds[axis+2]-patch.bounds[axis] < 2)
        throw std::runtime_error("Face history exceeds refinement-coordinate depth") ;
      const long long middle = patch.bounds[axis]+(patch.bounds[axis+2]-patch.bounds[axis])/2 ;
      for(int side=0; side<2; ++side) {
        FacePatch part = patch ;
        part.bounds[axis+(side == 0 ? 2 : 0)] = middle ;
        quad(children[side],part) ;
      }
    }
  }

  void FacePatchMap::add(const Cell* cell, int root) {
    for(int i=0; i<cell->numFace; ++i)
      add(cell->face[i],{BASE_FACE,root,i}) ;
    if(cell->child == 0) return ;
    std::vector<const DiamondCell*> children(cell->numNode) ;
    for(int i=0; i<cell->numNode; ++i) children[i] = cell->child[i] ;
    general(children,{CELL_FACE,root,GENERAL_CELL}) ;
  }

  void FacePatchMap::general(const std::vector<const DiamondCell*>& children,
                            const std::vector<long long>& parent) {
    // Register this split's surfaces before visiting any deeper split.
    for(size_t i=0; i<children.size(); ++i) {
      for(int f=0; f<2*children[i]->getNfold(); ++f) {
        const Face* face = children[i]->getFace(f) ;
        if(faces.count(face)) continue ;
        auto surface = parent ;
        surface.insert(surface.end(),{PATH_END,static_cast<long long>(i),f}) ;
        add(face,surface) ;
      }
    }
    for(size_t i=0; i<children.size(); ++i) {
      const DiamondCell* child = children[i] ;
      if(child->getChildCell() == 0) continue ;
      std::vector<const DiamondCell*> descendants(2*child->getNfold()+2) ;
      for(size_t j=0; j<descendants.size(); ++j) descendants[j] = child->getChildCell(j) ;
      auto path = parent ; path.push_back(i) ;
      general(descendants,path) ;
    }
  }

  void FacePatchMap::add(const HexCell* cell, int root) {
    hex(cell,root,{{0,0,0,extent,extent,extent}}) ;
  }

  void FacePatchMap::hex(const HexCell* cell, int root, const std::array<long long,6>& box) {
    for(int f=0; f<6; ++f) {
      const int axis = f/2 ;
      const long long plane = box[axis+(f%2 == 0 ? 3 : 0)] ;
      const int u = axis == 0 ? 1 : 0, v = axis == 2 ? 1 : 2 ;
      FacePatch patch = whole({CELL_FACE,root,HEX_CELL,axis,plane}) ;
      patch.bounds = {{box[u],box[v],box[u+3],box[v+3]}} ;
      quad(cell->getFace(f),patch) ;
    }
    const int code = cell->getMySplitCode() ;
    for(int i=0; i<cell->numChildren(); ++i) {
      auto childBox = box ;
      int bit = 0 ;
      for(int axis=2; axis>=0; --axis) {
        if(!(code & (4 >> axis))) continue ;
        const long long middle = box[axis]+(box[axis+3]-box[axis])/2 ;
        if(middle == box[axis]) throw std::runtime_error("Cell history exceeds refinement-coordinate depth") ;
        childBox[axis+((i & (1 << bit++)) ? 0 : 3)] = middle ;
      }
      hex(cell->getChildCell(i),root,childBox) ;
    }
  }

  struct FacePatchMap::PrismRegion {
    struct Side {
      std::vector<long long> surface ;
      long long lower, upper ;
      bool reversed ;
    } ;
    FacePatch end ;
    Corners points ;
    std::vector<Side> sides ;
    std::vector<long long> path ;
    long long lower = 0, upper = extent ;
  } ;

  void FacePatchMap::add(const Prism* cell, int root) {
    PrismRegion region ;
    region.end = whole({}) ;
    if(cell->getNfold() != 4) region.end.sector = -cell->getNfold() ;
    region.points = corners(region.end) ;
    region.path = {CELL_FACE,root,PRISM_CELL,PRISM_SIDE} ;
    for(int i=0; i<cell->getNfold(); ++i)
      region.sides.push_back({{BASE_FACE,root,i},0,extent,false}) ;
    prism(cell,root,region) ;
  }

  void FacePatchMap::prism(const Prism* cell, int root, const PrismRegion& region) {
    // Axial position and end-face subdivision are independent coordinates.
    for(int f=0; f<2; ++f) {
      FacePatch end = region.end ;
      end.surface = {CELL_FACE,root,PRISM_CELL,PRISM_END,f == 0 ? region.lower : region.upper} ;
      polygon(cell->getEndFace(f),end,region.points) ;
    }
    for(size_t f=0; f<region.sides.size(); ++f) {
      const auto& edge = region.sides[f] ;
      FacePatch side = whole(edge.surface) ;
      side.bounds = {{edge.lower,region.lower,edge.upper,region.upper}} ;
      quad(cell->getSideFace(f),side) ;
    }
    const int code = cell->getMySplitCode(), count = cell->getNfold() ;
    const long long middle = region.lower+(region.upper-region.lower)/2 ;
    for(int i=0; i<cell->numChildren(); ++i) {
      PrismRegion child = region ;
      if(code & 2) {
        require_split(region.end) ;
        const int sector = i%count, previous = (sector+count-1)%count ;
        child.path.push_back(sector) ;
        if(region.end.sector < 0) {
          child.end.sector = sector ; child.points = corners(child.end) ;
        } else {
          const auto& points = region.points ;
          child.points = {{points[sector],midpoint(points[sector],points[(sector+1)%4]),
            midpoint(points[0],points[2]),midpoint(points[previous],points[sector])}} ;
          child.end.bounds = {{std::min(child.points[0][0],child.points[2][0]),
            std::min(child.points[0][1],child.points[2][1]),
            std::max(child.points[0][0],child.points[2][0]),
            std::max(child.points[0][1],child.points[2][1])}} ;
        }
        auto first = region.sides[sector], last = region.sides[previous] ;
        if(first.upper-first.lower < 2 || last.upper-last.lower < 2)
          throw std::runtime_error("Prism history exceeds refinement-coordinate depth") ;
        const long long firstMiddle = first.lower+(first.upper-first.lower)/2 ;
        const long long lastMiddle = last.lower+(last.upper-last.lower)/2 ;
        if(first.reversed) first.lower = firstMiddle ; else first.upper = firstMiddle ;
        if(last.reversed) last.upper = lastMiddle ; else last.lower = lastMiddle ;
        auto firstSurface = region.path, lastSurface = region.path ;
        firstSurface.insert(firstSurface.end(),{PATH_END,sector}) ;
        lastSurface.insert(lastSurface.end(),{PATH_END,previous}) ;
        child.sides = {first,{firstSurface,0,extent,true},{lastSurface,0,extent,false},last} ;
      }
      if(code & 1) {
        if(middle == region.lower)
          throw std::runtime_error("Prism history exceeds refinement-coordinate depth") ;
        const bool top = code == 1 ? i == 1 : i >= count ;
        if(top) child.lower = middle ; else child.upper = middle ;
      }
      prism(cell->getChildCell(i),root,child) ;
    }
  }

  void FacePatchConverter::getState(long long* buffer, int& size) const {
    size = getSize() ;
    std::copy(value.bounds.begin(),value.bounds.end(),buffer) ;
    buffer[4] = value.sector ; buffer[5] = value.orientation ;
    std::copy(value.surface.begin(),value.surface.end(),buffer+6) ;
  }

  void FacePatchConverter::setState(const long long* buffer, int size) {
    if(size < 6) throw std::runtime_error("Incomplete face patch") ;
    std::copy(buffer,buffer+4,value.bounds.begin()) ;
    value.sector = buffer[4] ; value.orientation = buffer[5] ;
    value.surface.assign(buffer+6,buffer+size) ;
  }

  int FacePatchesConverter::getSize() const {
    size_t size = 0 ;
    for(const FacePatch& patch : value) size += 7+patch.surface.size() ;
    if(size > size_t(std::numeric_limits<int>::max()))
      throw std::runtime_error("Face patches exceed container size") ;
    return size ;
  }

  void FacePatchesConverter::getState(long long* buffer, int& size) const {
    size = getSize() ;
    for(FacePatch& patch : value) {
      FacePatchConverter converter(patch) ;
      int count = converter.getSize() ;
      *buffer++ = count ; converter.getState(buffer,count) ; buffer += count ;
    }
  }

  void FacePatchesConverter::setState(const long long* buffer, int size) {
    value.clear() ;
    while(size > 0) {
      const long long count = *buffer++ ; --size ;
      if(count < 6 || count > size) throw std::runtime_error("Invalid face-patch length") ;
      value.emplace_back() ;
      FacePatchConverter(value.back()).setState(buffer,count) ;
      buffer += count ; size -= count ;
    }
  }

  void HistoryFaceConverter::getState(long long* buffer, int& size) const {
    buffer[0] = value.face ;
    buffer[1] = value.left ;
    buffer[2] = value.right ;
    FacePatchConverter(value.patch).getState(buffer+3,size) ;
    size += 3 ;
  }

  void HistoryFaceConverter::setState(const long long* buffer, int size) {
    if(size < 9) throw std::runtime_error("Incomplete face record") ;
    value.face = buffer[0] ;
    value.left = buffer[1] ;
    value.right = buffer[2] ;
    FacePatchConverter(value.patch).setState(buffer+3,size-3) ;
  }

  std::ostream& operator<<(std::ostream& out, const HistoryFace& face) {
    return out << face.face << ' ' << face.left << ' ' << face.right << ' ' << face.patch ;
  }

  std::istream& operator>>(std::istream& in, HistoryFace& face) {
    return in >> face.face >> face.left >> face.right >> face.patch ;
  }

  std::ostream& operator<<(std::ostream& out, const FacePatch& patch) {
    out << patch.surface.size() ;
    for(long long value : patch.surface) out << ' ' << value ;
    for(long long value : patch.bounds) out << ' ' << value ;
    return out << ' ' << patch.sector << ' ' << patch.orientation ;
  }

  std::istream& operator>>(std::istream& in, FacePatch& patch) {
    size_t size ;
    if(!(in >> size)) return in ;
    patch.surface.resize(size) ;
    for(long long& value : patch.surface) in >> value ;
    for(long long& value : patch.bounds) in >> value ;
    return in >> patch.sector >> patch.orientation ;
  }
}
