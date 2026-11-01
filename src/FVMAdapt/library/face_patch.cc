// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/face_history_internal.h>
#include <algorithm>
#include <limits>
#include <cstring>
#include <stdexcept>

namespace Loci {
  void FacePatchConverter::getState(long long* buffer, int& size) const {
    size = getSize() ;
    std::copy(value.bounds.begin(),value.bounds.end(),buffer) ;
    buffer[4] = value.sector ; buffer[5] = value.orientation ;
    buffer[6] = value.surface.size() ; buffer[7] = value.polygon.size() ;
    static_assert(sizeof(double)==sizeof(long long),"face patch wire word size") ;
    std::memcpy(buffer+8,&value.area,sizeof(double)) ;
    buffer[9] = value.mappedArea ;
    std::copy(value.surface.begin(),value.surface.end(),buffer+10) ;
    buffer += 10+value.surface.size() ;
    for(const auto& p:value.polygon) for(double x:p) { std::memcpy(buffer++,&x,sizeof(double)) ; }
  }

  void FacePatchConverter::setState(const long long* buffer, int size) {
    if(size < 10 || buffer[6]<0 || buffer[7]<0 || buffer[6]>size-10 || buffer[7]>(size-10-buffer[6])/2 ||
       10+buffer[6]+2*buffer[7]!=size || (buffer[9]!=0 && buffer[9]!=1)) throw std::runtime_error("Incomplete face patch") ;
    std::copy(buffer,buffer+4,value.bounds.begin()) ;
    value.sector = buffer[4] ; value.orientation = buffer[5] ;
    value.mappedArea = buffer[9]!=0 ;
    value.surface.assign(buffer+10,buffer+10+buffer[6]) ; value.polygon.resize(buffer[7]) ;
    std::memcpy(&value.area,buffer+8,sizeof(double)) ; buffer+=10+buffer[6] ;
    for(auto& p:value.polygon) for(double& x:p) { std::memcpy(&x,buffer++,sizeof(double)) ; }
  }

  int FacePatchesConverter::getSize() const {
    size_t size = 0 ;
    for(FacePatch& patch : value) size += 1+FacePatchConverter(patch).getSize() ;
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
      if(count < 10 || count > size) throw std::runtime_error("Invalid face-patch length") ;
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
    if(size < 13) throw std::runtime_error("Incomplete face record") ;
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
    out << ' ' << patch.sector << ' ' << patch.orientation << ' ' << patch.area << ' ' << patch.mappedArea << ' ' << patch.polygon.size() ;
    for(const auto& p:patch.polygon) out << ' ' << p[0] << ' ' << p[1] ;
    return out ;
  }

  std::istream& operator>>(std::istream& in, FacePatch& patch) {
    size_t size ;
    if(!(in >> size)) return in ;
    patch.surface.resize(size) ;
    for(long long& value : patch.surface) in >> value ;
    for(long long& value : patch.bounds) in >> value ;
    in >> patch.sector >> patch.orientation >> patch.area >> patch.mappedArea >> size ; patch.polygon.resize(size) ;
    for(auto& p:patch.polygon) in >> p[0] >> p[1] ;
    return in ;
  }
}
