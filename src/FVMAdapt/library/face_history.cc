//#############################################################################
//#
//# Copyright 2026, Mississippi State University
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

#include "face_history.h"
#include "../face_history_internal.h"

#include <algorithm>
#include <limits>
#include <map>

namespace Loci {
  DatatypeP data_schema_traits<FaceSource>::get_type() {
    CompoundDatatypeP type = CompoundFactory(FaceSource()) ;
    LOCI_INSERT_TYPE(type,FaceSource,source) ;
    LOCI_INSERT_TYPE(type,FaceSource,orientation) ;
    LOCI_INSERT_TYPE(type,FaceSource,referenceFraction) ;
    return DatatypeP(type) ;
  }

  std::ostream& operator<<(std::ostream& out, const FaceSource& source) {
    return out << source.source << ' ' << source.orientation << ' ' << source.referenceFraction ;
  }

  std::istream& operator>>(std::istream& in, FaceSource& source) {
    return in >> source.source >> source.orientation >> source.referenceFraction ;
  }

  namespace face_history_detail {
    using Bounds = std::array<long long,4> ;

    bool intersects(const Bounds& a, const Bounds& b) {
      return a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3] ;
    }

    /// Bound the search on each refinement surface instead of comparing every face pair.
    class SurfaceIndex {
      struct Node {
        Bounds bounds ;
        size_t begin, end ;
        int left = -1, right = -1 ;
      } ;
      std::vector<size_t> faces ;
      std::vector<Node> nodes ;

      int build(const std::vector<HistoryFace>& source, size_t begin, size_t end) {
        Bounds bounds = source[faces[begin]].patch.bounds ;
        for(size_t i=begin+1; i<end; ++i) {
          const auto& b = source[faces[i]].patch.bounds ;
          for(int axis=0; axis<2; ++axis) {
            bounds[axis] = std::min(bounds[axis],b[axis]) ;
            bounds[axis+2] = std::max(bounds[axis+2],b[axis+2]) ;
          }
        }
        const int index = nodes.size() ;
        nodes.push_back({bounds,begin,end}) ;
        if(end-begin > 8) {
          const int axis = bounds[2]-bounds[0] >= bounds[3]-bounds[1] ? 0 : 1 ;
          const size_t middle = begin+(end-begin)/2 ;
          std::nth_element(faces.begin()+begin,faces.begin()+middle,faces.begin()+end,
            [&](size_t a, size_t b) {
              const auto& x = source[a].patch.bounds ;
              const auto& y = source[b].patch.bounds ;
              return x[axis]+(x[axis+2]-x[axis])/2 < y[axis]+(y[axis+2]-y[axis])/2 ;
            }) ;
          const int left = build(source,begin,middle), right = build(source,middle,end) ;
          nodes[index].left = left ; nodes[index].right = right ;
        }
        return index ;
      }

      void find(int index, const Bounds& bounds, std::vector<size_t>& result) const {
        const auto& node = nodes[index] ;
        if(!intersects(node.bounds,bounds)) return ;
        if(node.left < 0) {
          result.insert(result.end(),faces.begin()+node.begin,faces.begin()+node.end) ;
        } else {
          find(node.left,bounds,result) ; find(node.right,bounds,result) ;
        }
      }
    public:
      SurfaceIndex(const std::vector<HistoryFace>& source, std::vector<size_t> indices):
        faces(std::move(indices)) {
        if(!faces.empty()) build(source,0,faces.size()) ;
      }
      void find(const Bounds& bounds, std::vector<size_t>& result) const {
        if(!nodes.empty()) find(0,bounds,result) ;
      }
    } ;

    long double patchMeasure(const FacePatch& patch) {
      const auto& b = patch.bounds ;
      return static_cast<long double>(b[2]-b[0])*(b[3]-b[1])*
        (patch.sector < 0 ? -patch.sector : 1) ;
    }

    long double commonMeasure(const FacePatch& a, const FacePatch& b) {
      if(a.surface != b.surface) return 0 ;
      if(a.sector >= 0 && b.sector >= 0 && a.sector != b.sector) return 0 ;
      if(a.sector < 0 && b.sector >= -a.sector) return 0 ;
      if(b.sector < 0 && a.sector >= -b.sector) return 0 ;
      const long long x = std::min(a.bounds[2],b.bounds[2])-std::max(a.bounds[0],b.bounds[0]) ;
      const long long y = std::min(a.bounds[3],b.bounds[3])-std::max(a.bounds[1],b.bounds[1]) ;
      if(x <= 0 || y <= 0) return 0 ;
      const int sectors = a.sector < 0 && b.sector < 0 ? std::min(-a.sector,-b.sector) : 1 ;
      return static_cast<long double>(x)*y*sectors ;
    }

    bool validPatch(const FacePatch& patch) {
      return !patch.surface.empty() && patch.bounds[0] >= 0 && patch.bounds[1] >= 0 &&
        patch.bounds[2] > patch.bounds[0] && patch.bounds[3] > patch.bounds[1] &&
        patch.sector > std::numeric_limits<int>::min() &&
        patch.sector < std::numeric_limits<int>::max() &&
        (patch.orientation == 1 || patch.orientation == -1) ;
    }

    FaceSource contribution(const HistoryFace& old, const HistoryFace& next) {
      const auto& a = old.patch ;
      const auto& b = next.patch ;
      return {old.face,a.orientation*b.orientation,
        double(commonMeasure(a,b)/patchMeasure(a))} ;
    }

    /// Check rectangular coverage using integer endpoints, without area rounding.
    bool covers(const FacePatch& face, const std::vector<FacePatch>& parts,
                std::string* reason) {
      const auto fail = [reason](const char* message) {
        if(reason) *reason = message ;
        return false ;
      } ;
      const int first = face.sector < 0 ? 0 : face.sector ;
      const int last = face.sector < 0 ? -face.sector : face.sector+1 ;
      for(int sector=first; sector<last; ++sector) {
        struct Event { long long x, low, high ; int delta ; } ;
        std::vector<Event> events ;
        for(const FacePatch& part : parts) {
          if(part.sector >= 0 && part.sector != sector) continue ;
          if(part.sector < 0 && sector >= -part.sector) continue ;
          std::array<long long,4> b{{std::max(face.bounds[0],part.bounds[0]),
            std::max(face.bounds[1],part.bounds[1]),std::min(face.bounds[2],part.bounds[2]),
            std::min(face.bounds[3],part.bounds[3])}} ;
          if(b[2] <= b[0] || b[3] <= b[1]) continue ;
          events.push_back({b[0],b[1],b[3],1}) ;
          events.push_back({b[2],b[1],b[3],-1}) ;
        }
        std::sort(events.begin(),events.end(),[](const Event& a, const Event& b) {
          return a.x != b.x ? a.x < b.x : a.delta < b.delta ;
        }) ;
        std::map<long long,long long> active ;
        long long length = 0, last = face.bounds[0] ;
        for(size_t i=0; i<events.size();) {
          const long long x = events[i].x ;
          if(x != last && length != face.bounds[3]-face.bounds[1])
            return fail("incomplete patch coverage") ;
          last = x ;
          do {
            const auto& event = events[i++] ;
            if(event.delta < 0) {
              active.erase(event.low) ;
              length -= event.high-event.low ;
            } else {
              const auto next = active.lower_bound(event.low) ;
              if(next != active.end() && next->first < event.high)
                return fail("overlapping patches") ;
              if(next != active.begin() && std::prev(next)->second > event.low)
                return fail("overlapping patches") ;
              active.emplace(event.low,event.high) ;
              length += event.high-event.low ;
            }
          } while(i < events.size() && events[i].x == x) ;
        }
        if(last != face.bounds[2] || length != 0) return fail("incomplete patch coverage") ;
      }
      return true ;
    }

    void matchPatches(const std::vector<HistoryFace>& previous,
                     const std::vector<HistoryFace>& next,
                     const std::function<void(size_t,size_t)>& match) {
      std::map<std::vector<long long>,std::vector<size_t> > groups ;
      for(size_t i=0; i<previous.size(); ++i) groups[previous[i].patch.surface].push_back(i) ;
      std::map<std::vector<long long>,SurfaceIndex> search ;
      for(auto& group : groups)
        search.emplace(group.first,SurfaceIndex(previous,std::move(group.second))) ;
      for(size_t i=0; i<next.size(); ++i) {
        const auto surface = search.find(next[i].patch.surface) ;
        if(surface == search.end()) continue ;
        std::vector<size_t> candidates ;
        surface->second.find(next[i].patch.bounds,candidates) ;
        for(size_t j : candidates) {
          if(commonMeasure(previous[j].patch,next[i].patch) > 0) {
            match(j,i) ;
          }
        }
      }
    }
  }

}
