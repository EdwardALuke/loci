// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#include "../face_history_internal.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <set>

namespace Loci {
  namespace {
    using Word = long long ;
    using Buffer = std::vector<Word> ;
    using Messages = std::vector<Buffer> ;
    using Lookup = std::map<Word,Word> ;
    using namespace face_history_detail ;

    struct PrivateComm {
      MPI_Comm comm ;
      explicit PrivateComm(MPI_Comm input) { MPI_Comm_dup(input,&comm) ; }
      ~PrivateComm() { MPI_Comm_free(&comm) ; }
      PrivateComm(const PrivateComm&) = delete ;
      PrivateComm& operator=(const PrivateComm&) = delete ;
    } ;

    bool agree(bool valid, std::string& error, MPI_Comm comm) {
      int rank, bad ; MPI_Comm_rank(comm,&rank) ;
      int local = valid ? std::numeric_limits<int>::max() : rank ;
      MPI_Allreduce(&local,&bad,1,MPI_INT,MPI_MIN,comm) ;
      if(bad == std::numeric_limits<int>::max()) return true ;
      char message[256] = {} ;
      if(rank == bad) {
        if(error.empty()) error = "Invalid face-history input" ;
        std::strncpy(message,error.c_str(),sizeof(message)-1) ;
      }
      MPI_Bcast(message,sizeof(message),MPI_CHAR,bad,comm) ; error = message ;
      return false ;
    }

    /// Bound individual MPI messages without narrowing the full payload size.
    bool exchange(const Messages& send, Messages& receive, MPI_Comm comm,
                  std::string& error, FaceExchangeStats* stats = nullptr) {
      int p ; MPI_Comm_size(comm,&p) ;
      Buffer sizes(p), incoming(p) ;
      for(int i=0; i<p; ++i) sizes[i] = send[i].size() ;
      MPI_Alltoall(sizes.data(),1,MPI_LONG_LONG,incoming.data(),1,MPI_LONG_LONG,comm) ;
      bool valid = true ;
      try {
        receive.resize(p) ;
        for(int i=0; i<p; ++i) receive[i].resize(incoming[i]) ;
      } catch(const std::bad_alloc&) { valid = false ; error = "Face exchange allocation failed" ; }
      if(!agree(valid,error,comm)) return false ;
      const Word chunk = 65536 ;
      Word local = 0, rounds = 0 ;
      for(int i=0; i<p; ++i) {
        local = std::max(local,(sizes[i]+chunk-1)/chunk) ;
        if(stats) { stats->sentWords += sizes[i] ; stats->receivedWords += incoming[i] ; }
      }
      MPI_Allreduce(&local,&rounds,1,MPI_LONG_LONG,MPI_MAX,comm) ;
      for(Word round=0; round<rounds; ++round) {
        const Word offset = round*chunk ;
        std::vector<MPI_Request> requests ; requests.reserve(2*p) ;
        for(int i=0; i<p; ++i) {
          if(incoming[i] > offset) {
            requests.emplace_back() ;
            MPI_Irecv(receive[i].data()+offset,std::min(chunk,incoming[i]-offset),
                      MPI_LONG_LONG,i,0,comm,&requests.back()) ;
          }
          if(sizes[i] > offset) {
            requests.emplace_back() ;
            MPI_Isend(send[i].data()+offset,std::min(chunk,sizes[i]-offset),
                      MPI_LONG_LONG,i,0,comm,&requests.back()) ;
          }
        }
        MPI_Waitall(requests.size(),requests.data(),MPI_STATUSES_IGNORE) ;
      }
      return true ;
    }

    std::uint64_t mix(std::uint64_t value) {
      value ^= value >> 30 ; value *= UINT64_C(0xbf58476d1ce4e5b9) ;
      value ^= value >> 27 ; value *= UINT64_C(0x94d049bb133111eb) ;
      return value ^ (value >> 31) ;
    }

    Word surfaceNumber(const FacePatch& patch) {
      std::uint64_t value = 0 ;
      for(Word item : patch.surface) value = mix(value ^ mix(item)) ;
      return value & std::numeric_limits<Word>::max() ;
    }

    /// Find a face owner or a unique related cell (-1 means several cells).
    bool lookup(const std::vector<std::pair<Word,Word> >& entries,
                const std::set<Word>& needed, bool unique, Lookup& result,
                MPI_Comm comm, std::string& error, FaceExchangeStats* stats = nullptr) {
      int p ; MPI_Comm_size(comm,&p) ; Messages send(p), receive ;
      std::map<Word,Word> summary ;
      if(!unique) {
        for(const auto& entry : entries) {
          auto found = summary.emplace(entry.first,entry.second).first ;
          if(found->second != entry.second) found->second = -1 ;
        }
        for(const auto& entry : summary) {
          auto& out = send[mix(entry.first)%p] ;
          out.push_back(entry.first) ; out.push_back(entry.second) ;
        }
      } else {
        for(const auto& entry : entries) {
          auto& out = send[mix(entry.first)%p] ;
          out.push_back(entry.first) ; out.push_back(entry.second) ;
        }
      }
      if(!exchange(send,receive,comm,error,stats)) return false ;
      Lookup directory ; bool valid = true ;
      for(const auto& buffer : receive)
        for(size_t i=0; i<buffer.size(); i+=2) {
          const auto entry = directory.emplace(buffer[i],buffer[i+1]) ;
          if(!entry.second) {
            if(unique) valid = false ;
            else if(entry.first->second != buffer[i+1]) entry.first->second = -1 ;
          }
        }
      if(!valid) error = "Duplicate face identity" ;
      if(!agree(valid,error,comm)) return false ;
      send.assign(p,{}) ;
      for(Word key : needed) send[mix(key)%p].push_back(key) ;
      if(!exchange(send,receive,comm,error,stats)) return false ;
      send.assign(p,{}) ;
      for(int rank=0; rank<p; ++rank)
        for(Word key : receive[rank]) {
          const auto found = directory.find(key) ;
          if(found == directory.end()) { valid = false ; continue ; }
          send[rank].push_back(key) ; send[rank].push_back(found->second) ;
        }
      if(!valid) error = "Missing face owner or cell correspondence" ;
      if(!agree(valid,error,comm)) return false ;
      if(!exchange(send,receive,comm,error,stats)) return false ;
      result.clear() ;
      for(const auto& buffer : receive)
        for(size_t i=0; i<buffer.size(); i+=2) result[buffer[i]] = buffer[i+1] ;
      return true ;
    }

    /// Query cell pairs directly; do not import an entire descendant list.
    bool cellRelations(const std::vector<std::pair<Word,Word> >& entries,
                       const std::set<std::pair<Word,Word> >& needed,
                       std::set<std::pair<Word,Word> >& found,
                       MPI_Comm comm, std::string& error, FaceExchangeStats* stats) {
      int p ; MPI_Comm_size(comm,&p) ; Messages send(p), receive ;
      const auto owner = [p](const std::pair<Word,Word>& pair) {
        return mix(mix(pair.first)^std::uint64_t(pair.second))%p ;
      } ;
      for(const auto& pair : entries) {
        auto& buffer = send[owner(pair)] ; buffer.push_back(pair.first) ; buffer.push_back(pair.second) ;
      }
      if(!exchange(send,receive,comm,error,stats)) return false ;
      std::set<std::pair<Word,Word> > directory ; bool valid = true ;
      for(const auto& buffer : receive)
        for(size_t i=0; i<buffer.size(); i+=2)
          if(!directory.emplace(buffer[i],buffer[i+1]).second) valid = false ;
      if(!valid) error = "Duplicate cell correspondence" ;
      if(!agree(valid,error,comm)) return false ;
      send.assign(p,{}) ;
      for(const auto& pair : needed) {
        auto& buffer = send[owner(pair)] ; buffer.push_back(pair.first) ; buffer.push_back(pair.second) ;
      }
      if(!exchange(send,receive,comm,error,stats)) return false ;
      send.assign(p,{}) ;
      for(int rank=0; rank<p; ++rank)
        for(size_t i=0; i<receive[rank].size(); i+=2) {
          const std::pair<Word,Word> pair{receive[rank][i],receive[rank][i+1]} ;
          if(directory.count(pair)) { send[rank].push_back(pair.first) ; send[rank].push_back(pair.second) ; }
        }
      if(!exchange(send,receive,comm,error,stats)) return false ;
      for(const auto& buffer : receive)
        for(size_t i=0; i<buffer.size(); i+=2) found.emplace(buffer[i],buffer[i+1]) ;
      return true ;
    }

    struct Region {
      int first, last, axis = -1, left = -1, right = -1 ;
      Word cut = 0 ;
    } ;

    /// Split crowded refinement surfaces as well as unrelated base surfaces.
    class Regions {
      std::vector<Region> tree ;

      int locate(const std::array<Word,3>& point) const {
        int node = 0 ;
        while(tree[node].left >= 0) {
          const auto& split = tree[node] ;
          node = point[split.axis] <= split.cut ? split.left : split.right ;
        }
        return node ;
      }

      void visit(int node, const std::array<Word,3>& low, const std::array<Word,3>& high,
                 std::vector<int>& ranks) const {
        const auto& region = tree[node] ;
        if(region.left < 0) { ranks.push_back(region.first) ; return ; }
        if(low[region.axis] <= region.cut) visit(region.left,low,high,ranks) ;
        if(high[region.axis] > region.cut) visit(region.right,low,high,ranks) ;
      }
    public:
      void build(const std::vector<HistoryFace>& oldFaces,
                 const std::vector<HistoryFace>& newFaces, MPI_Comm comm) {
        int p ; MPI_Comm_size(comm,&p) ; tree.push_back({0,p}) ;
        if(p == 1) return ;
        std::vector<std::array<Word,3> > points ;
        for(const auto* faces : {&oldFaces,&newFaces})
          for(const auto& face : *faces) {
            const auto& b = face.patch.bounds ;
            points.push_back({{surfaceNumber(face.patch),b[0]+(b[2]-b[0]-1)/2,
                              b[1]+(b[3]-b[1]-1)/2}}) ;
          }
        std::vector<int> active{0} ;
        while(!active.empty()) {
          const size_t n = active.size() ;
          std::vector<int> slot(tree.size(),-1) ;
          for(size_t i=0; i<n; ++i) slot[active[i]] = i ;
          Buffer low(3*n,std::numeric_limits<Word>::max()), high(3*n,0), total(n,0) ;
          for(const auto& point : points) {
            const int i = slot[locate(point)] ; if(i < 0) continue ;
            ++total[i] ;
            for(int a=0; a<3; ++a) {
              low[3*i+a] = std::min(low[3*i+a],point[a]) ;
              high[3*i+a] = std::max(high[3*i+a],point[a]) ;
            }
          }
          MPI_Allreduce(MPI_IN_PLACE,low.data(),low.size(),MPI_LONG_LONG,MPI_MIN,comm) ;
          MPI_Allreduce(MPI_IN_PLACE,high.data(),high.size(),MPI_LONG_LONG,MPI_MAX,comm) ;
          MPI_Allreduce(MPI_IN_PLACE,total.data(),total.size(),MPI_LONG_LONG,MPI_SUM,comm) ;
          const int bins = 64 ; Buffer histogram(3*n*bins,0) ;
          for(const auto& point : points) {
            const int i = slot[locate(point)] ; if(i < 0) continue ;
            for(int a=0; a<3; ++a) {
              const auto span = static_cast<unsigned long long>(high[3*i+a]-low[3*i+a])+1 ;
              const int bin = static_cast<unsigned __int128>(point[a]-low[3*i+a])*bins/span ;
              ++histogram[(3*i+a)*bins+bin] ;
            }
          }
          MPI_Allreduce(MPI_IN_PLACE,histogram.data(),histogram.size(),MPI_LONG_LONG,MPI_SUM,comm) ;
          std::vector<int> next ;
          for(size_t i=0; i<n; ++i) {
            const int node = active[i] ;
            if(!total[i]) continue ;
            const int first = tree[node].first, last = tree[node].last, middle = (first+last)/2 ;
            const Word target = std::max<Word>(1,static_cast<unsigned __int128>(total[i])*(middle-first)/(last-first)) ;
            Word imbalance = std::numeric_limits<Word>::max(), cut = 0 ; int axis = -1 ;
            // Prefer separating surfaces on ties; otherwise split their crowded regions.
            for(int a=0; a<3; ++a) {
              if(low[3*i+a] == high[3*i+a]) continue ;
              const auto span = static_cast<unsigned long long>(high[3*i+a]-low[3*i+a])+1 ;
              Word sum = 0 ;
              for(int bin=0; bin<bins-1; ++bin) {
                sum += histogram[(3*i+a)*bins+bin] ;
                if(sum == 0 || sum == total[i]) continue ;
                const Word difference = sum > target ? sum-target : target-sum ;
                if(difference >= imbalance) continue ;
                imbalance = difference ; axis = a ;
                // ceil(span*(bin+1)/bins)-1 is the last integer in this bin.
                cut = low[3*i+a]+(static_cast<unsigned __int128>(span)*(bin+1)+bins-1)/bins-1 ;
              }
            }
            if(axis < 0) continue ;
            const int left = tree.size(), right = left+1 ;
            tree[node].axis = axis ; tree[node].cut = cut ;
            tree[node].left = left ; tree[node].right = right ;
            tree.push_back({first,middle}) ; tree.push_back({middle,last}) ;
            if(middle-first > 1) next.push_back(left) ;
            if(last-middle > 1) next.push_back(right) ;
          }
          active.swap(next) ;
        }
      }

      std::vector<int> owners(const FacePatch& patch) const {
        const Word surface = surfaceNumber(patch) ; const auto& b = patch.bounds ;
        std::vector<int> result ;
        visit(0,{{surface,b[0],b[1]}},{{surface,b[2]-1,b[3]-1}},result) ;
        return result ;
      }

      int owner(const FacePatch& a, const FacePatch& b) const {
        return tree[locate({{surfaceNumber(a),std::max(a.bounds[0],b.bounds[0]),
                            std::max(a.bounds[1],b.bounds[1])}})].first ;
      }
    } ;

    void pack(Buffer& out, const HistoryFace& face) {
      HistoryFaceConverter converter(const_cast<HistoryFace&>(face)) ;
      int size = converter.getSize() ;
      out.push_back(size) ;
      const size_t start = out.size() ;
      out.resize(start+size) ;
      converter.getState(out.data()+start,size) ;
    }

    HistoryFace unpack(const Buffer& buffer, size_t& i) {
      const Word size = buffer[i++] ;
      HistoryFace face ;
      HistoryFaceConverter(face).setState(buffer.data()+i,size) ;
      i += size ;
      return face ;
    }

  }

  bool buildFaceHistory(const std::vector<HistoryFace>& previous,
                        const std::vector<HistoryFace>& next,
                        const std::vector<std::pair<CellId,CellId> >& cell2parent,
                        FaceMatches& history, std::string& error, MPI_Comm comm,
                        FaceExchangeStats* stats) {
    history = FaceMatches() ;
    error.clear() ;
    PrivateComm privateComm(comm) ;
    comm = privateComm.comm ;
    if(stats) {
      *stats = FaceExchangeStats() ;
      stats->localFaces = previous.size()+next.size() ;
    }
    int rank, p ;
    MPI_Comm_rank(comm,&rank) ;
    MPI_Comm_size(comm,&p) ;
    bool valid = true ;
    const auto reject = [&](const char* mesh, FaceId face, const std::string& reason) {
      if(error.empty()) error = std::string(mesh)+" face "+std::to_string(face)+": "+reason ;
      valid = false ;
    } ;
    for(const auto* faces : {&previous,&next})
      for(const auto& face : *faces) valid = valid && validPatch(face.patch) && face.left >= 0 ;
    for(const auto& pair : cell2parent) valid = valid && pair.first >= 0 && pair.second >= 0 ;
    if(!agree(valid,error,comm)) return false ;
    Lookup ignored ;
    for(const auto* faces : {&previous,&next}) {
      std::vector<std::pair<Word,Word> > identities ;
      for(const auto& face : *faces) identities.emplace_back(face.face,rank) ;
      if(!lookup(identities,{},true,ignored,comm,error,stats)) return false ;
    }
    Lookup parents, children ; std::set<Word> needed ;
    for(const auto& face : next) { needed.insert(face.left) ; if(face.right >= 0) needed.insert(face.right) ; }
    if(!lookup(cell2parent,needed,false,parents,comm,error,stats)) return false ;
    std::vector<std::pair<Word,Word> > reverse ;
    for(const auto& pair : cell2parent) reverse.emplace_back(pair.second,pair.first) ;
    needed.clear() ;
    for(const auto& face : previous) { needed.insert(face.left) ; if(face.right >= 0) needed.insert(face.right) ; }
    if(!lookup(reverse,needed,false,children,comm,error,stats)) return false ;

    // Old and new meshes may have different owners. Match by refinement region,
    // then return contributions to the owning ranks using the Loci record converter.
    Regions regions ;
    regions.build(previous,next,comm) ;
    Messages send(p), receive ;
    for(int side=0; side<2; ++side) {
      const auto& faces = side ? next : previous ;
      for(size_t i=0; i<faces.size(); ++i)
        for(int owner : regions.owners(faces[i].patch)) {
          auto& out = send[owner] ;
          out.insert(out.end(),{side,rank,Word(i)}) ;
          pack(out,faces[i]) ;
        }
    }
    if(!exchange(send,receive,comm,error,stats)) return false ;
    std::vector<HistoryFace> faces[2] ;
    std::vector<std::pair<int,size_t> > locations[2] ;
    for(const auto& buffer : receive)
      for(size_t i=0; i<buffer.size();) {
        const int side = buffer[i++], owner = buffer[i++] ;
        const size_t index = buffer[i++] ;
        locations[side].emplace_back(owner,index) ;
        faces[side].push_back(unpack(buffer,i)) ;
      }
    receive.clear() ;
    send.assign(p,{}) ;
    if(stats) stats->matchingFaces = faces[0].size()+faces[1].size() ;
    matchPatches(faces[0],faces[1],[&](size_t previous, size_t next) {
      const auto& old = faces[0][previous] ;
      const auto& now = faces[1][next] ;
      if(regions.owner(old.patch,now.patch) != rank) return ;
      const auto a = locations[0][previous], b = locations[1][next] ;
      send[a.first].insert(send[a.first].end(),{0,Word(a.second)}) ;
      pack(send[a.first],now) ;
      send[b.first].insert(send[b.first].end(),{1,Word(b.second)}) ;
      pack(send[b.first],old) ;
      if(stats) ++stats->contributions ;
    }) ;
    for(int i=0; i<2; ++i) {
      std::vector<HistoryFace>().swap(faces[i]) ;
      std::vector<std::pair<int,size_t> >().swap(locations[i]) ;
    }
    if(!exchange(send,receive,comm,error,stats)) return false ;
    FaceMatches result ;
    result.sources.resize(next.size()) ;
    result.change.resize(next.size()) ;
    result.origin.assign(next.size(),-1) ;
    std::vector<std::vector<FacePatch> > coverage[2] ;
    coverage[0].resize(previous.size()) ;
    coverage[1].resize(next.size()) ;
    struct Neighbors {
      FaceId face ;
      CellId oldLeft, oldRight, newLeft, newRight ;
    } ;
    std::vector<Neighbors> neighbors ;
    std::set<std::pair<CellId,CellId> > neededPairs, foundPairs ;
    for(const auto& buffer : receive)
      for(size_t i=0; i<buffer.size();) {
        const int side = buffer[i++] ;
        const size_t index = buffer[i++] ;
        const auto other = unpack(buffer,i) ;
        coverage[side][index].push_back({{},other.patch.bounds,other.patch.sector,1}) ;
        if(side) {
          const auto& target = next[index] ;
          neighbors.push_back({target.face,other.left,other.right,target.left,target.right}) ;
          for(CellId current : {target.left,target.right})
            for(CellId previous : {other.left,other.right})
              if(current >= 0 && previous >= 0) neededPairs.emplace(current,previous) ;
          result.sources[index].push_back(contribution(other,target)) ;
        }
      }
    Messages().swap(send) ;
    Messages().swap(receive) ;
    if(!cellRelations(cell2parent,neededPairs,foundPairs,comm,error,stats)) return false ;
    for(const auto& cells : neighbors) {
      const CellId a = cells.oldLeft, b = cells.oldRight ;
      const CellId left = cells.newLeft, right = cells.newRight ;
      const bool adjacent = right < 0 ? b == right && foundPairs.count({left,a}) : b >= 0 &&
        ((foundPairs.count({left,a}) && foundPairs.count({right,b})) ||
         (foundPairs.count({left,b}) && foundPairs.count({right,a}))) ;
      if(!adjacent) reject("New",cells.face,"neighbors disagree with cell correspondence") ;
    }
    for(size_t i=0; i<next.size(); ++i) {
      auto& sources = result.sources[i] ;
      const auto& face = next[i] ;
      if(sources.empty()) {
        const auto& left = parents.at(face.left) ;
        const bool created = face.right >= 0 && left >= 0 && left == parents.at(face.right) ;
        if(!created) reject("New",face.face,"unmatched face is not inside one old cell") ;
        result.change[i] = FACE_CREATED ;
        if(created) result.origin[i] = left ;
      } else {
        std::string reason ;
        if(!covers(face.patch,coverage[1][i],&reason)) reject("New",face.face,reason) ;
        bool split = false ;
        for(const auto& patch : coverage[1][i])
          split = split || !covers(patch,{face.patch}) ;
        result.change[i] = sources.size() == 1 ? (split ? FACE_SPLIT : FACE_RETAINED) :
          (split ? FACE_RESHAPED : FACE_MERGED) ;
        std::sort(sources.begin(),sources.end(),[](const FaceSource& a, const FaceSource& b) {
          return a.source < b.source ;
        }) ;
      }
    }
    for(size_t i=0; i<previous.size(); ++i) {
      const auto& face = previous[i] ;
      if(coverage[0][i].empty()) {
        const auto& left = children.at(face.left) ;
        const bool removed = face.right >= 0 && left >= 0 && left == children.at(face.right) ;
        if(!removed) reject("Old",face.face,"unmatched face was not removed by coarsening") ;
        if(removed) {
          result.removed.push_back(face.face) ;
          result.removedInto.push_back(left) ;
        }
      } else {
        std::string reason ;
        if(!covers(face.patch,coverage[0][i],&reason)) reject("Old",face.face,reason) ;
      }
    }
    if(!agree(valid,error,comm)) return false ;
    history = std::move(result) ;
    return true ;
  }

  struct FaceDataExchange::Impl {
    MPI_Comm comm = MPI_COMM_NULL ;
    size_t owned = 0 ;
    std::vector<FaceId> sources ;
    std::vector<std::vector<size_t> > sendIndices, receiveIndices ;
    ~Impl() {
      int finalized ; MPI_Finalized(&finalized) ;
      if(!finalized && comm != MPI_COMM_NULL) MPI_Comm_free(&comm) ;
    }
  } ;

  bool FaceDataExchange::setup(const std::vector<FaceId>& owned, const std::vector<FaceId>& requested,
                               std::string& error, MPI_Comm comm) {
    impl.reset() ; error.clear() ;
    PrivateComm privateComm(comm) ; comm = privateComm.comm ;
    int rank, p ; MPI_Comm_rank(comm,&rank) ; MPI_Comm_size(comm,&p) ;
    std::vector<std::pair<Word,Word> > entries ; std::map<Word,size_t> local ;
    for(size_t i=0; i<owned.size(); ++i) { entries.emplace_back(owned[i],rank) ; local[owned[i]] = i ; }
    const std::set<Word> needed(requested.begin(),requested.end()) ;
    Lookup owners ;
    if(!lookup(entries,needed,true,owners,comm,error)) return false ;
    auto schedule = std::make_shared<Impl>() ; schedule->owned = owned.size() ;
    schedule->sources.assign(needed.begin(),needed.end()) ;
    schedule->sendIndices.resize(p) ; schedule->receiveIndices.resize(p) ;
    Messages send(p), receive ;
    for(size_t i=0; i<schedule->sources.size(); ++i) {
      const Word id = schedule->sources[i] ; const int owner = owners.at(id) ;
      send[owner].push_back(id) ; schedule->receiveIndices[owner].push_back(i) ;
    }
    if(!Loci::exchange(send,receive,comm,error)) return false ;
    for(int i=0; i<p; ++i) for(Word id : receive[i]) schedule->sendIndices[i].push_back(local.at(id)) ;
    MPI_Comm_dup(comm,&schedule->comm) ; impl = std::move(schedule) ; return true ;
  }

  const std::vector<FaceId>& FaceDataExchange::sources() const {
    static const std::vector<FaceId> empty ; return impl ? impl->sources : empty ;
  }


  bool FaceDataExchange::exchange(storeRepP input, storeRepP output,
                                  const entitySet& owned, int keyspace, int width,
                                  std::string& error) const {
    error.clear() ;
    if(!impl) { error = "Face data exchange has not been set up" ; return false ; }
    int widths[2] = {width,width} ;
    MPI_Allreduce(MPI_IN_PLACE,widths,1,MPI_INT,MPI_MIN,impl->comm) ;
    MPI_Allreduce(MPI_IN_PLACE,widths+1,1,MPI_INT,MPI_MAX,impl->comm) ;
    const bool valid = input.operator->() != output.operator->() && owned.size() == impl->owned &&
      (owned-input->domain()) == EMPTY && input->getDomainKeySpace() == keyspace &&
      width > 0 && widths[0] == widths[1] &&
      impl->sources.size() <= size_t(std::numeric_limits<int>::max()) ;
    if(!valid) error = "Old face field does not match its retained mesh or field width" ;
    if(!agree(valid,error,impl->comm)) return false ;

    // Donor values replace the output; its previous domain and width are unrelated.
    output->allocate(EMPTY) ;
    output->set_elem_size(width) ;
    output->allocate(impl->sources.empty() ? EMPTY : entitySet(interval(0,impl->sources.size()-1))) ;
    std::vector<Entity> entities ; entities.reserve(owned.size()) ;
    for(Entity entity : owned) entities.push_back(entity) ;
    const int p = impl->sendIndices.size() ;
    int peer = 0 ; size_t entry = 0 ;
    // Limit each packing round rather than retaining a second copy of the field.
    for(;;) {
      Messages send(p), receive ;
      size_t records = 0, packedBytes = 0 ;
      while(peer < p && records < 65536 && packedBytes < 1024*1024) {
        if(entry == impl->sendIndices[peer].size()) { ++peer ; entry = 0 ; continue ; }
        const Entity entity = entities[impl->sendIndices[peer][entry]] ;
        const entitySet domain = interval(entity,entity) ;
        int bytes = input->pack_size(domain), position = 0 ;
        auto& buffer = send[peer] ;
        buffer.push_back(entry++) ; buffer.push_back(bytes) ;
        const size_t start = buffer.size() ;
        buffer.resize(start+(size_t(bytes)+sizeof(Word)-1)/sizeof(Word)) ;
        input->pack(buffer.data()+start,position,bytes,domain) ;
        ++records ; packedBytes += 2*sizeof(Word)+size_t(bytes) ;
      }
      int active = records != 0 ;
      MPI_Allreduce(MPI_IN_PLACE,&active,1,MPI_INT,MPI_MAX,impl->comm) ;
      if(!active) break ;
      if(!Loci::exchange(send,receive,impl->comm,error)) return false ;
      for(int source=0; source<p; ++source) {
        auto& buffer = receive[source] ;
        for(size_t offset=0; offset<buffer.size();) {
          const size_t index = impl->receiveIndices[source][buffer[offset++]] ;
          int bytes = buffer[offset++], position = 0 ;
          output->unpack(buffer.data()+offset,position,bytes,sequence(interval(index,index))) ;
          offset += (size_t(bytes)+sizeof(Word)-1)/sizeof(Word) ;
        }
      }
    }
    return true ;
  }
}
