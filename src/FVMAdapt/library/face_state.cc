// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include <Loci.h>
#include <FVMAdapt/defines.h>
#include <FVMAdapt/dataxferDB.h>
#include <FVMAdapt/gridInterface.h>
#include "../face_history_internal.h"

namespace Loci {
  class FaceRecording {
  public:
    store<FacePatch> patches ;
  } ;
  namespace {
    bool collective(bool valid, std::string& error) {
      int result = valid ;
      MPI_Allreduce(MPI_IN_PLACE,&result,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      if(!result && error.empty()) error = "Face history input is invalid on another rank" ;
      return result ;
    }

    std::string baseName(const std::string& name) {
      char* path = realpath((name+".vog").c_str(),nullptr) ;
      if(!path) return name+".vog" ;
      std::string result(path) ; free(path) ; return result ;
    }

    std::array<unsigned long long,2> planStamp() {
      const_store<std::vector<char> > plan(DataXFER_DB.getItem("currentPlan")) ;
      std::array<unsigned long long,2> result{{0,0}} ;
      for(Entity cell : plan.domain()) {
        unsigned long long hash = 14695981039346656037ULL ;
        hash = (hash ^ static_cast<unsigned int>(cell))*1099511628211ULL ;
        for(unsigned char code : plan[cell]) hash = (hash ^ code)*1099511628211ULL ;
        result[0] ^= hash ; result[1] += hash ;
      }
      MPI_Allreduce(MPI_IN_PLACE,result.data(),1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,MPI_COMM_WORLD) ;
      MPI_Allreduce(MPI_IN_PLACE,result.data()+1,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD) ;
      return result ;
    }
  }

  struct FaceState::Impl {
    struct Face {
      FaceId id ;
      CellId left, right ;
      size_t surface ;
      std::array<long long,4> bounds ;
      int sector, orientation ;
    } ;
    std::string base ;
    int level = 0, keyspace = 0 ;
    CellId cellOffset = 0 ;
    std::array<unsigned long long,2> plan ;
    entitySet domain ;
    std::vector<std::vector<long long> > surfaces ;
    std::vector<Face> faces ;
    // Pending handoff only; installed states do not retain their ancestors.
    std::shared_ptr<const FaceState> previous ;
    std::vector<std::pair<CellId,CellId> > parents ;

    std::vector<HistoryFace> expand() const {
      std::vector<HistoryFace> result ; result.reserve(faces.size()) ;
      for(const auto& face : faces)
        result.push_back({face.id,face.left,face.right,
          {surfaces[face.surface],face.bounds,face.sector,face.orientation}}) ;
      return result ;
    }
  } ;

  struct FaceHistory::Impl {
    std::shared_ptr<const FaceState> previous ;
    multiStore<FaceSource> sources ;
    store<int> change ;
    store<CellId> origin ;
    std::vector<std::pair<FaceId,CellId> > removed ;
    FaceDataExchange exchange ;
  } ;

  struct FaceHistoryAccess {
    static void check(const RefinementOptions& options, const std::string& base, int level) {
      int presence[2] = {bool(options.previousFaceState),bool(options.previousFaceState)} ;
      MPI_Allreduce(MPI_IN_PLACE,presence,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      MPI_Allreduce(MPI_IN_PLACE,presence+1,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
      std::string error ;
      if(presence[0] != presence[1]) error = "Every rank must supply the previous face state" ;
      if(!collective(presence[0] == presence[1],error)) throw std::runtime_error(error) ;
      if(!presence[0]) return ;
      const auto& state = options.previousFaceState->impl ;
      bool valid = state && level == state->level && baseName(base) == state->base &&
        (level == 0 || DataXFER_DB.getItem("currentPlan") != 0) ;
      if(!valid) error = "Previous face state does not match the base mesh or adaptation level" ;
      if(!collective(valid,error)) throw std::runtime_error(error) ;
      valid = level == 0 || planStamp() == state->plan ;
      if(!valid) error = "Previous face state does not match the current refinement plan" ;
      if(!collective(valid,error)) throw std::runtime_error(error) ;
    }
    static store<FacePatch>& patches(refinedGridData& grid) {
      if(!grid.faceRecording) grid.faceRecording = std::make_shared<FaceRecording>() ;
      return grid.faceRecording->patches ;
    }
    static const store<FacePatch>& patches(const refinedGridData& grid) {
      return grid.faceRecording->patches ;
    }

    static std::shared_ptr<const FaceState> base(fact_db& facts, const std::string& name) {
      auto state = std::make_shared<FaceState::Impl>() ; state->base = baseName(name) ;
      const_multiMap f2n = facts.get_variable("face2node") ;
      const_Map cl = facts.get_variable("cl"), cr = facts.get_variable("cr"), ref = facts.get_variable("ref") ;
      const_store<std::string> tags = facts.get_variable("boundary_tags") ;
      constraint cells = facts.get_variable("geom_cells"), boundary = facts.get_variable("boundary_faces") ;
      const_store<vector3d<double> > pos = facts.get_variable("pos") ;
      state->domain = f2n.domain() ; state->keyspace = f2n.Rep()->getDomainKeySpace() ;
      dMap faceFile, cellFile ;
      if(MPI_processes == 1) {
        for(Entity face : state->domain) faceFile[face] = face ;
        for(Entity cell : *cells) cellFile[cell] = cell ;
      } else {
        auto dist = facts.get_distribute_info() ;
        faceFile = dist->g2fv[state->keyspace].Rep() ; cellFile = dist->g2fv[0].Rep() ;
        entitySet needed = MapRepP(cl.Rep())->image(state->domain)+
          MapRepP(cr.Rep())->image(state->domain-*boundary) ;
        storeRepP rep = cellFile.Rep() ;
        entitySet clones = needed-cellFile.domain() ;
        fill_clone(rep,clones,facts.get_init_ptn(0)) ;
      }
      long long nodeCount = pos.domain().size() ;
      MPI_Allreduce(MPI_IN_PLACE,&nodeCount,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD) ;
      int first[2] = {std::numeric_limits<int>::max(),std::numeric_limits<int>::max()} ;
      for(Entity face : state->domain) first[0] = std::min(first[0],faceFile[face]) ;
      for(Entity cell : *cells) first[1] = std::min(first[1],cellFile[cell]) ;
      MPI_Allreduce(MPI_IN_PLACE,first,2,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      state->cellOffset = first[1] ;
      const long long extent = 1LL << MAXLEVEL ;
      std::vector<HistoryFace> faces ; bool valid = true ; std::string error ;
      for(Entity face : state->domain) {
        CellId right ;
        if((*boundary).inSet(face)) {
          const auto& tag = tags[ref[face]] ;
          int id = 0 ; std::istringstream number(tag.size() >= 3 ? tag.substr(3) : "") ;
          valid = valid && tag.substr(0,3) == "BC_" && bool(number >> id) && id > 0 ;
          right = -id ;
        } else right = CellId(cellFile[cr[face]])-first[1] ;
        FacePatch patch{{0,nodeCount+FaceId(faceFile[face])-first[0]},{{0,0,extent,extent}},
                        f2n[face].size() == 4 ? 0 : -int(f2n[face].size()),1} ;
        faces.push_back({faceFile[face],CellId(cellFile[cl[face]])-first[1],right,std::move(patch)}) ;
      }
      if(!valid) error = "Base face history requires VOG boundary tags" ;
      if(!collective(valid,error)) throw std::runtime_error(error) ;
      record(*state,faces) ;
      auto result = std::shared_ptr<FaceState>(new FaceState) ; result->impl = state ;
      blackbox<FaceStateP> fact ; fact.set_entitySet(~EMPTY) ; *fact = result ;
      facts.create_fact("faceState",fact) ; return result ;
    }
    static void record(FaceState::Impl& state, const std::vector<HistoryFace>& faces) {
      std::map<std::vector<long long>,size_t> surfaces ;
      state.faces.reserve(faces.size()) ;
      for(const auto& face : faces) {
        auto entry = surfaces.emplace(face.patch.surface,surfaces.size()) ;
        if(entry.second) state.surfaces.push_back(face.patch.surface) ;
        state.faces.push_back({face.face,face.left,face.right,entry.first->second,
          face.patch.bounds,face.patch.sector,face.patch.orientation}) ;
      }
    }

    static void prepare(refinedGridData& grid, const std::string& base, int level,
                        const RefinementOptions& options, bool hasParents) {
      auto state = std::make_shared<FaceState::Impl>() ;
      state->base = baseName(base) ; state->level = level ;
      state->plan = planStamp() ;
      state->previous = options.previousFaceState ;
      state->domain = grid.local_faces[MPI_rank] ;
      state->cellOffset = std::numeric_limits<int>::max() ;
      for(const auto& cells : grid.local_cells)
        if(cells != EMPTY) state->cellOffset = std::min<CellId>(state->cellOffset,cells.Min()) ;

      std::string error ;
      int presence[2] = {bool(state->previous),bool(state->previous)} ;
      MPI_Allreduce(MPI_IN_PLACE,presence,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      MPI_Allreduce(MPI_IN_PLACE,presence+1,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
      bool valid = presence[0] == presence[1] && (!state->previous || (hasParents && state->previous->impl &&
        state->previous->impl->base == state->base && state->previous->impl->level+1 == level)) ;
      if(!valid) error = "Previous face state does not belong to this base mesh and adaptation level" ;
      if(!collective(valid,error)) throw std::runtime_error(error) ;

      std::vector<HistoryFace> faces ; faces.reserve(state->domain.size()) ;
      for(Entity face : state->domain)
        faces.push_back({face,grid.new_cl[face]-state->cellOffset,
          grid.new_cr[face] < 0 ? grid.new_cr[face] : grid.new_cr[face]-state->cellOffset,
          recordedFacePatches(grid)[face]}) ;
      record(*state,faces) ;
      if(state->previous) {
        const_store<std::pair<int,int> > pairs(DataXFER_DB.getItem("c2p")) ;
        int offsets[2] = {std::numeric_limits<int>::max(),std::numeric_limits<int>::max()} ;
        for(Entity pair : pairs.domain()) {
          offsets[0] = std::min(offsets[0],pairs[pair].first) ;
          offsets[1] = std::min(offsets[1],pairs[pair].second) ;
        }
        MPI_Allreduce(MPI_IN_PLACE,offsets,2,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
        for(Entity pair : pairs.domain())
          state->parents.emplace_back(CellId(pairs[pair].first)-offsets[0],
                                       CellId(pairs[pair].second)-offsets[1]) ;
      }
      auto result = std::shared_ptr<FaceState>(new FaceState) ; result->impl = state ;
      grid.faceState = result ;
    }

    static bool install(fact_db& facts, const refinedGridData& grid,
                        const Map& faceNumbers, std::string& error) {
      const auto& pending = *grid.faceState->impl ;
      const entitySet domain = faceNumbers.domain() ;
      const auto keyspace = faceNumbers.Rep()->getDomainKeySpace() ;

      // Move topology by face identity, independently of the new mesh partition.
      store<HistoryFace> records ;
      records.allocate(pending.domain) ;
      std::vector<FaceId> ids, requested ;
      auto entity = pending.domain.begin() ;
      for(const auto& face : pending.faces) {
        ids.push_back(face.id) ;
        records[*entity++] = {face.id,face.left,face.right,
          {pending.surfaces[face.surface],face.bounds,face.sector,face.orientation}} ;
      }
      for(Entity face : domain) requested.push_back(faceNumbers[face]) ;
      FaceDataExchange exchange ;
      if(!exchange.setup(ids,requested,error,MPI_COMM_WORLD)) return false ;
      store<HistoryFace> received ;
      if(!exchange.exchange(records.Rep(),received.Rep(),pending.domain,0,1,error)) return false ;
      records.allocate(EMPTY) ;

      std::vector<HistoryFace> next ; next.reserve(domain.size()) ;
      const auto& donors = exchange.sources() ;
      for(Entity face : domain) {
        const size_t index = std::lower_bound(donors.begin(),donors.end(),faceNumbers[face])-donors.begin() ;
        next.push_back(std::move(received[index])) ;
      }
      received.allocate(EMPTY) ;
      auto state = std::make_shared<FaceState::Impl>() ;
      state->base = pending.base ; state->level = pending.level ;
      state->plan = pending.plan ;
      state->cellOffset = pending.cellOffset ; state->domain = domain ; state->keyspace = keyspace ;
      record(*state,next) ;
      auto retained = std::shared_ptr<FaceState>(new FaceState) ; retained->impl = state ;

      std::shared_ptr<FaceHistory> history ;
      if(pending.previous) {
        const auto& previous = *pending.previous->impl ;
        FaceMatches matches ;
        if(!buildFaceHistory(previous.expand(),next,pending.parents,matches,error,MPI_COMM_WORLD)) return false ;
        auto data = std::make_shared<FaceHistory::Impl>() ; data->previous = pending.previous ;
        ids.clear() ;
        requested.clear() ;
        for(const auto& face : previous.faces) ids.push_back(face.id) ;
        for(const auto& row : matches.sources)
          for(const auto& source : row) requested.push_back(source.source) ;
        if(!data->exchange.setup(ids,requested,error,MPI_COMM_WORLD)) return false ;
        store<int> counts ;
        counts.allocate(domain) ;
        size_t row = 0 ;
        for(Entity face : domain) counts[face] = matches.sources[row++].size() ;
        data->sources.allocate(counts) ; data->change.allocate(domain) ; data->origin.allocate(domain) ;
        data->sources.Rep()->setDomainKeySpace(keyspace) ;
        data->change.Rep()->setDomainKeySpace(keyspace) ; data->origin.Rep()->setDomainKeySpace(keyspace) ;
        row = 0 ;
        for(Entity face : domain) {
          for(size_t i=0; i<matches.sources[row].size(); ++i) {
            data->sources[face][i] = matches.sources[row][i] ;
          }
          data->change[face] = matches.change[row] ;
          data->origin[face] = matches.origin[row] < 0 ? -1 : matches.origin[row]+previous.cellOffset ;
          ++row ;
        }
        for(size_t i=0; i<matches.removed.size(); ++i)
          data->removed.emplace_back(matches.removed[i],matches.removedInto[i]+pending.cellOffset) ;
        history = std::shared_ptr<FaceHistory>(new FaceHistory) ; history->impl = data ;
        facts.create_fact("faceSources",data->sources) ;
        facts.create_fact("faceChange",data->change) ; facts.create_fact("faceOriginCell",data->origin) ;
      }
      blackbox<std::shared_ptr<const FaceState> > stateFact ;
      stateFact.set_entitySet(~EMPTY) ; *stateFact = retained ; facts.create_fact("faceState",stateFact) ;
      if(history) {
        blackbox<std::shared_ptr<const FaceHistory> > historyFact ;
        historyFact.set_entitySet(~EMPTY) ; *historyFact = history ; facts.create_fact("faceHistory",historyFact) ;
      }
      return true ;
    }
  } ;

  void prepareFaceState(refinedGridData& grid, const std::string& base, int level,
                        const RefinementOptions& options, bool hasParents) {
    FaceHistoryAccess::prepare(grid,base,level,options,hasParents) ;
  }

  void checkFaceState(const RefinementOptions& options, const std::string& base, int level) {
    FaceHistoryAccess::check(options,base,level) ;
  }

  store<FacePatch>& recordedFacePatches(refinedGridData& grid) {
    return FaceHistoryAccess::patches(grid) ;
  }

  const store<FacePatch>& recordedFacePatches(const refinedGridData& grid) {
    return FaceHistoryAccess::patches(grid) ;
  }

  bool installFaceHistory(fact_db& facts, const refinedGridData& grid,
                          const Map& faceNumbers, std::string& error) {
    return FaceHistoryAccess::install(facts,grid,faceNumbers,error) ;
  }

  std::shared_ptr<const FaceState> getFaceState(fact_db& facts) {
    auto rep = facts.get_variable("faceState") ;
    if(rep == 0) return {} ;
    const_blackbox<std::shared_ptr<const FaceState> > state(rep) ; return *state ;
  }

  std::shared_ptr<const FaceState> recordBaseFaces(fact_db& facts, const std::string& casename) {
    return FaceHistoryAccess::base(facts,casename) ;
  }

  std::shared_ptr<const FaceHistory> getFaceHistory(fact_db& facts) {
    auto rep = facts.get_variable("faceHistory") ;
    if(rep == 0) return {} ;
    const_blackbox<std::shared_ptr<const FaceHistory> > history(rep) ; return *history ;
  }

  const_multiStore<FaceSource> FaceHistory::sources() const { return const_multiStore<FaceSource>(impl->sources.Rep()) ; }
  const_store<int> FaceHistory::change() const { return const_store<int>(impl->change.Rep()) ; }
  const_store<CellId> FaceHistory::origin() const { return const_store<CellId>(impl->origin.Rep()) ; }
  const std::vector<std::pair<FaceId,CellId> >& FaceHistory::removed() const { return impl->removed ; }

  int FaceHistory::donorIndex(FaceId source) const {
    const auto& donors = impl->exchange.sources() ;
    const auto found = std::lower_bound(donors.begin(),donors.end(),source) ;
    if(found == donors.end() || *found != source) throw std::out_of_range("Old face is not a local donor") ;
    return found-donors.begin() ;
  }

  bool FaceHistory::gather(const FaceState& previous, storeRepP input, storeRepP output,
                            int width, std::string& error) const {
    error.clear() ;
    const bool valid = impl && previous.impl && impl->previous->impl == previous.impl ;
    if(!valid) error = "Face field belongs to a different retained mesh" ;
    if(!collective(valid,error)) return false ;
    return impl->exchange.exchange(input,output,previous.impl->domain,previous.impl->keyspace,width,error) ;
  }
}
