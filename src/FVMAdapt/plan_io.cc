// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
// Collective I/O for source geometry and accepted refinement histories.
#include <FVMAdapt/parametric_plan.h>
#include <FVMAdapt/parametric_cell.h>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <unistd.h>

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr int planFormatVersion=2 ;
    constexpr int maxPlanSettingsBytes=4096 ;

    void broadcast(std::string &value) {
      int size=int(value.size()) ; MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD) ;
      value.resize(size) ; MPI_Bcast(value.data(),size,MPI_CHAR,0,MPI_COMM_WORLD) ;
    }
    std::string encodePlanSettings(const RefinementOptions &options,int roots) {
      std::ostringstream out ; out<<std::setprecision(17) ;
      out<<"FVMAdaptPlan "<<planFormatVersion<<' '<<(MPI_processes>1)<<' '<<roots<<'\n'
        <<options.splitFamily<<' '<<options.centerMethod<<' '<<options.extrusionAxis<<'\n'
        <<options.grading<<' '<<options.faceBalance<<' '<<options.directionalBalance<<'\n'
        <<options.automaticSplit<<' '<<options.splitAspectRatio<<' '<<options.minEdgeLength<<'\n'
        <<options.maxDepth<<' '<<options.maxNodes<<' '<<options.levelDifference<<' '<<options.maxRounds<<' '<<options.refinementLevels<<'\n'
        <<options.quality.mode<<' '<<options.quality.centroid<<' '<<options.quality.maxNonOrthogonality<<' '<<options.quality.maxFaceWarpage<<' '<<options.quality.maxChildVolumeRatio<<'\n' ;
      return out.str() ;
    }
    bool decodePlanSettings(const std::string &text,RefinementOptions &options,int &roots,std::string &error) {
      std::istringstream in(text) ; std::string magic ; int version=0 ; bool parallel=false ;
      if(!(in>>magic>>version>>parallel>>roots) || magic!="FVMAdaptPlan" || (version!=1 && version!=planFormatVersion)) {
        error="not a supported parametric plan; regenerate old plans with marker" ; return false ;
      }
      if(parallel!=(MPI_processes>1)) { error="plans cannot cross between serial and parallel execution" ; return false ; }
      if(!(in>>options.splitFamily>>options.centerMethod>>options.extrusionAxis
          >>options.grading>>options.faceBalance>>options.directionalBalance
          >>options.automaticSplit>>options.splitAspectRatio>>options.minEdgeLength
          >>options.maxDepth>>options.maxNodes>>options.levelDifference>>options.maxRounds>>options.refinementLevels)
          || roots<1) { error="invalid parametric plan settings" ; return false ; }
      if(version==planFormatVersion && !(in>>options.quality.mode>>options.quality.centroid>>options.quality.maxNonOrthogonality>>options.quality.maxFaceWarpage>>options.quality.maxChildVolumeRatio)) {
        error="invalid parametric quality settings" ; return false ;
      }
      if(!(in>>std::ws).eof()) { error="unexpected parametric plan settings" ; return false ; }
      try { options.quality.validate() ; } catch(const std::exception &e) { error=e.what() ; return false ; }
      return true ;
    }
    bool writeSettings(hid_t file,const std::string &text) {
      const hsize_t size=text.size() ; const auto space=H5Screate_simple(1,&size,nullptr) ;
      const auto data=H5Dcreate(file,"settings",H5T_NATIVE_CHAR,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
      const bool ok=data>=0 && H5Dwrite(data,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,text.data())>=0 ;
      if(data>=0) H5Dclose(data) ;
      H5Sclose(space) ; return ok ;
    }
    std::string readSettings(const std::string &name) {
      std::string text ; hid_t file=-1 ;
      H5E_BEGIN_TRY { file=H5Fopen(name.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ; } H5E_END_TRY ;
      if(file<0) return text ;
      if(H5Lexists(file,"settings",H5P_DEFAULT)>0 && H5Lexists(file,"parametricHistory",H5P_DEFAULT)>0) {
        const auto data=H5Dopen(file,"settings",H5P_DEFAULT),space=H5Dget_space(data) ;
        const auto size=H5Sget_simple_extent_npoints(space) ;
        if(size>0 && size<=maxPlanSettingsBytes && H5Sget_simple_extent_ndims(space)==1) {
          text.resize(size) ;
          if(H5Dread(data,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,text.data())<0) text.clear() ;
        }
        H5Sclose(space) ; H5Dclose(data) ;
      }
      H5Fclose(file) ; return text ;
    }
  }
  bool writeParametricPlan(const std::string &name,const RefinementState &state,const RefinementOptions &options,std::string &error) {
    error.clear() ;
    if(state.histories==0) { error="cannot write an uninitialized parametric plan" ; return false ; }
    int count=state.histories->domain().size(),roots=0 ; MPI_Allreduce(&count,&roots,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    std::string temporary ;
    if(MPI_rank==0) {
      temporary=name+".XXXXXX" ; const int fd=mkstemp(temporary.data()) ;
      if(fd<0) { error="cannot create plan output: "+name ; temporary.clear() ; } else close(fd) ;
    }
    broadcast(temporary) ; broadcast(error) ; if(!error.empty()) return false ;
    const auto file=hdf5CreateFile(temporary.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
    const bool owner=MPI_rank==0 || use_parallel_io ; hid_t group=0 ; int bad=0,anyBad=0 ;
    if(owner) {
      bad=file<0 || !writeSettings(file,encodePlanSettings(options,roots)) ;
      if(!bad) { group=H5Gcreate(file,"parametricHistory",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ; bad=group<0 ; }
    }
    MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(!anyBad) write_store(group,state.histories,state.histories->domain(),0,MPI_COMM_WORLD) ;
    if(owner && group>0) H5Gclose(group) ;
    hdf5CloseFile(file) ;
    if(MPI_rank==0) {
      if(anyBad || std::rename(temporary.c_str(),name.c_str())!=0) { error="cannot write plan: "+name ; std::remove(temporary.c_str()) ; }
    }
    broadcast(error) ; return error.empty() ;
  }
  bool readParametricPlan(const std::string &name,const std::string &mesh,RefinementState &state,
      RefinementOptions &options,std::string &error) {
    error.clear() ; std::string text ;
    if(MPI_rank==0) text=readSettings(name) ;
    broadcast(text) ;
    RefinementOptions controls ; int roots=0 ;
    if(!decodePlanSettings(text,controls,roots,error)) return false ;
    const auto file=hdf5OpenFile(name.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ; hid_t group=0 ;
    const bool owner=MPI_rank==0 || use_parallel_io ;
    if(owner) group=H5Gopen(file,"parametricHistory",H5P_DEFAULT) ;
    store<ParametricRecord> history ; storeRepP saved=history.Rep() ; int offset=0 ;
    read_store(group,saved,offset,MPI_COMM_WORLD) ; saved->shift(offset) ;
    if(owner) H5Gclose(group) ;
    hdf5CloseFile(file) ;
    int count=saved->domain().size(),total=0 ; MPI_Allreduce(&count,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    if(total!=roots) { error="parametric plan root count differs from its history" ; return false ; }
    RefinementState next ; next.histories=saved ; next.source=mesh ; state=std::move(next) ; options=controls ; return true ;
  }
}
}
