// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
// Adaptation queries and publication of accepted mesh state.
#include <FVMAdapt/parametric_plan.h>
#include "parametric_cell_data.h"
#include "cell_state_internal.h"
#include <FVMAdapt/defines.h>
#include <dataxferDB.h>
#include "face_history_helpers.h"
#include <iostream>

namespace Loci {
namespace FVMAdapt {
  namespace {
    /// Supply a named parameter as an input fact.
    template<class T> void setInputParameter(fact_db &facts,const std::string &name,const T &value) {
      param<T> p ;
      *p=value ;
      facts.create_fact(name,p) ;
    }

    /// Supply refinement controls as input facts.
    void installRefinementOptions(fact_db &facts,const RefinementOptions &options) {
      setInputParameter(facts,"parametricSplitFamily",options.splitFamily) ;
      setInputParameter(facts,"parametricRegionFile",options.regionFile) ;
      setInputParameter(facts,"parametricSourceFile",options.sourceFile) ;
      setInputParameter(facts,"parametricCenterMethod",options.centerMethod) ;
      setInputParameter(facts,"parametricExtrusionAxis",options.extrusionAxis) ;
      setInputParameter(facts,"parametricMaxDepth",options.maxDepth) ;
      setInputParameter(facts,"parametricMaxNodes",options.maxNodes) ;
      setInputParameter(facts,"parametricLevelDifference",options.levelDifference) ;
      setInputParameter(facts,"parametricMaxRounds",options.maxRounds) ;
      setInputParameter(facts,"parametricGrading",options.grading) ;
      setInputParameter(facts,"parametricFaceBalance",options.faceBalance) ;
      setInputParameter(facts,"parametricDirectionalBalance",options.directionalBalance) ;
      setInputParameter(facts,"parametricAutomaticSplit",options.automaticSplit) ;
      setInputParameter(facts,"parametricRefinementLevels",options.refinementLevels) ;
      setInputParameter(facts,"parametricSplitAspectRatio",options.splitAspectRatio) ;
      setInputParameter(facts,"parametricMinEdgeLength",options.minEdgeLength) ;
      setInputParameter(facts,"parametricQualityMode",options.quality.mode) ;
      setInputParameter(facts,"parametricQualityCentroid",options.quality.centroid) ;
      setInputParameter(facts,"parametricMaxNonOrthogonality",options.quality.maxNonOrthogonality) ;
      setInputParameter(facts,"parametricMaxFaceWarpage",options.quality.maxFaceWarpage) ;
      setInputParameter(facts,"parametricMaxChildVolumeRatio",options.quality.maxChildVolumeRatio) ;
      setInputParameter(facts,"recordFaceHistory",options.faceHistory) ;
    }

    /// Select mesh, history and donor outputs for an adaptation query.
    struct AdaptationOutput {
      bool historyOnly=false,transferData=true ;
      std::string tagFile,cellParents ;
      const RefinementState *parent=nullptr ; ///< Override donors on replay; otherwise use the input history.
    } ;

    /// Run adaptation and retain accepted results; rejection leaves grid and state unchanged.
    bool adapt(CPTR<refinedGridData> &grid,rule_db &rules,RefinementState &state,
        storeRepP tags,const std::string &mesh,const RefinementOptions &options,const AdaptationOutput &output,std::string &error,MeshQualityReport *qualityReport=nullptr) {
      const std::string casename=mesh.size()>4 && mesh.substr(mesh.size()-4)==".vog" ? mesh.substr(0,mesh.size()-4):mesh ;
      error.clear() ;
      if(qualityReport) *qualityReport=MeshQualityReport() ;
      try {
        options.quality.validate() ;
      } catch(const std::exception &e) {
        error=e.what() ;
        return false ;
      }
      if(state.histories!=0 && state.source!=mesh) {
        error="parametric state belongs to another source mesh" ;
        return false ;
      }
      if(options.tagLocation!="cells" && options.tagLocation!="nodes") {
        error="parametric tag location must be cells or nodes" ;
        return false ;
      }

      // The grid reader consumes its weight store. Keep accepted state available on rejection.
      store<int> weights ;
      storeRepP inputWeights=0 ;
      if(state.weights!=0) {
        const_store<int> saved=state.weights ;
        weights.allocate(saved.domain()) ;
        for(int cell:saved.domain()) weights[cell]=saved[cell] ;
        inputWeights=weights.Rep() ;
      }

      // Read the original mesh and restore accepted histories.
      fact_db facts ;
      if(!setupFVMGridWithWeightInStore(facts,mesh,inputWeights)) {
        error="unable to read parametric source grid" ;
        return false ;
      }
      createLowerUpper(facts) ;
      createEdgesPar(facts) ;
      constraint geometry=facts.get_variable("geom_cells") ;
      const entitySet cells=*geometry ;
      const auto checkCount=[&](storeRepP history) {
        if(history==0) return true ;
        int local[2]={int(cells.size()),int(history->domain().size())},global[2]={0,0} ;
        MPI_Allreduce(local,global,2,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
        return global[0]==global[1] ;
      } ;
      if(!checkCount(state.histories) || (output.parent && !checkCount(output.parent->histories))) {
        error="parametric history root count differs from source grid" ;
        return false ;
      }

      blackbox<storeRepP> savedHistory ;
      savedHistory.set_entitySet(~EMPTY) ;
      *savedHistory=state.histories ;
      facts.create_fact("parametricHistoryInput",savedHistory) ;
      blackbox<storeRepP> cubicInput ;
      cubicInput.set_entitySet(~EMPTY) ;
      *cubicInput=options.cubicEdges ;
      facts.create_fact("parametricCubicGeometryInput",cubicInput) ;
      installRefinementOptions(facts,options) ;
      setInputParameter(facts,"parametricHistoryReplay",tags==0 && output.tagFile.empty() && options.regionFile.empty() && options.sourceFile.empty() && state.histories!=0) ;
      const bool parentOverride=!output.cellParents.empty() && output.parent;
      if(!output.cellParents.empty())
        setInputParameter(facts,"cell2parent_file_par",output.cellParents) ;
      if(parentOverride) {
        blackbox<storeRepP> parent ;
        parent.set_entitySet(~EMPTY) ;
        *parent=output.parent->histories ;
        facts.create_fact("parametricParentHistoryInput",parent) ;
      }

      const std::string initialResults="initialParametricCell,initialParametricHistory,parametricTopology,parametricInitialError" ;
      std::string initialQuery=initialResults ;
      const std::string inputName=options.tagLocation=="nodes" ? "inputTagsData":"inputCellTagsData" ;
      if(!output.tagFile.empty()) {
        setInputParameter(facts,options.tagLocation=="nodes" ? "tagfile_par":"cell_tagfile_par",output.tagFile) ;
        initialQuery+=","+inputName ;
      }
      if(!makeQuery(rules,facts,initialQuery)) {
        error="parametric initialization query failed" ;
        return false ;
      }
      const_store<std::string> initialErrors=facts.get_variable("parametricInitialError") ;
      int bad=0,anyBad=0 ;
      for(int cell:cells) if(!initialErrors[cell].empty()) {
        bad=1 ;
        std::cerr<<initialErrors[cell]<<std::endl ;
      }
      MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
      if(anyBad) {
        error="invalid parametric source geometry or history" ;
        return false ;
      }

      // Explicit inputs let subsequent queries redistribute the validated source state.
      facts.create_fact("preparedParametricCell",facts.get_variable("initialParametricCell")) ;
      facts.create_fact("preparedParametricHistory",facts.get_variable("initialParametricHistory")) ;
      facts.create_fact("preparedParametricTopology",facts.get_variable("parametricTopology")) ;
      if(!output.tagFile.empty()) {
        blackbox<std::vector<char>> input=facts.get_variable(inputName) ;
        store<char> marking ;
        marking.allocate(input->empty() ? EMPTY:entitySet(interval(0,int(input->size())-1))) ;
        for(int i:marking.domain()) marking[i]=(*input)[i] ;
        tags=marking.Rep() ;
      }

      // Prepare replay donors from the validated source cells.
      if(parentOverride) {
        if(!makeQuery(rules,facts,"parametricDonorCell,parametricParentError")) {
          error="parametric parent history query failed" ;
          return false ;
        }
        const_store<std::string> parentErrors=facts.get_variable("parametricParentError") ;
        for(int cell:cells) if(!parentErrors[cell].empty()) {
          bad=1 ;
          std::cerr<<parentErrors[cell]<<std::endl ;
        }
        MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
        if(anyBad) {
          error="invalid parametric parent history" ;
          return false ;
        }
        facts.create_fact("preparedParametricDonorCell",facts.get_variable("parametricDonorCell")) ;
      }

      // Validate tags before constructing a candidate.
      const_store<ParametricCellData> initialized=facts.get_variable("preparedParametricCell") ;
      int oldCount=0,tagCount=0,totalOld=0,totalTags=0 ;
      if(tags!=0 && options.tagLocation=="nodes") {
        if(!makeQuery(rules,facts,"parametricInnerNodeCount,num_original_nodes")) {
          error="parametric node count query failed" ;
          return false ;
        }
        const_param<int> inner=facts.get_variable("parametricInnerNodeCount"),originals=facts.get_variable("num_original_nodes") ;
        if(MPI_rank==0) oldCount=*inner+*originals ;
      } else for(int cell:cells) oldCount+=initialized[cell].leafCount() ;
      if(tags!=0) {
        store<char> marking=tags ;
        tagCount=marking.domain().size() ;
        for(int cell:marking.domain()) if(marking[cell]<0 || marking[cell]>2) bad=1 ;
      }
      MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
      MPI_Allreduce(&oldCount,&totalOld,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
      MPI_Allreduce(&tagCount,&totalTags,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
      if(anyBad || (tags!=0 && totalTags!=totalOld)) {
        error="invalid parametric tags: old "+options.tagLocation+"="+std::to_string(totalOld)+", tags="+std::to_string(totalTags)+", invalid value="+std::to_string(anyBad) ;
        return false ;
      }

      blackbox<std::vector<char>> input ;
      input.set_entitySet(~EMPTY) ;
      if(tags!=0) {
        store<char> marking=tags ;
        for(int cell:marking.domain()) input->push_back(marking[cell]) ;
        facts.create_fact(inputName,input) ;
      } else {
        store<ParametricRecord> requests ;
        requests.allocate(cells) ;
        facts.create_fact("parametricRequests",requests) ;
      }

      std::string query="parametricAccepted,parametricQualityReport,parametricError,parametricOutputState" ;
      if(!output.historyOnly) {
        blackbox<std::vector<std::pair<std::string,entitySet>>> volumeTags ;
        volumeTags.set_entitySet(~EMPTY) ;
        if(MPI_rank==0 || use_parallel_io) {
          const auto file=hdf5OpenFile(mesh.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
          readVolTags(file,*volumeTags) ;
          hdf5CloseFile(file) ;
        }
        facts.create_fact("volTagInput",volumeTags) ;
        query+=",inner_nodes_cell,inner_nodes_face,inner_nodes_edge,fine_faces_cell,fine_faces,volTag_blackbox" ;
        if(output.transferData || options.faceHistory) query+=",cell2parentData" ;
        if(!output.cellParents.empty()) query+=",cell2parent_output" ;
        if(options.faceHistory)
          query+=",parametricCellPatchError,parametricFacePatchError,fine_face_patches_cell,fine_face_patches" ;
        if(options.cellState) query+=",fineCellDepth,fineCellChange,planRootFileNumber" ;
        if(options.edgeLengths) query+=",fineMaxEdgeLength,fineMaxEdgeLengthXY,balanced_num_fine_cells" ;
        if(options.cellState || options.edgeLengths) query+=",balanced_cell_offset" ;
      }

      // Output rules depend on acceptance; transfer data remains private to this query.
      if(!makeQuery(rules,facts,query)) {
        error="parametric adaptation query failed" ;
        return false ;
      }
      const_param<bool> accepted=facts.get_variable("parametricAccepted") ;
      const_param<MeshQualityReport> measured=facts.get_variable("parametricQualityReport") ;
      if(qualityReport) *qualityReport=*measured ;
      if(options.quality.mode!="off" && MPI_rank==0)
        std::cout<<"FVMAdapt quality: "<<measured->description(options.quality)<<std::endl ;
      if(!*accepted) {
        const_store<std::string> reasons=facts.get_variable("parametricError") ;
        for(int cell:cells) if(!reasons[cell].empty())
          std::cerr<<"Parametric cell "<<cell<<": "<<reasons[cell]<<std::endl ;
        error="parametric adaptation rejected; accepted mesh retained" ;
        return false ;
      }

      blackbox<RefinementState> next=facts.get_variable("parametricOutputState") ;
      if(output.historyOnly) {
        state=*next ;
        state.source=mesh ;
        return true ;
      }

      if(options.faceHistory) {
        bad=0 ;
        for(const char *name:{"parametricCellPatchError","parametricFacePatchError"}) {
          const_store<std::string> errors=facts.get_variable(name) ;
          for(int entity:errors.domain()) if(!errors[entity].empty()) {
            bad=1 ;
            std::cerr<<errors[entity]<<std::endl ;
          }
        }
        MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
        if(anyBad) {
          error="parametric face patches rejected; accepted mesh retained" ;
          return false ;
        }
      }

      storeRepP cellParents=0 ;
      if(output.transferData || options.faceHistory) {
        blackbox<storeRepP> pairs=facts.get_variable("cell2parentData") ;
        cellParents=*pairs ;
      }

      // Construct the accepted mesh before publishing state.
      CPTR<refinedGridData> candidate=new refinedGridData ;
      readBCfromVOG(mesh,candidate->boundary_ids) ;
      store<FineNodes> cellNodes=facts.get_variable("inner_nodes_cell"),faceNodes=facts.get_variable("inner_nodes_face"),edgeNodes=facts.get_variable("inner_nodes_edge") ;
      int nodeCount=0,faceCount=0,cellCount=0 ;
      createVOGNode(candidate->new_pos,cellNodes,faceNodes,edgeNodes,nodeCount,facts,candidate->local_nodes) ;
      store<FineFaces> cellFaces=facts.get_variable("fine_faces_cell"),faces=facts.get_variable("fine_faces") ;
      createVOGFace(nodeCount,cellFaces,faces,facts,faceCount,cellCount,candidate->new_cl,candidate->new_cr,
                    candidate->new_face2node,candidate->local_faces,candidate->local_cells,options.faceHistory ? &recordedFacePatches(*candidate):nullptr) ;
      blackbox<std::vector<std::pair<std::string,entitySet>>> outputTags=facts.get_variable("volTag_blackbox") ;
      candidate->volTags=*outputTags ;
      collectCellFields(*candidate,facts,options.cellState,options.edgeLengths) ;
      candidate->hasFaceHistory=options.faceHistory ;
      if(options.faceHistory) prepareParametricFaceState(*candidate,casename,state.faces,state.histories,next->histories,cellParents) ;

      if(output.transferData) {
        DataXFER_DB.insertItem("c2p",cellParents) ;
        DataXFER_DB.insertItem("parametricHistory",next->histories) ;
        DataXFER_DB.insertItem("cellweights",next->weights) ;
      }
      state=*next ;
      state.source=mesh ;
      grid=candidate ;
      return true ;
    }
  }

  /// Balance refinement requests and retain history without constructing the output mesh.
  bool markParametricPlan(rule_db &rules,RefinementState &state,const std::string &mesh,
      const RefinementOptions &options,const std::string &tagFile,std::string &error) {
    AdaptationOutput output ;
    output.historyOnly=true ;
    output.tagFile=tagFile ;
    auto controls=options ;
    controls.faceHistory=false ;
    CPTR<refinedGridData> unused ;
    return adapt(unused,rules,state,0,mesh,controls,output,error) ;
  }

  /// Reconstruct a saved refinement plan and write its mesh and optional donor pairs.
  bool writeParametricMesh(rule_db &rules,const RefinementState &state,const RefinementState *parent,
      const RefinementOptions &options,const std::string &output,const std::string &cellParents,std::string &error) {
    if(state.histories==0) {
      error="mesh reconstruction requires a saved plan" ;
      return false ;
    }
    if(parent && parent->source!=state.source) {
      error="parent plan belongs to another source mesh" ;
      return false ;
    }

    const RefinementState original ;
    AdaptationOutput files ;
    files.transferData=false ;
    files.parent=parent ? parent:&original ;
    files.cellParents=cellParents ;
    auto controls=options ;
    controls.regionFile.clear() ;
    controls.sourceFile.clear() ;
    controls.faceHistory=false ;

    auto replay=state ;
    CPTR<refinedGridData> grid ;
    if(!adapt(grid,rules,replay,0,state.source,controls,files,error)) return false ;
    writeVOG(output,grid->new_pos,grid->new_cl,grid->new_cr,grid->new_face2node,grid->boundary_ids,grid->volTags) ;
    return true ;
  }

  /// Apply file-based requests, write the adapted mesh and retain its history.
  bool refineParametricFile(rule_db &rules,RefinementState &state,const std::string &mesh,const RefinementOptions &options,
      const std::string &tagFile,const std::string &output,const std::string &cellParents,std::string &error) {
    AdaptationOutput files ;
    files.transferData=false ;
    files.tagFile=tagFile ;
    files.cellParents=cellParents ;
    auto controls=options ;
    controls.faceHistory=false ;
    CPTR<refinedGridData> grid ;
    if(!adapt(grid,rules,state,0,mesh,controls,files,error)) return false ;
    writeVOG(output,grid->new_pos,grid->new_cl,grid->new_cr,grid->new_face2node,grid->boundary_ids,grid->volTags) ;
    return true ;
  }
}

  /// Adapt the solver mesh and publish transfer data after acceptance.
  bool onlineRefineMesh(CPTR<refinedGridData> &grid,rule_db &rules,RefinementState &state,
      storeRepP tags,const std::string &casename,const RefinementOptions &options,std::string &error,FVMAdapt::MeshQualityReport *qualityReport) {
    return FVMAdapt::adapt(grid,rules,state,tags,casename+".vog",options,{},error,qualityReport) ;
  }

  /// Reconstruct an accepted mesh from saved history without applying new requests.
  bool initializeGridFromPlan(CPTR<refinedGridData> &grid,rule_db &rules,RefinementState &state,
      const std::string &casename,const RefinementOptions &options,std::string &error) {
    if(state.histories==0) {
      error="mesh reconstruction requires a saved plan" ;
      return false ;
    }

    auto candidate=state ;
    candidate.faces.reset() ;
    auto controls=options ;
    controls.regionFile.clear() ;
    controls.sourceFile.clear() ;

    // Restoring an accepted mesh does not apply proposal limits.
    if(controls.quality.mode=="reject") controls.quality.mode="report" ;
    if(!onlineRefineMesh(grid,rules,candidate,0,casename,controls,error)) return false ;
    state=std::move(candidate) ;
    return true ;
  }

  /// Record original face identities before the first adaptation.
  bool recordRefinementBaseFaces(fact_db &facts,rule_db &rules,RefinementState &state,
      const std::string &casename,const RefinementOptions &options,std::string &error) {
    if(!options.faceHistory || state.histories!=0) {
      error="base face recording requires a fresh face-enabled parametric state" ;
      return false ;
    }

    RefinementState candidate ;
    CPTR<refinedGridData> grid ;
    auto baseOptions=options ;
    baseOptions.regionFile.clear() ;
    baseOptions.sourceFile.clear() ;
    if(baseOptions.quality.mode=="reject") baseOptions.quality.mode="report" ;
    if(!onlineRefineMesh(grid,rules,candidate,0,casename,baseOptions,error)) return false ;
    candidate.faces=recordBaseParametricFaces(facts,*grid,casename,candidate.histories) ;
    state=std::move(candidate) ;

    DataXFER_DB.deleteItem("parametricHistory") ;
    DataXFER_DB.deleteItem("cellweights") ;
    DataXFER_DB.deleteItem("c2p") ;
    return true ;
  }

  /// Reattach recorded faces to the restored cell history.
  void restoreRefinementFaceState(refinedGridData& grid,const std::string& casename,const RefinementState& state) {
    if(grid.hasFaceHistory) prepareParametricFaceState(grid,casename,{},0,state.histories) ;
  }
}
