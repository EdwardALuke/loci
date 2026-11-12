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
    // Retain complete rule outputs as inputs for this call's remaining queries.
    // Removing their producers lets Loci communicate them as distributed inputs.
    void keepQueryResults(rule_db &rules,fact_db &facts,const std::string &query) {
      ruleSet producers ;
      for(const auto &name:variableSet(expression::create(query))) {
        facts.make_extensional_fact(name) ; producers+=rules.rules_by_target(name) ;
      }
      rules.remove_rules(producers) ;
    }
    template<class T> void parameter(fact_db &facts,const std::string &name,const T &value) {
      param<T> p ; *p=value ; facts.create_fact(name,p) ;
    }
  }
  namespace {
  struct PlanOutput {
    bool historyOnly=false,transferData=true ;
    std::string tagFile,cellParents ;
    const RefinementState *parent=nullptr ; ///< Override donors on replay; otherwise use the input history.
  } ;
  bool adapt(CPTR<refinedGridData> &grid,rule_db &rules,RefinementState &state,
      storeRepP tags,const std::string &mesh,const RefinementOptions &options,const PlanOutput &outputFiles,std::string &error,MeshQualityReport *qualityReport=nullptr) {
    rule_db queryRules=rules ;
    const std::string casename=mesh.size()>4 && mesh.substr(mesh.size()-4)==".vog" ? mesh.substr(0,mesh.size()-4):mesh ;
    error.clear() ; if(qualityReport) *qualityReport=MeshQualityReport() ;
    try { options.quality.validate() ; } catch(const std::exception &e) { error=e.what() ; return false ; }
    if(state.histories!=0 && state.source!=mesh) { error="parametric state belongs to another source mesh" ; return false ; }
    if(options.tagLocation!="cells" && options.tagLocation!="nodes") { error="parametric tag location must be cells or nodes" ; return false ; }
    // The grid reader consumes its weight store. Keep accepted state available on rejection.
    store<int> weights ; storeRepP inputWeights=0 ;
    if(state.weights!=0) {
      const_store<int> saved=state.weights ; weights.allocate(saved.domain()) ;
      for(int cell:saved.domain()) weights[cell]=saved[cell] ;
      inputWeights=weights.Rep() ;
    }
    fact_db facts ;
    if(!setupFVMGridWithWeightInStore(facts,mesh,inputWeights)) { error="unable to read parametric source grid" ; return false ; }
    createLowerUpper(facts) ; createEdgesPar(facts) ;
    constraint geometry=facts.get_variable("geom_cells") ; const entitySet cells=*geometry ;
    const auto checkCount=[&](storeRepP history) {
      if(history==0) return true ;
      int local[2]={int(cells.size()),int(history->domain().size())},global[2]={0,0} ;
      MPI_Allreduce(local,global,2,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ; return global[0]==global[1] ;
    } ;
    if(!checkCount(state.histories) || (outputFiles.parent && !checkCount(outputFiles.parent->histories))) {
      error="parametric history root count differs from source grid" ; return false ;
    }
    blackbox<storeRepP> savedHistory ; savedHistory.set_entitySet(~EMPTY) ; *savedHistory=state.histories ;
    facts.create_fact("parametricHistoryInput",savedHistory) ;
    blackbox<storeRepP> cubicInput ; cubicInput.set_entitySet(~EMPTY) ; *cubicInput=options.cubicEdges ;
    facts.create_fact("parametricCubicGeometryInput",cubicInput) ;
    parameter(facts,"parametricSplitFamily",options.splitFamily) ;
    parameter(facts,"parametricRegionFile",options.regionFile) ; parameter(facts,"parametricSourceFile",options.sourceFile) ;
    parameter(facts,"parametricCenterMethod",options.centerMethod) ; parameter(facts,"parametricExtrusionAxis",options.extrusionAxis) ;
    parameter(facts,"parametricMaxDepth",options.maxDepth) ; parameter(facts,"parametricMaxNodes",options.maxNodes) ;
    parameter(facts,"parametricLevelDifference",options.levelDifference) ; parameter(facts,"parametricMaxRounds",options.maxRounds) ;
    parameter(facts,"parametricGrading",options.grading) ;
    parameter(facts,"parametricFaceBalance",options.faceBalance) ;
    parameter(facts,"parametricDirectionalBalance",options.directionalBalance) ;
    parameter(facts,"parametricAutomaticSplit",options.automaticSplit) ;
    parameter(facts,"parametricRefinementLevels",options.refinementLevels) ;
    parameter(facts,"parametricSplitAspectRatio",options.splitAspectRatio) ; parameter(facts,"parametricMinEdgeLength",options.minEdgeLength) ;
    parameter(facts,"parametricQualityMode",options.quality.mode) ; parameter(facts,"parametricQualityCentroid",options.quality.centroid) ;
    parameter(facts,"parametricMaxNonOrthogonality",options.quality.maxNonOrthogonality) ;
    parameter(facts,"parametricMaxFaceWarpage",options.quality.maxFaceWarpage) ;
    parameter(facts,"parametricMaxChildVolumeRatio",options.quality.maxChildVolumeRatio) ;
    parameter(facts,"recordFaceHistory",options.faceHistory) ;
    parameter(facts,"parametricHistoryReplay",tags==0 && outputFiles.tagFile.empty() && options.regionFile.empty() && options.sourceFile.empty() && state.histories!=0) ;
    const std::string initialResults="initialParametricCell,initialParametricHistory,parametricTopology,parametricInitialError" ;
    std::string initialQuery=initialResults ;
    const std::string inputName=options.tagLocation=="nodes" ? "inputTagsData":"inputCellTagsData" ;
    if(!outputFiles.tagFile.empty()) {
      parameter(facts,options.tagLocation=="nodes" ? "tagfile_par":"cell_tagfile_par",outputFiles.tagFile) ;
      initialQuery+=","+inputName ;
    }
    if(!makeQuery(queryRules,facts,initialQuery)) { error="parametric initialization query failed" ; return false ; }
    const_store<std::string> initialErrors=facts.get_variable("parametricInitialError") ; int bad=0,anyBad=0 ;
    for(int cell:cells) if(!initialErrors[cell].empty()) { bad=1 ; std::cerr<<initialErrors[cell]<<std::endl ; }
    MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    if(anyBad) { error="invalid parametric source geometry or history" ; return false ; }
    keepQueryResults(queryRules,facts,initialResults) ;
    if(!outputFiles.tagFile.empty()) {
      blackbox<std::vector<char>> input=facts.get_variable(inputName) ;
      store<char> marking ; marking.allocate(input->empty() ? EMPTY:entitySet(interval(0,int(input->size())-1))) ;
      for(int i:marking.domain()) marking[i]=(*input)[i] ;
      tags=marking.Rep() ;
    }
    const_store<ParametricCellData> initialized=facts.get_variable("initialParametricCell") ;
    int oldCount=0,tagCount=0,totalOld=0,totalTags=0 ;
    if(tags!=0 && options.tagLocation=="nodes") {
      if(!makeQuery(queryRules,facts,"parametricInnerNodeCount,num_original_nodes")) { error="parametric node count query failed" ; return false ; }
      const_param<int> inner=facts.get_variable("parametricInnerNodeCount"),originals=facts.get_variable("num_original_nodes") ;
      if(MPI_rank==0) oldCount=*inner+*originals ;
    } else for(int cell:cells) oldCount+=initialized[cell].leafCount() ;
    if(tags!=0) {
      store<char> marking=tags ; tagCount=marking.domain().size() ;
      for(int cell:marking.domain()) if(marking[cell]<0 || marking[cell]>2) bad=1 ;
    }
    MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    MPI_Allreduce(&oldCount,&totalOld,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ; MPI_Allreduce(&tagCount,&totalTags,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    if(anyBad || (tags!=0 && totalTags!=totalOld)) { error="invalid parametric tags: old "+options.tagLocation+"="+std::to_string(totalOld)+", tags="+std::to_string(totalTags)+", invalid value="+std::to_string(anyBad) ; return false ; }
    blackbox<std::vector<char>> input ; input.set_entitySet(~EMPTY) ;
    if(tags!=0) {
      store<char> marking=tags ; for(int cell:marking.domain()) input->push_back(marking[cell]) ;
      facts.create_fact(inputName,input) ;
    } else {
      store<ParametricRecord> requests ; requests.allocate(cells) ; facts.create_fact("parametricRequests",requests) ;
    }
    // Check acceptance and face geometry before outputs publish transfer data.
    std::string query="candidateParametricCell,parametricBalanced,parametricAccepted,parametricQualityReport,parametricError,balancedParametricCell,balancedParametricHistory,balanced_num_fine_cells,fineCellDepth" ;
    if(options.faceHistory) query+=",parametricCellPatchError,parametricFacePatchError,fine_face_patches_cell,fine_face_patches" ;
    if(!makeQuery(queryRules,facts,query)) { error="parametric acceptance query failed" ; return false ; }
    const_param<bool> accepted=facts.get_variable("parametricAccepted") ;
    const_param<MeshQualityReport> measured=facts.get_variable("parametricQualityReport") ;
    if(qualityReport) *qualityReport=*measured ;
    if(options.quality.mode!="off" && MPI_rank==0) std::cout<<"FVMAdapt quality: "<<measured->description(options.quality)<<std::endl ;
    if(!*accepted) {
      const_store<std::string> reasons=facts.get_variable("parametricError") ;
      for(int cell:cells) if(!reasons[cell].empty()) std::cerr<<"Parametric cell "<<cell<<": "<<reasons[cell]<<std::endl ;
      error="parametric adaptation rejected; accepted mesh retained" ; return false ;
    }
    if(options.faceHistory) {
      bad=0 ;
      for(const char *name:{"parametricCellPatchError","parametricFacePatchError"}) {
        const_store<std::string> errors=facts.get_variable(name) ;
        for(int entity:errors.domain()) if(!errors[entity].empty()) { bad=1 ; std::cerr<<errors[entity]<<std::endl ; }
      }
      MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
      if(anyBad) { error="parametric face patches rejected; accepted mesh retained" ; return false ; }
    }
    keepQueryResults(queryRules,facts,query) ;
    // The accepted results replace the complete balancing iteration.
    ruleSet balancing ;
    for(const auto &rule:queryRules.all_rules())
      if(rule.source_time().level_name()=="parametricBalance" || rule.target_time().level_name()=="parametricBalance") balancing+=rule ;
    queryRules.remove_rules(balancing) ;
    if(outputFiles.historyOnly) {
      if(!makeQuery(queryRules,facts,"parametricOutputState")) { error="parametric plan output query failed" ; return false ; }
      blackbox<RefinementState> next=facts.get_variable("parametricOutputState") ;
      state=*next ; state.source=mesh ; return true ;
    }
    if(!outputFiles.cellParents.empty()) {
      parameter(facts,"cell2parent_file_par",outputFiles.cellParents) ;
    }
    if(!outputFiles.cellParents.empty() && outputFiles.parent) {
      blackbox<storeRepP> parent ; parent.set_entitySet(~EMPTY) ;
      *parent=outputFiles.parent->histories ;
      facts.create_fact("parametricParentHistoryInput",parent) ;
      const std::string parentQuery="parametricDonorCell,parametricParentError" ;
      if(!makeQuery(queryRules,facts,parentQuery)) { error="parametric parent history query failed" ; return false ; }
      const_store<std::string> errors=facts.get_variable("parametricParentError") ; bad=0 ;
      for(int cell:cells) if(!errors[cell].empty()) { bad=1 ; std::cerr<<errors[cell]<<std::endl ; }
      MPI_Allreduce(&bad,&anyBad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
      if(anyBad) { error="invalid parametric parent history" ; return false ; }
      keepQueryResults(queryRules,facts,parentQuery) ;
    }
    blackbox<std::vector<std::pair<std::string,entitySet>>> volumeTags ; volumeTags.set_entitySet(~EMPTY) ;
    if(MPI_rank==0 || use_parallel_io) {
      const auto file=hdf5OpenFile(mesh.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ; readVolTags(file,*volumeTags) ; hdf5CloseFile(file) ;
    }
    facts.create_fact("volTagInput",volumeTags) ;
    std::string output="inner_nodes_cell,inner_nodes_face,inner_nodes_edge,fine_faces_cell,fine_faces,volTag_blackbox,parametricOutputState" ;
    if(outputFiles.transferData) output+=",cell2parent_DB" ;
    if(!outputFiles.cellParents.empty()) output+=",cell2parent_output" ;
    if(options.faceHistory) output+=",fine_face_patches_cell,fine_face_patches" ;
    if(options.cellState) output+=",fineCellDepth,fineCellChange,planRootFileNumber" ;
    if(options.edgeLengths) output+=",fineMaxEdgeLength,fineMaxEdgeLengthXY,balanced_num_fine_cells" ;
    if(options.cellState || options.edgeLengths) output+=",balanced_cell_offset" ;
    if(!makeQuery(queryRules,facts,output)) {
      error="parametric mesh output query failed" ; return false ;
    }
    CPTR<refinedGridData> candidate=new refinedGridData ;
    readBCfromVOG(mesh,candidate->boundary_ids) ;
    store<FineNodes> cellNodes=facts.get_variable("inner_nodes_cell"),faceNodes=facts.get_variable("inner_nodes_face"),edgeNodes=facts.get_variable("inner_nodes_edge") ;
    int nodeCount=0,faceCount=0,cellCount=0 ;
    createVOGNode(candidate->new_pos,cellNodes,faceNodes,edgeNodes,nodeCount,facts,candidate->local_nodes) ;
    store<FineFaces> cellFaces=facts.get_variable("fine_faces_cell"),faces=facts.get_variable("fine_faces") ;
    createVOGFace(nodeCount,cellFaces,faces,facts,faceCount,cellCount,candidate->new_cl,candidate->new_cr,
                  candidate->new_face2node,candidate->local_faces,candidate->local_cells,options.faceHistory ? &recordedFacePatches(*candidate):nullptr) ;
    blackbox<std::vector<std::pair<std::string,entitySet>>> outputTags=facts.get_variable("volTag_blackbox") ; candidate->volTags=*outputTags ;
    collectCellFields(*candidate,facts,options.cellState,options.edgeLengths) ;
    blackbox<RefinementState> next=facts.get_variable("parametricOutputState") ;
    candidate->hasFaceHistory=options.faceHistory ;
    if(options.faceHistory) prepareParametricFaceState(*candidate,casename,state.faces,state.histories,next->histories) ;
    if(outputFiles.transferData) {
      DataXFER_DB.insertItem("parametricHistory",next->histories) ; DataXFER_DB.insertItem("cellweights",next->weights) ;
    }
    state=*next ; state.source=mesh ; grid=candidate ;
    return true ;
  }
  }
  bool markParametricPlan(rule_db &rules,RefinementState &state,const std::string &mesh,
      const RefinementOptions &options,const std::string &tagFile,std::string &error) {
    PlanOutput output ; output.historyOnly=true ; output.tagFile=tagFile ;
    auto controls=options ; controls.faceHistory=false ;
    CPTR<refinedGridData> unused ; return adapt(unused,rules,state,0,mesh,controls,output,error) ;
  }
  bool writeParametricMesh(rule_db &rules,const RefinementState &state,const RefinementState *parent,
      const RefinementOptions &options,const std::string &output,const std::string &cellParents,std::string &error) {
    if(state.histories==0) { error="mesh reconstruction requires a saved plan" ; return false ; }
    if(parent && parent->source!=state.source) { error="parent plan belongs to another source mesh" ; return false ; }
    const RefinementState original ;
    PlanOutput files ; files.transferData=false ; files.parent=parent ? parent:&original ; files.cellParents=cellParents ;
    auto controls=options ; controls.regionFile.clear() ; controls.sourceFile.clear() ; controls.faceHistory=false ;
    auto replay=state ; CPTR<refinedGridData> grid ;
    if(!adapt(grid,rules,replay,0,state.source,controls,files,error)) return false ;
    writeVOG(output,grid->new_pos,grid->new_cl,grid->new_cr,grid->new_face2node,grid->boundary_ids,grid->volTags) ;
    return true ;
  }
  bool refineParametricFile(rule_db &rules,RefinementState &state,const std::string &mesh,const RefinementOptions &options,
      const std::string &tagFile,const std::string &output,const std::string &cellParents,std::string &error) {
    PlanOutput files ; files.transferData=false ; files.tagFile=tagFile ; files.cellParents=cellParents ;
    auto controls=options ; controls.faceHistory=false ; CPTR<refinedGridData> grid ;
    if(!adapt(grid,rules,state,0,mesh,controls,files,error)) return false ;
    writeVOG(output,grid->new_pos,grid->new_cl,grid->new_cr,grid->new_face2node,grid->boundary_ids,grid->volTags) ;
    return true ;
  }
}
  bool onlineRefineMesh(CPTR<refinedGridData> &grid,rule_db &rules,RefinementState &state,
      storeRepP tags,const std::string &casename,const RefinementOptions &options,std::string &error,FVMAdapt::MeshQualityReport *qualityReport) {
    return FVMAdapt::adapt(grid,rules,state,tags,casename+".vog",options,{},error,qualityReport) ;
  }
  bool initializeGridFromPlan(CPTR<refinedGridData> &grid,rule_db &rules,RefinementState &state,
      const std::string &casename,const RefinementOptions &options,std::string &error) {
    if(state.histories==0) { error="mesh reconstruction requires a saved plan" ; return false ; }
    auto candidate=state ; candidate.faces.reset() ;
    auto controls=options ; controls.regionFile.clear() ; controls.sourceFile.clear() ;
    // Restoring an accepted mesh does not apply proposal limits.
    if(controls.quality.mode=="reject") controls.quality.mode="report" ;
    if(!onlineRefineMesh(grid,rules,candidate,0,casename,controls,error)) return false ;
    state=std::move(candidate) ; return true ;
  }
  bool recordRefinementBaseFaces(fact_db &facts,rule_db &rules,RefinementState &state,
      const std::string &casename,const RefinementOptions &options,std::string &error) {
    if(!options.faceHistory || state.histories!=0) { error="base face recording requires a fresh face-enabled parametric state" ; return false ; }
    RefinementState candidate ; CPTR<refinedGridData> grid ;
    auto baseOptions=options ; baseOptions.regionFile.clear() ; baseOptions.sourceFile.clear() ;
    if(baseOptions.quality.mode=="reject") baseOptions.quality.mode="report" ;
    if(!onlineRefineMesh(grid,rules,candidate,0,casename,baseOptions,error)) return false ;
    candidate.faces=recordBaseParametricFaces(facts,*grid,casename,candidate.histories) ; state=std::move(candidate) ;
    DataXFER_DB.deleteItem("parametricHistory") ; DataXFER_DB.deleteItem("cellweights") ; DataXFER_DB.deleteItem("c2p") ;
    return true ;
  }
  void restoreRefinementFaceState(refinedGridData& grid,const std::string& casename,const RefinementState& state) {
    if(grid.hasFaceHistory) prepareParametricFaceState(grid,casename,{},0,state.histories) ;
  }
}
