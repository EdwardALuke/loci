//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include <hdf5.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include "sciTypes.h"
#include "defines.h"
#include <LociGridReaders.h>
#include <Tools/tools.h>
#include <map> 
#include <distribute.h>
#include <distribute_container.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <sstream>

#include "FVMAdapt/defines.h"
#include "FVMAdapt/dataxferDB.h"
#include "FVMAdapt/gridInterface.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using Loci::storeRepP;
using Loci::constraint;
using std::vector;
using Loci::MPI_rank;
using Loci::MPI_processes;
using std::map ;
using std::ostringstream ;
using std::istringstream ;
///----------------------------------------------------------------------
// HACK WARNING
// 
// The FVMGridReader has been changed to support 4 gig entities.  This has
// an impact on how boundary faces are represented in the FVMadapt code.  We
// just copied the old code into here for now, but that means that adapt
// doesn't work with large interval sets.  This needs to be fixed... however
// we are reworking how the fact_db will be created so this can wait until that
// process is further along.

// RSM MOD 20181108
#ifdef LOCI_USE_METIS
#ifdef USE_SCOTCH
#include <scotch.h>
#endif
#include <metis.h>
#include <parmetis.h>

#ifdef USE_SCOTCH
typedef float metisreal_t ;
#else
#if REALTYPEWIDTH == 32
typedef float metisreal_t ;
#else
typedef double metisreal_t ;
#endif
#endif

#endif

namespace Loci {
  extern int metis_cpp_threshold ;
  
  storeRepP mapCellPartitionWeights(storeRepP wptr,
				    storeRepP c2p_ptr,
				    entitySet localcelldom) {
    if(wptr == 0)
      return wptr ;
    store<int> cell_weights_parent ;
    cell_weights_parent = wptr ;

    store<pair<int,int> > c2pset ;
    if(c2p_ptr == 0)
      return wptr ;
    c2pset = c2p_ptr ;
    entitySet domc2p = c2pset.domain() ;
    int r = Loci::MPI_rank ;
    int p = Loci::MPI_processes ;
     
    int minParent = std::numeric_limits<int>::max() ;
    int minTarget = std::numeric_limits<int>::max() ;
    FORALL(domc2p,ii) {
      minParent = min(minParent,c2pset[ii].second) ;
      minTarget = min(minTarget,c2pset[ii].first) ;
    } ENDFORALL ;
    
    int parentOffset = minParent ;
    int targetOffset = minTarget ;
    MPI_Allreduce(&minParent,&parentOffset,1,MPI_INT,MPI_MIN,
		  MPI_COMM_WORLD) ;
    MPI_Allreduce(&minTarget,&targetOffset,1,MPI_INT,MPI_MIN,
		  MPI_COMM_WORLD) ;
    int sz = domc2p.size() ;
    vector<pair<int,int> > p2c(sz) ;
    int cnt = 0 ;
    FORALL(domc2p,ii) {
      p2c[cnt].first = c2pset[ii].second-parentOffset ;
      p2c[cnt].second = c2pset[ii].first-targetOffset ;
      cnt++ ;;
    } ENDFORALL ;
    
    sort(p2c.begin(),p2c.end()) ;
    
    // Now distribute c2p to processors in parent file ordering
    entitySet cwdom = cell_weights_parent.domain() ;
    int cwsz = cwdom.size() ;
    vector<int> parentsizes(p) ;
    MPI_Allgather(&cwsz,1,MPI_INT,&parentsizes[0],1,MPI_INT,
		  MPI_COMM_WORLD) ;
    vector<int> parentoffsets(p+1,0) ;
    for(int i=0;i<p;++i)
      parentoffsets[i+1] = parentoffsets[i]+parentsizes[i] ;
    
    vector<pair<int,int> > splits(p-1) ;
    for(int i=0;i<p-1;++i)
      splits[i] = pair<int,int>(parentoffsets[i+1],-1) ;

    Loci::parSplitSort(p2c,splits,MPI_COMM_WORLD) ;

    int psz = p2c.size() ;

    vector<pair<int,int> > cell2weights(psz) ;
    entitySet::const_iterator ii = cwdom.begin() ;
    cnt = parentoffsets[r] ;
    int mincell = p2c[0].second ;
    int cnt2 = 0 ;
    for(int i=0;i<psz;) {
      while(i<psz && p2c[i].first == cnt) {
	mincell = min(mincell,p2c[i].second) ;
	cell2weights[cnt2] =pair<int,int>(p2c[i].second,
					  cell_weights_parent[*ii]) ;
	cnt2++ ;
	i++ ;
      }
      cnt++ ;
      ii++ ;
    }
    int mincell_global =mincell ;
    MPI_Allreduce(&mincell,&mincell_global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;

    int dmsz = localcelldom.size() ;
    vector<int> cellsizes(p) ;
    MPI_Allgather(&dmsz,1,MPI_INT,&cellsizes[0],1,MPI_INT,
		  MPI_COMM_WORLD) ;
    parentoffsets[0] = mincell_global ;
    for(int i=0;i<p;++i)
      parentoffsets[i+1] = parentoffsets[i]+cellsizes[i] ;
    
    for(int i=0;i<p-1;++i)
      splits[i] = pair<int,int>(parentoffsets[i+1],-1) ;


    sort(cell2weights.begin(),cell2weights.end()) ;
    Loci::parSplitSort(cell2weights,splits,MPI_COMM_WORLD) ;

    store<int> cellweightschild ;
    cellweightschild.allocate(localcelldom) ;
    cnt = 0;

    FORALL(localcelldom,ii) {
      cellweightschild[ii] = cell2weights[cnt].second ;
      cnt++ ;
    } ENDFORALL ;
    return cellweightschild.Rep() ;
  }


  void onlineRefineMesh(Loci::CPTR<refinedGridData> &gridDataP,
			rule_db &refmesh_rdb,
			int adaptmode,
			int level,
			storeRepP tags,
			string casename  ) {
    gridDataP = new(refinedGridData) ;
      
        
    Loci::debugout<< " running refinement in adaptcycle" << endl;
    
    // This is the name of the mesh file that you want to read in.
    string meshFile = casename+".vog";
           
    //otherwise, restart is false
    bool firstlevel = (level>0) ;
    // Setup the rule database.
    // Add all registered rules.  
    rule_db refine_rules;
    refine_rules = refmesh_rdb ;
      
    // Setup the refine_fact database.
    fact_db refine_facts;
      
    // Read in the mesh file.
    Loci::storeRepP cellwt = Loci::DataXFER_DB.getItem("cellweights") ;
    if(!Loci::setupFVMGridWithWeightInStore(refine_facts,meshFile, cellwt)) {
      std::cerr << "unable to read grid file '" << meshFile << "'" << std::endl ;
      Loci::Abort() ;
    }
    Loci::DataXFER_DB.deleteItem("cellweights") ;

      
    blackbox<vector<pair<string, entitySet> > > origVolTags ;
    origVolTags.Rep()->allocate(~EMPTY) ;
    if(MPI_rank==0){
      hid_t file_id = Loci::hdf5OpenFile(meshFile.c_str(),
					 H5F_ACC_RDONLY,H5P_DEFAULT) ;
    
      Loci::readVolTags(file_id,*origVolTags);
      Loci::hdf5CloseFile(file_id) ;
    }
    refine_facts.create_fact("volTagInput",origVolTags) ;
     
    //Setup Loci datastructures
    Loci::createLowerUpper(refine_facts) ;
    Loci::createEdgesPar(refine_facts) ;
    Loci::parallelClassifyCell(refine_facts);
       
      
    //this is a dummy parameter to trick Loci scheduler
    param<bool> beginWithMarker;
    *beginWithMarker = true;
    refine_facts.create_fact("beginWithMarker",beginWithMarker) ; 
      
    param<std::string> plan_outDB_par ;
    *plan_outDB_par = "nextPlan" ;
    refine_facts.create_fact("plan_outDB_par",plan_outDB_par) ;
  
    if(tags != 0) {
      store<char> tag_info ;
      tag_info = tags ;
      Loci::entitySet dom = tag_info.domain() ;
      vector<char> tagset(dom.size()) ;
      int cnt = 0 ;
      FORALL(dom,ii) {
	tagset[cnt] = tag_info[ii] ;
	cnt++ ;
      } ENDFORALL ;
      blackbox<vector<char> > inputTagsData ;
      inputTagsData.set_entitySet(~EMPTY) ;
      (*inputTagsData).swap(tagset) ;
      refine_facts.create_fact("inputCellTagsData",inputTagsData) ;
      Loci::DataXFER_DB.deleteItem("refineTag") ;
    }


    if(firstlevel){
      param<std::string> parent_planfile_par ;
      *parent_planfile_par = "currentPlan" ;
      refine_facts.create_fact("parent_planDB_par",parent_planfile_par) ;
    }
      
    gridDataP->boundary_ids.clear();
    Loci::readBCfromVOG(meshFile, gridDataP->boundary_ids);

      
                              
    //parameters to identify different options
    if(firstlevel){
      param<std::string> planfile_par ;
      *planfile_par = "currentPlan" ;
      refine_facts.create_fact("planDB_par",planfile_par) ;
    
      param<int> restart_tag_par;
      *restart_tag_par = 1;
      refine_facts.create_fact("restart_tag_par",restart_tag_par);
    }else{
      param<int> norestart_tag_par;
      *norestart_tag_par = 1;
      refine_facts.create_fact("norestart_tag_par",norestart_tag_par);
    }
  
    param<int> split_mode_par;
    *split_mode_par = adaptmode;
    refine_facts.create_fact("split_mode_par", split_mode_par);
        
    {
      param<std::string> cellweight_outDB_par;
      *cellweight_outDB_par = "cellweights";
      refine_facts.create_fact("cellweight_outDB_par",cellweight_outDB_par);
        
      if(!Loci::makeQuery(refine_rules,refine_facts,
			  "cellplan_output,cellweight_output,cell2parent_DB,inner_nodes,fine_faces,volTag_blackbox")) {
	std::cerr << "adapt query failed!" << std::endl;
	Loci::Abort();
      }
    }
      
    Loci::DataXFER_DB.deleteItem("currentPlan") ;
    Loci::storeRepP newPlan = Loci::DataXFER_DB.getItem("nextPlan") ;
    Loci::DataXFER_DB.insertItem("currentPlan",newPlan) ;
    Loci::DataXFER_DB.deleteItem("nextPlan") ;
	
    store<Loci::FineNodes> inner_nodes;
    inner_nodes = refine_facts.get_variable("inner_nodes");
        
       
    int num_nodes;
    createVOGNode(gridDataP->new_pos,
		  inner_nodes,
		  num_nodes,
		  refine_facts,
		  gridDataP->local_nodes
		  );
    if(Loci::MPI_rank ==0)
      cerr<< "num_nodes: " << num_nodes << " in adaptcycle" << endl;
        
    store<Loci::FineFaces> fine_faces;
    fine_faces = refine_facts.get_variable("fine_faces");
        
    int num_faces = 0;
    int num_cells = 0;
      
    createVOGFace( num_nodes,
		   fine_faces,
		   refine_facts,
		   num_faces,
		   num_cells,
		   gridDataP->new_cl,
		   gridDataP->new_cr,
		   gridDataP->new_face2node,
		   gridDataP->local_faces,
		   gridDataP->local_cells
		   );

    //update volume tags
    {
      blackbox<std::vector<std::pair<std::string, entitySet> > > volTag_blackbox;
      volTag_blackbox = refine_facts.get_variable("volTag_blackbox"); 
      gridDataP->volTags =  *volTag_blackbox;
    }

    if(Loci::MPI_rank==0){
      cerr<< "num_faces: " << num_faces << " in adaptcycle" <<endl;
      cerr<< "num_cells: " << num_cells << " in adaptcycle" << endl;
    }
  }

  void initializeGridFromPlan(Loci::CPTR<refinedGridData> &gridDataP,
			      int &level,
			      rule_db &refmesh_rdb,
			      string casename, string weightfile,
			      string restartplanfile) {
    string meshFile = casename+".vog";
    // Setup the rule database.
    // Add all registered rules.  
    rule_db refine_rules;
    refine_rules = refmesh_rdb ;
	  
    // Setup the refine_fact database.
    fact_db refine_facts;
	  
    if(Loci::MPI_rank == 0) 
      cout <<"reading in meshfile: " 
	   << meshFile << " for chem run" << std::endl;
	    
    // Read in the mesh file.
    if(!Loci::setupFVMGridWithWeightInFile(refine_facts,meshFile, weightfile)) {
      std::cerr << "unable to read mesh file '" << meshFile << "'" << std::endl ;
      Loci::Abort() ;
    }
          
    //Setup Loci datastructures
    Loci::createLowerUpper(refine_facts) ;
    Loci::createEdgesPar(refine_facts) ;
    Loci::parallelClassifyCell(refine_facts);
    blackbox<vector<pair<string, entitySet> > > origVolTags ;
    origVolTags.Rep()->allocate(~EMPTY) ;
    if(MPI_rank==0){
      hid_t file_id = Loci::hdf5OpenFile(meshFile.c_str(),
					 H5F_ACC_RDONLY,H5P_DEFAULT) ;
	    
      Loci::readVolTags(file_id,*origVolTags);
      Loci::hdf5CloseFile(file_id) ;
    }
    refine_facts.create_fact("volTagInput",origVolTags) ;
	  
    if(Loci::MPI_rank == 0) 
      cout <<"reading in cellplan file: " 
	   << restartplanfile << std::endl;
	  
    hid_t file_id = Loci::hdf5OpenFile(restartplanfile.c_str(),
				       H5F_ACC_RDONLY,H5P_DEFAULT) ;
	  
    store<std::vector<char> > cellPlan ;
    param<int> refineLevel ;
    hid_t group_id = 0 ;
    if(Loci::MPI_rank == 0)
      group_id = H5Gopen(file_id,"refineLevel",H5P_DEFAULT) ;
	    
    Loci::read_parameter(group_id,refineLevel.Rep(),MPI_COMM_WORLD) ;
    if(Loci::MPI_rank == 0) {
      H5Gclose(group_id) ;
      group_id = H5Gopen(file_id,"cellPlan",H5P_DEFAULT) ;
      cout << "restart refine level = " << *refineLevel << endl ;
    }
    level = *refineLevel ;
    int offset = 0; 
    Loci::read_store(group_id,cellPlan.Rep(),offset,MPI_COMM_WORLD) ;//read in file numbering store, simplePartition
    cellPlan.Rep()->shift(offset) ;
    if(Loci::MPI_rank == 0) {
      H5Gclose(group_id) ;
    }
    Loci::hdf5CloseFile(file_id) ;

    //finish reading currentPlan
    //inserting currentPlan into DataXFER_DB
    Loci::DataXFER_DB.insertItem("currentPlan",cellPlan.Rep()) ;
    param<std::string> balanced_planDB_par ;
    *balanced_planDB_par = "currentPlan" ;
    refine_facts.create_fact("balanced_planDB_par",balanced_planDB_par) ;
	      
    if(!Loci::makeQuery(refine_rules,refine_facts,
			"inner_nodes,fine_faces,volTag_blackbox")) {
      std::cerr << "adapt query failed!" << std::endl;
      Loci::Abort();
    }
       
    store<Loci::FineNodes> inner_nodes;
    inner_nodes = refine_facts.get_variable("inner_nodes");

	  
    gridDataP = new(refinedGridData) ;
    gridDataP->boundary_ids.clear();
    Loci::readBCfromVOG(meshFile, gridDataP->boundary_ids); 
       
    int num_nodes;
    createVOGNode(gridDataP->new_pos,
		  inner_nodes,
		  num_nodes,
		  refine_facts,
		  gridDataP->local_nodes
		  );
    if(Loci::MPI_rank ==0)cerr<< "num_nodes: " << num_nodes << " before chem run" <<  endl;
        
    store<Loci::FineFaces> fine_faces;
    fine_faces = refine_facts.get_variable("fine_faces");
	      
    int num_faces = 0;
    int num_cells = 0;
	      
    createVOGFace( num_nodes,
		   fine_faces,
		   refine_facts,
		   num_faces,
		   num_cells,
		   gridDataP->new_cl,
		   gridDataP->new_cr,
		   gridDataP->new_face2node,
		   gridDataP->local_faces,
		   gridDataP->local_cells
		   );

    if(Loci::MPI_rank ==0)cerr<< "num_faces: " << num_faces << " before chem run" <<  endl;
    if(Loci::MPI_rank ==0)cerr<< "num_cells: " << num_cells << " before chem run" <<  endl;
              
    {
      blackbox<std::vector<std::pair<std::string, entitySet> > > volTag_blackbox;
      volTag_blackbox = refine_facts.get_variable("volTag_blackbox"); 
      gridDataP->volTags =  *volTag_blackbox;
    }
  }

  extern bool use_simple_partition ;
  extern bool use_orb_partition ;
  extern bool use_sfc_partition ;
  extern bool load_cell_weights ;
  extern string cell_weight_file ;
  extern storeRepP cell_weight_store ; 
  bool redistribute_cell_weight(storeRepP old_store, storeRepP new_store);

  void fill_clone_proc( map<int,int> &mapdata, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  void redistribute_container(const vector<entitySet> &ptn,
                              const vector<entitySet> &ptn_t,
                              entitySet new_alloc,
                              storeRepP inRep,storeRepP outRep) ;
  extern void ORBPartition(const vector<vector3d<float> > &pnts,
                           vector<int> &procid,
                           MPI_Comm comm) ;
  void assignOwner(vector<pair<int,pair<int,int> > > &scratchPad,
                   vector<entitySet> ptn,
                   vector<entitySet> &out_ptn) ;
  void SFC_Partition_Mesh(const vector<entitySet> &local_nodes,
                          const vector<entitySet> &local_faces,
                          const vector<entitySet> &local_cells,
                          const store<vector3d<double> > &pos,
                          const Map &cl, const Map &cr,
                          const multiMap &face2node,
			  const store<string> &boundary_tags,
			  const store<int> &cell_weights,
                          vector<entitySet> &cell_ptn,
                          vector<entitySet> &face_ptn,
                          vector<entitySet> &node_ptn) ;
  Loci::storeRepP getCellPartitionWeights(entitySet localcelldom) {
    Loci::storeRepP wptr = Loci::DataXFER_DB.getItem("cellPartitionWeights") ;
    if(wptr == 0)
      return wptr ;
    store<int> cell_weights_parent ;
    cell_weights_parent = wptr ;
    store<pair<int,int> > c2pset ;
    Loci::storeRepP ptr = Loci::DataXFER_DB.getItem("c2p") ;
    if(ptr == 0)
      return ptr ;
    c2pset = ptr ;
    entitySet domc2p = c2pset.domain() ;
    int r = Loci::MPI_rank ;
    int p = Loci::MPI_processes ;
     
    int minParent = std::numeric_limits<int>::max() ;
    int minTarget = std::numeric_limits<int>::max() ;
    FORALL(domc2p,ii) {
      minParent = min(minParent,c2pset[ii].second) ;
      minTarget = min(minTarget,c2pset[ii].first) ;
    } ENDFORALL ;

    int parentOffset = minParent ;
    int targetOffset = minTarget ;
    MPI_Allreduce(&minParent,&parentOffset,1,MPI_INT,MPI_MIN,
		  MPI_COMM_WORLD) ;
    MPI_Allreduce(&minTarget,&targetOffset,1,MPI_INT,MPI_MIN,
		  MPI_COMM_WORLD) ;
    int sz = domc2p.size() ;
    vector<pair<int,int> > p2c(sz) ;
    int cnt = 0 ;
    FORALL(domc2p,ii) {
      p2c[cnt].first = c2pset[ii].second-parentOffset ;
      p2c[cnt].second = c2pset[ii].first-targetOffset ;
      cnt++ ;;
    } ENDFORALL ;

    sort(p2c.begin(),p2c.end()) ;
      
    // Now distribute c2p to processors in parent file ordering
    entitySet cwdom = cell_weights_parent.domain() ;
    int cwsz = cwdom.size() ;
    vector<int> parentsizes(p) ;
    MPI_Allgather(&cwsz,1,MPI_INT,&parentsizes[0],1,MPI_INT,
		  MPI_COMM_WORLD) ;
    vector<int> parentoffsets(p+1,0) ;
    for(int i=0;i<p;++i)
      parentoffsets[i+1] = parentoffsets[i]+parentsizes[i] ;
    
    vector<pair<int,int> > splits(p-1) ;
    for(int i=0;i<p-1;++i)
      splits[i] = pair<int,int>(parentoffsets[i+1],-1) ;

    Loci::parSplitSort(p2c,splits,MPI_COMM_WORLD) ;

    int psz = p2c.size() ;

    vector<pair<int,int> > cell2weights(psz) ;
    entitySet::const_iterator ii = cwdom.begin() ;
    cnt = parentoffsets[r] ;
    int mincell = p2c[0].second ;
    int cnt2 = 0 ;
    for(int i=0;i<psz;) {
      while(i<psz && p2c[i].first == cnt) {
	mincell = min(mincell,p2c[i].second) ;
	cell2weights[cnt2] =pair<int,int>(p2c[i].second,
					  cell_weights_parent[*ii]) ;
	cnt2++ ;
	i++ ;
      }
      cnt++ ;
      ii++ ;
    }
    int mincell_global =mincell ;
    MPI_Allreduce(&mincell,&mincell_global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;

    int dmsz = localcelldom.size() ;
    vector<int> cellsizes(p) ;
    MPI_Allgather(&dmsz,1,MPI_INT,&cellsizes[0],1,MPI_INT,
		  MPI_COMM_WORLD) ;
    parentoffsets[0] = mincell_global ;
    for(int i=0;i<p;++i)
      parentoffsets[i+1] = parentoffsets[i]+cellsizes[i] ;
    
    for(int i=0;i<p-1;++i)
      splits[i] = pair<int,int>(parentoffsets[i+1],-1) ;


    sort(cell2weights.begin(),cell2weights.end()) ;
    Loci::parSplitSort(cell2weights,splits,MPI_COMM_WORLD) ;

    store<int> cellweightschild ;
    cellweightschild.allocate(localcelldom) ;
    cnt = 0;

    FORALL(localcelldom,ii) {
      cellweightschild[ii] = cell2weights[cnt].second ;
      cnt++ ;
    } ENDFORALL ;
    return cellweightschild.Rep() ;
  }
    


#ifdef LOCI_USE_METIS
  vector<entitySet> AdaptMetisPartitionOfCells(const vector<entitySet> &local_cells,
                                               const Map &cl, const Map &cr,
                                               const store<string> &boundary_tags,
					       storeRepP cellwts) {
    entitySet dom = cl.domain() & cr.domain() ;
    entitySet refset = boundary_tags.domain() ;
    entitySet bcfaces = cr.preimage(refset).first ;
    dom -= bcfaces ;
    int cnt = dom.size() ;
    
    vector<pair<int,int> > rawMap(cnt*2) ;
    int j = 0 ;
    entitySet::const_iterator ei ;
    for(ei=dom.begin();ei!=dom.end();++ei) {
      rawMap[j++] = pair<int,int>(cl[*ei],cr[*ei]) ;
      rawMap[j++] = pair<int,int>(cr[*ei],cl[*ei]) ;
    }

    sort(rawMap.begin(),rawMap.end()) ;

    multiMap cell2cell ;
    entitySet all_cells ;
    for(int i=0;i<MPI_processes;++i)
      all_cells += local_cells[i] ;

    distributed_inverseMap(cell2cell,rawMap, all_cells,all_cells,local_cells) ;

    vector<pair<int,int> >().swap(rawMap) ; // Free up memory from rawMap
    int count = 0 ;
    idx_t size_map = local_cells[Loci::MPI_rank].size() ;
    vector<idx_t> size_adj(size_map) ;
    count = 0 ;
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) {
      size_adj[count] = cell2cell[*ei].size() ;
      ++count ;
    }

    vector<idx_t> part(size_map) ;
    vector<idx_t> xadj(size_map+1) ;
    idx_t edgecut ;
    vector<idx_t> vdist(Loci::MPI_processes + 1) ;
    int cmin = local_cells[0].Min();
    for(int i = 0; i < Loci::MPI_processes; i++) {
      cmin = min(local_cells[i].Min(), cmin);
    }
    
    // check local_cells to be consistent with parmetis partitioning
    for(int i=0;i<Loci::MPI_processes;++i) {
      if(local_cells[i].size() == 0 ||
         (int)local_cells[i].size() != local_cells[i].Max()-local_cells[i].Min()+1) {
        cerr << "invalid local cell set, p=" << i
             << ", local_cells=" << local_cells[i] << endl ;
      }
      if(i>0 && local_cells[i-1].Max()+1 != local_cells[i].Min()) {
        cerr << "gap between processors in cell numbering" << endl ;
      }
    }
    
      
    
    edgecut = 0 ;
    xadj[0] = 0 ;
    for(idx_t i = 0; i < size_map; ++i)
      xadj[i+1] = xadj[i] + size_adj[i] ;

    idx_t min_size = size_adj[0] ;
    for(idx_t i = 0; i < size_map; ++i)
      min_size = min(min_size,size_adj[i]) ;
    if(min_size == 0) 
      cerr << "cell with no adjacency!" << endl ;
    
    idx_t tot = xadj[size_map] ;
    vector<idx_t> adjncy(tot) ;
    count = 0 ;
    for(entitySet::const_iterator ei = local_cells[Loci::MPI_rank].begin(); ei != local_cells[Loci::MPI_rank].end(); ++ei) {
      size_t sz = cell2cell[*ei].size() ;
      for(size_t i = 0; i != sz; ++i)        {
	adjncy[count] = cell2cell[*ei][i] - cmin ;
	count ++ ;
      }
    }
    cell2cell.setRep(multiMap().Rep()) ;// Free up memory from multiMap

    vdist[0] = 0 ;
    for(int i = 1; i <= Loci::MPI_processes; ++i)
      vdist[i] = vdist[i-1] + local_cells[i-1].size() ;
    idx_t top = vdist[Loci::MPI_processes] ;
    
    bool trouble = false ;
    for(idx_t i=0;i<tot;++i)
      if(adjncy[i] >= top)
        trouble = true ;
    if(trouble)
      cerr << "adjacency list contains out of bounds reference" << endl ;
    
    MPI_Comm mc = MPI_COMM_WORLD ;
    idx_t nparts = Loci::MPI_processes ; // number of partitions
    idx_t wgtflag = 0 ;
    idx_t numflag = 0 ;
    idx_t options = 0 ;

    //    Loci::storeRepP ptr = getCellPartitionWeights(local_cells[MPI_rank]) ;
    Loci::storeRepP ptr = cellwts ;
    
    if(ptr != 0) {
      store<int> cell_weights ;
      cell_weights = ptr ;
        
      // compute necessary ParMETIS data-structure
      wgtflag = 2 ;           // weights on the vertices only
      idx_t ncon = 2 ;          // number of weights per vertex
      idx_t tpwgts_len = ncon*nparts ;
      vector<metisreal_t> tpwgts(tpwgts_len) ;
      
      for(idx_t i=0;i<tpwgts_len;++i)
	tpwgts[i] = 1.0 / double(nparts) ;
      
      vector<metisreal_t> ubvec(ncon) ;
      for(idx_t i=0;i<ncon;++i)
	ubvec[i] = 1.05 ;     // as recommended by the ParMETIS manual
      
      // now construct the vertex weights
      vector<idx_t> vwgt(ncon*size_map) ;
      int cnt = 0 ;
      for(entitySet::const_iterator
	    ei=local_cells[Loci::MPI_rank].begin();
	  ei!=local_cells[Loci::MPI_rank].end();++ei,cnt+=ncon) {
	// first weight for cell is 1 (the cell computation)
	vwgt[cnt] = 1 ;
	// the second weight is from the store cell_weights[*ei]
	vwgt[cnt+1] = cell_weights[*ei] ;
      }

      // now call the ParMETIS routine (V3)
      ParMETIS_V3_PartKway(&vdist[0],&xadj[0],&adjncy[0],&vwgt[0],NULL,
			   &wgtflag,&numflag,&ncon,&nparts,
			   &tpwgts[0],&ubvec[0],&options,&edgecut,&part[0],
			   &mc) ;
      
    } else {
      //      idx_t ncon = 1 ;
      idx_t ncon = 2 ;
      idx_t tpwgts_len = ncon*nparts ;
      vector<metisreal_t> tpwgts(tpwgts_len) ;
      for(idx_t i=0;i<tpwgts_len;++i)
        tpwgts[i] = 1.0 / double(nparts) ;

      vector<metisreal_t> ubvec(ncon) ;
      for(idx_t i=0;i<ncon;++i)
	ubvec[i] = 1.05 ;     // as recommended by the ParMETIS manual
      wgtflag = 2 ;
      // now construct the vertex weights
      vector<idx_t> vwgt(ncon*size_map) ;

      int cnt = 0 ;
      for(int i=0;i<size_map;++i) {
	// first weight for cell is 1 (the cell computation)
	vwgt[cnt] = 1 ;
	// the second weight is from the store cell_weights[*ei]
	vwgt[cnt+1] =  size_adj[i] ;
	cnt += ncon ;
      }
      ParMETIS_V3_PartKway(&vdist[0],&xadj[0],&adjncy[0],&vwgt[0],NULL,
                           &wgtflag,&numflag,&ncon,&nparts,
                           &tpwgts[0],&ubvec[0],&options,&edgecut,
                           &part[0],&mc) ;
    }

    if(Loci::MPI_rank == 0)
      Loci::debugout << " Parmetis Edge cut   " <<  edgecut << endl ;

    //find the partition ptn given by Metis
    vector<entitySet> ptn ;

    for(int i = 0; i < Loci::MPI_processes; ++i)
      ptn.push_back(EMPTY) ;
    cmin = local_cells[Loci::MPI_rank].Min() ;
    for(int i=0;i<size_map;++i) {
      ptn[part[i]] += i + cmin ;
    }
   
    return ptn;

  }
#endif /* LOCI_USE_METIS */

  
  void ORB_Partition_Mesh(const vector<entitySet> &local_nodes,
                          const vector<entitySet> &local_faces,
                          const vector<entitySet> &local_cells,
                          const store<vector3d<double> > &pos,
                          const Map &cl, const Map &cr,
                          const multiMap &face2node,
			  const store<string> &boundary_tags,
                          vector<entitySet> &cell_ptn,
                          vector<entitySet> &face_ptn,
                          vector<entitySet> &node_ptn) ;
  vector<entitySet> partitionFaces(vector<entitySet> cell_ptn, const Map &cl,
                                   const Map &cr,
				   const store<string> &boundary_tags) ;
  
  //Input: Mapping from faces to its right cells.
  //Output: Entities of boundary cells.
  entitySet getBoundaryCells(const MapRepP tmp_cr_rep) {
    entitySet cri = tmp_cr_rep->image(tmp_cr_rep->domain()) ;
    return(cri & interval(Loci::UNIVERSE_MIN,-1)) ;
  }

  void copyGridStructures( entitySet nodes, entitySet faces, entitySet cells,
			   const store<vector3d<double> > &t_pos,
			   const Map &tmp_cl, const Map &tmp_cr,
			   const multiMap &tmp_face2node,
			   store<vector3d<double> > &pos, Map &cl, Map &cr,
			   multiMap &face2node) {

    entitySet boundary_cells = getBoundaryCells(Loci::MapRepP(tmp_cr.Rep()));

    dMap identity_map;
    FORALL(nodes, ei) {
      identity_map[ei] = ei;
    } ENDFORALL ;
    FORALL(faces, ei) {
      identity_map[ei] = ei ;
    } ENDFORALL ;
    FORALL(cells, ei) {
      identity_map[ei] = ei ;
    } ENDFORALL ;
    FORALL(boundary_cells, ei) {
      identity_map[ei] = ei ;
    } ENDFORALL ;

    pos = t_pos.Rep()->remap(identity_map);
    cl = tmp_cl.Rep()->remap(identity_map);
    cr = tmp_cr.Rep()->remap(identity_map);
    face2node = MapRepP(tmp_face2node.Rep())->get_map();

  }


  vector<entitySet> transposePtn(const vector<entitySet> &ptn);
  vector<entitySet> partitionNodes(vector<entitySet> face_ptn, MapRepP face2node,entitySet old_node_dom);
  inline bool fieldSort(const std::pair<Entity,Entity> &p1,
                        const std::pair<Entity,Entity> &p2) {
    return p1.first < p2.first ;
  }

  
  void create_face_info(fact_db &facts);
  void create_ref(fact_db &facts) ;
  void create_ghost_cells(fact_db &facts);
  
  std::vector<int> simplePartitionVec(int mn, int mx, int p);
  vector<entitySet> simplePartition(int mn, int mx, MPI_Comm comm);
  vector<sequence> transposeSeq(const vector<sequence> sv);
}
void colorMatrix(Map &cl, Map &cr, multiMap &face2node) ;

namespace Loci{
  
  // Convert container from global numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator(comm)
  storeRepP Global2FileOrder(storeRepP sp, entitySet dom, int &offset,
                             fact_db::distribute_infoP dist, MPI_Comm comm) {
   
    // Now get global to file numbering
    dMap g2f ;
    g2f = dist->g2fv[0].Rep() ; // FIX THIS

    // Compute map from local numbering to file numbering
    Map newnum ;
    newnum.allocate(dom) ;
    FORALL(dom,i) {
      newnum[i] = g2f[i] ;
    } ENDFORALL ;

    int imx = std::numeric_limits<int>::min() ;
    int imn = std::numeric_limits<int>::max() ;

    // Find bounds in file numbering from this processor
    FORALL(dom,i) {
      imx = max(newnum[i],imx) ;
      imn = min(newnum[i],imn) ;
    } ENDFORALL ;
    
    // Find overall bounds
    imx = GLOBAL_MAX(imx) ;
    imn = GLOBAL_MIN(imn) ;

    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    // Get partitioning of file numbers across processors
    vector<entitySet> out_ptn = simplePartition(imn,imx,comm) ;

    // Now compute where to send data to put in file ordering
    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;

    // Loop over processors and compute sets of entities to send
    // To efficiently compute this mapping, first sort the transpose
    // of the newnum map to quickly find the set of entities to send
    // without searching entire newnum map for each processor
    vector<pair<int,int> > file2num(dom.size()) ;
    size_t cnt = 0 ;
    FORALL(dom,ii) {
      file2num[cnt].first = newnum[ii] ;
      file2num[cnt].second = ii ;
      cnt++ ;
    } ENDFORALL ;
    sort(file2num.begin(),file2num.end()) ;

    // Check each processor, find out which sets to send
    cnt = 0 ;
    for(int i=0;i<p;++i) {
      int mxi = out_ptn[i].Max() ;
      while(cnt < file2num.size() && file2num[cnt].first <= mxi) {
        send_sets[i] += file2num[cnt].second ;
        cnt++ ;
      }
      sequence s ;
      FORALL(send_sets[i],j) {
        s+= newnum[j] ;
      } ENDFORALL ;
      send_seqs[i] = s ;
    }

    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;


    // shift by the offset
    offset = out_ptn[prank].Min() ;
    for(int i=0;i<p;++i)
      recv_seqs[i] <<= offset ;

    // Compute allocation domain
    entitySet file_dom ;
    for(int i=0;i<p;++i)
      file_dom += entitySet(recv_seqs[i]) ;

    // allocate store over shifted domain
    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(file_dom) ;

    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
		  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
		  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    return qcol_rep ;
  } 


  
  std::vector<entitySet> getDist( Loci::entitySet &faces,
                                  Loci::entitySet &cells,
                                  Map &cl, Map &cr, multiMap &face2node) {
  
    // First establish current distribution of entities across processors
    std::vector<Loci::entitySet> ptn(Loci::MPI_processes) ; // entity Partition

    // Get entity distributions
  
    faces = face2node.domain() ;
    entitySet allFaces = Loci::all_collect_entitySet(faces) ;
    std::vector<int> facesizes(MPI_processes) ;
    int  size = faces.size() ;
    MPI_Allgather(&size,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    int  cnt = allFaces.Min() ;
    for(int i=0;i<MPI_processes;++i) {
      ptn[i] += interval(cnt,cnt+facesizes[i]-1) ;
      cnt += facesizes[i] ;
    }
    
    entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
    entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
    entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
    int mn = geom_cells.Min() ;
    int mx = geom_cells.Max() ;
    std:: vector<int> pl = Loci::simplePartitionVec(mn,mx,MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      ptn[i] += interval(pl[i],pl[i+1]-1) ;
    faces = allFaces ;
    cells = geom_cells ;
    return ptn ;
  }

  extern  bool useDomainKeySpaces  ;
  extern void remapGrid(vector<entitySet> &node_ptn,
                        vector<entitySet> &face_ptn,
                        vector<entitySet> &cell_ptn,
                        vector<entitySet> &node_ptn_t,
                        vector<entitySet> &face_ptn_t,
                        vector<entitySet> &cell_ptn_t,
                        store<vector3d<double> > &t_pos, Map &tmp_cl,
                        Map &tmp_cr, multiMap &tmp_face2node,
                        vector<entitySet> &bcsurf_ptn,
                        store<string> &tmp_boundary_names,
                        store<string> &tmp_boundary_tags,
                        entitySet nodes, entitySet faces, entitySet cells,
                        store<vector3d<double> > &pos, Map &cl, Map &cr,
                        multiMap &face2node,
                        store<string> &boundary_names, 
                        store<string> &boundary_tags, 
                        entitySet bcsurfset,
                        fact_db &facts) ;

  bool inputFVMGrid(fact_db &facts,
                    vector<entitySet>& local_nodes,
                    vector<entitySet>& local_faces,
                    vector<entitySet>& local_cells,
                    store<vector3d<double> >& t_pos,
                    Map& tmp_cl,
                    Map& tmp_cr,
                    multiMap& tmp_face2node,
                    vector<pair<int,string> >& boundary_ids,
                    vector<pair<string,entitySet> >& volTags,
                    storeRepP cellwts) {
    double t1 = MPI_Wtime() ;
    // Identify boundary tags
    if(Loci::MPI_processes == 1) {
      entitySet local_boundary_cells = getBoundaryCells(MapRepP(tmp_cr.Rep()));
      entitySet global_boundary_cells = all_collect_entitySet(local_boundary_cells) ;

      int npnts = local_nodes[0].size();
      int nfaces = local_faces[0].size();
      int ncells = local_cells[0].size();

      entitySet nodes = facts.get_distributed_alloc(npnts,0).first ; // FIX THIS
      entitySet faces = facts.get_distributed_alloc(nfaces,0).first ;
      entitySet cells = facts.get_distributed_alloc(ncells,0).first;

      store<vector3d<double> > pos ;
      Map cl ;
      Map cr ;
      multiMap face2node ;
      Loci::copyGridStructures(nodes, faces, cells,
                               t_pos, tmp_cl, tmp_cr, tmp_face2node,
                               pos, cl, cr, face2node);
      store<string> boundary_names ;
      store<string> boundary_tags ;
      boundary_names.allocate(global_boundary_cells) ;
      boundary_tags.allocate(global_boundary_cells) ;
      cout << " boundaries identified as:" ;
      
      FORALL(global_boundary_cells, bc) {
        char buf[512] ;
	bzero(buf,512) ;
        snprintf(buf,511,"BC_%d",-bc) ;
        boundary_tags[bc] = string(buf) ;
        boundary_names[bc] = string(buf) ;
	cout << " " << boundary_names[bc] ;
      } ENDFORALL ;

      for(size_t i=0;i<boundary_ids.size();++i) {
        int id = boundary_ids[i].first ;
        if(global_boundary_cells.inSet(-id))
          boundary_names[-id] = boundary_ids[i].second ;
      }

      FORALL(global_boundary_cells, bc) {
	cout << " " << boundary_names[bc] ;
      } ENDFORALL ;
      cout << endl ;
      
      facts.create_fact("cl", cl) ;
      facts.create_fact("cr", cr) ;
      facts.create_fact("pos", pos) ;
      facts.create_fact("face2node",face2node) ;
      facts.create_fact("boundary_names", boundary_names) ;
      facts.create_fact("boundary_tags", boundary_tags) ;

      int cells_base = local_cells[0].Min() ;
      for(size_t i=0;i<volTags.size();++i) {
        param<string> Tag ;
        *Tag = volTags[i].first ;
        Tag.set_entitySet(volTags[i].second >> cells_base) ;
        std::ostringstream oss ;
        oss << "volumeTag(" << volTags[i].first << ")" ;
        facts.create_fact(oss.str(),Tag) ;
      }
        
      return true ;
    }

    entitySet local_boundary_cells = getBoundaryCells(MapRepP(tmp_cr.Rep()));
    entitySet global_boundary_cells = all_collect_entitySet(local_boundary_cells) ;
    int maxc = 0 ;
    entitySet dom = tmp_cr.domain()&tmp_cl.domain() ;
    FORALL(dom,ii) {
      maxc = max(maxc,max(tmp_cr[ii],tmp_cl[ii])) ;
    } ENDFORALL ;
    int maxctot = maxc ;
    MPI_Allreduce(&maxc,&maxctot,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    entitySet boundary_taglist = global_boundary_cells ;
    int ref_base = maxctot+1 ;
    entitySet refSet = interval(ref_base,ref_base+(boundary_taglist.size()-1)) ;
    
    store<int> boundary_info ;
    boundary_info.allocate(refSet) ;
    map<int,int> boundary_remap ;
    int ct = 0 ;
    FORALL(boundary_taglist,ii) {
      boundary_remap[ii]=ct+ref_base ;
      boundary_info[ct+ref_base]=ii ;
      ct++ ;
    } ENDFORALL ;

    entitySet domcr = tmp_cr.domain() ;
    FORALL(domcr,fc) {
      if(tmp_cr[fc] < 0)
	tmp_cr[fc] = boundary_remap[tmp_cr[fc]] ;
    } ENDFORALL ;

    store<string> tmp_boundary_tags,tmp_boundary_names ;
    tmp_boundary_names.allocate(refSet) ;
    tmp_boundary_tags.allocate(refSet) ;

    map<int,int> bcmap ;
    FORALL(refSet, ii) {
      int bc = boundary_info[ii] ;
      char buf[128] ;
      bzero(buf,128) ;
      snprintf(buf,127,"BC_%d",-bc) ;
      bcmap[-bc] = ii ;
      tmp_boundary_tags[ii] = string(buf) ;
      tmp_boundary_names[ii] = string(buf) ;
    } ENDFORALL ;

    for(size_t i=0;i<boundary_ids.size();++i) {
      int id = boundary_ids[i].first ;
      map<int,int>::const_iterator mi = bcmap.find(id) ;
      if(mi != bcmap.end())
	tmp_boundary_names[mi->second] = boundary_ids[i].second ;
      else 
	debugout << "id " << id << " for boundary surface " << boundary_ids[i].second
                 << " not found!" << endl ;
    }

      
    if(Loci::MPI_rank == 0) {
      Loci::debugout << " boundaries identified as:" ;
      FORALL(refSet, bc) {
	debugout << " " << tmp_boundary_names[bc] ;
      } ENDFORALL ;
      Loci::debugout << endl ;
    }
    
    REPORTMEM() ;

    vector<entitySet> cell_ptn,face_ptn,node_ptn ;

    store<int> cell_weights ;
    if(cellwts != 0)
      cell_weights = cellwts ;
    else {
      Loci::storeRepP ptr = getCellPartitionWeights(local_cells[MPI_rank]) ;
      if(ptr != 0) {
	cell_weights = ptr ;
      }
    }

    enum {ORB=0,SFC=1,SIMPLE=2,GRAPH=3} partitioner_type = GRAPH ;
    if(use_orb_partition)
      partitioner_type = ORB ;
    if(use_sfc_partition)
      partitioner_type = SFC ;
    if(use_simple_partition) 
      partitioner_type = SIMPLE ;

    if(partitioner_type == GRAPH) {
      int lcpp = local_cells[MPI_rank].size() ;
      int mincpp = lcpp ;
      MPI_Allreduce(&lcpp,&mincpp,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      if(mincpp < metis_cpp_threshold) {
	partitioner_type = SFC ; // Space filling curve partitioner
	debugout << "switching from metis to space filling curve partitioner"
		 << " mincpp = " << mincpp << endl ;
      }
#ifndef LOCI_USE_METIS
      partitioner_type = SFC ; // Space filling curve partitioner if METIS
      // not available
#endif
    }
    switch(partitioner_type) {
    case ORB:
      ORB_Partition_Mesh(local_nodes, local_faces, local_cells,
			 t_pos, tmp_cl, tmp_cr, tmp_face2node,
			 tmp_boundary_tags,
			 cell_ptn,face_ptn,node_ptn) ;
      break ;
    case SFC:
      SFC_Partition_Mesh(local_nodes, local_faces, local_cells,
			 t_pos, tmp_cl, tmp_cr, tmp_face2node,
			 tmp_boundary_tags,
			 cell_weights,
			 cell_ptn,face_ptn,node_ptn) ;
      break ;
    case SIMPLE: // Simple partition
      {
	cell_ptn = vector<entitySet>(MPI_processes) ;
        cell_ptn[MPI_rank] = local_cells[MPI_rank] ;
	REPORTMEM() ;
	face_ptn = partitionFaces(cell_ptn,tmp_cl,tmp_cr,tmp_boundary_tags) ;
	REPORTMEM() ;

	node_ptn = partitionNodes(face_ptn,
				  MapRepP(tmp_face2node.Rep()),
                                t_pos.domain()) ;
	
      }
      break ;
#ifdef LOCI_USE_METIS
    case GRAPH: // METIS partition
      {
        cell_ptn = AdaptMetisPartitionOfCells(local_cells,tmp_cl,tmp_cr,tmp_boundary_tags,cellwts) ;
	//        cell_ptn = newMetisPartitionOfCells(local_cells,tmp_cl,tmp_cr,tmp_boundary_tags,cellwts) ;
	REPORTMEM() ;
	face_ptn = partitionFaces(cell_ptn,tmp_cl,tmp_cr,tmp_boundary_tags) ;
	REPORTMEM() ;

	node_ptn = partitionNodes(face_ptn,
				  MapRepP(tmp_face2node.Rep()),
				  t_pos.domain()) ;
      }
      break ;
#endif
    default: // SFC partition no weight balance
      {
	store<int> cell_weights ;

	SFC_Partition_Mesh(local_nodes, local_faces, local_cells,
			   t_pos, tmp_cl, tmp_cr, tmp_face2node,
			   tmp_boundary_tags,
			   cell_weights,
			   cell_ptn,face_ptn,node_ptn) ;
      }
      break ;
    }

    vector<entitySet> bcsurf_ptn(MPI_processes) ;
    entitySet refset = tmp_boundary_tags.domain() ;
      
    // round robin allocate boundary faces
    // NOTE!!!! This could be improved.
    int cnt = 0 ;
    //    int refsetsz = refset.size() ;
    FORALL(refset,ii) {
      if(cnt == MPI_rank)
	bcsurf_ptn[cnt] += ii ;
      cnt++ ;
      if(cnt == MPI_processes)
	cnt = 0 ;
    } ENDFORALL ;
    REPORTMEM() ;
      
    vector<entitySet> cell_ptn_t = transposePtn(cell_ptn) ;
    vector<entitySet> face_ptn_t = transposePtn(face_ptn) ;
    vector<entitySet> node_ptn_t = transposePtn(node_ptn) ;

    int newnodes = 0 ;
    for(int p=0;p<MPI_processes;++p)
      newnodes += node_ptn_t[p].size() ;

    vector<int> node_alloc(newnodes) ;
    int i=0;
    for(int p=0;p<MPI_processes;++p)
      FORALL(node_ptn_t[p], ni) {
        node_alloc[i++] = ni ;
      } ENDFORALL;

    entitySet nodes = facts.get_distributed_alloc(node_alloc,0).first ;// FIX THIS
    node_alloc.resize(0) ;

    int newfaces = 0 ;
    for(int p=0;p<MPI_processes;++p)
      newfaces += face_ptn_t[p].size() ;

    vector<int> face_alloc(newfaces) ;
    i = 0 ;
    for(int p=0;p<MPI_processes;++p)
      FORALL(face_ptn_t[p], ni) {
        face_alloc[i++] = ni ;
      }ENDFORALL;

    int fk = facts.getKeyDomain("Faces") ;

    if(!useDomainKeySpaces)
      fk = 0 ;

    entitySet faces = facts.get_distributed_alloc(face_alloc,fk).first ;
    face_alloc.resize(0) ;

    int newcells = 0 ;
    for(int p=0;p<MPI_processes;++p)
      newcells += cell_ptn_t[p].size() ;

    vector<int> cell_alloc(newcells) ;
    i = 0 ;
    for(int p=0;p<MPI_processes;++p)
      FORALL(cell_ptn_t[p], ni) {
        cell_alloc[i++] = ni ;
      }ENDFORALL;

    entitySet cells = facts.get_distributed_alloc(cell_alloc,0).first ;//Fix This

    Loci::debugout << "nodes = " << nodes << ", size= "
                   << nodes.size() << endl;
    Loci::debugout << "faces = " << faces << ", size = "
                   << faces.size() << endl ;
    Loci::debugout << "cells = " << cells << ", size = "
                   << cells.size() << endl ;

    vector<int> bcsurf_alloc(bcsurf_ptn[MPI_rank].size()) ;
    i=0 ;
    FORALL(bcsurf_ptn[MPI_rank],ii) {
      bcsurf_alloc[i++] = ii ;
    } ENDFORALL ;
    
    entitySet bcsurfset = facts.get_distributed_alloc(bcsurf_alloc,0).first ;// FIX THIS
    
    REPORTMEM() ;
    Map cl, cr ;
    multiMap face2node ;
    tmp_cl.Rep()->setDomainKeySpace(fk) ;
    tmp_cr.Rep()->setDomainKeySpace(fk) ;
    tmp_face2node.Rep()->setDomainKeySpace(fk) ;

    store<vector3d<double> > pos ;
    store<string> boundary_names,boundary_tags ;
    remapGrid(node_ptn, face_ptn, cell_ptn,
              node_ptn_t, face_ptn_t, cell_ptn_t,
              t_pos, tmp_cl, tmp_cr, tmp_face2node,
	      bcsurf_ptn,tmp_boundary_names,tmp_boundary_tags,
              nodes, faces, cells,
              pos, cl, cr, face2node,
	      boundary_names,boundary_tags, bcsurfset, facts);
    REPORTMEM() ;

    facts.create_fact("cl", cl) ;
    facts.create_fact("cr", cr) ;
    facts.create_fact("pos", pos) ;

    facts.create_fact("face2node",face2node) ;
    facts.create_fact("boundary_names", boundary_names) ;
    facts.create_fact("boundary_tags", boundary_tags) ;

    // update remap from global to file numbering for faces after sorting
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    dMap g2f ;
    g2f = df->g2fv[0].Rep() ; // FIX THIS

    int cells_base=local_cells[0].Min() ;
    for(size_t i=0;i<volTags.size();++i) {
      param<string> Tag ;
      *Tag = volTags[i].first ;

      // Map entitySet to new ordering
      entitySet inputTag = (volTags[i].second >> cells_base) ;
      // Map entitySet to new ordering
      entitySet tagset ;
      FORALL(cells,cc) {
        if(inputTag.inSet(g2f[cc])) 
          tagset += cc ;
      } ENDFORALL ;
      Tag.set_entitySet(tagset) ;
      ostringstream oss ;
      oss << "volumeTag(" << volTags[i].first << ")" ;
      facts.create_fact(oss.str(),Tag) ;
    }
    //	if(collect_perf_data)
    //		perfAnalysis->stop_timer(read_file_timer);
    double t2 = MPI_Wtime() ;
    debugout << "Time to process and partition adapted mesh is " << t2-t1
             << endl ;
    REPORTMEM() ;
    return true ;
    
  }

  
  bool setupFVMGridFromContainer(fact_db &facts,
                                 vector<entitySet>& local_nodes,
                                 vector<entitySet>& local_faces,
                                 vector<entitySet>& local_cells,
                                 store<vector3d<double> >& t_pos,
                                 Map& tmp_cl,
                                 Map& tmp_cr,
                                 multiMap& tmp_face2node,
                                 vector<pair<int,string> >& boundary_ids,
                                 vector<pair<string,entitySet> >& volTags,
				 storeRepP cellwts) {
       
    if(!inputFVMGrid(facts,
                     local_nodes,
                     local_faces,
                     local_cells,
                     t_pos,
                     tmp_cl,
                     tmp_cr,
                     tmp_face2node,
                     boundary_ids,
                     volTags,cellwts))
      return false ;
    REPORTMEM() ;

    create_face_info(facts) ;
    create_ref(facts) ;
    create_ghost_cells(facts) ;

    return true ;
  }
}


namespace Loci {
  //create a new store new_pos from pos and inner_nodes
  //re_number the nodes in pos and inner_nodes,
  //and then redistribute them across the processes
  void createVOGNode(store<vector3d<double> > &new_pos,
                     const store<Loci::FineNodes> &inner_nodes,
                     int& num_nodes,
                     fact_db & facts,//in global numbering
                     vector<entitySet>& nodes_ptn
                     ){
    //get store pos
    store<vector3d<double> > pos;
    pos =  facts.get_variable("pos");
       
    if(MPI_processes == 1){
      //firsr write out numNodes
      long  num_original_nodes  = pos.domain().size();
      long  num_inner_nodes  = 0;
      FORALL(inner_nodes.domain(), cc){
        num_inner_nodes += inner_nodes[cc].size();
      }ENDFORALL;
    
      int node_base = 0;
      long npnts = num_original_nodes + num_inner_nodes;
      entitySet new_domain= interval(node_base, node_base+npnts-1);
      new_pos.allocate(new_domain);
    
      entitySet::const_iterator nei = new_domain.begin();
      entitySet::const_iterator ei = pos.domain().begin();
    
      for(long count = 0; count < num_original_nodes; count++, nei++, ei++){
        new_pos[*nei] = pos[*ei];
      }

      Loci::constraint faces, geom_cells;
      faces = facts.get_variable("faces");
      geom_cells = facts.get_variable("geom_cells");
      entitySet local_edges = facts.get_variable("edge2node")->domain();

      FORALL(local_edges, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++, nei++){
          new_pos[*nei] = inner_nodes[cc][i];
        }
      }ENDFORALL;
      FORALL(*geom_cells, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++, nei++){
          new_pos[*nei] = inner_nodes[cc][i];
        }
      }ENDFORALL; 
      FORALL(*faces, cc){
        for(unsigned int i = 0; i < inner_nodes[cc].size(); i++, nei++){
          new_pos[*nei] = inner_nodes[cc][i];
        }
      }ENDFORALL;
      num_nodes = new_domain.size();
      nodes_ptn.resize(1);
      nodes_ptn[0] = new_domain;
      return;
    }

 
 
  
    store<vector3d<double> > pos_t; //temp container 
    vector<entitySet> temp_node_ptn;
    entitySet my_temp_nodes;
    {
      //reorder store first, from global to io entities
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      constraint  my_faces, my_geom_cells, my_edges; 
   
      my_faces = facts.get_variable("faces");
      my_geom_cells = facts.get_variable("geom_cells");
      my_edges =  facts.get_variable("edges");
    
      entitySet local_edges = *my_edges;
      entitySet local_faces =  *my_faces;
      entitySet local_cells =   *my_geom_cells;
      entitySet local_nodes =  pos.domain();

      //before write out,create stores for pos and inner_nodes which
      //are ordered across processors in the file numbering, the domain of this container
      //shifted by offset is the actual file numbering. offset will be modified after function call  
     
      int noffset, eoffset, coffset, foffset;
      noffset = 0;
      store<vector3d<double> > pos_io;
      pos_io = Loci::Global2FileOrder(pos.Rep(), local_nodes, noffset, dist, MPI_COMM_WORLD) ;
      entitySet file_nodes = pos_io.domain(); 
    
    
      eoffset = 0;
      store<Loci::FineNodes> edge_inner_nodes;
      edge_inner_nodes = Loci::Global2FileOrder(inner_nodes.Rep(),local_edges,eoffset,dist,MPI_COMM_WORLD) ;
      entitySet file_edges = edge_inner_nodes.domain();
      
      coffset= 0;
      // Create container 
      store<Loci::FineNodes> cell_inner_nodes;
      cell_inner_nodes = Loci::Global2FileOrder(inner_nodes.Rep(),local_cells,coffset,dist,MPI_COMM_WORLD) ;
      entitySet file_cells = cell_inner_nodes.domain();
  
      foffset= 0;
      // Create container
      store<Loci::FineNodes> face_inner_nodes;
      face_inner_nodes = Loci::Global2FileOrder(inner_nodes.Rep(),local_faces,foffset,dist,MPI_COMM_WORLD) ;
      entitySet file_faces = face_inner_nodes.domain();
      
      //Now allocate temp entitySet in file numbering
   
   
      //compute the size of pos
      int num_pos_nodes = file_nodes.size();
      std::vector<int> pos_sizes(Loci::MPI_processes) ;
      pos_sizes = Loci::all_collect_sizes(num_pos_nodes);
      
      //compute the size of inner_nodes
      int num_local_edge_nodes = 0;
      FORALL(file_edges, cc){
        num_local_edge_nodes += edge_inner_nodes[cc].size();
      }ENDFORALL;
  
      int num_local_cell_nodes = 0;
      FORALL(file_cells, cc){
        num_local_cell_nodes += cell_inner_nodes[cc].size();
      }ENDFORALL;
  
      int num_local_face_nodes = 0;
      FORALL(file_faces, cc){
        num_local_face_nodes += face_inner_nodes[cc].size();
      }ENDFORALL;
      std::vector<int> inner_edge_nodes_sizes(Loci::MPI_processes) ;
      inner_edge_nodes_sizes = Loci::all_collect_sizes(num_local_edge_nodes);
    
      std::vector<int> inner_cell_nodes_sizes(Loci::MPI_processes) ;
      inner_cell_nodes_sizes = Loci::all_collect_sizes(num_local_cell_nodes);

      std::vector<int> inner_face_nodes_sizes(Loci::MPI_processes) ;
      inner_face_nodes_sizes = Loci::all_collect_sizes(num_local_face_nodes);
    
      //set up the offset
      int my_id = Loci::MPI_rank;
      int num_procs = Loci::MPI_processes;
      noffset = 0;
      for(int i = 0; i < my_id; i++){
        noffset += pos_sizes[i];
      }
    
      eoffset = 0;
      for(int i = 0; i < num_procs; i++){
        eoffset += pos_sizes[i];
      }
      for(int i = 0; i < my_id; i++){
        eoffset += inner_edge_nodes_sizes[i];
      }
    
      coffset = 0;
      for(int i = 0; i < num_procs; i++){
        coffset += pos_sizes[i]+inner_edge_nodes_sizes[i];
      }
      for(int i = 0; i < my_id; i++){
        coffset += inner_cell_nodes_sizes[i];
      }
    
    
      foffset = 0;
      for(int i = 0; i < num_procs; i++){
        foffset += pos_sizes[i]+inner_edge_nodes_sizes[i]+inner_cell_nodes_sizes[i];
      }
      for(int i = 0; i < my_id; i++){
        foffset += inner_face_nodes_sizes[i];
      }
    

      num_nodes = 0;
      for(int i = 0; i < num_procs; i++){
        num_nodes += pos_sizes[i]+inner_edge_nodes_sizes[i]+
          inner_cell_nodes_sizes[i]+inner_face_nodes_sizes[i];
      }
      
      //create a new store pos_t in file numbering that contains all the nodes  
      entitySet my_original_nodes, my_inner_edge_nodes;
      entitySet my_inner_cell_nodes, my_inner_face_nodes;
      if(pos_sizes[my_id]>0) my_original_nodes = interval(noffset, noffset+pos_sizes[my_id]-1);
      if(inner_edge_nodes_sizes[my_id]>0)my_inner_edge_nodes = interval(eoffset, eoffset+inner_edge_nodes_sizes[my_id]-1);
      if(inner_cell_nodes_sizes[my_id]>0)my_inner_cell_nodes = interval(coffset, coffset+inner_cell_nodes_sizes[my_id]-1);
      if(inner_face_nodes_sizes[my_id]>0)my_inner_face_nodes = interval(foffset, foffset+inner_face_nodes_sizes[my_id]-1);
      my_temp_nodes = my_original_nodes + my_inner_edge_nodes + 
        my_inner_cell_nodes + my_inner_face_nodes;
      pos_t.allocate(my_temp_nodes);
    
      //fill container pos_t will pos_io
      entitySet::const_iterator ti = my_temp_nodes.begin();
      FORALL(file_nodes, ei){
        pos_t[*ti] = pos_io[ei];
        ti++;
      }ENDFORALL;
      //fill container pos_t will edge nodes
      FORALL(file_edges, ei){
        for(unsigned int i = 0; i < edge_inner_nodes[ei].size(); i++){
          pos_t[*ti] = edge_inner_nodes[ei][i];
          ti++;
        }
     
      }ENDFORALL;
      //fill container pos_t will cell nodes
      FORALL(file_cells, ei){
        for(unsigned int i = 0; i < cell_inner_nodes[ei].size(); i++){
          pos_t[*ti] = cell_inner_nodes[ei][i];
          ti++;
        }
      
      }ENDFORALL;
      //fill container pos_t will face nodes
      FORALL(file_faces, ei){
        for(unsigned int i = 0; i < face_inner_nodes[ei].size(); i++){
          pos_t[*ti] = face_inner_nodes[ei][i];
          ti++;
        }
      }ENDFORALL;
      if(ti != my_temp_nodes.end()){
        debugout<<"ERROR in createVOGNOde()" << endl;
      }
      temp_node_ptn = all_collect_vectors(my_temp_nodes);
    }
    //create new partition  
    nodes_ptn.resize(Loci::MPI_processes);
    // create node allocation
    long npnts = num_nodes ;
    int node_ivl = npnts / Loci::MPI_processes;
    int node_ivl_rem = npnts % Loci::MPI_processes ;
    int node_accum = 0 ;
    int nodes_base = 0 ;
    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int j = Loci::MPI_processes - i - 1 ;
      int node_accum_update = node_accum + node_ivl + ((j<node_ivl_rem)?1:0) ;
      if(i == Loci::MPI_processes-1) {
        nodes_ptn[i] = interval(nodes_base + node_accum,
                                nodes_base + npnts - 1) ;
      } else {
        nodes_ptn[i] = interval(nodes_base + node_accum,
                                nodes_base + node_accum_update - 1) ;
      }
      node_accum = node_accum_update ;
    }
    
    entitySet new_nodes = nodes_ptn[MPI_rank];
    new_pos.allocate(new_nodes) ;
  
    // Now compute where to send data
    int p = Loci::MPI_processes;
    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;
    
    // Check each processor, find out which sets to send
    for(int i=0;i<p;++i) {
      send_sets[i] = nodes_ptn[i]&my_temp_nodes;
      sequence s= send_sets[i] ;
      send_seqs[i] = s ;
    }
    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;
  
    storeRepP sp =pos_t.Rep();
    storeRepP qcol_rep =new_pos.Rep();
  
  
    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;
    
    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD ) ;
    
    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
		  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
                  MPI_COMM_WORLD ) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    pos_t.allocate(EMPTY) ;
  }


  //create cl, cr and face2node maps from fine_faces
  //re_number the nodes in pos and inner_nodes,
  //and then redistribute them across the processes
  void createVOGFace(int numNodes,
                     const store<Loci::FineFaces> &fine_faces,
                     fact_db & facts,
                     int& numFaces,
                     int& ncells,
                     Map& cl,
                     Map& cr,
                     multiMap& face2node,
                     vector<entitySet>& local_faces,
                     vector<entitySet>& local_cells
                     ){

    //compute numFaces
    constraint my_faces, my_geom_cells;
    my_faces = facts.get_variable("faces");
    my_geom_cells = facts.get_variable("geom_cells");
    entitySet dom = fine_faces.domain() & (*my_faces+*my_geom_cells);
  
    int local_num_face = 0;
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ei++){
      local_num_face += fine_faces[*ei].size();
    }
    numFaces = 0;
    std::vector<int> face_sizes= Loci::all_collect_sizes(local_num_face);
    for(int i =0; i < MPI_processes; i++) numFaces += face_sizes[i];
  
    //get face domains and allocate the maps 
    int face_min = numNodes;
    for(int i =0; i < MPI_rank; i++) face_min += face_sizes[i];
    int face_max = face_min + face_sizes[MPI_rank] -1;
    store<int> count;
    entitySet faces = interval(face_min, face_max);
    cl.allocate(faces);
    cr.allocate(faces);
    count.allocate(faces);
    local_faces.resize(MPI_processes);
    local_faces = all_collect_vectors(faces);
    //fill up the maps cl , cr and count 
    int cell_base = numNodes + numFaces;
    int cell_max = std::numeric_limits<int>::min();
    int cell_min = std::numeric_limits<int>::max();
    entitySet::const_iterator fid = faces.begin();
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ei++){
      for(unsigned int i = 0; i < fine_faces[*ei].size(); i++){
        cl[*fid] = fine_faces[*ei][i][0] + cell_base -1;// -1 finefaces cell index start at 1
        if(fine_faces[*ei][i][1]>=0) cr[*fid] = fine_faces[*ei][i][1] + cell_base -1;
        else  cr[*fid] = fine_faces[*ei][i][1];
        cell_max = max(cell_max,cl[*fid]);
        cell_max = max(cell_max,cr[*fid]);
        cell_min = min(cell_min,cl[*fid]);
        if (cr[*fid]>=0) cell_min = min(cell_min,cr[*fid]);
        count[*fid] = fine_faces[*ei][i].size()-2;
        fid++;
      }
    }
  
    //get cells distribution
    int global_max_cell = cell_max ;
    MPI_Allreduce(&cell_max,&global_max_cell,1,MPI_INT,MPI_MAX,
                  MPI_COMM_WORLD) ;
    int global_min_cell = cell_max ;
    MPI_Allreduce(&cell_min,&global_min_cell,1,MPI_INT,MPI_MIN,
                  MPI_COMM_WORLD) ;
  
  
    ncells = global_max_cell - global_min_cell +1;
  

    local_cells.resize(MPI_processes);
    int cell_ivl = ncells / MPI_processes;
    int cell_ivl_rem = ncells % MPI_processes ;
    int cell_accum = 0 ;

    for(int i = 0; i < Loci::MPI_processes; ++i) {
      int j = MPI_processes - i - 1 ;
      int cell_accum_update = cell_accum + cell_ivl + ((j<cell_ivl_rem)?1:0) ;
    
      if(i == MPI_processes-1) {
        local_cells[i] = interval(cell_base + cell_accum,
                                  cell_base + ncells-1) ;
      } else {
        local_cells[i] = interval(cell_base + cell_accum,
                                  cell_base + cell_accum_update - 1) ;
      }
      cell_accum = cell_accum_update ;
    }
    //fill up face2node  
    face2node.allocate(count);
    fid = faces.begin();
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ei++){
      for(unsigned int i = 0; i < fine_faces[*ei].size(); i++){
        for(int j = 0; j < count[*fid]; j++){
          //vog file node index start with 0
          // face2node[*fid][j] = fine_faces[*ei][i][j+2]-1;
          // if(fine_faces[*ei][i][j+2] < 0) cerr <<"WARNING: negative node index" << endl;
          face2node[*fid][j] = fine_faces[*ei][i][j+2];
        }
        fid++;
      }
    }
    colorMatrix(cl, cr, face2node);
  }
  
}

namespace Loci{
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm);
  
  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm);


  entitySet faceCluster(const multiMap &face2node,
                        const Map &cl, const Map &cr, entitySet faces,
                        vector<unsigned char> &cluster_info,
                        vector<unsigned short> &cluster_sizes) ;
  bool readBCfromVOG(string filename,
                     vector<pair<int,string> > &boundary_ids);
  bool readVolTags(hid_t input_fid,
                   vector<pair<string,Loci::entitySet> > &volDat);
  
  hid_t writeVOGOpen(string filename);
  void writeVOGSurf(hid_t file_id, std::vector<pair<int,string> > surface_ids);
  void writeVOGTag(hid_t output_fid,  vector<pair<string,entitySet> >& volTags);
  void writeVOGClose(hid_t file_id) ;
  void writeVOGNode(hid_t file_id,
                    Loci::storeRepP &pos,
                    const_store<Loci::FineNodes> &inner_nodes);


}



void colorMatrix(Map &cl, Map &cr, multiMap &face2node);
namespace Loci{
  
  //copied from ditribute_io.cc
  hid_t writeVOGOpen(string filename) {
    if(use_parallel_io){    
      hid_t file_id = 0;
      hid_t  acc_plist;
      // open collectively by all processor in MPI_COMM_WORLD,
      acc_plist = Loci::create_faccess_plist(MPI_COMM_WORLD,
                                             Loci::PHDF5_MPI_Info,
                                             Loci::hdf5_const::facc_type); 
      file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,acc_plist) ;
      H5Pclose(acc_plist);
      if(file_id == 0) {
        if(MPI_rank==0) cerr << "unable to open file " << filename << endl ;
        Loci::Abort() ;
      }
      return file_id ;
    }else{
      hid_t file_id = 0 ;
      if(MPI_rank==0) 
        file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
      return file_id ;
    }
  }
  
  
  //copied from FVMGridWriter.cc
  void writeVOGClose(hid_t file_id) {//parallel io included
    if(MPI_rank == 0 || use_parallel_io) H5Fclose(file_id) ;
  }

  //similar to the one in FVMGridWriter.cc, but no modofication to surface_ids
  void writeVOGSurf(hid_t file_id, std::vector<pair<int,string> > surface_ids) {
    hid_t group_id = 0 ;
    if(MPI_rank == 0 || use_parallel_io) {
      if(surface_ids.size() != 0) {
#ifdef H5_USE_16_API
        group_id = H5Gcreate(file_id,"surface_info",0) ;
#else
        group_id = H5Gcreate(file_id,"surface_info",
			     H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
        for(size_t i=0;i<surface_ids.size();++i) {
          hid_t bc_id = 0 ;
#ifdef H5_USE_16_API
          bc_id = H5Gcreate(group_id,surface_ids[i].second.c_str(),0) ;
#else
          bc_id = H5Gcreate(group_id,surface_ids[i].second.c_str(),
			    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
          hsize_t dims = 1 ;
          hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;

#ifdef H5_USE_16_API          
          hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                                   dataspace_id, H5P_DEFAULT) ;
#else
          hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                                   dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
#endif
          H5Awrite(att_id,H5T_NATIVE_INT,&surface_ids[i].first) ;
          H5Aclose(att_id) ;
          H5Gclose(bc_id) ;
        }
        H5Gclose(group_id) ;
      }
    }
  }

  //same as the function in FVMGridReader.cc
  unsigned long readAttributeLong(hid_t group, const char *name) {
    hid_t id_a = H5Aopen_name(group,name) ;
    unsigned long val = 0;
    H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
    H5Aclose(id_a) ;
    return val ;
  }
  //same as the function in FVMGridReader.cc
  bool readVolTags(hid_t input_fid,
                   vector<pair<string,Loci::entitySet> > &volDat) {
    using namespace Loci ;
    /* Save old error handler */
    H5E_auto_t old_func = 0;
    void *old_client_data = 0 ;
#ifdef H5_USE_16_API
    H5Eget_auto(&old_func, &old_client_data);
    H5Eset_auto(NULL, NULL);
#else
    H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);
    H5Eset_auto(H5E_DEFAULT,NULL, NULL);
#endif
    /* Turn off error handling */
    
    vector<pair<string,entitySet> > volTags ;
#ifdef H5_USE_16_API
    hid_t cell_info = H5Gopen(input_fid,"cell_info") ;
#else
    hid_t cell_info = H5Gopen(input_fid,"cell_info",H5P_DEFAULT) ;
#endif
    if(cell_info > 0) {
      vector<string> vol_tag ;
      vector<entitySet> vol_set ;
      vector<int> vol_id ;
      
      hsize_t num_tags = 0 ;
      H5Gget_num_objs(cell_info,&num_tags) ;
      for(hsize_t tg=0;tg<num_tags;++tg) {
        char buf[1024] ;
        memset(buf, '\0', 1024) ;
        H5Gget_objname_by_idx(cell_info,tg,buf,sizeof(buf)) ;
        buf[1023]='\0' ;
        
        string name = string(buf) ;
#ifdef H5_USE_16_API
        hid_t vt_g = H5Gopen(cell_info,buf) ;
#else
        hid_t vt_g = H5Gopen(cell_info,buf,H5P_DEFAULT) ;
#endif
        hid_t id_a = H5Aopen_name(vt_g,"Ident") ;
        int ident ;
        H5Aread(id_a,H5T_NATIVE_INT,&ident) ;
        H5Aclose(id_a) ;
        entitySet dom ;
        HDF5_ReadDomain(vt_g,dom) ;
        vol_tag.push_back(name) ;
        vol_set.push_back(dom) ;
        vol_id.push_back(ident) ;
        H5Gclose(vt_g) ;
      }
      int maxi =0 ;
      for(size_t i=0;i<vol_id.size();++i)
        maxi = max(vol_id[i],maxi) ;
      vector<pair<string,entitySet> > tmp(maxi+1) ;
      volTags.swap(tmp) ;
      for(size_t i=0;i<vol_id.size();++i) {
        volTags[vol_id[i]].first = vol_tag[i] ;
        volTags[vol_id[i]].second = vol_set[i] ;
      }
    } else {
#ifdef H5_USE_16_API
      hid_t file_info = H5Gopen(input_fid,"file_info") ;
#else
      hid_t file_info = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
#endif

      long numCells = readAttributeLong(file_info,"numCells") ;
      volTags.push_back(pair<string,entitySet>
                        (string("Main"),
                         entitySet(interval(0,numCells-1)))) ;
      H5Gclose(file_info) ;
    }
  
    /* Restore previous error handler */
#ifdef H5_USE_16_API
    H5Eset_auto(old_func, old_client_data);
#else
    H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
#endif
    volDat.swap(volTags) ;
    return true ;
  }

}
//same as the function in FVMGridWriter.cc
void writeVOGFace(hid_t file_id, Map &cl, Map &cr, multiMap &face2node) {
  // Compute cell set
  entitySet tmp_cells = cl.image(cl.domain())+cr.image(cr.domain()) ;
  entitySet loc_geom_cells = tmp_cells & interval(0,Loci::UNIVERSE_MAX) ;
  entitySet geom_cells = Loci::all_collect_entitySet(loc_geom_cells) ;
  
  Map tmp_cl, tmp_cr;
  multiMap tmp_face2node;

  
  long long local_num_faces = face2node.domain().size()  ;

  
  long long num_cells = geom_cells.size() ;
  long long num_faces = 0 ;

  // Reduce these variables
  MPI_Allreduce(&local_num_faces,&num_faces,1,MPI_LONG_LONG_INT,
                MPI_SUM,MPI_COMM_WORLD) ;

  hid_t group_id = 0 ;
  if(MPI_rank == 0 || Loci::use_parallel_io) {
#ifdef H5_USE_16_API
    group_id = H5Gopen(file_id,"file_info") ;
#else
    group_id = H5Gopen(file_id,"file_info",H5P_DEFAULT) ;
#endif

    if(MPI_rank == 0) std::cerr<< "num_cells = " << num_cells << endl
                               << "num_faces = " << num_faces << endl ;

    hsize_t dims = 1 ;
    hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
#ifdef H5_USE_16_API
    hid_t att_id = H5Acreate(group_id,"numFaces", H5T_STD_I64BE,
                             dataspace_id, H5P_DEFAULT) ;
#else
    hid_t att_id = H5Acreate(group_id,"numFaces", H5T_STD_I64BE,
                             dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    H5Awrite(att_id,H5T_NATIVE_LLONG,&num_faces) ;
    H5Aclose(att_id) ;
#ifdef H5_USE_16_API
    att_id = H5Acreate(group_id,"numCells", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT) ;
#else
    att_id = H5Acreate(group_id,"numCells", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    H5Awrite(att_id,H5T_NATIVE_LLONG,&num_cells) ;
    H5Aclose(att_id) ;
    H5Gclose(group_id) ;
#ifdef H5_USE_16_API
    group_id = H5Gcreate(file_id,"face_info",0) ;
#else
    group_id = H5Gcreate(file_id,"face_info",
			 H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
  }
  
  entitySet faces = face2node.domain() ;
  vector<pair<pair<int,int>, int> > f_ord(faces.size()) ;
  int i = 0 ;
  // For small number of cells, sort to keep bc groupings
  if(num_cells<100000) {
    FORALL(faces,fc) {
      f_ord[i].first.first = cr[fc] ;
      f_ord[i].first.second = cl[fc] ;
      f_ord[i].second = fc ;
      i++ ;
    } ENDFORALL ;
    std::sort(f_ord.begin(),f_ord.end()) ;
  } else {
    FORALL(faces,fc) {
      f_ord[i].first.first = cl[fc] ;
      f_ord[i].first.second = cr[fc] ;
      f_ord[i].second = fc ;
      i++ ;
    } ENDFORALL ;
  }

  i=0 ;
  store<int> count ;
  count.allocate(faces) ;
  FORALL(faces,fc) {
    int nfc = f_ord[i].second ;
    count[fc] = face2node[nfc].size() ;
    i++ ;
  } ENDFORALL ;
  tmp_face2node.allocate(count) ;
  tmp_cl.allocate(faces) ;
  tmp_cr.allocate(faces) ;
  i=0 ;
  
  int mc = (geom_cells).Min() ;
  // Nodes should be adjusted to start from zero also... for the general case
  FORALL(faces,fc) {
    int nfc = f_ord[i].second ;
    tmp_cl[fc] = cl[nfc]-mc ;
    tmp_cr[fc] = cr[nfc] ;
    if(tmp_cr[fc] >= 0)
      tmp_cr[fc] -= mc ;
    for(int j=0;j<count[fc];++j)
      tmp_face2node[fc][j] = face2node[nfc][j] ;
    i++ ;
  } ENDFORALL ;

  vector<unsigned char> cluster_info ;
  vector<unsigned short> cluster_sizes ;
  while(faces != EMPTY) {
    entitySet fcluster = Loci::faceCluster(tmp_face2node,tmp_cl,tmp_cr,faces,
                                           cluster_info,cluster_sizes) ;
    faces -= fcluster ;
  }

  Loci::writeUnorderedVector(group_id,"cluster_sizes",cluster_sizes) ;
  Loci::writeUnorderedVector(group_id,"cluster_info",cluster_info) ;
  
  
  if(MPI_rank == 0 || Loci::use_parallel_io) {
    H5Gclose(group_id) ;
  }
}


vector<pair<string,entitySet> > getVOGTagFromLocal(const vector<pair<string,entitySet> > &origVolTags,
                                                   // string outfile,
                                                   const_store<int> &cell_offset,
                                                   const_store<int> & num_fine_cells,
                                                   int num_original_nodes,
                                                   int num_original_faces) {
  

  
  vector<pair<string,entitySet> > volTags(origVolTags.size());

  //serial version
  if(Loci::MPI_processes == 1){
    for(unsigned int i = 0; i < origVolTags.size(); i++){
      string name = origVolTags[i].first;
      //transfer entitySet into vector of interval
      entitySet set=origVolTags[i].second;
      int sz = set.num_intervals() ;
      vector<interval> vlist(sz) ;
      for(int j=0;j<sz;++j)vlist[j] = set[j] ;
        
      entitySet new_set = EMPTY;
      for(int j=0;j<sz;++j){
        int start = vlist[j].first+num_original_nodes+num_original_faces;
        int end = vlist[j].second+num_original_nodes+num_original_faces ;
          
        start = cell_offset[start]; //the begin of interval changes to  cell_offset+1
        end = cell_offset[end]+num_fine_cells[end]-1;//the end of interval changes to cell_offset+num_fine_cells
        new_set += interval(start, end);
      }
      volTags[i] = make_pair(name, new_set);
      std::cout << "old tag:  "  << name << " " << set << std::endl;
      std::cout << "new tag:  " << name << " " << new_set << std::endl;
    }
   
    return volTags  ;
  }
       
    
  fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
  Loci::constraint  geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
  Loci::constraint my_entities;
  my_entities = dist->my_entities ;
  //don't know if it's necessray
  entitySet local_geom_cells = (*my_entities)&(*geom_cells);
    
  if(Loci::MPI_processes > 1){
     
    //trasnfer the store to file numbering 
    int offset = 0; 
    store<int> file_num_fine_cells;
    file_num_fine_cells = Loci::Local2FileOrder(num_fine_cells.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;
    offset = 0;
    store<int> file_cell_offset;
    file_cell_offset = Loci::Local2FileOrder(cell_offset.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;

    int coffset = 0; 
    std::vector<int> local_cell_sizes;
    int num_local_cells = file_cell_offset.domain().size();
    local_cell_sizes = Loci::all_collect_sizes(num_local_cells);
    for(int i = 0; i < MPI_rank; i++){
      coffset += local_cell_sizes[i];
    }  
     
      
    //process 0 broadcast the entitySet to all
    int buf_size = 0;
    if(MPI_rank==0){
      for(unsigned int i = 0; i < origVolTags.size(); i++){
        //transfer entitySet into vector of interval
        entitySet set=origVolTags[i].second;
        int sz = set.num_intervals() ;
        buf_size += sz;
      }
    }
    MPI_Bcast(&buf_size,1,MPI_INT, 0, MPI_COMM_WORLD) ;
    vector<int> buffer(2*buf_size);
    if(MPI_rank==0){
      int pnt = 0;
      for(unsigned int i = 0; i < origVolTags.size(); i++){
        //transfer entitySet into vector of interval
        entitySet set=origVolTags[i].second;
        int sz = set.num_intervals() ;
        for(int j=0;j<sz;++j){
          buffer[pnt++] =  set[j].first;
          buffer[pnt++] = set[j].second;
        }
        
      }
    }
    MPI_Bcast(&buffer[0],2*buf_size,MPI_INT,0,MPI_COMM_WORLD) ;
    //each process re-assign the intervals
    entitySet domain = file_cell_offset.domain();
     
    for (int i = 0; i < buf_size; i++){
      int start = buffer[2*i]-coffset;
      int end = buffer[2*i+1]-coffset;
      if(start >=0 && domain.inSet(start)) buffer[2*i] = file_cell_offset[start];
      else buffer[2*i] = 0;
      if(end >=0 && domain.inSet(end)) buffer[2*i+1] = file_cell_offset[end]+file_num_fine_cells[end]-1;
      else buffer[2*i+1] = 0;
       
    }
    //process 0 receive the values in buffer from all the other processes, add up the values from each process. 
    if(MPI_rank == 0){
       
      for(int pi = 1; pi <MPI_processes; pi++){
        //send a flag
        int flag = 0 ;
        MPI_Send(&flag,1,MPI_INT,pi,6,MPI_COMM_WORLD) ;
        vector<int> recv_buf(2*buf_size);
        MPI_Status mstat ;
        MPI_Recv(&recv_buf[0],sizeof(int)*2*buf_size,MPI_BYTE,pi,7,MPI_COMM_WORLD,
                 &mstat) ;
        for(int i = 0; i < 2*buf_size; i++)buffer[i] += recv_buf[i];
      }
    }else{
      int flag = 0;
      MPI_Status mstat ;
      MPI_Recv(&flag,1,MPI_INT,0,6,MPI_COMM_WORLD,&mstat) ;
      MPI_Send(&buffer[0],sizeof(int)*2*buf_size,MPI_BYTE,0,7,MPI_COMM_WORLD);
    }
    //process 0 re-organize the data in buffer into volTags and write it out
    if(MPI_rank == 0){
      int pnt = 0;
      for(unsigned int i = 0; i < origVolTags.size(); i++){
        //transfer entitySet into vector of interval
        string name = origVolTags[i].first;
        entitySet set=origVolTags[i].second;
        int sz = set.num_intervals() ;
        entitySet new_set = EMPTY;
        for(int j = 0; j < sz; j++){ 
          int start = buffer[pnt++]; 
          int end = buffer[pnt++];
          new_set += interval(start, end);
        }
        volTags[i] = make_pair(name, new_set);
        std::cout << "old tag:  "  << name << " " << set << std::endl;
        std::cout << "new tag:  " << name << " " << new_set << std::endl;
      }
    
    }

    //broadcast volTags to other processes
    int nvtags = volTags.size() ;
    MPI_Bcast(&nvtags,1,MPI_INT,0,MPI_COMM_WORLD) ;
    if(MPI_rank != 0) volTags.clear();
    for(int i=0;i<nvtags;++i) {
      int sz = 0 ;
      if(MPI_rank == 0) sz = volTags[i].first.size() ;
      MPI_Bcast(&sz,1,MPI_INT,0,MPI_COMM_WORLD) ;
      char *buf = new char[sz+1] ;
      buf[sz] = '\0' ;
      if(MPI_rank == 0) strcpy(buf,volTags[i].first.c_str()) ;
      MPI_Bcast(buf,sz,MPI_CHAR,0,MPI_COMM_WORLD) ;
      string name = string(buf) ;
      delete[] buf ;
      int nivals = 0 ;
      if(MPI_rank == 0) nivals = volTags[i].second.num_intervals() ;
      MPI_Bcast(&nivals,1,MPI_INT,0,MPI_COMM_WORLD) ;

      int *ibuf = new int[nivals*2] ;
      if(MPI_rank == 0) 
        for(int j=0;j<nivals;++j) {
          ibuf[j*2]= volTags[i].second[j].first ;
          ibuf[j*2+1]= volTags[i].second[j].second ;
        }
      MPI_Bcast(ibuf,nivals*2,MPI_INT,0,MPI_COMM_WORLD) ;
      entitySet set ;
      for(int j=0;j<nivals;++j) 
        set += interval(ibuf[j*2],ibuf[j*2+1]) ;
      if(MPI_rank != 0) 
        volTags.push_back(pair<string,entitySet>(name,set)) ;
      delete[] ibuf ;
    }
     
  }
  return volTags; 
}






