#ifndef GRID_INTERFACE_H
#define GRID_INTERFACE_H
#include <store_rep.h>

namespace Loci {
  void parallelClassifyCell(fact_db &facts) ;

  void createVOGNode(store<vector3d<double> > &new_pos,
                     const store<Loci::FineNodes> &inner_nodes_cell,
                     const store<Loci::FineNodes> &inner_nodes_face,
                     const store<Loci::FineNodes> &inner_nodes_edge,
                     int& num_nodes,
                     fact_db & facts,//in global numbering
                     vector<entitySet>& nodes_ptn) ;

  void createVOGFace(int numNodes,
                     const store<Loci::FineFaces> &fine_faces_cell,
                     const store<Loci::FineFaces> &fine_faces,
                     fact_db & facts,
                     int& numFaces,
                     int& ncells,
                     Map& cl,
                     Map& cr,
                     multiMap& face2node,
                     vector<entitySet>& local_faces,
                     vector<entitySet>& local_cells) ;


  // Map a store from the parent cells to the child cells of the adapted mesh
  storeRepP mapCellPartitionWeights(storeRepP wptr,
				    storeRepP c2p_ptr,
				    entitySet localcelldom) ;

  class refinedGridData : public Loci::CPTR_type {
  public:
    //global variables that pass data from refinemesh to chem
    Map new_cl;
    Map new_cr;
    multiMap new_face2node;
    store<vector3d<double> > new_pos;
    vector<entitySet> local_nodes;
    vector<entitySet> local_faces;
    vector<entitySet> local_cells;
    vector<pair<int,string> > boundary_ids;
    vector<pair<string,entitySet> > volTags;
  } ;

  void initializeGridFromPlan(Loci::CPTR<refinedGridData> &gridDataP,
			      int &level,
			      rule_db &refmesh_rdb,
			      string casename, string weightfile,
			      string restartplanfile) ;

  void onlineRefineMesh(Loci::CPTR<refinedGridData> &gridDataP,
			rule_db &refmesh_rdb,
			int adaptmode,
			int level,
			storeRepP tags,
			string casename  ) ;


  storeRepP getC2PGlobal(fact_db &facts) ;
  
  // Holds Data Structures for processing AMR interpolation
  struct AMRrefinementMapping {
    entitySet geom_cells_global ;
    // parent cell numbers 
    store<int> parent ;
    // local numbering of parent to child mapping in local numbering
    // for scattering communication
    multiStore<int> parent2child_l ;
    // List of refined parent cells
    entitySet refinedCells ;
    // Allocation set for gradient input gather
    entitySet gradAccessSet ;
    // Communication schedule needed to compute gradients of parent nodes
    gatherCommSchedule gradientComm ;
    // stencils for computing gradient of parent cell data
    multiStore<int> gradCellStencil ;
    // weights for computing graidents from stencil
    multiStore<vector3d<float> > stencilWeights ;
    // deltas used to compute stencil
    multiStore<vector3d<float> > grad_dvs ;
    // Vectors from parent centroid to cell centroids
    multiStore<vector3d<double> > child_dvs ;
    // Communication schedule gathering refined cells to parent processor
    gatherCommSchedule  refineCellComm ;
    // child -> local numbering of parent data
    vector<pair<int,int> > c2lp ;
    // Set of child cells that have direct mapping to parents
    entitySet directMapCells ;
    // Communicationschedule gathering parent cell data to children
    gatherCommSchedule directMapComm ;
    // Weights for volume weighted average for coarsened cells
    vector<float>  directWeights ;
  } ;

  void AMRprocessRefinementMapping(struct AMRrefinementMapping &mappings,
                                   const store<pair<int,int> > &c2pg,
                                   const store<float> &volw,
                                   multiStore<int> &gradCells,
                                   multiStore<vector3d<float> > &deltas,
                                   const_store<vector3d<double> > &cell_center,
                                   const_store<double> &vol,
                                   fact_db &facts) ;
  
  template<class T>
  void AMRinterpolateData(AMRrefinementMapping &mapping,
                          storeVec<T> &parent_data,
                          storeVec<T> &child_data) {
    //########################################################################
    //
    // First interpolate to cells that resulted from splitting parent cells
    //
    // Step 1, compute gradients using parent data for parents that form
    // split cells
    // Step 2, compute limiter for these gradients
    // Step 3, compute limited second order interpolation to child cells
    // Step 4, distribute second order reconstruction from parent to child cells
    // Step 5, collect direct mapped cells and use volume averaging to
    //         interpolate from many parents to child for coarsened cells
    // Note, if the variables given are conservative, the extrapolation is
    // conservative
    //
    const int vs = parent_data.vecSize() ;
    // gather parent data needed to compute parent gradients
    storeVec<T> grad_data ;
    mapping.gradientComm.gatherData(grad_data,parent_data) ;
    
    storeVec<T> refineCell_data ;
    refineCell_data.setVecSize(vs) ;
    int refcells = mapping.refinedCells.size() ;
#ifdef DEBUG
    Loci::debugout  << "#refcells interpolations = " << refcells << endl ;
#endif
    entitySet refineCellDomain ;
    if(refcells>0)
      refineCellDomain = interval(0,refcells-1) ;
    refineCell_data.allocate(refineCellDomain) ;
    int vsz = parent_data.vecSize() ;
    FORALL(mapping.parent.domain(),ii) {
      // get stencil size
      const int ssz = mapping.gradCellStencil[ii].getSize() ;
      int parent = mapping.parent[ii] ;
      for(int k=0;k<vsz;++k) { // loop over vector
        // Compute gradient
        vector3d<float> gradk(0,0,0) ;
        T Xcc = parent_data[parent][k] ;
        T max_val = Xcc ;
        T min_val = max_val ;
        for(int j=0;j<ssz;++j) {
          const int cc = mapping.gradCellStencil[ii][j] ;
          T val = grad_data[cc][k] ;
          max_val = max(max_val,val) ;
          min_val = min(min_val,val) ;
          gradk += mapping.stencilWeights[ii][j]*(val-Xcc) ;
        }
        T limi = 1.0 ;
        // Apply limiter to source points
        for(int j=0;j<ssz;++j) {
          const int cc = mapping.gradCellStencil[ii][j] ;
          T qdif = dot(gradk,mapping.grad_dvs[ii][j]) ;
          if(qdif > 0)
            limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
          if(qdif < 0)
            limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
        }
        // Apply limiter to target points
        for(int j=0;j<mapping.child_dvs[ii].getSize();++j) {
          T qdif =  dot(vector3d<double>(gradk),
                        mapping.child_dvs[ii][j]) ;
          if(qdif > 0)
            limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
          if(qdif < 0)
            limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
        }

        // reduce gradient by limiter factor
        gradk *= limi ;
        // compute child cell values
        for(int j=0;j<mapping.child_dvs[ii].getSize();++j) {
          T val = Xcc + dot(vector3d<double>(gradk),
                            mapping.child_dvs[ii][j]) ;
          refineCell_data[mapping.parent2child_l[ii][j]][k] = val ;
        }
      }
    } ENDFORALL ;

    mapping.refineCellComm.scatterData(refineCell_data,child_data) ;


    // direct map cells with averaging for coarse cells
    storeVec<T> mapCell_data ;
    mapping.directMapComm.gatherData(mapCell_data,parent_data) ;
    for(size_t i=0;i<mapping.c2lp.size();++i) {
      int child = mapping.c2lp[i].first ;
      float w = mapping.directWeights[i] ;
      for(int k=0;k<vs;++k)
        child_data[child][k] = w*mapCell_data[mapping.c2lp[i].second][k] ;
      size_t j = i+1 ;
      for(;j<mapping.c2lp.size()&&child == mapping.c2lp[j].first;++j) {
        float w = mapping.directWeights[j] ;
        for(int k=0;k<vs;++k)
          child_data[child][k] +=  w*mapCell_data[mapping.c2lp[j].second][k] ;
      }
      int cnt = j-i ;
      i+= cnt-1 ;
    }
  }
  
  template<class T>
  void AMRinterpolateData(AMRrefinementMapping &mapping,
                          store<T> &parent_data,
                          store<T> &child_data) {
    //########################################################################
    //
    // First interpolate to cells that resulted from splitting parent cells
    //
    // Step 1, compute gradients using parent data for parents that form
    // split cells
    // Step 2, compute limiter for these gradients
    // Step 3, compute limited second order interpolation to child cells
    // Step 4, distribute second order reconstruction from parent to child cells
    // Step 5, collect direct mapped cells and use volume averaging to
    //         interpolate from many parents to child for coarsened cells
    // Note, if the variables given are conservative, the extrapolation is
    // conservative
    //
    // gather parent data needed to compute parent gradients
    store<T> grad_data ;
    mapping.gradientComm.gatherData(grad_data,parent_data) ;
    
    store<T> refineCell_data ;
    int refcells = mapping.refinedCells.size() ;
#ifdef DEBUG
    Loci::debugout  << "#refcells interpolations = " << refcells << endl ;
#endif
    entitySet refineCellDomain ;
    if(refcells>0)
      refineCellDomain = interval(0,refcells-1) ;
    refineCell_data.allocate(refineCellDomain) ;
    FORALL(mapping.parent.domain(),ii) {
      // get stencil size
      const int ssz = mapping.gradCellStencil[ii].getSize() ;
      int parent = mapping.parent[ii] ;
      // Compute gradient
      vector3d<float> gradk(0,0,0) ;
      T Xcc = parent_data[parent] ;
      T max_val = Xcc ;
      T min_val = max_val ;
      for(int j=0;j<ssz;++j) {
        const int cc = mapping.gradCellStencil[ii][j] ;
        T val = grad_data[cc] ;
        max_val = max(max_val,val) ;
        min_val = min(min_val,val) ;
        gradk += mapping.stencilWeights[ii][j]*(val-Xcc) ;
      }
      T limi = 1.0 ;
      // Apply limiter to source points
      for(int j=0;j<ssz;++j) {
        const int cc = mapping.gradCellStencil[ii][j] ;
        T qdif = dot(gradk,mapping.grad_dvs[ii][j]) ;
        if(qdif > 0)
          limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
        if(qdif < 0)
          limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
      }
      // Apply limiter to target points
      for(int j=0;j<mapping.child_dvs[ii].getSize();++j) {
        T qdif =  dot(vector3d<double>(gradk),
                      mapping.child_dvs[ii][j]) ;
        if(qdif > 0)
          limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
        if(qdif < 0)
          limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
      }

      // reduce gradient by limiter factor
      gradk *= limi ;
      // compute child cell values
      for(int j=0;j<mapping.child_dvs[ii].getSize();++j) {
        T val = Xcc + dot(vector3d<double>(gradk),
                          mapping.child_dvs[ii][j]) ;
        refineCell_data[mapping.parent2child_l[ii][j]] = val ;
      }
    } ENDFORALL ;

    mapping.refineCellComm.scatterData(refineCell_data,child_data) ;


    // direct map cells with averaging for coarse cells
    store<T> mapCell_data ;
    mapping.directMapComm.gatherData(mapCell_data,parent_data) ;
    for(size_t i=0;i<mapping.c2lp.size();++i) {
      int child = mapping.c2lp[i].first ;
      float w = mapping.directWeights[i] ;
      child_data[child] = w*mapCell_data[mapping.c2lp[i].second] ;
      size_t j = i+1 ;
      for(;j<mapping.c2lp.size()&&child == mapping.c2lp[j].first;++j) {
        float w = mapping.directWeights[j] ;
        child_data[child] +=  w*mapCell_data[mapping.c2lp[j].second] ;
      }
      int cnt = j-i ;
      i+= cnt-1 ;
    }
  }
  

    template <class T> void AMRinterpolateStore(store<T> &qout,
                                              store<T> &qin,
                                              store<pair<int,int> > &c2pg,
                                              entitySet dom) {
    
    store<float> weights ;
    weights.allocate(dom) ;
    zeroStore(weights,dom) ;
    zeroStore(qout,dom) ;

    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;

    if(MPI_processes == 1) {
      entitySet parentSet ;
      entitySet targetSet ;
      entitySet domc2p = c2pg.domain() ;

      FORALL(domc2p,ii) {
	parentSet += c2pg[ii].second ;
	targetSet += c2pg[ii].first ;
      } ENDFORALL ;

      int parent_offset = qin.domain().Min() - parentSet.Min() ;
      int local_offset = dom.Min()-targetSet.Min() ;
      int sz = domc2p.size() ;
      vector<pair<int,int> > c2p(sz) ;
      int cnt = 0 ;
      FORALL(domc2p,ii) {
	c2p[cnt].first = c2pg[ii].second+parent_offset ;
	c2p[cnt].second = c2pg[ii].first+local_offset ;
	cnt++ ;;
      } ENDFORALL ;

      sort(c2p.begin(),c2p.end()) ;

      for(int ii=0;ii<sz;) {
	int pid = c2p[ii].first ;
	int j = 1 ;
	while(((ii+j) < sz) && (pid == c2p[ii+j].first))
	  j++ ;
	if(j == 1) {
	  int tid = c2p[ii].second ;
	  double w = 1.0 ;
	  qout[tid] += qin[pid] ;
	  weights[tid] += w ;
	} else {
	  for(int k=ii;k<ii+j;++k) {
	    int tid = c2p[k].second ;
	    qout[tid] += qin[pid] ;
	    weights[tid] += 1. ;
	  }
	}
	
	ii += j ;
      }
      FORALL(dom,ii) {
	if(weights[ii] > 0) {
	  double w = 1./double(weights[ii]) ;
          qout[ii] *= w ;
	}
      } ENDFORALL ;
    } else { // Parallel interpolation
      int p = Loci::MPI_processes ;
     
      entitySet domc2p = c2pg.domain() ;

      int minParent = std::numeric_limits<int>::max() ;
      if(domc2p != EMPTY) {
	int first = domc2p.Min() ;
	minParent = c2pg[first].second ;
	FORALL(domc2p,ii) {
	  minParent = min(minParent,c2pg[ii].second) ;
	} ENDFORALL ;
      }
      int parentOffset = minParent ;
      MPI_Allreduce(&minParent,&parentOffset,1,MPI_INT,MPI_MIN,
		    MPI_COMM_WORLD) ;
      int sz = domc2p.size() ;
      vector<pair<int,int> > p2c(sz) ;
      int cnt = 0 ;
      FORALL(domc2p,ii) {
	p2c[cnt].first = c2pg[ii].second-parentOffset ;
	p2c[cnt].second = c2pg[ii].first ;
	cnt++ ;;
      } ENDFORALL ;

      sort(p2c.begin(),p2c.end()) ;
      
      // Now distribute c2p to processors in parent file ordering
      int prevdomsz = qin.domain().size() ;
      vector<int> parentsizes(p) ;
      MPI_Allgather(&prevdomsz,1,MPI_INT,&parentsizes[0],1,MPI_INT,
		    MPI_COMM_WORLD) ;
      vector<int> parentoffsets(p+1,0) ;
      for(int i=0;i<p;++i)
	parentoffsets[i+1] = parentoffsets[i]+parentsizes[i] ;
      
      vector<pair<int,int> > splits(p-1) ;
      for(int i=0;i<p-1;++i)
	splits[i] = pair<int,int>(parentoffsets[i+1],-1) ;

      Loci::parSplitSort(p2c,splits,MPI_COMM_WORLD) ;

      int p2csz = p2c.size() ;
      // convert parent to local numbering and collect 
      // id's of target cells in new mesh
      vector<int> tlist(p2csz) ;
      for(int i =0;i<p2csz;++i) {
	tlist[i] = p2c[i].second ;
      }
      // identify list of target points to which we are interpolating
      sort(tlist.begin(),tlist.end()) ;
      vector<int>::iterator l = std::unique(tlist.begin(),tlist.end()) ;
      tlist.resize(std::distance(tlist.begin(),l)) ;
      // Grab data needed for interpolation from target mesh 
      // First reorder into file numbering the cellcenter and vol
      // data which will be used to compute the interpolation

      Map l2g ;
      l2g = dist->l2g.Rep() ;
      Loci::entitySet globaldom ;
      FORALL(dom,ii) {
	globaldom += l2g[ii] ;
      } ENDFORALL ;

      Map g2l ;
      g2l.allocate(globaldom) ;
      FORALL(dom,ii) {
	g2l[l2g[ii]] = ii ;
      } ENDFORALL ;
      
      vector<int> targetsizes(p) ;
      int tlsize = globaldom.size() ;
      MPI_Allgather(&tlsize,1,MPI_INT,&targetsizes[0],1,MPI_INT,
		    MPI_COMM_WORLD) ;
      vector<int> toffsets(p+1,0) ;
      int tmin = globaldom.Min() ;
      MPI_Allgather(&tmin,1,MPI_INT,&toffsets[0],1,MPI_INT,MPI_COMM_WORLD) ;
      toffsets[p] = toffsets[p-1]+targetsizes[p-1] ;

      // Now compute which processors to communicate with
      vector<int> tlistoffsets(p+1,0) ;
      int tsz = tlist.size() ;
      tlistoffsets[p] = tsz ;
      cnt = 1 ;
      for(int i=0;i<tsz;++i) {
	while(cnt != p && tlist[i] >= toffsets[cnt]) {
	  tlistoffsets[cnt] = i ;
	  cnt++ ;
	}
	if(cnt == p)
	  break ;
      }
      if(cnt != p) { // handle end case
	for(int i=cnt;i<p;++i)
	  tlistoffsets[i] = tsz ;
      }
      // Now send the requests
      vector<int> tlistsizes(p) ;
      for(int i=0;i<p;++i)
	tlistsizes[i] = tlistoffsets[i+1]-tlistoffsets[i] ;
      vector<int> rlistsizes(p) ;
      MPI_Alltoall(&tlistsizes[0],1,MPI_INT,
		   &rlistsizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
      vector<int> rlistoffsets(p+1,0) ;
      for(int i=0;i<p;++i)
	rlistoffsets[i+1] = rlistoffsets[i]+rlistsizes[i] ;
      
      vector<int> recvlists(rlistoffsets[p]) ;
      vector<MPI_Request> req_queue ;
      for(int i=0;i<p;++i) {
	if(rlistsizes[i] > 0) {
	  MPI_Request tmp ;
	  MPI_Irecv(&recvlists[rlistoffsets[i]],rlistsizes[i],MPI_INT,i,99,
		    MPI_COMM_WORLD,&tmp) ;
	  req_queue.push_back(tmp) ;
	}
      }
      for(int i=0;i<p;++i) {
	if(tlistsizes[i] > 0) {
	  MPI_Send(&tlist[tlistoffsets[i]],tlistsizes[i],MPI_INT,i,99,
		   MPI_COMM_WORLD) ;
	}
      }
      if(req_queue.size() > 0) {
	vector<MPI_Status> stat_queue(req_queue.size()) ;
	MPI_Waitall(req_queue.size(),&req_queue[0],&stat_queue[0]) ;
	req_queue.clear() ;
      }

      int rsz = recvlists.size() ;
      // correct indexes to local numbering
      
      // Now translate the p2c map target to the local numbering
      std::map<int,int> f2c ;
      for(int i=0;i<tsz;++i)
	f2c[tlist[i]] = i ;
      
      for(int i=0;i<p2csz;++i) {
	std::map<int,int>::const_iterator mi = f2c.find(p2c[i].second) ;
	if(mi == f2c.end()) {
	  cerr << "unable to remap p2c to local number" << endl ;
	  Loci::Abort() ;
	}
	
	p2c[i].second = mi->second ;
      }
      
      entitySet locdom = interval(0,tsz-1) ;
      store<T> qic_loc ;
      qic_loc.allocate(locdom) ;
      store<double> weights ;
      weights.allocate(locdom) ;

      for(int i=0;i<tsz;++i) {
	weights[i] = 0 ;
        setZero(qic_loc[i]) ;
      }
	
      for(int ii=0;ii<p2csz;) {
	int pid = p2c[ii].first ;

	int j = 1 ;
	while(((ii+j) < p2csz) && (pid == p2c[ii+j].first))
	  j++ ;
	if(j == 1) {
	  int tid = p2c[ii].second ;
	  double w = 1.0 ;

          qic_loc[tid] += qin[pid] ;
	  weights[tid] += w ;
	} else {
	  for(int k=ii;k<ii+j;++k) {
	    int tid = p2c[k].second ;
            qic_loc[tid] += qin[pid] ;
	    weights[tid] += 1. ;
	  }
	}
	
	ii += j ;
      }

      // Now return interpolated versions to new file numbering
      vector<T> senddata(tsz) ;
      vector<double> senddataw(tsz,0) ;
      vector<T> recvdata(rsz) ;
      vector<double> recvdataw(rsz,1) ;
      // Fill send data buffer
      for(int i=0;i<tsz;++i) {
        senddata[i] = qic_loc[i] ;
	senddataw[i] = weights[i] ;
      }
      for(int i=0;i<p;++i) {
	if(rlistsizes[i] > 0) {
	  MPI_Request tmp ;
	  MPI_Irecv(&recvdata[rlistoffsets[i]],rlistsizes[i]*sizeof(T),
		    MPI_BYTE,i,99,MPI_COMM_WORLD,&tmp) ;
	  req_queue.push_back(tmp) ;
	  MPI_Irecv(&recvdataw[rlistoffsets[i]],rlistsizes[i],
		    MPI_DOUBLE,i,98,MPI_COMM_WORLD,&tmp) ;
	  req_queue.push_back(tmp) ;
	}
      }
      for(int i=0;i<p;++i) {
	if(tlistsizes[i] > 0) {
	  MPI_Send(&senddata[tlistoffsets[i]],tlistsizes[i]*sizeof(T),
		   MPI_BYTE,i,99,MPI_COMM_WORLD) ;
	  MPI_Send(&senddataw[tlistoffsets[i]],tlistsizes[i],
		   MPI_DOUBLE,i,98,MPI_COMM_WORLD) ;
	}
      }
      if(req_queue.size() > 0) {
	vector<MPI_Status> stat_queue(req_queue.size()) ;
	MPI_Waitall(req_queue.size(),&req_queue[0],&stat_queue[0]) ;
	req_queue.clear() ;
      }

      store<double> q_ic_fn ;
      store<double> wt_fn ;

      entitySet dom2 = dom&qout.domain() ;
      WARN(dom!=dom2) ;
      wt_fn.allocate(dom) ;
      // sum contributions from processors
      FORALL(dom,ii) {
        setZero(qout[ii]) ;
	wt_fn[ii] = 0 ;
      } ENDFORALL ;
      for(int i=0;i<rsz;++i) {
	int id =g2l[recvlists[i]] ;

        qout[id] += recvdata[i] ;
	wt_fn[id] += recvdataw[i] ;
      }

      // re-normalize
      FORALL(dom2,ii) {
	double rw = 1./(wt_fn[ii]) ;
        qout[ii] *= rw;
      } ENDFORALL ;

    }
  }

  template <class T> void AMRinterpolateStore(storeVec<T> &qvecout,
                                              storeVec<T> &qvecin,
                                              store<pair<int,int> > &c2pg,
                                              entitySet dom, const int vs ) {

    cerr << "AMRinterpolateStoreVec" << endl ;
    store<float> weights ;
    weights.allocate(dom) ;
    zeroStore(weights,dom) ;
    zeroStore(qvecout,dom) ;

    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;

    if(MPI_processes == 1) {
      entitySet parentSet ;
      entitySet targetSet ;
      entitySet domc2p = c2pg.domain() ;

      FORALL(domc2p,ii) {
	parentSet += c2pg[ii].second ;
	targetSet += c2pg[ii].first ;
      } ENDFORALL ;

      int parent_offset = qvecin.domain().Min() - parentSet.Min() ;
      int local_offset = dom.Min()-targetSet.Min() ;
      int sz = domc2p.size() ;
      vector<pair<int,int> > c2p(sz) ;
      int cnt = 0 ;
      FORALL(domc2p,ii) {
	c2p[cnt].first = c2pg[ii].second+parent_offset ;
	c2p[cnt].second = c2pg[ii].first+local_offset ;
	cnt++ ;;
      } ENDFORALL ;

      sort(c2p.begin(),c2p.end()) ;

      for(int ii=0;ii<sz;) {
	int pid = c2p[ii].first ;
	int j = 1 ;
	while(((ii+j) < sz) && (pid == c2p[ii+j].first))
	  j++ ;
	if(j == 1) {
	  int tid = c2p[ii].second ;
	  double w = 1.0 ;
	  qvecout[tid] += qvecin[pid] ;
	  weights[tid] += w ;
	} else {
	  for(int k=ii;k<ii+j;++k) {
	    int tid = c2p[k].second ;
	    qvecout[tid] += qvecin[pid] ;
	    weights[tid] += 1. ;
	  }
	}
	
	ii += j ;
      }
      FORALL(dom,ii) {
	if(weights[ii] > 0) {
	  double w = 1./double(weights[ii]) ;
	  for(int j=0;j<vs;++j)
	    qvecout[ii][j] *= w ;
	}
      } ENDFORALL ;
    } else { // Parallel interpolation
      int p = Loci::MPI_processes ;
     
      entitySet domc2p = c2pg.domain() ;

      int minParent = std::numeric_limits<int>::max() ;
      if(domc2p != EMPTY) {
	int first = domc2p.Min() ;
	minParent = c2pg[first].second ;
	FORALL(domc2p,ii) {
	  minParent = min(minParent,c2pg[ii].second) ;
	} ENDFORALL ;
      }
      int parentOffset = minParent ;
      MPI_Allreduce(&minParent,&parentOffset,1,MPI_INT,MPI_MIN,
		    MPI_COMM_WORLD) ;
      int sz = domc2p.size() ;
      vector<pair<int,int> > p2c(sz) ;
      int cnt = 0 ;
      FORALL(domc2p,ii) {
	p2c[cnt].first = c2pg[ii].second-parentOffset ;
	p2c[cnt].second = c2pg[ii].first ;
	cnt++ ;;
      } ENDFORALL ;

      sort(p2c.begin(),p2c.end()) ;
      
      // Now distribute c2p to processors in parent file ordering
      int prevdomsz = qvecin.domain().size() ;
      vector<int> parentsizes(p) ;
      MPI_Allgather(&prevdomsz,1,MPI_INT,&parentsizes[0],1,MPI_INT,
		    MPI_COMM_WORLD) ;
      vector<int> parentoffsets(p+1,0) ;
      for(int i=0;i<p;++i)
	parentoffsets[i+1] = parentoffsets[i]+parentsizes[i] ;
      
      vector<pair<int,int> > splits(p-1) ;
      for(int i=0;i<p-1;++i)
	splits[i] = pair<int,int>(parentoffsets[i+1],-1) ;

      Loci::parSplitSort(p2c,splits,MPI_COMM_WORLD) ;

      int p2csz = p2c.size() ;
      // convert parent to local numbering and collect 
      // id's of target cells in new mesh
      vector<int> tlist(p2csz) ;
      for(int i =0;i<p2csz;++i) {
	tlist[i] = p2c[i].second ;
      }
      // identify list of target points to which we are interpolating
      sort(tlist.begin(),tlist.end()) ;
      vector<int>::iterator l = std::unique(tlist.begin(),tlist.end()) ;
      tlist.resize(std::distance(tlist.begin(),l)) ;
      // Grab data needed for interpolation from target mesh 
      // First reorder into file numbering the cellcenter and vol
      // data which will be used to compute the interpolation

      Map l2g ;
      l2g = dist->l2g.Rep() ;
      Loci::entitySet globaldom ;
      FORALL(dom,ii) {
	globaldom += l2g[ii] ;
      } ENDFORALL ;

      Map g2l ;
      g2l.allocate(globaldom) ;
      FORALL(dom,ii) {
	g2l[l2g[ii]] = ii ;
      } ENDFORALL ;
      
      vector<int> targetsizes(p) ;
      int tlsize = globaldom.size() ;
      MPI_Allgather(&tlsize,1,MPI_INT,&targetsizes[0],1,MPI_INT,
		    MPI_COMM_WORLD) ;
      vector<int> toffsets(p+1,0) ;
      int tmin = globaldom.Min() ;
      MPI_Allgather(&tmin,1,MPI_INT,&toffsets[0],1,MPI_INT,MPI_COMM_WORLD) ;
      toffsets[p] = toffsets[p-1]+targetsizes[p-1] ;

      // Now compute which processors to communicate with
      vector<int> tlistoffsets(p+1,0) ;
      int tsz = tlist.size() ;
      tlistoffsets[p] = tsz ;
      cnt = 1 ;
      for(int i=0;i<tsz;++i) {
	while(cnt != p && tlist[i] >= toffsets[cnt]) {
	  tlistoffsets[cnt] = i ;
	  cnt++ ;
	}
	if(cnt == p)
	  break ;
      }
      if(cnt != p) { // handle end case
	for(int i=cnt;i<p;++i)
	  tlistoffsets[i] = tsz ;
      }
      // Now send the requests
      vector<int> tlistsizes(p) ;
      for(int i=0;i<p;++i)
	tlistsizes[i] = tlistoffsets[i+1]-tlistoffsets[i] ;
      vector<int> rlistsizes(p) ;
      MPI_Alltoall(&tlistsizes[0],1,MPI_INT,
		   &rlistsizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
      vector<int> rlistoffsets(p+1,0) ;
      for(int i=0;i<p;++i)
	rlistoffsets[i+1] = rlistoffsets[i]+rlistsizes[i] ;
      
      vector<int> recvlists(rlistoffsets[p]) ;
      vector<MPI_Request> req_queue ;
      for(int i=0;i<p;++i) {
	if(rlistsizes[i] > 0) {
	  MPI_Request tmp ;
	  MPI_Irecv(&recvlists[rlistoffsets[i]],rlistsizes[i],MPI_INT,i,99,
		    MPI_COMM_WORLD,&tmp) ;
	  req_queue.push_back(tmp) ;
	}
      }
      for(int i=0;i<p;++i) {
	if(tlistsizes[i] > 0) {
	  MPI_Send(&tlist[tlistoffsets[i]],tlistsizes[i],MPI_INT,i,99,
		   MPI_COMM_WORLD) ;
	}
      }
      if(req_queue.size() > 0) {
	vector<MPI_Status> stat_queue(req_queue.size()) ;
	MPI_Waitall(req_queue.size(),&req_queue[0],&stat_queue[0]) ;
	req_queue.clear() ;
      }

      int rsz = recvlists.size() ;
      // correct indexes to local numbering
      
      // Now translate the p2c map target to the local numbering
      std::map<int,int> f2c ;
      for(int i=0;i<tsz;++i)
	f2c[tlist[i]] = i ;
      
      for(int i=0;i<p2csz;++i) {
	std::map<int,int>::const_iterator mi = f2c.find(p2c[i].second) ;
	if(mi == f2c.end()) {
	  cerr << "unable to remap p2c to local number" << endl ;
	  Loci::Abort() ;
	}
	
	p2c[i].second = mi->second ;
      }
      
      entitySet locdom = interval(0,tsz-1) ;
      storeVec<T> qic_loc ;
      qic_loc.allocate(locdom) ;
      qic_loc.setVecSize(vs) ;
      store<double> weights ;
      weights.allocate(locdom) ;

      for(int i=0;i<tsz;++i) {
	weights[i] = 0 ;
	for(int j=0;j<vs;++j)
	  setZero(qic_loc[i][j]) ;
      }
	
      for(int ii=0;ii<p2csz;) {
	int pid = p2c[ii].first ;

	int j = 1 ;
	while(((ii+j) < p2csz) && (pid == p2c[ii+j].first))
	  j++ ;
	if(j == 1) {
	  int tid = p2c[ii].second ;
	  double w = 1.0 ;

	  for(int kk=0;kk<vs;++kk)
	    qic_loc[tid][kk] += qvecin[pid][kk] ;
	  weights[tid] += w ;
	} else {
	  for(int k=ii;k<ii+j;++k) {
	    int tid = p2c[k].second ;
	    for(int kk=0;kk<vs;++kk)
	      qic_loc[tid][kk] += qvecin[pid][kk] ;
	    weights[tid] += 1. ;
	  }
	}
	
	ii += j ;
      }

      // Now return interpolated versions to new file numbering
      vector<T> senddata(tsz*vs) ;
      vector<double> senddataw(tsz,0) ;
      vector<T> recvdata(rsz*vs) ;
      vector<double> recvdataw(rsz,1) ;
      // Fill send data buffer
      for(int i=0;i<tsz;++i) {
	for(int j=0;j<vs;++j)
	  senddata[i*vs+j] = qic_loc[i][j] ;
	senddataw[i] = weights[i] ;
      }
      for(int i=0;i<p;++i) {
	if(rlistsizes[i] > 0) {
	  MPI_Request tmp ;
	  MPI_Irecv(&recvdata[rlistoffsets[i]*vs],rlistsizes[i]*vs*sizeof(T),
		    MPI_BYTE,i,99,MPI_COMM_WORLD,&tmp) ;
	  req_queue.push_back(tmp) ;
	  MPI_Irecv(&recvdataw[rlistoffsets[i]],rlistsizes[i],
		    MPI_DOUBLE,i,98,MPI_COMM_WORLD,&tmp) ;
	  req_queue.push_back(tmp) ;
	}
      }
      for(int i=0;i<p;++i) {
	if(tlistsizes[i] > 0) {
	  MPI_Send(&senddata[tlistoffsets[i]*vs],tlistsizes[i]*vs*sizeof(T),
		   MPI_BYTE,i,99,MPI_COMM_WORLD) ;
	  MPI_Send(&senddataw[tlistoffsets[i]],tlistsizes[i],
		   MPI_DOUBLE,i,98,MPI_COMM_WORLD) ;
	}
      }
      if(req_queue.size() > 0) {
	vector<MPI_Status> stat_queue(req_queue.size()) ;
	MPI_Waitall(req_queue.size(),&req_queue[0],&stat_queue[0]) ;
	req_queue.clear() ;
      }

      storeVec<double> q_ic_fn ;
      store<double> wt_fn ;

      entitySet dom2=dom & qvecout.domain() ;
      WARN(dom!=dom2) ;
      wt_fn.allocate(dom) ;
      // sum contributions from processors
      FORALL(dom,ii) {
	for(int i=0;i<vs;++i)
	  setZero(qvecout[ii][i]) ;
	wt_fn[ii] = 0 ;
      } ENDFORALL ;
      for(int i=0;i<rsz;++i) {
	int id =g2l[recvlists[i]] ;

	for(int k=0;k<vs;++k)
	  qvecout[id][k] += recvdata[i*vs+k] ;
	wt_fn[id] += recvdataw[i] ;
      }

      // re-normalize
      FORALL(dom2,ii) {
	double rw = 1./(wt_fn[ii]) ;
	for(int i=0;i<vs;++i)
	  qvecout[ii][i] *= rw;
      } ENDFORALL ;

    }
  }

}
#endif
