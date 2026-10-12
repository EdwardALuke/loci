//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
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
#ifndef GRID_INTERFACE_H
#define GRID_INTERFACE_H
#include <store_rep.h>
#include <memory>

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
    // global variables that pass data from refinemesh to a solver
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

  // Holds Data Structures for processning AMR interpolation
  class AMRrefinementMapping {
  public:
    entitySet geom_cells_local ;
    entitySet geom_cells_global ;
    Map l2g ;
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
    multiStore<vector3d<double> > stencilWeights ;
    // deltas used to compute stencil
    multiStore<vector3d<double> > grad_dvs ;
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
    vector<double>  directWeights ;
    void setupRefinementMapping(const store<pair<int,int> > &c2pg,
                                const store<double> &volw,
                                multiStore<int> &gradCells,
                                multiStore<vector3d<double> > &deltas,
                                const_store<vector3d<double> > &cell_center,
                                const_store<double> &vol,
                                fact_db &facts) ;
    template<class T>
    void interpolateData(storeVec<T> &parent_data,
                         storeVec<T> &child_data,
                         double limstop = 1.0) const {
      //########################################################################
      //
      // First interpolate to cells that resulted from splitting parent cells
      //
      // Step 1, compute gradients using parent data for parents that form
      // split cells
      // Step 2, compute limiter for these gradients
      // Step 3, compute limited second order interpolation to child cells
      // Step 4, distribute second order reconstruction from parent to child
      //         cells
      // Step 5, collect direct mapped cells and use volume averaging to
      //         interpolate from many parents to child for coarsened cells
      // Note, if the variables given are conservative, the extrapolation is
      // conservative
      //
      const int vs = parent_data.vecSize() ;
      // gather parent data needed to compute parent gradients
      storeVec<T> grad_data ;
      gradientComm.gatherData(grad_data,parent_data) ;

      storeVec<T> refineCell_data ;
      refineCell_data.setVecSize(vs) ;
      int refcells = refinedCells.size() ;
#ifdef DEBUG
      Loci::debugout  << "#refcells interpolations = " << refcells << endl ;
#endif
      entitySet refineCellDomain ;
      if(refcells>0)
        refineCellDomain = interval(0,refcells-1) ;
      refineCell_data.allocate(refineCellDomain) ;
      int vsz = parent_data.vecSize() ;
      FORALL(parent.domain(),ii) {
        // get stencil size
        const int ssz = gradCellStencil.vec_size(ii) ;
        int parent_id = parent[ii] ;
        for(int k=0;k<vsz;++k) { // loop over vector
          // Compute gradient
          vector3d<double> gradk(0,0,0) ;
          T Xcc = parent_data[parent_id][k] ;
          T max_val = Xcc ;
          T min_val = max_val ;
          for(int j=0;j<ssz;++j) {
            const int cc = gradCellStencil[ii][j] ;
            T val = grad_data[cc][k] ;
            max_val = max(max_val,val) ;
            min_val = min(min_val,val) ;
            gradk += stencilWeights[ii][j]*(val-Xcc) ;
          }
          T limi = limstop ;
          // Apply limiter to source points
          for(int j=0;j<ssz;++j) {
            T qdif = dot(gradk,grad_dvs[ii][j]) ;
            if(qdif > 0)
              limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
            if(qdif < 0)
              limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
          }
          // Apply limiter to target points
          for(int j=0;j<child_dvs.vec_size(ii);++j) {
            T qdif =  dot(vector3d<double>(gradk),
                          child_dvs[ii][j]) ;
            if(qdif > 0)
              limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
            if(qdif < 0)
              limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
          }

          // reduce gradient by limiter factor
          gradk *= limi ;
          // compute child cell values
          for(int j=0;j<child_dvs.vec_size(ii);++j) {
            T val = Xcc + dot(vector3d<double>(gradk),
                              child_dvs[ii][j]) ;
            refineCell_data[parent2child_l[ii][j]][k] = val ;
          }
        }
      } ENDFORALL ;

      storeVec<T> child_tmp ;
      child_tmp.setVecSize(child_data.vecSize()) ;
      child_tmp.allocate(geom_cells_global) ;

      refineCellComm.scatterData(refineCell_data,child_tmp) ;


      // direct map cells with averaging for coarse cells
      storeVec<T> mapCell_data ;
      directMapComm.gatherData(mapCell_data,parent_data) ;
      for(size_t i=0;i<c2lp.size();++i) {
        int child = c2lp[i].first ;
        double w = directWeights[i] ;
        for(int k=0;k<vs;++k)
          child_tmp[child][k] = w*mapCell_data[c2lp[i].second][k] ;
        size_t j = i+1 ;
        for(;j<c2lp.size()&&child == c2lp[j].first;++j) {
          double w = directWeights[j] ;
          for(int k=0;k<vs;++k)
            child_tmp[child][k] +=  w*mapCell_data[c2lp[j].second][k] ;
        }
        int cnt = j-i ;
        i+= cnt-1 ;
      }
      FORALL(geom_cells_local,ii) {
        child_data[ii] = child_tmp[l2g[ii]] ;
      } ENDFORALL ;
    }

    template<class T>
    void interpolateData(store<T> &parent_data,
                         store<T> &child_data,
                         double limstop = 1.0) const {
      //########################################################################
      //
      // First interpolate to cells that resulted from splitting parent cells
      //
      // Step 1, compute gradients using parent data for parents that form
      // split cells
      // Step 2, compute limiter for these gradients
      // Step 3, compute limited second order interpolation to child cells
      // Step 4, distribute second order reconstruction from parent to child
      //         cells
      // Step 5, collect direct mapped cells and use volume averaging to
      //         interpolate from many parents to child for coarsened cells
      // Note, if the variables given are conservative, the extrapolation is
      // conservative
      //
      // gather parent data needed to compute parent gradients
      store<T> grad_data ;
      gradientComm.gatherData(grad_data,parent_data) ;

      store<T> refineCell_data ;
      int refcells = refinedCells.size() ;
#ifdef DEBUG
      Loci::debugout  << "#refcells interpolations = " << refcells << endl ;
#endif
      entitySet refineCellDomain ;
      if(refcells>0)
        refineCellDomain = interval(0,refcells-1) ;
      refineCell_data.allocate(refineCellDomain) ;
      FORALL(parent.domain(),ii) {
        // get stencil size
        const int ssz = gradCellStencil.vec_size(ii) ;
        // Compute gradient
        vector3d<double> gradk(0,0,0) ;
        T Xcc = parent_data[parent[ii]] ;
        T max_val = Xcc ;
        T min_val = max_val ;
        for(int j=0;j<ssz;++j) {
          const int cc = gradCellStencil[ii][j] ;
          T val = grad_data[cc] ;
          max_val = max(max_val,val) ;
          min_val = min(min_val,val) ;
          gradk += stencilWeights[ii][j]*(val-Xcc) ;
        }

        T limi = limstop ;
        // Apply limiter to source points
        for(int j=0;j<ssz;++j) {
          T qdif = dot(gradk,grad_dvs[ii][j]) ;
          if(qdif > 0)
            limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
          if(qdif < 0)
            limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
        }
        // Apply limiter to target points
        for(int j=0;j<child_dvs.vec_size(ii);++j) {
          T qdif =  dot(vector3d<double>(gradk),
                        child_dvs[ii][j]) ;
          if(qdif > 0)
            limi = min(limi,(max_val-Xcc)/(qdif+1e-100)) ;
          if(qdif < 0)
            limi = min(limi,(min_val-Xcc)/(qdif-1e-100)) ;
        }

        // reduce gradient by limiter factor
        gradk *= limi ;
        // compute child cell values
        for(int j=0;j<child_dvs.vec_size(ii);++j) {
          T val = Xcc + dot(vector3d<double>(gradk),
                            child_dvs[ii][j]) ;
          refineCell_data[parent2child_l[ii][j]] = val ;
        }
      } ENDFORALL ;

      store<T> child_tmp ;
      child_tmp.allocate(geom_cells_global) ;
      refineCellComm.scatterData(refineCell_data,child_tmp) ;


      // direct map cells with averaging for coarse cells
      store<T> mapCell_data ;
      directMapComm.gatherData(mapCell_data,parent_data) ;
      for(size_t i=0;i<c2lp.size();++i) {
        int child = c2lp[i].first ;
        double w = directWeights[i] ;
        child_tmp[child] = w*mapCell_data[c2lp[i].second] ;
        size_t j = i+1 ;
        for(;j<c2lp.size()&&child == c2lp[j].first;++j) {
          double w = directWeights[j] ;
          child_tmp[child] +=  w*mapCell_data[c2lp[j].second] ;
        }
        int cnt = j-i ;
        i+= cnt-1 ;
      }
      FORALL(geom_cells_local,ii) {
        child_data[ii] = child_tmp[l2g[ii]] ;
      } ENDFORALL ;
    }

    template<class T>
    void interpolateData(store<vector3d<T> > &parent_data,
                         store<vector3d<T> > &child_data,
                         double limstop = 1.0) const {
      //########################################################################
      //
      // First interpolate to cells that resulted from splitting parent cells
      //
      // Step 1, compute gradients using parent data for parents that form
      // split cells
      // Step 2, compute limiter for these gradients
      // Step 3, compute limited second order interpolation to child cells
      // Step 4, distribute second order reconstruction from parent to child
      //         cells
      // Step 5, collect direct mapped cells and use volume averaging to
      //         interpolate from many parents to child for coarsened cells
      // Note, if the variables given are conservative, the extrapolation is
      // conservative
      //
      // gather parent data needed to compute parent gradients
      store<vector3d<T> > grad_data ;
      gradientComm.gatherData(grad_data,parent_data) ;

      store<vector3d<T> > refineCell_data ;
      int refcells = refinedCells.size() ;
#ifdef DEBUG
      Loci::debugout  << "#refcells interpolations = " << refcells << endl ;
#endif
      entitySet refineCellDomain ;
      if(refcells>0)
        refineCellDomain = interval(0,refcells-1) ;
      refineCell_data.allocate(refineCellDomain) ;
      FORALL(parent.domain(),ii) {
        // get stencil size
        const int ssz = gradCellStencil.vec_size(ii) ;
        // Compute gradient
        vector3d<double> gradkx(0,0,0) ;
        vector3d<double> gradky(0,0,0) ;
        vector3d<double> gradkz(0,0,0) ;
        vector3d<T> Xcc = parent_data[parent[ii]] ;
        vector3d<T> max_val = Xcc ;
        vector3d<T> min_val = max_val ;
        for(int j=0;j<ssz;++j) {
          const int cc = gradCellStencil[ii][j] ;
          vector3d<T> val = grad_data[cc] ;
          max_val.x = max(max_val.x,val.x) ;
          max_val.y = max(max_val.y,val.y) ;
          max_val.z = max(max_val.z,val.z) ;
          min_val.x = min(min_val.x,val.x) ;
          min_val.y = min(min_val.y,val.y) ;
          min_val.z = min(min_val.z,val.z) ;
          gradkx += stencilWeights[ii][j]*(val.x-Xcc.x) ;
          gradky += stencilWeights[ii][j]*(val.y-Xcc.y) ;
          gradkz += stencilWeights[ii][j]*(val.z-Xcc.z) ;
        }

        T limi = limstop ;
        // Apply limiter to source points
        for(int j=0;j<ssz;++j) {
          T qdifx = dot(gradkx,grad_dvs[ii][j]) ;
          if(qdifx > 0)
            limi = min(limi,(max_val.x-Xcc.x)/(qdifx+1e-100)) ;
          if(qdifx < 0)
            limi = min(limi,(min_val.x-Xcc.x)/(qdifx-1e-100)) ;
          T qdify = dot(gradky,grad_dvs[ii][j]) ;
          if(qdify > 0)
            limi = min(limi,(max_val.y-Xcc.y)/(qdify+1e-100)) ;
          if(qdify < 0)
            limi = min(limi,(min_val.y-Xcc.y)/(qdify-1e-100)) ;
          T qdifz = dot(gradkz,grad_dvs[ii][j]) ;
          if(qdifz > 0)
            limi = min(limi,(max_val.z-Xcc.z)/(qdifz+1e-100)) ;
          if(qdifz < 0)
            limi = min(limi,(min_val.z-Xcc.z)/(qdifz-1e-100)) ;
        }
        // Apply limiter to target points
        for(int j=0;j<child_dvs.vec_size(ii);++j) {
          T qdifx =  dot(vector3d<double>(gradkx),
                        child_dvs[ii][j]) ;
          if(qdifx > 0)
            limi = min(limi,(max_val.x-Xcc.x)/(qdifx+1e-100)) ;
          if(qdifx < 0)
            limi = min(limi,(min_val.x-Xcc.x)/(qdifx-1e-100)) ;
          T qdify =  dot(vector3d<double>(gradky),
                        child_dvs[ii][j]) ;
          if(qdify > 0)
            limi = min(limi,(max_val.y-Xcc.y)/(qdify+1e-100)) ;
          if(qdify < 0)
            limi = min(limi,(min_val.y-Xcc.y)/(qdify-1e-100)) ;
          T qdifz =  dot(vector3d<double>(gradkz),
                        child_dvs[ii][j]) ;
          if(qdifz > 0)
            limi = min(limi,(max_val.z-Xcc.z)/(qdifz+1e-100)) ;
          if(qdifz < 0)
            limi = min(limi,(min_val.z-Xcc.z)/(qdifz-1e-100)) ;
        }

        // reduce gradient by limiter factor
        gradkx *= limi ;
        gradky *= limi ;
        gradkz *= limi ;
        // compute child cell values
        for(int j=0;j<child_dvs.vec_size(ii);++j) {
          T valx = Xcc.x + dot(vector3d<double>(gradkx),
                               child_dvs[ii][j]) ;
          T valy = Xcc.y + dot(vector3d<double>(gradky),
                               child_dvs[ii][j]) ;
          T valz = Xcc.z + dot(vector3d<double>(gradkz),
                               child_dvs[ii][j]) ;
          refineCell_data[parent2child_l[ii][j]] = vector3d<T>(valx,valy,valz) ;
        }
      } ENDFORALL ;

      store<vector3d<T> > child_tmp ;
      child_tmp.allocate(geom_cells_global) ;
      refineCellComm.scatterData(refineCell_data,child_tmp) ;


      // direct map cells with averaging for coarse cells
      store<vector3d<T> > mapCell_data ;
      directMapComm.gatherData(mapCell_data,parent_data) ;
      for(size_t i=0;i<c2lp.size();++i) {
        int child = c2lp[i].first ;
        double w = directWeights[i] ;
        child_tmp[child] = w*mapCell_data[c2lp[i].second] ;
        size_t j = i+1 ;
        for(;j<c2lp.size()&&child == c2lp[j].first;++j) {
          double w = directWeights[j] ;
          child_tmp[child] +=  w*mapCell_data[c2lp[j].second] ;
        }
        int cnt = j-i ;
        i+= cnt-1 ;
      }
      FORALL(geom_cells_local,ii) {
        child_data[ii] = child_tmp[l2g[ii]] ;
      } ENDFORALL ;
    }
  } ;

  class AMRinterpolation {
  private:
    std::shared_ptr<AMRrefinementMapping> interpolator ;
  public:
    void setupRefinementMapping(const store<pair<int,int> > &c2pg,
                                const store<double> &volw,
                                multiStore<int> &gradCells,
                                multiStore<vector3d<double> > &deltas,
                                const_store<vector3d<double> > &cell_center,
                                const_store<double> &vol,
                                fact_db &facts) {
      interpolator = std::make_shared<AMRrefinementMapping>() ;
      interpolator->setupRefinementMapping(c2pg,volw,gradCells,deltas,
                                           cell_center,vol,facts) ;
    }
    entitySet getOutputAllocate() const { return interpolator->geom_cells_global; }
    template<class T>
    void interpolateData(storeVec<T> &parent_data,
                         storeVec<T> &child_data,
                         double limstop=1.0) const {
      interpolator->interpolateData(parent_data,child_data,limstop) ;
    }
    template<class T>
    void interpolateData(store<T> &parent_data,
                         store<T> &child_data,
                         double limstop=1.0) const {
      interpolator->interpolateData(parent_data,child_data,limstop) ;
    }
    template<class T>
    void interpolateData(store<vector3d<T> > &parent_data,
                         store<vector3d<T> > &child_data,
                         double limstop=1.0) const {
      interpolator->interpolateData(parent_data,child_data,limstop) ;
    }
  } ;

}
#endif
