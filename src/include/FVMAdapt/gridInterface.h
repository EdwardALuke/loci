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

}
#endif
