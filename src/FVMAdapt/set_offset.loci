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
//////////////////////////////////////////////////////////////////////////////////////////////////////
//                       set_offset.cc
// 
//
//////////////////////////////////////////////////////////////////////////////////////////////////////


#include <Loci.h>
#include <stdio.h>
#include "mpi.h"
#include "defines.h"
#include "dataxferDB.h"

#include <iostream>
#include <fstream>
using Loci::storeRepP;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
typedef Loci::vector3d<double> vect3d;
namespace Loci{
  std::vector<int> all_collect_sizes(int size);
  // Convert container from local numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm);

  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm);
  
}





class get_node_offset : public pointwise_rule {
  const_param<int> num_original_nodes;
  const_store<int> num_inner_nodes;
  store<int> node_offset;
  
public:
  get_node_offset(){
    name_store("num_inner_nodes", num_inner_nodes);
    name_store("num_original_nodes", num_original_nodes);
    name_store("node_offset", node_offset);
    input("num_inner_nodes");
    input("num_original_nodes");
    output("node_offset");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {

    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
   
  
    Loci::constraint edges, geom_cells, faces;
    Loci::storeRepP e2n = Loci::exec_current_fact_db->get_variable("edge2node");
    *edges = e2n->domain();
    faces = Loci::exec_current_fact_db->get_variable("faces");
    geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");


 

    if(num_procs ==  1 || dist == 0){
      
      // int offset = *num_original_nodes + 1; //in cobalt, index of node start with 1
      int offset = *num_original_nodes ; //in vog, index of node start with 0
      FORALL(*edges, ei){
        node_offset[ei] = offset;
        offset += num_inner_nodes[ei];
      }ENDFORALL;
      //write out cell_nodes first, then write face_nodes
      FORALL(*geom_cells, ei){
        node_offset[ei] = offset;
        offset += num_inner_nodes[ei];
      }ENDFORALL;
      
      FORALL(*faces, ei){
        node_offset[ei] = offset;
        offset += num_inner_nodes[ei];
      }ENDFORALL;
      return;
    }
    
    Loci::constraint my_entities ; 
    my_entities = dist->my_entities ;

    entitySet local_nodes, local_edges, local_faces, local_geom_cells;
    //  local_nodes = *my_entities & nodes ;
    //don't know if it's necessray
    local_edges = (*my_entities) & (*edges) ;
    local_faces = (*my_entities) & (*faces);
    local_geom_cells = (*my_entities)&(*geom_cells);

    //each process computes its node  offset

    int noffset = *num_original_nodes ;

    int offset = 0;
    int num_local_nodes = 0;
    std::vector<int> local_nodes_sizes;
    storeRepP localVar=node_offset.Rep();
    //put into block to make the memory for store<int> freeed after it's used 
    {
      store<int> edge_num_inner_nodes;
      edge_num_inner_nodes = Loci::Local2FileOrder(num_inner_nodes.Rep(),local_edges,offset,dist,MPI_COMM_WORLD) ;
      FORALL(edge_num_inner_nodes.domain(), ee){
        num_local_nodes += edge_num_inner_nodes[ee];
      }ENDFORALL;
      local_nodes_sizes = Loci::all_collect_sizes(num_local_nodes);
      
    
      for(int i = 0; i < my_id; i++){
        noffset += local_nodes_sizes[i];
      }
      
      //compute the store values
      store<int> edge_file_offset;
      edge_file_offset.allocate(edge_num_inner_nodes.domain());
      FORALL(edge_num_inner_nodes.domain(), ei){
        edge_file_offset[ei] = noffset;
        noffset += edge_num_inner_nodes[ei];
      }ENDFORALL;
      Loci::File2LocalOrder(localVar, local_edges,
                            edge_file_offset.Rep(), offset,
                            dist,
                            MPI_COMM_WORLD);
      
    }
    //finish with edge nodes
    int total_num_edge_nodes = 0;
    for(int i = 0; i < num_procs; i++){
      total_num_edge_nodes += local_nodes_sizes[i];
    } 

    noffset = *num_original_nodes + total_num_edge_nodes;

    
   
    
    //repeat for geom_cells

 
    
    
    offset= 0;
    // Create container vardist that is ordered across processors in the
    // file numbering, the domain of this container shifted by offset
    // is the actual file numbering.
    {
      store<int> cell_num_inner_nodes;
      cell_num_inner_nodes = Loci::Local2FileOrder(num_inner_nodes.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;
   
      num_local_nodes = 0;
      FORALL(cell_num_inner_nodes.domain(), ei){
        num_local_nodes += cell_num_inner_nodes[ei];
      }ENDFORALL;
      local_nodes_sizes = Loci::all_collect_sizes(num_local_nodes);
      for(int i = 0; i < my_id; i++){
        noffset += local_nodes_sizes[i];
      }
      
      //compute the store values
      store<int> cell_file_offset;
      cell_file_offset.allocate(cell_num_inner_nodes.domain());
      FORALL(cell_num_inner_nodes.domain(), ei){
        cell_file_offset[ei] = noffset;
        noffset += cell_num_inner_nodes[ei];
      }ENDFORALL;
      //File2Local use the offset value set by Local2File    
      Loci::File2LocalOrder(localVar, local_geom_cells,
                            cell_file_offset.Rep(), offset,
                            dist,
                            MPI_COMM_WORLD);
      //finish with cell nodes
    }
    //update noffset
    int total_num_cell_nodes = 0;
    for(int i = 0; i < num_procs; i++){
      total_num_cell_nodes += local_nodes_sizes[i];
    } 
    noffset = *num_original_nodes + total_num_edge_nodes + total_num_cell_nodes;

    
    
     
    //repeat for faces
    offset = 0; 
    // Create container vardist that is ordered across processors in the
    // file numbering, the domain of this container shifted by offset
    // is the actual file numbering.
    {
      store<int> face_num_inner_nodes;
      face_num_inner_nodes = Loci::Local2FileOrder(num_inner_nodes.Rep(),local_faces, offset,dist,MPI_COMM_WORLD) ;
          
      num_local_nodes = 0;
      FORALL(face_num_inner_nodes.domain(), ei){
        num_local_nodes += face_num_inner_nodes[ei];
      }ENDFORALL;
      
      local_nodes_sizes = Loci::all_collect_sizes(num_local_nodes);
      for(int i = 0; i < my_id; i++){
        noffset += local_nodes_sizes[i];
      }
      
      //compute the store values
      store<int> face_file_offset;
      face_file_offset.allocate(face_num_inner_nodes.domain());
      FORALL(face_num_inner_nodes.domain(), ei){
        face_file_offset[ei] = noffset;
        noffset += face_num_inner_nodes[ei];
      }ENDFORALL;
      
      Loci::File2LocalOrder(localVar, local_faces,
                            face_file_offset.Rep(), offset,
                            dist,
                            MPI_COMM_WORLD);
      
    } 
  }
} ;
register_rule<get_node_offset> register_get_node_offset;


class get_cell_offset : public pointwise_rule {
  const_store<int> num_fine_cells;
  store<int> cell_offset;
    
public:
  get_cell_offset(){
    name_store("num_fine_cells", num_fine_cells);
    name_store("cell_offset", cell_offset);
    input("num_fine_cells");
    output("cell_offset");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
    
    Loci::constraint  geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
       
    if(num_procs ==  1){
      int offset = 0;
      FORALL(*geom_cells, ei){
        cell_offset[ei] = offset;
        offset += num_fine_cells[ei];
      }ENDFORALL;
      
      return;
    }
    
    Loci::constraint my_entities ; 
    my_entities = dist->my_entities ;
   
    //don't know if it's necessray
    entitySet local_geom_cells = (*my_entities)&(*geom_cells);

    //each process computes its cell  offset

    int coffset = 0 ;

    int offset = 0;
    int num_local_cells = 0;
    std::vector<int> local_cells_sizes;
    storeRepP localVar=cell_offset.Rep();
    //put into block to make the memory for store<int> freeed after it's used 
    {
      store<int> file_num_fine_cells;
      file_num_fine_cells = Loci::Local2FileOrder(num_fine_cells.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;
      FORALL(file_num_fine_cells.domain(), ee){
        num_local_cells += file_num_fine_cells[ee];
      }ENDFORALL;
      local_cells_sizes = Loci::all_collect_sizes(num_local_cells);
      
    
      for(int i = 0; i < my_id; i++){
        coffset += local_cells_sizes[i];
      }
      
      //compute the store values
      store<int> cell_file_offset;
      cell_file_offset.allocate(file_num_fine_cells.domain());
      FORALL(file_num_fine_cells.domain(), ei){
        cell_file_offset[ei] = coffset;
        coffset += file_num_fine_cells[ei];
      }ENDFORALL;
      Loci::File2LocalOrder(localVar, local_geom_cells,
                            cell_file_offset.Rep(), offset,
                            dist,
                            MPI_COMM_WORLD);
      
    }
    
    
   
  }
    
} ;
register_rule<get_cell_offset> register_get_cell_offset;



class get_parent_cell_offset : public pointwise_rule {
  const_store<int> num_fine_cells;
  store<int> cell_offset;
    
public:
  get_parent_cell_offset(){
    name_store("parent_num_fine_cells", num_fine_cells);
    name_store("parent_cell_offset", cell_offset);
    input("parent_num_fine_cells");
    output("parent_cell_offset");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq) {
    int num_procs = Loci::MPI_processes;
    int my_id = Loci::MPI_rank;
    fact_db::distribute_infoP dist = Loci::exec_current_fact_db->get_distribute_info() ;
    
    Loci::constraint  geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
       
    if(num_procs ==  1){
      int offset = 0;
      FORALL(*geom_cells, ei){
        cell_offset[ei] = offset;
        offset += num_fine_cells[ei];
      }ENDFORALL;
      
      return;
    }
    
    Loci::constraint my_entities ; 
    my_entities = dist->my_entities ;
   
    //don't know if it's necessray
    entitySet local_geom_cells = (*my_entities)&(*geom_cells);

    //each process computes its cell  offset

    int coffset = 0 ;

    int offset = 0;
    int num_local_cells = 0;
    std::vector<int> local_cells_sizes;
    storeRepP localVar=cell_offset.Rep();
    //put into block to make the memory for store<int> freeed after it's used 
    {
      store<int> file_num_fine_cells;
      file_num_fine_cells = Loci::Local2FileOrder(num_fine_cells.Rep(),local_geom_cells,offset,dist,MPI_COMM_WORLD) ;
      FORALL(file_num_fine_cells.domain(), ee){
        num_local_cells += file_num_fine_cells[ee];
      }ENDFORALL;
      local_cells_sizes = Loci::all_collect_sizes(num_local_cells);
      
    
      for(int i = 0; i < my_id; i++){
        coffset += local_cells_sizes[i];
      }
      
      //compute the store values
      store<int> cell_file_offset;
      cell_file_offset.allocate(file_num_fine_cells.domain());
      FORALL(file_num_fine_cells.domain(), ei){
        cell_file_offset[ei] = coffset;
        coffset += file_num_fine_cells[ei];
      }ENDFORALL;
      Loci::File2LocalOrder(localVar, local_geom_cells,
                            cell_file_offset.Rep(), offset,
                            dist,
                            MPI_COMM_WORLD);
      
   
    }
  }
    
} ;
register_rule<get_parent_cell_offset> register_get_parent_cell_offset;


class init_num_original_faces : public unit_rule{
  param<int> num_original_faces;
public:
  init_num_original_faces(){
    name_store("num_original_faces", num_original_faces);
    output("num_original_faces");
    constraint("UNIVERSE");
    
  }
  //parameter, no loop, 
  virtual void compute(const sequence &seq){
    
    
    *num_original_faces = 0;
  }
}; 
register_rule<init_num_original_faces> register_init_num_original_faces;

class apply_num_original_faces : public apply_rule<param<int>, Loci::Summation<int> >{
  param<int> num_original_faces;
  
public:
  apply_num_original_faces(){
    name_store("num_original_faces", num_original_faces);
    input("num_original_faces");
    output("num_original_faces");
    constraint("faces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
    *num_original_faces +=1;
  }
}; 
register_rule<apply_num_original_faces> register_apply_num_original_faces;

class get_cell_parent : public pointwise_rule {
  const_store<int> num_fine_cells;
  const_store<int> cell_offset;
  const_store<int>  cell_l2f;
  const_param<int> num_original_nodes;
  const_param<int> num_original_faces;
  
  store<std::vector<pair<int32, int32> > > cell2parent;
  
public:
  get_cell_parent(){
    name_store("balanced_num_fine_cells", num_fine_cells);
    name_store("balanced_cell_offset", cell_offset);
    name_store("fileNumber(geom_cells)", cell_l2f);
    name_store("cell2parent", cell2parent);
    name_store("num_original_nodes", num_original_nodes);
    name_store("num_original_faces", num_original_faces);
   
    input("balanced_num_fine_cells");
    input("balanced_cell_offset");
    input("fileNumber(geom_cells)");
    input("num_original_nodes");
    input("num_original_faces");
    
   
    output("cell2parent");
    constraint("geom_cells");
  }
  virtual void compute(const sequence &seq) {
    if(seq.size()!=0){
       
      do_loop(seq, this);
    }
  }
  void calculate(Entity cc){
    std::vector<pair<int32, int32> > c2p(num_fine_cells[cc]);
    for(int i = 0; i < num_fine_cells[cc]; i++){
      
      int child_index = cell_offset[cc]+i+1;//local cellindex start with 1
      int parent_index = cell_l2f[cc]-*num_original_nodes-*num_original_faces+1;
      
      
      c2p[i] = make_pair(child_index, parent_index);
     
    }
    c2p.swap(cell2parent[cc]);
  }   
} ;
register_rule<get_cell_parent> register_get_cell_parent;

class get_cell2parent : public pointwise_rule {
  const_store<int> cell_offset;
  const_store<int> parent_cell_offset;
  const_store<vector<pair<int32, int32> > > indexMap;
  store<vector<pair<int32, int32> > > cell2parent;
public:
  get_cell2parent(){
    name_store("balanced_cell_offset", cell_offset);
    name_store("parent_cell_offset", parent_cell_offset);
    name_store("indexMap", indexMap);
    name_store("priority::restart::cell2parent", cell2parent);
    input("balanced_cell_offset");
    input("parent_cell_offset");
    input("indexMap");
    output("priority::restart::cell2parent");
    constraint("geom_cells");
  }
  virtual void compute(const sequence &seq) {
    
    if(seq.size()!=0){
      do_loop(seq, this);
    }
  }
  void calculate(Entity cc){
    vector<pair<int32, int32> > c2p(indexMap[cc].size());

    for(unsigned int i = 0; i < indexMap[cc].size(); i++){
      c2p[i]=  make_pair(indexMap[cc][i].first+cell_offset[cc],
                         indexMap[cc][i].second+parent_cell_offset[cc]);
    }
    c2p.swap(cell2parent[cc]);
  }
} ;
register_rule<get_cell2parent> register_get_cell2parent;


class write_cell2parent : public pointwise_rule{
  const_param<std::string> c2pfile_par;
  const_store<std::vector<pair<int32, int32> > > cell2parent;
  store<bool> c2p_output;
public:
  write_cell2parent(){
    name_store("cell2parent_file_par", c2pfile_par);
    name_store("cell2parent", cell2parent);
    name_store("cell2parent_output", c2p_output);
    input("cell2parent_file_par");
    input("cell2parent");
    output("cell2parent_output");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    std::ofstream outFile;
    int nprocs = Loci::MPI_processes;
    
    //process 0 open c2p File
    if(Loci::MPI_rank == 0){
      outFile.open((*c2pfile_par).c_str());
      if(!outFile){
        cerr <<"can not open " << *c2pfile_par << " for output" << endl;
        Loci::Abort();
      }
    }
    
    entitySet dom = entitySet(seq);
    
    //serial version
    if(nprocs == 1){
      FORALL(dom, cc){
        for(unsigned int i = 0; i < cell2parent[cc].size(); i++){
          outFile << cell2parent[cc][i].first << ' '<<cell2parent[cc][i].second<<endl;
        }
      }ENDFORALL;
      //close c2p file
      outFile.close();
      return;
    }
    //parallel version
    
   
    //first create store in the order of file numbering      
   
   
   
    int my_size = 0;
    FORALL(dom, cc){
      my_size += 2*cell2parent[cc].size();
    }ENDFORALL;
    
    int* size_buf = new int[nprocs];
    
    //compute buf_size on each process
    unsigned int  buf_size = 0;
    MPI_Allreduce(&my_size, &buf_size, 1, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);
    //process 0 find out size of buffer for each process
     
    MPI_Gather(&my_size, 1, MPI_INT, size_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
    int32 *buf = new int32[buf_size]; 
    
    if(Loci::MPI_rank == 0){
      
      //process 0 write its local cell2parent 
      FORALL(dom, cc){
        for(unsigned int i = 0; i < cell2parent[cc].size(); i++){
          outFile << cell2parent[cc][i].first << ' '<<cell2parent[cc][i].second<<endl;
        }
      }ENDFORALL;
      
      //process 0 recv the buf and write it out
      for(int i = 1; i < nprocs; i++){
        MPI_Status status;
        MPI_Recv(buf, buf_size, MPI_INT, i, 20, MPI_COMM_WORLD, &status);
        int recv_size = size_buf[i]/2; 
        for(int j = 0; j < recv_size; j++)  
          {
            outFile << buf[2*j] << ' '<<buf[2*j+1]<<endl;
          }
       
      }
      
    }else{ //other processes send buf to process 0
   
      int ptr = 0;
      FORALL(dom, cc){
        for(unsigned int i = 0; i < cell2parent[cc].size(); i++){
          buf[ptr++] = cell2parent[cc][i].first;
          buf[ptr++] = cell2parent[cc][i].second;
        }
      }ENDFORALL;
      int send_size = ptr;
      MPI_Send(buf, send_size, MPI_INT, 0, 20, MPI_COMM_WORLD);
    }
    
   
               
    delete [] buf;
    delete [] size_buf;
    //process 0 close file
    if(Loci::MPI_rank == 0){
      outFile.close();
      //  cout << "Finish reading  posTag " << endl;
    }
    
  } 

};
register_rule<write_cell2parent> register_write_cell2parent;
  

class save_cell2parent : public pointwise_rule{
  const_store<std::vector<pair<int32, int32> > > cell2parent;
  store<bool> c2p_output;
public:
  save_cell2parent(){
    name_store("cell2parent", cell2parent);
    name_store("cell2parent_DB", c2p_output);
    input("cell2parent");
    output("cell2parent_DB");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    entitySet dom = entitySet(seq);
    int cnt = 0 ;
    FORALL(dom, cc){
      cnt += cell2parent[cc].size() ;
    }ENDFORALL;
    store<pair<int,int> > c2pset ;
    entitySet domset = interval(0,cnt-1) ;
    c2pset.allocate(domset) ;
    
    cnt = 0 ;
    FORALL(dom, cc){
      for(unsigned int i = 0; i < cell2parent[cc].size(); i++){
	c2pset[cnt] = cell2parent[cc][i] ;
	cnt++ ;
      }
    } ENDFORALL ;
    
    Loci::DataXFER_DB.insertItem("c2p",c2pset.Rep()) ;
  } 

};
register_rule<save_cell2parent> register_save_cell2parent;
  





