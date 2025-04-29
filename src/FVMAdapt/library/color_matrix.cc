//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
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
#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include "sciTypes.h"
#include "defines.h"
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using Loci::storeRepP;
using Loci::constraint;
using Loci::MPI_processes;
using Loci::MPI_rank;
using Loci::entitySet;
using Loci::UNIVERSE_MAX;
using Loci::UNIVERSE_MIN;

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
  std::vector<int> simplePartitionVec(int mn, int mx, int p);


  //writeVOGNode is different from the one from FVMGridWriter.cc
  namespace pio{ 
    void writeVOGNodeS(hid_t file_id,
                       Loci::storeRepP &pos,
                       const_store<Loci::FineNodes> &inner_nodes){ //serial io

      hid_t group_id = 0 ;
  
      if(MPI_processes == 1){
        //firsr write out numNodes
        long long num_original_nodes  = pos->domain().size();
        long long num_inner_nodes  = 0;
        FORALL(inner_nodes.domain(), cc){
          num_inner_nodes += inner_nodes[cc].size();
        }ENDFORALL;
      
        long long array_size = num_original_nodes + num_inner_nodes;
      
        group_id = H5Gcreate(file_id,"file_info",
                             H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

        cout << "num_nodes = " << array_size << endl ;
      
        hsize_t dims = 1 ;
        hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;

        hid_t att_id = H5Acreate(group_id,"numNodes", H5T_STD_I64BE,
                                 dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
        H5Awrite(att_id,H5T_NATIVE_LLONG,&array_size) ;
        H5Aclose(att_id) ;
        H5Gclose(group_id) ;

        if(array_size == 0)
          return ;


        //prepare to write positions
        group_id = H5Gcreate(file_id,"node_info",
                             H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
      
        //create dataspace and dataset
        int rank = 1 ;
        hsize_t dimension = array_size ;
        hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
        typedef data_schema_traits<vect3d> traits_type ;
        Loci::DatatypeP dp = traits_type::get_type() ;

        hsize_t start = 0 ;
        hsize_t stride = 1 ;
      
        hid_t dataset = H5Dcreate(group_id,"positions",dp->get_hdf5_type(),
                                  dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
      
        //first write out pos
        hsize_t count = num_original_nodes ;
        if(count != 0) {
          std::vector<vect3d> v_pos(count);
          //put pos_io in a vector
          store<vect3d>pos_io;
          pos_io = pos;
          int index = 0;
          FORALL(pos_io.domain(),cc){
            v_pos[index++] = pos_io[cc];
          }ENDFORALL;
        
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &v_pos[0]) ;
          H5Sclose(memspace) ;
        }
      
        //put inner_nodes in a vector
        Loci::constraint faces, geom_cells;
        faces = Loci::exec_current_fact_db->get_variable("faces");
        geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
        entitySet local_edges = Loci::exec_current_fact_db->get_variable("edge2node")->domain();
    
   
      
        //next, write out inner_nodes
        start += num_original_nodes  ;
      
        count = num_inner_nodes;
        if(num_inner_nodes != 0) {
          std::vector<vect3d> v_nodes(num_inner_nodes);
          //put inner_nodes in a vector

          long index = 0;
          FORALL(local_edges, cc){
            for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
              v_nodes[index++] = inner_nodes[cc][i];
            }
          }ENDFORALL;
          FORALL(*geom_cells, cc){
            for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
              v_nodes[index++] = inner_nodes[cc][i];
            }
          }ENDFORALL; 
          FORALL(*faces, cc){
            for(unsigned int i = 0; i < inner_nodes[cc].size(); i++){
              v_nodes[index++] = inner_nodes[cc][i];
            }
          }ENDFORALL;
                
        
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &v_nodes[0]) ;
          H5Sclose(memspace) ;
        }

        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        H5Gclose(group_id) ;
        return;
      }//end of if(Loci::MPI_Processes==1)
  
      //reorder store first, from local to io entities
      fact_db::distribute_infoP dist =  Loci::exec_current_fact_db->get_distribute_info() ;
      constraint  my_faces, my_geom_cells; 
      entitySet my_entities = dist->my_entities ;
      my_faces = Loci::exec_current_fact_db->get_variable("faces");
      my_geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
  
  
  
      entitySet local_edges =my_entities &
        (Loci::exec_current_fact_db->get_variable("edge2node"))->domain(); 
      entitySet local_faces = my_entities & *my_faces;
      entitySet local_cells =  my_entities & *my_geom_cells;
      entitySet local_nodes = my_entities & pos->domain();
  
  

      //before write out,create stores for pos, inner_nodes of cells and faces which
      //are ordered across processors in the file numbering, the domain of this container
      //shifted by offset is the actual file numbering. offset will be modified after function call  
 
      int offset = 0;
      store<vect3d> pos_io;
      pos_io = Loci::Local2FileOrder(pos, local_nodes, offset, dist, MPI_COMM_WORLD) ;
      entitySet file_nodes = pos_io.domain();
  
      offset = 0;
      store<Loci::FineNodes> edge_inner_nodes;
      edge_inner_nodes = Loci::Local2FileOrder(inner_nodes.Rep(),local_edges,offset,dist,MPI_COMM_WORLD) ;
      entitySet file_edges = edge_inner_nodes.domain();
  
      offset= 0;
      // Create container vardist that
      store<Loci::FineNodes> cell_inner_nodes;
      cell_inner_nodes = Loci::Local2FileOrder(inner_nodes.Rep(),local_cells,offset,dist,MPI_COMM_WORLD) ;
      entitySet file_cells = cell_inner_nodes.domain();
  
      offset= 0;
      // Create container vardist that
      store<Loci::FineNodes> face_inner_nodes;
      face_inner_nodes = Loci::Local2FileOrder(inner_nodes.Rep(),local_faces,offset,dist,MPI_COMM_WORLD) ;
      entitySet file_faces = face_inner_nodes.domain();
  
  
      //compute the size of pos
      int local_pos_size = file_nodes.size();
      std::vector<int> pos_sizes(Loci::MPI_processes) ;
      MPI_Gather(&local_pos_size,1,MPI_INT,
                 &pos_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;

  
  
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
  
  
      std::vector<int> inner_edge_nodes_sizes(Loci::MPI_processes);
      MPI_Gather(&num_local_edge_nodes,1,MPI_INT,
                 &inner_edge_nodes_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;
  
  
      std::vector<int> inner_cell_nodes_sizes(Loci::MPI_processes);
      MPI_Gather(&num_local_cell_nodes,1,MPI_INT,
                 &inner_cell_nodes_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;
  
      std::vector<int> inner_face_nodes_sizes(Loci::MPI_processes);
      MPI_Gather(&num_local_face_nodes,1,MPI_INT,
                 &inner_face_nodes_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD) ;
  




  
      if(Loci::MPI_rank == 0) {
        //compute array_size
        hsize_t array_size = 0 ;
        for(int i=0;i<MPI_processes;++i)
          array_size += (pos_sizes[i]+ inner_edge_nodes_sizes[i]+
                         inner_cell_nodes_sizes[i] + inner_face_nodes_sizes[i] );
        //first write out numNodes
    
        group_id = H5Gcreate(file_id,"file_info",
                             H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    
        cout << "num_nodes = " << array_size << endl ;
    
        hsize_t dims = 1 ;
        hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
        hid_t att_id = H5Acreate(group_id,"numNodes", H5T_STD_I64BE,
                                 dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
        H5Awrite(att_id,H5T_NATIVE_LLONG,&array_size) ;
        H5Aclose(att_id) ;
        H5Gclose(group_id) ;
    
        if(array_size == 0)
          return ;

        group_id = H5Gcreate(file_id,"node_info",
                             H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
        //create dataspace and dataset
        int rank = 1 ;
        hsize_t dimension = array_size ;
        hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
        typedef data_schema_traits<vect3d> traits_type ;
        Loci::DatatypeP dp = traits_type::get_type() ;

        hsize_t start = 0 ;
        hsize_t stride = 1 ;
    
        hid_t dataset = H5Dcreate(group_id,"positions",dp->get_hdf5_type(),
                                  dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

    
    
        //first write out pos
        hsize_t count = pos_sizes[0] ;
        if(count != 0) {
          std::vector<vect3d> v_pos(count);
          //put pos_io in a vector
          int index = 0;
          FORALL(file_nodes,cc){
            v_pos[index++] = pos_io[cc];
          }ENDFORALL;
      
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &v_pos[0]) ;
          H5Sclose(memspace) ;
        }
        for(int i=1;i<MPI_processes;++i) {
          start += pos_sizes[i-1] ;
          if(pos_sizes[i] == 0)
            continue ;
          int flag = 0 ;
          MPI_Send(&flag,1,MPI_INT,i,0,MPI_COMM_WORLD) ;
          std::vector<vect3d> rv(pos_sizes[i]) ;
          MPI_Status mstat ;
          MPI_Recv(&rv[0],sizeof(vect3d)*pos_sizes[i],MPI_BYTE,i,1,MPI_COMM_WORLD,
                   &mstat) ;
          count = pos_sizes[i] ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &rv[0]) ;
          H5Sclose(memspace) ;
        }
        //next, write out inner_nodes
        //local edge nodes
        start += pos_sizes[MPI_processes-1] ;
        count = num_local_edge_nodes;
        if(num_local_edge_nodes != 0) {
          std::vector<vect3d> v_nodes(num_local_edge_nodes);
          //put inner_nodes in a vector
      
          long index = 0;
          FORALL(file_edges, cc){
            for(unsigned int i = 0; i < edge_inner_nodes[cc].size(); i++){
              v_nodes[index++] = edge_inner_nodes[cc][i];
            }
          }ENDFORALL;
        
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &v_nodes[0]) ;
          H5Sclose(memspace) ;
        }
        //edge nodes from the other processor
    
        for(int i=1;i<MPI_processes;++i) {
          start += inner_edge_nodes_sizes[i-1] ;
          if(inner_edge_nodes_sizes[i] == 0)
            continue ;
          int flag = 0 ;
          MPI_Send(&flag,1,MPI_INT,i,2,MPI_COMM_WORLD) ;
          std::vector<vector3d<double> > rv(inner_edge_nodes_sizes[i]) ;
          MPI_Status mstat ;
          MPI_Recv(&rv[0],sizeof(vector3d<double> )*inner_edge_nodes_sizes[i],MPI_BYTE,i,3,MPI_COMM_WORLD,
                   &mstat) ;
          count = inner_edge_nodes_sizes[i] ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &rv[0]) ;
          H5Sclose(memspace) ;
        }
    
        //local cell nodes
        start += inner_edge_nodes_sizes[MPI_processes-1] ;
        count = num_local_cell_nodes;
        if(num_local_cell_nodes != 0) {
          std::vector<vector3d<double> > v_nodes(num_local_cell_nodes);
          //put inner_nodes in a vector
      
          long index = 0;
          FORALL(file_cells, cc){
            for(unsigned int i = 0; i < cell_inner_nodes[cc].size(); i++){
              v_nodes[index++] = cell_inner_nodes[cc][i];
            }
          }ENDFORALL;
      
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &v_nodes[0]) ;
          H5Sclose(memspace) ;
        }
        //cell nodes from the other processor
    
        for(int i=1;i<MPI_processes;++i) {
          start += inner_cell_nodes_sizes[i-1] ;
          if(inner_cell_nodes_sizes[i] == 0)
            continue ;
          int flag = 0 ;
          MPI_Send(&flag,1,MPI_INT,i,4,MPI_COMM_WORLD) ;
          std::vector<vector3d<double> > rv(inner_cell_nodes_sizes[i]) ;
          MPI_Status mstat ;
          MPI_Recv(&rv[0],sizeof(vector3d<double> )*inner_cell_nodes_sizes[i],MPI_BYTE,i,5,MPI_COMM_WORLD,
                   &mstat) ;
          count = inner_cell_nodes_sizes[i] ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &rv[0]) ;
          H5Sclose(memspace) ;
        }

        //local face nodes
        start += inner_cell_nodes_sizes[MPI_processes-1] ;
        count = num_local_face_nodes;
        if(num_local_face_nodes != 0) {
          std::vector<vector3d<double> > v_nodes(num_local_face_nodes);
          //put inner_nodes in a vector
      
          long index = 0;
          FORALL(file_faces, cc){
            for(unsigned int i = 0; i < face_inner_nodes[cc].size(); i++){
              v_nodes[index++] = face_inner_nodes[cc][i];
            }
          }ENDFORALL;
      
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &v_nodes[0]) ;
          H5Sclose(memspace) ;
        }
        //face nodes from the other processor
    
        for(int i=1;i<MPI_processes;++i) {
          start += inner_face_nodes_sizes[i-1] ;
          if(inner_face_nodes_sizes[i] == 0)
            continue ;
          int flag = 0 ;
          MPI_Send(&flag,1,MPI_INT,i,6,MPI_COMM_WORLD) ;
          std::vector<vector3d<double> > rv(inner_face_nodes_sizes[i]) ;
          MPI_Status mstat ;
          MPI_Recv(&rv[0],sizeof(vector3d<double> )*inner_face_nodes_sizes[i],MPI_BYTE,i,7,MPI_COMM_WORLD,
                   &mstat) ;
          count = inner_face_nodes_sizes[i] ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                   H5P_DEFAULT, &rv[0]) ;
          H5Sclose(memspace) ;
        }
        start += inner_face_nodes_sizes[MPI_processes-1] ;
        cout << "nodes written " << start << endl;
    
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        //all the other processors
      } else {
        if(local_pos_size != 0){
          std::vector<vector3d<double> > v_pos(local_pos_size);
          //put pos_io in a vector
          int index = 0;
          FORALL(file_nodes,cc){
            v_pos[index++] = pos_io[cc];
          }ENDFORALL;
        
          int flag = 0;
          MPI_Status mstat ;
          MPI_Recv(&flag,1,MPI_INT,0,0,MPI_COMM_WORLD,&mstat) ;
          MPI_Send(&v_pos[0],sizeof(vector3d<double> )*local_pos_size,MPI_BYTE,0,1,MPI_COMM_WORLD) ;
        }
        if(num_local_edge_nodes != 0){
          std::vector<vector3d<double> > v_nodes(num_local_edge_nodes);
          //put inner_nodes in a vector
          int index = 0;
          FORALL(file_edges, cc){
            for(unsigned int i = 0; i < edge_inner_nodes[cc].size(); i++){
              v_nodes[index++] = edge_inner_nodes[cc][i];
            }
          }ENDFORALL;
          int flag = 0;
          MPI_Status mstat ;
          MPI_Recv(&flag,1,MPI_INT,0,2,MPI_COMM_WORLD,&mstat) ;
          MPI_Send(&v_nodes[0],sizeof(vector3d<double> )*num_local_edge_nodes,MPI_BYTE,0,3,MPI_COMM_WORLD) ;
        }
    
        if(num_local_cell_nodes != 0){
          std::vector<vector3d<double> > v_nodes(num_local_cell_nodes);
          //put inner_nodes in a vector
          int index = 0;
          FORALL(file_cells, cc){
            for(unsigned int i = 0; i < cell_inner_nodes[cc].size(); i++){
              v_nodes[index++] = cell_inner_nodes[cc][i];
            }
          }ENDFORALL;
          int flag = 0;
          MPI_Status mstat ;
          MPI_Recv(&flag,1,MPI_INT,0,4,MPI_COMM_WORLD,&mstat) ;
          MPI_Send(&v_nodes[0],sizeof(vector3d<double> )*num_local_cell_nodes,MPI_BYTE,0,5,MPI_COMM_WORLD) ;
        }
    
        if(num_local_face_nodes != 0){
          std::vector<vector3d<double> > v_nodes(num_local_face_nodes);
          //put inner_nodes in a vector
          int index = 0;
          FORALL(file_faces, cc){
            for(unsigned int i = 0; i < face_inner_nodes[cc].size(); i++){
              v_nodes[index++] = face_inner_nodes[cc][i];
            }
          }ENDFORALL;
          int flag = 0;
          MPI_Status mstat ;
          MPI_Recv(&flag,1,MPI_INT,0,6,MPI_COMM_WORLD,&mstat) ;
          MPI_Send(&v_nodes[0],sizeof(vector3d<double> )*num_local_face_nodes,MPI_BYTE,0,7,MPI_COMM_WORLD) ;
        }
      }
  
      if(Loci::MPI_rank == 0) H5Gclose(group_id) ;
  
    }


    void writeVOGNodeP(hid_t file_id,
                       Loci::storeRepP &pos,
                       const_store<Loci::FineNodes> &inner_nodes){ //parallel io
  
#ifndef H5_HAVE_PARALLEL
  
      writeVOGNodeS( file_id,
                     pos,
                     inner_nodes);
#else
   
      hid_t group_id = 0 ;  
      //reorder store first, from local to io entities
      fact_db::distribute_infoP dist =  Loci::exec_current_fact_db->get_distribute_info() ;
      constraint  my_faces, my_geom_cells; 
      entitySet my_entities = dist->my_entities ;
      my_faces = Loci::exec_current_fact_db->get_variable("faces");
      my_geom_cells = Loci::exec_current_fact_db->get_variable("geom_cells");
  
  
  
      entitySet local_edges =my_entities &
        (Loci::exec_current_fact_db->get_variable("edge2node"))->domain(); 
      entitySet local_faces = my_entities & *my_faces;
      entitySet local_cells =  my_entities & *my_geom_cells;
      entitySet local_nodes = my_entities & pos->domain();
  
  

      //before write out,create stores for pos, inner_nodes of cells and faces which
      //are ordered across processors in the file numbering, the domain of this container
      //shifted by offset is the actual file numbering. offset will be modified after function call  
 
      int offset = 0;
      store<vect3d> pos_io;
      pos_io = Loci::Local2FileOrder(pos, local_nodes, offset, dist, MPI_COMM_WORLD) ;
      entitySet file_nodes = pos_io.domain();
  
      offset = 0;
      store<Loci::FineNodes> edge_inner_nodes;
      edge_inner_nodes = Loci::Local2FileOrder(inner_nodes.Rep(),local_edges,offset,dist,MPI_COMM_WORLD) ;
      entitySet file_edges = edge_inner_nodes.domain();
  
      offset= 0;
      // Create container vardist that
      store<Loci::FineNodes> cell_inner_nodes;
      cell_inner_nodes = Loci::Local2FileOrder(inner_nodes.Rep(),local_cells,offset,dist,MPI_COMM_WORLD) ;
      entitySet file_cells = cell_inner_nodes.domain();
  
      offset= 0;
      // Create container vardist that
      store<Loci::FineNodes> face_inner_nodes;
      face_inner_nodes = Loci::Local2FileOrder(inner_nodes.Rep(),local_faces,offset,dist,MPI_COMM_WORLD) ;
      entitySet file_faces = face_inner_nodes.domain();
  
  
      //compute the size of pos
      int local_pos_size = file_nodes.size();
      std::vector<int> pos_sizes(Loci::MPI_processes) ;
      MPI_Allgather(&local_pos_size,1,MPI_INT,
                    &pos_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  
  
  
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
  
  
      std::vector<int> inner_edge_nodes_sizes(Loci::MPI_processes);
      MPI_Allgather(&num_local_edge_nodes,1,MPI_INT,
                    &inner_edge_nodes_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  
  
      std::vector<int> inner_cell_nodes_sizes(Loci::MPI_processes);
      MPI_Allgather(&num_local_cell_nodes,1,MPI_INT,
                    &inner_cell_nodes_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  
      std::vector<int> inner_face_nodes_sizes(Loci::MPI_processes);
      MPI_Allgather(&num_local_face_nodes,1,MPI_INT,
                    &inner_face_nodes_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  




  

      //compute array_size
      long long array_size = 0 ;
      for(int i=0;i<MPI_processes;++i)
        array_size += (pos_sizes[i]+ inner_edge_nodes_sizes[i]+
                       inner_cell_nodes_sizes[i] + inner_face_nodes_sizes[i] );
  
 
  

  
      //first write out numNodes
     
      group_id = H5Gcreate(file_id,"file_info",
                           H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    
      if(MPI_rank == 0) cout << "num_nodes = " << array_size << endl ;
    
      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
      hid_t att_id = H5Acreate(group_id,"numNodes", H5T_STD_I64BE,
                               dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
      H5Awrite(att_id,H5T_NATIVE_LLONG,&array_size) ;
      H5Aclose(att_id) ;
      H5Gclose(group_id) ;
    
      if(array_size == 0)
        return ;

      group_id = H5Gcreate(file_id,"node_info",
                           H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
      //create dataspace and dataset
      int rank = 1 ;
      hsize_t stride = 1 ;
      hsize_t dimension = array_size ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
      typedef data_schema_traits<vect3d> traits_type ;
      Loci::DatatypeP dp = traits_type::get_type() ;


   
    
      hid_t dataset = H5Dcreate(group_id,"positions",dp->get_hdf5_type(),
                                dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
 
      hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
      std::vector<hsize_t> pdispls(Loci::MPI_processes) ; //the start point of each process in prime_comm
      pdispls[0] = 0 ;
      for(int i = 1; i < Loci::MPI_processes; i++) {
        pdispls[i] = pdispls[i-1]+pos_sizes[i-1] ;
      }

    
      //first write out pos
      hsize_t count = pos_sizes[Loci::MPI_rank] ;
  
      if(count != 0) {
        std::vector<vect3d> v_pos(count);
        //put pos_io in a vector
        int index = 0;
        FORALL(file_nodes,cc){
          v_pos[index++] = pos_io[cc];
        }ENDFORALL;
      
        hsize_t start = pdispls[Loci::MPI_rank] ;
    
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;

       
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 xfer_plist, &v_pos[0]) ;
     
        H5Sclose(memspace) ;
      }
  
      //next, write out inner_nodes
      //local edge nodes
  
      pdispls[0] =  pdispls[MPI_processes-1] + pos_sizes[MPI_processes-1] ;
      for(int i = 1; i < Loci::MPI_processes; i++) {
        pdispls[i] = pdispls[i-1]+inner_edge_nodes_sizes[i-1];
      }
      count = num_local_edge_nodes;

    
      if(num_local_edge_nodes != 0) {
        std::vector<vect3d> v_nodes(num_local_edge_nodes);
        //put inner_nodes in a vector
      
        long index = 0;
        FORALL(file_edges, cc){
          for(unsigned int i = 0; i < edge_inner_nodes[cc].size(); i++){
            v_nodes[index++] = edge_inner_nodes[cc][i];
          }
        }ENDFORALL;
    
        hsize_t start = pdispls[Loci::MPI_rank] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
  
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 xfer_plist, &v_nodes[0]) ;
     
        H5Sclose(memspace) ;
      }
   
      //local cell nodes
      pdispls[0] =  pdispls[MPI_processes-1] + inner_edge_nodes_sizes[MPI_processes-1] ;
      for(int i = 1; i < Loci::MPI_processes; i++) {
        pdispls[i] = pdispls[i-1]+inner_cell_nodes_sizes[i-1];
      }
      count = num_local_cell_nodes;

  
    
   
      if(num_local_cell_nodes != 0) {
        std::vector<vector3d<double> > v_nodes(num_local_cell_nodes);
        //put inner_nodes in a vector
      
        long index = 0;
        FORALL(file_cells, cc){
          for(unsigned int i = 0; i < cell_inner_nodes[cc].size(); i++){
            v_nodes[index++] = cell_inner_nodes[cc][i];
          }
        }ENDFORALL;

        hsize_t start = pdispls[Loci::MPI_rank] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 xfer_plist, &v_nodes[0]) ;
        H5Sclose(memspace) ;
      }
    
   

      //local face nodes

      pdispls[0] =  pdispls[MPI_processes-1] + inner_cell_nodes_sizes[MPI_processes-1] ;
      for(int i = 1; i < Loci::MPI_processes; i++) {
        pdispls[i] = pdispls[i-1]+inner_face_nodes_sizes[i-1];
      }
      count = num_local_face_nodes;
      if(num_local_face_nodes != 0) {
        std::vector<vector3d<double> > v_nodes(num_local_face_nodes);
        //put inner_nodes in a vector
      
        long index = 0;
        FORALL(file_faces, cc){
          for(unsigned int i = 0; i < face_inner_nodes[cc].size(); i++){
            v_nodes[index++] = face_inner_nodes[cc][i];
          }
        }ENDFORALL;
        hsize_t start = pdispls[Loci::MPI_rank] ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
        H5Dwrite(dataset,dp->get_hdf5_type(),memspace,dataspace,
                 xfer_plist, &v_nodes[0]) ;
        H5Sclose(memspace) ;
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Pclose(xfer_plist) ;
      H5Gclose(group_id) ;
#endif
    }
  }

  
  void writeVOGNode(hid_t file_id,
                    Loci::storeRepP &pos,
                    const_store<Loci::FineNodes> &inner_nodes){
    if(use_parallel_io)pio::writeVOGNodeP( file_id,pos,inner_nodes);
    else pio::writeVOGNodeS( file_id,pos,inner_nodes);
  }

  
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

void colorMatrix(Map &cl, Map &cr, multiMap &face2node) {
    
  entitySet  faces,cells ;
  std::vector<entitySet> ptn = getDist(faces,cells,
                                       cl,cr,face2node);
  entitySet loc_faces = faces & ptn[Loci::MPI_rank] ;
  entitySet geom_cells = cells & ptn[Loci::MPI_rank] ;
  entitySet negs = interval(UNIVERSE_MIN,-1) ;
  entitySet boundary_faces = cr.preimage(negs).first ;
  entitySet interior_faces = loc_faces - boundary_faces ;

  using std::pair ;
  std:: vector<pair<Entity,Entity> > cellmap(interior_faces.size()*2) ;
  int cnt = 0 ;
  FORALL(interior_faces,fc) {
    cellmap[cnt++] = pair<Entity,Entity>(cl[fc],cr[fc]) ;
    cellmap[cnt++] = pair<Entity,Entity>(cr[fc],cl[fc]) ;
  } ENDFORALL ;
  multiMap c2c ;
  Loci::distributed_inverseMap(c2c,cellmap,cells,cells,ptn) ;
  int ncells = cells.size() ;

  store<int> ctmp ;
  ctmp.allocate(geom_cells) ;
  FORALL(geom_cells,cc) {
    ctmp[cc] = -1 ;
  } ENDFORALL ;

  int col = ncells*Loci::MPI_rank ;
    
  std:: vector<int> visited ;
  entitySet left_out = geom_cells ;
  int lo_p = geom_cells.Min() ;
  while(left_out != EMPTY) {
    std:: vector<int> work ;
    work.push_back(left_out.Min()) ;
    while(work.size() != 0) {
      std::vector<int> working ;
      for(size_t i=0;i<work.size();++i) {
        int cc = work[i] ;
        if(ctmp[cc] == -1) {
          ctmp[cc] = col++ ;
          visited.push_back(cc) ;
          for(const int *pi = c2c.begin(cc);pi!=c2c.end(cc);++pi)
            if(geom_cells.inSet(*pi) && ctmp[*pi] == -1) {
              working.push_back(*pi) ;
            }
        }
      }
      work.swap(working) ;
    }
    left_out = EMPTY ;
    entitySet candidates = geom_cells & interval(lo_p,UNIVERSE_MAX) ;
    FORALL(candidates,cc) {
      if(ctmp[cc] == -1) {
        left_out += cc ;
        break ;
      }
    } ENDFORALL ;
    if(left_out != EMPTY)
      lo_p = left_out.Min() ;
  }

  dstore<int> color ;
  FORALL(geom_cells,cc) {
    color[cc] = ctmp[cc];
  } ENDFORALL ;

  entitySet clone_cells = cl.image(interior_faces)
    + cr.image(interior_faces) ;
  clone_cells -= geom_cells ;
  Loci::storeRepP cp_sp = color.Rep() ;
  Loci::fill_clone(cp_sp, clone_cells, ptn) ;

  FORALL(interior_faces,fc) {
    int color_l = color[cl[fc]] ;
    int color_r = color[cr[fc]] ;
    if(color_l == color_r) 
      cerr << "color equal == " << color_l << endl ;
    if((color_l == -1) || (color_r == -1))
      cerr << "matrix coloring internal error" << endl ;
                                                              
    if(color_l > color_r) {
      // change face orientation to match matrix coloring
      std::swap(cl[fc],cr[fc]) ;
      int i = 0 ;
      int j = face2node[fc].size() - 1;
      while(i < j) {
        std::swap(face2node[fc][i],face2node[fc][j]) ;
        i++ ;
        j-- ;
      } 
    }
  } ENDFORALL ;
    

}
