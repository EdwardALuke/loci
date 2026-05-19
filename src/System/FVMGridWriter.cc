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
#include <store.h>
#include <DStore.h>
#include <Map.h>
#include <DMap.h>
#include <multiMap.h>
#include <DMultiMap.h>
#include <constraint.h>
#include <distribute.h>
#include <distribute_container.h>
#include <parameter.h>
#include <fact_db.h>
#include <Loci_types.h>
#include <LociGridReaders.h>
#include "loci_globs.h"
#include <sstream>
#include <distribute_io.h>
#include <Tools/tools.h>


#include <map>

#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
using std::cout ;
using std::endl ;
using std::cerr ;
using std::ifstream ;
using std::ios ;

namespace Loci {
  // Utility routines for creating face cluster in vog file output
  void writeUnsignedVal(vector<unsigned char>& cluster, unsigned long long val) {
    do {
      unsigned char byte = val & 0x7f ;
      val = val >> 7 ;
      if(val != 0) {
        byte |= 0x80 ;
      }
      cluster.push_back(byte) ;
    } while(val != 0) ;
  }


  void writeSignedVal(vector<unsigned char> &cluster, long long val) {
    bool sign = false ;
    if(val < 0) {
      sign = true ;
      val = -val ;
    }
    unsigned char byte = val & 0x3f ;
    if(sign)
      byte |= 0x40 ;
    val = val >> 6 ;
    if(val != 0)
      byte |= 0x80 ;
    cluster.push_back(byte) ;
    if((byte & 0x80) == 0x80)
      writeUnsignedVal(cluster,val) ;
  }

  // Writes the table for node or cell id lookup by writing offsets
  // using adjustable precision integers
  void writeTable(vector<unsigned char> &cluster, entitySet set) {
    entitySet::const_iterator ei = set.begin() ;
    unsigned char sz = set.size() ;
    cluster.push_back(sz) ;
    writeSignedVal(cluster,*ei) ;
    long long last = *ei ;
    for(++ei;ei!=set.end();++ei) {
      unsigned long diff = *ei - last ;
      last = *ei ;
      writeUnsignedVal(cluster,diff) ;
    }
  }

  // Encode a face cluster for a vog file
  vector<unsigned char>
  encode_face_cluster(const multiMap &face2node,
                      const Map &cl, const Map &cr,
                      entitySet fcluster,
                      entitySet nodeSet,
                      entitySet cellSet) {
    vector<unsigned char> cluster ;



    std::map<int,int> node2local ;
    std::map<int,int> cell2local ;
    int cnt = 0 ;
    entitySet::const_iterator ei ;
    for(ei = nodeSet.begin();ei!=nodeSet.end();++ei) {
      node2local[*ei] = cnt++ ;
    }
    cnt = 0 ;
    for(ei = cellSet.begin();ei!=cellSet.end();++ei) {
      cell2local[*ei] = cnt++ ;
    }

    // Sort faces according to number of nodes
    vector<pair<int,Entity> > face_order(fcluster.size()) ;
    cnt = 0 ;
    for(ei=fcluster.begin();ei!=fcluster.end();++ei)
      face_order[cnt++] = pair<int,Entity>(face2node[*ei].size(),*ei) ;
    std::sort(face_order.begin(),face_order.end()) ;
    vector<pair<int,int> > rll ;
    int lsz = face_order[0].first ;
    cnt = 0 ;
    for(size_t i=0;i<face_order.size();++i) {
      if(lsz!=face_order[i].first) {
        while(cnt > 255) {
          rll.push_back(pair<int,int>(lsz,255)) ;
          cnt -= 255 ;
        }
        rll.push_back(pair<int,int>(lsz,cnt)) ;
        cnt = 0 ;
        lsz = face_order[i].first ;
      }
      cnt++ ;
    }
    while(cnt > 255) {
      rll.push_back(pair<int,int>(lsz,255)) ;
      cnt -= 255 ;
    }
    rll.push_back(pair<int,int>(lsz,cnt)) ;

    // Now write out the faces for each size category
    cnt = 0 ;
    for(size_t i=0;i<rll.size();++i) {
      cluster.push_back(rll[i].first) ;
      cluster.push_back(rll[i].second) ;
      int nds = rll[i].first ;
      for(int k=0;k<rll[i].second;++k) {
        int fc = face_order[cnt].second ;
        cnt++ ;
      
        for(int j=0;j<nds;++j)
          cluster.push_back(node2local[face2node[fc][j]]) ;
        cluster.push_back(cell2local[cl[fc]]) ;
        cluster.push_back(cell2local[cr[fc]]) ;
      }
    }
    // A zero face size marks end of cluster
    cluster.push_back(0) ;

    writeTable(cluster,nodeSet) ;
    writeTable(cluster,cellSet) ;
    // Cluster finished,return ;
    return cluster ;
  }
  
  entitySet faceCluster(const multiMap &face2node,
                        const Map &cl, const Map &cr, entitySet faces,
                        vector<unsigned char> &cluster_info,
                        vector<unsigned short> &cluster_sizes) {
    entitySet faceSet ;
    entitySet nodeSet ;
    entitySet cellSet ;
    entitySet fcluster ;
    
    int nnodes = 0 ;
    int ncells = 0 ;
    entitySet::const_iterator ei ;
    for(ei = faces.begin();ei!=faces.end();++ei) {
      int ncells_local = 0 ;
      int nnodes_local = 0 ;
      Entity fc = *ei ;
      if(!cellSet.inSet(cl[fc]))
        ncells_local++ ;
      if(!cellSet.inSet(cr[fc]))
        ncells_local++ ;
      int sz = face2node[fc].size() ;
      for(int i=0;i<sz;++i)
        if(!nodeSet.inSet(face2node[fc][i]))
          nnodes_local++ ;
      if(nnodes +nnodes_local >256 ||
         ncells +ncells_local > 256)
        break ;
      cellSet += cl[fc] ;
      cellSet += cr[fc] ;
      for(int i=0;i<sz;++i)
        nodeSet += face2node[fc][i] ;
      nnodes = nodeSet.size() ;
      ncells = cellSet.size() ;
      fcluster += fc ;
    }
    vector<unsigned char> cluster =
      encode_face_cluster(face2node,cl,cr, fcluster, nodeSet, cellSet) ;

    int cluster_size = cluster.size() ;
    cluster_sizes.push_back(cluster_size) ;
    for(int i=0;i<cluster_size;++i)
      cluster_info.push_back(cluster[i]) ;

    return fcluster ;
  }



  
  void writeVOGSurf(hid_t file_id, vector<pair<int,string> > surface_ids) {
    /*
      parallel io included, all attribute and group creations need to be called collectively by all ranks, including H5Awrite
    */
    hid_t group_id = 0 ;
    if(surface_ids.size() == 0) return;
    if(MPI_rank == 0 || use_parallel_io) {
      
      vector<pair<int,string> > surface_ids_mod= surface_ids ; 
      for(size_t i=0;i<surface_ids.size();++i) {
        string from = surface_ids_mod[i].second ;
        if((from.size()==0) || !((from[0] >='a' && from[0] <='z') ||
                                 (from[0] >='A' && from[0] <='Z'))) {
          from = "_"+surface_ids_mod[i].second ;
        }
        for(size_t j=0;j<from.size();++j) {
          if(!((from[j] >='a' && from[j] <= 'z') ||
               (from[j] >='A' && from[j] <= 'Z') ||
               (from[j] >='0' && from[j] <= '9'))) {
            from[j] = '_' ;
          }
        }
        surface_ids_mod[i].second = from ;
      }
      
      group_id = H5Gcreate(file_id,"surface_info",
                           H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

      for(size_t i=0;i<surface_ids.size();++i) {
        hid_t bc_id = 0 ;

        bc_id = H5Gcreate(group_id,surface_ids_mod[i].second.c_str(),
			  H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT) ;
        hsize_t dims = 1 ;
        hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
          
        hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                                 dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
         
        H5Awrite(att_id,H5T_NATIVE_INT,&surface_ids_mod[i].first) ;
        H5Aclose(att_id) ;
        H5Gclose(bc_id) ;
      }
      H5Gclose(group_id) ;
    }
    
  }
  
  // originally, no MPI_rank ==0,  can not run in parallel. Modified so that writeVOG can run in parallel io
  void writeVOGTag(hid_t output_fid,  vector<pair<string,entitySet> >& volTags){
    if(MPI_rank == 0 || (use_parallel_io && output_fid > 0)){

      hid_t cell_info = H5Gcreate(output_fid,"cell_info", 
                                  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

      for(size_t i=0;i<volTags.size();++i) {

        hid_t vol_id = H5Gcreate(cell_info,volTags[i].first.c_str(),
                                 H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
        hsize_t dims = 1 ;
        hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
      
        hid_t att_id = H5Acreate(vol_id,"Ident", H5T_NATIVE_INT,
                                 dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;

        int num = int(i) ;
        H5Awrite(att_id,H5T_NATIVE_INT,&num) ;
        H5Aclose(att_id) ;
        
        Loci::HDF5_WriteDomain(vol_id,volTags[i].second,
                               get_exec_comm()) ;
        
      
        H5Gclose(vol_id) ;
      }
      H5Gclose(cell_info) ;
    }
  }
  
  void writeVOGNode(hid_t file_id, store<vector3d<double> > &pos) {//parallel io included
    hid_t group_id = 0 ;
    if(MPI_rank == 0 || use_parallel_io){

      group_id = H5Gcreate(file_id,"node_info",H5P_DEFAULT,
                           H5P_DEFAULT,H5P_DEFAULT) ;

    }
    // Write out node info
    entitySet nodes = pos.domain() ;
    vector<vector3d<double> > vpos(nodes.size()) ;
    int cnt = 0 ;
    entitySet::const_iterator ei ;
    for(ei=nodes.begin();ei!=nodes.end();++ei)
      vpos[cnt++] = pos[*ei] ;

    // if(use_parallel_io)
    //   writeUnorderedVectorP(group_id,"positions",vpos) ;
    // else
    writeUnorderedVector(group_id,"positions",vpos) ;

    if(MPI_rank == 0 || use_parallel_io) H5Gclose(group_id) ;

    long long local_num_nodes = pos.domain().size() ;
    long long num_nodes = 0 ;
    MPI_Allreduce(&local_num_nodes,&num_nodes,1,MPI_LONG_LONG_INT,
                  MPI_SUM, get_exec_comm()) ;

    if(MPI_rank == 0 || use_parallel_io) {
      group_id = H5Gcreate(file_id,"file_info",
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    
      if(MPI_rank == 0) cerr << "num_nodes = " << num_nodes << endl ;

      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
     
      hid_t att_id = H5Acreate(group_id,"numNodes", H5T_STD_I64BE,
                               dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;

      H5Awrite(att_id,H5T_NATIVE_LLONG,&num_nodes) ;
      H5Aclose(att_id) ;
      H5Gclose(group_id) ;
    }
  }

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
                  MPI_SUM, get_exec_comm()) ;

    hid_t group_id = 0 ;
    if(MPI_rank == 0 || use_parallel_io) {
      group_id = H5Gopen(file_id,"file_info",H5P_DEFAULT) ;

      if(MPI_rank == 0)cerr << "num_cells = " << num_cells << endl
                            << "num_faces = " << num_faces << endl ;

      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;

      hid_t att_id = H5Acreate(group_id,"numFaces", H5T_STD_I64BE,
                               dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;

      H5Awrite(att_id,H5T_NATIVE_LLONG,&num_faces) ;
      H5Aclose(att_id) ;

      att_id = H5Acreate(group_id,"numCells", H5T_STD_I64BE,
                         dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;

      H5Awrite(att_id,H5T_NATIVE_LLONG,&num_cells) ;
      H5Aclose(att_id) ;
      H5Gclose(group_id) ;

      group_id = H5Gcreate(file_id,"face_info",
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

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
      entitySet fcluster = faceCluster(tmp_face2node,tmp_cl,tmp_cr,faces,
                                       cluster_info,cluster_sizes) ;
      faces -= fcluster ;
    }
    
    
    Loci::writeUnorderedVector(group_id,"cluster_sizes",cluster_sizes) ;
    Loci::writeUnorderedVector(group_id,"cluster_info",cluster_info) ;
      
    
    if(MPI_rank == 0 || use_parallel_io) {
      H5Gclose(group_id) ;
    }
  }

  void writeVOGClose(hid_t file_id) {//parallel io included
    if(MPI_rank == 0 || use_parallel_io) H5Fclose(file_id) ;
  }
  



  void writeVOG(string filename,store<vector3d<double> > &pos,
                Map &cl, Map &cr, multiMap &face2node,
                vector<pair<int,string> > surface_ids) {
    // write grid file
    hid_t file_id = 0;
   
    file_id = writeVOGOpen(filename) ;
    
    writeVOGSurf(file_id,surface_ids) ;
    writeVOGNode(file_id,pos) ;
    writeVOGFace(file_id,cl,cr,face2node) ;
    writeVOGClose(file_id) ; //parallel io included
  }

  void writeVOG(string filename,store<vector3d<double> > &pos,
                Map &cl, Map &cr, multiMap &face2node,
                vector<pair<int,string> >& surface_ids,
                vector<pair<string,entitySet> >& volTags ) {
    // write grid file
    hid_t file_id = writeVOGOpen(filename) ;
   
    writeVOGTag(file_id, volTags) ;
    writeVOGSurf(file_id,surface_ids) ;
    writeVOGNode(file_id,pos) ;
    writeVOGFace(file_id,cl,cr,face2node) ;
    writeVOGClose(file_id) ; //parallel io included
  }
  

}
