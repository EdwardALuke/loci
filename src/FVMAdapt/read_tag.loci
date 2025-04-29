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
//*****************************************************************************
// this file read in the refinement plan 
//moified 03/12/07, read in cell tags before face tags
//*****************************************************************************
#include <iostream>
#include <fstream>
#include <string>
#include <Loci.h>
#include <vector>
#include <Tools/tools.h>
#include "sciTypes.h"
#include "defines.h"
#include <hdf5.h>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using Loci::storeRepP ;


//functions from distribute_io
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

// Read in text version of tags file
vector<char> read_tags_txt(string filename,MPI_Comm &comm) {
  int p = 1 ;
  MPI_Comm_size(comm, &p) ;
  int r = 0 ;
  MPI_Comm_rank(comm, &r) ;
  if(r == 0) { // read in the file here
    // First scan the file to find how many items are in there
    ifstream inFileScan;
    inFileScan.open(filename.c_str(),std::ios::in);
    if(!inFileScan){
      cerr <<"cannot open '" << filename << "' for input" << endl;
      Loci::Abort();
    }

    int cnt = 0 ;
    for(;;) {
      int tmp =0;
      inFileScan >> tmp ;
      if(inFileScan.eof() || inFileScan.fail())
        break ;
      cnt++ ;
    }
    inFileScan.close();
    vector<int> proc_alloc(p,cnt/p) ;
    for(int i=0;i<p-1;++i)
      cnt -= proc_alloc[i] ;
    proc_alloc[p-1] = cnt ;
    int size ;
    MPI_Scatter(&proc_alloc[0],1,MPI_INT,&size,1,MPI_INT,0,comm) ;
    ifstream inFile ;
    inFile.open(filename.c_str(),std::ios::in) ;
    vector<char> retval(size) ;
    for(int i=0;i<size;++i) {
      int tmp ;
      inFile >> tmp ;
      retval[i] = ((0==tmp)?0:((tmp>0)?1:2)) ;
    }
    for(int k=1;k<p;++k) {
      vector<char> sendval(proc_alloc[k]) ;
      for(int i=0;i<size;++i) {
        int tmp ;
        inFile >> tmp ;
        sendval[i] = ((0==tmp)?0:((tmp>0)?1:2)) ;
      }
      MPI_Send(&sendval[0],proc_alloc[k],MPI_CHAR,k,111,comm) ;
    }
    return retval ;
  } else {
    int size ;
    MPI_Scatter(0,0,MPI_INT,&size,1,MPI_INT,0,comm) ;
    vector<char> retval(size) ;
    MPI_Status status ;
    MPI_Recv(&retval[0],size,MPI_CHAR,0,111,comm,&status) ;
    return retval ;
  }
}

// read in hdf5 version of tags file
vector<char> read_tags_hdf5(string filename, string varname, MPI_Comm &comm) {
  int p = 1 ;
  MPI_Comm_size(comm, &p) ;
  int r = 0 ;
  MPI_Comm_rank(comm, &r) ;
  if(r == 0) { // read in the file here
    hid_t file_id = 0;
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group_id = 0;
    group_id = H5Gopen(file_id, varname.c_str(),H5P_DEFAULT) ;
    //process 0 read in its local tag
    hid_t dataset =  H5Dopen(group_id, "data",H5P_DEFAULT) ;
    hid_t dataspace = H5Dget_space(dataset) ;
    // Get size of tag file
    hsize_t size = 0 ;
    H5Sget_simple_extent_dims(dataspace,&size,NULL) ;
    int cnt = size ;
    vector<int> proc_alloc(p,cnt/p) ;
    for(int i=0;i<p-1;++i)
      cnt -= proc_alloc[i] ;
    proc_alloc[p-1] = cnt ;
    int lsize ;
    MPI_Scatter(&proc_alloc[0],1,MPI_INT,&lsize,1,MPI_INT,0,comm) ;
    
    hsize_t start = 0 ;
    hsize_t stride = 1 ;
    hsize_t count = 0 ;
    int rank = 1 ;
    hsize_t dimension = proc_alloc[0]; // size of local data
    count = dimension ;
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    vector<float> buf(proc_alloc[0]) ;
    H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
             H5P_DEFAULT, &buf[0]);
    H5Sclose(memspace) ;
    vector<char> retval(proc_alloc[0]) ;
    for(int i=0;i<proc_alloc[0];++i) {
      float tmp = buf[i] ;
      retval[i] = ((0.0==tmp)?0:((tmp>0.0)?1:2)) ;
    }
    start += proc_alloc[0] ;
    for(int k=1;k<p;++k) {
      int send_size = proc_alloc[k] ;
      vector<char> sendval(send_size) ;
      dimension = send_size ;
      count = dimension ;
      
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      vector<float> buf(send_size) ;
      H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
               H5P_DEFAULT, &buf[0]);
      H5Sclose(memspace) ;
      start += count;

      for(int i=0;i<send_size;++i) {
        float tmp = buf[i] ;
        sendval[i] = ((0.0==tmp)?0:((tmp>0.0)?1:2)) ;
      }
      MPI_Send(&sendval[0],proc_alloc[k],MPI_CHAR,k,111,comm) ;
    }

    H5Dclose(dataset) ;
    H5Sclose(dataspace) ;
    H5Gclose(group_id) ;
    H5Fclose(file_id);
    return retval ;
  } else {
    int size ;
    MPI_Scatter(0,0,MPI_INT,&size,1,MPI_INT,0,comm) ;
    vector<char> retval(size) ;
    MPI_Status status ;
    MPI_Recv(&retval[0],size,MPI_CHAR,0,111,comm,&status) ;
    return retval ;
  }
}


// Add requests for node to datamap, datamap here is a file2local map
//but file number is shifted by offset so that it will start with 0
void add_to_datamap(entitySet dom, vector<pair<int,int> > &datamap, int offset,
                    fact_db::distribute_infoP dist) {
  if(dist == 0) {
    int cnt = 0 ;
    FORALL(dom,nd) {
      pair<int,int> item(cnt,nd) ;
      datamap.push_back(item) ;
      cnt++ ;
    } ENDFORALL ;
  } else {
    dMap g2f ;
    g2f = dist->g2fv[0].Rep() ; // FIX THIS
    Map l2g ;
    l2g = dist->l2g.Rep() ;
    FORALL(dom,nd) {
      pair<int,int> item(g2f[l2g[nd]]-offset,nd) ;
      datamap.push_back(item) ;
    } ENDFORALL ;
  }      
}

// Take datamap input and return data from the tags file for the
// nodes requested along with the index of the item that made the
// request
void communicateFileData(vector<char> &filedata,
                         vector<pair<int,int> > &datamap,
                         vector<pair<int,char> > &returndata,
                         MPI_Comm &comm) {
  int p = 1 ;
  MPI_Comm_size(comm, &p) ;
  int r = 0 ;
  MPI_Comm_rank(comm, &r) ;
  
  int dmsz = datamap.size() ;
  if(p == 1) {
    // On one processor then we just need to copy data from filedata
    // to tag using datamap
    {
      vector<pair<int,char> > tmp(dmsz) ;
      returndata.swap(tmp) ;
    }
    for(int i=0;i<dmsz;++i) {
      returndata[i].first = datamap[i].second ;
      returndata[i].second = filedata[datamap[i].first] ;
    }
    return ;
  }
  sort(datamap.begin(),datamap.end()) ;
 
  vector<int> sizes(p) ;
  int lsize = filedata.size() ;
  MPI_Allgather(&lsize,1,MPI_INT,&sizes[0],1,MPI_INT,comm) ;
  vector<int> offsets(p+1) ;
  offsets[0] = 0 ;
  for(int i=0;i<p;++i)
    offsets[i+1] = offsets[i]+sizes[i] ;
  vector<int> datamapoffsets(p+1,dmsz) ;
  datamapoffsets[0] =0 ;

  int cnt=1 ;
  for(int i=0;cnt < p && i<=dmsz;++i) {
    if(i == dmsz)
      datamapoffsets[cnt] = i ;
    else
      while(datamap[i].first >= offsets[cnt]) {
        datamapoffsets[cnt] = i ;
        if(cnt == p)
          break ;
        cnt++ ;
      }
  }
  datamapoffsets[p] = dmsz ;

  // Compute data requests for file data by processor
  vector<int> reqsizes(p) ;
  for(int i=0;i<p;++i)
    reqsizes[i] = datamapoffsets[i+1]-datamapoffsets[i] ;
  // convert requested item to local numbering on that processor.
  vector<int> reqitem(dmsz) ;
  for(int i=0;i<p;++i)
    for(int j=datamapoffsets[i];j<datamapoffsets[i+1];++j) {
      reqitem[j] = datamap[j].first-offsets[i] ;
    }
  
  // send request sizes
  vector<int> sendsizes(p) ;
  MPI_Alltoall(&reqsizes[0],1,MPI_INT,&sendsizes[0],1,MPI_INT, comm) ;
  
  // Now send item requests
  int nreq = sendsizes[0] ;
  for(int i=1;i<p;++i)
    nreq+=sendsizes[i] ;

  vector<int> recvreq(nreq) ;
  vector<int> sendoffsets(p+1) ;
  sendoffsets[0] = 0 ;
  for(int i=0;i<p;++i)
    sendoffsets[i+1] = sendoffsets[i]+sendsizes[i] ;
  MPI_Alltoallv(&reqitem[0],&reqsizes[0],&datamapoffsets[0],MPI_INT,
                &recvreq[0],&sendsizes[0],&sendoffsets[0],MPI_INT,
                comm) ;
  vector<char> sendinfo(nreq) ;
  for(int i=0;i<nreq;++i) {
    sendinfo[i] = filedata[recvreq[i]] ;
  }
  vector<char> recvinfo(dmsz) ;
  // then return requested items
  MPI_Alltoallv(&sendinfo[0],&sendsizes[0],&sendoffsets[0],MPI_CHAR,
                &recvinfo[0],&reqsizes[0],&datamapoffsets[0],MPI_CHAR,
                comm) ;

  // then store requested item in output
  {
    vector<pair<int,char> > tmp(dmsz) ;
    returndata.swap(tmp) ;
  }
  for(int i=0;i<dmsz;++i) {
    returndata[i].first = datamap[i].second ;
    returndata[i].second = recvinfo[i] ;
  }
    
}

class input_TagsData : public singleton_rule {
  const_param<std::string> tagfile_par ;
  blackbox<vector<char> > inputTagsData ;
public:
  input_TagsData() {
    name_store("tagfile_par",tagfile_par) ;
    name_store("inputTagsData",inputTagsData) ;
    input("tagfile_par") ;
    output("inputTagsData") ;
    disable_threading() ;
  }
  virtual void compute(const sequence &seq) {
    string filename = *tagfile_par;
    MPI_Comm comm = MPI_COMM_WORLD ;    

    size_t p = filename.find("ref_sca");
    if(p != string::npos){//read hdf5 file
      *inputTagsData = read_tags_hdf5(filename,string("ref"),comm) ;
    } else {
      *inputTagsData = read_tags_txt(filename,comm) ;
    }
  }
} ;

register_rule<input_TagsData> register_input_TagsData ;
  
class get_postag : public pointwise_rule{
  const_blackbox<vector<char> > inputTagsData ;
  store<char> posTag;
  
public:
  get_postag(){
    name_store("inputTagsData", inputTagsData);
    name_store("posTag", posTag);
    input("inputTagsData");
    output("posTag");
    constraint("pos");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    MPI_Comm comm = MPI_COMM_WORLD ;    
    
    vector<char> filedata  = *inputTagsData ;
    
    // Now determine what data we need to collect from filedata
    fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
    entitySet dom = entitySet(seq);

    vector<pair<int,int> > datamap ;
    add_to_datamap(dom,datamap,0,dist) ;
    
   
    vector<pair<int,char> > returndata ;
    communicateFileData(filedata,datamap,returndata,comm) ;
    
    for(size_t i=0;i<returndata.size();++i) {
      posTag[returndata[i].first] = returndata[i].second ;
    }

    
  }

};
register_rule<get_postag> register_get_postag;


  
int getFileNumberOffset(const entitySet& locdom,MPI_Comm &comm){
  // Now determine what data we need to collect from filedata
  fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
  int imn;
  if(dist == 0) {
    return locdom.Min();
  } else {
    entitySet dom = ~EMPTY;
    dom = locdom&(dist->my_entities) ;
    Map l2g ;
    l2g = dist->l2g.Rep() ;
    store<unsigned char> key_domain ;
    key_domain = dist->key_domain.Rep() ;
    // // Compute domain in global numbering
    //     entitySet dom_global = l2g.image(dom) ;
     
    //     // This shouldn't happen
    //     FATAL(dom.size() != dom_global.size()) ;

    // Now get global to file numbering
    //    dMap g2f ;
    //    g2f = dist->g2f.Rep() ;

    // Compute map from local numbering to file numbering
    entitySet filedom = EMPTY;
    FORALL(dom,i) {
      int kd = key_domain[i] ;
      filedom += dist->g2fv[kd][l2g[i]] ;
    } ENDFORALL ;

    // int imx = std::numeric_limits<int>::min() ;
    // int imn = std::numeric_limits<int>::max() ;

    // // Find bounds in file numbering from this processor
    //     FORALL(dom,i) {
    //       imx = max(newnum[i],imx) ;
    //       imn = min(newnum[i],imn) ;
    //     } ENDFORALL ;

    imn = filedom.Min();
    // Find overall bounds
    //imx = GLOBAL_MAX(imx) ;
    int local_min = imn;
    MPI_Allreduce(&local_min, &imn, 1, MPI_INT, MPI_MIN,comm);  
    return imn;
  }
  return imn;
}

  

// Return node location for each entity based on how many nodes are
// created through splits input in numsplits, numsplits is accessed
// over locdom in the local numbering
int
getNodeOffsets(store<int> &nodeloc, const_store<int> &numsplits,
               entitySet locdom,
               fact_db::distribute_infoP dist,
               MPI_Comm &comm) {
  entitySet my_entities = ~EMPTY ;
  if(dist != 0)
    my_entities = (dist->my_entities) ;
  locdom = locdom & my_entities ;
  // Number of splits in file numbering
  store<int> fo_numsplits ;
  int offset = 0 ;
  if(dist == 0) {
    fo_numsplits.allocate(locdom) ;
    FORALL(locdom,ii) {
      fo_numsplits[ii] = numsplits[ii] ;
    } ENDFORALL ;
  } else {
    fo_numsplits = Local2FileOrder(numsplits.Rep(),locdom,offset,dist,comm) ;
  }
  int psplits = 0 ;
  Loci::entitySet filedom = fo_numsplits.domain() ;
  FORALL(filedom,pt) {
    psplits += fo_numsplits[pt] ;
  } ENDFORALL ;
  int p = 1 ;
  MPI_Comm_size(comm, &p) ;
  int r = 0 ;
  MPI_Comm_rank(comm, &r) ;
  vector<int> gather_sizes(p) ;
  MPI_Allgather(&psplits,1,MPI_INT,&gather_sizes[0],1,MPI_INT,comm) ;
  vector<int> fo_offsets(p) ;
  fo_offsets[0] = 0 ;
  for(int i=1;i<p;++i)
    fo_offsets[i] = fo_offsets[i-1]+gather_sizes[i-1] ;
  int tnodes = fo_offsets[p-1]+gather_sizes[p-1] ;
  int tot = fo_offsets[r] ;
  FORALL(filedom,pt) {
    int sloc = fo_numsplits[pt] ;
    fo_numsplits[pt] = tot ;
    tot += sloc ;
  } ENDFORALL ;

  if(dist == 0)
    nodeloc = fo_numsplits.Rep() ;
  else {
    nodeloc.allocate(locdom) ;
    storeRepP localVar=nodeloc.Rep();
    File2LocalOrder(localVar, locdom,
                    fo_numsplits.Rep(), offset,
                    dist,
                    MPI_COMM_WORLD);
  }
  return tnodes ;
}
  
    
class get_nodetag : public pointwise_rule{
  const_blackbox<vector<char> > inputTagsData;
  const_param<int> num_original_nodes;
  const_store<int> num_inner_nodes;
  store<std::vector<char> > nodeTag;
public:
  get_nodetag(){
    name_store("inputTagsData", inputTagsData);
    name_store("num_original_nodes", num_original_nodes);
    name_store("nodeTag", nodeTag);
    name_store("num_inner_nodes", num_inner_nodes);
    input("inputTagsData");
    input("num_original_nodes");
    input("num_inner_nodes");
    output("nodeTag");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    vector<char> filedata  = *inputTagsData ;

    MPI_Comm comm = MPI_COMM_WORLD ;    
    
    // Now determine what data we need to collect from filedata
    fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
    entitySet my_entities ;
    my_entities = ~EMPTY ;

    if(dist != 0)
      my_entities = (dist->my_entities) ;

    int start = *num_original_nodes; //skip the posTag part

    Loci::constraint geom_cells, faces;
    Loci::storeRepP e2n = (Loci::exec_current_fact_db)->get_variable("edge2node");
    faces = (Loci::exec_current_fact_db)->get_variable("faces");
    geom_cells = (Loci::exec_current_fact_db)->get_variable("geom_cells");

    entitySet local_edges, local_faces, local_cells;
    //don't know if it's necessray
    local_edges = (my_entities) & (e2n->domain()) ;
    local_faces = (my_entities) & (*faces);
    local_cells = (my_entities)&(*geom_cells);
    
    // Gather data maps for nodes created at edges, then cells, then faces
    // Note, since each of these are given a file number starting from zero
    // and we are using local2file to get a consistent numbering we need
    // to process each of these steps one by one.
    vector<pair<int,int> > datamap ;
    vector<pair<int,int> > locations ; ;
    store<int> nodeloc ; // location of nodes

    // get offsets from the beginning of the edge node for each referenced
    // node
    int nenodes =
      getNodeOffsets(nodeloc,num_inner_nodes,local_edges,dist,comm) ;
    
    // allocate nodeTag
    int cnt = 0 ;
    FORALL(local_edges, cc){
      if(num_inner_nodes[cc] != 0){
        int nn = num_inner_nodes[cc] ;
        std::vector<char>(nn).swap(nodeTag[cc]);
        pair<int,int> loc(cc,0) ;
        pair<int,int> data(start+nodeloc[cc],cnt)  ;
        for(int i=0;i<nn;++i) {
          loc.second = i ;
          data.second = cnt ;
          datamap.push_back(data) ;
          locations.push_back(loc) ;
          data.first++ ;
          cnt++ ;
          nodeTag[cc][i] = 5 ;
        }
      }else{
        std::vector<char>(1).swap(nodeTag[cc]);//to avoid default allocated size for vector
        nodeTag[cc].clear();
      }
    } ENDFORALL; 
    start += nenodes ;
    nodeloc.allocate(EMPTY) ;
    
    // Now process cell created nodes
    int ncnodes =
      getNodeOffsets(nodeloc,num_inner_nodes,local_cells,dist,comm) ;
    
    
    FORALL(local_cells, cc){
      if(num_inner_nodes[cc] != 0){
        std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
        int nn = num_inner_nodes[cc] ;

        std::vector<char>(nn).swap(nodeTag[cc]);
        pair<int,int> loc(cc,0) ;
        pair<int,int> data(start+nodeloc[cc],cnt)  ;
        for(int i=0;i<nn;++i) {
          loc.second = i ;
          data.second = cnt ;
          datamap.push_back(data) ;
          locations.push_back(loc) ;
          data.first++ ;
          cnt++ ;
          nodeTag[cc][i] = 5 ;
        }
      }else{
        std::vector<char>(1).swap(nodeTag[cc]);
        nodeTag[cc].clear();
      }
    } ENDFORALL;
    
    start += ncnodes ;
    
    nodeloc.allocate(EMPTY) ;
    // now process face allocated nodes    
    int nfnodes =
      getNodeOffsets(nodeloc,num_inner_nodes,local_faces,dist,comm) ;
    

    FORALL(local_faces, cc){
      if(num_inner_nodes[cc] != 0){
        std::vector<char>(int(num_inner_nodes[cc])).swap(nodeTag[cc]);
        int nn = num_inner_nodes[cc] ;

        std::vector<char>(nn).swap(nodeTag[cc]);
        pair<int,int> loc(cc,0) ;
        pair<int,int> data(start+nodeloc[cc],cnt) ;
        for(int i=0;i<nn;++i) {
          loc.second = i ;
          data.second = cnt ;
          datamap.push_back(data) ;
          locations.push_back(loc) ;
          data.first++ ;
          cnt++ ;
          nodeTag[cc][i] = 5 ;
        }
      }else{
        std::vector<char>(1).swap(nodeTag[cc]);
        nodeTag[cc].clear();
      }
    } ENDFORALL;
    
    start += nfnodes ;

    // Sanity Check on file size
    int local_size = filedata.size() ;
    int global_size = 0 ;
    MPI_Allreduce(&local_size,&global_size,1,MPI_INT,MPI_SUM,comm) ;
    
    // Check to see if the tag file contained enough tags
    if(start != global_size) {
      if(Loci::MPI_rank == 0) {
        cerr << "Refinement tag data size mismatch! Looking for " << start << " tags, but found " << filedata.size() << "!" << endl ;
      }
    }

    // Communicate tag file information as requested in datamap
    vector<pair<int,char> > returndata ;
    communicateFileData(filedata,datamap,returndata,comm) ;

    // Store node tags in the nodeTag datastructure
    for(size_t i=0;i<returndata.size();++i) {
      int locid = returndata[i].first ;
      int loc = locations[locid].first ;
      int off = locations[locid].second ;
      nodeTag[loc][off] = returndata[i].second ;
    }
#ifdef DEBUG
    entitySet dom2 = entitySet(seq) ;
    
    FORALL(dom2,ii) {
      int sz = nodeTag[ii].size() ;
      for(int i=0;i<sz;++i)
        if(nodeTag[ii][i] < 0 || nodeTag[ii][i] > 2) {
          cerr <<  "invalid nodeTag " << i << ' ' << int(nodeTag[ii][i]) << endl ;
        }
    } ENDFORALL ;
#endif
    
  }
     
  
  
};

register_rule<get_nodetag> register_get_nodetag;



class get_fineCellTag : public pointwise_rule{
  const_blackbox<vector<char> > inputTagsData;
  const_store<int> num_fine_cells;
  store<std::vector<char> > fineCellTag;
public:
  get_fineCellTag(){
    name_store("inputCellTagsData", inputTagsData);
    name_store("fineCellTag", fineCellTag);
    name_store("num_fine_cells", num_fine_cells);
    input("inputCellTagsData");
    input("num_fine_cells");
    output("fineCellTag");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    vector<char> filedata  = *inputTagsData ;

    MPI_Comm comm = MPI_COMM_WORLD ;    
    
    // Now determine what data we need to collect from filedata
    fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
    entitySet my_entities ;
    my_entities = ~EMPTY ;

    if(dist != 0)
      my_entities = (dist->my_entities) ;

  

    Loci::constraint geom_cells;
    geom_cells = (Loci::exec_current_fact_db)->get_variable("geom_cells");

    entitySet  local_cells;
    //don't know if it's necessray
    local_cells = (my_entities)&(*geom_cells);
    
    // Gather data maps for nodes created at edges, then cells, then faces
    // Note, since each of these are given a file number starting from zero
    // and we are using local2file to get a consistent numbering we need
    // to process each of these steps one by one.
    vector<pair<int,int> > datamap ;
    vector<pair<int,int> > locations ; ;
    store<int> nodeloc ; // location of nodes

   
    int start = 0;
    int cnt = 0 ;
    // Now process cell created nodes
    int ncnodes =
      getNodeOffsets(nodeloc,num_fine_cells,local_cells,dist,comm) ;
    
    
    FORALL(local_cells, cc){
      if(num_fine_cells[cc] != 0){
        std::vector<char>(int(num_fine_cells[cc])).swap(fineCellTag[cc]);
        int nn = num_fine_cells[cc] ;

        std::vector<char>(nn).swap(fineCellTag[cc]);
        pair<int,int> loc(cc,0) ;
        pair<int,int> data(start+nodeloc[cc],cnt)  ;
        for(int i=0;i<nn;++i) {
          loc.second = i ;
          data.second = cnt ;
          datamap.push_back(data) ;
          locations.push_back(loc) ;
          data.first++ ;
          cnt++ ;
          fineCellTag[cc][i] = 5 ;
        }
      }else{
        std::cerr <<"ERROR:: num_fine_cells is 0"<< std::endl;
        std::vector<char>(1).swap(fineCellTag[cc]);
        fineCellTag[cc].clear();
      }
    } ENDFORALL;
    
    start += ncnodes ;
    
    nodeloc.allocate(EMPTY) ;
   

    // Sanity Check on file size
    int local_size = filedata.size() ;
    int global_size = 0 ;
    MPI_Allreduce(&local_size,&global_size,1,MPI_INT,MPI_SUM,comm) ;
    
    // Check to see if the tag file contained enough tags
    if(start != global_size) {
      if(Loci::MPI_rank == 0) {
        cerr << "Refinement tag data size mismatch! Looking for " << start << " tags, but found " << filedata.size() << "!" << endl ;
      }
    }

    // Communicate tag file information as requested in datamap
    vector<pair<int,char> > returndata ;
    communicateFileData(filedata,datamap,returndata,comm) ;

    // Store node tags in the fineCellTag datastructure
    for(size_t i=0;i<returndata.size();++i) {
      int locid = returndata[i].first ;
      int loc = locations[locid].first ;
      int off = locations[locid].second ;
      fineCellTag[loc][off] = returndata[i].second ;
    }
#ifdef DEBUG
    entitySet dom2 = entitySet(seq) ;
    
    FORALL(dom2,ii) {
      int sz = fineCellTag[ii].size() ;
      for(int i=0;i<sz;++i)
        if(fineCellTag[ii][i] < 0 || fineCellTag[ii][i] > 2) {
          cerr <<  "invalid fineCellTag " << i << ' ' << int(fineCellTag[ii][i]) << endl ;
        }
    } ENDFORALL ;
#endif
    
  }
     
  
  
};

register_rule<get_fineCellTag> register_get_fineCellTag;
                                                

class input_CellTagsData : public singleton_rule {
  const_param<std::string> tagfile_par ;
  blackbox<vector<char> > inputTagsData ;
public:
  input_CellTagsData() {
    name_store("cell_tagfile_par",tagfile_par) ;
    name_store("inputCellTagsData",inputTagsData) ;
    input("cell_tagfile_par") ;
    output("inputCellTagsData") ;
    disable_threading() ;
  }
  virtual void compute(const sequence &seq) {
    string filename = *tagfile_par;
    MPI_Comm comm = MPI_COMM_WORLD ;    

    size_t p = filename.find("ref_sca");
    if(p != string::npos){//read hdf5 file
      *inputTagsData = read_tags_hdf5(filename,string("ref"),comm) ;
    } else {
      *inputTagsData = read_tags_txt(filename,comm) ;
    }
  }
} ;

register_rule<input_CellTagsData> register_input_CellTagsData ;

class get_celltag : public pointwise_rule{
  const_blackbox<vector<char> > inputTagsData ;
  store<char> cellTag;
  
public:
  get_celltag(){
    name_store("inputCellTagsData", inputTagsData);
    name_store("cellTag", cellTag);
    input("inputCellTagsData");
    output("cellTag");
    constraint("geom_cells");
    disable_threading();
  }
  virtual void compute(const sequence &seq){
    MPI_Comm comm = MPI_COMM_WORLD ;    
    
    vector<char> filedata  = *inputTagsData ;
    
    // Now determine what data we need to collect from filedata
    fact_db::distribute_infoP dist = (Loci::exec_current_fact_db)->get_distribute_info() ;
    entitySet dom = entitySet(seq);
    if(dist!=0)dom = dom & dist->my_entities;
    
    int offset= getFileNumberOffset(dom,comm);
   
    vector<pair<int,int> > datamap ;
    add_to_datamap(dom,datamap,offset,dist) ;
    
    
    vector<pair<int,char> > returndata ;
    communicateFileData(filedata,datamap,returndata,comm) ;
    
    
    for(size_t i=0;i<returndata.size();++i) {
      cellTag[returndata[i].first] = returndata[i].second ;
    }
        
  }

};
register_rule<get_celltag> register_get_celltag;
