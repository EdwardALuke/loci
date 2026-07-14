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
#include <Loci.h> 
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::ofstream ;
using std::ios ;
using std::sort ;
using std::unique ;
#include "extract.h"

// Extract a surface and write it out as an ascii file

struct surf_info {
  string name ;
  vector<pair<int,int> > edge_list ;
} ;
  
void get_surf(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries) {
  FATAL(Loci::MPI_processes != 1) ;
  store<vector3d<double> > pos ;
  string posname = getPosFile(output_dir,iteration,casename) ;
  hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT, MPI_COMM_WORLD) ;
  if(file_id < 0) {
    cerr << "unable to get grid positions for iteration " << iteration
         << endl ;
    cerr << "does file '" << posname << "' exist?" << endl ;
    Loci::Abort() ;
    exit(-1) ;
  }

  fact_db facts ;
  Loci::readContainer(file_id,"pos",pos.Rep(),EMPTY,facts) ;
  Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

  string gridtopo = getTopoFileName(output_dir, casename, iteration) ;

  file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  
#ifdef H5_USE_16_API
  hid_t bndg = H5Gopen(file_id,"boundaries") ;
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT) ;
#endif
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  vector<string> processed_bcs ;

  vector<int> f2size ;
  vector<int> fsurfid ;
  vector<int> f2node ;

  vector<surf_info> bc_data ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;
#ifdef H5_USE_16_API
    hid_t bcg = H5Gopen(bndg,buf) ;
#else
    hid_t bcg = H5Gopen(bndg,buf,H5P_DEFAULT) ;
#endif
    
    int nquads = sizeElementType(bcg,"quads") ;
    int ntrias = sizeElementType(bcg,"triangles") ;
    int ngeneral = sizeElementType(bcg,"nside_sizes") ;
    
    vector<Array<int,3> > trias(ntrias) ;
    readElementType(bcg,"triangles",trias) ;
    vector<Array<int,4> > quads(nquads) ;
    readElementType(bcg,"quads",quads) ;
    vector<int> nside_sizes(ngeneral) ;
    readElementType(bcg,"nside_sizes",nside_sizes) ;
    int nside_nodes_size = sizeElementType(bcg,"nside_nodes") ;
    vector<int> nside_nodes(nside_nodes_size) ;
    readElementType(bcg,"nside_nodes",nside_nodes) ;

    int bc_to_extract = false;
    for(size_t k=0;k<boundaries.size();++k)
      if(boundaries[k] == buf)
        bc_to_extract = true ;
    if(!bc_to_extract) { // collect boundary edges
      vector<pair<int,int> > edges ;
      for(int i=0;i<ntrias;++i) {
        pair<int,int> e ;
        e.first = min(trias[i][0],trias[i][1]) ;
        e.second = max(trias[i][0],trias[i][1]) ;
        edges.push_back(e) ;
        e.first = min(trias[i][2],trias[i][1]) ;
        e.second = max(trias[i][2],trias[i][1]) ;
        edges.push_back(e) ;
        e.first = min(trias[i][2],trias[i][0]) ;
        e.second = max(trias[i][2],trias[i][0]) ;
        edges.push_back(e) ;
      }
      for(int i=0;i<nquads;++i) {
        pair<int,int> e ;
        e.first = min(quads[i][0],quads[i][1]) ;
        e.second = max(quads[i][0],quads[i][1]) ;
        edges.push_back(e) ;
        e.first = min(quads[i][2],quads[i][1]) ;
        e.second = max(quads[i][2],quads[i][1]) ;
        edges.push_back(e) ;
        e.first = min(quads[i][2],quads[i][3]) ;
        e.second = max(quads[i][2],quads[i][3]) ;
        edges.push_back(e) ;
        e.first = min(quads[i][0],quads[i][3]) ;
        e.second = max(quads[i][0],quads[i][3]) ;
        edges.push_back(e) ;
      }
      int off = 0 ;
      for(int i=0;i<ngeneral;++i) {
        int sz = nside_sizes[i] ;
        int n1 = nside_nodes[off] ;
        int n2 = nside_nodes[off+sz-1] ;
        pair<int,int> e ;
        e.first = min(n1,n2) ;
        e.second = max(n1,n2) ;
        edges.push_back(e) ;
        for(int j=1;j<sz;++j) {
          int n1 = nside_nodes[off+j-1] ;
          int n2 = nside_nodes[off+j] ;
          e.first = min(n1,n2) ;
          e.second = max(n1,n2) ;
          edges.push_back(e) ;
        }
        off += sz ;
      }
      sort(edges.begin(),edges.end()) ;
      bc_data.push_back(surf_info()) ;
      bc_data.back().name = string(buf) ;
      int esz = edges.size() ;
      for(int i=0;i<esz;++i) {
        if((i+1 == esz) || (edges[i] != edges[i+1]))
          bc_data.back().edge_list.push_back(edges[i]) ;
        else 
          i++ ;
      }
    } else {
      processed_bcs.push_back(string(buf)) ;
      cout << "processing bc: " << buf << endl ;

      for(int i=0;i<ntrias;++i) {
        f2size.push_back(3) ;
        fsurfid.push_back(bc) ;
        f2node.push_back(trias[i][0]) ;
        f2node.push_back(trias[i][1]) ;
        f2node.push_back(trias[i][2]) ;
      }
      for(int i=0;i<nquads;++i) {
        f2size.push_back(4) ;
        fsurfid.push_back(bc) ;
        f2node.push_back(quads[i][0]) ;
        f2node.push_back(quads[i][1]) ;
        f2node.push_back(quads[i][2]) ;
        f2node.push_back(quads[i][3]) ;
      }
      int off = 0 ;
      for(int i=0;i<ngeneral;++i) {
        int sz = nside_sizes[i] ;
        f2size.push_back(sz) ;
        fsurfid.push_back(bc) ;
        for(int j=0;j<sz;++j)
          f2node.push_back(nside_nodes[off+j]) ;
        off += sz ;
      }
    }
  }
  if(processed_bcs.size() != boundaries.size()) {
    cerr << "Warning: Not all boundaries were found in the input file!"
         << endl
         << "         Check -bc flags!"
         << endl ;
  }
  vector<int> node_set = f2node ;

  if(node_set.size() == 0) {
    cerr << "nothing to process in surf extract!  Exiting."
         << endl ;
    Loci::Abort() ;
    exit(-1) ;
  }
  sort(node_set.begin(),node_set.end()) ;
  node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;

  int mx = node_set[0] ;
  int mn = node_set[0] ;
  for(size_t i=0;i<node_set.size();++i) {
    mx = max(mx,node_set[i]) ;
    mn = min(mn,node_set[i]) ;
  }
  
  entitySet nset = interval(mn,mx);

  store<int> nmap ;
  
  nmap.allocate(nset) ;
  for(int i=mn;i<=mx;++i)
    nmap[i] = -1 ;
  
  for(size_t i=0;i<node_set.size();++i) {
    nmap[node_set[i]] = i+1 ;
  }

  int nbc = 0 ;
  for(size_t i=0;i<bc_data.size();++i) {
    vector<pair<int,int> > elist ;
    for(size_t j=0;j<bc_data[i].edge_list.size();++j) {
      int n1 = bc_data[i].edge_list[j].first ;
      int n2 = bc_data[i].edge_list[j].second ;
      if((n1 >= mn && n1 <= mx && nmap[n1] >=0) &&
         (n2 >= mn && n2 <= mx && nmap[n2] >=0) ) {
        pair<int,int> e(nmap[n1],nmap[n2]) ;
        elist.push_back(e) ;
      }
    }
    bc_data[i].edge_list.swap(elist) ;
    if(bc_data[i].edge_list.size() > 0)
      nbc++ ;
  }
    
  for(size_t i=0;i<f2node.size();++i) {
    if(nmap[f2node[i]] < 0)
      cerr << "nmap invalid for f2node, i = " << i << endl ;
    f2node[i] = nmap[f2node[i]] ;
  }


  string filename = processed_bcs[0] + ".gsurf" ;
  ofstream sfile(filename.c_str(),ios::out) ;
  sfile.precision(15) ;
  sfile << "# number of nodes" << endl ;
  sfile << node_set.size() << endl ;
  sfile << "# Node positions" << endl ;
  int minpos = pos.domain().Min() ;
  for(size_t i=0;i<node_set.size();++i) {
    int nd = node_set[i]-1+minpos ;
    sfile << pos[nd].x << ' ' << pos[nd].y << ' ' << pos[nd].z << endl ;
  }
  sfile << "# number of faces" << endl ;
  sfile << f2size.size() << endl ;
  sfile << "# Surface face Definition (#pnts nodes)" << endl ;
  int off = 0 ;
  for(size_t i=0;i<f2size.size();++i) {
    int sz = f2size[i] ;
    sfile << sz  ;
    for(int j=0;j<sz;++j)
      sfile << ' ' << f2node[off+j] ;
    sfile << endl ;
    off += sz ;
  }

  sfile << "# number of boundary edge tags" << endl;
  sfile << nbc << endl ;
  for(size_t i=0;i<bc_data.size();++i) {
    if(bc_data[i].edge_list.size() != 0) {
      sfile << bc_data[i].name << endl ;
      sfile << bc_data[i].edge_list.size() << endl ;
      for(size_t j=0;j<bc_data[i].edge_list.size();++j) {
        sfile << bc_data[i].edge_list[j].first << ' ' 
              <<  bc_data[i].edge_list[j].second << endl ;
      }
    }
  }
  
  sfile.close() ;

  // write out proper solid mesh file
  string solid_filename = casename + ".surf" ;
  if(processed_bcs.size() == 1)
    solid_filename = processed_bcs[0] + ".surf" ;

  ofstream ssfile(solid_filename.c_str(),ios::out) ;
  ssfile.precision(15) ;

  int ntri = 0 ;
  int nqua = 0 ;
  
  for(size_t i=0;i<f2size.size();++i) {
    int sz = f2size[i] ;
    if(sz == 3)
      ntri++ ;
    else if(sz == 4)
      nqua++ ;
    else
      ntri += sz-2 ; // Convert general faces to triangles
  }
  ssfile << ntri << ' ' << nqua << ' ' 
         << node_set.size() << endl ;
  
  minpos = pos.domain().Min() ;
  double normal_spacing = 0 ;
  for(size_t i=0;i<node_set.size();++i) {
    int nd = node_set[i]-1+minpos ;
    ssfile << pos[nd].x << ' ' << pos[nd].y << ' ' << pos[nd].z
           << ' '<< normal_spacing << endl ;
  }
  // output triangles
  off = 0 ;
  for(size_t i=0;i<f2size.size();++i) {
    int sz = f2size[i] ;
    int sid = fsurfid[i]+1 ;
    if(sz == 3) {
      ssfile << f2node[off+0] << ' '
             << f2node[off+1] << ' '
             << f2node[off+2] << ' ' 
             << sid << " 0 1" << endl ;
    } else if(sz > 4) { // convert general face to triangles
      for(int i=1;i<sz-1;++i) {
        ssfile << f2node[off+0] << ' '
               << f2node[off+i] << ' '
               << f2node[off+i+1] << ' ' 
               << sid << " 0 1" << endl ;
      }
    }
    off += sz ;
  }
  // output quads
  off = 0 ;
  for(size_t i=0;i<f2size.size();++i) {
    int sz = f2size[i] ;
    int sid = fsurfid[i]+1 ;
    if(sz == 4) {
      ssfile << f2node[off+0] << ' '
             << f2node[off+1] << ' '
             << f2node[off+2] << ' '
             << f2node[off+3] << ' '
             << sid << " 0 1" << endl ;
    }
    off += sz ;
  }
}
