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

void get_2dgv(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries,
              int view) {
  FATAL(Loci::MPI_processes != 1) ;
  store<vector3d<double> > pos ;
  string posname = getPosFile(output_dir,iteration,casename) ;
  hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to get grid positions for iteration " << iteration
         << endl ;
    cerr << "does file '" << posname << "' exist?" << endl ;
    Loci::Abort() ;
    exit(-1) ;
  }

  fact_db facts ;
  readData(file_id,"pos",pos.Rep(),EMPTY,facts) ;
  Loci::hdf5CloseFile(file_id) ;

  int npnts = pos.domain().size() ;

  
  string iblankname = "output/grid_iblank." + iteration + "_" + casename ;
  store<unsigned char> iblank ;
  entitySet pdom = interval(1,npnts) ;
  iblank.allocate(pdom) ;
  struct stat tmpstat ;
  if(stat(iblankname.c_str(),&tmpstat)== 0) {
    hid_t file_id = Loci::hdf5OpenFile(iblankname.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to get iblank info for iteration " << iteration
           << endl ;
      cerr << "is file '" << iblankname << "' corrupted?" << endl ;
      Loci::Abort() ;
      exit(-1) ;
    }

    store<unsigned char> iblank_tmp ;
    readData(file_id,"iblank",iblank_tmp.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;
    entitySet dom = iblank_tmp.domain() ;
    int cnt = 1 ;
    FORALL(dom,nd) {
      iblank[cnt++] = iblank_tmp[nd] ;
    } ENDFORALL ;
  } else {
    for(int i=1;i<=npnts;++i)
      iblank[i] = 0 ;
  }


  vector<Array<int,3> > edges ;

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

  int nelem = 1 ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;

    int bc_to_extract = false;
    for(size_t k=0;k<boundaries.size();++k)
      if(boundaries[k] == buf)
        bc_to_extract = true ;
    if(!bc_to_extract)
      continue ;
    processed_bcs.push_back(string(buf)) ;
    cout << "processing bc: " << buf << endl ;
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

    for(int i=0;i<ntrias;++i) {
      if(iblank[trias[i][0]] < 2 || iblank[trias[i][1]] < 2 ||
         iblank[trias[i][2]] < 2) {
        Array<int,3> e ;
        e[0] = min(trias[i][0],trias[i][1]) ;
        e[1] = max(trias[i][0],trias[i][1]) ;
        e[2] = nelem ;
        edges.push_back(e) ;
        e[0] = min(trias[i][2],trias[i][1]) ;
        e[1] = max(trias[i][2],trias[i][1]) ;
        edges.push_back(e) ;
        e[0] = min(trias[i][2],trias[i][0]) ;
        e[1] = max(trias[i][2],trias[i][0]) ;
        edges.push_back(e) ;
        nelem++ ;
      }
    }
    for(int i=0;i<nquads;++i) {
      if(iblank[quads[i][0]] < 2 || iblank[quads[i][1]] < 2 ||
         iblank[quads[i][2]] < 2 || iblank[quads[i][3]] < 2) {
        Array<int,3> e ;
        e[0] = min(quads[i][0],quads[i][1]) ;
        e[1] = max(quads[i][0],quads[i][1]) ;
        e[2] = nelem ;
        edges.push_back(e) ;
        e[0] = min(quads[i][2],quads[i][1]) ;
        e[1] = max(quads[i][2],quads[i][1]) ;
        edges.push_back(e) ;
        e[0] = min(quads[i][2],quads[i][3]) ;
        e[1] = max(quads[i][2],quads[i][3]) ;
        edges.push_back(e) ;
        e[0] = min(quads[i][0],quads[i][3]) ;
        e[1] = max(quads[i][0],quads[i][3]) ;
        edges.push_back(e) ;
        nelem++ ;
      }
    }
    int off = 0 ;
    for(int i=0;i<ngeneral;++i) {
      int sz = nside_sizes[i] ;
      int n1 = nside_nodes[off] ;
      int n2 = nside_nodes[off+sz-1] ;
      bool blanked = true ;
      for(int j=0;j<sz;++j) 
        if(iblank[nside_nodes[off+j]] < 2)
          blanked = false ;
      if(!blanked) {
        Array<int,3> e ;
        e[0] = min(n1,n2) ;
        e[1] = max(n1,n2) ;
        e[2] = nelem ;
        edges.push_back(e) ;
        for(int j=1;j<sz;++j) {
          int n1 = nside_nodes[off+j-1] ;
          int n2 = nside_nodes[off+j] ;
          e[0] = min(n1,n2) ;
          e[1] = max(n1,n2) ;
          e[2] = nelem ;
          edges.push_back(e) ;
        }
        nelem++ ;
      }
      off += sz ;
    }
  }
  sort(edges.begin(),edges.end()) ;
  vector<int> node_set(edges.size()*2) ;
  for(size_t i=0;i<edges.size();++i) {
    node_set[i*2] = edges[i][0] ;
    node_set[i*2+1] = edges[i][1] ;
  }
  if(processed_bcs.size() != boundaries.size()) {
    cerr << "Warning: Not all boundaries were found in the input file!"
         << endl
         << "         Check -bc flags!"
         << endl ;
  }
  if(node_set.size() == 0) {
    cerr << "nothing to process in 2dgv extract!  Exiting."
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
  int minpos = pos.domain().Min() ;
  for(size_t i=0;i<node_set.size();++i) {
    nmap[node_set[i]] = i+1 ;
  }

  string vname = variables[0] ;
  string outname = vname + '.' + iteration ;
  std::ofstream out(outname.c_str(),ios::out) ;
  out.precision(16) ;
  out << "general " << endl ;
  out << node_set.size() << " 1" << endl ;
  for(size_t i=0;i<node_set.size();++i) {
    int nd = node_set[i]-1+minpos ;
    double x = pos[nd].x ;
    double y = pos[nd].y ;
    double z = pos[nd].z ;

    switch(view) {
    case VIEWXR:
      out << x << ' ' << sqrt(y*y+z*z) << endl ;
      break ;
    case VIEWXY:
      out << x << ' ' << y << endl ;
      break ;
    case VIEWXZ:
      out << x << ' ' << z << endl ;
      break ;
    case VIEWYZ:
      out << y << ' ' << z << endl ;
      break ;
    default:
      cerr << "wrong output coordinate system!" << endl ;
      Loci::Abort() ;
      exit(EXIT_FAILURE) ;
    }
  }

  vector<Array<int,4> > edge_f ;
  for(size_t i=0;i<edges.size();) {
    Array<int,4> e ;
    e[0] = nmap[edges[i][0]] ;
    e[1] = nmap[edges[i][1]] ;
    e[2] = edges[i][2] ;
    e[3] = -1 ;
    ++i ;
    if(i<edges.size()             &&
       edges[i][0]==edges[i-1][0] &&
       edges[i][1]==edges[i-1][1] ) {
      e[3] = edges[i][2] ;
      ++i ;
    }
    edge_f.push_back(e) ;
  }    

  out << edge_f.size() << " 1 " << nelem-1 << " 1" << endl ;
  for(size_t i=0;i<edge_f.size();++i)
    out << edge_f[i][0] << ' '
        << edge_f[i][1] << ' '
        << edge_f[i][2] << ' '
        << edge_f[i][3] << endl ;

  const string var_name(variables[0]) ;
  const int vt = variable_types[0] ;
  const string filename(variable_filenames[0]) ;
  switch(vt) {
  case NODAL_SCALAR:
    {
      hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                         H5F_ACC_RDONLY,
                                         H5P_DEFAULT) ;
      if(file_id < 0) {
        cerr << "unable to open file '" << filename << "'!" << endl ;
        return ;
      }

      fact_db facts ;
      store<float> scalar ;
      readData(file_id,var_name,scalar.Rep(),EMPTY,facts) ;
      entitySet dom = scalar.domain() ;

      int min_val= dom.Min() ;
      for(size_t i=0;i<node_set.size();++i) {
        int nd = node_set[i]-1+min_val ;
        out << scalar[nd] << endl ;
      }
      Loci::hdf5CloseFile(file_id) ;
    }
    break;
  case NODAL_DERIVED:
    {
      vector<float> dval(npnts) ;
      getDerivedVar(dval,var_name,casename,iteration) ;

      for(size_t i=0;i<node_set.size();++i) {
        int nd = node_set[i]-1 ;
        out << dval[nd] << endl ;
      }
    }
    break;
  case NODAL_VECTOR:
    {
      hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                         H5F_ACC_RDONLY,
                                         H5P_DEFAULT) ;
      if(file_id < 0) {
        cerr << "unable to open file '" << filename << "'!" << endl ;
        Loci::Abort() ;
        exit(-1) ;
      }
      
      fact_db facts ;
      store<vector3d<float> > vec ;
      readData(file_id,var_name,vec.Rep(),EMPTY,facts) ;
      entitySet dom = vec.domain() ;

      int min_val= dom.Min() ;
      for(size_t i=0;i<node_set.size();++i) {
        int nd = node_set[i]-1+min_val ;
        out << norm(vec[nd]) << endl ;
      }
      Loci::hdf5CloseFile(file_id) ;
    }
    break;
  case NODAL_MASSFRACTION:
    {
      string filename = "output/mix." + iteration + "_" + casename ;
    
      hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                         H5F_ACC_RDONLY,
                                         H5P_DEFAULT) ;
      if(file_id < 0) {
        cerr << "unable to open file '" << filename << "'!" << endl ;
        Loci::Abort() ;
        exit(-1) ;
      }

      fact_db facts ;
      storeVec<float> mix ;
      readData(file_id,"mixture",mix.Rep(),EMPTY,facts) ;
      param<string> species_names ;
      readData(file_id,"species_names",species_names.Rep(),EMPTY,facts) ;
      Loci::hdf5CloseFile(file_id) ;
      
      map<string,int> smap ;
      std::istringstream iss(*species_names) ;
      for(int i=0;i<mix.vecSize();++i) {
        string s ;
        iss >> s ;
        smap[s] = i ;
      }
      
      entitySet dom = mix.domain() ;

      vector<float> vec(npnts) ;
    
      string sp = string(var_name.c_str()+1) ;
      map<string,int>::const_iterator mi = smap.find(sp) ;
      if(mi == smap.end()) {
        cerr << "warning, species " << sp << " does not exist in dataset!"
             << endl ;
      } else {
        const int ind = mi->second ;
        int c = 0 ;
        FORALL(dom,nd) {
          vec[c++] = mix[nd][ind] ;
        } ENDFORALL ;
        
        for(size_t i=0;i<node_set.size();++i) {
          int nd = node_set[i]-1 ;
          out << vec[nd] << endl ;
        }
      }
    }
    break ;
  default:
    cerr << "unable to process variable " << var_name << endl ;
    cerr << "this variable type can not be plotted by 2dgv" << endl ;
    break ;
  }
}
