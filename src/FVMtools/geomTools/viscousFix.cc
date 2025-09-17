//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
using std::string ;
using std::vector ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::ofstream ;
using std::ifstream ;
using std::ios ;
using std::sort ;
using std::unique ;

void dump_scalar(store<float> &c2n, string sname, string iter, string casename) {
  string filename = "output/" + sname + "_sca." + iter + "_" + casename ;
  
  if(Loci::MPI_rank == 0)
    cout << "writing file " << filename << endl ;
  
  fact_db facts ;

  hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT) ;
  
  Loci::writeContainer(file_id,sname,c2n.Rep(),facts) ;

  Loci::hdf5CloseFile(file_id) ;
}

template<class T> void readElementType(hid_t group_id, const char *element_name,
                                       vector<T> &v) {
  if(v.size() > 0) {
#ifdef H5_USE_16_API
    hid_t dataset = H5Dopen(group_id,element_name) ;
#else
    hid_t dataset = H5Dopen(group_id,element_name,H5P_DEFAULT) ;
#endif

    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    
    H5Dread(dataset,dp->get_hdf5_type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,&v[0]) ;
    H5Dclose(dataset) ;
  }
}

int  sizeElementType(hid_t group_id, const char *element_name) {
#ifdef H5_USE_16_API
  hid_t dataset = H5Dopen(group_id,element_name) ;
#else
  hid_t dataset = H5Dopen(group_id,element_name,H5P_DEFAULT) ;
#endif
  if(dataset < 0) {
#ifdef H5_USE_16_API
    H5Eclear() ;
#else
    H5Eclear(H5E_DEFAULT) ;
#endif
    return 0 ;
  }
  hid_t dspace = H5Dget_space(dataset) ;

  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;
  
  
  H5Dclose(dataset) ;
  return int(size) ;
  
}

void Usage() {
  cerr << "Usage:" << endl
       << "  viscousFix <arguments> <casename> <iteration>" << endl
       << "  where arguments are: " << endl
       << "  -bc <boundary> : identify viscous boundaries" << endl
       << "  -delta <distance>: control refinement off boundary to delta" << endl
       << "  -tag <file>: input tag filename" << endl
       << "  -o <file>: output tag filename" << endl ;
  exit(-1) ;
}

string getPosFile(string output_dir,string iteration, string casename) {
  string posname = output_dir+"/grid_pos." + iteration + "_" + casename ;
  struct stat tmpstat ;
  if(stat(posname.c_str(),&tmpstat) != 0) {
    posname = output_dir+"/grid_pos." + casename ;
  } else if(tmpstat.st_size == 0) {
    posname = output_dir+"/grid_pos." + casename ;
  }
  return posname ;
}

int main(int ac, char *av[]) {
  vector<string> bclist ;
  double delta = 0 ;
  string casename ;
  string iteration ;
  string input_file = "refine.dat" ;
  string output_file = "crefine.dat" ;
    
  bool found_casename =false ;
  bool found_iteration = false ;

  Loci::Init(&ac,&av) ;
  for(int i=1;i<ac;++i) {
    if(av[i][0] == '-') {
      if(!strcmp(av[i],"-bc")) {
        i++ ;
        string v(av[i]) ;
        if(av[i][0] >= '0' && av[i][0] <= '9')
          v = "BC_"+v ;
        bclist.push_back(v) ;
      } else if(!strcmp(av[i],"-delta")) {
        i++ ;
        delta = atof(av[i]) ;
      } else if(!strcmp(av[i],"-tag")) {
        i++ ;
        input_file = av[i] ;
      } else if(!strcmp(av[i],"-o")) {
        i++ ;
        output_file = av[i] ;
      } else {
        cerr << "bad option " << av[i] << endl ;
        Usage() ;
      }
    } else {
      if(!found_casename) {
        casename = av[i] ;
        found_casename = true ;
      } else if(!found_iteration) {
        iteration = av[i] ;
        found_iteration = true ;
      } else {
        cerr << "extra argument '" << av[i] << "' is not understood" << endl ;
        Usage() ;
      }
    }
  }
  if(!found_iteration) {
    Usage() ;
  }


  string gridtopo = "output/" + casename +".topo" ;


#ifdef H5_USE_16_API
  H5Eset_auto(NULL,NULL) ;
#else
  H5Eset_auto(H5E_DEFAULT,NULL,NULL) ;
#endif
  hid_t file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to open file '" << gridtopo << "'" << endl ;
    exit(-1) ;
  }
                   

#ifdef H5_USE_16_API  
  hid_t bndg = H5Gopen(file_id,"boundaries") ;
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT) ;
#endif
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  vector<string> processed_bcs ;

  vector<int> node_set ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;

    int bc_to_extract = false;
    for(size_t k=0;k<bclist.size();++k)
      if(bclist[k] == buf)
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

    if(ntrias > 0) {    
      vector<Array<int,3> > trias(ntrias) ;
      readElementType(bcg,"triangles",trias) ;
      for(int i=0;i<ntrias;++i) {
        node_set.push_back(trias[i][0]) ;
        node_set.push_back(trias[i][1]) ;
        node_set.push_back(trias[i][2]) ;
      }
      sort(node_set.begin(),node_set.end()) ;
      node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;
    }
    if(nquads > 0) {
      vector<Array<int,4> > quads(nquads) ;
      readElementType(bcg,"quads",quads) ;
      for(int i=0;i<nquads;++i) {
        node_set.push_back(quads[i][0]) ;
        node_set.push_back(quads[i][1]) ;
        node_set.push_back(quads[i][2]) ;
        node_set.push_back(quads[i][3]) ;
      }
      sort(node_set.begin(),node_set.end()) ;
      node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;
    }
    if(ngeneral > 0) {
      vector<int> nside_sizes(ngeneral) ;
      readElementType(bcg,"nside_sizes",nside_sizes) ;
      int nside_nodes_size = sizeElementType(bcg,"nside_nodes") ;
      vector<int> nside_nodes(nside_nodes_size) ;
      readElementType(bcg,"nside_nodes",nside_nodes) ;
      sort(nside_nodes.begin(),nside_nodes.end()) ;
      node_set.push_back(nside_nodes[0]) ;
      for(int i=1;i<nside_nodes_size;++i)
        if(nside_nodes[i-1] != nside_nodes[i])
          node_set.push_back(nside_nodes[i]) ;
    }
  }
  H5Fclose(file_id) ;
  sort(node_set.begin(),node_set.end()) ;
  node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;

  for(size_t i=0;i<node_set.size();++i)
    node_set[i] -= 1 ;
  
  if(node_set.size() == 0) {
    cerr << "no boundary nodes found!" << endl ;
    exit(-1) ;
  }
  cout << "processing "<< node_set.size() << " boundary nodes. " << endl ;

  store<vector3d<double> > pos ;
  string posname = getPosFile("output",iteration, casename) ;

  file_id = Loci::hdf5OpenFile(posname.c_str(),
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
  Loci::readContainer(file_id,"pos",pos.Rep(),EMPTY,facts) ;
  Loci::hdf5CloseFile(file_id) ;
  
  cout << "processing " << pos.domain().size() << " grid nodes" << endl ;

  vector<Loci::kdTree::coord3d> pnts(node_set.size()) ;
  vector<int> pnt_id(node_set.size()) ;
  for(size_t i=0;i<node_set.size();++i) {
    pnt_id[i] = i ;
    pnts[i][0] = pos[node_set[i]].x ;
    pnts[i][1] = pos[node_set[i]].y ;
    pnts[i][2] = pos[node_set[i]].z ;
  }
  Loci::kdTree::kd_tree nns(pnts,pnt_id) ;
  double ds = delta*delta ;

  vector<vector<int> > pnt_set(node_set.size()) ;


  entitySet pdom = pos.domain() ;
  FORALL(pdom,nd) {
    double rmin = ds ;
    Loci::kdTree::coord3d pt ;
    pt[0]= pos[nd].x ;
    pt[1]= pos[nd].y ;
    pt[2]= pos[nd].z ;
    int p = nns.find_closest(pt,rmin) ;
    if(p >= 0)
      pnt_set[p].push_back(nd) ;
  } ENDFORALL ;

  ifstream infile(input_file.c_str(),ios::in) ;

  store<unsigned char> refine ;
  refine.allocate(pdom) ;
  FORALL(pdom,nd) {
    int v ;
    infile >> v ;
    refine[nd] = 0 ;
    if(v != 0)
      refine[nd] = 1 ;
  } ENDFORALL ;

  int newmarks = 0 ;
  for(size_t i=0;i<node_set.size();++i) {
    int tmp =0;
    for(size_t j=0;j<pnt_set[i].size();++j)
      tmp += refine[pnt_set[i][j]] ;
    if(tmp != 0)
      for(size_t j=0;j<pnt_set[i].size();++j) {
        if(refine[pnt_set[i][j]] == 0)
          newmarks++ ;
        refine[pnt_set[i][j]] = 1 ;
      }
  }
  cout << "marking " << newmarks << " additional nodes near boundary." << endl ;
  ofstream ofile(output_file.c_str(),ios::out) ;

  FORALL(pdom,nd) {
   ofile << int(refine[nd]) << endl ;
  } ENDFORALL ;

  store<float> adapt ;
  adapt.allocate(pdom) ;
  
  FORALL(pdom,nd) {
    adapt[nd] = int(refine[nd]) ;
  } ENDFORALL ;

  dump_scalar(adapt,"reffix",iteration,casename) ;
  
  Loci::Finalize() ;
}
