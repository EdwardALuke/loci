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
#include <stdio.h>
#include <strings.h>

#include <Loci>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

using std::map ;
using std::vector ;
using std::string ;
using std::cout ;
using std::cerr ;
using std::endl ;
using Loci::vector3d ;

struct file_info {
  unsigned long numNodes ;
  unsigned long numFaces ;
  unsigned long numCells ;
  vector<short> cluster_sizes ;
} ;

unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

namespace Loci {
  extern int getClusterNumFaces(unsigned char *cluster) ;
  extern int fillClusterFaceSizes(unsigned char *cluster, int *sizes) ;
  int fillFaceInfo(unsigned char *cluster, multiMap &face2node,
                   Map &cl, Map &cr, int face_base) ;
  vector<unsigned char>
  encode_face_cluster(const multiMap &face2node,
                      const Map &cl, const Map &cr,
                      entitySet fcluster,
                      entitySet nodeSet,
                      entitySet cellSet) ;
}


void Usage() {
  cout << "Utility for extracting surfaces from vog grids." << endl
       << endl ;
  cout << "Usage: " << endl
       << "  vog2surf <grid> [ -o <ofile>.surf ] [ -surface <ofile>.surface ]" << endl
       << "    where <grid> is the input file name sans the '.vog' postfix" <<endl
       << "    <ofile> is the output filename, can be .surf format or .surface format " << endl
       << "    .surf is the solidMesh format, .surface will add the node index map from surface mesh " << endl
       << "     to volume mesh at the end of surf file. " <<endl 
       <<  endl ;
  cout << "  Outputs: <grid>.surf/<ofile>.surf/<ofile>.surface  and <grid>.names" << endl ;
  exit(-1) ;
}

struct surface_info {
  string name ;
  int id ;
  vector<Loci::Array<int,3> > trias ;
  vector<Loci::Array<int,4> > quads ;
  vector<vector<int> > gen_faces ;
} ;


void readSurfaces(string filename,
		  vector<string> &surfnames,
		  vector<surface_info> &surf_list,
		  vector<vector3d<double> > &pos,
                  vector<int>& node_map) {

  surf_list.clear() ;
  pos.clear() ;

  map<int,int> surf_lookup ;

  // read in boundary names.
  vector<pair<int,string> > boundary_ids ;
  Loci::readBCfromVOG(filename,boundary_ids) ;
  map<int,string> surf_id ;
  for(size_t i=0;i<boundary_ids.size();++i)
    surf_id[boundary_ids[i].first] = boundary_ids[i].second ;

  if(surfnames.size() == 0) {
    for(size_t i=0;i<boundary_ids.size();++i)
      surfnames.push_back(boundary_ids[i].second) ;
  }
  
  
  hid_t input_fid ; 
  input_fid = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(input_fid <= 0) {
    cerr << "unable to open file '" << filename << "'"<< endl ;
    Usage() ;
  }

  hid_t face_g = H5Gopen(input_fid,"face_info",H5P_DEFAULT) ;
  // Read cluster sizes
  hid_t dataset = H5Dopen(face_g,"cluster_sizes",H5P_DEFAULT) ;
  
  hid_t dspace = H5Dget_space(dataset) ;
  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;
  vector<unsigned short> csizes(size) ;
  hsize_t dimension = size ;
  hsize_t stride = 1 ;
  hsize_t start = 0 ;
  hsize_t count = size ;
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
  int rank = 1 ;
  hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
  hid_t err = H5Dread(dataset,H5T_NATIVE_USHORT,memspace,dspace,H5P_DEFAULT,&csizes[0]) ;
  if(err < 0) {
    cerr << "unable to read cluster sizes from file '" <<
      filename << "'" << endl ;
    exit(-1) ;
  }
  H5Dclose(dataset) ;
  H5Sclose(memspace) ;

  // Read in clusters and transform
  dataset = H5Dopen(face_g,"cluster_info",H5P_DEFAULT) ;
  dspace = H5Dget_space(dataset) ;
  start = 0 ;
  for(size_t c=0;c<size;++c) { // Loop over clusters
    count = csizes[c] ;
    vector<unsigned char> cluster(count) ;
    H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
    dimension = count ;
    memspace = H5Screate_simple(rank,&dimension,NULL) ;
    err = H5Dread(dataset,H5T_NATIVE_UCHAR,memspace,dspace,H5P_DEFAULT,&cluster[0]) ;
    if(err < 0) {
      cerr << "unable to read cluster from file '" << filename << "'" << endl ;
    }
    start += count ;
    // Now scan cluster into local buffers

    int nfaces = Loci::getClusterNumFaces(&cluster[0]) ;
    Loci::entitySet fclust = Loci::interval(0,nfaces-1) ;
    Loci::store<int> fcnts ;
    fcnts.allocate(fclust) ;
    Loci::fillClusterFaceSizes(&cluster[0],&fcnts[0]) ;
    Loci::multiMap face2node ;
    Loci::Map cl,cr ;
    face2node.allocate(fcnts) ;
    cl.allocate(fclust) ;
    cr.allocate(fclust) ;
    Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0) ;
    // Now loop over faces in cluster to determine if any are boundary 
    // faces
    for(int f=0;f<nfaces;++f) {
      if(cr[f] < 0) { // boundary face 
	// Now check to see if we have encountered it before
	map<int,int>::const_iterator mi = surf_lookup.find(-cr[f]) ;
	int bc ;
	if(mi == surf_lookup.end() ) { // First time to see this boundary tag
	  // Get name of boundary
	  string bc_name ;
	  map<int,string>::const_iterator si = surf_id.find(-cr[f]) ;
	  if(si == surf_id.end()) {
	    char buf[128] ;
	    bzero(buf,128) ;
	    snprintf(buf,127,"BC_%d",-cr[f]) ;
	    bc_name = buf ;
	  } else {
	    bc_name = si->second ;
	  }
	  int bc_id = surf_list.size() ;
	  surf_list.push_back(surface_info()) ;
	  surf_list[bc_id].name = bc_name ;
          surf_list[bc_id].id = -cr[f] ;
	  surf_lookup[-cr[f]] = bc_id ;
	  bc = bc_id ;
	} else {
	  bc = mi->second ;
	}
	int fsz = face2node[f].size() ;
	if(fsz == 3) {
	  Loci::Array<int,3> tri ;
	  tri[0] = face2node[f][0] ;
	  tri[1] = face2node[f][1] ;
	  tri[2] = face2node[f][2] ;
	  surf_list[bc].trias.push_back(tri) ;
	} else if(fsz == 4) {
	  Loci::Array<int,4> qua ;
	  qua[0] = face2node[f][0] ;
	  qua[1] = face2node[f][1] ;
	  qua[2] = face2node[f][2] ;
	  qua[3] = face2node[f][3] ;
	  surf_list[bc].quads.push_back(qua) ;
	} else {
	  vector<int> tmp(fsz) ;
	  for(int i=0;i<fsz;++i)
	    tmp[i] = face2node[f][i] ;
	  surf_list[bc].gen_faces.push_back(tmp) ;
	}
      }
    }
  }

  // read in positions
  hid_t fi = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
  unsigned long numNodes = readAttributeLong(fi,"numNodes") ;
    
  H5Gclose(fi) ;

  count = numNodes ;

    hsize_t lstart = 0 ;
      
  // Read in pos data from file i
  vector<Loci::vector3d<double> > pos_dat(numNodes) ;
  hid_t node_g = H5Gopen(input_fid,"node_info",H5P_DEFAULT) ;
  dataset = H5Dopen(node_g,"positions",H5P_DEFAULT) ;
  dspace = H5Dget_space(dataset) ;
      
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
  rank = 1 ;
  dimension = count ;
  memspace = H5Screate_simple(rank,&dimension,NULL) ;
  typedef Loci::data_schema_traits<Loci::vector3d<double> > traits_type ;
  Loci::DatatypeP dp = traits_type::get_type() ;
  hid_t datatype = dp->get_hdf5_type() ;
  err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
		      &pos_dat[0]) ;
  if(err < 0) {
    cerr << "unable to read positions from '" << filename << "'" << endl ;
    exit(-1) ;
  }
  H5Sclose(dspace) ;
  H5Dclose(dataset) ;
  H5Gclose(node_g) ;

  // Now mark all nodes that the faces access

  vector<int> used(numNodes) ;
  for(size_t i=0;i<numNodes;++i)
    used[i] = 0 ;
  
  int ssz = surf_list.size() ;

  
  for(int i=0;i<ssz;++i) {
    for(size_t j=0;j<surf_list[i].trias.size();++j) {
      used[surf_list[i].trias[j][0]] = 1 ;
      used[surf_list[i].trias[j][1]] = 1 ;
      used[surf_list[i].trias[j][2]] = 1 ;
    }
    for(size_t j=0;j<surf_list[i].quads.size();++j) {
      used[surf_list[i].quads[j][0]] = 1 ;
      used[surf_list[i].quads[j][1]] = 1 ;
      used[surf_list[i].quads[j][2]] = 1 ;
      used[surf_list[i].quads[j][3]] = 1 ;
    }
    for(size_t j=0;j<surf_list[i].gen_faces.size();++j) {
      for(size_t k=0;k<surf_list[i].gen_faces[j].size();++k)
	used[surf_list[i].gen_faces[j][k]] = 1 ;
    }
  }

  // Count nodes in the surface mesh
  int nnodes = 0 ;
  for(size_t i=0;i<numNodes;++i)
    if(used[i] != 0) {
      used[i] += nnodes++ ;
    }


  vector<vector3d<double> > ptmp(nnodes) ;
   node_map.resize(nnodes);
  
  for(size_t i=0;i<numNodes;++i)
    {
      if(used[i] != 0){
      ptmp[used[i]-1] = pos_dat[i] ;
      node_map[used[i]-1] = i;
    }
    }
  pos.swap(ptmp) ;

  for(int i=0;i<ssz;++i) {
    for(size_t j=0;j<surf_list[i].trias.size();++j) {
      surf_list[i].trias[j][0] = used[surf_list[i].trias[j][0]] ;
      surf_list[i].trias[j][1] = used[surf_list[i].trias[j][1]] ;
      surf_list[i].trias[j][2] = used[surf_list[i].trias[j][2]] ;
    }
    for(size_t j=0;j<surf_list[i].quads.size();++j) {
      surf_list[i].quads[j][0] = used[surf_list[i].quads[j][0]] ;
      surf_list[i].quads[j][1] = used[surf_list[i].quads[j][1]] ;
      surf_list[i].quads[j][2] = used[surf_list[i].quads[j][2]] ;
      surf_list[i].quads[j][3] = used[surf_list[i].quads[j][3]] ;
    }
    for(size_t j=0;j<surf_list[i].gen_faces.size();++j) {
      for(size_t k=0;k<surf_list[i].gen_faces[j].size();++k)
	surf_list[i].gen_faces[j][k] = used[surf_list[i].gen_faces[j][k]] ;
    }
  }

#ifdef DEBUG
  cout << "nnodes = "<< nnodes << endl ;



  for(int i=0;i<ssz;++i) {
    cout << "surf = " << surf_list[i].name << endl ;
    cout << "ntrias=" << surf_list[i].trias.size() 
    << ",nquads=" << surf_list[i].quads.size()
    << ",ngenfs=" << surf_list[i].gen_faces.size()
         << endl;
  }
#endif

}



void writeSurf(string filename_surf, string filename_name,
               vector<surface_info> &tmp_surf,
               vector<vector3d<double> > &tmp_p) {
  
  std::ofstream ofile(filename_surf.c_str(),std::ios::out) ;
  
  if(ofile.fail()) {
    cerr << "unable to open output file '" << filename_surf << "'" << endl ;
    Usage() ;
  }
  std::ofstream nfile(filename_name.c_str(),std::ios::out) ;
  if(nfile.fail()) {
    cerr << "unable to open output file '" << filename_name << "'" << endl ;
    Usage() ;
  }


  int ntri=0,nqua=0,ngen=0 ;
  for(size_t i=0;i<tmp_surf.size();++i) {
    ntri += tmp_surf[i].trias.size() ;
    nqua += tmp_surf[i].quads.size() ;
    ngen += tmp_surf[i].gen_faces.size() ;
  }

  ofile << ntri << ' ' << nqua << ' ' << tmp_p.size() << endl ;
  ofile.precision(14) ;
  double normal_spacing = 0 ;
  for(size_t i=0;i<tmp_p.size();++i) {
    ofile << tmp_p[i].x << ' ' << tmp_p[i].y << ' ' << tmp_p[i].z
	  << ' ' << normal_spacing << endl ;
  }
  
  // output triangle faces
  for(size_t i=0;i<tmp_surf.size();++i) {  
    for(size_t j=0;j<tmp_surf[i].trias.size();++j)
      ofile << tmp_surf[i].trias[j][0] << ' '
	    << tmp_surf[i].trias[j][1] << ' '
	    << tmp_surf[i].trias[j][2] << ' '
	    << tmp_surf[i].id << ' '
	    << 0 << ' '  // reconnection flag
	    << 0 << endl ; // bc flag
  }
  // output quad faces
  for(size_t i=0;i<tmp_surf.size();++i) {  
    for(size_t j=0;j<tmp_surf[i].quads.size();++j)
      ofile << tmp_surf[i].quads[j][0] << ' '
	    << tmp_surf[i].quads[j][1] << ' '
	    << tmp_surf[i].quads[j][2] << ' '
	    << tmp_surf[i].quads[j][3] << ' '
	    << tmp_surf[i].id << ' '
	    << 0 << ' '  // reconnection flag
	    << 0 << endl ; // bc flag
  }
  // Now write out general faces
  if(ngen > 0) {
    ofile << ngen << endl ;
    for(size_t i=0;i<tmp_surf.size();++i) {  
      for(size_t j=0;j<tmp_surf[i].gen_faces.size();++j) {
        size_t nf = tmp_surf[i].gen_faces[j].size() ;
        ofile << nf ;
        for(size_t k=0;k<nf;++k)
          ofile << ' ' << tmp_surf[i].gen_faces[j][k] ;
        ofile << endl ;
      }
    }
  }

  for(size_t i=0;i<tmp_surf.size();++i) {
    nfile << tmp_surf[i].id << ' ' << tmp_surf[i].name << endl ;
  }
  nfile.close() ;
  ofile.close();
}

void writeSurfaces(string filename,
                   vector<surface_info> &surf_list,
                   vector<vector3d<double> > &pos,
                   vector<int> &node_map) {
  std::ofstream ofile(filename.c_str(),std::ios::out) ;
  
  if(ofile.fail()) {
    cerr << "unable to open output file '" << filename << "'" << endl ;
    Usage() ;
  }
  
  size_t  npos = pos.size();
  
  ofile << npos << endl ;
  ofile.precision(14) ;

  //output pos
  for(size_t i=0;i<npos;++i) {
    ofile << pos[i].x << ' ' << pos[i].y << ' ' << pos[i].z<< endl ;
  }
  //output surf_list
  size_t nsurf = surf_list.size();
  ofile << nsurf << endl;

  for(size_t i = 0; i < nsurf; i++){
    ofile<<surf_list[i].name<<endl;
  }
  for(size_t i = 0; i < nsurf; i++){ 
    size_t ntris=surf_list[i].trias.size();
    size_t nquads = surf_list[i].quads.size();
    size_t ngens = surf_list[i].gen_faces.size();
    size_t id = surf_list[i].id;
    ofile << id<<' ' <<ntris<<' '<<nquads<<' ' <<ngens << endl;
    for(size_t j =0; j < ntris; j++){
      ofile << surf_list[i].trias[j][0] << ' '
	    << surf_list[i].trias[j][1] << ' '
            << surf_list[i].trias[j][2] << endl;
    }
    
    for(size_t j=0;j<nquads;++j){
      
     ofile << surf_list[i].quads[j][0] << ' '
           << surf_list[i].quads[j][1] << ' '
           << surf_list[i].quads[j][2] << ' '
           << surf_list[i].quads[j][3] << endl;
    }

    for(size_t j=0;j<ngens;++j){
      size_t nf = surf_list[i].gen_faces[j].size() ;
      ofile << nf ;
      for(size_t k=0;k<nf;++k)
        ofile << ' ' << surf_list[i].gen_faces[j][k] ;
      ofile << endl ;
    }
  }

  for(size_t i=0;i<npos;++i) {
    ofile << node_map[i] << endl;
  }
  ofile.close();  
  
}


// int main(int ac, char *av[]) {
//   using Loci::entitySet ;
//   using Loci::vector3d ;
//   Loci::Init(&ac, &av) ;
//   if(Loci::MPI_processes > 1) {
//     cerr << "vog2surf is not parallel! Run on only one processor!" << endl ;
//     Loci::Abort() ;
//   }
//   string surface_file = "";
//   string surf_file = "";
//   string input_file = "";
  
//   for(int i=1;i<ac;++i) {
//     string opt = av[i] ;
//     bool parsed = false ;
//    if(opt == "-o") {
//       i++ ;
//       if(i<ac) {
//         surf_file = string(av[i]) ;
//         parsed = true ;
//       }
//     } else if(opt == "-surface") {
//       i++ ;
//       if(i<ac) {
//         surface_file = string(av[i]) ;
//         parsed = true ;
//       }
//     } else
//       if(input_file == "") {
// 	input_file = string(av[i]) ;
// 	parsed = true ;
//       } else
// 	parsed = false ;
    
//     if(!parsed) {
//       cerr << "unable to parse command line argument '" << av[i] << "'" << endl ;
//       Usage() ;
//     }
//   }

 

//   string name_file = input_file + ".names" ;
//   string file_input = input_file + ".vog";
//   if(surface_file=="" && surf_file=="")surf_file=input_file+".surf";
  
  
  
// #define DEBUG
// #ifndef DEBUG
//   /* Save old error handler */
//   herr_t (*old_func)(void*) = 0;
//   void *old_client_data = 0 ;
//   H5Eget_auto(&old_func, &old_client_data);
  
//   /* Turn off error handling */
//   H5Eset_auto(NULL, NULL);
// #endif
  
//   vector<surface_info> tmp_surf ;
//   vector<vector3d<double> > tmp_p ;
//   readSurfaces(file_input,tmp_surf,tmp_p) ;


  
//   if(surface_file!="")writeSurfaces(surface_file, tmp_surf, tmp_p);
//   if(surf_file!="")writeSurf(surf_file,name_file, tmp_surf, tmp_p);

  

//   Loci::Finalize();
// }
int main(int ac, char *av[]) {
  using Loci::entitySet ;
  using Loci::vector3d ;
  vector<string> bclist ;
  Loci::Init(&ac, &av) ;
  if(Loci::MPI_processes > 1) {
    cerr << "vog2surf is not parallel! Run on only one processor!" << endl ;
    Loci::Abort() ;
  }
  if(ac < 2){
    Usage();
    Loci::Abort() ;
  }
  string surface_file = "";
  string surf_file = "";
  string input_file = "";
  
  for(int i=1;i<ac;++i) {
    string opt = av[i] ;
    bool parsed = false ;
   if(opt == "-o") {
      i++ ;
      if(i<ac) {
        surf_file = string(av[i]) ;
        parsed = true ;
      }
    } else if(opt == "-surface") {
      i++ ;
      if(i<ac) {
        surface_file = string(av[i]) ;
        parsed = true ;
      }
    } else
      if(input_file == "") {
	input_file = string(av[i]) ;
	parsed = true ;
      } else
	parsed = false ;
    
    if(!parsed) {
      cerr << "unable to parse command line argument '" << av[i] << "'" << endl ;
      Usage() ;
    }
  }

 

  string name_file = input_file + ".names" ;
  string file_input = input_file + ".vog";
  if(surface_file=="" && surf_file=="")surf_file=input_file+".surf";
  
  
  
#define H5DEBUG
#ifndef H5DEBUG
  /* Save old error handler */
  herr_t (*old_func)(void*) = 0;
  void *old_client_data = 0 ;
  H5Eget_auto(&old_func, &old_client_data);
  
  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);
#endif
  
  vector<surface_info> tmp_surf ;
  vector<vector3d<double> > tmp_p ;
  vector<int> tmp_map;
  readSurfaces(file_input,bclist,tmp_surf,tmp_p,tmp_map) ;


  
  if(surface_file!="")writeSurfaces(surface_file, tmp_surf, tmp_p, tmp_map);
  if(surf_file!="")writeSurf(surf_file,name_file, tmp_surf, tmp_p);

  

  Loci::Finalize();
}
