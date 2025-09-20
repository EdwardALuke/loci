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
#include <Tools/fpe.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

using std::istringstream ;
using std::map ;
using std::vector ;
using std::string ;
using std::cerr ;
using std::cout ;
using std::endl ;
using std::ifstream ;
using std::ios ;

unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

bool readVolTags(hid_t input_fid,
                 vector<pair<string,Loci::entitySet> > &volDat) {
  using namespace Loci ;
  /* Save old error handler */
  H5E_auto_t old_func = 0 ;
  void *old_client_data = 0 ;
#ifdef H5_USE_16_API
      H5Eget_auto(&old_func, &old_client_data);
      /* Turn off error handling */
      H5Eset_auto(NULL, NULL);
#else
      H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);
      /* Turn off error handling */
      H5Eset_auto(H5E_DEFAULT,NULL, NULL);
#endif


  vector<pair<string,entitySet> > volTags ;
#ifdef H5_USE_16_API
  hid_t cell_info = H5Gopen(input_fid,"cell_info") ;
#else
  hid_t cell_info = H5Gopen(input_fid,"cell_info",H5P_DEFAULT) ;
#endif
  if(cell_info > 0) {
    vector<string> vol_tag ;
    vector<entitySet> vol_set ;
    vector<int> vol_id ;
    
    hsize_t num_tags = 0 ;
    H5Gget_num_objs(cell_info,&num_tags) ;
    for(hsize_t tg=0;tg<num_tags;++tg) {
      char buf[1024] ;
      memset(buf, '\0', 1024) ;
      H5Gget_objname_by_idx(cell_info,tg,buf,sizeof(buf)) ;
      buf[1023]='\0' ;
      
      string name = string(buf) ;
#ifdef H5_USE_16_API
      hid_t vt_g = H5Gopen(cell_info,buf) ;
#else
      hid_t vt_g = H5Gopen(cell_info,buf,H5P_DEFAULT) ;
#endif
      hid_t id_a = H5Aopen_name(vt_g,"Ident") ;
      int ident ;
      H5Aread(id_a,H5T_NATIVE_INT,&ident) ;
      H5Aclose(id_a) ;
      entitySet dom ;
      HDF5_ReadDomain(vt_g,dom) ;
      vol_tag.push_back(name) ;
      vol_set.push_back(dom) ;
      vol_id.push_back(ident) ;
      H5Gclose(vt_g) ;
    }
    int maxi =0 ;
    for(size_t i=0;i<vol_id.size();++i)
      maxi = max(vol_id[i],maxi) ;
    vector<pair<string,entitySet> > tmp(maxi+1) ;
    volTags.swap(tmp) ;
    for(size_t i=0;i<vol_id.size();++i) {
      volTags[vol_id[i]].first = vol_tag[i] ;
      volTags[vol_id[i]].second = vol_set[i] ;
    }
  } else {
#ifdef H5_USE_16_API
    hid_t file_info = H5Gopen(input_fid,"file_info") ;
#else
    hid_t file_info = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
#endif
    long numCells = readAttributeLong(file_info,"numCells") ;
    volTags.push_back(pair<string,entitySet>
                      (string("Main"),
                       entitySet(interval(0,numCells-1)))) ;
    H5Gclose(file_info) ;
  }
  
  /* Restore previous error handler */
#ifdef H5_USE_16_API
  H5Eset_auto(old_func, old_client_data);
#else
  H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
#endif
  volDat.swap(volTags) ;
  return true ;
}

void writeVogGrid(string infile,string outfile) {
  if(Loci::MPI_rank==0) {
    hid_t output_fid = 0 ;
    output_fid = H5Fcreate(outfile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
    if(output_fid <= 0) {
      cerr << "unable to open output file '" << outfile << "'" << endl ;
      Loci::Abort() ;
    }
    hid_t input_fid = 0 ;
    input_fid = H5Fopen(infile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;

    // Copy file_info
#ifdef H5_USE_16_API
    hid_t fi = H5Gopen(input_fid,"file_info") ;
#else
    hid_t fi = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
#endif
    unsigned long numNodes = readAttributeLong(fi,"numNodes") ;
    unsigned long numFaces = readAttributeLong(fi,"numFaces") ;
    unsigned long numCells = readAttributeLong(fi,"numCells") ;
    H5Gclose(fi) ;

#ifdef H5_USE_16_API
    fi = H5Gcreate(output_fid,"file_info",0) ;
#else
    fi = H5Gcreate(output_fid,"file_info",
		   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    hsize_t dims = 1 ;
    hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
#ifdef H5_USE_16_API
    hid_t att_id = H5Acreate(fi,"numNodes",H5T_STD_I64BE,
                             dataspace_id,H5P_DEFAULT) ;
#else
    hid_t att_id = H5Acreate(fi,"numNodes",H5T_STD_I64BE,
                             dataspace_id,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    H5Awrite(att_id,H5T_NATIVE_LLONG,&numNodes) ;
    H5Aclose(att_id) ;

#ifdef H5_USE_16_API
    att_id = H5Acreate(fi,"numFaces", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT) ;
#else
    att_id = H5Acreate(fi,"numFaces", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    H5Awrite(att_id,H5T_NATIVE_LLONG,&numFaces) ;
    H5Aclose(att_id) ;
#ifdef H5_USE_16_API
    att_id = H5Acreate(fi,"numCells", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT) ;
#else
    att_id = H5Acreate(fi,"numCells", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    H5Awrite(att_id,H5T_NATIVE_LLONG,&numCells) ;
    H5Aclose(att_id) ;
    H5Gclose(fi) ;


    {
      // Copy Face Information
#ifdef H5_USE_16_API
      hid_t face_g = H5Gopen(input_fid,"face_info") ;
#else
      hid_t face_g = H5Gopen(input_fid,"face_info",H5P_DEFAULT) ;
#endif

      // Read cluster sizes
#ifdef H5_USE_16_API
      hid_t dataset = H5Dopen(face_g,"cluster_sizes") ;
#else
      hid_t dataset = H5Dopen(face_g,"cluster_sizes",H5P_DEFAULT) ;
#endif
      hid_t dspace = H5Dget_space(dataset) ;
      hsize_t size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      vector<unsigned short> cluster_sizes(size) ;
      hsize_t dimension = size ;
      hsize_t stride = 1 ;
#ifdef H5_INTERFACE_1_6_4
      hsize_t start = 0 ;
      hsize_t lstart = 0 ;
#else
      hssize_t start = 0 ;
      hssize_t lstart = 0 ;
#endif
      hsize_t count = size ;
      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
      int rank = 1 ;
      hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
      hid_t err = H5Dread(dataset,H5T_NATIVE_USHORT,memspace,dspace,H5P_DEFAULT,&cluster_sizes[0]) ;
      if(err < 0) {
        cerr << "unable to read cluster sizes from vog file" << endl; 
      }
      H5Dclose(dataset) ;
      H5Sclose(memspace) ;

#ifdef H5_USE_16_API
      hid_t face_id = H5Gcreate(output_fid,"face_info",0) ;
#else
      hid_t face_id = H5Gcreate(output_fid,"face_info",
				H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
      // Write cluster sizes
      rank = 1 ;
      dimension = cluster_sizes.size() ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
      start = 0 ;
      lstart = 0 ;

      stride = 1;
#ifdef H5_USE_16_API
      dataset = H5Dcreate(face_id,"cluster_sizes",H5T_NATIVE_USHORT,
                          dataspace,H5P_DEFAULT) ;
#else
      dataset = H5Dcreate(face_id,"cluster_sizes",H5T_NATIVE_USHORT,
                          dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
      memspace = H5Screate_simple(rank,&dimension,NULL) ;
      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&start,&stride,&dimension,NULL) ;
      err = H5Dwrite(dataset,H5T_NATIVE_USHORT,memspace,dataspace,H5P_DEFAULT,
                     &cluster_sizes[0]) ;
      if(err<0) {
        cerr << "unable to write cluster sizes to '" << outfile << "'" << endl ;
        exit(-1) ;
      }
      H5Sclose(memspace) ;
      H5Sclose(dataspace) ;
      H5Dclose(dataset) ;

#ifdef H5_USE_16_API
      dataset = H5Dopen(face_g,"cluster_info") ;
#else
      dataset = H5Dopen(face_g,"cluster_info",H5P_DEFAULT) ;
#endif
      dspace = H5Dget_space(dataset) ;
      dspace = H5Dget_space(dataset) ;
      size = 0 ;
      H5Sget_simple_extent_dims(dspace,&size,NULL) ;
      vector<unsigned char> cluster_info(size) ;
      dimension = size ;
      stride = 1 ;
      start = 0 ;
      lstart = 0 ;

      count = size ;
      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
      rank = 1 ;
      memspace = H5Screate_simple(rank,&dimension,NULL) ;
      err = H5Dread(dataset,H5T_NATIVE_UCHAR,memspace,dspace,H5P_DEFAULT,&cluster_info[0]) ;
      if(err < 0) {
        cerr << "unable to read cluster info from vog file" << endl; 
      }
      H5Dclose(dataset) ;
      H5Sclose(memspace) ;

      // Write cluster info
      rank = 1 ;
      dimension = cluster_info.size() ;
      dataspace = H5Screate_simple(rank,&dimension,NULL) ;
      start = 0 ;
      lstart = 0 ;

      stride = 1;
#ifdef H5_USE_16_API
      dataset = H5Dcreate(face_id,"cluster_info",H5T_NATIVE_UCHAR,
                          dataspace,H5P_DEFAULT) ;
#else
      dataset = H5Dcreate(face_id,"cluster_info",H5T_NATIVE_UCHAR,
                          dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
      memspace = H5Screate_simple(rank,&dimension,NULL) ;
      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&start,&stride,&dimension,NULL) ;
      err = H5Dwrite(dataset,H5T_NATIVE_UCHAR,memspace,dataspace,H5P_DEFAULT,
                     &cluster_info[0]) ;
      if(err<0) {
        cerr << "unable to write cluster sizes to '" << outfile << "'" << endl ;
        exit(-1) ;
      }
      H5Sclose(memspace) ;
      H5Sclose(dataspace) ;
      H5Dclose(dataset) ;
    
      H5Gclose(face_id) ;
      H5Gclose(face_g) ;
    }
    // Copy Boundary Information
    vector<pair<int,string> > boundary_ids ;
    Loci::readBCfromVOG(infile,boundary_ids) ;
    
#ifdef H5_USE_16_API
    hid_t surf_id = H5Gcreate(output_fid,"surface_info",0) ;
#else
    hid_t surf_id = H5Gcreate(output_fid,"surface_info",
			      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    vector<pair<int,string> >::const_iterator mi ;
    
    for(mi=boundary_ids.begin();mi!=boundary_ids.end();++mi) {
#ifdef H5_USE_16_API
      hid_t bc_id = H5Gcreate(surf_id,mi->second.c_str(),0) ;
#else
      hid_t bc_id = H5Gcreate(surf_id,mi->second.c_str(),
			      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
      hsize_t dims = 1 ;
      hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
      
#ifdef H5_USE_16_API
      hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                               dataspace_id, H5P_DEFAULT) ;
#else
      hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                               dataspace_id, H5P_DEFAULT, H5P_DEFAULT) ;
#endif
      H5Awrite(att_id,H5T_NATIVE_INT,&(mi->first)) ;
      H5Aclose(att_id) ;
      H5Gclose(bc_id) ;
    }
    H5Gclose(surf_id) ;


    // Load in new positions
    hid_t pos_id = H5Fopen("newpos.dat",H5F_ACC_RDONLY,H5P_DEFAULT) ;
    // Read in pos data from file i
    vector<Loci::vector3d<double> > pos_dat(numNodes) ;

#ifdef H5_USE_16_API
    hid_t pos_g = H5Gopen(pos_id,"pos") ;
#else
    hid_t pos_g = H5Gopen(pos_id,"pos",H5P_DEFAULT) ;
#endif

    int rank = 1 ;
    hsize_t dimension = numNodes ;

#ifdef H5_INTERFACE_1_6_4
    hsize_t start = 0 ;
    hsize_t lstart = 0 ;
#else
    hssize_t start = 0 ;
    hssize_t lstart = 0 ;
#endif

    hsize_t stride = 1 ;
    typedef Loci::data_schema_traits<Loci::vector3d<double> > traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;

#ifdef H5_USE_16_API
    hid_t dataset = H5Dopen(pos_g,"data") ;
#else
    hid_t dataset = H5Dopen(pos_g,"data",H5P_DEFAULT) ;
#endif
    hid_t dspace = H5Dget_space(dataset) ;
      
    hsize_t count = numNodes ;
    H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
    rank = 1 ;
    dimension = numNodes ;
    hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
                  &pos_dat[0]) ;
    if(err < 0) {
      cerr << "unable to read positions from '" << "newpos.dat" << "'" << endl ;
      exit(-1) ;
    }
    H5Sclose(dspace) ;
    H5Dclose(dataset) ;
    H5Fclose(pos_id) ;

    // Write out node information
#ifdef H5_USE_16_API
    hid_t node_id = H5Gcreate(output_fid,"node_info",0) ;
#else
    hid_t node_id = H5Gcreate(output_fid,"node_info",
			      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif

    rank = 1 ;
    dimension = numNodes ;
    hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

    start = 0 ;
    lstart = 0 ;

    stride = 1 ;
#ifdef H5_USE_16_API
    hid_t dataseto = H5Dcreate(node_id,"positions",dp->get_hdf5_type(),
                               dataspace,H5P_DEFAULT) ;
#else
    hid_t dataseto = H5Dcreate(node_id,"positions",dp->get_hdf5_type(),
                               dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
    count = numNodes ;
    // Now write out this file's part of the positions
    H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, &start,&stride,&count,NULL)  ;
    dimension = count ;
    rank = 1 ;
    memspace = H5Screate_simple(rank,&dimension,NULL) ;
    datatype = dp->get_hdf5_type() ;
    err = H5Dwrite(dataseto,datatype,memspace,dataspace,
                   H5P_DEFAULT,&pos_dat[0]) ;
    if(err < 0) {
      cerr << "unable to write positions to '" << outfile << "'" << endl ;
      exit(-1) ;
    }
    H5Sclose(memspace) ;
    H5Sclose(dataspace) ;
    H5Dclose(dataseto) ;
  }
}

void Usage() {
  cerr << "Usage:  adjustpos <options>" << endl 
       << " were <options> are:" << endl
       << "  -i <input grid>: specify input grid filename"<< endl
       << "  -o <output grid>: specify output grid filename" << endl
       << "  -g <geometry>: specify geometry file" << endl
       << "  -sym <bc>: specify symmetry plane bc's" << endl
       << "  -stiffness <val>: specify stiffness value (0-1)." << endl
       << endl ;
}

int main(int ac, char *av[]) {
  string query = "newPos" ;

  Loci::Init(&ac, &av) ;

  string infile ;
  string outfile ;
  string geomfile ;
  vector<string> sym_list ;
  double epsilon = -1 ;

  bool error = false ;

  for(int i=1;i<ac;++i) {
    if(!strcmp(av[i],"-q")) {
      i++ ;
      query = av[i] ;
    } else if(!strcmp(av[i],"-i")) {
      i++ ;
      infile = av[i] ;
    } else if(!strcmp(av[i],"-o")) {
      i++ ;
      outfile = av[i] ;
    } else if(!strcmp(av[i],"-g")) {
      i++ ;
      geomfile = av[i] ;
    } else if(!strcmp(av[i],"-sym")) {
      i++ ;
      string v(av[i]) ;
      if(av[i][0] >= '0' && av[i][0] <= '9')
        v = "BC_"+v ;
      sym_list.push_back(v) ;
    } else if(!strcmp(av[i],"-stiffness")) {
      i++ ;
      epsilon = atof(av[i]) ;
    } else {
      cerr << "bad option " << av[i] << endl ;
      error = true ;
    }
  }

  if(infile == "") {
    cerr << "did not specify input grid file with '-i'" << endl ;
    error = true ;
  }

  if(outfile == "") {
    cerr << "did not specify output grid file with '-o' " << endl ;
    error = true ;
  }

  if(geomfile == "") {
    cerr << "did not specify geometry file with '-g'" << endl ;
    error = true ;
  }

  if(error) {
    Usage() ;
    return -1 ;
  }
  
  fact_db facts ;

  rule_db rdb ;
  rdb.add_rules(global_rule_list) ;
  Loci::load_module("fvm",rdb) ;

 

  param<std::string> geomfile_par ;
  *geomfile_par = geomfile;
  facts.create_fact("geomfile_par",geomfile_par) ;
    
  if(!Loci::setupFVMGrid(facts,infile)) {
    cerr << "unable to read file " << infile << endl ;
    Loci::Abort() ;
  }

  store<string> boundary_names ;
  boundary_names = facts.get_variable("boundary_names") ;

  map<string,string> binfo ;
  entitySet dom = boundary_names.domain() ;
  FORALL(dom,cc) {
    binfo[boundary_names[cc]] = string("geometry") ;
  } ENDFORALL ;
  for(size_t i = 0; i<sym_list.size();++i)
    binfo[sym_list[i]] = string("symmetry") ;

  try {
    string var_str ;
    var_str = "{ boundary_conditions: <" ;
    map<string,string>::const_iterator si ;
    for(si = binfo.begin();si!=binfo.end();) {
      var_str += si->first + "=" + si->second ;
      ++si ;
      if(si != binfo.end())
        var_str += "," ;
    }
    var_str += "> } " ;
    //    cout << "var_str="<<var_str <<endl ;
    istringstream ifile(var_str) ;

    facts.read_vars(ifile,rdb) ;
  } catch(const Loci::BasicException &err) {
    err.Print(cerr) ;
    cerr << "aborted in vars setup" << endl ;
    Loci::Abort() ;
  }

  if(epsilon > 0) {
    param<double> stiffness ;
    *stiffness = epsilon ;
    facts.update_fact("stiffness",stiffness) ;
  }
  
  setupBoundaryConditions(facts) ;

  constraint bc1 ;
  bc1 = facts.get_fact("geometry_BC") ;
  multiMap face2node ;
  face2node = facts.get_fact("face2node") ;
  Loci::MapRepP mp = Loci::MapRepP(face2node.Rep()) ;
  entitySet bcnodes = mp->image(*bc1) ;
  bcnodes = Loci::all_collect_entitySet(bcnodes) ;
  constraint boundary_nodes ;
  *boundary_nodes = bcnodes ;
  facts.create_fact("boundary_nodes",boundary_nodes) ;

  store<vector3d<double> > pos ;
  pos = facts.get_fact("pos") ;
  bcnodes &= pos.domain() ;
  vector<Loci::kdTree::coord3d> geom_pts(bcnodes.size()) ;
  vector<int> geom_ids(bcnodes.size()) ;
  int cnt = 0 ;
  FORALL(bcnodes,nd) {
    geom_pts[cnt][0] = pos[nd].x ;
    geom_pts[cnt][1] = pos[nd].y ;
    geom_pts[cnt][2] = pos[nd].z ;
    geom_ids[cnt] = nd ;
    cnt++ ;
  } ENDFORALL ;

  entitySet nodes = pos.domain() ;
  vector<Loci::kdTree::coord3d> node_pts(nodes.size()) ;
  vector<int> closest(nodes.size(),-1) ;
  cnt = 0 ;
  FORALL(nodes,nd) {
    node_pts[cnt][0] = pos[nd].x ;
    node_pts[cnt][1] = pos[nd].y ;
    node_pts[cnt][2] = pos[nd].z ;
    cnt++ ;
  } ENDFORALL ;

  Loci::parallelNearestNeighbors(geom_pts,geom_ids,node_pts,closest,
                                 MPI_COMM_WORLD) ;

  Map node2surface ;
  node2surface.allocate(nodes) ;
  cnt = 0 ;
  FORALL(nodes,nd) {
    node2surface[nd] = closest[cnt] ;
    
    cnt++ ;
  } ENDFORALL ;

  
  facts.create_fact("node2surface",node2surface) ;

  if(!Loci::makeQuery(rdb,facts,query)) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
  
  if(query == "newPos") {
    writeVogGrid(infile,outfile) ;
  }

  Loci::Finalize() ;
}
