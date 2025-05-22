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
#include <stdio.h>
#include <strings.h>
#include <strings.h>
#include <Loci.h>
#include "vogtools.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <limits>

using std::map ;
using std::vector ;
using std::string ;
using std::cout ;
using std::cerr ;
using std::endl ;
using std::ifstream;
using Loci::vector3d ;
using Loci::entitySet;
typedef vector3d<double> vect3d;
using std::ofstream;

const int NEW_BC_ID = std::numeric_limits<int>::min();


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
  H5E_auto_t old_func = 0;
  void *old_client_data = 0 ;

  H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);
  /* Turn off error handling */
  H5Eset_auto(H5E_DEFAULT,NULL, NULL);

  vector<pair<string,entitySet> > volTags ;

  hid_t cell_info = H5Gopen(input_fid,"cell_info",H5P_DEFAULT) ;
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

      hid_t vt_g = H5Gopen(cell_info,buf,H5P_DEFAULT) ;
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

    hid_t file_info = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
    long numCells = readAttributeLong(file_info,"numCells") ;
    volTags.push_back(pair<string,entitySet>
                      (string("Main"),
                       entitySet(interval(0,numCells-1)))) ;
    H5Gclose(file_info) ;
  }
  
  /* Restore previous error handler */
  H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
  volDat.swap(volTags) ;
  return true ;
}

namespace Loci {
  extern int getClusterNumFaces(unsigned char *cluster) ;
  extern entitySet faceCluster(const multiMap &face2node,
                               const Map &cl, const Map &cr, entitySet faces,
                               vector<unsigned char> &cluster_info,
                               vector<unsigned short> &cluster_sizes);
  extern int fillClusterFaceSizes(unsigned char *cluster, int *sizes) ;
  int fillFaceInfo(unsigned char *cluster, multiMap &face2node,
                   Map &cl, Map &cr, int face_base) ;

}


//if a node is outside the sphere
inline bool mark_node(const vect3d& p, const vect3d& center, double r){
  return (dot(p-center, p-center) >= r*r);
}

//if a face is outside the sphere:
//if all the nodes of a face are outside the sphere, then the face is outside the sphere
inline bool mark_face(const vector<bool>& nodestat, const Loci::multiMap& face2node, Loci::Entity f){
  
  int fsz = face2node[f].size() ;
  bool isOut = true;
  for(int i = 0; i < fsz; i++){
    if(!nodestat[face2node[f][i]]){
      isOut = false;
      break;
    }
  }
  return isOut;
}



//read in the positions of the nodes 
bool readNodes(hid_t input_fid, 
               vector<vect3d>& pos_dat //pos of all nodes
               ) {
      
  // read in positions

  hid_t fid = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;

  unsigned long numNodes = readAttributeLong(fid,"numNodes") ;
  H5Gclose(fid) ;
  
  hsize_t  count = numNodes ;

  hsize_t lstart = 0 ;
  // Read in pos data from file i
  pos_dat.resize(numNodes) ;


  hid_t node_g = H5Gopen(input_fid,"node_info",H5P_DEFAULT) ;
  hid_t dataset = H5Dopen(node_g,"positions",H5P_DEFAULT) ;
  hid_t dspace = H5Dget_space(dataset) ;
  
  hsize_t stride = 1 ;
   
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
  int rank = 1 ;
  hsize_t dimension = count ;
  hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
  typedef Loci::data_schema_traits<vect3d > traits_type ;
  Loci::DatatypeP dp = traits_type::get_type() ;
  hid_t datatype = dp->get_hdf5_type() ;
  hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
		      &pos_dat[0]) ;
  if(err < 0) {
    cerr << "unable to read positions "<< endl ;
    exit(-1) ;
  }
  H5Sclose(dspace) ;
  H5Dclose(dataset) ;
  H5Gclose(node_g) ;
  return true;
}

bool readTags(string tagfile, vector<bool>& tags){
  
  ifstream inFile;
  inFile.open(tagfile.c_str());
  if(!inFile){
    cerr <<"can not open " << tagfile<< " for input" << endl;
    return false;
  }
  short int  tag; 
  for(unsigned int i = 0; i < tags.size(); i++){
    if(inFile >> tag) tags[i] = (tag != 0);
    else tags[i] = false;
  }
  return true;
}
// bool writeTags(string tagfile, const vector<bool>& tags){
  
//   ofstream outFile;
//   outFile.open(tagfile.c_str());
//   if(!outFile){
//     cerr <<"can not open " << tagfile<< " for output" << endl;
//     return false;
//   }
   
//   for(unsigned int i = 0; i < tags.size(); i++){
//     if(tags[i]) outFile << "1 ";
//     else outFile<<"0 ";
//   }
//   return true;
// }

bool createMaps(hid_t input_fid,
                const vector<bool> nodestat, //the IN/OUT state of all nodes in the original grid
                vector<int>& nodeMap, //the index of nodes
                vector<int>& cellMap, //the index of cells
                unsigned long &newNumNodes, 
                unsigned long &newNumFaces,
                unsigned long &newNumCells,
                vector<vector<int> >& face_data, //the faces of new grids
                Loci::entitySet& bc_ids //the boundary ids of the new grid
                ){
  
  hid_t fid = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
  unsigned long numNodes = readAttributeLong(fid,"numNodes");
  unsigned long numCells = readAttributeLong(fid,"numCells");
  H5Gclose(fid) ;
  
  //initialize maps to -2
  {
    vector<int>(numNodes, -2).swap(nodeMap);
    vector<int>(numCells, -2).swap(cellMap);
  }

    //first time read through face_info to select cells to keep
   
    hid_t face_g  = H5Gopen(input_fid,"face_info",H5P_DEFAULT) ;
    hid_t  dataset = H5Dopen(face_g,"cluster_sizes",H5P_DEFAULT) ;
    hid_t dspace = H5Dget_space(dataset) ;
    hsize_t stride = 1 ;
    hsize_t size = 0 ;
    H5Sget_simple_extent_dims(dspace,&size,NULL) ;
    vector<unsigned short> csizes(size) ;
    hsize_t  dimension = size ;
 

    hsize_t start = 0 ;
    hsize_t  count = size ;
    H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
    int rank = 1 ;
    hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
    hid_t err = H5Dread(dataset,H5T_NATIVE_USHORT,memspace,dspace,H5P_DEFAULT,&csizes[0]) ;
    if(err < 0) {
      cerr << "unable to read cluster sizes from file" << endl ;
      exit(-1) ;
    }
    H5Dclose(dataset) ;
    H5Sclose(memspace) ;
  
  {
    // Read in clusters 

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
        cerr << "unable to read cluster from file " << endl ;
      }
      start += count ;
      // Now scan cluster into local buffers

      int nfaces = Loci::getClusterNumFaces(&cluster[0]) ;
      Loci::entitySet fclust = Loci::interval(0,nfaces-1) ;
      Loci::store<int> fcnts ;
      fcnts.allocate(fclust) ;
      Loci::fillClusterFaceSizes(&cluster[0],&fcnts[0]) ;
      multiMap face2node ;
      Map cl, cr ;
      face2node.allocate(fcnts) ;
      cl.allocate(fclust) ;
      cr.allocate(fclust) ;
      Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0);
      H5Sclose(memspace) ;
    
      // Now loop over faces in cluster to select cells to keep
      for(int f=0;f<nfaces;++f) {
        bool isOut = mark_face(nodestat, face2node, f);
        if(!isOut){//if any node of the face is in, both cl and cr is kept
          cellMap[cl[f]] = -1;
          if(cr[f]>=0) cellMap[cr[f]] = -1;
        }
      }//end loop over faces
    }//end loop over cluster
  
  
  }

  
  //second time read through face_info to select nodes and faces to keep
 

  {
    // Read in clusters 
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
        cerr << "unable to read cluster from file " << endl ;
      }
      start += count ;
      // Now scan cluster into local buffers

      int nfaces = Loci::getClusterNumFaces(&cluster[0]) ;
      Loci::entitySet fclust = Loci::interval(0,nfaces-1) ;
      Loci::store<int> fcnts ;
      fcnts.allocate(fclust) ;
      Loci::fillClusterFaceSizes(&cluster[0],&fcnts[0]) ;
      multiMap face2node ;
      Map cl,cr ;
      face2node.allocate(fcnts) ;
      cl.allocate(fclust) ;
      cr.allocate(fclust) ;
      Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0);
      H5Sclose(memspace) ;
    
      // Now loop over faces in cluster to select faces and nodes to keep
      for(int f=0;f<nfaces;++f) {
        if(cr[f] < 0 && cellMap[cl[f]] == -1){ //boundary face in old grid
          bc_ids += -cr[f];
        
          int sz = face2node.num_elems(f); 
          vector<int> tmpVec(sz+2);
          tmpVec[0] = cl[f];
          tmpVec[1] = cr[f];
          for(int i = 0; i < sz; i++){
            tmpVec[i+2] = face2node[f][i];
            nodeMap[face2node[f][i]] = -1;
          }
          face_data.push_back(tmpVec);
        }else if(cr[f] >= 0 &&(cellMap[cl[f]] == -1 || cellMap[cr[f]]==-1)){//inner face in old grid
          if(cellMap[cl[f]] == -1 && cellMap[cr[f]] == -1){//still inner face in new grid
            int sz = face2node.num_elems(f); 
            vector<int> tmpVec(sz+2);
            tmpVec[0] = cl[f];
            tmpVec[1] = cr[f];
            for(int i = 0; i < sz; i++){
              tmpVec[i+2] = face2node[f][i];
              nodeMap[face2node[f][i]] = -1;
            }
            face_data.push_back(tmpVec);
          }else if(cellMap[cl[f]] == -1 && cellMap[cr[f]] == -2){//boundary face in new grid, no flipping
          
            int sz = face2node.num_elems(f); 
            vector<int> tmpVec(sz+2);
            tmpVec[0] = cl[f];
            tmpVec[1] = NEW_BC_ID;
            for(int i = 0; i < sz; i++){
              tmpVec[i+2] = face2node[f][i];
              nodeMap[face2node[f][i]] = -1;
            }
            face_data.push_back(tmpVec);
          }else if(cellMap[cl[f]] == -2 && cellMap[cr[f]] == -1){//boundary face in new grid,  flipping
            int sz = face2node.num_elems(f);
            for(int n=0;n<sz/2;++n) {
              std::swap(face2node[f][n],face2node[f][sz-1-n]) ;
            }
          
            vector<int> tmpVec(sz+2);
            tmpVec[0] = cr[f];
            tmpVec[1] = NEW_BC_ID;
            for(int i = 0; i < sz; i++){
              tmpVec[i+2] = face2node[f][i];
              nodeMap[face2node[f][i]] = -1;
            }
            face_data.push_back(tmpVec);
          }
        }
    
      }//end loop over faces
    }//end loop over cluster
    H5Gclose(face_g) ;
  }
   
  // Count nodes 
  int nnodes = 0 ;
  for(unsigned long i = 0; i < numNodes; i++){
    if(nodeMap[i]==-1) nodeMap[i] = nnodes++;
  }
  newNumNodes = nnodes;
  
  int ncells = 0;
  for(unsigned long i = 0; i < numCells; i++){
    if(cellMap[i]==-1) cellMap[i] = ncells++;
  }
  newNumCells = ncells;
  newNumFaces = face_data.size();

  
  if(nnodes==0 ||ncells==0||newNumFaces==0)return false;  
  else return true; 
}



void Usage() {
  cout << "******************************************************************************" << endl ;

  cout << "vogcut is a utility for extracting  a sphere region from a grid." << endl
       << endl ;
  cout << "Usage: " << endl
       << "  vogcut <options>" << endl
       << "where <options> may be given by:" << endl 
       <<  endl ;
  cout << "  -g <inputfile>           | Specify input grid filename " << endl
       << "  -o <outputfile>          | Specify output grid filename" << endl
       << "  -not                     | complement the geometry defn" << endl 
       << "  -tag <tagfile>           | Specify nodal cut tag filename " <<endl
       << "  -sphere <args>           | Specify sphere geometry" << endl
       << "  -cylinder <args>         | Specify cylinder geometry" << endl
       << "  -cone <args>             | Specify cone geometry" << endl
       << "  -box <args>              | Specify box geometry" << endl << endl ;
  cout << "geometry arguments are as follows:" << endl << endl
       << " for -sphere:" << endl
       << "      -center <x> <y> <z> " << endl
       << "      -radius <r> << endl " << endl << endl
       << " for -cylinder:" << endl 
       << "      -pt1 <x> <y> <z> " << endl
       << "      -pt2 <x> <y> <z> " << endl
       << "      -radius <r> << endl " << endl << endl
       << " for -cone:" << endl 
       << "      -pt1 <x> <y> <z> " << endl
       << "      -pt2 <x> <y> <z> " << endl
       << "      -r1 <r> << endl " << endl 
       << "      -r2 <r> << endl " << endl << endl
       << " for -box:" << endl 
       << "      -pt1 <x> <y> <z> " << endl
       << "      -pt2 <x> <y> <z> " << endl << endl ;

  
  cout << "Note multiple geometries can be specified and the union of all of the" << endl 
       << " geometric entities will define the region to be cut." <<endl
       << " For more complex cases a nodal tag file can be used to control the mesh cut." << endl << endl ;
  
  cout << "example1: vogcut -g grid.vog -o grid_cut.vog -sphere -center 0 0 0 -radius 1"<<endl;
  cout << "example2: vogcut -g grid.vog -o grid_cut.vog -tag grid.tag"<<endl<<endl;
  cout << "******************************************************************************" << endl ;
}


typedef enum {SPHERE,CYLINDER,CONE,BOX} geoType ;
struct geoInfo {
  geoType type ;
  vect3d pt1, pt2 ;
  double r1, r2 ;
  void setSphere(vect3d pt, double r) {
    type = SPHERE ;
    pt1 = pt ;
    pt2 = pt ;
    r1 = r ;
    r2 = r ;
  }
  void setCylinder(vect3d pti1, vect3d pti2, double r) {
    type = CYLINDER ;
    pt1 = pti1 ;
    pt2 = pti2 ;
    r1 = r ;
    r2 = r ;
  }
  void setCone(vect3d pti1, vect3d pti2, double ri1, double ri2 ) {
    type = CONE ;
    pt1 = pti1 ;
    pt2 = pti2 ;
    r1 = ri1 ;
    r2 = ri2 ;
  }
  void setBox(vect3d pti1, vect3d pti2) {
    type = BOX ;
    pt1 = pti1 ;
    pt2 = pti2 ;
    r1 = 0 ;
    r2 = 0 ;
  }    
  geoInfo() {
    type = SPHERE ;
    pt1=vect3d(0,0,0) ;
    pt2=pt1 ;
    r1 = 0 ;
    r2 = 0 ;
  }

} ;

bool inGeo(const geoInfo &geo, const vect3d pt) {
  vect3d ptaxis(0,0,0),n(0,0,0) ;
  double r = 0,t = 0,alen = 1 ;
  switch(geo.type) {
  case SPHERE:
    return (norm(pt-geo.pt1) <= geo.r1) ;
  case CYLINDER:
  case CONE:
    n = geo.pt2-geo.pt1 ;
    alen = norm(n) ;
    n *= 1./alen ;
    // check endcaps
    if(dot(n,pt-geo.pt1) < 0 || dot(n,pt-geo.pt2)>0)
      return false ;
    t = dot(pt-geo.pt1,n)/alen ;
    r = geo.r1*(1.-t)+geo.r2*t ;
    ptaxis = geo.pt1*(1.-t)+geo.pt2*t ;
    return (norm(pt-ptaxis) <= r ) ;
  case BOX:
    return (pt.x <= max(geo.pt1.x,geo.pt2.x) &&
            pt.x >= min(geo.pt1.x,geo.pt2.x) &&
            pt.y <= max(geo.pt1.y,geo.pt2.y) &&
            pt.y >= min(geo.pt1.y,geo.pt2.y) &&
            pt.z <= max(geo.pt1.z,geo.pt2.z) &&
            pt.z >= min(geo.pt1.z,geo.pt2.z)) ;
  }
  return false ;
}


int main(int ac, char *av[]) {
  using Loci::entitySet ;
  using Loci::vector3d ;
  Loci::Init(&ac, &av) ;
  if(Loci::MPI_processes > 1) {
    cerr << "vogcut is not parallel! Run on only one processor!" << endl ;
    Loci::Abort() ;
  }
  string output_file="";
  string input_file="" ;
  string tag_file="";
  vector<geoInfo> vecGeoInfo ;
  
  bool tag_opt = false;
  bool complement_geom = false ;
  for(int i=1;i<ac;++i) {
    string opt = av[i] ;
    bool parsed = false ;
    if(opt == "-g") {
      i++ ;
      if(i<ac) {
        input_file = string(av[i]) ;
        parsed = true ;
      }
    }
    
    if(opt == "-not") {
      complement_geom = true ;
      parsed = true ;
    }
    if(opt == "-sphere") {
      bool haspt = false ;
      bool hasr = false ;
      vect3d pt(0,0,0) ;
      double r=0 ;
      while(i+1 < ac) {
        string arg = av[i+1] ;
        if(arg == "-center") {
          if(i+4 < ac) {
            pt.x = atof(av[i+2]) ;
            pt.y = atof(av[i+3]) ;
            pt.z = atof(av[i+4]) ;
            haspt = true ;
            i = i+4 ;
          } else {
            cerr << "unable to parse center after '-sphere'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else if(arg == "-radius") {
          if(i+2 < ac) {
            r = atof(av[i+2]) ;
            hasr = true ;
            i = i+2 ;
          } else {
            cerr << "unable to parse radius after '-sphere'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else {
          break ;
        }
      }
      if(haspt && hasr) {
        geoInfo sphere ;
        sphere.setSphere(pt,r);
        vecGeoInfo.push_back(sphere) ;
      } else {
        cerr << "expected -center and -radius after -sphere" << endl ;
        Usage() ;
        exit(-1) ;
      }
      parsed = true ;
    }
    if (opt == "-cylinder") {
      bool haspt1 = false,haspt2 = false ; 
      bool hasr = false ;
      vect3d pt1(0,0,0) ,pt2(0,0,0) ;
      double r=0 ;
      while(i+1 < ac) {
        string arg = av[i+1] ;
        if(arg == "-pt1") {
          if(i+4 < ac) {
            pt1.x = atof(av[i+2]) ;
            pt1.y = atof(av[i+3]) ;
            pt1.z = atof(av[i+4]) ;
            haspt1 = true ;
            i = i+4 ;
          } else {
            cerr << "unable to parse pt1 after '-cylinder'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else if(arg == "-pt2") {
          if(i+4 < ac) {
            pt2.x = atof(av[i+2]) ;
            pt2.y = atof(av[i+3]) ;
            pt2.z = atof(av[i+4]) ;
            haspt2 = true ;
            i = i+4 ;
          } else {
            cerr << "unable to parse pt2 after '-cylinder'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else if(arg == "-radius") {
          if(i+2 < ac) {
            r = atof(av[i+2]) ;
            hasr = true ;
            i = i+2 ;
          } else {
            cerr << "unable to parse radius after '-cylinder'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else {
          break ;
        }
      }
      if(haspt1 && haspt2 && hasr) {
        geoInfo cylinder ;
        cylinder.setCylinder(pt1,pt2,r);
        vecGeoInfo.push_back(cylinder) ;
      } else {
        cerr << "expected -pt1, -pt2,  and -radius after -cylinder" << endl ;
        Usage() ;
        exit(-1) ;
      }
      parsed = true ;
    }
    if (opt == "-cone") {
      bool haspt1 = false,haspt2 = false ; 
      bool hasr1,hasr2 = false ;
      vect3d pt1(0,0,0) ,pt2(0,0,0) ;
      double r1=0,r2=0 ;
      while(i+1 < ac) {
        string arg = av[i+1] ;
        if(arg == "-pt1") {
          if(i+4 < ac) {
            pt1.x = atof(av[i+2]) ;
            pt1.y = atof(av[i+3]) ;
            pt1.z = atof(av[i+4]) ;
            haspt1 = true ;
            i = i+4 ;
          } else {
            cerr << "unable to parse pt1 after '-cone'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else if(arg == "-pt2") {
          if(i+4 < ac) {
            pt2.x = atof(av[i+2]) ;
            pt2.y = atof(av[i+3]) ;
            pt2.z = atof(av[i+4]) ;
            haspt2 = true ;
            i = i+4 ;
          } else {
            cerr << "unable to parse pt2 after '-cone'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else if(arg == "-r1") {
          if(i+2 < ac) {
            r1 = atof(av[i+2]) ;
            hasr1 = true ;
            i = i+2 ;
          } else {
            cerr << "unable to parse r1 after '-cylinder'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else if(arg == "-r2") {
          if(i+2 < ac) {
            r2 = atof(av[i+2]) ;
            hasr2 = true ;
            i = i+2 ;
          } else {
            cerr << "unable to parse r2 after '-cylinder'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else {
          break ;
        }
      }
      if(haspt1 && haspt2 && hasr1 && hasr2) {
        geoInfo cone ;
        cone.setCone(pt1,pt2,r1,r2);
        vecGeoInfo.push_back(cone) ;
      } else {
        cerr << "expected -pt1, -pt2, -r1, and -r2 after -cone" << endl ;
        Usage() ;
        exit(-1) ;
      }
      parsed = true ;
    }
    if (opt == "-box") {
      bool haspt1 = false,haspt2 = false ; 
      vect3d pt1(0,0,0) ,pt2(0,0,0) ;
      while(i+1 < ac) {
        string arg = av[i+1] ;
        if(arg == "-pt1") {
          if(i+4 < ac) {
            pt1.x = atof(av[i+2]) ;
            pt1.y = atof(av[i+3]) ;
            pt1.z = atof(av[i+4]) ;
            haspt1 = true ;
            i = i+4 ;
          } else {
            cerr << "unable to parse pt1 after '-box'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else if(arg == "-pt2") {
          if(i+4 < ac) {
            pt2.x = atof(av[i+2]) ;
            pt2.y = atof(av[i+3]) ;
            pt2.z = atof(av[i+4]) ;
            haspt2 = true ;
            i = i+4 ;
          } else {
            cerr << "unable to parse pt2 after '-box'" << endl ;
            Usage() ;
            exit(-1) ;
          }
        } else {
          break ;
        }
      }
      if(haspt1 && haspt2) {
        geoInfo box ;
        box.setBox(pt1,pt2);
        vecGeoInfo.push_back(box) ;
      } else {
        cerr << "expected -pt1 and -pt2 -box" << endl ;
        Usage() ;
        exit(-1) ;
      }
      parsed = true ;
    }
    
    if(opt == "-tag") {
      i++ ;
      if(i<ac) {
        tag_file = string(av[i]);
        parsed = true ;
        tag_opt = true;
      }
    }  
           
    if(opt == "-o") {
      i++ ;
      if(i<ac) {
        output_file = string(av[i]) ;
        parsed = true ;
      }
    }
    
    if(!parsed) {
      cerr << "unable to parse command line argument '" << av[i] << "'" << endl ;
      Usage() ;
      exit(-1) ;
    }
  }

   
  if(input_file.empty()) {
    cerr << "no input grid files to process!" << endl ;
    Usage() ;
    exit(-1) ;
  }


  

  bool fail = false ;

#define DEBUG_HDF5
#ifndef DEBUG_HDF5
  /* Save old error handler */
  H5E_auto_t old_func = 0;
  void *old_client_data = 0 ;
  H5Eget_auto(H5E_DEFAULT,&old_func, &old_client_data);
  /* Turn off error handling */
  H5Eset_auto(H5E_DEFAULT,NULL, NULL);
#endif

  {
    size_t p = input_file.find('.');
    if(p == string::npos)input_file = input_file+".vog";
  }
    
  if(output_file==""){
    size_t p = input_file.find('.');
    output_file = input_file.substr(0, p)+"_cut.vog";
  }
  
  
  hid_t input_fid= H5Fopen(input_file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(input_fid <= 0) {
    cerr << "unable to open file '" << input_file << "'"<< endl ;
    exit(-1) ;
    fail = true ;
  }
  if(fail)
    exit(-1) ; 
  
  vector<vect3d> pos; //position of nodes in input grid
  //read in positions
  readNodes(input_fid, 
            pos);
 
  //mark the nodes
  unsigned long inNumNodes = pos.size();
  // note nodestat is true for nodes we wish to CUT out of the grid
  vector<bool> nodestat(inNumNodes, true);
  if(tag_opt)readTags(tag_file, nodestat);

  if(vecGeoInfo.size() > 0) {
    for(unsigned long i = 0; i < inNumNodes; i++) {
      bool inGeom = false ;
      for(size_t j=0;j<vecGeoInfo.size();++j) 
        inGeom = inGeom || inGeo(vecGeoInfo[j],pos[i]) ;
      bool cut = !inGeom ;
      
      nodestat[i] = nodestat[i] & cut ;
      if(complement_geom) {
	nodestat[i] = !nodestat[i] ;
      }
    }
  }
  
  vector<int> nodeMap;
  vector<int> cellMap;
  unsigned long newNumCells, newNumFaces, newNumNodes;
  vector<vector<int> > face_data;
  Loci::entitySet bc_ids;
  if(!createMaps(input_fid,
                 nodestat,
                 nodeMap,
                 cellMap,
                 newNumNodes,
                 newNumFaces,
                 newNumCells,
                 face_data,
                 bc_ids
                 )){
    cerr<<"ERROR: nothing to cut!"<< endl;
    return -1;
  }

   
  
  //get new positions
  store<vector3d<double> > act_pos_dat;
  entitySet act_nodes = interval(0, newNumNodes-1);
  act_pos_dat.allocate(act_nodes) ;
  size_t local_count = 0; 
  for(size_t ii=0;ii<pos.size();++ii){
    if(nodeMap[ii]!= -2){
      act_pos_dat[local_count++] = pos[ii] ;
    }
  } 
    
  //define the new boundary id
    
  int new_bc_id = 1;
  if(bc_ids.size()>0)new_bc_id = bc_ids.Max()+1;
    
  //create maps
  multiMap face2node;
  Map cl, cr;
  store<int> fcount;
  entitySet faces = interval(0, newNumFaces-1);
  cl.allocate(faces);
  cr.allocate(faces);
  fcount.allocate(faces);
  bool new_bc_id_used = false;
    
  for(unsigned long f = 0; f < newNumFaces; f++){
    fcount[f] = face_data[f].size()-2;
    cl[f] = cellMap[face_data[f][0]];
    if(face_data[f][1]>=0) cr[f] =cellMap[face_data[f][1]];
    else cr[f] = face_data[f][1];
    if(cr[f]==NEW_BC_ID){
      cr[f] = -new_bc_id;
      new_bc_id_used = true;
    }
  }
  face2node.allocate(fcount);
    
  for(unsigned long f = 0; f < newNumFaces; f++){
    for(int i = 0; i < fcount[f]; i++){
      face2node[f][i] = nodeMap[face_data[f][i+2]];
    }
  }
  VOG::colorMatrix(act_pos_dat,cl,cr,face2node) ;
    
  // Setup volume tags
    
  vector<pair<string,entitySet> > volTags ;
  vector<pair<string,entitySet> > volDat ;
  readVolTags(input_fid,volDat) ;
  if(!volDat.empty()){
    for(size_t j=0;j<volDat.size();++j){
      string name = volDat[j].first;
      entitySet domain = volDat[j].second;
      entitySet image;
      for(entitySet::const_iterator ei = domain.begin(); ei != domain.end(); ei++){
        if(cellMap[*ei] >=0)image += cellMap[*ei];
      }
      if(image != EMPTY) volTags.push_back(make_pair(volDat[j].first, image));
    }
        

    cout << "volume tags =" << endl  ;
    for(size_t i=0;i<volTags.size();++i)
      cout << i << " - " << volTags[i].first << volTags[i].second << endl ;
  }
  //read in boundary names 
  vector<pair<int,string> > boundary_ids;
  Loci::readBCfromVOG(input_file,boundary_ids) ;
  std::map<int,string> bc_map;
  for(unsigned int i = 0; i < boundary_ids.size(); i++){
    bc_map[boundary_ids[i].first] = boundary_ids[i].second;
  }
        
  vector<pair<int,string> > surf_ids ;
  if(new_bc_id !=1){
    bool bc_map_empty = (boundary_ids.size()==0);
    for(entitySet::const_iterator ei=bc_ids.begin(); ei != bc_ids.end(); ei++){
        
      if(!bc_map_empty){
        int id = *ei ;
        std::map<int, string>::const_iterator mi= bc_map.find(id) ;
        if(mi != bc_map.end()){
          surf_ids.push_back(pair<int,string>(*ei,mi->second)) ;
          
        }else{
          int id = *ei ;
          char bcname[128] ;
	  bzero(bcname,128) ;
          snprintf(bcname,127,"BC_%d",id) ;
          string bcstr(bcname) ;
          surf_ids.push_back(pair<int,string>(*ei,bcstr)) ;
        }
      }else{
        int id = *ei ;
        char bcname[128] ;
	bzero(bcname,128) ;
        snprintf(bcname,127,"BC_%d",id) ;
        string bcstr(bcname) ;
        surf_ids.push_back(pair<int,string>(*ei,bcstr)) ;
      }
    }
  }
 
  if(new_bc_id_used){     
    char bcname[128] ;
    bzero(bcname,128) ;
    snprintf(bcname,127,"BC_%d",new_bc_id) ;
    string bcstr(bcname) ;
    surf_ids.push_back(pair<int,string>(new_bc_id, bcstr)) ; 
  }
  Loci::writeVOG(output_file, act_pos_dat, cl, cr, face2node,surf_ids, volTags) ;
  // Close everything up
  H5Fclose(input_fid) ;
  Loci::Finalize() ;
  return 0 ;
}
    
