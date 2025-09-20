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
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <strings.h>
#include <Loci>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <kd_tree.h>
using namespace Loci::kdTree;

using std::map ;
using std::vector ;
using std::string ;
using std::cout ;
using std::cerr ;
using std::endl ;
using Loci::vector3d ;
using Loci::entitySet;
typedef vector3d<double> vect3d;

FILE *mytmpfile() {
  FILE *file = 0 ;
  char filename[128] ;
  pid_t pid = getpid() ;
  
  bzero(filename,128) ;
  snprintf(filename,127,"vogmerge.tmp.%d",pid) ;
  file = fopen(filename,"w+") ;
  if(!file) {
    return file ;
  }
  remove(filename) ;
  return file ;
}
    
  
  
const double NORMALIZE_ZERO_THRESHOLD = 0.0 ;
inline  void normalize(vect3d& v) {
  if( (fabs(v.x) <= NORMALIZE_ZERO_THRESHOLD) &&
      (fabs(v.y) <= NORMALIZE_ZERO_THRESHOLD) &&
      (fabs(v.z) <= NORMALIZE_ZERO_THRESHOLD)) {
    cerr << "WARNING: normalizing zero vector, nothing done!" << endl ;
    return ;
  }
  double t = sqrt(v.x*v.x + v.y*v.y + v.z*v.z) ;
  v.x /= t ;
  v.y /= t ;
  v.z /= t ;
}


struct file_info {
  //parameters in the original grid file
  unsigned long numNodes ; 
  unsigned long numFaces ; 
  unsigned long numCells ;
  //actual values after nodes and faces are deleted
  unsigned long numActualNodes; 
  unsigned long numActualFaces;
  //computed from actual values
  unsigned long nodeOffset; 
  unsigned long cellOffset; 
  vector<short> cluster_sizes ;
  //number of nodes/faces on the boundaries that will be glued 
  unsigned long numBdNodes;
  unsigned long numBdFaces;
} ;

unsigned long readAttributeLong(hid_t group, const char *name) {
  hid_t id_a = H5Aopen_name(group,name) ;
  unsigned long val = 0;
  H5Aread(id_a,H5T_NATIVE_ULONG,&val) ;
  H5Aclose(id_a) ;
  return val ;
}

struct affineMapping {
  double M[4][4] ;
  double determinant() { return
      (M[0][0]*M[1][1]*M[2][2]+
       M[1][0]*M[2][1]*M[0][2]+
       M[2][0]*M[0][1]*M[1][2]) -
      (M[0][0]*M[2][1]*M[1][2]+
       M[1][0]*M[0][1]*M[2][2]+
       M[2][0]*M[1][1]*M[0][2]) ;
  }

  bool leftHanded() {return (determinant() < 0) ; }
  affineMapping() {
    for(int i=0;i<4;++i) {
      for(int j=0;j<4;++j)
        M[i][j] = 0 ;
      M[i][i] = 1 ;
    }
  }
  void Combine(affineMapping a) {
    double Mtmp[4][4] ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        Mtmp[i][j] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j) {
        double mtmp = 0 ;
        for(int k=0;k<4;++k)
          mtmp += a.M[i][k]*M[k][j] ;
        Mtmp[i][j] = mtmp ;
      }
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        M[i][j] = Mtmp[i][j] ;
  }
  void translate(vect3d tv) {
    affineMapping tmp ;
    tmp.M[0][3] = tv.x ;
    tmp.M[1][3] = tv.y ;
    tmp.M[2][3] = tv.z ;
    Combine(tmp) ;
  }
  void scale(vect3d tv) {
    affineMapping tmp ;
    tmp.M[0][0] = tv.x ;
    tmp.M[1][1] = tv.y ;
    tmp.M[2][2] = tv.z ;
    Combine(tmp) ;
  }
  void rotateX(double theta) {
    double th = -theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[1][1] =  cth ;
    tmp.M[1][2] =  sth ;
    tmp.M[2][1] = -sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateY(double theta) {
    double th = -theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][2] = -sth ;
    tmp.M[2][0] =  sth ;
    tmp.M[2][2] =  cth ;
    Combine(tmp) ;
  }
  void rotateZ(double theta) {
    double th = -theta*2.*M_PI/360. ;
    double sth = sin(th) ;
    double cth = cos(th) ;
    affineMapping tmp ;

    tmp.M[0][0] =  cth ;
    tmp.M[0][1] =  sth ;
    tmp.M[1][0] = -sth ;
    tmp.M[1][1] =  cth ;
    Combine(tmp) ;
  }
  //mirror with a 3d plane with origin p and  normal n  
  void mirror(vect3d p, vect3d n){
   
    double k = dot(p, n);
    affineMapping tmp ;
    tmp.M[0][0] = 1-2.*n.x*n.x;
    tmp.M[1][1] = 1-2.*n.y*n.y;
    tmp.M[2][2] = 1-2.*n.z*n.z;
    
    tmp.M[0][1]  = tmp.M[1][0] = -2.*n.x*n.y;
    tmp.M[0][2] = tmp.M[2][0] = -2.*n.x*n.z;
    tmp.M[0][3] = 2.*n.x*k;
    tmp.M[3][0] = 0;

    tmp.M[1][2] = tmp.M[2][1] = -2.*n.y*n.z;
    tmp.M[1][3] = 2.*n.y*k;
    tmp.M[2][3] = 2.*n.z*k;
    Combine(tmp) ;
  }
    
    
  vect3d Map(vect3d v) {
    double tmp[4] ;
    tmp[0] = v.x ;
    tmp[1] = v.y ;
    tmp[2] = v.z ;
    tmp[3] = 1. ;
    double res[4] ;
    for(int i=0;i<4;++i)
      res[i] = 0 ;
    for(int i=0;i<4;++i)
      for(int j=0;j<4;++j)
        res[i] += M[i][j]*tmp[j] ;
    vect3d r(res[0],res[1],res[2]) ;
    return r ;
  }
   
} ;



namespace Loci {
  extern int getClusterNumFaces(unsigned char *cluster) ;
  extern entitySet faceCluster(const multiMap &face2node,
                               const Map &cl, const Map &cr, entitySet faces,
                               vector<unsigned char> &cluster_info,
                               vector<unsigned short> &cluster_sizes);
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



struct node_info{//the info of the nodes on the boundaries that need to be glued together
  int block_id; //which grid
  int local_index; //the index in this grid file
  int global_index; // the index in glued grid file
  int mapped; //0: initialized, 1: mapped, -1: deleted
  node_info(int i, int j):block_id(i), local_index(j), global_index(-1), mapped(-1){};
};

// Information about boundary faces.  Only faces on glued boundaries are stored.
struct face_info{
  int block_id; //which grid
  int face_id; // the local index in this grid file
  int cl;
  int cr; 
  int mapped; //0: initialized, 1: mapped, -1: deleted
  face_info(int i, int j, int k, int l):block_id(i), face_id(j), cl(k), cr(l),mapped(0){};
};

//read in the nodes and faces on the boundaries that need to be glued
bool readVolBC(int fi,//file id, i.e., block_id
               hid_t input_fid, 
               const vector<int>& bc_id, // the ids of all boundaries that need to be glued  
               affineMapping& gridXform , 
               bool applyMap,
               vector<node_info>& node_map, //info of all nodes on the boundaries that need gluing
               vector<coord3d>& pos, //pos of all nodes on the boundaries that need gluing
	       vector<double> &pos_len, // edge length per point
               double& min_len, //minimum edge length
               vector<vector<int> >& f2n, //f2n of all faces on the boundaries that need gluing
               vector<face_info>& face_map, // info of all faces on the boundaries that need gluing
               vector<file_info>& fileData 
               ) {
  
  int f2n_start = f2n.size();
  
    
  // read in positions
  hid_t fid = H5Gopen(input_fid,"file_info",H5P_DEFAULT) ;
  unsigned long numNodes = readAttributeLong(fid,"numNodes") ;
  H5Gclose(fid) ;
  
  hsize_t  count = numNodes ;

  hsize_t lstart = 0 ;

  // Read in pos data from file i
  vector<vect3d> pos_dat(numNodes) ;
  vector<double> node_len(numNodes,1e30) ;

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
  

  


  hid_t face_g = H5Gopen(input_fid,"face_info",H5P_DEFAULT) ;
  // Read cluster sizes
  dataset = H5Dopen(face_g,"cluster_sizes",H5P_DEFAULT) ;
  
  dspace = H5Dget_space(dataset) ;
  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;
  vector<unsigned short> csizes(size) ;
  dimension = size ;
 
  hsize_t start = 0 ;

  count = size ;
  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
  rank = 1 ;
  memspace = H5Screate_simple(rank,&dimension,NULL) ;
  err = H5Dread(dataset,H5T_NATIVE_USHORT,memspace,dspace,H5P_DEFAULT,&csizes[0]) ;
  if(err < 0) {
    cerr << "unable to read cluster sizes from file" << endl ;
    exit(-1) ;
  }
  H5Dclose(dataset) ;
  H5Sclose(memspace) ;

  // Read in clusters and transform
  dataset = H5Dopen(face_g,"cluster_info",H5P_DEFAULT) ;
  dspace = H5Dget_space(dataset) ;
  start = 0 ;

  
  

  entitySet nodeSet;
  unsigned long face_count = 0;
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
    Loci::multiMap face2node ;
    Loci::Map cl,cr ;
    face2node.allocate(fcnts) ;
    cl.allocate(fclust) ;
    cr.allocate(fclust) ;
    Loci::fillFaceInfo(&cluster[0],face2node,cl,cr,0);
    if(gridXform.leftHanded()) {
      // If coordinate system changes handedness, we need to flip face
      // ordering to get correct normal orientation
      for(int f=0;f<nfaces;++f) {
        int fsz = face2node[f].size() ;
        for(int n=0;n<fsz/2;++n) {
          std::swap(face2node[f][n],face2node[f][fsz-1-n]) ;
        }
      }
    }
       

    
    // Now loop over faces in cluster to determine if any are boundary 
    // faces
    for(int f=0;f<nfaces;++f) {
      
      if(find(bc_id.begin(), bc_id.end(), -cr[f])!= bc_id.end()) { // boundary face 
        int fsz = face2node[f].size() ;
        vector<int> tmpVec(fsz);
        for(int i =0; i < fsz; i++){
          nodeSet += face2node[f][i];
          tmpVec[i] = face2node[f][i];
	  int n1 = face2node[f][i] ;
	  int n2 = face2node[f][(i+1)%fsz] ;
          double edgelen = norm(pos_dat[n1]- pos_dat[n2]);
	  node_len[n1] = min(node_len[n1],edgelen) ;
	  node_len[n2] = min(node_len[n2],edgelen) ;
          min_len = min(min_len,
                        edgelen);
          
        }
        f2n.push_back(tmpVec);
        face_map.push_back(face_info(fi, face_count,  cl[f], cr[f]));
      }
      face_count++;
    }//end loop over faces
  }//end loop over cluster

  //update fileData
  int f2n_end = f2n.size();
  fileData[fi].numBdNodes=nodeSet.size();
  fileData[fi].numBdFaces=f2n.size()-f2n_start;
  
  // Now mark all nodes that the faces access
  vector<int> used(numNodes, 0) ;
   
  // Count nodes in the surface mesh
  int nnodes = 0 ;
  for(entitySet::const_iterator ei = nodeSet.begin(); ei != nodeSet.end(); ei++){
    used[*ei] =1+ nnodes++ ;
  }
  

  //copy used nodes to pos, and index them
  for(size_t i=0;i<numNodes;++i){
    if(used[i] != 0){
      vect3d p = pos_dat[i];
      if(applyMap) p = gridXform.Map(p);
      node_map.push_back(node_info(fi, i));
      pos.push_back(coord3d(p.x, p.y, p.z));
      pos_len.push_back(node_len[i]) ;
      used[i] = pos.size();
    }
  }

  //switch f2n to index in pos
  for(int i = f2n_start; i < f2n_end; i++){
    for(size_t j = 0; j <f2n[i].size(); j++){
      f2n[i][j] = used[f2n[i][j]]-1;//start with 0
    }
  }
 
  
  
#ifdef DEBUG
  cout << "nnodes = "<< nnodes << endl ;
#endif
  return true; 
}



    

void createNodeMap(const vector<coord3d>& pos,
                   const vector<double> &min_len,
                   std::map<int, vector<int> >& nodeID)
{
  int sz =pos.size();

  // Sparse matrix using entitySets  
  vector<entitySet> matrix(sz)  ;
  
  //each point get an id
  vector<int> pntid(sz);
  for(int i=0;i<sz;++i)pntid[i] = i ;

  //build search tree
  Loci::kdTree::kd_tree search_tree(pos,pntid);
  
  for(int i = 0; i < sz; i++){ 
    
    vect3d center = vect3d(pos[i][0],pos[i][1], pos[i][2]);
    kd_tree::bounds box;
    box.maxc[0] = center.x + min_len[i];
    box.maxc[1] = center.y+ min_len[i];
    box.maxc[2] = center.z + min_len[i];
    
    box.minc[0] = center.x - min_len[i];
    box.minc[1] = center.y - min_len[i];
    box.minc[2] = center.z - min_len[i];
    
    
    // Find all of the points within a given bounding box
    std::vector< kd_tree::coord_info> found_pts;
    search_tree.find_box(found_pts,  box);
      
    for(unsigned int pi = 0; pi < found_pts.size(); pi++){
      vect3d p = vect3d(found_pts[pi].coords[0], found_pts[pi].coords[1], found_pts[pi].coords[2]);
      if(norm(p - center) <= min_len[i] && i != found_pts[pi].id){
	matrix[i] += found_pts[pi].id ;
	matrix[found_pts[pi].id] += i  ;
      }   
    }
    
  }

  // transitive closure algorithm here 
  bool updated = true;
  while(updated){
    updated = false;
    for(int i = 0; i < sz; i++) {
      for(entitySet::const_iterator ei=matrix[i].begin();
	  ei != matrix[i].end(); ++ei) {
	int j = *ei ;
	for(entitySet::const_iterator ek=matrix[j].begin();
	    ek!= matrix[j].end();++ek) {
	  int k = *ek ;
	  if(i != k) {
	    if(!matrix[i].inSet(k)){
	      updated = true;
	    }
	    matrix[i] += k ;
	  }
	}
      }
    }
  }
  
  for(int i = 0; i < sz; i++){
    vector<int> eqPoints;
    for(entitySet::const_iterator ei=matrix[i].begin();
	ei != matrix[i].end(); ++ei) {
      int j = *ei ;
      eqPoints.push_back(j);
    }
    if(eqPoints.size() >0) nodeID[i] = eqPoints;
  }
}

//find deleted nodes,set the global index in nodeMap, and update fileData
void setNodeIndex(const std::map<int, vector<int> >& nodeID,
                  vector<node_info>& nodeMap,
                  vector<file_info>& fileData){

  //first compute numDelNodes
  int num_inputs = fileData.size();
 
  vector<unsigned long> numDelNodes(num_inputs);
  for(int i = 0; i < num_inputs; i++){
    numDelNodes[i] = 0;
  }
  for(std::map<int, vector<int> >::const_iterator ei = nodeID.begin();
      ei != nodeID.end(); ei++){
    if(ei->second.size()> 0){
      if(nodeMap[ei->first].block_id > nodeMap[ei->second[0]].block_id){
        numDelNodes[nodeMap[ei->first].block_id]++;
      }else if(nodeMap[ei->first].block_id == nodeMap[ei->second[0]].block_id){
        if(nodeMap[ei->first].local_index > nodeMap[ei->second[0]].local_index){
          numDelNodes[nodeMap[ei->first].block_id]++;
        }
      }
    }
  }
  for(int i = 0; i < num_inputs; i++){
    fileData[i].numActualNodes = fileData[i].numNodes - numDelNodes[i];
  }
  
  //compute nodeOffset in fileData
  unsigned long noffset = 0;
  for(int i =0; i < num_inputs; i++){
    fileData[i].nodeOffset = noffset;
    noffset += fileData[i].numActualNodes;
  }
  
  //compute global index
  int numUnmatchedNodes = 0;
  unsigned long nodemap_offset = 0;
  for(int fi = 0; fi < num_inputs; fi++){
    //mark deleted points in nodeMap and in used
    vector<int> used(fileData[fi].numNodes, 0); 
   
    for(unsigned int i = nodemap_offset; i < nodemap_offset+ fileData[fi].numBdNodes; i++){
      std::map<int, vector<int> >::const_iterator pnt = nodeID.find(i);
      if(pnt == nodeID.end() || (pnt->second).size()==0){//this node is not matched
        nodeMap[i].mapped = 1;
        numUnmatchedNodes++;
      }else if(nodeMap[i].block_id < nodeMap[(pnt->second)[0]].block_id){//this node need to be kept
        nodeMap[i].mapped = 1; 
      }else if(nodeMap[i].block_id == nodeMap[(pnt->second)[0]].block_id){
        if(nodeMap[pnt->first].local_index < nodeMap[(pnt->second)[0]].local_index) nodeMap[i].mapped = 1; 
        else{//this node need to be deleted
          nodeMap[i].mapped = -1;
          used[nodeMap[i].local_index] = -1;
        }
      }else{ //this node need to be deleted
        nodeMap[i].mapped = -1;
        used[nodeMap[i].local_index] = -1;
      }
    }
    
    //index the  points that are not deleted, start from 1
    unsigned long count = 1;
    for(unsigned int i = 0; i < used.size(); i++){
      if(used[i] != -1) used[i] = count++;
    }

   
    
    //set the global index in nodeMap 
    for(unsigned int i = nodemap_offset; i < nodemap_offset+ fileData[fi].numBdNodes; i++){
      
      if(nodeMap[i].mapped ==1)
        nodeMap[i].global_index = used[nodeMap[i].local_index]-1 + fileData[nodeMap[i].block_id].nodeOffset;
      else{
        std::map<int, vector<int> >::const_iterator pnt = nodeID.find(i);
        nodeMap[i].global_index = nodeMap[(pnt->second)[0]].global_index;
      }
    }
    nodemap_offset += fileData[fi].numBdNodes;
  }

  
  if(numUnmatchedNodes > 0){
    cerr<< "Warning: number of unmatched nodes on glued boundaries: " << numUnmatchedNodes << endl;
  }
}
      

//if two faces are equivalent                  
bool is_equivalent(const vector<int>& f2n_1,
                   const vector<int>& f2n_2){
  
  int size = f2n_1.size();
  if(size != (int)f2n_2.size() || size==0)return false;
  int offset = -1;
  for(int i = 0; i < size; i++){
    if(f2n_2[i] == f2n_1[0]){
      offset = i;
      break;
    }
  }
  if(offset == -1)return false;
  for(int i =0; i < size; i++){
    if(f2n_1[i] != f2n_2[(offset+size-i)%size])return false;
  }
  return true;
}
 
// match faces on glued boundaries,
//find deleted faces,
//set cl,cr, and mapped in faceMap
//update fileData
void createFaceMap(const vector<vector<int> >& face2node,
                   vector<face_info>& faceMap,
                   const Loci::multiMap& node2face,
                   const vector<node_info>& node_map,
                   const std::map<int, vector<int> >& nodeID,
                   vector<file_info>& fileData){
  
  //initialize number of deleted faces in each block
  int num_inputs = fileData.size();
  vector<int> numDelFaces(num_inputs);
  for(int i = 0; i < num_inputs; i++){
    numDelFaces[i] = 0;
  }
  
  //update cellOffset in fileData
  unsigned long coffset  = 0;
  for(int i = 0; i < num_inputs; i++){
    fileData[i].cellOffset = coffset;
    coffset += fileData[i].numCells;
  }
  
  //first build a face map faceID, between temporary global index 
  map<int, int> faceID;
  
  //for each face f in face2node
  for(unsigned int f = 0; f <faceMap.size(); f++){
    if(faceMap[f].mapped == 0){ //if not mapped yet
    
      //get its f2n in global indexes
      vector<int> f2n_1(face2node[f].size());
      for(unsigned int n=0; n < face2node[f].size(); n++){
        f2n_1[n] = node_map[face2node[f][n]].global_index;
      }

      
      //pick the first node of face f
      //for each face fi that is associated with the equivalence of first node
      int theNode = face2node[f][0];
      std::map<int, vector<int> >:: const_iterator pnt = nodeID.find(theNode);
      if(pnt==nodeID.end()){
        faceMap[f].mapped = 1;
        faceMap[f].cl =  faceMap[f].cl + fileData[faceMap[f].block_id].cellOffset;
        continue; //no matched nodes
      }
      vector<int>  firstNodes = vector<int>(pnt->second);
      bool found = false;
      for(unsigned int ni = 0; ni < firstNodes.size(); ni++){
       
        for(int fi = 0; fi < node2face.num_elems(firstNodes[ni]); fi++){
          int theface = node2face[firstNodes[ni]][fi];
          if(faceMap[theface].mapped==0 && face2node[theface].size()==face2node[f].size()){
            vector<int> f2n_2(face2node[f].size());
            for(unsigned int n=0; n < face2node[theface].size(); n++){
              f2n_2[n] = node_map[face2node[theface][n]].global_index;
            }
            found = is_equivalent(f2n_1, f2n_2);
            if(found){
              faceID[f]=theface;
              faceID[theface] = f; //mutual
              if(faceMap[f].block_id < faceMap[theface].block_id){
                faceMap[f].mapped = 1;
                faceMap[theface].mapped = -1;
                numDelFaces[faceMap[theface].block_id]++;
                faceMap[f].cl =  faceMap[f].cl + fileData[faceMap[f].block_id].cellOffset;
                faceMap[f].cr =  faceMap[theface].cl + fileData[faceMap[theface].block_id].cellOffset;
              }else if(faceMap[f].block_id == faceMap[theface].block_id){
                if(faceMap[f].face_id < faceMap[theface].face_id){
                  faceMap[f].mapped = 1;
		  faceMap[theface].mapped = -1;
		  numDelFaces[faceMap[theface].block_id]++;
		  faceMap[f].cl =  faceMap[f].cl + fileData[faceMap[f].block_id].cellOffset;
		  faceMap[f].cr =  faceMap[theface].cl + fileData[faceMap[theface].block_id].cellOffset;
                }else{
                  faceMap[theface].mapped = 1;
                  faceMap[f].mapped = -1;
                  numDelFaces[faceMap[f].block_id]++;
                  faceMap[theface].cl =  faceMap[theface].cl + fileData[faceMap[theface].block_id].cellOffset;
                  faceMap[theface].cr =  faceMap[f].cl + fileData[faceMap[f].block_id].cellOffset;
                }
                
                
              }else{

                faceMap[theface].mapped = 1;
                faceMap[f].mapped = -1;
                numDelFaces[faceMap[f].block_id]++;
                faceMap[theface].cl =  faceMap[theface].cl + fileData[faceMap[theface].block_id].cellOffset;
                faceMap[theface].cr =  faceMap[f].cl + fileData[faceMap[f].block_id].cellOffset;
              }
              break;
            }
          }
        }
        if(found)break;
      }
      if(faceMap[f].mapped==0){//the face has no equivalence
        faceMap[f].mapped = 1;
        faceMap[f].cl =  faceMap[f].cl + fileData[faceMap[f].block_id].cellOffset;
        
      }
    }
    
  }
  //update fileData
  for(int i = 0; i < num_inputs; i++)fileData[i].numActualFaces = fileData[i].numFaces - numDelFaces[i];
}
  
  


                
void inverseF2N( const vector<vector<int> >& face2node,
                 int numBdNodes,
                 Loci::multiMap& node2face){
  using namespace Loci;
  Loci::entitySet faces = Loci::interval(0, face2node.size()-1);
  Loci::store<int> fcount;
  fcount.allocate(faces);
  for(size_t i = 0; i < face2node.size(); i++){
    fcount[i] = face2node[i].size();
  }
  Loci::multiMap f2n;
  f2n.allocate(fcount);
  for(size_t i = 0; i < face2node.size(); i++){
    for( int j = 0; j < fcount[i]; j++){
      f2n[i][j] = face2node[i][j];
    }
  }
  Loci::entitySet image = Loci::interval(0, numBdNodes-1);
  Loci::entitySet preimage = f2n.domain();
  Loci::inverseMap(node2face, f2n, image, preimage);
}



//read in the nodes and faces on the boundaries that need to be glued
//find matched nodes and faces
//set global indexes for  boudnary nodes
//set cl and cr for boundary faces
//update fileData
void create_equivalence(const vector<hid_t>& input_fids,
                        const vector<vector<int> >& glue_info,
                        vector<affineMapping>& gridXforms,
                        vector<file_info>& fileData,
                        const vector<bool>& applyMap,
                        double tol,
                        vector<node_info>& nodeMap, //info of all nodes on the boundaries that need to be glued
                        vector<face_info>& faceMap//info of all faces on the boundaries that need to be glued
                        ){
  
  double   rmin =  std::numeric_limits<float>::max() ;
 
  vector<coord3d> pos; //coord of all nodes on the boundaries that need to be glued
  vector<double> pos_len ; // edge length associated with each pos
 
  //f2n of all faces on the boundaries that need to be glued, indexes of nodes are temperary global, i.e.,
  //the indexes in pos
  vector<vector<int> > f2n; 
  for(unsigned int fi = 0; fi < input_fids.size(); fi++){
    readVolBC(fi,
              input_fids[fi] ,
              glue_info[fi],
              gridXforms[fi],
              applyMap[fi],
              nodeMap, 
              pos,
	      pos_len,
              rmin,
              f2n,
              faceMap,
              fileData
              ) ;                   
  }
  
  double factor = 0.1 ;
  if(tol > 0 && tol < 1.0)
    factor = tol ;
  for(size_t i=0;i<pos_len.size();++i)
    pos_len[i] *= factor ;

  map<int, vector<int> > nodeID;
  createNodeMap(pos,
                pos_len,
                nodeID);
  setNodeIndex(nodeID, nodeMap, fileData);
 
  Loci::multiMap node2face;
  int totalNumBdNodes = pos.size();
  inverseF2N(f2n,
             totalNumBdNodes,
             node2face);

  createFaceMap(f2n,
                faceMap,
                node2face,
                nodeMap,
                nodeID,
                fileData);

  
  
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


void splitStr(string str, vector<string>& args){
 
  size_t found = str.find_first_of(',');
  size_t start = 0;
  if(found == string::npos){
    args.push_back(str);
    return;
  }else{
    while (found!=string::npos){
     
      args.push_back(str.substr(start, found-start));
      start = found+1;
      found=str.find_first_of(',', start);
    }
    args.push_back(str.substr(start)); 
  }
}




void Usage() {
  cout << "******************************************************************************" << endl ;

  cout << "vogmerge is a utility for merging/gluing vog grids into a single vog file." << endl
       << endl ;
  cout << "Usage: " << endl
       << "  vogglue <options>" << endl
       << "where <options> may be given by:" << endl 
       <<  endl ;
  cout << "  -g <gridfile>           | Specify input grid filename" << endl
       << "  -o <gridfile>           | Specify output grid filename" << endl
       << "  -xshift <value>         | translate grid in x-coordinate" << endl
       << "  -yshift <value>         | translate grid in y-coordinate" << endl
       << "  -zshift <value>         | translate grid in z-coordinate" << endl
       << "  -xscale <value>         | scale grid in x-coordinate" << endl
       << "  -yscale <value>         | scale grid in y-coordinate" << endl
       << "  -zscale <value>         | scale grid in z-coordinate" << endl
       << "  -scale  <value>         | scale grid in all coordinates" << endl
       << "  -xrotate <value>        | rotate scale about x-axis (degrees)" << endl
       << "  -yrotate <value>        | rotate scale about y-axis (degrees)" << endl
       << "  -zrotate <value>        | rotate scale about z-axis (degrees)" << endl
       << "  -bc <oldname>,<newname> | rename boundary surface" << endl
       << "  -tag <name>             | specify volume tag for input grid" << endl
       << "  -mirror <v1>,...,<v6>   | mirror the grid about a 3d plane " << endl
       << "                          | with origin (v1, v2, v3) and normal (v4, v5, v6) "<< endl
       << "  -mirrorx                | mirror the grid about x=0 plane " <<endl
       << "  -mirrory                | mirror the grid about y=0 plane " <<endl
       << "  -mirrorz                | mirror the grid about z=0 plane " <<endl
       << "  -tol <value>            | specify the tolerance for matching points, if not "<<endl
       << "                          | specified, the 10% of min_edge_len on glued boundares is used" << endl
       << "  -glue <name1>,...       | specify the names of boundaries need to be glued" <<endl
       << endl ;
  cout << "******************************************************************************" << endl ;
}
  
int main(int ac, char *av[]) {
  using Loci::entitySet ;
  using Loci::vector3d ;
  Loci::Init(&ac, &av) ;
  if(Loci::MPI_processes > 1) {
    cerr << "vogmerge is not parallel! Run on only one processor!" << endl ;
    Loci::Abort() ;
  }
  string output_file = "glued.vog" ;
  vector<string> input_files ;
  vector<map<string,string> > bc_rename ;
  vector<vector<string> > bc_glue;
  vector<vector<pair<string,entitySet> > > volTag ;
  vector<bool> applyMap ;
  vector<affineMapping> gridXform ;
  double tol = -1.0;
  

  for(int i=1;i<ac;++i) {
    string opt = av[i] ;
    bool parsed = false ;
    if(opt == "-g") {
      i++ ;
      if(i<ac) {
        input_files.push_back(string(av[i])) ;
        bc_rename.push_back(map<string,string>()) ;
        bc_glue.push_back(vector<string>());
        volTag.push_back(vector<pair<string,entitySet> >()) ;
        applyMap.push_back(false) ;
        gridXform.push_back(affineMapping()) ;
        parsed = true ;
      }
    }
    if(input_files.size() > 0) {
      if(opt == "-xshift") {
        i++ ;
        if(i<ac) {
          vect3d tv(atof(av[i]),0,0) ;
          gridXform.back().translate(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-yshift") {
        i++ ;
        if(i<ac) {
          vect3d tv(0,atof(av[i]),0) ;
          gridXform.back().translate(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-zshift") {
        i++ ;
        if(i<ac) {
          vect3d tv(0,0,atof(av[i])) ;
          gridXform.back().translate(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-xrotate") {
        i++ ;
        if(i<ac) {
          gridXform.back().rotateX(atof(av[i])) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-yrotate") {
        i++ ;
        if(i<ac) {
          gridXform.back().rotateY(atof(av[i])) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-zrotate") {
        i++ ;
        if(i<ac) {
          gridXform.back().rotateZ(atof(av[i])) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-xscale") {
        i++ ;
        if(i<ac) {
          vect3d tv(atof(av[i]),1.,1.) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-yscale") {
        i++ ;
        if(i<ac) {
          vect3d tv(1.,atof(av[i]),1.) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-zscale") {
        i++ ;
        if(i<ac) {
          vect3d tv(1.,1.,atof(av[i])) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }
      if(opt == "-scale") {
        i++ ;
        if(i<ac) {
          double v = atof(av[i]) ;
          vect3d tv(v,v,v) ;
          gridXform.back().scale(tv) ;
          applyMap.back() = true ;
          parsed = true ;
        }
      }      
      if(opt == "-tol") {
        i++ ;
        if(i<ac) {
          tol = atof(av[i]) ;
          parsed = true ;
        }
      }      
      if(opt== "-bc") {
        i++ ;
        if(i<ac) {
          char *p = index(av[i],',') ;
          if(p==0) {
            cerr << "-bc option, comma missing in argument" << endl ;
            Usage() ;
            exit(-1) ;
          }
          *p = '\0' ;
          
          string first(av[i]), second(p+1) ;
          bc_rename.back()[first]=second ;
          parsed = true ;
        }
      }

      if(opt== "-glue") {
        i++ ;
        if(i<ac) {
          splitStr(av[i], bc_glue.back());
          parsed = true ;
        }
      }

      if(opt== "-mirror") {
        i++ ;
        if(i<ac) {
         
          vector<string> PN;
          splitStr(av[i], PN);
          if(PN.size()!=6){
            cerr << "ERROR: -mirror '" << av[i] <<"' is invalid!" << endl 
                 << "****** should be followed by 6 double values separated by comma to specify the " << endl
                 << "****** point and normal of the plane" << endl; 
          }else{
            vect3d p = vect3d(atof(PN[0].c_str()), atof(PN[1].c_str()), atof(PN[2].c_str()));
            vect3d n = vect3d(atof(PN[3].c_str()), atof(PN[4].c_str()), atof(PN[5].c_str()));
            normalize(n);
            gridXform.back().mirror(p, n) ;
            applyMap.back() = true ;
            parsed = true;
          }
        }
      }

      if(opt== "-mirrorx") {
        vect3d p = vect3d(0, 0, 0);
        vect3d n = vect3d(1, 0, 0);
        gridXform.back().mirror(p, n) ;
        applyMap.back() = true ;
        parsed = true;
      }
    
      if(opt== "-mirrory") {
        vect3d p = vect3d(0, 0, 0);
        vect3d n = vect3d(0, 1, 0);
        gridXform.back().mirror(p, n) ;
        applyMap.back() = true ;
        parsed = true;
      }
    
       
      if(opt== "-mirrorz") {
        vect3d p = vect3d(0, 0, 0);
        vect3d n = vect3d(0, 0, 1);
        gridXform.back().mirror(p, n) ;
        applyMap.back() = true ;
        parsed = true;
      }
        
      if(opt=="-tag") {
        i++ ;
        if(i<ac) {
          vector<pair<string,entitySet> > vp ;
          string tagname = av[i] ;
          bool validtag = true ;
          if(tagname.size() == 0)
            validtag = false ;
          if(validtag & !isalpha(tagname[0]))
            validtag = false ;
          if(validtag) {
            for(size_t i=1;i<tagname.size();++i)
              if(!isalnum(tagname[i]) && tagname[i] != '_')
                validtag = false ;
          }
          if(validtag) {
            vp.push_back(pair<string,entitySet>(tagname,Loci::EMPTY)) ;
            volTag.back() = vp ;
            parsed = true ;
          } else {
            cerr << "ERROR: -tag '" << tagname <<"' is invalid!" << endl 
                 << "****** Volume tag names must start with a letter and must only contain " << endl 
                 << "****** characters that are letters, numbers, or underscore." << endl;
          }
        }
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

  int num_inputs = input_files.size() ;
  if(num_inputs == 0) {
    cerr << "no input grid files to process!" << endl ;
    Usage() ;
    exit(-1) ;
  }


  vector<hid_t> input_fid(num_inputs) ;

  bool fail = false ;
  
  hid_t output_fid = 0 ;
  output_fid = H5Fcreate(output_file.c_str(),
                         H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT) ;
  if(output_fid <= 0) {
    cerr << "unable to open output file '" << output_file << "'" << endl ;
    fail = true ;
  }
#define DEBUG_TEST
#ifndef DEBUG_TEST
  /* Save old error handler */
  herr_t (*old_func)(void*) = 0;
  void *old_client_data = 0 ;
  H5Eget_auto(&old_func, &old_client_data);
  
  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);
#endif


  
  for(int i=0;i<num_inputs;++i) 
    input_fid[i] = H5Fopen(input_files[i].c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  for(int i=0;i<num_inputs;++i)  
    if(input_fid[i] <= 0) {
      cerr << "unable to open file '" << input_files[i] << "'"<< endl ;
      fail = true ;
    }
  if(fail)
    exit(-1) ;
  
  vector<file_info> fileData(num_inputs) ;
  for(int i=0;i<num_inputs;++i) {
    hid_t fi = H5Gopen(input_fid[i],"file_info",H5P_DEFAULT) ;

    fileData[i].numNodes = readAttributeLong(fi,"numNodes") ;
    fileData[i].numFaces = readAttributeLong(fi,"numFaces") ;
    fileData[i].numCells = readAttributeLong(fi,"numCells") ;

   
    H5Gclose(fi) ;
  }
  
  
  // Setup volume tags
  for(int i=0;i<num_inputs;++i) {
    entitySet compSet = Loci::interval(0,fileData[i].numCells-1) ;
    if(volTag[i].size() == 0) {
      vector<pair<string,entitySet> > volDat ;
      readVolTags(input_fid[i],volDat) ;
      for(size_t j=0;j<volDat.size();++j)
        volTag[i].push_back(volDat[j]) ;
    } else {
      if(volTag[i].size() != 1)
        cerr << "Grid '" << input_files[i] << "' has more than one volume tag!" << endl ;
      volTag[i][0].second = compSet ;
    }
  }
  
  
  vector<vector<pair<int,string> > > boundary_ids(num_inputs);
  for(int i=0;i<num_inputs;++i) {
    Loci::readBCfromVOG(input_files[i],boundary_ids[i]) ;
    int sz = boundary_ids[i].size() ;

   
    map<string,string>::const_iterator mi ;
    for(int bc=0;bc<sz;++bc) {
      if((mi = bc_rename[i].find(boundary_ids[i][bc].second)) != bc_rename[i].end())
        boundary_ids[i][bc].second = mi->second ;
    }
  }
  
 
  
  //before merge, don't know how many boundaries left
  map<string,int> newbcs ;
  int bcid = 1 ;
  for(int i=0;i<num_inputs;++i) {
    int bndsi = boundary_ids[i].size() ;
    for(int j=0;j<bndsi;++j) {
      string bnd_name = boundary_ids[i][j].second ;
      if(newbcs.find(bnd_name) == newbcs.end()) {
        newbcs[bnd_name] = bcid++ ;
      }
    }
  }

  
  bool needGlue = false;
  //process glue-surf
  vector<vector<int> > glue_info(num_inputs); 
  for(int i=0;i<num_inputs;++i) {
    if(bc_glue[i].size() > 0){
      needGlue = true;
      glue_info[i].resize(bc_glue[i].size());
      for(unsigned int j = 0; j < glue_info[i].size(); j++){
        int bc_id = -1;
        int sz = boundary_ids[i].size();
        for(int bc=0;bc<sz;++bc) {
          if(boundary_ids[i][bc].second == bc_glue[i][j]){
            bc_id = boundary_ids[i][bc].first;
            break;
	  }
        }
        if(bc_id==-1 && bc_glue[i][j].substr(0,3)=="BC_"){
          bc_id = atoi(bc_glue[i][j].substr(3).c_str());
        }
        
        if(bc_id >= 0)glue_info[i][j]=bc_id;
        else{
          cerr<<"Warning: can not find boundary " << bc_glue[i][j] << " in grid " << input_files[i] << endl;
          
        }
      }
    }
  }

 
  vector<node_info> nodeMap;
  vector<face_info> faceMap;
  

  if(needGlue)create_equivalence(input_fid,
                                 glue_info,
                                 gridXform,
                                 fileData,
                                 applyMap,
                                 tol,
                                 nodeMap,
                                 faceMap);     
  
 
  unsigned long numNodes = 0 ;
  unsigned long numFaces = 0 ;
  unsigned long numCells = 0 ;
  if(needGlue){
    for(int i=0;i<num_inputs;++i) {
      numNodes += fileData[i].numActualNodes ;
      numFaces += fileData[i].numActualFaces ;
      numCells += fileData[i].numCells ;
    }
  }else{
    for(int i=0;i<num_inputs;++i) {
      numNodes += fileData[i].numNodes ;
      numFaces += fileData[i].numFaces ;
      numCells += fileData[i].numCells ;
    }
  }
  //compute nodeOffset and cellOffset
  if(!needGlue) {
    unsigned long noffset = 0;
    unsigned long coffset = 0;
    for(int i=0;i<num_inputs;++i) {
      fileData[i].nodeOffset = noffset;
      fileData[i].cellOffset = coffset;
      noffset += fileData[i].numNodes;
      coffset += fileData[i].numCells;
    }
  }

    
  cout << "numNodes: " << numNodes<< endl;
  cout << "numFaces: " << numFaces<< endl;
  cout << "numCells: " << numCells<< endl;
  // Begin writing output file

  {
    // Write out file info
    hid_t file_id = H5Gcreate(output_fid,"file_info",
			      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

    hsize_t dims = 1 ;
    hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;

    hid_t att_id = H5Acreate(file_id,"numNodes",H5T_STD_I64BE,
                             dataspace_id,H5P_DEFAULT,H5P_DEFAULT) ;

    H5Awrite(att_id,H5T_NATIVE_LLONG,&numNodes) ;
    H5Aclose(att_id) ;

    att_id = H5Acreate(file_id,"numFaces", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
    H5Awrite(att_id,H5T_NATIVE_LLONG,&numFaces) ;
    H5Aclose(att_id) ;

    att_id = H5Acreate(file_id,"numCells", H5T_STD_I64BE,
                       dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;

    H5Awrite(att_id,H5T_NATIVE_LLONG,&numCells) ;
    H5Aclose(att_id) ;
    H5Gclose(file_id) ;
  }
  
  {
    // Write out node information
    hid_t node_id = H5Gcreate(output_fid,"node_info",
			      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

    int rank = 1 ;
    hsize_t dimension = numNodes ;
    hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;


    hsize_t start = 0 ;
    hsize_t lstart = 0 ;
    hsize_t stride = 1 ;
    typedef Loci::data_schema_traits<vect3d > traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    hid_t dataseto = H5Dcreate(node_id,"positions",dp->get_hdf5_type(),
                               dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    // Write out node info
    unsigned long nodeMap_start = 0; //the start index in nodeMap
    for(int i=0;i<num_inputs;++i) {
      
      
      hsize_t count = fileData[i].numNodes ;
      
      // Read in pos data from file i
      vector<vect3d > pos_dat(fileData[i].numNodes) ;

      hid_t node_g = H5Gopen(input_fid[i],"node_info",H5P_DEFAULT) ;
      hid_t dataset = H5Dopen(node_g,"positions",H5P_DEFAULT) ;
      hid_t dspace = H5Dget_space(dataset) ;
      
      H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&lstart,&stride,&count,NULL) ;
      int rank = 1 ;
      hsize_t dimension = count ;
      hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t err = H5Dread(dataset,datatype,memspace,dspace,H5P_DEFAULT,
                          &pos_dat[0]) ;
      if(err < 0) {
        cerr << "unable to read positions from '" << input_files[i] << "'" << endl ;
        exit(-1) ;
      }
      H5Sclose(dspace) ;
      H5Dclose(dataset) ;
      H5Gclose(node_g) ;

      if(applyMap[i]) {
        for(size_t ii=0;ii<pos_dat.size();++ii)
          pos_dat[ii] = gridXform[i].Map(pos_dat[ii]) ;
      }
      if(needGlue){
        //mark the nodes deleted
        vector<int> used(fileData[i].numNodes, 0); 
        for(unsigned ii = 0; ii < fileData[i].numBdNodes; ii++){
          if(nodeMap[nodeMap_start+ii].mapped==-1)used[nodeMap[nodeMap_start+ii].local_index] = -1;
        }
        
        //take out deleted nodes, only write out the nodes left
        count = fileData[i].numActualNodes;
        vector<vect3d > act_pos_dat(fileData[i].numActualNodes) ;
        size_t local_count = 0; 
        for(size_t ii=0;ii<pos_dat.size();++ii){
          if(used[ii]!= -1){
            act_pos_dat[local_count++] = pos_dat[ii] ;
          }
        }
      
        // Now write out this file's part of the positions
     
        H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, &start,&stride,&count,NULL)  ;
      
        dimension = count ;
        rank = 1 ;
      
        memspace = H5Screate_simple(rank,&dimension,NULL) ;
        err = H5Dwrite(dataseto,datatype,memspace,dataspace,
                       H5P_DEFAULT,&act_pos_dat[0]) ;
        if(err < 0) {
          cerr << "unable to write positions to '" << output_file << "'" << endl ;
          exit(-1) ;
        }
      
        start += count ;
      
        H5Sclose(memspace) ;
        nodeMap_start += fileData[i].numBdNodes;
      }else{
        
	// Now write out this file's part of the positions
        H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, &start,&stride,&count,NULL)  ;
        dimension = count ;
        rank = 1 ;
        err = H5Dwrite(dataseto,datatype,memspace,dataspace,
                       H5P_DEFAULT,&pos_dat[0]) ;
        if(err < 0) {
          cerr << "unable to write positions to '" << output_file << "'" << endl ;
          exit(-1) ;
        }
        
        start += count ;
      
        H5Sclose(memspace) ;
      }
    
    

        
    }//for(int i =0; i < num_inputs...
    
    H5Sclose(dataspace) ;
    H5Dclose(dataseto) ;
    H5Gclose(node_id) ;
  }
 
  {
    // Write Face Info
    
    // First read in clusters from each file and write modified
    // clusters to a scratch file
    FILE *scratch = mytmpfile() ;
    
    if(scratch == 0) {
      perror("tmpfile") ;
      cerr << "failed to open tmpfile, check setting of TMPDIR environment variable" << endl ;
      Loci::Abort() ;
    }
    
    // Total output file cluster sizes
    vector<unsigned short> cluster_sizes ;

   
    // Loop over files
    unsigned long nodeMap_start = 0;
    unsigned long faceMap_start = 0;
    for(int i=0;i<num_inputs;++i) {
      
     

      unsigned long face_id_start = 0;
      // Generate bc remap ;
      map<int,int> bcr ;
      int bndsi = boundary_ids[i].size() ;
      for(int j=0;j<bndsi;++j) {
        string bnd_name = boundary_ids[i][j].second ;
        int bnd_tag = boundary_ids[i][j].first ;
        bcr[-bnd_tag] = -newbcs[bnd_name] ;
        
      }        

      
      hid_t face_g = H5Gopen(input_fid[i],"face_info",H5P_DEFAULT) ;
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
          input_files[i] << "'" << endl ;
        exit(-1) ;
      }
      H5Dclose(dataset) ;
      H5Sclose(memspace) ;

      // Read in clusters and transform
      dataset = H5Dopen(face_g,"cluster_info",H5P_DEFAULT) ;
      dspace = H5Dget_space(dataset) ;
      start = 0 ;
      
      if(needGlue){
       
        //first create an index map used for node
        vector<int> used(fileData[i].numNodes, 0);
        {
          for(unsigned long ii = 0; ii < fileData[i].numBdNodes; ii++){
            if(nodeMap[nodeMap_start+ii].mapped==-1)used[nodeMap[nodeMap_start+ii].local_index] = -1;
          }
          unsigned long local_count = 0;
          for(unsigned int ii = 0; ii<used.size(); ii++){
            if(used[ii]!=-1) used[ii] = fileData[i].nodeOffset + local_count++;
          }
          for(unsigned long ii = 0; ii < fileData[i].numBdNodes; ii++){
            if(nodeMap[nodeMap_start+ii].mapped==-1)used[nodeMap[nodeMap_start+ii].local_index] = nodeMap[nodeMap_start+ii].global_index;
          }
        }

       

        
        //create an index map fused for face
        map<int, int> CR;
        vector<int> fused(fileData[i].numFaces, 0);
        {
        
          for(unsigned long ii = 0; ii < fileData[i].numBdFaces; ii++){
            
	    if(faceMap[faceMap_start+ii].mapped==-1)fused[faceMap[faceMap_start+ii].face_id] = -1;
	    else CR[faceMap[faceMap_start+ii].face_id]=faceMap[faceMap_start+ii].cr;
          
          }
          
       
        }
        


        
	for(size_t c=0;c<size;++c) { // Loop over clusters
	  count = csizes[c] ;
	  vector<unsigned char> cluster(count) ;
	  H5Sselect_hyperslab(dspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
	  dimension = count ;
	  memspace = H5Screate_simple(rank,&dimension,NULL) ;
	  err = H5Dread(dataset,H5T_NATIVE_UCHAR,memspace,dspace,H5P_DEFAULT,&cluster[0]) ;
	  if(err < 0) {
	    cerr << "unable to read cluster from file '" << input_files[i]
		 << "'" << endl ;
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
	  if(gridXform[i].leftHanded()) {
	    // If coordinate system changes handedness, we need to flip face
	    // ordering to get correct normal orientation
	    for(int f=0;f<nfaces;++f) {
	      int fsz = face2node[f].size() ;
	      for(int n=0;n<fsz/2;++n) {
		std::swap(face2node[f][n],face2node[f][fsz-1-n]) ;
	      }
	    }
	  }
        
        
	  Loci::entitySet act_fclust;
	  Loci::store<int> act_fcnts ;
	  Loci::multiMap act_face2node ;
	  Loci::Map act_cl,act_cr ;
	  for(int f=0;f<nfaces;++f) {
          
	    if(fused[face_id_start+f] != -1) act_fclust += f;
	  }
       

       
	  if(act_fclust.size() ==0){
	    face_id_start += nfaces;
	    continue;
	  }
	  act_fcnts.allocate(act_fclust);
	  for(entitySet::const_iterator ei = act_fclust.begin(); 
	      ei != act_fclust.end(); ei++){
	    act_fcnts[*ei] = fcnts[*ei];
	  }
	  act_cl.allocate(act_fclust);
	  act_cr.allocate(act_fclust);
	  act_face2node.allocate(act_fcnts);
        
       
	  for(entitySet::const_iterator  f=act_fclust.begin();
	      f != act_fclust.end(); f++) {
	    if(find( glue_info[i].begin(),  glue_info[i].end(), -cr[*f])!=  glue_info[i].end()){//boundary face
            
	      int fsz = face2node[*f].size() ;
	      for(int n=0;n<fsz;++n){
		act_face2node[*f][n] = used[face2node[*f][n]];
	      }
            
	      act_cl[*f] = cl[*f] + fileData[i].cellOffset ;
	      act_cr[*f] = CR[face_id_start+int(*f)];
	      if(act_cr[*f]<0){
		map<int,int>::const_iterator mi ;
		if((mi=bcr.find(cr[*f])) != bcr.end()) {
		  act_cr[*f] = mi->second ;
		} else{
		  int id = -act_cr[*f] ;
		  char bcname[128] ;
		  bzero(bcname,128) ;
		  snprintf(bcname,127,"BC_%d",id) ;
		  string bcstr(bcname) ;
              
		  map<string,string>::const_iterator mis ;
              
		  if((mis = bc_rename[i].find(bcstr)) != bc_rename[i].end())
		    bcstr = mis->second ;
		  if(newbcs.find(bcstr) == newbcs.end())
		    newbcs[bcstr] = bcid++ ;
		  bcr[-id] = -newbcs[bcstr] ;
		  act_cr[*f] = bcr[-id] ;
               
		}
	      }
	    }else{
	      int fsz = face2node[*f].size() ;
	      for(int n=0;n<fsz;++n){
		act_face2node[*f][n] = used[face2node[*f][n]];
             
	      }
           
	      act_cl[*f] = cl[*f] + fileData[i].cellOffset ;
	      if(cr[*f] >=0)
		act_cr[*f] = cr[*f] +  fileData[i].cellOffset ;
	      else {
		map<int,int>::const_iterator mi ;
		if((mi=bcr.find(cr[*f])) != bcr.end()) {
		  act_cr[*f] = mi->second ;
		} else {
		  int id = -cr[*f] ;
		  char bcname[128] ;
		  bzero(bcname,128) ;
		  snprintf(bcname,127,"BC_%d",id) ;
		  string bcstr(bcname) ;

		  map<string,string>::const_iterator mis ;
                
		  if((mis = bc_rename[i].find(bcstr)) != bc_rename[i].end())
		    bcstr = mis->second ;
		  if(newbcs.find(bcstr) == newbcs.end())
		    newbcs[bcstr] = bcid++ ;
		  bcr[-id] = -newbcs[bcstr] ;
		  act_cr[*f] = bcr[cr[*f]] ;
               
		}
              
              
	      }
          
                                     
          
	    }
          
	  }//for(f = ....
       
        
	  entitySet cellSet = act_cl.image(act_fclust) +act_cr.image(act_fclust) ;
	  entitySet nodeSet = Loci::MapRepP(act_face2node.Rep())->image(act_fclust) ;
	  if(cellSet.size() > 256 || nodeSet.size() > 256) {
         
         
	    while(act_fclust.size() != 0) {
	      vector<unsigned char> cluster_out ;
	      entitySet fcluster =Loci::faceCluster(act_face2node,act_cl,act_cr,act_fclust,
						    cluster_out,cluster_sizes) ;
	      size_t swrite = fwrite(&cluster_out[0],cluster_out.size(),1,scratch) ;
	      if(swrite != 1) {
		perror("fwrite") ;
		cerr << "write to temporary file failed"<< endl ;
		fclose(scratch) ;
		Loci::Abort() ;
	      }
           
	      act_fclust -= fcluster ;
           
	    }
         
         

	  }else{
         
	    vector<unsigned char> clusterout =
	      Loci::encode_face_cluster(act_face2node,act_cl,act_cr,act_fclust,nodeSet,cellSet) ;
	    // Write cluster to tmp file
	    size_t swrite = fwrite(&clusterout[0],clusterout.size(),1,scratch) ;
	    if(swrite != 1) {
	      perror("fwrite") ;
	      cerr << "write to temporary file failed"<< endl ;
	      fclose(scratch) ;
	      Loci::Abort() ;
	    }
	    unsigned short clsz = clusterout.size() ;
	    cluster_sizes.push_back(clsz) ;
	  }
	  face_id_start += nfaces;

       
	}//for(size_t c = 0; c < csizes.size()...
      }else{

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
            cerr << "unable to read cluster from file '" << input_files[i]
                 << "'" << endl ;
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
	  if(gridXform[i].leftHanded()) {
	    // If coordinate system changes handedness, we need to flip face
	    // ordering to get correct normal orientation
	    for(int f=0;f<nfaces;++f) {
	      int fsz = face2node[f].size() ;
	      for(int n=0;n<fsz/2;++n) {
		std::swap(face2node[f][n],face2node[f][fsz-1-n]) ;
	      }
	    }
	  }
	  for(int f=0;f<nfaces;++f) {
	    int fsz = face2node[f].size() ;
	    for(int n=0;n<fsz;++n)
	      face2node[f][n] += fileData[i].nodeOffset ;
	    cl[f] += fileData[i].cellOffset ;
	    if(cr[f] >=0)
	      cr[f] += fileData[i].cellOffset ;
	    else {
	      map<int,int>::const_iterator mi ;
	      if((mi=bcr.find(cr[f])) != bcr.end()) {
		cr[f] = mi->second ;
	      } else {
		int id = -cr[f] ;
		char bcname[128] ;
		bzero(bcname,128) ;
		snprintf(bcname,127,"BC_%d",id) ;
		string bcstr(bcname) ;

		map<string,string>::const_iterator mis ;
              
		if((mis = bc_rename[i].find(bcstr)) != bc_rename[i].end())
		  bcstr = mis->second ;
		if(newbcs.find(bcstr) == newbcs.end())
		  newbcs[bcstr] = bcid++ ;
		bcr[-id] = -newbcs[bcstr] ;
		cr[f] = bcr[cr[f]] ;
	      }
	    }
	  }
	  entitySet cellSet = cl.image(fclust) +cr.image(fclust) ;
	  entitySet nodeSet = Loci::MapRepP(face2node.Rep())->image(fclust) ;
	  if(cellSet.size() > 256 || nodeSet.size() > 256) {
	    cerr << "problem encoding cluster" << endl ;
	  }

	  vector<unsigned char> clusterout =
	    Loci::encode_face_cluster(face2node,cl,cr,fclust,nodeSet,cellSet) ;
	  // Write cluster to tmp file
	  size_t swrite = fwrite(&clusterout[0],clusterout.size(),1,scratch) ;
	  if(swrite != 1) {
	    perror("fwrite") ;
	    cerr << "write to temporary file failed"<< endl ;
	    fclose(scratch) ;

	    Loci::Abort() ;
	  }
	  unsigned short clsz = clusterout.size() ;
	  cluster_sizes.push_back(clsz) ;
        }
      }
      
      H5Gclose(face_g) ;
      if(needGlue){
        faceMap_start += fileData[i].numBdFaces;
        nodeMap_start += fileData[i].numBdNodes;
      }  
    }//for(int i = 0; i < num_inputs...

    rewind(scratch) ;

    hid_t face_id = H5Gcreate(output_fid,"face_info",
			      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

    // Write cluster sizes
    int rank = 1 ;
    hsize_t dimension = cluster_sizes.size() ;
    hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
    hsize_t start = 0 ;
    hsize_t stride = 1;

    hid_t dataset = H5Dcreate(face_id,"cluster_sizes",H5T_NATIVE_USHORT,
                              dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    hid_t memspace = H5Screate_simple(rank,&dimension,NULL) ;
    H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&start,&stride,&dimension,NULL) ;
    hid_t err = H5Dwrite(dataset,H5T_NATIVE_USHORT,memspace,dataspace,H5P_DEFAULT,
                         &cluster_sizes[0]) ;
    if(err<0) {
      cerr << "unable to write cluster sizes to '" << output_file << "'" << endl ;
      exit(-1) ;
    }
    H5Sclose(memspace) ;
    H5Sclose(dataspace) ;
    H5Dclose(dataset) ;

    // Write cluster_info

    long cluster_info_size = 0 ;
    for(size_t s=0;s<cluster_sizes.size();++s) 
      cluster_info_size += cluster_sizes[s] ;

    rank = 1 ;
    dimension = cluster_info_size ;
    dataspace = H5Screate_simple(rank,&dimension,NULL) ;
    start = 0 ;
    stride = 1 ;
    dataset = H5Dcreate(face_id,"cluster_info",H5T_NATIVE_UCHAR,dataspace,
                        H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    
    vector<unsigned char > data(8096) ;
    while(cluster_info_size > 0) {
      int sz = 8096 ;
      if(cluster_info_size < 8096)
        sz = cluster_info_size ;

      size_t rsz = fread(&data[0],sz,1,scratch) ;
      if(rsz != 1) {
        perror("fread") ;
        cerr << "read failed when reading temporary scratch file" << endl ;
        fclose(scratch) ;
        Loci::Abort() ;
      }

      hsize_t dim = sz ;
      hsize_t count = sz ;
      memspace = H5Screate_simple(rank,&dim,NULL) ;
      H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,&start,&stride,&count,NULL) ;
      err = H5Dwrite(dataset,H5T_NATIVE_UCHAR,memspace,dataspace,H5P_DEFAULT,
                     &data[0]) ;
      if(err < 0) {
        cerr << "unable to write cluster_info" << endl ;
        exit(-1) ;
      }
      H5Sclose(memspace) ;
      cluster_info_size -= sz ;
      start += sz ;
    }

    H5Sclose(dataspace) ;
    H5Dclose(dataset) ;
    H5Gclose(face_id) ;
    fclose(scratch) ;
  }

  cout << "boundary conditions in output are:" << endl ;
  map<string,int>::const_iterator mi ;
  for(mi=newbcs.begin();mi!=newbcs.end();++mi) 
    cout << mi->first <<  ' ' << mi->second << endl ;

  // Write out new boundary conditions
  hid_t surf_id = H5Gcreate(output_fid,"surface_info",
			    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
  for(mi=newbcs.begin();mi!=newbcs.end();++mi) {
    hid_t bc_id = H5Gcreate(surf_id,mi->first.c_str(),
			    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    hsize_t dims = 1 ;
    hid_t dataspace_id = H5Screate_simple(1,&dims,NULL) ;
    
    hid_t att_id = H5Acreate(bc_id,"Ident", H5T_NATIVE_INT,
                             dataspace_id, H5P_DEFAULT,H5P_DEFAULT) ;
    H5Awrite(att_id,H5T_NATIVE_INT,&(mi->second)) ;
    H5Aclose(att_id) ;
    H5Gclose(bc_id) ;
  }
  H5Gclose(surf_id) ;


  // Compute Volume Tags
  map<string,entitySet> volTagMap ;
  int cellcnt = 0 ;
  for(int i=0;i<num_inputs;++i) {
    for(size_t j=0;j<volTag[i].size();++j) {
      string tag = volTag[i][j].first ;
      entitySet tagSet = (volTag[i][j].second >> cellcnt) ;
      volTagMap[tag] += tagSet ;
    }
    cellcnt+= fileData[i].numCells ;
  }
  vector<pair<string,entitySet> > volTags ;

  if(volTagMap["Main"] != Loci::EMPTY)
    volTags.push_back(pair<string,entitySet>(string("Main"),volTagMap["Main"])) ;
  
                    
  map<string,entitySet>::const_iterator vtmi ;
  for(vtmi=volTagMap.begin();vtmi!=volTagMap.end();++vtmi) {
    if(vtmi->first != "Main")
      volTags.push_back(pair<string,entitySet>(vtmi->first,vtmi->second)) ;
  }

  cout << "volume tags =" << endl  ;
  for(size_t i=0;i<volTags.size();++i)
    cout << i << " - " << volTags[i].first << volTags[i].second << endl ;


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
    Loci::HDF5_WriteDomain(vol_id,volTags[i].second, MPI_COMM_WORLD) ;
    H5Gclose(vol_id) ;
  }
  H5Gclose(cell_info) ;
  
  // Close everything up
  for(int i=0;i<num_inputs;++i) 
    H5Fclose(input_fid[i]) ;
  H5Fclose(output_fid) ;
  Loci::Finalize() ;
}
