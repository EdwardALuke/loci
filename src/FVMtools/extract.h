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
#ifndef EXTRACT_H
#define EXTRACT_H
#include <sstream>
#include <string>
#include <stdlib.h>

using std::list;

enum var_type {NODAL_SCALAR,NODAL_VECTOR,NODAL_DERIVED,NODAL_MASSFRACTION, BOUNDARY_SCALAR,BOUNDARY_VECTOR, BOUNDARY_DERIVED_SCALAR, BOUNDARY_DERIVED_VECTOR, PARTICLE_SCALAR, PARTICLE_VECTOR, UNDEFINED} ;
enum view_type {VIEWXY=0,VIEWYZ=1,VIEWXZ=2,VIEWXR=3}  ;


extern string getTopoFileName(string output_dir, string casename, string iteration) ;

// convert a string to an integer
inline int str2int(string s) {
  std::stringstream ss ;
  ss << s ;
  int ret ;
  ss >> ret ;
  return ret ;
}

// determine whether range [b,e) contains only digits
// [b,e) must be a valid range and contains characters
template <class ForwardIter> inline bool
alldigits(ForwardIter b, ForwardIter e) {
  for(;b!=e;++b)
    if(!isdigit(*b))
      return false ;
                                                                                
  return true ;
}
// determine whether s contains a valid integer (including the sign)
inline bool
valid_int(const std::string& s) {
  if(s.empty())
    return false ;
  if( (s[0] == '-') || (s[0] == '+')) {
    if(s.size() > 1)
      return alldigits(s.begin()+1,s.end()) ;
    else
      return false ;
  }else
    return alldigits(s.begin(),s.end()) ;
}



struct affineMapping {
  float M[4][4] ;
  affineMapping() ;
  void Combine(affineMapping a) ;
  void translate(vector3d<float> tv) ;
  void rotateX(float theta) ;
  void rotateY(float theta) ;
  void rotateZ(float theta) ;
  vector3d<float> MapNode(vector3d<float> v) ;
};

void get_2dgv(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries,
              int view) ;

void get_surf(string casename, string iteration,
              vector<string> variables,
              vector<int> variable_types,
              vector<string> variable_filenames,
              vector<string> boundaries) ;

void process_ascii_nodal(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames) ;

void process_ascii_bndry(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames,
                         vector<string> boundaries) ;

void process_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter) ;

void combine_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter,
		  bool do_favre) ;

size_t  sizeElementType(hid_t group_id, const char *element_name) ;

template<class T> void readElementType(hid_t group_id, const char *element_name,
                                       vector<T> &v) {
  if(v.size() > 0) {
    hid_t dataset = H5Dopen(group_id,element_name,H5P_DEFAULT) ;
    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;
    
    H5Dread(dataset,dp->get_hdf5_type(),H5S_ALL,H5S_ALL,H5P_DEFAULT,&v[0]) ;
    H5Dclose(dataset) ;
  }
}

template<class T>
void readElementTypeBlock(hid_t group_id,
                          const char *element_name,
                          vector<T> &v,
                          size_t start_elem,
                          int block_size) {
  if(block_size > 0) {
    hid_t dataset = H5Dopen(group_id,element_name,H5P_DEFAULT) ;

    typedef data_schema_traits<T> traits_type ;
    Loci::DatatypeP dp = traits_type::get_type() ;

    hsize_t count = block_size ;
    hsize_t start = start_elem ;
    hid_t dataspace = H5Dget_space(dataset) ;
    hsize_t stride = 1 ;
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
    int rank = 1 ;
    hsize_t dimension = count ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    hid_t datatype = dp->get_hdf5_type() ;
    H5Dread(dataset,datatype,memspace,dataspace,H5P_DEFAULT,&v[0]) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    H5Sclose(dataspace) ;
    H5Dclose(dataset) ;
  }
}


template<class T> void writeElementType(hid_t group_id,
                                        const char *element_name,
                                        std::vector<T> &v) {
  hsize_t array_size = v.size() ;
  if(array_size == 0)
    return ;
  int rank = 1 ;
  hsize_t dimension = array_size ;

  hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

  typedef data_schema_traits<T> traits_type ;
  Loci::DatatypeP dp = traits_type::get_type() ;
  
  hsize_t start = 0 ;
  hsize_t stride = 1 ;
  hsize_t count = v.size() ;
  hid_t datatype = dp->get_hdf5_type() ;

  hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                            dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

  if(count != 0) {
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                        &start, &stride, &count, NULL) ;
    hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
    H5Dwrite(dataset,datatype,memspace,dataspace,
             H5P_DEFAULT, &v[0]) ;
    H5Sclose(memspace) ;
  }
  H5Dclose(dataset) ;
  H5Sclose(dataspace) ;
  H5Tclose(datatype) ;
}
  
  
void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration) ;

namespace Loci {
  inline bool operator<(const Array<int,3> &a1,
                        const Array<int,3> &a2) {
    return a1[0]<a2[0] || (a1[0] == a2[0] && a1[1] < a2[1]) ||
      (a1[0] == a2[0] && a1[1] == a2[1] && a1[2] < a2[2]) ;
  }
}

extern string output_dir ;

string getPosFile(string output_dir,string iteration, string casename) ;


class particlePartBase : public Loci::CPTR_type {
protected:
  bool error ;
  string partName ;
  size_t numParticles ;
public:
  bool fail() const { return error ; }
  string getPartName() const { return partName ; }
  size_t getNumParticles() const { return numParticles ; }
  virtual bool hasScalarVar(string var) const = 0 ;
  virtual bool hasVectorVar(string var) const = 0 ;
  virtual std::vector<string> getScalarVars() const = 0 ;
  virtual std::vector<string> getVectorVars() const = 0 ;
  virtual void getParticlePositions(vector<vector3d<float> > &ppos) const = 0;
  virtual void getParticleScalar(string varname, vector<float> &val) const = 0 ;
  virtual void getParticleVector(string varname, 
				 vector<vector3d<float> > &val) const  = 0 ;
} ;

typedef Loci::CPTR<particlePartBase> particlePartP ;

class particlePart: public particlePartBase {
  string directory ;
  string posfile ;
  map<string,string> scalarVars ;
  map<string,string> vectorVars ;
  int stride_size ;
public:
  particlePart() { error = true ; }
  particlePart(string output_dir, string iteration, string casename,
	       vector<string> vars, int maxparticles) ;
  virtual bool hasScalarVar(string var) const ;
  virtual bool hasVectorVar(string var) const ;
  virtual std::vector<string> getScalarVars() const ;
  virtual std::vector<string> getVectorVars() const ;
  virtual void getParticlePositions(vector<vector3d<float> > &ppos) const ;
  virtual void getParticleScalar(string varname, vector<float> &val) const ;
  virtual void getParticleVector(string varname, 
				 vector<vector3d<float> > &val) const ;
} ;

// Create abstraction for parts
class volumePartBase : public Loci::CPTR_type {
protected:
  bool error ;
  string partName ;
  size_t nnodes ;
  size_t ntets, nhexs, nprsm, npyrm, ngenc ; 
  size_t ntetsIblank, nhexsIblank,nprsmIblank,npyrmIblank,ngencIblank ;
public:
  bool fail() const { return error ; }
  string getPartName() const { return partName ; }
  size_t getNumNodes() const { return nnodes ; }
  
  size_t getNumTets() const { return ntets ; }
  size_t getNumHexs() const { return nhexs ; }
  size_t getNumPrsm() const { return nprsm ; }
  size_t getNumPyrm() const { return npyrm ; }
  size_t getNumGenc() const { return ngenc ; }
  size_t getNumTetsIblank() const { return ntetsIblank ; }
  size_t getNumHexsIblank() const { return nhexsIblank ; }
  size_t getNumPrsmIblank() const { return nprsmIblank ; }
  size_t getNumPyrmIblank() const { return npyrmIblank ; }
  size_t getNumGencIblank() const { return ngencIblank ; }

  virtual bool hasNodalScalarVar(string var) const = 0 ;
  virtual bool hasNodalVectorVar(string var) const = 0 ;
  virtual std::vector<string> getNodalScalarVars() const = 0 ;
  virtual std::vector<string> getNodalVectorVars() const = 0 ;
  virtual void getPos(vector<vector3d<float> > &pos) const = 0 ;
  virtual void getPos(vector<vector3d<double> > &pos) const = 0 ;
  virtual void getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const = 0 ;
  virtual void getTetIds(vector<int> &tetids, size_t start, size_t size) const = 0 ;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const = 0 ;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const = 0 ;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const = 0 ;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const = 0 ;
  virtual void getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const = 0 ;
  virtual void getHexIds(vector<int> &hexids, size_t start, size_t size) const = 0 ;
  virtual void getGenCell(vector<int> &genCellNfaces, 
			  vector<int> &genCellNsides,
			  vector<int> &genCellNodes) const = 0;
  virtual void getGenIds(vector<int> &genids) const = 0 ;
  
  virtual void getNodalScalar(string varname, vector<float> &vals) const = 0 ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const = 0 ;
  virtual void getNodalIblank(vector<unsigned char> &blank) const = 0 ;
} ;

typedef Loci::CPTR<volumePartBase> volumePartP ;

class volumePartDerivedVars : public volumePartBase {
  volumePartP shadowPart ;
  enum derivedVar_t { VAR_M, VAR_P, VAR_logp, VAR_U, VAR_0, VAR_1, VAR_2,
		      VAR_X, VAR_Y, VAR_Z } ;
  map<string,derivedVar_t> derivedVars ;
  float Pambient ;
  void processDerivedVars(const vector<string> &vars) ;
public:
  volumePartDerivedVars() {error = true ;}
  volumePartDerivedVars(volumePartP part,
			string output_dir, string iteration, string casename,
			vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual std::vector<string> getNodalScalarVars() const ;
  virtual std::vector<string> getNodalVectorVars() const ;
  virtual void getPos(vector<vector3d<float> > &val) const ;
  virtual void getPos(vector<vector3d<double> > &pos) const;
  virtual void getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const ;
  virtual void getTetIds(vector<int> &tetids, size_t start, size_t size) const ;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const ;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const ;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const ;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const ;
  virtual void getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const ;
  virtual void getHexIds(vector<int> &hexids, size_t start, size_t size) const ;
  virtual void getGenCell(vector<int> &genCellNfaces, 
			  vector<int> &genCellNsides,
			  vector<int> &genCellNodes) const ;
  virtual void getGenIds(vector<int> &genids) const ;
  
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getNodalIblank(vector<unsigned char> &blank) const ;
} ;

class volumePart : public volumePartBase {
  string directory ;
  string topoFile ;
  string posFile ;
  // maps from variables to file name
  map<string,string> nodalScalarVars ;
  map<string,string> nodalVectorVars ;
  bool has_iblank ;
  store<unsigned char> iblank ;
  size_t ntets_orig, nhexs_orig, nprsm_orig, npyrm_orig, ngenc_orig ; 
  Loci::entitySet tetsIblanked, hexsIblanked, prsmIblanked, pyrmIblanked ;
  Loci::entitySet gencIblanked ;
public:
  volumePart() {error = true ;}
  volumePart(string output_dir, string iteration, string casename,
             vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual std::vector<string> getNodalScalarVars() const ;
  virtual std::vector<string> getNodalVectorVars() const ;
  virtual void getPos(vector<vector3d<float> > &val) const ;
  virtual void getPos(vector<vector3d<double> > &pos) const;
  virtual void getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const ;
  virtual void getTetIds(vector<int> &tetids, size_t start, size_t size) const ;
  virtual void getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const ;
  virtual void getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const ;
  virtual void getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const ;
  virtual void getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const ;
  virtual void getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const ;
  virtual void getHexIds(vector<int> &hexids, size_t start, size_t size) const ;
  virtual void getGenCell(vector<int> &genCellNfaces, 
			  vector<int> &genCellNsides,
			  vector<int> &genCellNodes) const ;
  virtual void getGenIds(vector<int> &genids) const ;
  
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getNodalIblank(vector<unsigned char> &blank) const ;
} ;

class surfacePartBase : public Loci::CPTR_type {
protected:
  bool error ;
  string partName ;
  int nnodes, nquads, ntrias, ngenf ;
public:
  bool fail() const { return error ; }
  string getPartName() const { return partName ; }
  int getNumNodes() const { return nnodes ; }
  int getNumQuads() const { return nquads ; }
  int getNumTrias() const { return ntrias ; }
  int getNumGenfc() const { return ngenf ; }
  virtual bool hasNodalScalarVar(string var) const = 0 ;
  virtual bool hasNodalVectorVar(string var) const = 0 ;
  virtual bool hasElementScalarVar(string var) const = 0 ;
  virtual bool hasElementVectorVar(string var) const = 0 ;
  virtual std::vector<string> getNodalScalarVars() const = 0 ;
  virtual std::vector<string> getNodalVectorVars() const = 0 ;
  virtual std::vector<string> getElementScalarVars() const = 0 ;
  virtual std::vector<string> getElementVectorVars() const = 0 ;
  virtual void getQuads(vector<Array<int,4> > &quads) const = 0 ;
  virtual void getTrias(vector<Array<int,3> > &trias) const = 0 ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const =0;
  virtual void getQuadsIds(vector<int> &quads_ids) const = 0 ;
  virtual void getTriasIds(vector<int> &trias_ids) const = 0 ;
  virtual void getGenfIds(vector<int> &genface_ids) const =0;
  
  virtual void getPos(vector<vector3d<float> > &pos) const = 0;
  virtual void getPos(vector<vector3d<double> > &pos) const = 0 ;
  virtual void getNodalScalar(string varname, vector<float> &vals) const = 0 ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const = 0;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const =0 ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const = 0 ;
} ;
typedef Loci::CPTR<surfacePartBase> surfacePartP ;

class surfacePartDerivedVars : public surfacePartBase {
  surfacePartP shadowPart ;
  enum derivedVar_t { VAR_M, VAR_P, VAR_logp, VAR_U, VAR_0, VAR_1, VAR_2,
		      VAR_X, VAR_Y, VAR_Z } ;
  map<string,derivedVar_t> derivedVars ;
  float Pambient ;
  void processDerivedVars(const vector<string> &vars) ;
public:
  surfacePartDerivedVars() { error = true ; shadowPart=0;} 
  surfacePartDerivedVars(surfacePartP part, string output_dir,
			 string casename ,
			 string iteration, 
			 vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual bool hasElementScalarVar(string var) const ;
  virtual bool hasElementVectorVar(string var) const ;
  virtual vector<string> getNodalScalarVars() const ;
  virtual vector<string> getNodalVectorVars() const ;
  virtual vector<string> getElementScalarVars() const ;
  virtual vector<string> getElementVectorVars() const ;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const ;
  virtual void getTrias(vector<Array<int,3> > &trias) const ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const ;
  virtual void getQuadsIds(vector<int> &quads_ids) const  ;
  virtual void getTriasIds(vector<int> &trias_ids) const ;
  virtual void getGenfIds(vector<int> &genface_ids) const ;
  virtual void getPos(vector<vector3d<float> > &pos) const ;
  virtual void getPos(vector<vector3d<double> > &pos) const;
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const ;
} ;
    
      
class surfacePart : public surfacePartBase {
  string directory ;
  string topoFile ;
  string posFile ;
  // maps from variables to file name
  map<string,string> nodalScalarVars ;
  map<string,string> nodalVectorVars ;
  map<string,string> elementScalarVars ;
  map<string,string> elementVectorVars ;

  vector<int> quad_ord ;
  vector<int> tri_ord ;
  vector<int> gen_ord ;
  entitySet quadSet,triSet,genSet ;
public:
  surfacePart() {error = true ;}
  surfacePart(string name, string directory, string iteration,
	      vector<string> vars) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual bool hasElementScalarVar(string var) const ;
  virtual bool hasElementVectorVar(string var) const ;
  virtual vector<string> getNodalScalarVars() const ;
  virtual vector<string> getNodalVectorVars() const ;
  virtual vector<string> getElementScalarVars() const ;
  virtual vector<string> getElementVectorVars() const ;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const ;
  virtual void getTrias(vector<Array<int,3> > &trias) const ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const ;
  virtual void getQuadsIds(vector<int> &quads_ids) const  ;
  virtual void getTriasIds(vector<int> &trias_ids) const ;
  virtual void getGenfIds(vector<int> &genface_ids) const ;
  
  virtual void getPos(vector<vector3d<float> > &pos) const ;
  virtual void getPos(vector<vector3d<double> > &pos) const;
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const ;
} ;

class surfacePartCopy : public surfacePartBase {
  vector<Array<int,3> > trifaces ;
  vector<Array<int,4> > quadfaces ;
  vector<int> nfacenodes, gennodes ;
  vector<int> triaIds;
  vector<int> quadIds;
  vector<int> genIds;
  vector<int> nodemap ;
  vector<vector3d<double> > pos ;
  map<string,vector<float> > nodalScalars ;
  map<string,vector<vector3d<float> > > nodalVectors ;
  map<string,Array<vector<float>,3> > elementScalars ;
  map<string,Array<vector<vector3d<float> >,3> > elementVectors ;
public:
  surfacePartCopy() {error = true ;}
  surfacePartCopy(string name, vector<Array<int,3> > &triangles,
                  vector<int> &tria_ids,
                  vector<Array<int,4> > &quads,
                  vector<int>& quads_ids,
                  vector<int> &genface2n, vector<int> &gnodes,
                  vector<int>&gen_ids) ;
  void registerPos(const vector<vector3d<double> > &pos) ;
  void registerNodalScalar(string name,const vector<float> &val) ;
  void registerNodalVector(string name,const vector<vector3d<float> > &val) ;
  void registerElementScalar(string name, 
			     const vector<float> &qval,
			     const vector<float> &tval,
			     const vector<float> &gval) ;
  void registerElementVector(string name, 
			     const vector<vector3d<float> > &qval,
			     const vector<vector3d<float> > &tval,
			     const vector<vector3d<float> > &gval) ;
  virtual bool hasNodalScalarVar(string var) const ;
  virtual bool hasNodalVectorVar(string var) const ;
  virtual bool hasElementScalarVar(string var) const ;
  virtual bool hasElementVectorVar(string var) const ;
  virtual vector<string> getNodalScalarVars() const ;
  virtual vector<string> getNodalVectorVars() const ;
  virtual vector<string> getElementScalarVars() const ;
  virtual vector<string> getElementVectorVars() const ;
  
  virtual void getQuads(vector<Array<int,4> > &quads) const ;
  virtual void getTrias(vector<Array<int,3> > &trias) const ;
  virtual void getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const ;
  virtual void getQuadsIds(vector<int> &quads_ids) const;  
  virtual void getTriasIds(vector<int> &trias_ids) const; 
  virtual void getGenfIds(vector<int> &genface_ids) const; 
  
  virtual void getPos(vector<vector3d<float> > &pos) const ;
  virtual void getPos(vector<vector3d<double> > &pos) const ;
  virtual void getNodalScalar(string varname, vector<float> &vals) const ;
  virtual void getNodalVector(string varname, vector<vector3d<float> > &vals) const ;
  virtual void getElementScalar(string varname, vector<float> &qvals,
				vector<float> &tvals, 
				vector<float> &gvals) const ;
  virtual void getElementVector(string varname, vector<vector3d<float> > &qvals,
				vector<vector3d<float> > &tvals,
				vector<vector3d<float> > &gvals) const ;
} ;

class postProcessorConvert : public Loci::CPTR_type {
protected:
  vector<surfacePartP> surfacePartList ;
  vector<volumePartP> volumePartList ;
  vector<particlePartP> particlePartList ;
public:
  void addSurfaceParts(const vector<surfacePartP> &list) {
    for(size_t i=0;i<list.size();++i)
      surfacePartList.push_back(list[i]) ;
  }
  void addVolumePart(volumePartP volpart) {
    volumePartList.push_back(volpart) ;
  }
  void addParticlePart(particlePartP particles) {
    particlePartList.push_back(particles) ;
  }
  virtual bool processesVolumeElements() const = 0 ;
  virtual bool processesSurfaceElements() const = 0 ;
  virtual bool processesParticleElements() const = 0 ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const = 0 ;
} ;

typedef Loci::CPTR<postProcessorConvert> postProcessorP ;

class ensightPartConverter : public postProcessorConvert {
  bool id_required;
public:
  ensightPartConverter(bool input) {id_required = input; } ;
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

//#ifdef HAVE_CGNS
class cgnsPartConverter : public postProcessorConvert {
  bool id_required;
public:
  cgnsPartConverter(bool input) {id_required = input; } ;
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;
//#endif

class tecplotPartConverter : public postProcessorConvert {
public:
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class vtkPartConverter : public postProcessorConvert {
  bool bit64 ;
public:
  vtkPartConverter() {bit64 = false; } 
  vtkPartConverter(bool input) {bit64 = input; } ;
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class vtkSurfacePartConverter : public postProcessorConvert {
  bool bit64 ;
public:
  vtkSurfacePartConverter() {bit64=false; }
  vtkSurfacePartConverter(bool input) { bit64=input; }
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class fieldViewPartConverter : public postProcessorConvert {
public:
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

class cuttingPlanePartConverter : public postProcessorConvert 
{
  affineMapping transformMatrix ;
  float xShift, yShift, zShift ;
public:
  cuttingPlanePartConverter(const affineMapping &m,
			    float xs,float ys,float zs) {
    transformMatrix = m ;
    xShift=xs ;
    yShift=ys ;
    zShift=zs ;
  }
  virtual bool processesVolumeElements() const ;
  virtual bool processesSurfaceElements() const ;
  virtual bool processesParticleElements() const ;

  virtual void exportPostProcessorFiles(string casename, string iteration) const ;
} ;

extern void readData(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet, fact_db &facts) ;

#endif
