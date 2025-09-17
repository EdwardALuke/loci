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
#ifndef DEFINES_H
#define DEFINES_H
#include <Loci.h>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>

using std::cerr;
using std::endl;
using std::ostream;
using std::istream;


//some machine, long int and int have the same size, both are 32 bits
typedef unsigned long int64;
typedef unsigned int int32;
typedef Loci::vector3d<double> vect3d;

//the maximum value of levels allowed, if any level of refinement
// plan is greater than MAXLEVEL, the integer coordinates will
//exceed integer limit
const int MAXLEVEL = sizeof(int64)*8 - 2;
//const int MAXLEVEL = 30;
//normalize a vector
const double NORMALIZE_ZERO_THRESHOLD = 1e-10 ;
//0: simple(node average), 1: wireframe for face, area_weighted for cell,
//2: extract
const int CENTROID = 1;

inline void normalize(vect3d& v) {
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

//this function is used in build_general_cell with quadface. 
inline bool int_equal(const vect3d& v1, const vect3d& v2) {
  return (int64(v1.x) == int64(v2.x) && int64(v1.y) == int64(v2.y));
}





//to reduce the memory allocation of a vector to its size
template<class T> inline void reduce_vector(std::vector<T>& v1){
  std::vector<T> tmpVec = v1;
  std::swap(v1, tmpVec);
}

//calculate the unweighted mass center of a set of points
inline vect3d point_center(const std::vector<vect3d>& fnodes){
  int numNodes = fnodes.size();
  vect3d center = vect3d(0.0, 0.0, 0.0);
  for(int i = 0; i<numNodes; i++){
    center += fnodes[i];
  }
  center = center/double(numNodes);
  return center;
 }


inline vect3d weighted_center(const std::vector<vect3d>& fnodes, const std::vector<double>& weights){
  vect3d nodesum = vect3d(0.0, 0.0, 0.0);
  double lensum = 0.0;
  int numNodes = fnodes.size();
  for(int i = 0; i<numNodes; i++){
    nodesum += weights[i]*fnodes[i];
    lensum += weights[i];
  }

  return nodesum/lensum;
}
    
//the indexes of two neighbors of each face
struct NeibIndex{
  int c1;
  int c2;
  NeibIndex(int ii, int jj):c1(ii), c2(jj){}
  NeibIndex():c1(0), c2(0){}
};

namespace Loci {  
  struct SetLong{
  
    std::set<unsigned long> aset;
    friend ostream& operator << (ostream &, const SetLong &);
    friend istream& operator >> (istream &, SetLong &);
  };


  class SetLong_SchemaConverter ;
  template<>
  struct data_schema_traits<SetLong> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef unsigned long Converter_Base_Type ;
    typedef SetLong_SchemaConverter Converter_Type ;
  };

  class SetLong_SchemaConverter {
    SetLong& RefObj ;
  public:
    explicit SetLong_SchemaConverter(SetLong& new_obj): RefObj(new_obj) {}
    int getSize() const {
      return RefObj.aset.size () ;
    }
    void getState(unsigned long* buf, int& size) {
      size = getSize() ;
      int ii = 0;         
      for(std::set<unsigned long>::const_iterator ci=RefObj.aset.begin();
          ci!=RefObj.aset.end();++ci)
        buf[ii++] = *ci ;
    }
    void setState(unsigned long* buf, int size) {
      RefObj.aset.clear() ;
      for(int i=0;i<size;++i)
        RefObj.aset.insert(buf[i]) ;
    }
  };

  inline std::ostream& operator << (std::ostream &s, const SetLong &obj){
    std::set<unsigned long>::const_iterator ci;
    std::set<unsigned long>::const_iterator begin = obj.aset.begin();
    std::set<unsigned long>::const_iterator end = obj.aset.end();

    s << obj.aset.size() << '{' ;
    for(ci = begin; ci!=end; ci++) s << *ci << ' ';
    s << '}';
    return s;
    
  }
  
  inline std::istream& operator >> (std::istream &s,  SetLong &obj){
    
    obj.aset.clear();
    int size;
    s >> size;
    char token;
    s >> token;
    if(token != '{') {
      cerr << "parse error in SetLong istream" << endl;
      exit(0);
    }
    unsigned long value;
    for(int i = 0; i < size; i++){
      s>>value;
      obj.aset.insert(value);
    }
    if(size != 0) s>>token;
    s>>token;
    
    if(token != '}') {
      cerr << "parse error in SetLong istream" << endl;
      exit(0);
    } 
  
    return s;
  }
  
  struct SetLongUnion {
    void operator()(SetLong &f1, const SetLong &f2) {
      SetLong result;
      std::set_union(f1.aset.begin(), f1.aset.end(), f2.aset.begin(), f2.aset.end(),
                std::inserter(result.aset, result.aset.begin())); 
      f1 = result;
    }
  };
}

namespace Loci {
  
  typedef std::vector<vect3d> FineNodes;
  
  class FineNodes_SchemaConverter ;
  template<>
  struct data_schema_traits<FineNodes> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef vect3d Converter_Base_Type ;
    typedef FineNodes_SchemaConverter Converter_Type ;
  } ;
  
  class FineNodes_SchemaConverter {
    FineNodes& RefObj ;
  public:
    explicit FineNodes_SchemaConverter(FineNodes& new_obj): RefObj(new_obj) {}
    int getSize() const {
      return RefObj.size () ;
    }
    void getState(vect3d* buf, int& size) {
      size = getSize() ;
      for(int ii = 0; ii < size; ii++)
        buf[ii] = RefObj[ii] ;
    }
    void setState(vect3d* buf, int size) {
      RefObj.resize(size); ;
      for(int i=0;i<size;++i)
        RefObj[i] = buf[i] ;
    }
  };
  
  inline std::istream&
  operator >> (std::istream& s, std::vector<vect3d>& sl) {
    
    int sz ;
    s >> sz ;
    sl.resize(sz);
    char tok;
    s >> tok;
    if(tok != '{'){
      cerr << " parse error in std::vector<vect3d> istream operator" << endl;
      exit(0);
    }
      
    for(int i=0;i<sz;++i) {
      s >> sl[i] ;
    }
    if(sz != 0) s>>tok;
    s>>tok;
    if(tok != '}'){
      cerr << " parse error in std::vector<vect3d> istream operator" << endl;
      exit(0);
    }
    
    return s ;
  }

  inline std::ostream&
  operator << (std::ostream& s, const std::vector<vect3d>& sl) {
    int size = sl.size();
    s << size <<'{'  ;
    
    for(int i = 0; i < size; ++i) {
      s << sl[i] <<' ';
    }
    s << '}' ;
    return s ;
  }
  
}


namespace Loci {

  typedef std::vector<std::vector<int> > FineFaces;
  
  class FineFaces_SchemaConverter ;
  template<>
  struct data_schema_traits<FineFaces> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef int Converter_Base_Type ;
    typedef FineFaces_SchemaConverter Converter_Type ;
  } ;

  class FineFaces_SchemaConverter {
    FineFaces& RefObj ;
  public:
    explicit FineFaces_SchemaConverter(FineFaces& new_obj): RefObj(new_obj) {}
    int getSize() const {
      int vec_size = RefObj.size();
      if(vec_size == 0) return 0;
      int total_size = 0;
      for(int i = 0; i < vec_size; i++) total_size += RefObj[i].size();
      return total_size + vec_size + 1;
    }
    void getState(int* buf, int& size) {
      size = getSize() ;
      if(size ==0) return;
      unsigned ii = 0;
      buf[ii++] = RefObj.size();
      for(unsigned int i = 0; i < RefObj.size(); i++){
        buf[ii++] = RefObj[i].size();
        for(unsigned int j = 0; j < RefObj[i].size(); j++)
          buf[ii++] = RefObj[i][j];
      }
    }
    void setState(int* buf, int size) {
      if(size == 0){
        RefObj.clear();
        return;
      }
    
      unsigned int ii = 0;
      RefObj.resize(buf[ii++]);
      for(unsigned int i=0; i<RefObj.size(); i++){
        RefObj[i].resize(buf[ii++]);
        for(unsigned int j = 0; j <RefObj[i].size(); j++)
          RefObj[i][j] = buf[ii++];
      }
    }
  };

  inline std::istream&
  operator >> (std::istream& s, std::vector<std::vector<int> >& sl) {
    
    int vec_size;
    s >> vec_size;
    sl.resize(vec_size);
    char tok;
    s>>tok;
     if(tok != '{'){
      cerr << " parse error in std::vector<vector<int> > istream operator" << endl;
      exit(0);
     }
     
     int size;
     for(int i = 0; i < vec_size; i++){
       s >> size;
       s>>tok;
       if(tok != '{'){
         cerr << " parse error in std::vector<vector<int> > istream operator" << endl;
         exit(0);
       }
     
       sl[i].resize(size);
       for(int j = 0; j <size; j++){
         s >> sl[i][j];
       }
       if(size != 0) s>>tok;
       s>>tok;
       if(tok != '}'){
         cerr << " parse error in std::vector<vector<int> > istream operator" << endl;
         exit(0);
       } 
     }
     if(vec_size!= 0) s >> tok;
     s>>tok;
     if(tok != '}'){
       cerr << " parse error in std::vector<vector<int> > istream operator" << endl;
       exit(0);
     } 
     return s ;
  }
  
  inline std::ostream&
  operator << (std::ostream& s, const std::vector<std::vector<int> >& sl) {
    s <<sl.size();
    s <<'{';
    for(unsigned int i = 0; i < sl.size(); i++)
      {
        s <<sl[i].size()<<'{' ;
        for(unsigned int j = 0; j <sl[i].size(); j++)
          {
            s << sl[i][j]<<' ' ;
          }
        s << '}';
      }
    s << '}';
    return s;
  }

}










struct logicalAnd {
  void operator()(bool &f1, const bool &f2) {
    f1 = f1 && f2;
    
  }
} ;  

//for a quadface, the two edge that need apply
struct TwoEdge{
  TwoEdge(){};
  
  pair<int64, int64> e0;
  pair<int64, int64> e3;
};



namespace Loci{
  storeRepP my_get_node_remap(entitySet nodes, entitySet faces, entitySet edges);
}


void colorMatrix(Map &cl, Map &cr, multiMap &face2node);


namespace Loci {
  extern hid_t writeVOGOpen(string filename);
  extern  void writeVOGSurf(hid_t file_id, std::vector<pair<int,string> > surface_ids);
  extern  void writeVOGClose(hid_t file_id) ;
  extern bool readVolTags(hid_t input_fid,
                   std::vector<pair<string,Loci::entitySet> > &volDat);

  extern void writeVOGNode(hid_t file_id,
                           Loci::storeRepP &pos,
                           const_store<Loci::FineNodes> &inner_nodes_cell,
                           const_store<Loci::FineNodes> &inner_nodes_face,
                           const_store<Loci::FineNodes> &inner_nodes_edge
                           );
  extern void writeVOGFace(hid_t file_id, Map &cl, Map &cr, multiMap &face2node) ;
  extern  unsigned long readAttributeLong(hid_t group, const char *name);
  
  bool setupFVMGridFromContainer(fact_db &facts,
                                 std::vector<entitySet>& local_nodes,
                                 std::vector<entitySet>& local_faces,
                                 std::vector<entitySet>& local_cells,
                                 store<vector3d<double> >& t_pos,
                                 Map& tmp_cl,
                                 Map& tmp_cr,
                                 multiMap& tmp_face2node,
                                 std::vector<pair<int,string> >& boundary_ids,
                                 std::vector<pair<string,entitySet> >& volTags,
				 Loci::storeRepP cellwts = 0) ;
  
  
  inline std::ostream &operator <<(std::ostream &s, const std::vector<std::pair<int32,int32> > &v) {
    s << v.size() << endl ;
    for(size_t i=0;i<v.size();++i) {
      s << v[i].first << ' ' << v[i].second << endl ;
    }
    return s ;
  }
  
  inline std::istream &operator >>(std::istream &s, std::vector<std::pair<int32,int32> > &v) {
    size_t sz ;
    s >> sz ;
    std::vector<std::pair<int32,int32> > tmp(sz) ;
    v.swap(tmp) ;
    for(size_t i=0;i<sz;++i) {
      s >> v[i].first >> v[i].second ;
    }
    return s ;
  }
}

#endif
