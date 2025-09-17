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
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::ofstream ;
using std::ios ;
using std::list ;

#include "extract.h"

namespace Loci {
  //Define struct Area which data members of normal vector and area of the area
  struct Aread {
    vector3d<double> n ;  //normal vector of the face
    double sada ; //area of the face
  } ;

  //Overload ostream and istream (Input/Output) operators for struct Area
  inline std::ostream & operator<<(std::ostream &s, const Aread &v)
  {
    s << v.n << ' ' << v.sada << ' ' ;
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, Aread &v)
  {
    s >> v.n >> v.sada  ;
    return s ;
  }

  template<> struct data_schema_traits<Loci::Aread> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(Loci::Aread()) ;
      LOCI_INSERT_TYPE(ct,Loci::Aread,n) ;
      LOCI_INSERT_TYPE(ct,Loci::Aread,sada) ;
      return DatatypeP(ct) ;
    }
  } ;
}
  


void process_ascii_nodal(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames) {
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


  Loci::hdf5CloseFile(file_id) ;

  list<vector<float> > values ;
  for(size_t i=0;i<variables.size();++i) {
    string var_name = variables[i] ;
    string filename = variable_filenames[i] ;
    switch(variable_types[i]) {
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
        Loci::hdf5CloseFile(file_id) ;

        values.push_back(vector<float>(0)) ;
        vector<float> val(npnts) ;
        int c = 0 ;
        FORALL(dom,nd) {
          val[c++] = scalar[nd] ;
        } ENDFORALL ;
        values.back().swap(val) ;
      }      
      break;
    case NODAL_DERIVED:
      {
        vector<float> dval(npnts) ;
        getDerivedVar(dval,var_name,casename,iteration) ;

        values.push_back(vector<float>(0)) ;
        values.back().swap(dval) ;
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
        Loci::hdf5CloseFile(file_id) ;
        values.push_back(vector<float>(npnts)) ;
        int c = 0 ;
        FORALL(dom,nd) {
          values.back()[c++] = vec[nd].x ;
        } ENDFORALL ;
        values.push_back(vector<float>(npnts)) ;
        c = 0 ;
        FORALL(dom,nd) {
          values.back()[c++] = vec[nd].y ;
        } ENDFORALL ;
        values.push_back(vector<float>(npnts)) ;
        c = 0 ;
        FORALL(dom,nd) {
          values.back()[c++] = vec[nd].z ;
        } ENDFORALL ;
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
        
          values.push_back(vector<float>(0)) ;
          values.back().swap(vec) ;
        }
      }
      break ;
    default:
      cerr << "unable to process variable " << var_name << "!"<< endl ;
      cerr << "this variable is ignored!" << endl ;
      break ;
    }
  }
  for(int i=0;i<npnts;++i) {
    list<vector<float> >::const_iterator li= values.begin() ;
    cout << (*li)[i] ;
    ++li ;
    for(;li!=values.end();++li) 
      cout << ' ' << (*li)[i] ;
    cout << endl ;
  }

}

void process_ascii_bndry(string casename, string iteration,
                         vector<string> variables,
                         vector<int> variable_types,
                         vector<string> variable_filenames,
                         vector<string> boundaries) {

  string gridtopo = getTopoFileName(output_dir, casename, iteration) ;


  hid_t file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;

#ifdef H5_USE_16_API
  hid_t bndg = H5Gopen(file_id,"boundaries") ;
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT) ;
#endif
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  vector<string> processed_bcs ;

  vector<int> elem_ids ;
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
#ifdef H5_USE_16_API
    hid_t bcg = H5Gopen(bndg,buf) ;
#else
    hid_t bcg = H5Gopen(bndg,buf,H5P_DEFAULT) ;
#endif
    int nquads = sizeElementType(bcg,"quads") ;
    int ntrias = sizeElementType(bcg,"triangles") ;
    int ngeneral = sizeElementType(bcg,"nside_sizes") ;
    vector<int > trias_id(ntrias) ;
    readElementType(bcg,"triangles_id",trias_id) ;
    vector<int > quads_id(nquads) ;
    readElementType(bcg,"quads_id",quads_id) ;
    vector<int > nside_id(ngeneral) ;
    readElementType(bcg,"nside_id",nside_id) ;
    for(int i=0;i<ntrias;++i)
      elem_ids.push_back(trias_id[i]) ;
    for(int i=0;i<nquads;++i)
      elem_ids.push_back(quads_id[i]) ;
    for(int i=0;i<ngeneral;++i)
      elem_ids.push_back(nside_id[i]) ;
    H5Gclose(bcg) ;
  }
  H5Gclose(bndg) ;
  H5Fclose(file_id) ;

  map<int,int> nmap ;
  int nelem = elem_ids.size() ;
  for(int i=0;i<nelem;++i)
    nmap[elem_ids[i]] = i ;

  bool needxyz = false ;
  bool needa = false ;
  for(size_t i=0;i<variables.size();++i) {
    if(variables[i] == "x" || variables[i] == "y" || variables[i] == "z")
      needxyz = true ;
    if(variables[i] == "n" || variables[i] == "area")
      needa = true ;
  }
  int nalloc = 0 ;
  if(needxyz || needa)
    nalloc = nelem ;
  vector<float> xc(nalloc),yc(nalloc),zc(nalloc) ;
  vector<float> nx(nalloc),ny(nalloc),nz(nalloc),area(nalloc) ;
  if(needxyz || needa) {
    string filename = "output/bc_geom." + iteration + "_" + casename ;

    hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      Loci::Abort() ;
      exit(-1) ;
    }
        
#ifdef H5_USE_16_API
    hid_t di = H5Gopen(file_id,"dataInfo") ;
#else
    hid_t di = H5Gopen(file_id,"dataInfo",H5P_DEFAULT) ;
#endif

    if(di < 0) {
      cerr << "unable to open group dataInfo in file " << filename << endl ;
    }
    int nbel = sizeElementType(di,"entityIds") ;

    vector<int> elemIds(nbel) ;
    readElementType(di,"entityIds",elemIds) ;
    
    H5Gclose(di) ;
    if(needxyz) {
      vector<vector3d<float> > facecenter(nbel) ;
      readElementType(file_id,"facecenter",facecenter) ;

      map<int,int>:: const_iterator mi ;
      for(int j=0;j<nbel;++j) {
        if((mi=nmap.find(elemIds[j])) != nmap.end()) {
          xc[mi->second] = facecenter[j].x ;
          yc[mi->second] = facecenter[j].y ;
          zc[mi->second] = facecenter[j].z ;
        }
      }
    }
    if(needa) {
      vector<Loci::Aread> a(nbel) ;
      readElementType(file_id,"area",a) ;

      map<int,int>:: const_iterator mi ;
      for(int j=0;j<nbel;++j) {
        if((mi=nmap.find(elemIds[j])) != nmap.end()) {
          nx[mi->second] = a[j].n.x ;
          ny[mi->second] = a[j].n.y ;
          nz[mi->second] = a[j].n.z ;
          area[mi->second] = a[j].sada ;
        }
      }
    }      
    H5Fclose(file_id) ;      
  }

  list< vector<float> > values ;
  for(size_t i=0;i<variables.size();++i) {
    switch(variable_types[i]) {
    case BOUNDARY_SCALAR:
      {
        const string filename(variable_filenames[i]) ;
        hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                           H5F_ACC_RDONLY,
                                           H5P_DEFAULT) ;
        if(file_id < 0) {
          cerr << "unable to open file '" << filename << "'!" << endl ;
          continue ;
        }
        
#ifdef H5_USE_16_API
        hid_t di = H5Gopen(file_id,"dataInfo") ;
#else
        hid_t di = H5Gopen(file_id,"dataInfo",H5P_DEFAULT) ;
#endif
        int nbel = sizeElementType(di,"entityIds") ;
        
        vector<int> elemIds(nbel) ;
        readElementType(di,"entityIds",elemIds) ;
        
        H5Gclose(di) ;
        vector<float> var(nbel) ;
        readElementType(file_id,variables[i].c_str(),var) ;
        H5Fclose(file_id) ;
        vector<float> val(nelem) ;
        int cnt = 0 ;
        map<int,int>:: const_iterator mi ;
        for(int j=0;j<nbel;++j) {
          if((mi=nmap.find(elemIds[j])) != nmap.end()) {
            val[mi->second] = var[j] ;
            cnt++ ;
          }
        }
        if(cnt == nelem) {
          values.push_back(vector<float>(0)) ;
          values.back().swap(val) ;
        } else {
          cerr << "Unable to extract values for all boundary elements for variable " << variables[i] << "..." << endl << "Skipping." << endl ;
        }
      }
      break;
    case BOUNDARY_VECTOR:
      {
        const string filename(variable_filenames[i]) ;
        hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                           H5F_ACC_RDONLY,
                                           H5P_DEFAULT) ;
        if(file_id < 0) {
          cerr << "unable to open file '" << filename << "'!" << endl ;
          continue ;
        }
        
#ifdef H5_USE_16_API
        hid_t di = H5Gopen(file_id,"dataInfo") ;
#else
        hid_t di = H5Gopen(file_id,"dataInfo",H5P_DEFAULT) ;
#endif
        int nbel = sizeElementType(di,"entityIds") ;
        
        vector<int> elemIds(nbel) ;
        readElementType(di,"entityIds",elemIds) ;
        
        H5Gclose(di) ;
        vector<vector3d<float> > var(nbel) ;
        readElementType(file_id,variables[i].c_str(),var) ;
        H5Fclose(file_id) ;
        vector<float> valx(nelem) ;
        vector<float> valy(nelem) ;
        vector<float> valz(nelem) ;
        int cnt = 0 ;
        map<int,int>:: const_iterator mi ;
        for(int j=0;j<nbel;++j) {
          if((mi=nmap.find(elemIds[j])) != nmap.end()) {
            valx[mi->second] = var[j].x ;
            valy[mi->second] = var[j].y ;
            valz[mi->second] = var[j].z ;
            cnt++ ;
          }
        }
        if(cnt == nelem) {
          values.push_back(vector<float>(0)) ;
          values.back().swap(valx) ;
          values.push_back(vector<float>(0)) ;
          values.back().swap(valy) ;
          values.push_back(vector<float>(0)) ;
          values.back().swap(valz) ;
        } else {
          cerr << "Unable to extract values for all boundary elements for variable " << variables[i] << "..." << endl << "Skipping." << endl ;
        }
      }
      break ;
    case BOUNDARY_DERIVED_SCALAR:
      if(variables[i] == "x") {
        values.push_back(xc) ;
      } else if(variables[i] == "y") {
        values.push_back(yc) ;
      } else if(variables[i] == "z") {
        values.push_back(zc) ;
      } else if(variables[i] == "area") {
        values.push_back(area) ;
      } else
        cerr << "unknown derived boundary variable " << variables[i] <<endl ;
                  
      break ;
    case BOUNDARY_DERIVED_VECTOR:
      if(variables[i] == "n") {
        values.push_back(nx) ;
        values.push_back(ny) ;
        values.push_back(nz) ;
      } else {
        cerr << "unknown derived boundary variable " << variables[i] << endl ;
      }
      break ;
    default:
      cerr << "variable " << variables[i] << "not a boundary variable, ignored!" ;
      break ;
    }
  }
  if(values.size() == 0) {
    cerr << "no values to extract!" << endl ;
    Loci::Abort() ;
    exit(-1) ;
  }
  for(int i=0;i<nelem;++i) {
    list<vector<float> >::const_iterator li= values.begin() ;
    cout << (*li)[i] ;
    ++li ;
    for(;li!=values.end();++li) 
      cout << ' ' << (*li)[i] ;
    cout << endl ;
  }
  
}

