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
#include <strings.h>
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


void process_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter) {

  string postfix ="" ;
  string iter_part = "" ;
  size_t i = 0;
  while(iteration[i] >= '0' && iteration[i] <= '9' && i<iteration.size()) {
    iter_part += iteration[i] ;
    i++ ;
  }
  while(i<iteration.size()) {
    postfix += iteration[i] ;
    i++ ;
  }
  int start_iter = atoi(iter_part.c_str()) ;
  
  for(size_t i=0;i<variables.size();++i) {
    string var_name = variables[i] ;
    cout << "processing variable: " << var_name << endl ;
    string var = var_name ;
    string filename = variable_filenames[i] ;
    switch(variable_types[i]) {
    case NODAL_SCALAR:
      {
        store<double> mean ;
        store<double> M2 ;
        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          n = n + 1.0 ;
          string filename = output_dir+'/' + var + "_sca." +
            string(buf) + postfix + "_" + casename ;

          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT, MPI_COMM_WORLD) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            continue ;
          }

          fact_db facts ;
          store<float> scalar ;
          readData(file_id,var_name,scalar.Rep(),EMPTY,facts) ;
          entitySet dom = scalar.domain() ;
          Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
          if(it == start_iter) {
            mean.allocate(dom) ;
            M2.allocate(dom) ;
            FORALL(dom,ii) {
              mean[ii] = 0. ;
              M2[ii] = 0. ;
            } ENDFORALL ;
          }

          FORALL(dom,nd) {
            double delta = scalar[nd] - mean[nd] ;
            mean[nd] += delta/n ;
            M2[nd] += delta*(scalar[nd]-mean[nd]) ;
          } ENDFORALL ;
        }
        store<float> m ;
        store<float> v ;
        entitySet dom = mean.domain() ;
        m.allocate(dom) ;
        v.allocate(dom) ;
        FORALL(dom,ii) {
          m[ii] = mean[ii] ;
          v[ii] = M2[ii]/(n) ;
        } ENDFORALL ;
        string filename = output_dir+'/' + var + "Mean_sca."
          + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"Mean" ;
	cout << "Writing Variables: " << sname ;
        fact_db facts ;
        Loci::writeContainer(file_id,sname,m.Rep(),facts) ;

        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "Var_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Var" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,v.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	cout << endl ;
      }
      break;
    case NODAL_DERIVED:
      {
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
        readData(file_id,"pos",pos.Rep(),EMPTY,facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
        int npnts = pos.domain().size() ;
        entitySet dom = pos.domain() ;

        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        store<double> mean ;
        store<double> M2 ;
        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          
          n = n + 1.0 ;
          string filename = output_dir+'/' + var + "_sca." +
            string(buf) + postfix + "_" + casename ;

          string iter_it = string(buf)+postfix ;
          vector<float> dval(npnts) ;
          getDerivedVar(dval,var_name,casename,iter_it) ;
          Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
          if(it == start_iter) {
            mean.allocate(dom) ;
            M2.allocate(dom) ;
            FORALL(dom,ii) {
              mean[ii] = 0. ;
              M2[ii] = 0. ;
            } ENDFORALL ;
          }

          int cnt = 0 ;
          FORALL(dom,nd) {
            double delta = dval[cnt] - mean[nd] ;
            mean[nd] += delta/n ;
            M2[nd] += delta*(dval[cnt]-mean[nd]) ;
            cnt++ ;
          } ENDFORALL ;
        }
        store<float> m ;
        store<float> v ;
        m.allocate(dom) ;
        v.allocate(dom) ;
        FORALL(dom,ii) {
          m[ii] = mean[ii] ;
          v[ii] = M2[ii]/(n) ;
        } ENDFORALL ;
        string filename = output_dir+'/' + var + "Mean_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"Mean" ;
	cout << "Writing Variables: " << sname ;
        Loci::writeContainer(file_id,sname,m.Rep(),facts) ;

        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "Var_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Var" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,v.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	cout << endl ;
      }
      break;
    case NODAL_VECTOR:
      {
        store<vector3d<double> > mean ;
        store<vector3d<double> > M2 ;
	store<double> Muv,Muw,Mvw ;

        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          
          n = n + 1.0 ;
          string filename = output_dir+'/' + var + "_vec." +
            string(buf) + postfix + "_" + casename ;

          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT, MPI_COMM_WORLD) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            continue ;
          }

          fact_db facts ;
          store<vector3d<float> > vect ;
          readData(file_id,var_name,vect.Rep(),EMPTY,facts) ;
          entitySet dom = vect.domain() ;
          Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
          if(it == start_iter) {
            mean.allocate(dom) ;
            M2.allocate(dom) ;
	    Muv.allocate(dom) ;
	    Muw.allocate(dom) ;
	    Mvw.allocate(dom) ;
            FORALL(dom,ii) {
              mean[ii] = vector3d<double>(0.,0.,0.) ;
              M2[ii] = vector3d<double>(0.,0.,0.) ;
	      Muv[ii] = 0 ;
	      Muw[ii] = 0 ;
	      Mvw[ii] = 0 ;
            } ENDFORALL ;
          }

          FORALL(dom,nd) {
            vector3d<double> delta = vector3d<double>(vect[nd].x - mean[nd].x,
                                                      vect[nd].y - mean[nd].y,
                                                      vect[nd].z - mean[nd].z);
	    Muv[nd] += (n-1)*(vect[nd].x-mean[nd].x)*
	      (vect[nd].y-mean[nd].y)/n ;
	    Muw[nd] += (n-1)*(vect[nd].x-mean[nd].x)*
	      (vect[nd].z-mean[nd].z)/n ;
	    Mvw[nd] += (n-1)*(vect[nd].y-mean[nd].y)*
	      (vect[nd].z-mean[nd].z)/n ;
            mean[nd] += (1/n)*delta ;
            M2[nd].x += delta.x*(vect[nd].x-mean[nd].x) ;
            M2[nd].y += delta.y*(vect[nd].y-mean[nd].y) ;
            M2[nd].z += delta.z*(vect[nd].z-mean[nd].z) ;
          } ENDFORALL ;
        }
        store<vector3d<float> > m ;
        store<vector3d<float> > v ;
	store<float> cuv,cuw,cvw ;
        entitySet dom = mean.domain() ;
        m.allocate(dom) ;
        FORALL(dom,ii) {
          m[ii].x = mean[ii].x ;
          m[ii].y = mean[ii].y ;
          m[ii].z = mean[ii].z ;
	} ENDFORALL ;
	mean.allocate(EMPTY) ;
        v.allocate(dom) ;
        FORALL(dom,ii) {
          v[ii].x = M2[ii].x/(n) ;
          v[ii].y = M2[ii].y/(n) ;
          v[ii].z = M2[ii].z/(n) ;
	} ENDFORALL ;
	M2.allocate(EMPTY) ;
	cuv.allocate(dom) ;
	cuw.allocate(dom) ;
	cvw.allocate(dom) ;
        FORALL(dom,ii) {
	  cuv[ii] = Muv[ii]/n ;
	  cuw[ii] = Muw[ii]/n ;
	  cvw[ii] = Mvw[ii]/n ;
        } ENDFORALL ;
	Muv.allocate(EMPTY);
	Muw.allocate(EMPTY);
	Mvw.allocate(EMPTY);

        string filename = output_dir+'/' + var + "Mean_vec."
          + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"Mean" ;
	cout << "Writing Variables: " << sname ;
        fact_db facts ;
        Loci::writeContainer(file_id,sname,m.Rep(),facts) ;

        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "Var_vec."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Var" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,v.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "Cuv_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Cuv" ;
	cout << ", " << sname ;

        Loci::writeContainer(file_id,sname,cuv.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
        filename = output_dir+'/' + var + "Cuw_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Cuw" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,cuw.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
        Loci::writeContainer(file_id,sname,cuv.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "Cvw_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Cvw" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,cvw.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	cout << endl ;
      }
      break;
    case NODAL_MASSFRACTION:
      {
        store<double> mean ;
        store<double> M2 ;
        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          
          n = n + 1.0 ;
          string int_str = string(buf) + postfix ;
          store<float> scalar ;
          string filename = "output/mix." + int_str + "_" + casename ;
    
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT, MPI_COMM_WORLD) ;
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
          Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
      
          map<string,int> smap ;
          std::istringstream iss(*species_names) ;
          for(int i=0;i<mix.vecSize();++i) {
            string s ;
            iss >> s ;
            smap[s] = i ;
          }
      
          entitySet dom = mix.domain() ;

          string sp = string(var_name.c_str()+1) ;
          map<string,int>::const_iterator mi = smap.find(sp) ;
          if(mi == smap.end()) {
            cerr << "warning, species " << sp << " does not exist in dataset!"
                 << endl ;
          } else {
            const int ind = mi->second ;
            int c = 0 ;
            scalar.allocate(dom) ;
            FORALL(dom,nd) {
              scalar[c++] = mix[nd][ind] ;
            } ENDFORALL ;
          }
          dom = scalar.domain() ;
          if(it == start_iter) {
            mean.allocate(dom) ;
            M2.allocate(dom) ;
            FORALL(dom,ii) {
              mean[ii] = 0. ;
              M2[ii] = 0. ;
            } ENDFORALL ;
          }
        
          FORALL(dom,nd) {
            double delta = scalar[nd] - mean[nd] ;
            mean[nd] += delta/n ;
            M2[nd] += delta*(scalar[nd]-mean[nd]) ;
          } ENDFORALL ;
        }
        store<float> m ;
        store<float> v ;
        entitySet dom = mean.domain() ;
        m.allocate(dom) ;
        v.allocate(dom) ;
        FORALL(dom,ii) {
          m[ii] = mean[ii] ;
          v[ii] = M2[ii]/(n) ;
        } ENDFORALL ;
        string filename = output_dir+'/' + var + "Mean_sca."
          + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"Mean" ;
	cout << "Writing Variables: " << sname ;
        fact_db facts ;
        Loci::writeContainer(file_id,sname,m.Rep(),facts) ;

        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "Var_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Var" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,v.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	cout << endl ;
      }
      break ;
    case BOUNDARY_SCALAR:
      {
        vector<double> mean ;
        vector<double> M2 ;
        vector<int> ids ;
        map<int,int> id2local ;
        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          
          n = n + 1.0 ;
          string filename = output_dir+'/' + var + "_bnd." +
            string(buf) + postfix + "_" + casename ;
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT, MPI_COMM_WORLD) ;
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
          
          if(it == start_iter) {
            vector<int> elemIds(nbel) ;
            readElementType(di,"entityIds",elemIds) ;
            ids.swap(elemIds) ;
            for(size_t i=0;i<ids.size();++i)
              id2local[ids[i]] = i ;
            vector<double> tmp1(nbel,0.0) ;
            mean.swap(tmp1) ;
            vector<double> tmp2(nbel,0.0) ;
            M2.swap(tmp2) ;
          }
        
          H5Gclose(di) ;
          vector<float> bvar(nbel) ;
          readElementType(file_id,variables[i].c_str(),bvar) ;
#ifdef H5_USE_16_API
          di = H5Gopen(file_id,"dataInfo") ;
#else
          di = H5Gopen(file_id,"dataInfo",H5P_DEFAULT) ;
#endif
          vector<int> eid(nbel,-1) ;
          readElementType(di,"entityIds",eid) ;
          H5Fclose(file_id) ;
        
          for(int i=0;i<nbel;++i) {
            map<int,int>::const_iterator mi ;
            if((mi=id2local.find(eid[i])) == id2local.end()) {
              cerr << "problems with ids in boundaries not matching" << endl ;
              cerr << "i=" << i << ",eid=" << eid[i] << endl ;
              exit(-1) ;
            }
            int lid = mi->second ;
            double delta = bvar[i] - mean[lid] ;
            mean[lid] += delta/n ;
            M2[lid] += delta*(bvar[i]-mean[lid]) ;
          }
        }
        vector<float> m(mean.size()) ;
        vector<float> v(mean.size()) ;
        for(size_t i=0;i<mean.size();++i) {
          m[i] = mean[i] ;
          v[i] = M2[i]/(n) ;
        }
        string filename = output_dir+'/' + var + "Mean_bnd."
          + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"Mean" ;
	cout << "Writing Variables: " << sname ;
#ifdef H5_USE_16_API
        hid_t group_id = H5Gcreate(file_id,"dataInfo",0) ;
#else
        hid_t group_id = H5Gcreate(file_id,"dataInfo",
				   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
        writeElementType(group_id,"entityIds",ids) ;
        H5Gclose(group_id) ;
        writeElementType(file_id,sname.c_str(),m) ;
        H5Fclose(file_id) ;

        filename = output_dir+'/' + var + "Var_bnd."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Var" ;
	cout << ", " << sname ;
#ifdef H5_USE_16_API
        group_id = H5Gcreate(file_id,"dataInfo",0) ;
#else
        group_id = H5Gcreate(file_id,"dataInfo",
			     H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
        writeElementType(group_id,"entityIds",ids) ;
        H5Gclose(group_id) ;
        writeElementType(file_id,sname.c_str(),v) ;
        H5Fclose(file_id) ;
	cout << endl ;
      }
      break;
    case BOUNDARY_VECTOR:
      {
        vector<vector3d<double> > mean ;
        vector<vector3d<double> > M2 ;
        vector<int> ids ;
        map<int,int> id2local ;
        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          
          n = n + 1.0 ;
          string filename = output_dir+'/' + var + "_bndvec." +
            string(buf) + postfix + "_" + casename ;
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT, MPI_COMM_WORLD) ;
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

          if(it == start_iter) {
            vector<int> elemIds(nbel) ;
            readElementType(di,"entityIds",elemIds) ;
            ids.swap(elemIds) ;
            for(size_t i=0;i<ids.size();++i)
              id2local[ids[i]] = i ;
            vector<vector3d<double> > tmp1(nbel,vector3d<double>(0.,0.,0.)) ;
            mean.swap(tmp1) ;
            vector<vector3d<double> > tmp2(nbel,vector3d<double>(0.,0.,0.)) ;
            M2.swap(tmp2) ;
          }
        
          H5Gclose(di) ;
          vector<vector3d<float> > bvar(nbel) ;
          readElementType(file_id,variables[i].c_str(),bvar) ;
#ifdef H5_USE_16_API
          di = H5Gopen(file_id,"dataInfo") ;
#else
          di = H5Gopen(file_id,"dataInfo",H5P_DEFAULT) ;
#endif
          vector<int> eid(nbel,-1) ;
          readElementType(di,"entityIds",eid) ;
          H5Fclose(file_id) ;

          for(int i=0;i<nbel;++i) {
            map<int,int>::const_iterator mi ;
            if((mi=id2local.find(eid[i])) == id2local.end()) {
              cerr << "problems with ids in boundaries not matching" << endl ;
              exit(-1) ;
            }
            int lid = mi->second ;
            
            vector3d<double> delta = vector3d<double>(bvar[i].x - mean[lid].x,
                                                      bvar[i].y - mean[lid].y,
                                                      bvar[i].z - mean[lid].z) ;
            mean[lid] += (1./n)*delta ;
            M2[lid].x += delta.x*(bvar[i].x-mean[lid].x) ;
            M2[lid].y += delta.y*(bvar[i].y-mean[lid].y) ;
            M2[lid].z += delta.z*(bvar[i].z-mean[lid].z) ;
          }
        }
        vector<vector3d<float> > m(mean.size()) ;
        vector<vector3d<float> > v(mean.size()) ;
        for(size_t i=0;i<mean.size();++i) {
          m[i] = vector3d<float>(mean[i].x, mean[i].y, mean[i].z) ;
          v[i] = vector3d<float>(M2[i].x,M2[i].y,M2[i].z)/(n) ;
        }
        string filename = output_dir+'/' + var + "Mean_bndvec."
          + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"Mean" ;
	cout << "Writing Variables: " << sname ;
#ifdef H5_USE_16_API
        hid_t group_id = H5Gcreate(file_id,"dataInfo",0) ;
#else
        hid_t group_id = H5Gcreate(file_id,"dataInfo",
				   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
        writeElementType(group_id,"entityIds",ids) ;
        H5Gclose(group_id) ;
        writeElementType(file_id,sname.c_str(),m) ;
        H5Fclose(file_id) ;

        filename = output_dir+'/' + var + "Var_bndvec."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"Var" ;
	cout << ", " << sname ;

#ifdef H5_USE_16_API
        group_id = H5Gcreate(file_id,"dataInfo",0) ;
#else
        group_id = H5Gcreate(file_id,"dataInfo",
			     H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
#endif
        writeElementType(group_id,"entityIds",ids) ;
        H5Gclose(group_id) ;
        writeElementType(file_id,sname.c_str(),v) ;
        H5Fclose(file_id) ;
	cout << endl ;
      }
      break ;
    case BOUNDARY_DERIVED_SCALAR:
    case BOUNDARY_DERIVED_VECTOR:
      cerr << "variable " << variables[i] << " ignored for Mean operation."
           << endl ;
      break ;
    default:
      cerr << "unable to process variable " << var_name << "!"<< endl ;
      cerr << "this variable is ignored!" << endl ;
      break ;
    }
  }
}

void combine_mean(string casename, string iteration,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames,
                  int end_iter, int inc_iter,
		  bool do_favre) {

  string postfix ="" ;
  string iter_part = "" ;
  size_t i = 0;
  while(iteration[i] >= '0' && iteration[i] <= '9' && i<iteration.size()) {
    iter_part += iteration[i] ;
    i++ ;
  }
  while(i<iteration.size()) {
    postfix += iteration[i] ;
    i++ ;
  }
  int start_iter = atoi(iter_part.c_str()) ;
  
  for(size_t i=0;i<variables.size();++i) {
    int ssz = variables[i].size()  ;
    string var_name ;
    for(int j=0;j<ssz-4;++j)
      var_name += variables[i][j] ;
    
    string rem ;
    for(int j=max(ssz-4,0);j<ssz;++j)
      rem += variables[i][j] ;
    if(rem != "mean") {
      cerr << "error variable " << variables[i] << " should end with 'mean' for combine to work properly!  Ignoring this variable!" << endl ;
      continue ;
    }
    
    cout << "processing variable: " << var_name << endl ;
    string var = var_name ;
    string filename = variable_filenames[i] ;
    switch(variable_types[i]) {
    case NODAL_SCALAR:
      {
	store<double> rbar ;
        store<double> rxmean ;
	store<double> rxvar ;
        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          n = n + 1.0 ;
          fact_db facts ;
	  store<float> scalar_mean ;
	  { 
	    string vname = var+"mean" ;
	    string filename = output_dir+'/'  + vname+ "_sca." +
	      string(buf) + postfix + "_" + casename ;
	    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
						H5F_ACC_RDONLY,
						H5P_DEFAULT, MPI_COMM_WORLD) ;
	    readData(file_id,vname,scalar_mean.Rep(),EMPTY,facts) ;
	    Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	  }
	  entitySet dom = scalar_mean.domain() ;

	  store<float> scalar_var ;
	  { 
	    string vname = var+"var" ;
	    string filename = output_dir+'/'  + vname+ "_sca." +
	      string(buf) + postfix + "_" + casename ;
	    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
						H5F_ACC_RDONLY,
						H5P_DEFAULT, MPI_COMM_WORLD) ;
	    readData(file_id,vname,scalar_var.Rep(),EMPTY,facts) ;
	    Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	  }

	  store<float> rentry ;
	  if(do_favre) { 
	    string vname = "rmean" ;
	    string filename = output_dir+'/'  + vname+ "_sca." +
	      string(buf) + postfix + "_" + casename ;
	    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
						H5F_ACC_RDONLY,
						H5P_DEFAULT, MPI_COMM_WORLD) ;
	    readData(file_id,vname,rentry.Rep(),EMPTY,facts) ;
	    Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	  } else {
	    rentry.allocate(dom) ;
	    FORALL(dom,ii) {
	      rentry[ii] = 1.0 ;
	    } ENDFORALL ;
	  }

          if(it == start_iter) {
	    rbar.allocate(dom) ;
	    rxmean.allocate(dom) ;
	    rxvar.allocate(dom) ;
            FORALL(dom,ii) {
              rbar[ii] = 0. ;
              rxmean[ii] = 0. ;
	      rxvar[ii] = 0. ;
            } ENDFORALL ;
          }

          FORALL(dom,nd) {
            double rdelta = rentry[nd] - rbar[nd] ;
	    double rbar_old = rbar[nd] ;
	    rbar[nd] += rdelta/n ;
	    double rxmean_old = rxmean[nd] ;
	    double rxdelta = rentry[nd]*scalar_mean[nd] - rxmean[nd] ;
	    rxmean[nd] += rxdelta/n ;
	    if(it!=start_iter) {
	      double xtilde_new = rxmean[nd]/rbar[nd] ;
	      // Correct variances for change in mean
	      double eps1 = (xtilde_new-rxmean_old/rbar_old) ;
	      double rxvarcor = rxvar[nd] + rbar_old*eps1*eps1 ;
	      double eps2 = (xtilde_new-scalar_mean[nd]) ;
	      double rvarcor = scalar_var[nd]*rentry[nd] + rentry[nd]*eps2*eps2 ;
	      // Update corrected variances using recursive delta form
	      rxvar[nd] = rxvarcor + (rvarcor-rxvarcor)/n ;
	    } else {
	      rxvar[nd] = rentry[nd]*scalar_var[nd] ;
	    }
          } ENDFORALL ;
        }
        store<float> m ;
        store<float> v ;
        entitySet dom = rxmean.domain() ;
        m.allocate(dom) ;
        v.allocate(dom) ;
        FORALL(dom,ii) {
	  cout << "rbar=" << rbar[ii] << endl ;
          m[ii] = rxmean[ii]/rbar[ii] ;
          v[ii] = rxvar[ii]/rbar[ii] ;
        } ENDFORALL ;
        string filename = output_dir+'/' + var + "CMean_sca."
          + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"CMean" ;
	cout << "Writing Variables: " << sname ;
        fact_db facts ;
        Loci::writeContainer(file_id,sname,m.Rep(),facts) ;

        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "CVar_sca."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"CVar" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,v.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	cout << endl ;
      }
      break;
    case NODAL_DERIVED:
      cerr << "variable " << variables[i] << " not supported NODAL_DERIVED variable for combine" << endl ;
      break ;
    case NODAL_VECTOR:
      {
	store<double> rbar ;
        store<vector3d<double> > rxmean ;
	store<vector3d<double> > rxvar ;
        double n = 0;
        for(int it = start_iter;it<=end_iter;it+=inc_iter) {
          char buf[128] ;
	  bzero(buf,128) ;
          snprintf(buf,127,"%d",it) ;
          cout << "iteration: " << buf << postfix << '\r' ;
	  cout.flush() ;
          n = n + 1.0 ;
          fact_db facts ;
	  store<vector3d<float> > vect_mean ;
	  { 
	    string vname = var+"mean" ;
	    string filename = output_dir+'/'  + vname+ "_vec." +
	      string(buf) + postfix + "_" + casename ;
	    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
						H5F_ACC_RDONLY,
						H5P_DEFAULT, MPI_COMM_WORLD) ;
	    readData(file_id,vname,vect_mean.Rep(),EMPTY,facts) ;
	    Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	  }
	  entitySet dom = vect_mean.domain() ;

	  store<vector3d<float> > vect_var ;
	  { 
	    string vname = var+"var" ;
	    string filename = output_dir+'/'  + vname+ "_vec." +
	      string(buf) + postfix + "_" + casename ;
	    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
						H5F_ACC_RDONLY,
						H5P_DEFAULT, MPI_COMM_WORLD) ;
	    readData(file_id,vname,vect_var.Rep(),EMPTY,facts) ;
	    Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	  }

	  store<float> rentry ;
	  if(do_favre) { 
	    string vname = "rmean" ;
	    string filename = output_dir+'/'  + vname+ "_sca." +
	      string(buf) + postfix + "_" + casename ;
	    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
						H5F_ACC_RDONLY,
						H5P_DEFAULT, MPI_COMM_WORLD) ;
	    readData(file_id,vname,rentry.Rep(),EMPTY,facts) ;
	    Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	  } else {
	    rentry.allocate(dom) ;
	    FORALL(dom,ii) {
	      rentry[ii] = 1.0 ;
	    } ENDFORALL ;
	  }

          if(it == start_iter) {
	    rbar.allocate(dom) ;
	    rxmean.allocate(dom) ;
	    rxvar.allocate(dom) ;
            FORALL(dom,ii) {
              rbar[ii] = 0. ;
              rxmean[ii] = vector3d<double>(0.,0.,0.) ;
	      rxvar[ii] = vector3d<double>(0.,0.,0.) ;
            } ENDFORALL ;
          }

          FORALL(dom,nd) {
            double rdelta = rentry[nd] - rbar[nd] ;
	    double rbar_old = rbar[nd] ;
	    rbar[nd] += rdelta/n ;
	    vector3d<double> rxmean_old = rxmean[nd] ;
	    vector3d<double> rxdelta = double(rentry[nd])*vector3d<double>(vect_mean[nd]) - rxmean[nd] ;
	    rxmean[nd] += rxdelta/n ;
	    if(it!=start_iter) {
	      vector3d<double> xtilde_new = rxmean[nd]/rbar[nd] ;
	      // Correct variances for change in mean
	      vector3d<double> eps1 = (xtilde_new-rxmean_old/rbar_old) ;
	      vector3d<double> eps1sq = vector3d<double>(eps1.x*eps1.x,
							 eps1.y*eps1.y,
							 eps1.z*eps1.z) ;
	      vector3d<double> rxvarcor = rxvar[nd] + rbar_old*eps1sq ;
	      vector3d<double> eps2 = (xtilde_new-vector3d<double>(vect_mean[nd])) ;
	      vector3d<double> eps2sq = vector3d<double>(eps2.x*eps2.x,
							 eps2.y*eps2.y,
							 eps2.z*eps2.z) ;
	      vector3d<double> rvarcor = vector3d<double>(vect_var[nd]*rentry[nd]) + rentry[nd]*eps2sq ;
	    // Update corrected variances using recursive delta form
	      rxvar[nd] = rxvarcor + (rvarcor-rxvarcor)/n ;
	    } else {
	      rxvar[nd] = rentry[nd]*vect_var[nd] ;
	    }
	  } ENDFORALL ;
        }
        store<vector3d<float> > m ;
        store<vector3d<float> > v ;
        entitySet dom = rxmean.domain() ;
        m.allocate(dom) ;
        v.allocate(dom) ;
        FORALL(dom,ii) {
          m[ii] = rxmean[ii]/rbar[ii] ;
          v[ii] = rxvar[ii]/rbar[ii] ;
        } ENDFORALL ;
        string filename = output_dir+'/' + var + "CMean_vec."
          + iteration + "_" + casename ;
        hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        string sname = var+"CMean" ;
	cout << "Writing Variables: " << sname ;
        fact_db facts ;
        Loci::writeContainer(file_id,sname,m.Rep(),facts) ;

        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;

        filename = output_dir+'/' + var + "CVar_vec."
          + iteration + "_" + casename ;
        file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT, MPI_COMM_WORLD) ;
        sname = var+"CVar" ;
	cout << ", " << sname ;
        Loci::writeContainer(file_id,sname,v.Rep(),facts) ;
        Loci::hdf5CloseFile(file_id, MPI_COMM_WORLD) ;
	cout << endl ;
      }
      break;
    case NODAL_MASSFRACTION:
      cerr << "variable " << variables[i] << "not supporting nodal mass fraction for combine" << endl ;
      break ;
    case BOUNDARY_SCALAR:
      cerr << "variable " << variables[i] << " not supprting boundary scalar for combine!" << endl ;
      break;
    case BOUNDARY_VECTOR:
      cerr << "varaible " << variables[i] << " not supporting bouddary vector for combine!" << endl ;
      break ;
    case BOUNDARY_DERIVED_SCALAR:
    case BOUNDARY_DERIVED_VECTOR:
      cerr << "variable " << variables[i] << " ignored for Mean operation."
           << endl ;
      break ;
    default:
      cerr << "unable to process variable " << var_name << "!"<< endl ;
      cerr << "this variable is ignored!" << endl ;
      break ;
    }
  }
}

