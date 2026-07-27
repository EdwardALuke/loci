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
#include <set> 
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
using std::ifstream ;

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "extract.h"

#pragma  GCC diagnostic ignored "-Wunused-variable"

string output_dir ;
#define MAX_NAME 1024
string getVarNameFromFile(hid_t file_id, string varName) {
  string name = varName ;
  hid_t grp = H5Gopen(file_id,"/", H5P_DEFAULT);
  hsize_t nobj ;
  H5Gget_num_objs(grp,&nobj) ;
  if(nobj == 1) {
    char memb_name[MAX_NAME] ;
    ssize_t len = H5Gget_objname_by_idx(grp,(hsize_t)0,
					memb_name,(size_t)MAX_NAME) ;
     name = string(memb_name) ;
  }
  H5Gclose(grp) ;
  return name ;
}


void readData(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet, fact_db &facts) {
  hid_t grp = H5Gopen(file_id,"/", H5P_DEFAULT);
  hsize_t nobj ;
  H5Gget_num_objs(grp,&nobj) ;
  if(nobj == 1) {
    char memb_name[MAX_NAME] ;
    ssize_t len = H5Gget_objname_by_idx(grp,(hsize_t)0,
					memb_name,(size_t)MAX_NAME) ;
    string name(memb_name) ;

#ifdef VERBOSE
    if(name != vname) {
      cerr << "NOTE: reading dataset '" << name << "' instead of '" << vname << "'" << endl ;
    }
#endif
    H5Gclose(grp) ;
    Loci::readContainer(file_id,name,var,readSet,facts) ;
  } else {
    H5Gclose(grp) ;
    Loci::readContainer(file_id,vname,var,readSet,facts) ;
  }
}
    

void Usage(int ac, char *av[]) {
  cerr << av[0] << ": Incorrect Usage" << endl ;
  cout << "Usage:" << endl;
  cout << av[0] << " <package> [package options] <case_name> <time_step> [variable(s)]" << endl ;
  cout << endl ;
  cout << "where <package> may be:" << endl
       << "-2d :  extract for the 2dgv plotting package" << endl
       << "-fv :  extract for the FieldView post-processing package" << endl
       << "-en :  extract for the Ensight post-processing package" << endl
       << "-en_with_id :  extract for the Ensight post-processing package with node id and element id" << endl
    //#ifdef HAVE_CGNS
       << "-cgns :  extract for the CGNS post-processing package" << endl
    //#endif
       << "-tec:  extract for the TecPlot post-procesing package" << endl
       << "-vtk:   extract for the Paraview post-procesing package" << endl
       << "-vtk64: extract for the Paraview post-procesing package (for large cases, must use >= Paraview 3.98)" << endl
       << "-vtk_surf:  extract boundary surface mesh for the Paraview post-procesing package" << endl
       << "-vtk_surf64:  extract boundary surface mesh for the Paraview post-procesing package (for large cases, must use >= Paraview 3.98)" << endl
       << "-ascii: extract to an ascii file" << endl
       << "-surf: extract boundary surface mesh" << endl
       << "-cut:  extract a cutting plane for the 2dgv plotting package" << endl
       << "-mean: generate mean and variance from a family of ouput variables" << endl 
       << "-combine: combine mean and variance from online averaging ouput" << endl 
       << "-fcombine: combine favre mean and variance from online averaging ouput" << endl 
       << endl ;
  cout << "If no variables are specified, extract scans for all available variables"
       << " for the requested case and time step." << endl ;
  cout << "Variables are defined by the solver, but typically include: " << endl
       << "r     - nodal density" << endl 
       << "p     - nodal log10 pressure" << endl
       << "P     - nodal absolute pressure" << endl
       << "pg    - nodal gage pressure" << endl 
       << "u     - nodal velocity magnitude" << endl
       << "m     - nodal mach number" << endl
       << "t     - nodal temperature" << endl
       << "a     - nodal soundspeed" << endl
    //       << "G - extract Gamma" << endl
    //       << "R - extract Gas Constant R~" << endl
       << "f<species> - nodal mass fractions for <species>" << endl
       << "v     - nodal velocity vector"<<endl
       << "x     - x coordinate" << endl 
       << "y     - y coordinate" << endl 
       << "z     - z coordinate" << endl
       << "Boundary Variables:" << endl 
       << "qdot  - wall boundary heat flux" << endl
       << "yplus - wall boundary y plus" << endl
       << "tau   - viscous wall shear stress vector" << endl
       << "tw    - viscous wall boundary temperature" << endl
       << "pw    - viscous wall boundary pressure" << endl
       << "n     - boundary normal (-ascii only)" << endl
       << "area  - boundary area   (-ascii only)" << endl
       << endl ;

  cout << "extra options for particle extraction" << endl ;
  cout << "  -mp <n> : maximum particles to extract" << endl ;
  cout << "          : default (if not specified) is all particles available"
       << endl ;

  cout << "extra options for the 2dgv postprocessing package" << endl 
       << "  -bc <boundary_tag>  : specify which boundary to extract for (multiple -bc's"
       << endl 
       << "             are merged)" << endl 
       << "  -xy : project boundary onto z=0 plane (default)" << endl
       << "  -yz : project boundary onto x=0 plane" << endl
       << "  -xz : project boundary onto y=0 plane" << endl
       << "  -xr : project boundary onto x,radius plane (for axisymmetric grids)"
       << endl << endl ;
  cout << "extra options for orienting cutting plane" << endl
       << "  -xy : originate transformations from z=0 plane (default)" << endl
       << "  -yz : originate transformations from x=0 plane" << endl
       << "  -xz : originate transformations from y=0 plane" << endl
       << "  -Rx <amount> : rotate cutting plane about x-axis" << endl
       << "  -Ry <amount> : rotate cutting plane about y-axis" << endl
       << "  -Rz <amount> : rotate cutting plane about z-axis" << endl
       << "  -Sx <amount> : translate cutting plane along x-axis" << endl
       << "  -Sy <amount> : translate cutting plane along y-axis" << endl
       << "  -Sz <amount> : translate cutting plane along z-axis" << endl << endl;
  cout << "extra options for averaging feature '-mean'" << endl 
       << "  -end <value> : ending iteration number for averaging" << endl
       << "  -inc <value> : value to increment between iterations for averaging" << endl 
       << endl ;
  cout << "extra options for controlling extract" << endl
       << "   -dir <directory> : change extract directory from default 'output'" << endl
       << "   -skipPartDirs    : do not process the data in the part sub-directories, if present" << endl
       << "                    : default is to process the data in the part sub-directories, if present"

       << endl << endl ;
  cout << "example:  to extract OH species from time step 50 of ssme simulation for visualization with 2dgv use:" << endl
       << av[0] << " -2d -bc 1 -xr ssme 50 fOH" << endl ;
  cout << "example: to extract an ascii table of boundary heat flux and x locations:"<<endl
       << av[0] << " -ascii -bc 4 nozzle 0 x qdot" << endl ;
  cout << "example: to extract a cutting plane with various transformations:" << endl
       << av[0] << " -cut -xz -Sy 1.5 -Rx 30 -Rz -15 nozzle 0 P t" << endl;
  cout << "example: to extract 5000 particles with associated temperature"
       << " from time step 100 for visualization with Ensight:" << endl
       << av[0] << " -en combustor 100 particle_temp -mp 5000" << endl ;

  cout << "example: to compute mean and variance values for velocity and temperature" 
       << " for iterations 1000-3000 output every 100 iterations run:"
       << endl
       << av[0] << " -mean -end 3000 -inc 100 nozzle 1000 t v" << endl
       << "NOTE: outputs variabes tMean, tVar, vMean, vVar, vCuv, vCuw, vCvw"
       << "      at iteration 1000." << endl ;
  exit(-1) ;
}

size_t  sizeElementType(hid_t group_id, const char *element_name) {
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
  return size ;
  
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

string getTopoFileName(string output_dir, string casename, string iteration) {
  string gridtopo = output_dir+"/" + casename +".topo" ;
  string toponamefile = output_dir + "/topo_file." + iteration + "_" + casename ;
  struct stat tmpstat ;
  if(stat(toponamefile.c_str(),&tmpstat)== 0) {
    ifstream tinput(toponamefile.c_str()) ;
    string name  ;
    tinput >> name ;
    name = output_dir + "/" + name ;
    if(stat(name.c_str(),&tmpstat)==0)
      gridtopo=name ;
  }
  return gridtopo ;
}

volumePart::volumePart(string out_dir, string iteration, string casename,
		       vector<string> vars) {
  error = true ;
  partName = "Volume" ;

  // Check number of nodes
  //-------------------------------------------------------------------------
  posFile = getPosFile(out_dir,iteration,casename) ;
  hid_t file_id = H5Fopen(posFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0)
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"pos") ;
#else
  hid_t elg = H5Gopen(file_id,"pos",H5P_DEFAULT) ;
#endif
  if(elg < 0) {
    H5Fclose(file_id) ;
    return ;
  }
  nnodes = sizeElementType(elg,"data") ;
  H5Gclose(elg) ;
  H5Fclose(file_id) ;
    

  // Check for iblank information
  //-------------------------------------------------------------------------

  string iblankname = out_dir+"/grid_iblank." + iteration + "_" + casename ;
  struct stat tmpstat ;
  has_iblank = false ;
  if(stat(iblankname.c_str(),&tmpstat)== 0) {
    file_id = Loci::hdf5OpenFile(iblankname.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT) ;
    if(file_id < 0) {
      return ;
    }
    fact_db facts ;
    store<unsigned char> iblank_tmp ;
    readData(file_id,"iblank",iblank_tmp.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;
    entitySet pdom = interval(1,nnodes) ;
    iblank.allocate(pdom) ;
    entitySet dom = iblank_tmp.domain() ;
    int cnt = 1 ;
    FORALL(dom,nd) {
      iblank[cnt++] = iblank_tmp[nd] ;
    } ENDFORALL ;
    has_iblank = true ;
  } 

  // Check for element types in topo file
  //-------------------------------------------------------------------------
  topoFile = getTopoFileName(out_dir, casename, iteration) ;
  file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  elg = H5Gopen(file_id,"elements") ;
#else
  elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;

  ntets = sizeElementType(elg,"tetrahedra") ;
  nhexs = sizeElementType(elg,"hexahedra") ;
  nprsm = sizeElementType(elg,"prism") ;
  npyrm = sizeElementType(elg,"pyramid") ;
  ngenc = sizeElementType(elg,"GeneralCellNfaces") ;
  
  size_t ntets_b = ntets ;
  size_t nhexs_b = nhexs ;
  size_t nprsm_b = nprsm ;
  size_t npyrm_b = npyrm ;
  size_t ngenc_b = ngenc ;
  const int block_size=65536 ; // Size of blocking factor
  if(has_iblank) {
    // need to adjust the number of elements based on iblanking.
    if(ntets > 0) {
      int nblocks = (ntets-1)/block_size+1 ;
      int remain = ntets ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        vector<Array<int,4> > tets(size) ;
        readElementTypeBlock(elg,"tetrahedra",tets,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<4;++j)
            if(iblank[tets[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
	  else
	    tetsIblanked += start-size+i ;
        }
      }
      WARN(remain != 0) ;
      ntets_b = cnt ;
      if(ntets-ntets_b > 0)
        cout << ntets-ntets_b << " tetrahedra iblanked" << endl ;
    }
    if(nhexs > 0) {
      int nblocks = (nhexs-1)/block_size+1 ;
      int remain = nhexs ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        
        vector<Array<int,8> > hexs(size) ;
        readElementTypeBlock(elg,"hexahedra",hexs,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<8;++j)
            if(iblank[hexs[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
	  else 
	    hexsIblanked += start-size+i ;
        }
      }
      WARN(remain != 0) ;
      nhexs_b = cnt ;
      if(nhexs-nhexs_b > 0)
        cout << nhexs-nhexs_b << " hexahedra iblanked" << endl ;
    }
    if(nprsm > 0) {
      int nblocks = (nprsm-1)/block_size+1 ;
      int remain = nprsm ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        vector<Array<int,6> > prsm(size) ;
        readElementTypeBlock(elg,"prism",prsm,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<6;++j)
            if(iblank[prsm[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
	  else
	    prsmIblanked += start-size+i ;
        }
      }
      WARN(remain != 0) ;
      nprsm_b = cnt ;
      if(nprsm-nprsm_b > 0)
        cout << nprsm-nprsm_b << " prisms iblanked" << endl ;
        
    }
    if(npyrm > 0) {
      int nblocks = (npyrm-1)/block_size+1 ;
      int remain = npyrm ;
      int start = 0 ;
      int cnt = 0 ;
      for(int b=0;b<nblocks;++b) {
        int size = min(block_size,remain) ;
        vector<Array<int,5> > pyrm(size) ;
        readElementTypeBlock(elg,"pyramid",pyrm,start,size) ;
        remain -= size ;
        start += size ;
        for(int i=0;i<size;++i) {
          bool blank = true ;
          for(int j=0;j<5;++j)
            if(iblank[pyrm[i][j]] < 2)
              blank = false ;
          if(!blank)
            cnt++ ;
	  else 
	    pyrmIblanked += start-size+i ;
        }
      }
      WARN(remain != 0) ;
      npyrm_b = cnt ;
      if(npyrm-npyrm_b > 0)
        cout << npyrm-npyrm_b << " pyramids iblanked" << endl ;
    }
    if(ngenc > 0) {
      vector<int> GeneralCellNfaces(ngenc) ;
      readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces) ;
      size_t nside = sizeElementType(elg,"GeneralCellNsides") ;
      vector<int> GeneralCellNsides(nside) ;
      readElementType(elg,"GeneralCellNsides",GeneralCellNsides) ;
      size_t nnodes = sizeElementType(elg,"GeneralCellNodes") ;
      vector<int> GeneralCellNodes(nnodes) ;
      readElementType(elg,"GeneralCellNodes",GeneralCellNodes) ;
      int cnt1 = 0 ;
      int cnt2 = 0 ;
      int cnt = 0 ;
      for(size_t i=0;i<ngenc;++i) {
        bool blank = true ;
        int nf = GeneralCellNfaces[i] ;
        for(int f=0;f<nf;++f) {
          int fs = GeneralCellNsides[cnt1++] ;
          for(int n=0;n<fs;++n) {
            int nd = GeneralCellNodes[cnt2++] ;
            if(iblank[nd] < 2)
              blank = false ;
          }
        }
        if(!blank)
          cnt++ ;
        else
          gencIblanked += (int)i ;
      }
      ngenc_b = cnt ;
      if(ngenc-ngenc_b > 0)
        cout << ngenc-ngenc_b << " general cells iblanked" << endl ;
    }
    H5Gclose(elg) ;
    Loci::hdf5CloseFile(file_id) ;
     
  }
  ntets_orig = ntets ;
  nhexs_orig = nhexs ;
  nprsm_orig = nprsm ;
  npyrm_orig = npyrm ;
  ngenc_orig = ngenc ;

  ntets = ntets_b ;
  nhexs = nhexs_b ;
  nprsm = nprsm_b ;
  npyrm = npyrm_b ;
  ngenc = ngenc_b ;
  ntetsIblank = ntets_orig-ntets ;
  nhexsIblank = nhexs_orig-nhexs ;
  nprsmIblank = nprsm_orig-nprsm ;
  npyrmIblank = npyrm_orig-npyrm ;
  ngencIblank = ngenc_orig-ngenc ;
  
  // Now check for variables
  for(size_t i=0;i<vars.size();++i) {
    string varname = vars[i] ;
    string filename = out_dir+'/' + varname + "_sca." + iteration + "_" + casename ;
    struct stat tmpstat ;
    if(stat(filename.c_str(),&tmpstat) == 0) {
      nodalScalarVars[varname] = filename ;
    } else {
      filename = out_dir+'/' + varname + "_vec." + iteration + "_" + casename ;
      if(stat(filename.c_str(),&tmpstat) == 0) {
        nodalVectorVars[varname] = filename ;
      } else {
      }
    }
  }
  
  error = false ;
}
bool volumePart::hasNodalScalarVar(string var) const {
  map<string,string>::const_iterator mi=nodalScalarVars.find(var) ;
  return (mi != nodalScalarVars.end()) ;
}
bool volumePart::hasNodalVectorVar(string var) const {
  map<string,string>::const_iterator mi=nodalVectorVars.find(var) ;
  return (mi != nodalVectorVars.end()) ;
}
std::vector<string> volumePart::getNodalScalarVars() const  {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=nodalScalarVars.begin();mi!=nodalScalarVars.end();++mi)
    tmp.push_back(mi->first) ;

  return tmp ;
}

std::vector<string> volumePart::getNodalVectorVars() const {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=nodalVectorVars.begin();mi!=nodalVectorVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}

void volumePart::getPos(vector<vector3d<float> > &pos) const {
  pos.clear() ;
  string filename = posFile ;
  hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0)
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"pos") ;
#else
  hid_t elg = H5Gopen(file_id,"pos",H5P_DEFAULT) ;
#endif
  if(elg < 0) {
    H5Fclose(file_id) ;
    return ;
  }
  size_t nsz = sizeElementType(elg,"data") ;
  if(nsz != nnodes) {
    H5Gclose(elg) ;
    H5Fclose(file_id) ;
    return ;
  }
    
  vector<vector3d<float> > tmp(nsz) ;
  pos.swap(tmp) ;
  readElementType(elg,"data",pos) ;
  H5Gclose(elg) ;
  H5Fclose(file_id) ;
}


void volumePart::getPos(vector<vector3d<double> > &pos) const {
  pos.clear() ;
  string filename = posFile ;
  hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0)
    return ;

#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"pos") ;
#else
  hid_t elg = H5Gopen(file_id,"pos",H5P_DEFAULT) ;
#endif
  
  if(elg < 0) {
    H5Fclose(file_id) ;
    return ;
  }
  size_t nsz = sizeElementType(elg,"data") ;
  if(nsz != nnodes) {
    H5Gclose(elg) ;
    H5Fclose(file_id) ;
    return ;
  }
    
  vector<vector3d<double> > tmp(nsz) ;
  pos.swap(tmp) ;
  readElementType(elg,"data",pos) ;
  H5Gclose(elg) ;
  H5Fclose(file_id) ;
}





void volumePart::getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const {
  tets.clear() ;
  if(ntets <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,ntets_orig-start) ;
  vector<Array<int,4> > tets_local(lsize) ;
  readElementTypeBlock(elg,"tetrahedra",tets_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = tetsIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    tets.swap(tets_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<Array<int,4> > tets_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    tets_new[cnt] = tets_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  tets.swap(tets_new) ;
}

void volumePart::getTetIds(vector<int> &tetids, size_t start, size_t size) const {
  tetids.clear() ;
  if(ntets <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,ntets_orig-start) ;
  vector<int > tets_local(lsize) ;
  readElementTypeBlock(elg,"tetrahedra_ids",tets_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = tetsIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    tetids.swap(tets_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<int > tets_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    tets_new[cnt] = tets_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  tetids.swap(tets_new) ;
}

void volumePart::getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const {
  pyrms.clear() ;
  if(npyrm <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,npyrm_orig-start) ;
  vector<Array<int,5> > pyrms_local(lsize) ;
  readElementTypeBlock(elg,"pyramid",pyrms_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = pyrmIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    pyrms.swap(pyrms_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<Array<int,5> > pyrms_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    pyrms_new[cnt] = pyrms_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  pyrms.swap(pyrms_new) ;
}

void volumePart::getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const {
  pyrmids.clear() ;
  if(npyrm <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,npyrm_orig-start) ;
  vector<int > pyrms_local(lsize) ;
  readElementTypeBlock(elg,"pyramid_ids",pyrms_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = pyrmIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    pyrmids.swap(pyrms_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<int > pyrms_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    pyrms_new[cnt] = pyrms_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  pyrmids.swap(pyrms_new) ;
}
void volumePart::getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const {
  prsms.clear() ;
  if(nprsm <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,nprsm_orig-start) ;
  vector<Array<int,6> > prsms_local(lsize) ;
  readElementTypeBlock(elg,"prism",prsms_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = prsmIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    prsms.swap(prsms_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<Array<int,6> > prsms_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    prsms_new[cnt] = prsms_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  prsms.swap(prsms_new) ;
}
void volumePart::getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const {
  prsmids.clear() ;
  if(nprsm <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,nprsm_orig-start) ;
  vector<int > prsms_local(lsize) ;
  readElementTypeBlock(elg,"prism_ids",prsms_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = prsmIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    prsmids.swap(prsms_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<int > prsms_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    prsms_new[cnt] = prsms_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  prsmids.swap(prsms_new) ;
}

void volumePart::getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const {
  hexs.clear() ;
  if(nhexs <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,nhexs_orig-start) ;
  vector<Array<int,8> > hexs_local(lsize) ;
  readElementTypeBlock(elg,"hexahedra",hexs_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = hexsIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    hexs.swap(hexs_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<Array<int,8> > hexs_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    hexs_new[cnt] = hexs_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  hexs.swap(hexs_new) ;
}
void volumePart::getHexIds(vector<int> &hexids, size_t start, size_t size) const {
  hexids.clear() ;
  if(nhexs <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = min(size,nhexs_orig-start) ;
  vector<int > hexs_local(lsize) ;
  readElementTypeBlock(elg,"hexahedra_ids",hexs_local,start,lsize) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = hexsIblanked & interval(start,start+lsize-1) ;
  if(iblank==EMPTY) {
    hexids.swap(hexs_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<int > hexs_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank<<start)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    hexs_new[cnt] = hexs_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  hexids.swap(hexs_new) ;
}
void volumePart::getGenCell(vector<int> &genCellNfaces, 
			    vector<int> &genCellNsides,
			    vector<int> &genCellNodes) const {
  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;

  vector<int> GeneralCellNfaces(ngenc_orig) ;
  readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces) ;
  size_t nside = sizeElementType(elg,"GeneralCellNsides") ;
  vector<int> GeneralCellNsides(nside) ;
  readElementType(elg,"GeneralCellNsides",GeneralCellNsides) ;
  size_t nnodes = sizeElementType(elg,"GeneralCellNodes") ;
  vector<int> GeneralCellNodes(nnodes) ;
  readElementType(elg,"GeneralCellNodes",GeneralCellNodes) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;
  // If no general cells iblanked, then return
  if(ngenc_orig == ngenc) {
    genCellNfaces.swap(GeneralCellNfaces) ;
    genCellNsides.swap(GeneralCellNsides) ;
    genCellNodes.swap(GeneralCellNodes) ;
    return ;
  }
  int currentFaceOffset = 0 ;
  int currentNodeOffset = 0 ;
  int skip_cells = 0;
  int skip_faces = 0 ;
  int skip_nodes = 0 ;
  for(size_t i=0;i<ngenc_orig;++i) {
    bool blank = gencIblanked.inSet(i) ;
    // nf is number of faces for this general cell
    int nf = GeneralCellNfaces[i] ;
    // nn is the number of nodes firthis general cell
    int nn = 0 ;
    for(int f=0;f<nf;++f) {
      nn += GeneralCellNsides[currentFaceOffset+f] ;
    }
    if(blank) {
      skip_cells += 1 ;
      skip_faces += nf ;
      skip_nodes += nn ;
    } else {
      if(skip_cells > 0) {
	GeneralCellNfaces[i-skip_cells] = GeneralCellNfaces[i] ;
	for(int f=0;f<nf;++f)
	  GeneralCellNsides[currentFaceOffset-skip_faces+f] =
	    GeneralCellNsides[currentFaceOffset+f] ;
	for(int n=0;n<nn;++n)
	  GeneralCellNodes[currentNodeOffset-skip_nodes+n] =
	    GeneralCellNodes[currentNodeOffset+n] ;
      }
    }
    currentFaceOffset += nf ;
    currentNodeOffset += nn ;
  }

  GeneralCellNfaces.resize(GeneralCellNfaces.size()-skip_cells) ;
  GeneralCellNsides.resize(GeneralCellNsides.size()-skip_faces) ;
  GeneralCellNodes.resize(GeneralCellNodes.size()-skip_nodes) ;
  genCellNfaces.swap(GeneralCellNfaces) ;
  genCellNsides.swap(GeneralCellNsides) ;
  genCellNodes.swap(GeneralCellNodes) ;
}
void volumePart::getGenIds(vector<int> &genids) const {
  genids.clear() ;
  if(ngenc <=0)
    return ;

  hid_t file_id = H5Fopen(topoFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,"elements") ;
#else
  hid_t elg = H5Gopen(file_id,"elements",H5P_DEFAULT) ;
#endif
  if(elg < 0) 
    return ;
  int lsize = ngenc_orig ;
  vector<int > gens_local(lsize) ;
  readElementType(elg,"GeneralCell_ids",gens_local) ;
  H5Gclose(elg) ;
  Loci::hdf5CloseFile(file_id) ;

  entitySet iblank = gencIblanked & interval(0,lsize-1) ;
  if(iblank==EMPTY) {
    genids.swap(gens_local) ;
    return ;
  }
  // Remove iblanked cells
  vector<int > gens_new(lsize-iblank.size()) ;
  int cnt = 0 ;
  entitySet dom = (~(iblank)) & interval(0,lsize-1) ;
  FORALL(dom,cp) {
    gens_new[cnt] = gens_local[cp] ;
    cnt++ ;
  }ENDFORALL ;
  genids.swap(gens_new) ;
}
  
void volumePart::getNodalScalar(string varname, vector<float> &vals) const {
  vals.clear() ;
  map<string,string>::const_iterator mi = nodalScalarVars.find(varname) ;
  if(mi == nodalScalarVars.end())
    return ;
  string filename = mi->second;
  hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0)
    return ;
  string vname = getVarNameFromFile(file_id,varname) ;
#ifdef VERBOSE
  if(vname != varname) {
    cerr << "reading var '" << vname << "' from file '" << filename << "'" << endl ;
  }
#endif
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,vname.c_str()) ;
#else
  hid_t elg = H5Gopen(file_id,vname.c_str(),H5P_DEFAULT) ;

#endif
  if(elg < 0)
    return ;
  int nsz = sizeElementType(elg,"data") ;
  vector<float> tmp(nsz) ;
  vals.swap(tmp) ;
  readElementType(elg,"data",vals) ;
  H5Gclose(elg) ;
  H5Fclose(file_id) ;
}

void volumePart::getNodalVector(string varname, vector<vector3d<float> > &vals) const {
  vals.clear() ;
  map<string,string>::const_iterator mi = nodalVectorVars.find(varname) ;
  if(mi == nodalVectorVars.end())
    return ;
  string filename = mi->second;
  hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0)
    return ;
  string vname = getVarNameFromFile(file_id,varname) ;
#ifdef VERBOSE
  if(vname != varname) {
    cerr << "reading var '" << vname << "' from file '" << filename << "'" << endl ;
  }
#endif
#ifdef H5_USE_16_API
  hid_t elg = H5Gopen(file_id,vname.c_str()) ;
#else
  hid_t elg = H5Gopen(file_id,vname.c_str(),H5P_DEFAULT) ;
#endif
  if(elg < 0)
    return ;
  int nsz = sizeElementType(elg,"data") ;
  vector<vector3d<float> > tmp(nsz) ;
  vals.swap(tmp) ;
  readElementType(elg,"data",vals) ;
  H5Gclose(elg) ;
  H5Fclose(file_id) ;
}

void volumePart::getNodalIblank(vector<unsigned char> &blank) const {
  if(has_iblank) {
    Loci::entitySet idom = iblank.domain() ;
    vector<unsigned char> tmp(idom.size()) ;
    int cnt = 0 ;
    FORALL(idom,ii) {
      tmp[cnt] = iblank[ii] ;
      cnt++ ;
    } ENDFORALL ;
    blank.swap(tmp) ;
  } else {
    blank.clear() ;
  } 
}
//----
void volumePartDerivedVars::processDerivedVars(const vector<string> &vars) {
  for(size_t i=0;i<vars.size();++i) {
    if(vars[i] == "m" && !shadowPart->hasNodalScalarVar("m")) {
      if(shadowPart->hasNodalScalarVar("a") && 
	 shadowPart->hasNodalVectorVar("v"))
	derivedVars["m"] = VAR_M ;
    }
    if(vars[i] == "P" && !shadowPart->hasNodalScalarVar("P")) {
      if(shadowPart->hasNodalScalarVar("pg")) {
	derivedVars["P"] = VAR_P ;
      }
    }
    if(vars[i] == "p" && !shadowPart->hasNodalScalarVar("p")) {
      if(shadowPart->hasNodalScalarVar("pg")) {
	derivedVars["p"] = VAR_logp ;
      }
    }
    if(vars[i] == "u" && !shadowPart->hasNodalScalarVar("u")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["u"] = VAR_U ;
      }
    }
    if(vars[i] == "0" && !shadowPart->hasNodalScalarVar("0")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["0"] = VAR_0 ;
      }
    }
    if(vars[i] == "1" && !shadowPart->hasNodalScalarVar("1")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["1"] = VAR_1 ;
      }
    }
    if(vars[i] == "2" && !shadowPart->hasNodalScalarVar("2")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["2"] = VAR_2 ;
      }
    }
    if(vars[i] == "x")
      derivedVars["x"] = VAR_X ;
    if(vars[i] == "y")
      derivedVars["y"] = VAR_Y ;
    if(vars[i] == "z")
      derivedVars["z"] = VAR_Z ;
  }
}
volumePartDerivedVars::volumePartDerivedVars(volumePartP part,
					     string output_dir, 
					     string casename, string iteration,
					     vector<string> vars) {
  error = part->fail() ;
  partName = part->getPartName() ;
  nnodes = part->getNumNodes() ;
  ntets = part->getNumTets() ;
  nhexs = part->getNumHexs() ;
  nprsm = part->getNumPrsm() ;
  npyrm = part->getNumPyrm() ;
  ngenc = part->getNumGenc() ;
  ntetsIblank = part->getNumTetsIblank() ;
  nhexsIblank = part->getNumHexsIblank() ;
  nprsmIblank = part->getNumPrsmIblank() ;
  npyrmIblank = part->getNumPyrmIblank() ;
  ngencIblank = part->getNumGencIblank() ;

  shadowPart = part ;

  string filename = output_dir+"/Pambient_par." + iteration +"_" + casename ;
  struct stat tmpstat ;
  if(stat(filename.c_str(),&tmpstat) == 0) {
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				       H5F_ACC_RDONLY,
				       H5P_DEFAULT) ;
    Pambient = 0 ;
    if(file_id >= 0) {
      fact_db facts ;
      param<float> Pamb ;
      readData(file_id,"Pambient",Pamb.Rep(),EMPTY,facts) ;
      Loci::hdf5CloseFile(file_id) ;
      Pambient = *Pamb ;
    } else {
      cerr << "Unable to open file " << filename << endl ;
    }
  }
	  
  processDerivedVars(vars) ;
}

bool volumePartDerivedVars::hasNodalScalarVar(string var) const {
  map<string,derivedVar_t>::const_iterator mi=derivedVars.find(var) ;
  if(mi==derivedVars.end())
    return shadowPart->hasNodalScalarVar(var) ;
  else
    return true ;
}
bool volumePartDerivedVars::hasNodalVectorVar(string var) const {
  return shadowPart->hasNodalVectorVar(var) ;
}
std::vector<string> volumePartDerivedVars::getNodalScalarVars() const  {
  
  vector<string> tmp = shadowPart->getNodalScalarVars() ;
  map<string,derivedVar_t>::const_iterator mi ;
  for(mi=derivedVars.begin();mi!=derivedVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}

std::vector<string> volumePartDerivedVars::getNodalVectorVars() const {
  return shadowPart->getNodalVectorVars() ;
}

void volumePartDerivedVars::getPos(vector<vector3d<float> > &pos) const {
  shadowPart->getPos(pos) ;
}
void volumePartDerivedVars::getPos(vector<vector3d<double> > &pos) const {
  shadowPart->getPos(pos) ;
}
void volumePartDerivedVars::getTetBlock(vector<Array<int,4> > &tets, size_t start, size_t size) const {
  shadowPart->getTetBlock(tets,start,size) ;
}

void volumePartDerivedVars::getTetIds(vector<int> &tetids, size_t start, size_t size) const {
  shadowPart->getTetIds(tetids,start,size) ;
}

void volumePartDerivedVars::getPyrmBlock(vector<Array<int,5> > &pyrms, size_t start, size_t size) const {
  shadowPart->getPyrmBlock(pyrms,start,size) ;
}

void volumePartDerivedVars::getPyrmIds(vector<int> &pyrmids, size_t start, size_t size) const {
  shadowPart->getPyrmIds(pyrmids,start,size) ;
}
void volumePartDerivedVars::getPrsmBlock(vector<Array<int,6> > &prsms, size_t start, size_t size) const {
  shadowPart->getPrsmBlock(prsms,start,size) ;
}
void volumePartDerivedVars::getPrsmIds(vector<int> &prsmids, size_t start, size_t size) const {
  shadowPart->getPrsmIds(prsmids,start,size) ;
}
 
void volumePartDerivedVars::getHexBlock(vector<Array<int,8> > &hexs, size_t start, size_t size) const {
  shadowPart->getHexBlock(hexs,start,size) ;
}
void volumePartDerivedVars::getHexIds(vector<int> &hexids, size_t start, size_t size) const {
  shadowPart->getHexIds(hexids,start,size) ;
}
void volumePartDerivedVars::getGenCell(vector<int> &genCellNfaces, 
                                       vector<int> &genCellNsides,
                                       vector<int> &genCellNodes) const {
  shadowPart->getGenCell(genCellNfaces,genCellNsides,genCellNodes) ;
}
void volumePartDerivedVars::getGenIds(vector<int> &genids) const {
  shadowPart->getGenIds(genids) ;
}
  
void volumePartDerivedVars::getNodalScalar(string varname, vector<float> &vals) const {
  map<string,derivedVar_t>::const_iterator mi=derivedVars.find(varname) ;
  if(mi==derivedVars.end())
    shadowPart->getNodalScalar(varname,vals) ;
  else {
    derivedVar_t vartype = mi->second ;
    switch(vartype) {
    case VAR_M: 
      {
	vector<float> a ;
	vector<vector3d<float> > v ;
	shadowPart->getNodalScalar("a",a) ;
	shadowPart->getNodalVector("v",v) ;
	vector<float> m(a.size()) ;
	for(size_t i=0;i<a.size();++i)
	  m[i] = norm(v[i])/a[i] ;
	vals.swap(m) ;
      }
      break ;
    case VAR_P:
    case VAR_logp:
      {
	vector<float> pg ;
	shadowPart->getNodalScalar("pg",pg) ;
	vector<float> P(pg.size()) ;
	for(size_t i=0;i<P.size();++i) 
	  P[i] = (vartype==VAR_logp)?log10(max(pg[i]+Pambient,1e-30f)):
	    (pg[i]+Pambient) ;
	vals.swap(P) ;
      }
      break ;
    case VAR_U:
    case VAR_0:
    case VAR_1:
    case VAR_2:
      {
	vector<vector3d<float> > v ;
	shadowPart->getNodalVector("v",v) ;
	vector<float> tmp(v.size()) ;
	for(size_t i=0;i<v.size();++i) {
	  switch(vartype) {
	  case VAR_U:
	    tmp[i] = norm(v[i]) ;
	    break ;
	  case VAR_0:
	    tmp[i] = v[i].x ;
	    break ;
	  case VAR_1:
	    tmp[i] = v[i].y ;
	    break ;
	  case VAR_2:
	    tmp[i] = v[i].z ;
	    break ;
	  default:
	    tmp[i] = 0 ;
	  }
	}
	vals.swap(tmp) ;
      }
      break ;
    case VAR_X:
    case VAR_Y:
    case VAR_Z:
      {
	vector<vector3d<float> > pos ;
	shadowPart->getPos(pos) ;
	vector<float> tmp(pos.size()) ;
	for(size_t i=0;i<pos.size();++i) {
	  switch(vartype) {
	  case VAR_X:
	    tmp[i] = pos[i].x ;
	    break ;
	  case VAR_Y:
	    tmp[i] = pos[i].y ;
	    break ;
	  case VAR_Z:
	    tmp[i] =pos[i].z ;
	    break ;
	  default:
	    tmp[i] = 0 ;
	  }
	}
	vals.swap(tmp) ;
      }
      break ;
    }
  }
}

void volumePartDerivedVars::getNodalVector(string varname, vector<vector3d<float> > &vals) const {
  shadowPart->getNodalVector(varname,vals) ;
}

void volumePartDerivedVars::getNodalIblank(vector<unsigned char> &blank) const {
  shadowPart->getNodalIblank(blank) ;
}

//----
surfacePart::surfacePart(string name, string dir, string iteration,
			 vector<string> vars) {
  partName = name ;
  directory = dir ;
  nnodes = 0 ;
  nquads = 0 ;
  ntrias = 0 ;
  ngenf = 0 ;
  error = true  ;
  string topo_link = dir + "/topo_file."+iteration ;
  ifstream topo_links(topo_link.c_str(),ios::in) ;
  //  if(topo_links.fail()) cerr << "topo_links fail, " << topo_link << endl ;
  if(topo_links.fail()) return ;

  string topolink ;
  topo_links>> topolink ;
  topo_links.close();
  
  posFile = dir + "/pos." + iteration ;
  hid_t file_id = Loci::hdf5OpenFile(posFile.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;
  if(file_id < 0) cerr << posFile << " fail" << endl ;
  if(file_id < 0) return ;
  
  nnodes = sizeElementType(file_id,"data") ;
  Loci::hdf5CloseFile(file_id) ;
  
  topoFile = dir + "/" + topolink ;
  file_id = Loci::hdf5OpenFile(topoFile.c_str(),
                               H5F_ACC_RDONLY,
                               H5P_DEFAULT) ;
  if(file_id < 0) cerr << topoFile << " fail" << endl ;
  if(file_id < 0) return ;

  ngenf = sizeElementType(file_id,"nside_sizes") ;
  nquads = sizeElementType(file_id,"quads") ;
  ntrias = sizeElementType(file_id,"triangles") ;
  
  if(nquads > 0)
    quadSet = interval(0,nquads-1) ;
  if(ntrias > 0)
    triSet = interval(0,ntrias-1) ;
  if(ngenf > 0)
    genSet = interval(0,ngenf-1) ;
  
  Loci::hdf5CloseFile(file_id) ;
  bool has_element_data = false ;
  for(size_t i=0;i<vars.size();++i) {
    string varname = vars[i] ;
    // try scalar
    string svar = dir+"/" + varname+"_sca."+iteration ;
    file_id = Loci::hdf5OpenFile(svar.c_str(),
				 H5F_ACC_RDONLY,
				 H5P_DEFAULT) ;
    bool found_var = false ;
    if(file_id >= 0) {
      int nsz = sizeElementType(file_id,"data") ;
      if(nsz == nnodes) {
	nodalScalarVars[varname] = svar ;
	found_var = true ;
      }
      Loci::hdf5CloseFile(file_id) ;
    }

    if(!found_var) {
      svar = dir+"/" + varname+"_vec."+iteration ;
      file_id = Loci::hdf5OpenFile(svar.c_str(),
				   H5F_ACC_RDONLY,
				   H5P_DEFAULT) ;
      if(file_id >= 0) {
	int nsz = sizeElementType(file_id,"data") ;
	if(nsz == nnodes) {
	  nodalVectorVars[varname] = svar ;
	  found_var = true ;
	}
	Loci::hdf5CloseFile(file_id) ;
      }
    }
    if(!found_var) {
      svar = dir+"/" + varname+"_bsca."+iteration ;
      file_id = Loci::hdf5OpenFile(svar.c_str(),
				   H5F_ACC_RDONLY,
				   H5P_DEFAULT) ;
      if(file_id >= 0) {
	int nsz = sizeElementType(file_id,"data") ;
	if(nsz == (nquads+ntrias+ngenf)) {
	  elementScalarVars[varname] = svar ;
	  found_var = true ;
          has_element_data = true ;
	}
	Loci::hdf5CloseFile(file_id) ;
      }
    }
    if(!found_var) {
      svar = dir+"/" + varname+"_bvec."+iteration ;
      file_id = Loci::hdf5OpenFile(svar.c_str(),
				   H5F_ACC_RDONLY,
				   H5P_DEFAULT) ;
      if(file_id >= 0) {
	int nsz = sizeElementType(file_id,"data") ;
	if(nsz == (nquads+ntrias+ngenf)) {
	  elementVectorVars[varname] = svar ;
	  found_var = true ;
          has_element_data = true ;
	}
	Loci::hdf5CloseFile(file_id) ;
      }
    }
  }
  vector<unsigned char> iblank ;
  string iblank_file = dir +"/iblank."+iteration ;
  file_id = Loci::hdf5OpenFile(iblank_file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if(file_id >=0) {
    vector<unsigned char> tmp(nquads+ntrias+ngenf) ;
    iblank.swap(tmp) ;
    readElementType(file_id,"data",iblank) ;
    Loci::hdf5CloseFile(file_id) ;
    for(size_t i=0;i<iblank.size();++i)
      has_element_data = (iblank[i] > 1) || has_element_data ;
  }
  if(has_element_data) {
    file_id = Loci::hdf5OpenFile(topoFile.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT) ;
    if(nquads > 0) {
      vector<int> tmp(nquads) ;
      readElementType(file_id,"quads_ord",tmp) ;
      quad_ord.swap(tmp) ;
    }

    if(ntrias > 0) {
      vector<int> tmp(ntrias) ;
      readElementType(file_id,"triangles_ord",tmp) ;
      tri_ord.swap(tmp) ;
    }

    if(ngenf > 0) {
      vector<int> tmp(ngenf) ;
      readElementType(file_id,"nside_ord",tmp) ;
      gen_ord.swap(tmp) ;
    }
    Loci::hdf5CloseFile(file_id) ;
    if(iblank.size() > 0) {
      // compute iblanked set
      quadSet = EMPTY ;
      for(int i=0;i<nquads;++i)
        if(iblank[quad_ord[i]] < 2)
          quadSet += i ;
      nquads = quadSet.size() ;
      triSet = EMPTY ;
      for(int i=0;i<ntrias;++i)
        if(iblank[tri_ord[i]] < 2)
          triSet += i ;
      ntrias = triSet.size() ;
      genSet = EMPTY ;
      for(int i=0;i<ngenf;++i)
        if(iblank[gen_ord[i]] < 2)
          genSet += i ;
      ngenf = genSet.size() ;
    }
  }
    
  error = false ;
}

bool surfacePart::hasNodalScalarVar(string var) const {
  map<string,string>::const_iterator mi=nodalScalarVars.find(var) ;
  return (mi != nodalScalarVars.end()) ;
}
bool surfacePart::hasNodalVectorVar(string var) const {
  map<string,string>::const_iterator mi=nodalVectorVars.find(var) ;
  return (mi != nodalVectorVars.end()) ;
}
bool surfacePart::hasElementScalarVar(string var) const {
  map<string,string>::const_iterator mi=elementScalarVars.find(var) ;
  return (mi != elementScalarVars.end()) ;
}
bool surfacePart::hasElementVectorVar(string var) const {
  map<string,string>::const_iterator mi=elementVectorVars.find(var) ;
  return (mi != elementVectorVars.end()) ;
}

vector<string> surfacePart::getNodalScalarVars() const {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=nodalScalarVars.begin();mi!=nodalScalarVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}

vector<string> surfacePart::getNodalVectorVars() const {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=nodalVectorVars.begin();mi!=nodalVectorVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}
  
vector<string> surfacePart::getElementScalarVars() const {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=elementScalarVars.begin();mi!=elementScalarVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}

vector<string> surfacePart::getElementVectorVars() const {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=elementVectorVars.begin();mi!=elementVectorVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}
  
void surfacePart::getQuads(vector<Array<int,4> > &quads) const {
  quads.clear() ;
  if(nquads > 0) {
    hid_t file_id = Loci::hdf5OpenFile(topoFile.c_str(),
				       H5F_ACC_RDONLY,
				       H5P_DEFAULT) ;
    if(file_id < 0) return ;
    int nq = sizeElementType(file_id,"quads") ;

    vector<Array<int,4> > tmp(nq) ;
    readElementType(file_id,"quads",tmp) ;
    quads.resize(nquads) ;
    int cnt = 0 ;
    FORALL(quadSet,ii) {
      quads[cnt] = tmp[ii] ;
      cnt++ ;
    } ENDFORALL ;
    Loci::hdf5CloseFile(file_id) ;
  }
}
void surfacePart::getQuadsIds(vector<int> &quads_ids) const {
  cout << "start getQuadsIds " << endl;
  quads_ids.clear();
  FORALL(quadSet,ii) {
    // quads_ids.push_back(quad_ord[ii]) ;
    quads_ids.push_back(ii);
  } ENDFORALL ;
  cout << "start getQuadsIds " << endl;
}

 
  
void surfacePart::getTrias(vector<Array<int,3> > &trias) const {
  trias.clear() ;
  if(ntrias > 0) {
    hid_t file_id = Loci::hdf5OpenFile(topoFile.c_str(),
				       H5F_ACC_RDONLY,
				       H5P_DEFAULT) ;
    if(file_id < 0) return ;
    int nt = sizeElementType(file_id,"triangles") ;
    vector<Array<int,3> > tmp(nt) ;
    readElementType(file_id,"triangles",tmp) ;
    trias.resize(ntrias) ;
    int cnt = 0 ;
    FORALL(triSet,ii) {
      trias[cnt] = tmp[ii] ;
      cnt++ ;
    } ENDFORALL ;
    Loci::hdf5CloseFile(file_id) ;
  }
}
void  surfacePart::getTriasIds(vector<int> &trias_ids) const{
  cout << "start getTriasIds " << endl;
  trias_ids.clear();
  FORALL(triSet,ii) {
    //  trias_ids.push_back(tri_ord[ii]) ;
    trias_ids.push_back(ii);
  } ENDFORALL ;
  cout << "end getTriasIds " << endl;
}

void surfacePart::getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const {
  numGenFnodes.clear() ;
  genNodes.clear() ;
  if(ngenf > 0) {
    hid_t file_id = Loci::hdf5OpenFile(topoFile.c_str(),
				       H5F_ACC_RDONLY,
				       H5P_DEFAULT) ;
    if(file_id < 0) return ;

    int ng = sizeElementType(file_id,"nside_sizes") ;
    vector<int> tmp(ng) ;
    readElementType(file_id,"nside_sizes",tmp) ;

    int nsz = sizeElementType(file_id,"nside_nodes") ;
    vector<int> tmp2(nsz) ;
    readElementType(file_id,"nside_nodes",tmp2) ;

    vector<int> sum(ng) ;
    sum[0] = 0 ;
    for(int i=1;i<ng;++i)
      sum[i] = sum[i-1]+tmp[i-1] ;
    FORALL(genSet,ii) {
      numGenFnodes.push_back(tmp[ii]) ;
      for(int i=0;i<tmp[ii];++i)
        genNodes.push_back(tmp2[sum[ii]+i]) ;
    } ENDFORALL ;
    Loci::hdf5CloseFile(file_id) ;
  }
}

void  surfacePart::getGenfIds(vector<int> &genface_ids) const{
  cout << "start getGenIds " << endl;
  genface_ids.clear();
  FORALL(genSet,ii) {
    //genface_ids.push_back(gen_ord[ii]) ;
    genface_ids.push_back(ii);
  } ENDFORALL ;
  cout << "end getGenIds " << endl;
}

void surfacePart::getPos(vector<vector3d<float> > &pos) const {
  vector<vector3d<float> > tmp(nnodes) ;
  pos.swap(tmp) ;
  hid_t file_id = Loci::hdf5OpenFile(posFile.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;

  if(file_id < 0) return ;
  readElementType(file_id,"data",pos) ;
  Loci::hdf5CloseFile(file_id) ;
}

void surfacePart::getPos(vector<vector3d<double> > &pos) const {
  vector<vector3d<double> > tmp(nnodes) ;
  pos.swap(tmp) ;
  hid_t file_id = Loci::hdf5OpenFile(posFile.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;

  if(file_id < 0) return ;
  readElementType(file_id,"data",pos) ;
  Loci::hdf5CloseFile(file_id) ;
}

void surfacePart::getNodalScalar(string varname,
				 vector<float> &vals) const {
  vector<float> tmp(nnodes) ;
  vals.swap(tmp) ;
  map<string,string>::const_iterator mi = nodalScalarVars.find(varname) ;
  if(mi == nodalScalarVars.end())
    return ;
  string filename = mi->second ;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;

  if(file_id < 0) return ;
  readElementType(file_id,"data",vals) ;
  Loci::hdf5CloseFile(file_id) ;
  
}
void surfacePart::getNodalVector(string varname,
				 vector<vector3d<float> > &vals) const {
  vector<vector3d<float> > tmp(nnodes) ;
  vals.swap(tmp) ;
  map<string,string>::const_iterator mi = nodalVectorVars.find(varname) ;
  if(mi == nodalScalarVars.end())
    return ;
  string filename = mi->second ;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;

  if(file_id < 0) return ;
  readElementType(file_id,"data",vals) ;
  Loci::hdf5CloseFile(file_id) ;
}

void surfacePart::getElementScalar(string varname,
                                   vector<float> &qvals,
                                   vector<float> &tvals,
                                   vector<float> &gvals) const {
  { vector<float> tmpq(nquads) ;  qvals.swap(tmpq) ; }
  { vector<float> tmpt(ntrias) ;  tvals.swap(tmpt) ; }
  { vector<float> tmpg(ngenf) ;  gvals.swap(tmpg) ; }
  
  map<string,string>::const_iterator mi = elementScalarVars.find(varname) ;
  if(mi == elementScalarVars.end())
    return ;
  string filename = mi->second ;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;

  if(file_id < 0) return ;
  vector<float> vals(quad_ord.size()+tri_ord.size()+gen_ord.size()) ;
  readElementType(file_id,"data",vals) ;
  Loci::hdf5CloseFile(file_id) ;
  int i=0 ;
  FORALL(quadSet,ii) {
    qvals[i] = vals[quad_ord[ii]] ;
    i++ ;
  } ENDFORALL ;
  i=0; 
  FORALL(triSet,ii) {
    tvals[i] = vals[tri_ord[ii]] ;
    i++ ;
  } ENDFORALL ;
  i=0;
  FORALL(genSet,ii) {
    gvals[i] = vals[gen_ord[ii]] ;
    i++ ;
  } ENDFORALL ;
}

void surfacePart::getElementVector(string varname,
                                   vector<vector3d<float> > &qvals,
                                   vector<vector3d<float> > &tvals,
                                   vector<vector3d<float> > &gvals) const {
  { vector<vector3d<float> > tmpq(nquads) ;  qvals.swap(tmpq) ; }
  { vector<vector3d<float> > tmpt(ntrias) ;  tvals.swap(tmpt) ; }
  { vector<vector3d<float> > tmpg(ngenf) ;  gvals.swap(tmpg) ; }
  
  map<string,string>::const_iterator mi = elementVectorVars.find(varname) ;
  if(mi == elementVectorVars.end())
    return ;
  string filename = mi->second ;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;

  if(file_id < 0) return ;
  vector<vector3d<float> > vals(quad_ord.size()+tri_ord.size()+gen_ord.size()) ;
  readElementType(file_id,"data",vals) ;
  Loci::hdf5CloseFile(file_id) ;
  int i=0 ;
  FORALL(quadSet,ii) {
    qvals[i] = vals[quad_ord[ii]] ;
    i++ ;
  } ENDFORALL ;
  i=0; 
  FORALL(triSet,ii) {
    tvals[i] = vals[tri_ord[ii]] ;
    i++ ;
  } ENDFORALL ;
  i=0;
  FORALL(genSet,ii) {
    gvals[i] = vals[gen_ord[ii]] ;
    i++ ;
  } ENDFORALL ;
}
//----
void surfacePartDerivedVars::processDerivedVars(const vector<string> &vars) {
  for(size_t i=0;i<vars.size();++i) {
    if(vars[i] == "m" && !shadowPart->hasNodalScalarVar("m")) {
      if(shadowPart->hasNodalScalarVar("a") && 
	 shadowPart->hasNodalVectorVar("v"))
	derivedVars["m"] = VAR_M ;
    }
    if(vars[i] == "P" && !shadowPart->hasNodalScalarVar("P")) {
      if(shadowPart->hasNodalScalarVar("pg")) {
	derivedVars["P"] = VAR_P ;
      }
    }
    if(vars[i] == "p" && !shadowPart->hasNodalScalarVar("p")) {
      if(shadowPart->hasNodalScalarVar("pg")) {
	derivedVars["p"] = VAR_logp ;
      }
    }
    if(vars[i] == "u" && !shadowPart->hasNodalScalarVar("u")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["u"] = VAR_U ;
      }
    }
    if(vars[i] == "0" && !shadowPart->hasNodalScalarVar("0")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["0"] = VAR_0 ;
      }
    }
    if(vars[i] == "1" && !shadowPart->hasNodalScalarVar("1")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["1"] = VAR_1 ;
      }
    }
    if(vars[i] == "2" && !shadowPart->hasNodalScalarVar("2")) {
      if(shadowPart->hasNodalVectorVar("v")) {
	derivedVars["2"] = VAR_2 ;
      }
    }
    if(vars[i] == "x")
      derivedVars["x"] = VAR_X ;
    if(vars[i] == "y")
      derivedVars["y"] = VAR_Y ;
    if(vars[i] == "z")
      derivedVars["z"] = VAR_Z ;
  }
}

surfacePartDerivedVars::surfacePartDerivedVars(surfacePartP part, 
					       string output_dir,
					       string casename ,
					       string iteration, 
					       vector<string> vars) {
  error = part->fail() ;
  partName = part->getPartName() ;
  nnodes = part->getNumNodes() ;
  nquads = part->getNumQuads() ;
  ntrias = part->getNumTrias() ;
  ngenf = part->getNumGenfc() ;
  shadowPart = part ;

  string filename = output_dir+"/Pambient_par." + iteration +"_" + casename ;
  struct stat tmpstat ;
  if(stat(filename.c_str(),&tmpstat) == 0) {
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;
    Pambient = 0 ;
    if(file_id >= 0) {
      fact_db facts ;
      param<float> Pamb ;
      readData(file_id,"Pambient",Pamb.Rep(),EMPTY,facts) ;
      Loci::hdf5CloseFile(file_id) ;
      Pambient = *Pamb ;
    } else { 
      cerr << "unable to open " << filename << endl ;
    }
    processDerivedVars(vars) ;
  }
}

bool surfacePartDerivedVars::hasNodalScalarVar(string var) const {
  map<string,derivedVar_t>::const_iterator mi=derivedVars.find(var) ;
  if(mi==derivedVars.end())
    return shadowPart->hasNodalScalarVar(var) ;
  else
    return true ;
}
bool surfacePartDerivedVars::hasNodalVectorVar(string var) const {
  return shadowPart->hasNodalVectorVar(var) ;
}
bool surfacePartDerivedVars::hasElementScalarVar(string var) const {
  return shadowPart->hasElementScalarVar(var) ;
}
bool surfacePartDerivedVars::hasElementVectorVar(string var) const {
  return shadowPart->hasElementVectorVar(var)  ;
}

vector<string> surfacePartDerivedVars::getNodalScalarVars() const {
  vector<string> tmp = shadowPart->getNodalScalarVars() ;
  map<string,derivedVar_t>::const_iterator mi ;
  for(mi=derivedVars.begin();mi!=derivedVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}

vector<string> surfacePartDerivedVars::getNodalVectorVars() const {
  return shadowPart->getNodalVectorVars() ;
}
  
vector<string> surfacePartDerivedVars::getElementScalarVars() const {
  return shadowPart->getElementScalarVars() ;
}

vector<string> surfacePartDerivedVars::getElementVectorVars() const {
  return shadowPart->getElementVectorVars() ;
}
  
void surfacePartDerivedVars::getQuads(vector<Array<int,4> > &quads) const {
  shadowPart->getQuads(quads) ;
}

void surfacePartDerivedVars::getTrias(vector<Array<int,3> > &trias) const {
  shadowPart->getTrias(trias) ;
}

void surfacePartDerivedVars::getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const {
  shadowPart->getGenf(numGenFnodes,genNodes) ;
}
void  surfacePartDerivedVars::getQuadsIds(vector<int> &quads_ids) const{
  shadowPart->getQuadsIds(quads_ids);
}
void  surfacePartDerivedVars::getTriasIds(vector<int> &trias_ids) const{
  shadowPart->getTriasIds(trias_ids); 
}
void surfacePartDerivedVars::getGenfIds(vector<int> &genface_ids) const{
  shadowPart->getGenfIds(genface_ids);
}
void surfacePartDerivedVars::getPos(vector<vector3d<float> > &pos) const {
  shadowPart->getPos(pos) ;
}
void surfacePartDerivedVars::getPos(vector<vector3d<double> > &pos) const {
  shadowPart->getPos(pos) ;
}
void surfacePartDerivedVars::getNodalScalar(string varname,
					    vector<float> &vals) const {
  map<string,derivedVar_t>::const_iterator mi=derivedVars.find(varname) ;
  if(mi==derivedVars.end())
    shadowPart->getNodalScalar(varname,vals) ;
  else {
    derivedVar_t vartype = mi->second ;
    switch(vartype) {
    case VAR_M: 
      {
	vector<float> a ;
	vector<vector3d<float> > v ;
	shadowPart->getNodalScalar("a",a) ;
	shadowPart->getNodalVector("v",v) ;
	vector<float> m(a.size()) ;
	for(size_t i=0;i<a.size();++i)
	  m[i] = norm(v[i])/a[i] ;
	vals.swap(m) ;
      }
      break ;
    case VAR_P:
    case VAR_logp:
      {
	vector<float> pg ;
	shadowPart->getNodalScalar("pg",pg) ;
	vector<float> P(pg.size()) ;
	for(size_t i=0;i<P.size();++i) 
	  P[i] = (vartype==VAR_logp)?log10(max(pg[i]+Pambient,1e-30f)):
	    (pg[i]+Pambient) ;
	vals.swap(P) ;
      }
      break ;
    case VAR_U:
    case VAR_0:
    case VAR_1:
    case VAR_2:
      {
	vector<vector3d<float> > v ;
	shadowPart->getNodalVector("v",v) ;
	vector<float> tmp(v.size()) ;
	for(size_t i=0;i<v.size();++i) {
	  switch(vartype) {
	  case VAR_U:
	    tmp[i] = norm(v[i]) ;
	    break ;
	  case VAR_0:
	    tmp[i] = v[i].x ;
	    break ;
	  case VAR_1:
	    tmp[i] = v[i].y ;
	    break ;
	  case VAR_2:
	    tmp[i] = v[i].z ;
	    break ;
	  default:
	    tmp[i] = 0 ;
	  }
	}
	vals.swap(tmp) ;
      }
      break ;
    case VAR_X:
    case VAR_Y:
    case VAR_Z:
      {
	vector<vector3d<float> > pos ;
	shadowPart->getPos(pos) ;
	vector<float> tmp(pos.size()) ;
	for(size_t i=0;i<pos.size();++i) {
	  switch(vartype) {
	  case VAR_X:
	    tmp[i] = pos[i].x ;
	    break ;
	  case VAR_Y:
	    tmp[i] = pos[i].y ;
	    break ;
	  case VAR_Z:
	    tmp[i] =pos[i].z ;
	    break ;
	  default:
	    tmp[i] = 0 ;
	  }
	}
	vals.swap(tmp) ;
      }
      break ;
    }
  }
}
void surfacePartDerivedVars::getNodalVector(string varname,
                                            vector<vector3d<float> > &vals) const {
  shadowPart->getNodalVector(varname,vals) ;
}

void surfacePartDerivedVars::getElementScalar(string varname,
                                              vector<float> &qvals,
                                              vector<float> &tvals,
                                              vector<float> &gvals) const {
  shadowPart->getElementScalar(varname,qvals,tvals,gvals) ;
}

void surfacePartDerivedVars::getElementVector(string varname,
                                              vector<vector3d<float> > &qvals,
                                              vector<vector3d<float> > &tvals,
                                              vector<vector3d<float> > &gvals) const {
  shadowPart->getElementVector(varname,qvals,tvals,gvals) ;
}

//----
surfacePartCopy::surfacePartCopy(string name,
                                 vector<Array<int,3> > &triangles,
                                 vector<int> &tria_ids,
                                 vector<Array<int,4> > &quads,
                                 vector<int> &quad_ids,
                                 vector<int> &genface2n,
                                 vector<int> &gnodes,
                                 vector<int> &gen_ids) {
  error = false ;

  vector<int> node_set ;
  for(size_t i=0;i<gnodes.size();++i)
    node_set.push_back(gnodes[i]) ;
        
  for(size_t i=0;i<triangles.size();++i) {
    node_set.push_back(triangles[i][0]) ;
    node_set.push_back(triangles[i][1]) ;
    node_set.push_back(triangles[i][2]) ;
  }
  for(size_t i=0;i<quads.size();++i) {
    node_set.push_back(quads[i][0]) ;
    node_set.push_back(quads[i][1]) ;
    node_set.push_back(quads[i][2]) ;
    node_set.push_back(quads[i][3]) ;
  }
  sort(node_set.begin(),node_set.end()) ;
  node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;
  
  map<int,int> nmap ;
  for(size_t i=0;i<node_set.size();++i) {
    nmap[node_set[i]] = i+1 ;
  }
  nodemap = node_set ;
  trifaces = triangles ;
  quadfaces = quads ;
  nfacenodes = genface2n ;
  gennodes = gnodes ;
  triaIds = tria_ids;
  quadIds = quad_ids;
  genIds = gen_ids;
  
  for(size_t i=0;i<trifaces.size();++i) {
    trifaces[i][0] = nmap[trifaces[i][0]] ;
    trifaces[i][1] = nmap[trifaces[i][1]] ;
    trifaces[i][2] = nmap[trifaces[i][2]] ;
  }
  for(size_t i=0;i<quadfaces.size();++i) {
    quadfaces[i][0] = nmap[quadfaces[i][0]] ;
    quadfaces[i][1] = nmap[quadfaces[i][1]] ;
    quadfaces[i][2] = nmap[quadfaces[i][2]] ;
    quadfaces[i][3] = nmap[quadfaces[i][3]] ;
  }
  for(size_t i=0;i<gennodes.size();++i)
    gennodes[i] = nmap[gennodes[i]] ;
  
  partName = name ;
  nnodes = nodemap.size() ;
  nquads = quadfaces.size() ;
  ntrias = trifaces.size() ;
  ngenf = nfacenodes.size() ;
}

void surfacePartCopy::registerPos(const vector<vector3d<double> > &posvol) {
  vector<vector3d<double> > tmp(nodemap.size()) ;
  pos.swap(tmp) ;
  for(size_t i=0;i<nodemap.size();++i)
    pos[i] = posvol[nodemap[i]-1] ;
}

void surfacePartCopy::registerNodalScalar(string name, const vector<float> &val) {
  vector<float> tmp(nodemap.size()) ;
  for(size_t i=0;i<nodemap.size();++i)
    tmp[i] = val[nodemap[i]-1] ;
  nodalScalars[name] = tmp ;
}

void surfacePartCopy::registerNodalVector(string name, const vector<vector3d<float> > &val) {
  vector<vector3d<float> > tmp(nodemap.size())
    ;
  for(size_t i=0;i<nodemap.size();++i)
    tmp[i] = val[nodemap[i]-1] ;
  nodalVectors[name] = tmp ;
}

void surfacePartCopy::
registerElementScalar(string name, 
		      const vector<float> &qval,
		      const vector<float> &tval,
		      const vector<float> &gval)  {
  Array<vector<float>,3> &tmp = elementScalars[name] ;
  tmp[0] = qval ;
  tmp[1] = tval ;
  tmp[2] = gval ;
}

void surfacePartCopy::
registerElementVector(string name, 
		      const vector<vector3d<float> > &qval,
		      const vector<vector3d<float> > &tval,
		      const vector<vector3d<float> > &gval)  {
  Array<vector<vector3d<float> >,3> &tmp = elementVectors[name] ;
  tmp[0] = qval ;
  tmp[1] = tval ;
  tmp[2] = gval ;
}


bool surfacePartCopy::hasNodalScalarVar(string var) const {
  map<string,vector<float> >::const_iterator mi ;
  mi = nodalScalars.find(var) ;
  return !(mi == nodalScalars.end()) ;
}
bool surfacePartCopy::hasNodalVectorVar(string var) const {
  map<string,vector<vector3d<float> > >::const_iterator mi ;
  mi = nodalVectors.find(var) ;
  return !(mi == nodalVectors.end()) ;
}
bool surfacePartCopy::hasElementScalarVar(string var) const {
  map<string,Array<vector<float>,3> >::const_iterator mi ;
  mi = elementScalars.find(var) ;
  return mi != elementScalars.end() ;
}
bool surfacePartCopy::hasElementVectorVar(string var) const {
  map<string,Array<vector<vector3d<float> >,3> >::const_iterator mi ;
  mi = elementVectors.find(var) ;
  return mi != elementVectors.end() ;
}

vector<string> surfacePartCopy::getNodalScalarVars() const {
  vector<string> tmp ;
  map<string,vector<float> >::const_iterator mi ;
  for(mi=nodalScalars.begin();mi!=nodalScalars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}

vector<string> surfacePartCopy::getNodalVectorVars() const {
  vector<string> tmp ;
  map<string,vector<vector3d<float> > >::const_iterator mi ;
  for(mi=nodalVectors.begin();mi!=nodalVectors.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}
  
vector<string> surfacePartCopy::getElementScalarVars() const {
  vector<string> tmp ;
  map<string,Array<vector<float>,3> >::const_iterator mi ;
  for(mi=elementScalars.begin();mi!=elementScalars.end();++mi) 
    tmp.push_back(mi->first) ;
  return tmp ;
}

vector<string> surfacePartCopy::getElementVectorVars() const {
  vector<string> tmp ;
  map<string,Array<vector<vector3d<float> >,3> >::const_iterator mi ;
  for(mi=elementVectors.begin();mi!=elementVectors.end();++mi) 
    tmp.push_back(mi->first) ;
  return tmp ;
}
  
void surfacePartCopy::getQuads(vector<Array<int,4> > &quads) const {
  quads = quadfaces ;
}
void surfacePartCopy::getQuadsIds(vector<int> &quads_ids) const{
  quads_ids = quadIds;
}
void surfacePartCopy::getTrias(vector<Array<int,3> > &trias) const {
  trias = trifaces ;
}
void surfacePartCopy::getTriasIds(vector<int> &trias_ids) const{
  trias_ids = triaIds;
}
void surfacePartCopy::getGenf(vector<int> &numGenFnodes, vector<int> &genNodes) const {
  numGenFnodes = nfacenodes ;
  genNodes = gennodes ;
}
void surfacePartCopy::getGenfIds(vector<int> &genface_ids) const{
  genface_ids = genIds;
}

void surfacePartCopy::getPos(vector<vector3d<float> > &pos_out) const {
  pos_out.resize(pos.size());
  for(size_t i = 0; i < pos.size(); i++){
    pos_out[i].x = pos[i].x;
    pos_out[i].y = pos[i].y;
    pos_out[i].z = pos[i].z;
  }
  
}

void surfacePartCopy::getPos(vector<vector3d<double> > &pos_out) const {
  pos_out = pos;
}

void surfacePartCopy::getNodalScalar(string varname,
                                     vector<float> &vals) const {
  map<string,vector<float> >::const_iterator mi ;
  mi = nodalScalars.find(varname) ;
  if(!(mi == nodalScalars.end()))
    vals = mi->second ;
}
void surfacePartCopy::getNodalVector(string varname,
                                     vector<vector3d<float> > &vals) const {
  map<string,vector<vector3d<float> > >::const_iterator mi ;
  mi = nodalVectors.find(varname) ;
  if(!(mi == nodalVectors.end()))
    vals = mi->second ;
}

void surfacePartCopy::getElementScalar(string name,
                                       vector<float> &qvals,
                                       vector<float> &tvals,
                                       vector<float> &gvals) const {
  map<string,Array<vector<float>,3> >::const_iterator mi ;
  mi = elementScalars.find(name) ;
  qvals = mi->second[0] ;
  tvals = mi->second[1] ;
  gvals = mi->second[2] ;
}

void surfacePartCopy::getElementVector(string name,
                                       vector<vector3d<float> > &qvals,
                                       vector<vector3d<float> > &tvals,
                                       vector<vector3d<float> > &gvals) const {
  map<string,Array<vector<vector3d<float> >,3> >::const_iterator mi ;
  mi = elementVectors.find(name) ;
  qvals = mi->second[0] ;
  tvals = mi->second[1] ;
  gvals = mi->second[2] ;
}

particlePart::particlePart(string output_dir, string iteration, string casename,
			   vector<string> vars,
			   int maxparticles) {
  error = true ;
  partName = "Particles" ;
  directory = output_dir ;
  posfile = output_dir + "/particle_pos."+iteration + "_" + casename ;
  
  struct stat tmpstat ;
  if(stat(posfile.c_str(),&tmpstat) != 0) {
    return ;
  }
  hid_t file_id = Loci::hdf5OpenFile(posfile.c_str(),
				     H5F_ACC_RDONLY,
				     H5P_DEFAULT) ;
  if(file_id < 0) 
    return ;
  numParticles = sizeElementType(file_id, "particle position") ;
  H5Fclose(file_id) ;
  stride_size = 1 ;
  if(maxparticles > 0) {
    stride_size = numParticles/maxparticles ;
    stride_size = max(stride_size,1) ;
    int num_blocks = numParticles/stride_size ;
    numParticles = num_blocks ; //*stride_size ;
  }  
  for(size_t i=0;i<vars.size();++i) {
    string varname = vars[i] ;
    string scalarfile = output_dir+"/"+varname+"_ptsca."
      +iteration+"_"+casename ;
    if(stat(scalarfile.c_str(),&tmpstat) == 0) {
      scalarVars[varname] = scalarfile ;
    } else {
      string vectorfile = output_dir+"/"+varname+"_ptvec."
	+iteration+"_"+casename ;
      if(stat(vectorfile.c_str(),&tmpstat) == 0) 
	vectorVars[varname] = vectorfile ;
    }
  }
}

bool particlePart::hasScalarVar(string var) const {
  map<string,string>::const_iterator mi ;
  mi = scalarVars.find(var) ;
  return !(mi == scalarVars.end()) ;
}
bool particlePart::hasVectorVar(string var) const {
  map<string,string>::const_iterator mi ;
  mi = vectorVars.find(var) ;
  return !(mi == vectorVars.end()) ;
}
std::vector<string> particlePart::getScalarVars() const {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=scalarVars.begin();mi!=scalarVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}
std::vector<string> particlePart::getVectorVars() const {
  vector<string> tmp ;
  map<string,string>::const_iterator mi ;
  for(mi=vectorVars.begin();mi!=vectorVars.end();++mi)
    tmp.push_back(mi->first) ;
  return tmp ;
}
void particlePart::getParticlePositions(vector<vector3d<float> > &ppos) const {
  hid_t file_id = Loci::hdf5OpenFile(posfile.c_str(),
				     H5F_ACC_RDONLY, H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to open file '" << posfile << "'!" << endl ;
    return ;
  }
  size_t np = sizeElementType(file_id, "particle position") ;
  vector<vector3d<float> > tmp(np) ;
  readElementType(file_id, "particle position", tmp) ;
  Loci::hdf5CloseFile(file_id) ;
  
  if(stride_size == 1)
    ppos.swap(tmp) ;
  else {
    vector<vector3d<float> > cpy(numParticles) ;
    for(size_t i=0;i<numParticles;++i)
      cpy[i] = tmp[i*stride_size] ;
    ppos.swap(cpy) ;
  }
}
void particlePart::getParticleScalar(string varname, vector<float> &val) const {
  map<string,string>::const_iterator mi ;
  mi = scalarVars.find(varname) ;
  string filename = mi->second ;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				     H5F_ACC_RDONLY, H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to open file '" << filename << "'!" << endl ;
  }
  size_t np = sizeElementType(file_id, varname.c_str()) ;
  vector<float> scalar(np) ;
  readElementType(file_id, varname.c_str(), scalar) ;
  Loci::hdf5CloseFile(file_id) ;
  
  if(stride_size == 1)
    val.swap(scalar) ;
  else {
    vector<float > cpy(numParticles) ;
    for(size_t i=0;i<numParticles;++i)
      cpy[i] = scalar[i*stride_size] ;
    val.swap(cpy) ;
  }
}

void particlePart::getParticleVector(string varname, 
				     vector<vector3d<float> > &val) const {
  map<string,string>::const_iterator mi ;
  mi = vectorVars.find(varname) ;
  string filename = mi->second ;
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
				     H5F_ACC_RDONLY, H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to open file '" << filename << "'!" << endl ;
    return ;
  }
  size_t np = sizeElementType(file_id, varname.c_str()) ;
  vector<vector3d<float> > tmp(np) ;
  readElementType(file_id, varname.c_str(), tmp) ;
  Loci::hdf5CloseFile(file_id) ;
  
  if(stride_size == 1)
    val.swap(tmp) ;
  else {
    vector<vector3d<float> > cpy(numParticles) ;
    for(size_t i=0;i<numParticles;++i)
      cpy[i] = tmp[i*stride_size] ;
    val.swap(cpy) ;
  }
}


void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration) {
  if(var_name == "m") {
    string filename = output_dir+"/a_sca."+iteration + "_" + casename ;
    
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    fact_db facts ;
    store<float> soundSpeed ;
    readData(file_id,"a",soundSpeed.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    filename = output_dir+"/v_vec." + iteration +"_" + casename ;
    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    readData(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    FORALL(dom,nd) {
      float m = norm(u[nd])/soundSpeed[nd] ;
      dval[c++] = m ;
    } ENDFORALL ;
  } else if(var_name == "p" || var_name == "P") {
    string filename = output_dir+"/pg_sca."+iteration + "_" + casename ;
    
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    fact_db facts ;
    store<float> pg ;
    readData(file_id,"pg",pg.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    filename = output_dir+"/Pambient_par." + iteration +"_" + casename ;

    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    param<float> Pambient ;
    readData(file_id,"Pambient",Pambient.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = pg.domain() ;
    bool log = (var_name == "p") ;
    int c = 0 ;
    FORALL(dom,nd) {
      float p = pg[nd]+Pambient[nd] ;
      if(log)
        dval[c++] = log10(p) ;
      else
        dval[c++] = p ;
    } ENDFORALL ;
  } else if (var_name == "u") {
    fact_db facts ;
    string filename = output_dir+"/v_vec." + iteration +"_" + casename ;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    readData(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    FORALL(dom,nd) {
      float m = norm(u[nd]) ;
      dval[c++] = m ;
    } ENDFORALL ;
  } else if(var_name == "x" || var_name =="y" || var_name == "z") {
    store<vector3d<float> > pos ;
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
    entitySet dom = pos.domain() ;
    int c = 0 ;
    if(var_name == "x") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].x ;
      } ENDFORALL ;
    }
    if(var_name == "y") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].y ;
      } ENDFORALL ;
    }
    if(var_name == "z") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].z ;
      } ENDFORALL ;
    }
  } else if(var_name == "0" || var_name =="1" || var_name == "2") {
    fact_db facts ;
    string filename = output_dir+"/v_vec." + iteration +"_" + casename ;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    readData(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    if(var_name == "0") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].x ;
      } ENDFORALL ;
    }
    if(var_name == "1") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].y ;
      } ENDFORALL ;
    }
    if(var_name == "2") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].z ;
      } ENDFORALL ;
    }
  } else {
    cerr << "don't know how to get derived variable " << var_name << endl ;
  }
}




vector<string> volumeSurfaceNames(string output_dir, string iteration,
				  string casename) {
  string gridtopo = getTopoFileName(output_dir, casename, iteration) ;
  hid_t file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
#ifdef H5_USE_16_API
  hid_t bndg = H5Gopen(file_id,"boundaries") ;
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT) ;
#endif
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  vector<string>  bc_names ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;
    bc_names.push_back(string(buf)) ;
  }
  H5Gclose(bndg) ;
  H5Fclose(file_id) ;
  return bc_names ;
}

void extractVolumeSurfaces(vector<surfacePartP> &volSurface,
                           volumePartP vp,
                           string output_dir,
                           string iteration,
                           string casename,
			   vector<string> varlist) {
  string gridtopo = getTopoFileName(output_dir, casename, iteration) ;

  cout << "extracting topology from '" << gridtopo << "'" << endl;

  hid_t file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;

  map<string,int> elementScalars ;
  map<string,int> elementVectors ;
  vector<string> variable_file(varlist.size()) ;
  vector<map<int,int> > entityMap(varlist.size()) ;
  vector<vector<int> > entityIds(varlist.size()) ;
  vector<vector<float> > scalarElementVars(varlist.size()) ;
  vector<vector<vector3d<float> > > vectorElementVars(varlist.size()) ;

  for(size_t i = 0;i<varlist.size();++i) {
    string var = varlist[i] ;
    string filename = output_dir+'/' + var + "_bnd." + iteration + "_" + casename ;
    struct stat tmpstat ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      elementScalars[var] = i ;
      variable_file[i] = filename ;
      continue ;
    }
    filename = output_dir+'/' + var + "_bndvec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      elementVectors[var] = i ;
      variable_file[i] = filename ;
      continue ;
    }
  }

  map<string,int>::const_iterator mi ;
  for(mi=elementScalars.begin();mi!=elementScalars.end();++mi) {
    int id = mi->second ;
    string filename(variable_file[id]) ;
    string varname = varlist[id] ;
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
    size_t nbel = sizeElementType(di,"entityIds") ;
    
    vector<int> elemIds(nbel) ;
    readElementType(di,"entityIds",elemIds) ;
          
    H5Gclose(di) ;
    vector<float> var(nbel) ;
    readElementType(file_id,varname.c_str(),var) ;
    
    entityIds[id] = elemIds ;
    for(size_t i=0;i<nbel;++i) {
      entityMap[id][elemIds[i]] = int(i) ;
    }
    scalarElementVars[id] = var ;
    H5Fclose(file_id) ;
  }
  
  for(mi=elementVectors.begin();mi!=elementVectors.end();++mi) {
    int id = mi->second ;
    string filename(variable_file[id]) ;
    string varname = varlist[id] ;
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
    size_t nbel = sizeElementType(di,"entityIds") ;
    
    vector<int> elemIds(nbel) ;
    readElementType(di,"entityIds",elemIds) ;
          
    H5Gclose(di) ;
    vector<vector3d<float> > var(nbel) ;
    readElementType(file_id,varname.c_str(),var) ;
    
    entityIds[id] = elemIds ;
    for(size_t i=0;i<nbel;++i) {
      entityMap[id][elemIds[i]] = int(i) ;
    }
    vectorElementVars[id] = var ;
    H5Fclose(file_id) ;
  }

#ifdef H5_USE_16_API
  hid_t bndg = H5Gopen(file_id,"boundaries") ;
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT) ;
#endif
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  vector<string>  bc_names ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;
    bc_names.push_back(string(buf)) ;
  }
  vector<surfacePartCopy * > surfaceWork(num_bcs) ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
#ifdef H5_USE_16_API
    hid_t bcg = H5Gopen(bndg,bc_names[bc].c_str()) ;
#else
    hid_t bcg = H5Gopen(bndg,bc_names[bc].c_str(),H5P_DEFAULT) ;
#endif
    
    size_t nquads = sizeElementType(bcg,"quads") ;
    size_t ntrias = sizeElementType(bcg,"triangles") ;
    size_t ngeneral = sizeElementType(bcg,"nside_sizes") ;
        
    vector<Array<int,3> > trias(ntrias) ;
    readElementType(bcg,"triangles",trias) ;
    vector<Array<int,4> > quads(nquads) ;
    readElementType(bcg,"quads",quads) ;
    
    vector<int> nside_sizes(ngeneral) ;
    readElementType(bcg,"nside_sizes",nside_sizes) ;
    size_t nside_nodes_size = sizeElementType(bcg,"nside_nodes") ;
    vector<int> nside_nodes(nside_nodes_size) ;
    readElementType(bcg,"nside_nodes",nside_nodes) ;

    vector<int > trias_id(ntrias) ;
    readElementType(bcg,"triangles_id",trias_id) ;
    vector<int > quads_id(nquads) ;
    readElementType(bcg,"quads_id",quads_id) ;
    vector<int > nside_id(ngeneral) ;
    readElementType(bcg,"nside_id",nside_id) ;

    vector<unsigned char> iblank ;
    vp->getNodalIblank(iblank) ;

    if(iblank.size() > 0) {
      int cnt = 0 ;
      for(size_t i=0;i<ntrias;++i) {
	bool blank = true ;
	for(int j=0;j<3;++j)
	  if(iblank[trias[i][j]-1] < 2)
	    blank = false ;
	if(blank)
	  cnt++ ;
	else 
	  if(cnt != 0) { // If there are some blanked copy into place
	    trias[i-cnt]=trias[i] ;
	    trias_id[i-cnt] = trias_id[i] ;
	  }
      }
      if(cnt > 0) {
	size_t newsz = trias.size()-cnt ;
	trias.resize(newsz) ;
	trias_id.resize(newsz) ;
      }
      cnt = 0 ;
      for(size_t i=0;i<nquads;++i) {
	bool blank = true ;
	for(int j=0;j<4;++j)
	  if(iblank[quads[i][j]-1] < 2)
	    blank = false ;
	if(blank)
	  cnt++ ;
	else 
	  if(cnt != 0) { // If there are some blanked copy into place
	    quads[i-cnt]=quads[i] ;
	    quads_id[i-cnt] = quads_id[i] ;
	  }
      }
      if(cnt > 0) {
	size_t newsz = quads.size()-cnt ;
	quads.resize(newsz) ;
	quads_id.resize(newsz) ;
      }
      cnt  = 0 ;
      int cnt2 = 0 ;
      int nside_off = 0 ;
      for(size_t i=0;i<ngeneral;++i) {
	bool blank = true ;
	for(int j=0;j<nside_sizes[i];++j) {
	  if(iblank[nside_nodes[nside_off+j]-1] < 2)
	    blank = false ;
	}
	if(blank) {
	  cnt++ ;
	  cnt2 += nside_sizes[i] ;
	} else {
	  if(cnt != 0) {
	    nside_sizes[i-cnt] = nside_sizes[i] ;
	    nside_id[i-cnt] = nside_id[i] ;
	    for(int j=0;j<nside_sizes[i];++j)
	      nside_nodes[nside_off-cnt2+j] = nside_nodes[nside_off+j] ;
	  }
	}
	nside_off += nside_sizes[i] ;
      }
      if(cnt > 0) {
	size_t newsz = nside_sizes.size()-cnt ;
	nside_sizes.resize(newsz) ;
	nside_id.resize(newsz) ;
	size_t newsz2 = nside_nodes.size()-cnt2 ;
	nside_nodes.resize(newsz2) ;
      }
    }
    surfaceWork[bc] = 
      new surfacePartCopy(bc_names[bc],trias,trias_id, quads,quads_id, nside_sizes,nside_nodes, nside_id) ;

    for(mi=elementScalars.begin();mi!=elementScalars.end();++mi) {
      int id = mi->second ;
      string varname = varlist[id] ;
      vector<float> qvals(quads.size()) ;
      vector<float> tvals(trias.size()) ;
      vector<float> gvals(nside_sizes.size()) ;
      bool valid = true ;
      for(size_t i=0;i<quads.size();++i) {
	map<int,int>::const_iterator ii = entityMap[id].find(quads_id[i]) ;
	if(ii==entityMap[id].end()) {
	  valid = false ;
	  break ;
	} else {
	  qvals[i] = scalarElementVars[id][ii->second] ;
	}
      }
      for(size_t i=0;i<trias.size();++i) {
	map<int,int>::const_iterator ii = entityMap[id].find(trias_id[i]) ;
	if(ii==entityMap[id].end()) {
	  valid = false ;
	  break ;
	} else {
	  tvals[i] = scalarElementVars[id][ii->second] ;
	}
      }
      for(size_t i=0;i<nside_sizes.size();++i) {
	map<int,int>::const_iterator ii = entityMap[id].find(nside_id[i]) ;
	if(ii==entityMap[id].end()) {
	  valid = false ;
	  break ;
	} else {
	  gvals[i] = scalarElementVars[id][ii->second] ;
	}
      }
      if(valid)
	surfaceWork[bc]->registerElementScalar(varname,qvals,tvals,gvals) ;
    }

    for(mi=elementVectors.begin();mi!=elementVectors.end();++mi) {
      int id = mi->second ;
      string varname = varlist[id] ;
      vector<vector3d<float> > qvals(quads.size()) ;
      vector<vector3d<float> > tvals(trias.size()) ;
      vector<vector3d<float> > gvals(nside_sizes.size()) ;
      bool valid = true ;
      for(size_t i=0;i<quads.size();++i) {
	map<int,int>::const_iterator ii = entityMap[id].find(quads_id[i]) ;
	if(ii==entityMap[id].end()) {
	  valid = false ;
	} else {
	  qvals[i] = vectorElementVars[id][ii->second] ;
	}
      }
      for(size_t i=0;i<trias.size();++i) {
	map<int,int>::const_iterator ii = entityMap[id].find(trias_id[i]) ;
	if(ii==entityMap[id].end()) {
	  valid = false ;
	} else {
	  tvals[i] = vectorElementVars[id][ii->second] ;
	}
      }
      for(size_t i=0;i<nside_sizes.size();++i) {
	map<int,int>::const_iterator ii = entityMap[id].find(nside_id[i]) ;
	if(ii==entityMap[id].end()) {
	  valid = false ;
	} else {
	  gvals[i] = vectorElementVars[id][ii->second] ;
	}
      }
      
      if(valid)
	surfaceWork[bc]->registerElementVector(varname,qvals,tvals,gvals) ;
    }
    
  }
  H5Gclose(bndg) ;
  H5Fclose(file_id) ;
  
  {
    vector<vector3d<double> > posvol ;
    vp->getPos(posvol) ;
    for(hsize_t bc=0;bc<num_bcs;++bc) {
      surfaceWork[bc]->registerPos(posvol) ;
    }
  }
  {
    vector<string> vars = vp->getNodalScalarVars() ;
    for(size_t i=0;i<vars.size();++i) {
      vector<float> val ;
      vp->getNodalScalar(vars[i],val) ;
      for(hsize_t bc=0;bc<num_bcs;++bc) {
        surfaceWork[bc]->registerNodalScalar(vars[i],val) ;
      }
    }
  }
  {
    vector<string> vars = vp->getNodalVectorVars() ;
    for(size_t i=0;i<vars.size();++i) {
      vector<vector3d<float> > val ;
      vp->getNodalVector(vars[i],val) ;
      for(hsize_t bc=0;bc<num_bcs;++bc) {
        surfaceWork[bc]->registerNodalVector(vars[i],val) ;
      }
    }
  }
  volSurface = vector<surfacePartP>(num_bcs) ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    volSurface[bc] = surfaceWork[bc] ;
  }
  
}





namespace Loci {
  void disableDebugDir() ;
}
int main(int ac, char *av[]) {
  output_dir = "output" ;
  Loci::disableDebugDir() ;
  Loci::Init(&ac,&av) ;

  enum {ASCII,TWODGV,ENSIGHT,CGNS,FIELDVIEW,TECPLOT,VTK,VTK_SURFACE,VTK64,VTK_SURFACE64,CUTTINGPLANE, SURFACE, MEAN, COMBINE, FCOMBINE, NONE} plot_type = NONE ;

  string casename ;
  bool found_casename = false ;
  bool found_iteration = false ;
  bool skip_part_dirs = false;
  string iteration ;
  vector<string> variables ;
  vector<string> boundaries ;
  float xShift, yShift, zShift, temp;

  xShift = yShift = zShift = 0.0;
  int view = VIEWXY ;
  affineMapping transformMatrix;

  // record the maximum particle number to extract
  // a value < 0 means that there is no maximum particle
  // number limit, i.e., all particles are to be extracted
  // default is to extract all particles.  users can use
  // command line switch "-mp <n>" to set the maximum number
  // of particles to be extracted. if the requested particle
  // number is larger than the available particle number, then
  // all particles will be extracted.
  int max_particles = -1 ;

  int end_iter = -1 ;
  int inc_iter = -1 ;
  bool id_required = false;//ensight has the option to display node and element ids

  vector<string> partlist ;
  bool novolume = false ;
  for(int i=1;i<ac;++i) {
    if(av[i][0] == '-') {
      if(!strcmp(av[i],"-ascii"))
        plot_type = ASCII ;
      else if(!strcmp(av[i],"-surf"))
        plot_type = SURFACE ;
      else if(!strcmp(av[i],"-mean"))
        plot_type = MEAN ;
      else if(!strcmp(av[i],"-combine"))
        plot_type = COMBINE ;
      else if(!strcmp(av[i],"-fcombine"))
        plot_type = FCOMBINE ;
      else if(!strcmp(av[i],"-2d"))
        plot_type = TWODGV ;
      else if(!strcmp(av[i],"-en"))
        plot_type = ENSIGHT ;
      else if(!strcmp(av[i],"-en_with_id"))
        {
          plot_type = ENSIGHT ;
          id_required = true;
        }
      else if(!strcmp(av[i],"-cgns"))
        {
          plot_type = CGNS ;
          id_required = true;
        }
      else if(!strcmp(av[i],"-fv"))
        plot_type = FIELDVIEW ;
      else if(!strcmp(av[i],"-tec")) {
        plot_type = TECPLOT ;
#ifndef USE_NATIVE_TECPLOT
	cerr << "Note, This compiled version is using the older ASCII tecplot format." << endl 
	     << "If you are using a recent version you can configure Loci to use" << endl 
	     << "the native binary tecplot format that will be more effective for use with" << endl
	     << "tecplot360. " << endl ;
	  
#endif
      } 
      else if(!strcmp(av[i],"-vtk"))
        plot_type = VTK ;
      else if(!strcmp(av[i],"-vtk_surf"))
        plot_type = VTK_SURFACE ;
      else if(!strcmp(av[i],"-vtk64"))
        plot_type = VTK64 ;
      else if(!strcmp(av[i],"-vtk_surf64"))
        plot_type = VTK_SURFACE64 ;
      else if(!strcmp(av[i],"-cut"))
	plot_type = CUTTINGPLANE ;
      else if(!strcmp(av[i],"-Sx")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> xShift).fail())
          Usage(ac, av) ;
      }
      else if(!strcmp(av[i],"-Sy")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> yShift).fail())
          Usage(ac, av) ;
      }
      else if(!strcmp(av[i],"-Sz")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> zShift).fail())
          Usage(ac, av) ;
      }
      else if(!strcmp(av[i],"-Rx")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> temp).fail())
          Usage(ac, av) ;
        transformMatrix.rotateX(-temp) ;
      }
      else if(!strcmp(av[i],"-Ry")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> temp).fail())
          Usage(ac, av) ;
	transformMatrix.rotateY(-temp) ;
      }
      else if(!strcmp(av[i],"-Rz")) {
	i++ ;
	std::istringstream iss(av[i]) ;
	if ((iss >> std::dec >> temp).fail())
          Usage(ac, av) ;
	transformMatrix.rotateZ(-temp) ;
      }
      else if(!strcmp(av[i],"-xy")) 
        view=VIEWXY ;
      else if(!strcmp(av[i],"-yz")) {
        view=VIEWYZ ;
	transformMatrix.rotateY(90.0) ;
	transformMatrix.rotateZ(90.0) ;
      }
      else if(!strcmp(av[i],"-xz")) {
        view=VIEWXZ ;
	transformMatrix.rotateX(90.0) ;
      }
      else if(!strcmp(av[i],"-xr")) 
        view=VIEWXR ;
      else if(!strcmp(av[i],"-novolume"))
	novolume = true ;
      else if(!strcmp(av[i],"-bc")) {
        i++ ;
        string v(av[i]) ;
        if(av[i][0] >='0' && av[i][0] <= '9')
          v = "BC_"+v ;
        boundaries.push_back(v) ;
	partlist.push_back(v) ;
      } else if(!strcmp(av[i],"-part")) {
        i++ ;
        string v(av[i]) ;
        partlist.push_back(v) ;
      } else if(!strcmp(av[i],"-mp")) {
        // get the number of particles
        ++i ;
        string n(av[i]) ;
        if(!valid_int(n)) {
          cerr << "argument followed option '-mp' is not an integer,"
               << " used default value" << endl ;
        } else {
          max_particles = str2int(n) ;
        }
        
      } else if(!strcmp(av[i],"-dir")) {
        ++i ;
        output_dir = string(av[i]) ;
      } else if(!strcmp(av[i],"-inc")) {
        ++i ;
        inc_iter = atoi(av[i]) ;
      } else if(!strcmp(av[i],"-end")) {
        ++i ;
        end_iter = atoi(av[i]) ;
      } else if(!strcmp(av[i],"-skipPartDirs") ||
		!strcmp(av[i],"-onlyvolume") ) {
        skip_part_dirs = true ;
      } else {
        cerr << "unknown option " << av[i] << endl ;
        Usage(ac,av) ;
      }
      
    } else {
      if(found_iteration)
        variables.push_back(string(av[i])) ;
      else if(found_casename) {
        iteration = string(av[i]) ;
        found_iteration = true ;
      } else {
        casename = string(av[i]) ;
        found_casename = true ;
      }
    }
  }
  if(boundaries.size() > 0)
    novolume = true ;
  if(plot_type == NONE) {
    Usage(ac,av) ;
  }

  // if output directory doesn't exist, create one
  struct stat statbuf ;
  if(stat(output_dir.c_str(),&statbuf))
    mkdir(output_dir.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file '" << output_dir
           <<"' should be a directory!, rename 'output' and start again."
           << endl ;
      Loci::Abort() ;
    }
  
  if(partlist.size() == 0) { // scan for parts
    std::set<string> partfind ;
    DIR *dp = opendir(output_dir.c_str()) ;
    // Look in output directory and find all variables
    if(dp == 0) {
      cerr << "unable to open directory '" << output_dir << "'" << endl ;
      exit(-1) ;
    }
    dirent *entry = readdir(dp) ;
  
    while(entry != 0) {
      string filename = entry->d_name ;
      string postfix ;
      string vname ;
      string vtype ;
      string searchHeader = casename+"_SURF." ;
      string filesub = filename.substr(0,searchHeader.size()) ;
      if(searchHeader == filesub && !skip_part_dirs) {
	size_t len = filename.size()-searchHeader.size() ;
	string partname = filename.substr(searchHeader.size(),len) ;
	string filecheck = output_dir + "/" + filename + "/topo_file." + iteration ;
	struct stat tmpstat ;
	if(stat(filecheck.c_str(),&tmpstat)== 0) {
	  partfind.insert(partname) ;
	  partlist.push_back(partname) ;
	}
      }
      entry = readdir(dp) ;
    }
    closedir(dp) ;

    vector<string> vsurfs = volumeSurfaceNames(output_dir,iteration,
					       casename) ;
    for(size_t i=0;i<vsurfs.size();++i) {
      if(partfind.find(vsurfs[i]) == partfind.end()) {
	partlist.push_back(vsurfs[i]) ;
      }
    }
    
  }    
      
  if(variables.size() == 0) {
    std::set<string> varset ;
    DIR *dp = opendir(output_dir.c_str()) ;
    // Look in output directory and find all variables
    if(dp == 0) {
      cerr << "unable to open directory '" << output_dir << "'" << endl ;
      exit(-1) ;
    }
    dirent *entry = readdir(dp) ;
  
    string tail = iteration + "_" + casename ;
    while(entry != 0) {
      string filename = entry->d_name ;
      string postfix ;
      string vname ;
      string vtype ;
      int nsz = filename.size() ;
      int dot = -1 ;
      for(int i=nsz-1;i>=0;--i)
        if(filename[i] == '.') {
          dot = i ;
          break ;
        }
      for(int i=dot+1;i<nsz;++i)
        postfix += filename[i] ;
      int und = -1 ;
      if(dot > 0) {
        for(int i=dot-1;i>=0;--i)
          if(filename[i] == '_') {
            und = i ;
            break ;
          }
        if(und > 0) {
          for(int i=und+1;i<dot;++i)
            vtype += filename[i] ;
          for(int i=0;i<und;++i)
            vname += filename[i] ;
        } else
          for(int i=0;i<dot;++i)
            vname += filename[i] ;
        
        
      }
      //add derived variables  
      if(dot>0 && und>0 && postfix == tail) {
        if(vtype == "sca" || vtype == "vec" || vtype == "bnd" ||
           vtype == "bndvec" || vtype == "ptsca" || vtype == "ptvec") {
          varset.insert(vname) ;
        }
        if(vtype == "sca" && vname == "pg") {
          varset.insert("P") ;
        }
        if(vtype == "sca" && vname == "a") {
          varset.insert("m") ;
        }
      }
       
      entry = readdir(dp) ;
    }
    closedir(dp) ;
    
    if(partlist.size() > 0) { // Now check each part for variables
      for(size_t i=0;i<partlist.size();++i) {
	string dirname = output_dir+"/"+casename+"_SURF."+partlist[i] ;
	DIR *dp = opendir(dirname.c_str()) ;
	// Look in output directory and find all variables
	if(dp == 0 || skip_part_dirs) {
	  continue ;
	}
	dirent *entry = readdir(dp) ;
  	for(;entry != 0;entry=readdir(dp)) {
	  string filename = entry->d_name ;	
	  int fsz = filename.size() ;
	  int isz = iteration.size() ;
	  if(fsz <= isz)
	    continue ;
	  string fiter = filename.substr(fsz-(isz+1),isz+1) ;

	  if(fiter != string("."+iteration)) 
	    continue ;

	  string remainder  = filename.substr(0,fsz-(isz+1)) ;
	  int remsz = remainder.size() ;
	  if(remsz <= 4)
	    continue ;
	  string postfix = remainder.substr(remsz-4,4) ;
	  if(postfix == "_sca" || postfix == "_vec") {
	    string vname = remainder.substr(0,remsz-4) ;
	    varset.insert(vname) ;
	    continue ;
	  }
	  if(remsz <= 5)
	    continue ;
	  postfix = remainder.substr(remsz-5,5) ;
	  if(postfix == "_bsca" || postfix == "_bvec") {
	    string vname = remainder.substr(0,remsz-5) ;
	    varset.insert(vname) ;
	  }
	}
	closedir(dp) ;
      }
    }
    
    std::set<string>::const_iterator vi ;
    for(vi=varset.begin();vi!=varset.end();++vi)
      variables.push_back(*vi) ;

    if(variables.size() == 0) {
      variables.push_back("x") ;
      variables.push_back("y") ;
      variables.push_back("z") ;
    }
    cout << "extracting variables:" ;
    for(size_t i=0;i<variables.size();++i) {
      cout << ' ' << variables[i] ;
    }
    cout << endl ;
      
  }

  //find out variable types and variable files
  vector<int> variable_type(variables.size()) ;
  vector<string> variable_file(variables.size()) ;

  bool particle_info_requested = false ;
  
  for(size_t i=0;i<variables.size();++i) {
    const string var(variables[i]) ;
    string filename = output_dir+'/' + var + "_hdf5." + iteration ;
    struct stat tmpstat ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }
      
    filename = output_dir+'/' + var + "_sca." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }

    filename = output_dir+'/' + var + "_vec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_VECTOR ;
      variable_file[i] = filename ;
      continue ;
    }
    filename = output_dir+'/' + var + "_bnd." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = BOUNDARY_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }
    filename = output_dir+'/' + var + "_bndvec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = BOUNDARY_VECTOR ;
      variable_file[i] = filename ;
      continue ;
    }

    filename = output_dir+'/' + var + "_ptsca." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)==0) {
      variable_type[i] = PARTICLE_SCALAR ;
      variable_file[i] = filename ;
      particle_info_requested = true ;
      continue ;
    }

    filename = output_dir+'/' + var + "_ptvec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)==0) {
      variable_type[i] = PARTICLE_VECTOR ;
      variable_file[i] = filename ;
      particle_info_requested = true ;
      continue ;
    }

    if(plot_type == ASCII && boundaries.size() > 0) {
      if(var == "n") {
        variable_type[i] = BOUNDARY_DERIVED_VECTOR ;
        continue ;
      }
      if(var == "area" || var == "x" || var == "y" || var == "z") {
        variable_type[i] = BOUNDARY_DERIVED_SCALAR ;
        continue ;
      }
    }
    if(var == "m") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "p") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "P") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "u") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "0" || var == "1" || var == "2") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "x" || var == "y" || var == "z") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
      
    // ... other derived variables here

    if(partlist.size() == 0) 
      cerr << "Warning, variable '" << var << "' is unknown and will not be processed." << endl ;
    variable_type[i] = UNDEFINED ;
  }


  // we will first check to see if particle position is present
  // in case of any particle information extraction
  if(particle_info_requested) {
    string filename = output_dir +"/particle_pos." + iteration + "_" + casename ;
    struct stat tmpstat ;
    if(stat(filename.c_str(),&tmpstat)!=0) {
      cerr << "Warning: particle geometry '" << filename << "' must"
           << " be presented for any"
           << " particle related variable extraction." << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
  }
#ifdef H5_USE_16_API
  H5Eset_auto(NULL,NULL) ;
#else
  H5Eset_auto(H5E_DEFAULT,NULL,NULL) ;
#endif

  string filename = getTopoFileName(output_dir, casename, iteration) ;
  struct stat tmpstat ;
  string posfile = getPosFile(output_dir,iteration,casename) ;

  bool timesyncerror = false ;
  if(stat(filename.c_str(),&tmpstat)==0) {
    struct stat gridstat ;
    string gridfile = casename + ".vog" ;
    if(stat(gridfile.c_str(),&gridstat)==0) {
      if(gridstat.st_mtime > tmpstat.st_mtime)
        timesyncerror = true ;
    }
  }

  if(timesyncerror) {
    cerr << "WARNING!!!:  grid file newer than topology file in output directory!  "
         << endl
         << "             You are not extracting the present state of the mesh!"
         << endl
         << "             Rerun chem or vogcheck to regenerate mesh topology file in "
         << endl
         << "             output directory."
         << endl ;
  }

  if(plot_type == ASCII) {
    if(boundaries.size() == 0) {
      process_ascii_nodal(casename,iteration,
                          variables,variable_type,variable_file) ;
      Loci::Finalize() ;
      exit(0) ;
    } else {
      process_ascii_bndry(casename,iteration,
                          variables,variable_type,variable_file,
                          boundaries) ;
      Loci::Finalize() ;
      exit(0) ;
    }
  }

  if(plot_type == MEAN) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
	   << "       and option -inc to specify iteration increment value for iterations" << endl
	   << "       to specify which files to average!" << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
	
    process_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == COMBINE) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
	   << "       and option -inc to specify iteration increment value for iterations" << endl
	   << "       to specify which files to average!" << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
	
    combine_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter,false) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == FCOMBINE) {
    if(end_iter<0 ||inc_iter< 0) {
      cerr << "ERROR: Must use option -end to specify ending iteration for average" << endl
	   << "       and option -inc to specify iteration increment value for iterations" << endl
	   << "       to specify which files to average!" << endl ;
      Loci::Finalize() ;
      exit(-1) ;
    }
	
    combine_mean(casename,iteration,variables,variable_type,
                 variable_file,end_iter,inc_iter,true) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == TWODGV) {
    if(variables.size() != 1) {
      cerr << "2dgv extract can only extract one variable at a time."
           << endl ;
      Usage(ac,av) ;
    }
    if(boundaries.size() == 0) {
      cerr << "2dgv extract must have the projected boundaries identified using the '-bc' flag" << endl ;
      Usage(ac,av) ;
    }
    get_2dgv(casename,iteration,variables,variable_type,variable_file,
             boundaries,view) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  if(plot_type == SURFACE) {
    if(boundaries.size() ==0) {
      cerr << "'extract -surf' must have one boundary surface identified using the '-bc' flag" << endl ;
      Usage(ac,av) ;
    }
    if(iteration == "") {
      cerr << "'extract -surf' must specify iteration to extract from"
           << endl ;
      Usage(ac,av) ;
    }
    get_surf(casename,iteration,variables,variable_type,variable_file,
             boundaries) ;
    Loci::Finalize() ;
    exit(0) ;
  }
  
  // New grid topology processor
  postProcessorP postprocessor = 0 ;
  
  switch(plot_type) {
  case ENSIGHT:
    postprocessor = new ensightPartConverter(id_required) ;
    break ;
  case CGNS:
    postprocessor = new cgnsPartConverter(id_required) ;
    break ;  
  case FIELDVIEW:
    postprocessor = new fieldViewPartConverter ;
    break ;
  case TECPLOT:
    postprocessor = new tecplotPartConverter ;
    break ;
  case VTK:
    postprocessor = new vtkPartConverter(false) ;
    break ;
  case VTK_SURFACE:
    postprocessor = new vtkSurfacePartConverter(false) ;
    break ;
  case VTK64:
    postprocessor = new vtkPartConverter(true) ;
    break ;
  case VTK_SURFACE64:
    postprocessor = new vtkSurfacePartConverter(true) ;
    break ;
  case CUTTINGPLANE:
    postprocessor = new cuttingPlanePartConverter(transformMatrix, -xShift, -yShift, -zShift) ;
    break ;
  default:
    cerr << "Unknown export method!" << endl ;
    break ;
  }

  if(postprocessor != 0) {
#ifdef H5_USE_16_API
    H5Eset_auto(NULL,NULL) ;
#else
    H5Eset_auto(H5E_DEFAULT,NULL,NULL) ;
#endif
    vector<surfacePartP> parts ;


    // Check for derived variable input requirements
    std::set<string> varset ;
    for(size_t i=0;i<variables.size();++i)
      varset.insert(variables[i]) ;
    if(varset.find("m") != varset.end()) {
      varset.insert("a") ;
      varset.insert("v") ;
    }
    if(varset.find("P") != varset.end()) {
      varset.insert("pg") ;
    }
    if(varset.find("p") != varset.end()) {
      varset.insert("pg") ;
    }
    if(varset.find("u") != varset.end()) {
      varset.insert("v") ;
    }
    if(varset.find("0") != varset.end()) {
      varset.insert("v") ;
    }
    if(varset.find("1") != varset.end()) {
      varset.insert("v") ;
    }
    if(varset.find("2") != varset.end()) {
      varset.insert("v") ;
    }
    vector<string> varlist ;
    std::set<string>::const_iterator vi ;
    for(vi=varset.begin();vi!=varset.end();++vi)
      varlist.push_back(*vi) ;
    variables.swap(varlist) ;
    
    if(postprocessor->processesSurfaceElements()) {
      std::set<string> volsearch ;
      for(size_t i=0;i<partlist.size();++i) {
	string name = partlist[i] ;
	string dir = output_dir + "/" + casename + "_SURF." + name ;
        if (skip_part_dirs) {
	  volsearch.insert(name) ;
        } else {
	  surfacePartP sp = new surfacePart(name,dir,iteration,variables) ;
	  if(sp->fail()) {
	    volsearch.insert(name) ;
	  } else {
            cout << "part: " << name << endl ;
	    parts.push_back(sp) ;
	  }
        }
      }
      if(!volsearch.empty()) {
	vector<surfacePartP> volSurface ;

	volumePartP vp = 
	  new volumePart(output_dir,iteration,casename,variables) ;	
	if(!vp->fail())
	  extractVolumeSurfaces(volSurface,vp,output_dir,iteration,
				casename,variables) ;
	
	for(size_t i=0;i<volSurface.size();++i) {
	  if(volsearch.find(volSurface[i]->getPartName()) != volsearch.end()) {
	    cout << "part: " << volSurface[i]->getPartName() << endl ;
	    parts.push_back(volSurface[i]) ;
	  }
	}
      }
    }
    volumePartP vp = 0 ;
    if(!novolume && postprocessor->processesVolumeElements()) {
      string testfile = getTopoFileName(output_dir, casename, iteration) ;
      struct stat tmpstat ;
      cout << "checking " << testfile << endl ;
      if(stat(testfile.c_str(),&tmpstat)==0) {
	// topo file exists, so there is a volume grid
	cout << "creating volume part" << endl;
	vp = new volumePart(output_dir,iteration,casename,variables) ;
	if(vp->fail()) {
	  vp = 0 ;
	} 
      }
    }
    particlePartP pp = 0 ;
    if(postprocessor->processesParticleElements()) {
      string testfile = output_dir + "/particle_pos."+iteration + "_" + casename ;
      struct stat tmpstat ;
      cout << "checking " << testfile << endl ;
      if(stat(testfile.c_str(),&tmpstat)==0) {
	pp = new particlePart(output_dir,iteration,casename,variables,max_particles) ;
      }
    }
    
    if(parts.size() > 0) {
      vector<surfacePartP> modparts(parts.size()) ;
      for(size_t i=0;i<parts.size();++i)
	modparts[i] = new surfacePartDerivedVars(parts[i],
                                                 output_dir,casename,
                                                 iteration, variables) ;
      postprocessor->addSurfaceParts(modparts) ;
    }
    
    if(vp!=0) {
      volumePartP vpn = new volumePartDerivedVars(vp,
						  output_dir,casename,
						  iteration,variables) ;
      postprocessor->addVolumePart(vpn) ;
    }
    if(pp!=0)
      postprocessor->addParticlePart(pp) ;
    postprocessor->exportPostProcessorFiles(casename,iteration) ;
    Loci::Finalize() ;
    exit(0) ;
  } 


  if(timesyncerror) {
    cerr << "WARNING!!!:  grid file newer than topology file in output directory!  "
         << endl
         << "             You are not extracting the present state of the mesh!"
         << endl
         << "             Rerun chem or vogcheck to regenerate mesh topology file in "
         << endl
         << "             output directory."
         << endl ;
  }
  Loci::Finalize() ;
}

