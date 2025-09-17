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
#include <Loci.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <set> 
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::set ;
using std::ofstream ;
using std::ios ;

#include "extract.h"

#include <sys/types.h>
#include <sys/stat.h>


bool ensightPartConverter::processesVolumeElements() const {
  return true ; 
}
bool ensightPartConverter::processesSurfaceElements() const {
  return true ; 
}
bool ensightPartConverter::processesParticleElements() const {
  return true ; 
}

void ensightPartConverter::exportPostProcessorFiles(string casename,
						    string iteration) const {
  string dirname = casename + "_ensight."+iteration ;
  struct stat statbuf ;
  if(stat(dirname.c_str(),&statbuf))
    mkdir(dirname.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file " << dirname << " should be a directory!, rename 'ensight' and start again."
           << endl ;
      exit(-1) ;
    }
  
  string geo_filename =  casename + ".geo" ;
  string case_filename = dirname + "/" + casename + ".case" ;
  set<string> nodal_scalars ;
  set<string> nodal_vectors ;
  set<string> element_scalars ;
  set<string> element_vectors ;
  set<string> particle_scalars ;
  set<string> particle_vectors ;
  for(size_t i=0;i<volumePartList.size();++i) {
    vector<string> nscalars = volumePartList[i]->getNodalScalarVars() ;
    for(size_t j=0;j<nscalars.size();++j) 
      nodal_scalars.insert(nscalars[j]) ;
    vector<string> nvectors = volumePartList[i]->getNodalVectorVars() ;
    for(size_t j=0;j<nvectors.size();++j) 
      nodal_vectors.insert(nvectors[j]) ;
  }    
  for(size_t i =0;i<surfacePartList.size();++i) {
    vector<string> nscalars = surfacePartList[i]->getNodalScalarVars() ;
    for(size_t j=0;j<nscalars.size();++j) 
      nodal_scalars.insert(nscalars[j]) ;
    vector<string> nvectors = surfacePartList[i]->getNodalVectorVars() ;
    for(size_t j=0;j<nvectors.size();++j) 
      nodal_vectors.insert(nvectors[j]) ;
    vector<string> escalars = surfacePartList[i]->getElementScalarVars() ;
    for(size_t j=0;j<escalars.size();++j) 
      element_scalars.insert(escalars[j]) ;
    vector<string> evectors = surfacePartList[i]->getElementVectorVars() ;
    for(size_t j=0;j<evectors.size();++j) 
      element_vectors.insert(evectors[j]) ;
  }
  for(size_t i=0;i<particlePartList.size();++i) {
    vector<string> nscalars = particlePartList[i]->getScalarVars() ;
    for(size_t j=0;j<nscalars.size();++j) 
      particle_scalars.insert(nscalars[j]) ;
    vector<string> nvectors = particlePartList[i]->getVectorVars() ;
    for(size_t j=0;j<nvectors.size();++j) 
      particle_vectors.insert(nvectors[j]) ;
  }    


  //write out case file
  ofstream of(case_filename.c_str(),ios::out) ;
  of << "FORMAT" << endl ;
  of << "type:  ensight gold" << endl ;
  of << "GEOMETRY" << endl ;
  of << "model:  " << geo_filename << endl ;
  string particle_geo_file ;
  if(particlePartList.size() > 0) {
    string pgeo_file = casename + "_particles.geo" ;
    of << "measured:  " << pgeo_file << endl ;
    particle_geo_file = dirname + "/" + pgeo_file ;
  }
  
  geo_filename = dirname + "/"+geo_filename ;
  
  if(nodal_scalars.size()+nodal_vectors.size()+
     element_scalars.size()+element_vectors.size()+
     particle_scalars.size()+particle_vectors.size()> 0) {
    of << "VARIABLE" << endl ;
    set<string>::const_iterator si ;
    for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {
      of << "scalar per node:\t " << *si << '\t' << *si << endl ;
    }
    for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
      of << "vector per node:\t " << *si << '\t' << *si << endl ;
    }
    for(si=element_scalars.begin();si!=element_scalars.end();++si) {
      of << "scalar per element:\t " << *si << '\t' << *si << endl ;
    }
    for(si=element_vectors.begin();si!=element_vectors.end();++si) {
      of << "vector per element:\t " << *si << '\t' << *si << endl ;
    }
    for(si=particle_scalars.begin();si!=particle_scalars.end();++si) {
      of << "scalar per measured node:\t " << *si << '\t' << *si << endl ;
    }
    for(si=particle_vectors.begin();si!=particle_vectors.end();++si) {
      of << "vector per measured node:\t " << *si << '\t' << *si << endl ;
    }
  }

  of.close() ;
  
  FILE *OFP ;
  OFP = fopen(geo_filename.c_str(),"wb") ;
  if(OFP == NULL) {
    cerr << "unable to open file '" << geo_filename << "' for writing!" ;
    exit(-1) ;
  }
  char tmp_buf[80] ;
  
  memset(tmp_buf, '\0', 80) ; 
  snprintf(tmp_buf,80, "C Binary") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf,80, "Ensight model geometry description") ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "Grid file used is %s", casename.c_str()) ;
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  snprintf(tmp_buf, 80, "node id off") ;  
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;
  memset(tmp_buf, '\0', 80) ;
  if(id_required)snprintf(tmp_buf,80, "element id given") ;
  else snprintf(tmp_buf, 80,"element id off") ;  
  fwrite(tmp_buf, sizeof(char), 80, OFP) ;

  int part_id = 0 ;
  // Volume Parts Output
  // -----------------------------------------------------------------------
  vector<int> vpartnums ;
  for(size_t i=0;i<volumePartList.size();++i) { 
    part_id++ ;
    vpartnums.push_back(part_id) ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf, 80, "part") ;
    fwrite(tmp_buf,sizeof(char), 80, OFP) ;
    fwrite(&part_id,sizeof(int),1,OFP) ;
    memset(tmp_buf, '\0', 80) ;
    string name = volumePartList[i]->getPartName();
    snprintf(tmp_buf,80, "%s", name.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf, 80, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    size_t pnts = volumePartList[i]->getNumNodes() ;
    vector<vector3d<float> > pos ;
    volumePartList[i]->getPos(pos) ;

    int pts = pnts ;
    fwrite(&pts,sizeof(int),1,OFP) ;
    for(size_t j=0;j<pnts;++j) {
      float x = pos[j].x ;
      fwrite(&x,sizeof(float),1,OFP) ;
    }
    for(size_t j=0;j<pnts;++j) {
      float y = pos[j].y ;
      fwrite(&y,sizeof(float),1,OFP) ;
    }
    for(size_t j=0;j<pnts;++j) {
      float z = pos[j].z ;
      fwrite(&z,sizeof(float),1,OFP) ;
    }
    const int block_size=65536 ; // Size of blocking factor
    
    if(volumePartList[i]->getNumTets() > 0) { // write out tets
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "tetra4") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = volumePartList[i]->getNumTets() ;
      fwrite(&tot, sizeof(int), 1, OFP) ;
      
      if(id_required){//write fake ids
        vector<int> ids(tot, 1);
        fwrite(&ids[0],sizeof(int),tot,OFP); 
      }
      
      int start = 0 ;
      int top = tot + volumePartList[i]->getNumTetsIblank() ;
      while(start < top) {
	vector<Array<int,4> > tets ;
	volumePartList[i]->getTetBlock(tets,start,block_size) ;
	start += block_size ;
	if(tets.size() > 0) {
	  int ntets = tets.size() ;
	  fwrite(&tets[0],sizeof(Array<int,4>),ntets,OFP) ;
	}
      }
    }

    if(volumePartList[i]->getNumHexs() > 0) { // write out hexs
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80, "hexa8") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = volumePartList[i]->getNumHexs() ;
      fwrite(&tot, sizeof(int), 1, OFP) ;
      
      if(id_required){//write fake ids
        vector<int> ids(tot, 1);
        fwrite(&ids[0],sizeof(int),tot,OFP); 
      }
      
      int start = 0 ;
      int top = tot+volumePartList[i]->getNumHexsIblank() ;
      while(start < top) {
	vector<Array<int,8> > hexs ;
	volumePartList[i]->getHexBlock(hexs,start,block_size) ;
	start += block_size ;
	if(hexs.size() > 0) {
	  int nhexs = hexs.size() ;
	  fwrite(&hexs[0],sizeof(Array<int,8>),nhexs,OFP) ;
	}
      }
    }

    if(volumePartList[i]->getNumPrsm() > 0) { // write out Prisms
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80, "penta6") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = volumePartList[i]->getNumPrsm() ;
      fwrite(&tot, sizeof(int), 1, OFP) ;

      if(id_required){//write fake ids
        vector<int> ids(tot, 1);
        fwrite(&ids[0],sizeof(int),tot,OFP); 
      }
      
      int start = 0 ;
      int top = tot + volumePartList[i]->getNumPrsmIblank() ;
      while(start < top) {
	vector<Array<int,6> > prsm ;
	volumePartList[i]->getPrsmBlock(prsm,start,block_size) ;
	start += block_size ;
	if(prsm.size() > 0) {
	  int nprsm = prsm.size() ;
	  fwrite(&prsm[0],sizeof(Array<int,6>),nprsm,OFP) ;
	}
      }
    }

    if(volumePartList[i]->getNumPyrm() > 0) { // write out Pyramids
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf,80, "pyramid5") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int tot = volumePartList[i]->getNumPyrm() ;
      fwrite(&tot, sizeof(int), 1, OFP) ;

      if(id_required){//write fake ids
        vector<int> ids(tot, 1);
        fwrite(&ids[0],sizeof(int),tot,OFP); 
      }
      
      int start = 0 ;
      int top = tot+volumePartList[i]->getNumPyrmIblank() ;
      while(start < top) {
	vector<Array<int,5> > pyrm ;
	volumePartList[i]->getPyrmBlock(pyrm,start,block_size) ;
	start += block_size ;
	if(pyrm.size() > 0) {
	  int npyrm = pyrm.size() ;
	  fwrite(&pyrm[0],sizeof(Array<int,5>),npyrm,OFP) ;
	}
      }
    }
    
    if(volumePartList[i]->getNumGenc() > 0) { // write out general cells
      vector<int> genCellNfaces, genCellNsides,genCellNodes ;
      volumePartList[i]->getGenCell(genCellNfaces,genCellNsides,genCellNodes) ;
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "nfaced") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int nnf = genCellNfaces.size() ;
      fwrite(&nnf, sizeof(int), 1, OFP) ;

      if(id_required){//write fake ids
        vector<int> ids(nnf, 1);
        fwrite(&ids[0],sizeof(int),nnf,OFP); 
      }
      
      fwrite(&genCellNfaces[0], sizeof(int), nnf, OFP) ;
      int nnsides = genCellNsides.size() ;
      fwrite(&genCellNsides[0], sizeof(int),nnsides,OFP) ;
      int nnodes = genCellNodes.size() ;
      fwrite(&genCellNodes[0],sizeof(int),nnodes,OFP) ;
    }
  }
  vector<int> partnums(surfacePartList.size()) ;
  // Loop over parts, write out their geometry
  for(size_t i =0;i<surfacePartList.size();++i) {
    part_id++ ;
    partnums[i] = part_id ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "part") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    fwrite(&part_id, sizeof(int), 1, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    string name = surfacePartList[i]->getPartName();
    snprintf(tmp_buf,80, "%s", name.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, OFP) ;
    int npt = surfacePartList[i]->getNumNodes() ;
    fwrite(&npt, sizeof(int), 1, OFP) ;
    vector<vector3d<float> > pos ;
    surfacePartList[i]->getPos(pos) ;
    for(int j=0;j<npt;++j) {
      float x = pos[j].x ;
      fwrite(&x,sizeof(float),1,OFP) ;
    }
    for(int j=0;j<npt;++j) {
      float y = pos[j].y ;
      fwrite(&y,sizeof(float),1,OFP) ;
    }
    for(int j=0;j<npt;++j) {
      float z = pos[j].z ;
      fwrite(&z,sizeof(float),1,OFP) ;
    }
    int nquads = surfacePartList[i]->getNumQuads() ;
    if(nquads > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "quad4") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int nq = nquads ;
      fwrite(&nq, sizeof(int),1,OFP) ;
      if(id_required){
        vector<int> quads_ids;
        surfacePartList[i]->getQuadsIds(quads_ids);
        fwrite(&quads_ids[0], sizeof(int),nquads,OFP) ;
      }
      vector<Array<int,4> > quads ;
      surfacePartList[i]->getQuads(quads) ;
      fwrite(&quads[0],sizeof(Array<int,4>),nquads,OFP) ;
    }
    int ntrias = surfacePartList[i]->getNumTrias() ;
    if(ntrias > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "tria3") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int nt = ntrias ;
      fwrite(&nt, sizeof(int),1,OFP) ;
      if(id_required){
        vector<int> trias_ids;
        surfacePartList[i]->getTriasIds(trias_ids);
        fwrite(&trias_ids[0], sizeof(int),ntrias,OFP) ;
      }
      vector<Array<int,3> > trias ;
      surfacePartList[i]->getTrias(trias) ;
      fwrite(&trias[0],sizeof(Array<int,3>),ntrias,OFP) ;
    }    
    int ngeneral = surfacePartList[i]->getNumGenfc() ;
    if(ngeneral > 0) {
      char tmp_buf[80] ;
      memset(tmp_buf, '\0', 80) ;
      snprintf(tmp_buf, 80, "nsided") ;
      fwrite(tmp_buf, sizeof(char), 80, OFP) ;
      int ngen = ngeneral ;
      fwrite(&ngen, sizeof(int),1,OFP) ;
      if(id_required){
        vector<int> nside_ids;
        surfacePartList[i]->getGenfIds(nside_ids);
        fwrite(&nside_ids[0], sizeof(int),ngeneral,OFP) ;
      }
      vector<int> nside_sizes,nside_nodes ;
      surfacePartList[i]->getGenf(nside_sizes,nside_nodes) ;
      
      int tot = 0 ;
      for(int j=0;j<ngeneral;++j)
	tot += nside_sizes[j] ;
      int nside_nodes_size = nside_nodes.size() ;
      if(nside_nodes_size != tot) {
	cerr << "mismatch in node size and faces size " << nside_nodes_size
	     << " was " << tot << endl ;
      }
      fwrite(&nside_sizes[0],sizeof(int),ngeneral,OFP) ;
      fwrite(&nside_nodes[0],sizeof(int),nside_nodes_size,OFP) ;
    }
  }

  fclose(OFP) ;
  
  // Finished writing out the the geo file, now write out the variables
  set<string>::const_iterator si ;
  // write out nodal scalars
  for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;

      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    // Loop over volume parts and write out variables for each part if they
    // exist
    for(size_t i =0;i<volumePartList.size();++i) {
      if(volumePartList[i]->hasNodalScalarVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = vpartnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "coordinates") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	vector<float> vals ;
	volumePartList[i]->getNodalScalar(varname,vals) ;
	fwrite(&vals[0],sizeof(float),vals.size(),FP) ;
      }
    }
    // Loop over surface parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      if(surfacePartList[i]->hasNodalScalarVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "coordinates") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	vector<float> vals ;
	surfacePartList[i]->getNodalScalar(varname,vals) ;
	fwrite(&vals[0],sizeof(float),vals.size(),FP) ;
      }
    }
    fclose(FP) ;
  }
  // write out nodal vector variables
  for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;

      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    // Loop over volume parts and write out variables for each part if they
    // exist
    for(size_t i =0;i<volumePartList.size();++i) {
      if(volumePartList[i]->hasNodalVectorVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = vpartnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "coordinates") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	vector<vector3d<float> > vals ;
	volumePartList[i]->getNodalVector(varname,vals) ;
        int nvals = vals.size() ;
        for(int j=0;j<nvals;++j)
          fwrite(&vals[j].x,sizeof(float),1,FP) ;
        for(int j=0;j<nvals;++j)
          fwrite(&vals[j].y,sizeof(float),1,FP) ;
        for(int j=0;j<nvals;++j)
          fwrite(&vals[j].z,sizeof(float),1,FP) ;
      }
    }
    // Loop over parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      memset(tmp_buf, '\0', 80) ;
      if(surfacePartList[i]->hasNodalVectorVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "coordinates") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	vector<vector3d<float> > vals ;
	surfacePartList[i]->getNodalVector(varname,vals) ;
	int nvals = vals.size() ;
	for(int i=0;i<nvals;++i)
	  fwrite(&vals[i].x,sizeof(float),1,FP) ;

	for(int i=0;i<nvals;++i)
	  fwrite(&vals[i].y,sizeof(float),1,FP) ;

	for(int i=0;i<nvals;++i)
	  fwrite(&vals[i].z,sizeof(float),1,FP) ;
      }
    }
    fclose(FP) ;
  }
  // write out element scalars
  for(si=element_scalars.begin();si!=element_scalars.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;
      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    // Loop over parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      if(surfacePartList[i]->hasElementScalarVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
        vector<float> qvals, tvals, gvals ;
        surfacePartList[i]->getElementScalar(varname,qvals,tvals,gvals) ;

        if(qvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "quad4") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&qvals[0],sizeof(float),qvals.size(),FP) ;
        }
        if(tvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "tria3") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&tvals[0],sizeof(float),tvals.size(),FP) ;
        }
        if(gvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "nsided") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          fwrite(&gvals[0],sizeof(float),gvals.size(),FP) ;
        }
      }
    }
    fclose(FP) ;
  }
  // write out element vectors
  for(si=element_vectors.begin();si!=element_vectors.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;
      continue ;
    }
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"variable : %s",varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    
    // Loop over parts and write out variables for each part if they 
    // exist ;
    for(size_t i =0;i<surfacePartList.size();++i) {
      if(surfacePartList[i]->hasElementVectorVar(varname)) {
	memset(tmp_buf, '\0', 80) ;
	snprintf(tmp_buf,80, "part") ;
	fwrite(tmp_buf, sizeof(char), 80, FP) ;
	int tmp = partnums[i] ;
	fwrite(&tmp, sizeof(int), 1, FP) ;
        vector<vector3d<float> > qvals, tvals, gvals ;
        surfacePartList[i]->getElementVector(varname,qvals,tvals,gvals) ;
        
        if(qvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "quad4") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          for(size_t i=0;i<qvals.size();++i)
            fwrite(&qvals[i].x,sizeof(float),1,FP) ;
          for(size_t i=0;i<qvals.size();++i)
            fwrite(&qvals[i].y,sizeof(float),1,FP) ;
          for(size_t i=0;i<qvals.size();++i)
            fwrite(&qvals[i].z,sizeof(float),1,FP) ;
        }
        if(tvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "tria3") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          for(size_t i=0;i<tvals.size();++i)
            fwrite(&tvals[i].x,sizeof(float),1,FP) ;
          for(size_t i=0;i<tvals.size();++i)
            fwrite(&tvals[i].y,sizeof(float),1,FP) ;
          for(size_t i=0;i<tvals.size();++i)
            fwrite(&tvals[i].z,sizeof(float),1,FP) ;
        }
        if(gvals.size() > 0) {
          memset(tmp_buf, '\0', 80) ;
          snprintf(tmp_buf,80, "nsided") ;
          fwrite(tmp_buf, sizeof(char), 80, FP) ;
          for(size_t i=0;i<gvals.size();++i)
            fwrite(&gvals[i].x,sizeof(float),1,FP) ;
          for(size_t i=0;i<gvals.size();++i)
            fwrite(&gvals[i].y,sizeof(float),1,FP) ;
          for(size_t i=0;i<gvals.size();++i)
            fwrite(&gvals[i].z,sizeof(float),1,FP) ;
        }
      }
    }
    fclose(FP) ;
  }
  if(particlePartList.size() > 0) {
    FILE *FP = 0 ;
    FP = fopen(particle_geo_file.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << particle_geo_file
	   << "' for writing particle geometry info!" << endl ;
      return ;
    }
    int npnts = 0 ;
    for(size_t i=0;i<particlePartList.size();++i)
      npnts += particlePartList[i]->getNumParticles() ;

    char tmp_buf[80] ;
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80,"C Binary") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "particle positions") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;

    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf,80, "particle coordinates") ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;

    // number of points
    fwrite(&npnts, sizeof(int), 1, FP) ;
    // point ids
    for(int i=1;i<npnts+1;++i)
      fwrite(&i, sizeof(int), 1, FP) ;
    

    for(size_t i=0;i<particlePartList.size();++i) {
      int np = particlePartList[i]->getNumParticles() ;
      vector<vector3d<float> > ppos ;
      particlePartList[i]->getParticlePositions(ppos) ;
      for(int k=0;k<np;++k) {
	float x = ppos[k].x ;
	float y = ppos[k].y ;
	float z = ppos[k].z ;
	fwrite(&x, sizeof(float), 1, FP) ;
	fwrite(&y, sizeof(float), 1, FP) ;
	fwrite(&z, sizeof(float), 1, FP) ;
      }
    }
    fclose(FP) ;
    
  }

  for(si=particle_scalars.begin();si!=particle_scalars.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;
      continue ;
    }
    char tmp_buf[80] ;
    
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf, 80, "Per particle scalar: %s", varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    for(size_t i=0;i<particlePartList.size();++i) {
      int np = particlePartList[i]->getNumParticles() ;
      vector<float> val ;
      particlePartList[i]->getParticleScalar(varname,val) ;
      fwrite(&val[0],sizeof(float),np,FP) ;
    }
    fclose(FP) ;
  }

  for(si=particle_vectors.begin();si!=particle_vectors.end();++si) {
    string varname = *si ;
    string filename = dirname + "/" + varname ;
    FILE *FP = 0 ;
    FP = fopen(filename.c_str(), "wb") ;
    if(FP==0) {
      cerr << "can't open file '" << filename << "' for writing variable info!"
	   << endl ;
      continue ;
    }
    char tmp_buf[80] ;
    
    memset(tmp_buf, '\0', 80) ;
    snprintf(tmp_buf, 80, "Per particle vector: %s", varname.c_str()) ;
    fwrite(tmp_buf, sizeof(char), 80, FP) ;
    for(size_t i=0;i<particlePartList.size();++i) {
      int np = particlePartList[i]->getNumParticles() ;
      vector<vector3d<float> > val ;
      particlePartList[i]->getParticleVector(varname,val) ;
      for(int k=0;k<np;++k) {
	float x = val[k].x ;
	float y = val[k].y ;
	float z = val[k].z ;
    
	fwrite(&x, sizeof(float), 1, FP) ;
	fwrite(&y, sizeof(float), 1, FP) ;
	fwrite(&z, sizeof(float), 1, FP) ;
      }
    }
    fclose(FP) ;
  }

}
// ... the end ...


