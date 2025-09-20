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

#include "extract.h"

#include <sys/types.h>
#include <sys/stat.h>

using namespace std ;

    
bool vtkSurfacePartConverter::processesVolumeElements() const {
  return true ;
}
bool vtkSurfacePartConverter::processesSurfaceElements() const {
  return true ;
}
bool vtkSurfacePartConverter::processesParticleElements() const {
  return true ;
}

void vtkSurfacePartConverter::exportPostProcessorFiles(string casename, string iteration) const {

  FILE *fid;
  string dirname = casename+"_vtk."+iteration ;
  struct stat statbuf ;
  if(stat(dirname.c_str(),&statbuf))
    mkdir(dirname.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file '" << dirname << "' should be a directory!, rename and start again."
           << endl ;
      exit(-1) ;
    }
  
  int nsurfaces = surfacePartList.size() ;
  for(int i=0;i<nsurfaces;++i) {
    string surfName = surfacePartList[i]->getPartName() ;
    string filename = dirname + "/" + surfName+"_Surf.vtu" ;
    if(bit64)
      filename = dirname+ "/" + surfName+"_Surf64.vtu" ;
  
    long long unsigned int npnts, ncells;
    
    long long unsigned int elem_offset;
    
    vector<int> elem_ids ;
    vector<int> elem_conn ;
    vector<int> elem_offsets ;
    
    
    vector<unsigned char> elem_types ;
    
    vector<int> data_size ;
    vector<float> elem_data ;
    vector<string> data_names ;
    vector<float> point_data ;
    vector<string>  point_data_names ;
    vector<int> point_data_size ;
    vector<float>  position;
    int int_size;
    
    npnts = 0; ncells = 0; elem_offset = 0;
    int_size = bit64 ? sizeof(long long unsigned int) : sizeof(int);
  

    set<string> element_scalars ;
    set<string> element_vectors ;
    vector<string> escalars = surfacePartList[i]->getElementScalarVars() ;
    for(size_t j=0;j<escalars.size();++j) 
      element_scalars.insert(escalars[j]) ;
    vector<string> evectors = surfacePartList[i]->getElementVectorVars() ;
    for(size_t j=0;j<evectors.size();++j) 
      element_vectors.insert(evectors[j]) ;
    set<string> point_scalars ;
    set<string> point_vectors ;
    vector<string> pscalars = surfacePartList[i]->getNodalScalarVars() ;
    vector<string> pvectors = surfacePartList[i]->getNodalVectorVars() ;
    for(size_t j=0;j<pscalars.size();++j)
      point_scalars.insert(pscalars[j]) ;
    for(size_t j=0;j<pvectors.size();++j)
      point_vectors.insert(pvectors[j]) ;

  
    fid = fopen(filename.c_str(),"w");
    fprintf(fid,"<?xml version='1.0'?>\n");
    if (bit64) fprintf(fid,"<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>\n");
    else fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
    fprintf(fid,"  <UnstructuredGrid>\n");
  
    int tot = 0;
    size_t npt = surfacePartList[i]->getNumNodes() ;
    tot += npt;
    npnts = tot;
  
    //write boundary mesh

    int face_id = 0;
    //write out quads
    int nquads = surfacePartList[i]->getNumQuads() ;
    if(nquads > 0) {
      vector<Array<int,4> > quads ;
      surfacePartList[i]->getQuads(quads) ;
      
      int type = 9;
      unsigned char char_type = (unsigned char)type;
      for (int fi=0;fi<nquads;fi++) {
	elem_ids.push_back(face_id++);
	for (int j=0;j<4;j++) {
	  int id = quads[fi][j] -1;
	  elem_conn.push_back(id);
	}
	elem_offset += 4;
	elem_offsets.push_back(elem_offset);
	elem_types.push_back(char_type);
      }
    }
    //write out trias
    int ntrias = surfacePartList[i]->getNumTrias() ;
    if(ntrias > 0) {
      vector<Array<int,3> > trias ;
      surfacePartList[i]->getTrias(trias) ; 
      int type = 5;
      unsigned char char_type = (unsigned char)type;
      for (int fi=0;fi<ntrias;fi++) {
	elem_ids.push_back(face_id++);
	for (int j=0;j<3;j++) {
	  int id = trias[fi][j] -1;
	  elem_conn.push_back(id);
	}
	elem_offset += 3;
	elem_offsets.push_back(elem_offset);
	elem_types.push_back(char_type);
      }
    }
    //write out general faces
    int ngeneral = surfacePartList[i]->getNumGenfc() ;
    if(ngeneral > 0) {
      vector<int> nside_sizes,nside_nodes ;
      surfacePartList[i]->getGenf(nside_sizes,nside_nodes) ;
      {
	int tot = 0 ;
	for(int j=0;j<ngeneral;++j)
	  tot += nside_sizes[j] ;
	int nside_nodes_size = nside_nodes.size() ;
	if(nside_nodes_size != tot) {
	  cerr << "mismatch in node size and faces size " << nside_nodes_size
	       << " was " << tot << endl ;
	}
      }
      int type = 7,cnt=0;
      unsigned char char_type = (unsigned char)type;
      for (int fi=0;fi<ngeneral;fi++) {
	
	elem_ids.push_back(face_id++);
	for (int j=0;j<nside_sizes[fi];j++) {
	  int id = nside_nodes[cnt++]-1;
	  elem_conn.push_back(id);
	}
	elem_offset += nside_sizes[fi];
	elem_offsets.push_back(elem_offset);
	elem_types.push_back(char_type);
      }
    }

    set<string>::const_iterator si ;
    // write out point scalars
    for( si=point_scalars.begin();si!=point_scalars.end();++si){
      string varname = *si ;
      point_data_names.push_back(varname);
      int size = surfacePartList[i]->getNumNodes() ;
      point_data_size.push_back(size);    
      vector<float> val ;
      surfacePartList[i]->getNodalScalar(varname,val);
      for(int j=0;j<size;++j) 
	point_data.push_back(val[j]) ;
    }

    // write put point vectors
    for( si=point_vectors.begin();si!=point_vectors.end();++si){
      string varname = *si ;
      point_data_names.push_back(varname);
      int size = surfacePartList[i]->getNumNodes() ;
      point_data_size.push_back(size*3);    
      vector<vector3d<float> > val ;
      surfacePartList[i]->getNodalVector(varname,val) ;
      for(int j=0;j<size;++j) {
	point_data.push_back(val[j].x) ;
	point_data.push_back(val[j].y) ;
	point_data.push_back(val[j].z) ;
      }
    }

    // write out element scalars
    for( si=element_scalars.begin();si!=element_scalars.end();++si){
      string varname = *si ;
      data_names.push_back(varname);
      int size = elem_ids.size();
      data_size.push_back(size);    
      vector<float> valout(size, 0);  
      if(surfacePartList[i]->hasElementScalarVar(varname)) {
        vector<float> qvals, tvals, gvals ;
        surfacePartList[i]->getElementScalar(varname,qvals,tvals,gvals) ;
        int nqval = qvals.size();
        int ntval = tvals.size();
        int ngval = gvals.size();
        for (int fi=0;fi<nqval;fi++) valout[fi] = qvals[fi];
        for (int fi=0;fi<ntval;fi++) valout[fi+nqval] = tvals[fi];
        for (int fi=0;fi<ngval;fi++) valout[fi+nqval+ntval] = gvals[fi];
      }
      for (int i=0;i<size;i++) elem_data.push_back(valout[i]);
    }

    //output boundary vector
    for(si=element_vectors.begin();si!=element_vectors.end();++si) {
      string varname = *si ;
      data_names.push_back(varname);
      int size = 3 * elem_ids.size();
      data_size.push_back(size);
      
      vector<float> xvalout(elem_ids.size(),0);  
      vector<float> yvalout(elem_ids.size(),0);  
      vector<float> zvalout(elem_ids.size(),0);
      if(surfacePartList[i]->hasElementVectorVar(varname)) {
        vector<vector3d<float> > qvals, tvals, gvals ;
        surfacePartList[i]->getElementVector(varname,qvals,tvals,gvals) ;
        int nqval = qvals.size();
        int ntval = tvals.size();
        int ngval = gvals.size();
        
        for (int fi=0;fi<nqval;fi++) xvalout[fi] = qvals[fi].x;
        for (int fi=0;fi<nqval;fi++) yvalout[fi] = qvals[fi].y;
        for (int fi=0;fi<nqval;fi++) zvalout[fi] = qvals[fi].z;

        for (int fi=0;fi<ntval;fi++) xvalout[fi+nqval] = tvals[fi].x;
        for (int fi=0;fi<ntval;fi++) yvalout[fi+nqval] = tvals[fi].y;
        for (int fi=0;fi<ntval;fi++) zvalout[fi+nqval] = tvals[fi].z;

        for (int fi=0;fi<ngval;fi++) xvalout[fi+nqval+ntval] = gvals[fi].x;
        for (int fi=0;fi<ngval;fi++) yvalout[fi+nqval+ntval] = gvals[fi].y;
        for (int fi=0;fi<ngval;fi++) zvalout[fi+nqval+ntval] = gvals[fi].z;
      }
      for (int i=0;i<(int)elem_ids.size();i++) {
	elem_data.push_back(xvalout[i]);
	elem_data.push_back(yvalout[i]);
	elem_data.push_back(zvalout[i]);
      }
    }


    //create mesh positions
    {
      position.resize(3*npnts);
      vector<vector3d<float> > pos ;
      surfacePartList[i]->getPos(pos) ;
      int npt  = pos.size();
      for (int j=0;j<npt;j++) {
        position[3*j+0] = pos[j].x; 
        position[3*j+1] = pos[j].y; 
        position[3*j+2] = pos[j].z; 
      }
    }

    ncells = elem_types.size();
    
    long long unsigned int Offset = 0;
    fprintf(fid,"    <Piece NumberOfPoints='%llu' NumberOfCells='%llu'>\n",npnts,ncells);
    fprintf(fid,"      <Points>\n");
    fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%llu'/>\n",Offset);
    fprintf(fid,"      </Points>\n");
    fprintf(fid,"      <Cells>\n");
    Offset += 3 * npnts * sizeof(float) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += (long long unsigned int)elem_conn.size() * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += ncells * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += ncells * sizeof (unsigned char) + int_size ;
    fprintf(fid,"      </Cells>\n");


    fprintf(fid,"      <PointData>\n");
    for (long long unsigned int i=0;i<(long long unsigned int)point_data_names.size();i++) { 
      int comp=-1;
      if      ( (int) point_data_size[i] ==     (int) npnts) comp = 1;
      else if ( (int) point_data_size[i] == 3 * (int) npnts) comp = 3;
      else { cout << "Wrong size" << endl; exit(1); }
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",point_data_names[i].c_str(),comp,Offset) ;
      Offset += comp * npnts * sizeof (float) + int_size; 
    }
    fprintf(fid,"      </PointData>\n");


    fprintf(fid,"      <CellData>\n");
    
    for (long long unsigned int i=0;i<(long long unsigned int)data_names.size();i++) { 
      int comp=-1;
      if      ( (int) data_size[i] ==     (int) ncells) comp = 1;
      else if ( (int) data_size[i] == 3 * (int) ncells) comp = 3;
      else { cout << "Wrong size" << endl; exit(1); }
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",data_names[i].c_str(),comp,Offset) ;
      Offset += comp * ncells * sizeof (float) + int_size; 
    }
    fprintf(fid,"      </CellData>\n");
    fprintf(fid,"    </Piece>\n");
    fprintf(fid,"  </UnstructuredGrid>\n");
    fprintf(fid,"  <AppendedData encoding='raw'>\n");
    fprintf(fid,"_");

    long long unsigned int Scalar = (long long unsigned int) npnts * sizeof (float), Vector = 3 * Scalar;
    long long unsigned int Cells = (long long unsigned int) elem_types.size() * sizeof(int);
    long long unsigned int CellChars = (long long unsigned int) elem_types.size() * sizeof(unsigned char);
    long long unsigned int Conn = (long long unsigned int) elem_conn.size() * sizeof(int);
    fwrite((const char *) (&Vector), int_size, 1, fid) ;
    fwrite((const char *) (&position[0]), sizeof (float), 3 * (long long unsigned int)npnts, fid) ;
    fwrite((const char *) (&Conn), int_size, 1, fid) ;
    fwrite((const char *) (&elem_conn[0]), sizeof (int), (long long unsigned int) elem_conn.size(), fid) ;
    fwrite((const char *) (&Cells), int_size, 1, fid) ;
    fwrite((const char *) (&elem_offsets[0]), sizeof (int), (long long unsigned int) elem_offsets.size(), fid) ;
    fwrite((const char *) (&CellChars), int_size, 1, fid) ;
    fwrite((const char *) (&elem_types[0]), sizeof (unsigned char), (long long unsigned int) elem_types.size(), fid) ;

    long long unsigned int curr_loc = 0;
    for (size_t j=0;j<point_data_size.size();j++) {
      long long unsigned int ndata = point_data_size[j] ;
      long long unsigned int ndata_size = ndata * sizeof (float);
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      fwrite((const char *) (&point_data[curr_loc]), sizeof (float), ndata, fid) ;
      curr_loc += ndata;
    }

    curr_loc = 0;
    for (size_t j=0;j<data_size.size();j++) {
      long long unsigned int ndata = data_size[j] ;
      long long unsigned int ndata_size = ndata * sizeof (float);
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      fwrite((const char *) (&elem_data[curr_loc]), sizeof (float), ndata, fid) ;
      curr_loc += ndata;
    }
    fprintf(fid,"  </AppendedData>\n");
    fprintf(fid,"</VTKFile>\n");

    fclose(fid);
  }
  if(particlePartList.size() > 0) {
    int npnts = 0 ;
    for(size_t i=0;i<particlePartList.size();++i)
      npnts += particlePartList[i]->getNumParticles() ;
    unsigned long long int nP = npnts ;
    set<string> particle_scalars ;
    set<string> particle_vectors ;
    for(size_t i=0;i<particlePartList.size();++i) {
      vector<string> nscalars = particlePartList[i]->getScalarVars() ;
      for(size_t j=0;j<nscalars.size();++j) {
	particle_scalars.insert(nscalars[j]) ;
      }
      vector<string> nvectors = particlePartList[i]->getVectorVars() ;
      for(size_t j=0;j<nvectors.size();++j) {
	particle_vectors.insert(nvectors[j]) ;
      }    
    }


    string filename = dirname + "/vtk_particle_"+casename+"_"+iteration+".vtu" ;
    FILE *fid = fopen(filename.c_str(),"w");
    fprintf(fid,"<?xml version='1.0'?>\n");
    fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
    fprintf(fid,"  <UnstructuredGrid>\n");
    unsigned long long Offset = 0;
    fprintf(fid,"    <Piece NumberOfPoints='%llu' NumberOfCells='%llu'>\n",nP,nP);
    fprintf(fid,"      <Points>\n");
    fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%llu'/>\n",Offset);
    fprintf(fid,"      </Points>\n");
    fprintf(fid,"      <Cells>\n");
    int int_size = sizeof(unsigned int);
    long long unsigned int node_size = (long long unsigned int)nP;
    Offset += 3 * nP * sizeof(float) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += node_size * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += nP * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += nP * sizeof (unsigned char) + int_size ;
    fprintf(fid,"      </Cells>\n");
    fprintf(fid,"      <PointData>\n");
    set<string>::const_iterator si ;
    for(si=particle_scalars.begin();si!=particle_scalars.end();++si) {
      string name = *si ;
      int comp=1 ;
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",name.c_str(),comp,Offset) ;
      Offset += comp * nP * sizeof (float) + int_size; 
    }

    for(si=particle_vectors.begin();si!=particle_vectors.end();++si) {
      string name = *si ;
      int comp=3 ;
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",name.c_str(),comp,Offset) ;
      Offset += comp * nP * sizeof (float) + int_size; 
    }
    fprintf(fid,"      </PointData>\n");
    fprintf(fid,"    </Piece>\n");
    fprintf(fid,"  </UnstructuredGrid>\n");
    fprintf(fid,"  <AppendedData encoding='raw'>\n");
    fprintf(fid,"_");
    long long unsigned int Scalar = node_size * sizeof (float), Vector = 3 * Scalar;
    long long unsigned int Cells = node_size * sizeof(int);
    long long unsigned int CellChars = node_size * sizeof(unsigned char);
    long long unsigned int Conn = node_size * sizeof(int);
    fwrite((const char *) (&Vector), int_size, 1, fid) ;
    for(size_t j=0;j<particlePartList.size();++j) {
      vector<vector3d<float> > ppos ;
      particlePartList[j]->getParticlePositions(ppos) ;
      fwrite((const char *)(&ppos[0]),sizeof(float),3*ppos.size(),fid) ;
    }
    
    //    fwrite((const char *) (&position[0]), sizeof (float), 3*nP, fid) ;
    fwrite((const char *) (&Conn), int_size, 1, fid) ;
    for(unsigned long long int i=0;i<nP;++i) 
      fwrite((const char *) (&i), sizeof (int), 1, fid) ;
    
    fwrite((const char *) (&Cells), int_size, 1, fid) ;
    for(unsigned long long int i=0;i<nP;++i) {
      int off = i+1;
      fwrite((const char *) (&off), sizeof (int), 1, fid) ;
    }
    fwrite((const char *) (&CellChars), int_size, 1, fid) ;
    unsigned char type = 1;
    
    for(unsigned long long int i=0;i<nP;++i) {
      fwrite((const char *) (&type), sizeof (unsigned char), 1, fid) ;
    }
  
    for(si=particle_scalars.begin();si!=particle_scalars.end();++si) {
      string name = *si ;
      int comp=1 ;
      long long unsigned int ndata_size = nP*sizeof(float)*comp ;
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      for(size_t i=0;i<particlePartList.size();++i) {
	int np = particlePartList[i]->getNumParticles() ;
	vector<float> val(np,0) ;
	particlePartList[i]->getParticleScalar(name,val) ;
	fwrite(&val[0],sizeof(float),np*comp,fid) ;
      }
    }
    for(si=particle_vectors.begin();si!=particle_vectors.end();++si) {
      string name = *si ;
      int comp=3 ;
      long long unsigned int ndata_size = nP*sizeof(float)*comp ;
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      for(size_t i=0;i<particlePartList.size();++i) {
	int np = particlePartList[i]->getNumParticles() ;
	vector<vector3d<float> > val(np,vector3d<float>(0,0,0)) ;
	particlePartList[i]->getParticleVector(name,val) ;
	fwrite(&val[0],sizeof(float),np*comp,fid) ;
      }
    }

    fprintf(fid,"  </AppendedData>\n");
    fprintf(fid,"</VTKFile>\n");
    fclose(fid);
  }
}

bool vtkPartConverter::processesVolumeElements() const {
  return true ;
}
bool vtkPartConverter::processesSurfaceElements() const {
  return true ;
}
bool vtkPartConverter::processesParticleElements() const {
  return true ;
}

void vtkPartConverter::exportPostProcessorFiles(string casename, string iteration) const {

  string dirname = casename+"_vtk."+iteration ;
  struct stat statbuf ;
  if(stat(dirname.c_str(),&statbuf))
    mkdir(dirname.c_str(),0755) ;
  else
    if(!S_ISDIR(statbuf.st_mode)) {
      cerr << "file '" << dirname << "' should be a directory!, rename and start again."
           << endl ;
      exit(-1) ;
    }

  if(volumePartList.size() > 0) {
    //open
    long long unsigned int npnts=0 ;
    long long unsigned int ntets=0, nprsm=0, npyrm=0, nhexs=0, ngen=0 ;
    long long unsigned int face_off = 0 ;
    long long unsigned int Offset = 0 ;
    long long unsigned int off = 0 ;
    long long unsigned int data_count = 0 ;
    string filename ;
    int nvars = 0 ;
    string boundary_name ;

    vector<float> data_store ;
    vector<int> data_size;
    vector<int> conn ;
    vector<int> cell_offsets ;
    vector<int> cell_faces ;
    vector<int> face_offsets ;
    vector<unsigned char> cell_types ;
 
    int int_size =bit64 ? sizeof(long long unsigned int) : sizeof(int); 



    //open file and write out header
 
    if (bit64) filename = dirname + "/Volume64.vtu" ;
    else filename = dirname + "/Volume.vtu" ;
  
    set<string> nodal_scalars ;
    set<string> nodal_vectors ;
 
    for(size_t i=0;i<volumePartList.size();++i) {
     
      vector<string> nscalars = volumePartList[i]->getNodalScalarVars() ;
      for(size_t j=0;j<nscalars.size();++j) 
	nodal_scalars.insert(nscalars[j]) ;
      vector<string> nvectors = volumePartList[i]->getNodalVectorVars() ;
      for(size_t j=0;j<nvectors.size();++j) 
	nodal_vectors.insert(nvectors[j]) ;
    
      npnts += volumePartList[i]->getNumNodes();
      ntets += volumePartList[i]->getNumTets();
      nprsm += volumePartList[i]->getNumPrsm();
      nhexs += volumePartList[i]->getNumHexs();
      npyrm += volumePartList[i]->getNumPyrm();
      ngen += volumePartList[i]->getNumGenc();
    }    
 
    long long unsigned int ncell = ntets+nprsm+npyrm+nhexs+ngen;
   
    
    nvars = nodal_scalars.size() + 3*nodal_vectors.size();
    data_store.resize(nvars * npnts);
    FILE * fid = fopen(filename.c_str(), "w"); 
    fprintf(fid,"<?xml version='1.0'?>\n");
    if (bit64) fprintf(fid,"<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>\n");
    else fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
    fprintf(fid,"  <UnstructuredGrid>\n");
    fprintf(fid,"    <Piece NumberOfPoints='%llu' NumberOfCells='%llu'>\n",npnts,ncell);
    fprintf(fid,"      <Points>\n");
    fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%llu'/>\n",Offset);
    fprintf(fid,"      </Points>\n");
    Offset += 3 * npnts * sizeof(float) + int_size ;
    fclose(fid);
  
    //create mesh positions
    vector<float> pos(3*npnts);
    int nparts = volumePartList.size();
    vector<int> node_offset(nparts);
  
    {
      vector<vector3d<float> > position(npnts);
      size_t tot = 0;
      for(size_t i=0;i<volumePartList.size();++i) {
	vector<vector3d<float> > part_pos;
	volumePartList[i]->getPos(part_pos) ;
	for(size_t j = 0; j<part_pos.size(); j++){
	  position[j+tot] = part_pos[j];
	}
	node_offset[i] = tot;
	tot += part_pos.size();
      }
    
      for(long long unsigned int i=0;i<npnts;++i) {
	long long unsigned int j=3*i;
	pos[j]   = position[i].x ;
	pos[j+1] = position[i].y ;
	pos[j+2] = position[i].z ;
      }
    }

    //write out tets
    const int block_size=65536 ; // Size of blocking factor
    for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumTets() > 0) { 
	int tot = volumePartList[i]->getNumTets() ;
	int start = 0 ;
	int top = tot + volumePartList[i]->getNumTetsIblank() ;
	while(start < top) {
	  vector<Array<int,4> > cells ;
	  volumePartList[i]->getTetBlock(cells,start,block_size) ;
	  size_t ncells = cells.size();
	  start += block_size ;
	  if(cells.size() > 0) {
	    for(size_t j=0;j<cells.size();++j) {
	      for(int k = 0; k < 4; k++){
		cells[j][k] += node_offset[i];
	      }
	    }
            
         
	    unsigned int type = 10 ;
	    unsigned char ctype = (unsigned char)type ;
	    int nnpc = 4;

	    const int nf=4,nnpf[nf]={3,3,3,3};
	    vector< vector<int> > ff(nf);
	    for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
	    ff[0][0] = 0; ff[0][1] = 1; ff[0][2] = 3;
	    ff[1][0] = 1; ff[1][1] = 2; ff[1][2] = 3;
	    ff[2][0] = 2; ff[2][1] = 0; ff[2][2] = 3;
	    ff[3][0] = 0; ff[3][1] = 2; ff[3][2] = 1;
          
	    if (ngen) {
	      for(size_t i=0;i<ncells;++i) {
		for (int j=0;j<nnpc;j++) { 
		  int nd = cells[i][j]-1 ;
		  conn.push_back(nd) ;
		  off++ ;
		}
		cell_types.push_back(ctype);
		cell_offsets.push_back(off);
              
		cell_faces.push_back(nf); face_off++;
		for (int j=0;j<nf;j++) {
		  cell_faces.push_back(nnpf[j]); face_off++;
		  for (int k=0;k<nnpf[j];k++) {
		    int nn = cells[i][ff[j][k]] - 1;
		    cell_faces.push_back(nn); face_off++;
		  }
		}
		face_offsets.push_back(face_off) ;
	      }
	    }
	    else {
	      for(size_t i=0;i<ncells;++i) {
		for (int j=0;j<nnpc;j++) { 
		  int nd = cells[i][j]-1 ;
		  conn.push_back(nd) ;
		  off++ ;
		}
		cell_types.push_back(ctype);
		cell_offsets.push_back(off);
	      }
	    }
	  }
	}
      }
    }

    // write out pyrms
  
    for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumPyrm() > 0) { 
	int tot = volumePartList[i]->getNumPyrm() ;
	int start = 0 ;
	int top = tot + volumePartList[i]->getNumPyrmIblank() ;
	while(start < top) {
	  vector<Array<int,5> > cells ;
	  volumePartList[i]->getPyrmBlock(cells,start,block_size) ;
	  size_t ncells = cells.size();
	  start += block_size ;
            
	  for(size_t j=0;j<cells.size();++j) {
	    for(int k = 0; k < 5; k++){
	      cells[j][k] += node_offset[i];
	    }
	  }
            
	  unsigned int type = 14 ;
	  unsigned char ctype = (unsigned char)type ;
	  int nnpc = 5;

	  const int nf=5,nnpf[nf]={4,3,3,3,3};
	  vector< vector<int> > ff(nf);
	  for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
	  ff[0][0] = 0; ff[0][1] = 3; ff[0][2] = 2; ff[0][3] = 1;
	  ff[1][0] = 0; ff[1][1] = 1; ff[1][2] = 4;
	  ff[2][0] = 1; ff[2][1] = 2; ff[2][2] = 4;
	  ff[3][0] = 2; ff[3][1] = 3; ff[3][2] = 4;
	  ff[4][0] = 3; ff[4][1] = 0; ff[4][2] = 4; 
        
	  if (ngen) {
	    for(size_t i=0;i<ncells;++i) {
	      for (int j=0;j<nnpc;j++) { 
		int nd = cells[i][j]-1 ;
		conn.push_back(nd) ;
		off++ ;
	      }
	      cell_types.push_back(ctype);
	      cell_offsets.push_back(off);
      
	      cell_faces.push_back(nf); face_off++;
	      for (int j=0;j<nf;j++) {
		cell_faces.push_back(nnpf[j]); face_off++;
		for (int k=0;k<nnpf[j];k++) {
		  int nn = cells[i][ff[j][k]] - 1;
		  cell_faces.push_back(nn); face_off++;
		}
	      }
	      face_offsets.push_back(face_off) ;
	    }
	  }
	  else {
	    for(size_t i=0;i<ncells;++i) {
	      for (int j=0;j<nnpc;j++) { 
		int nd = cells[i][j]-1 ;
		conn.push_back(nd) ;
		off++ ;
	      }
	      cell_types.push_back(ctype);
	      cell_offsets.push_back(off);
	    }
	  }
	}
      }
    }
  

    // write out prsms
  
    for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumPrsm() > 0) { 
	int tot = volumePartList[i]->getNumPrsm() ;
	int start = 0 ;
	int top = tot + volumePartList[i]->getNumPrsmIblank() ;
	while(start < top) {
	  vector<Array<int,6> > cells ;
	  volumePartList[i]->getPrsmBlock(cells,start,block_size) ;
	  size_t ncells = cells.size();
	  start += block_size ;
	  for(size_t j=0;j<cells.size();++j) {
	    for(int k = 0; k < 6; k++){
	      cells[j][k] += node_offset[i];
	    }
	  }
            
            
	  unsigned int type = 13 ;
	  unsigned char ctype = (unsigned char)type ;
	  int nnpc = 6;

	  const int nf=5,nnpf[nf]={4,4,4,3,3};
	  vector< vector<int> > ff(nf);
	  for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
	  ff[0][0] = 0; ff[0][1] = 2; ff[0][2] = 5; ff[0][3] = 3;
	  ff[1][0] = 2; ff[1][1] = 1; ff[1][2] = 4; ff[1][3] = 5;
	  ff[2][0] = 0; ff[2][1] = 3; ff[2][2] = 4; ff[2][3] = 1;
	  ff[3][0] = 0; ff[3][1] = 1; ff[3][2] = 2;
	  ff[4][0] = 3; ff[4][1] = 5; ff[4][2] = 4; 

	  if (ngen) {
	    for(size_t i=0;i<ncells;++i) {
	      for (int j=0;j<nnpc;j++) { 
		int nd = cells[i][j]-1 ;
		conn.push_back(nd) ;
		off++ ;
	      }
	      cell_types.push_back(ctype);
	      cell_offsets.push_back(off);
      
	      cell_faces.push_back(nf); face_off++;
	      for (int j=0;j<nf;j++) {
		cell_faces.push_back(nnpf[j]); face_off++;
		for (int k=0;k<nnpf[j];k++) {
		  int nn = cells[i][ff[j][k]] - 1;
		  cell_faces.push_back(nn); face_off++;
		}
	      }
	      face_offsets.push_back(face_off) ;
	    }
	  }
	  else {
	    for(size_t i=0;i<ncells;++i) {
	      for (int j=0;j<nnpc;j++) { 
		int nd = cells[i][j]-1 ;
		conn.push_back(nd) ;
		off++ ;
	      }
	      cell_types.push_back(ctype);
	      cell_offsets.push_back(off);
	    }
	  }
	}
      }
    }
  

    // write out hexs
  
    for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumHexs() > 0) { 
	int tot = volumePartList[i]->getNumHexs() ;
	int start = 0 ;
	int top = tot + volumePartList[i]->getNumHexsIblank() ;
	while(start < top) {
	  vector<Array<int,8> > cells ;
	  volumePartList[i]->getHexBlock(cells,start,block_size) ;
	  size_t ncells = cells.size();
	  start += block_size ;
	  for(size_t j=0;j<cells.size();++j) {
	    for(int k = 0; k < 8; k++){
	      cells[j][k] += node_offset[i];
	    }
	  }
            
	  unsigned int type = 12 ;
	  unsigned char ctype = (unsigned char)type ;
	  int nnpc = 8;
        
	  const int nf=6,nnpf[nf]={4,4,4,4,4,4};
	  vector< vector<int> > ff(nf);
	  for (int i=0;i<nf;i++) ff[i].resize(nnpf[i]);
	  ff[0][0] = 0; ff[0][1] = 1; ff[0][2] = 5; ff[0][3] = 4;
	  ff[1][0] = 1; ff[1][1] = 2; ff[1][2] = 6; ff[1][3] = 5;
	  ff[2][0] = 2; ff[2][1] = 3; ff[2][2] = 7; ff[2][3] = 6;
	  ff[3][0] = 3; ff[3][1] = 0; ff[3][2] = 4; ff[3][3] = 7;
	  ff[4][0] = 4; ff[4][1] = 5; ff[4][2] = 6; ff[4][3] = 7;
	  ff[5][0] = 0; ff[5][1] = 3; ff[5][2] = 2; ff[5][3] = 1;

	  if (ngen) {
	    for(size_t i=0;i<ncells;++i) {
	      for (int j=0;j<nnpc;j++) { 
		int nd = cells[i][j]-1 ;
		conn.push_back(nd) ;
		off++ ;
	      }
	      cell_types.push_back(ctype);
	      cell_offsets.push_back(off);
      
	      cell_faces.push_back(nf); face_off++;
	      for (int j=0;j<nf;j++) {
		cell_faces.push_back(nnpf[j]); face_off++;
		for (int k=0;k<nnpf[j];k++) {
		  int nn = cells[i][ff[j][k]] - 1;
		  cell_faces.push_back(nn); face_off++;
		}
	      }
	      face_offsets.push_back(face_off) ;
	    }
	  }
	  else {
	    for(size_t i=0;i<ncells;++i) {
	      for (int j=0;j<nnpc;j++) { 
		int nd = cells[i][j]-1 ;
		conn.push_back(nd) ;
		off++ ;
	      }
	      cell_types.push_back(ctype);
	      cell_offsets.push_back(off);
	    }
	  }
	}
      }
    }
  

    //write out genc
  
    for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumGenc() > 0) {
	vector<int> nfaces;
	vector<int> nsides;
	vector<int> nodes;
	volumePartList[i]->getGenCell(nfaces,nsides,nodes) ;
	size_t nnfaces = nfaces.size();
	size_t nnsides = nsides.size();
	size_t nnodes = nodes.size();
          
	for(size_t j=0;j<nodes.size();++j) {
	  nodes[j] += node_offset[i];
	}
           
            
	unsigned int type = 42 ;
	unsigned char ctype = (unsigned char)type ;
	int face = 0, node = 0 ;
	vector<int> tmp;
        
	for(size_t i=0;i<nnfaces;++i) {
	  int tsz=0,nf = nfaces[i] ;
	  for (int j=0;j<nf;j++) tsz += nsides[face++] ; 
	  tmp.resize(tsz) ; face-=nf ; tsz=0 ;
	  cell_faces.push_back(nf) ; face_off++ ;
	  for (int j=0;j<nf;j++) { 
	    int ne = nsides[face++] ; 
	    cell_faces.push_back(ne) ; face_off++ ;
	    for (int k=0;k<ne;k++) { 
	      int nd = nodes[node++]-1 ;
	      cell_faces.push_back(nd) ; face_off++ ;
	      tmp[tsz++] = nd ;
	    }
	  }
          
	  sort(tmp.begin(),tmp.end()) ;
	  tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end()) ;
          
	  for (size_t j=0;j<tmp.size();j++) conn.push_back(tmp[j]) ;  
	  off += (long long unsigned int) tmp.size() ;
          
	  cell_types.push_back(ctype) ;
	  cell_offsets.push_back(off) ;
	  face_offsets.push_back(face_off) ;
	}
        
	if (nnodes  != size_t(node)) cerr << "Problem in write_general_cell" << endl ;
	if (nnsides != size_t(face)) cerr << "Problem in write_general_cell" << endl ;
      }
    }
    
  
    //close mesh elements
    {
      unsigned int ncells = ntets+nprsm+npyrm+nhexs+ngen;
      FILE * fid = fopen(filename.c_str(),"a");
      fprintf(fid,"      <Cells>\n");
      fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
      Offset += off * sizeof (int) + int_size ;
      fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
      Offset += ncells * sizeof (int) + int_size ;
      fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
      Offset += ncells * sizeof (unsigned char) + int_size ;
      if (ngen) { // include faces and faceoffsets
	fprintf(fid,"        <DataArray type='Int32' Name='faces' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
	Offset += (long long unsigned int) cell_faces.size() * sizeof (int) + int_size ;
	fprintf(fid,"        <DataArray type='Int32' Name='faceoffsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
	Offset += (long long unsigned int) face_offsets.size() * sizeof (int) + int_size ;
      }
      fprintf(fid,"      </Cells>\n");
      fprintf(fid,"      <PointData>\n");
      fclose(fid);
    }
  
    set<string>::const_iterator si ;
    // write out nodal scalars
    for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {                                         
      string varname = *si ;
      vector<float> val_out(npnts, 0.0);
      for(size_t i =0;i<volumePartList.size();++i) {
	if(volumePartList[i]->hasNodalScalarVar(varname)) {
	  vector<float> val ;
	  volumePartList[i]->getNodalScalar(varname,val);
	  for(long long unsigned int j = 0; j < (long long unsigned int)val.size(); j++) val_out[j+node_offset[i]] = val[j] ;
	}
      }
      for(long long unsigned int j = 0; j < (long long unsigned int)npnts; j++) data_store[data_count++] = val_out[j] ; 
    
      FILE * fid = fopen(filename.c_str(),"a") ;
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='1' format='appended' offset='%llu'/>\n",varname.c_str(),Offset) ;
      Offset += (long long unsigned int) npnts * sizeof(float) + int_size ;
      long long unsigned int ndata = (long long unsigned int) npnts;
      data_size.push_back(ndata);
      fclose(fid);
    }

    //output nodal vector
    for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
      string varname = *si ;
      vector<vector3d<float> > val_out(npnts, vector3d<float>(0.0, 0.0, 0.0));
      for(size_t i =0;i<volumePartList.size();++i) {
	if(volumePartList[i]->hasNodalVectorVar(varname)) {
	  vector<vector3d<float> > val ;
	  volumePartList[i]->getNodalVector(varname,val) ;
	  for(long long unsigned int j = 0; j < (long long unsigned int)val.size(); j++) val_out[j+node_offset[i]] = val[j] ;
         
	}
      }
      for(long long unsigned int j = 0; j < (long long unsigned int) npnts; j++) { 
	data_store[data_count++] = val_out[j].x ;
	data_store[data_count++] = val_out[j].y ;
	data_store[data_count++] = val_out[j].z ;
      }
      FILE * fid = fopen(filename.c_str(),"a");
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='3' format='appended' offset='%llu'/>\n",varname.c_str(),Offset);
      Offset += 3 * (long long unsigned) npnts * sizeof(float) + int_size ;
      long long unsigned int vec_size = 3 * (long long unsigned int)npnts;
      data_size.push_back(vec_size);
      fclose(fid);
    }


    //close
    {
 
      FILE * fid = fopen(filename.c_str(),"a");
      fprintf(fid,"      </PointData>\n");
      fprintf(fid,"    </Piece>\n");
      fprintf(fid,"  </UnstructuredGrid>\n");
      fprintf(fid,"  <AppendedData encoding='raw'>\n");
      fprintf(fid,"_");
      long long unsigned int Scalar = npnts * sizeof (float), Vector = 3 * Scalar;
      long long unsigned int Cells = (long long unsigned int) cell_types.size() * sizeof(int);
      long long unsigned int CellChars = (long long unsigned int) cell_types.size() * sizeof(unsigned char);
      long long unsigned int Conn = (long long unsigned int) conn.size() * sizeof(int);
      if (!bit64) {
	if ( Scalar > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
	if ( Vector > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
	if ( Cells > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
	if ( Conn > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
      }
      fwrite((const char *) (&Vector), int_size, 1, fid) ;
      fwrite((const char *) (&pos[0]), sizeof (float), 3 * npnts, fid) ;

      fwrite((const char *) (&Conn), int_size, 1, fid) ;
      fwrite((const char *) &conn[0], sizeof (int), (long long unsigned int) conn.size(), fid) ;
      vector<int>().swap(conn) ;
      fwrite((const char *) (&Cells), int_size, 1, fid) ;
      fwrite((const char *) &cell_offsets[0], sizeof (int), (long long unsigned int) cell_offsets.size(), fid) ;
      vector<int>().swap(cell_offsets) ;
      fwrite((const char *) (&CellChars), int_size, 1, fid) ;
      fwrite((const char *) &cell_types[0], sizeof (unsigned char), (long long unsigned int) cell_types.size(), fid) ;
      vector<unsigned char>().swap(cell_types) ;
      if (ngen) {
	long long unsigned int Faces = (long long unsigned int) cell_faces.size() * sizeof(int) ;
	long long unsigned int FaceConn = (long long unsigned int) face_offsets.size() * sizeof(int) ;
	if (!bit64) {
	  if ( Faces > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
	  if ( FaceConn > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
	}
	fwrite((const char *) (&Faces), int_size, 1, fid) ;
	fwrite((const char *) &cell_faces[0], sizeof (int), (long long unsigned int) cell_faces.size(), fid) ;
	fwrite((const char *) (&FaceConn), int_size, 1, fid) ;
	fwrite((const char *) &face_offsets[0], sizeof (int), (long long unsigned int) face_offsets.size(), fid) ;
      }
      long long unsigned int curr_loc = 0;
      for (long long unsigned int i=0;i<(long long unsigned int)data_size.size();i++) {
	long long unsigned int ndata = data_size[i], ndata_size = ndata * sizeof (float); 
	if (!bit64) {
	  if ( ndata_size > std::numeric_limits<unsigned int>::max()) { cerr << "Dataset too large. Must use -vtk64 and Paraview 3.98." << endl; exit(1); }
	}
	fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
	fwrite((const char *) (&data_store[curr_loc]), sizeof (float), ndata, fid) ;
	curr_loc += ndata;
      }
      vector<int>().swap(data_size);
      fprintf(fid,"  </AppendedData>\n");
      fprintf(fid,"</VTKFile>\n");
      fclose(fid);
    }
  }

  int nsurfaces = surfacePartList.size() ;
  for(int i=0;i<nsurfaces;++i) {
    string surfName = surfacePartList[i]->getPartName() ;
    string filename = dirname + "/" + surfName+"_Surf.vtu" ;
    if(bit64)
      filename = dirname+ "/" + surfName+"_Surf64.vtu" ;
  
    long long unsigned int npnts, ncells;
    
    long long unsigned int elem_offset;
    
    vector<int> elem_ids ;
    vector<int> elem_conn ;
    vector<int> elem_offsets ;
    
    
    vector<unsigned char> elem_types ;
    
    vector<int> data_size ;
    vector<float> elem_data ;
    vector<string> data_names ;
    vector<float> point_data ;
    vector<string>  point_data_names ;
    vector<int> point_data_size ;
    vector<float>  position;
    int int_size;
    
    npnts = 0; ncells = 0; elem_offset = 0;
    int_size = bit64 ? sizeof(long long unsigned int) : sizeof(int);
  

    set<string> element_scalars ;
    set<string> element_vectors ;
    vector<string> escalars = surfacePartList[i]->getElementScalarVars() ;
    for(size_t j=0;j<escalars.size();++j) 
      element_scalars.insert(escalars[j]) ;
    vector<string> evectors = surfacePartList[i]->getElementVectorVars() ;
    for(size_t j=0;j<evectors.size();++j) 
      element_vectors.insert(evectors[j]) ;
    set<string> point_scalars ;
    set<string> point_vectors ;
    vector<string> pscalars = surfacePartList[i]->getNodalScalarVars() ;
    vector<string> pvectors = surfacePartList[i]->getNodalVectorVars() ;
    for(size_t j=0;j<pscalars.size();++j)
      point_scalars.insert(pscalars[j]) ;
    for(size_t j=0;j<pvectors.size();++j)
      point_vectors.insert(pvectors[j]) ;

  
    FILE *fid = fopen(filename.c_str(),"w");
    fprintf(fid,"<?xml version='1.0'?>\n");
    if (bit64) fprintf(fid,"<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>\n");
    else fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
    fprintf(fid,"  <UnstructuredGrid>\n");
  
    int tot = 0;
    size_t npt = surfacePartList[i]->getNumNodes() ;
    tot += npt;
    npnts = tot;
  
    //write boundary mesh

    int face_id = 0;
    //write out quads
    int nquads = surfacePartList[i]->getNumQuads() ;
    if(nquads > 0) {
      vector<Array<int,4> > quads ;
      surfacePartList[i]->getQuads(quads) ;
      
      int type = 9;
      unsigned char char_type = (unsigned char)type;
      for (int fi=0;fi<nquads;fi++) {
	elem_ids.push_back(face_id++);
	for (int j=0;j<4;j++) {
	  int id = quads[fi][j] -1;
	  elem_conn.push_back(id);
	}
	elem_offset += 4;
	elem_offsets.push_back(elem_offset);
	elem_types.push_back(char_type);
      }
    }
    //write out trias
    int ntrias = surfacePartList[i]->getNumTrias() ;
    if(ntrias > 0) {
      vector<Array<int,3> > trias ;
      surfacePartList[i]->getTrias(trias) ; 
      int type = 5;
      unsigned char char_type = (unsigned char)type;
      for (int fi=0;fi<ntrias;fi++) {
	elem_ids.push_back(face_id++);
	for (int j=0;j<3;j++) {
	  int id = trias[fi][j] -1;
	  elem_conn.push_back(id);
	}
	elem_offset += 3;
	elem_offsets.push_back(elem_offset);
	elem_types.push_back(char_type);
      }
    }
    //write out general faces
    int ngeneral = surfacePartList[i]->getNumGenfc() ;
    if(ngeneral > 0) {
      vector<int> nside_sizes,nside_nodes ;
      surfacePartList[i]->getGenf(nside_sizes,nside_nodes) ;
      {
	int tot = 0 ;
	for(int j=0;j<ngeneral;++j)
	  tot += nside_sizes[j] ;
	int nside_nodes_size = nside_nodes.size() ;
	if(nside_nodes_size != tot) {
	  cerr << "mismatch in node size and faces size " << nside_nodes_size
	       << " was " << tot << endl ;
	}
      }
      int type = 7,cnt=0;
      unsigned char char_type = (unsigned char)type;
      for (int fi=0;fi<ngeneral;fi++) {
	
	elem_ids.push_back(face_id++);
	for (int j=0;j<nside_sizes[fi];j++) {
	  int id = nside_nodes[cnt++]-1;
	  elem_conn.push_back(id);
	}
	elem_offset += nside_sizes[fi];
	elem_offsets.push_back(elem_offset);
	elem_types.push_back(char_type);
      }
    }

    set<string>::const_iterator si ;
    // write out point scalars
    for( si=point_scalars.begin();si!=point_scalars.end();++si){
      string varname = *si ;
      point_data_names.push_back(varname);
      int size = surfacePartList[i]->getNumNodes() ;
      point_data_size.push_back(size);    
      vector<float> val ;
      surfacePartList[i]->getNodalScalar(varname,val);
      for(int j=0;j<size;++j) 
	point_data.push_back(val[j]) ;
    }

    // write put point vectors
    for( si=point_vectors.begin();si!=point_vectors.end();++si){
      string varname = *si ;
      point_data_names.push_back(varname);
      int size = surfacePartList[i]->getNumNodes() ;
      point_data_size.push_back(size*3);    
      vector<vector3d<float> > val ;
      surfacePartList[i]->getNodalVector(varname,val) ;
      for(int j=0;j<size;++j) {
	point_data.push_back(val[j].x) ;
	point_data.push_back(val[j].y) ;
	point_data.push_back(val[j].z) ;
      }
    }

    // write out element scalars
    for( si=element_scalars.begin();si!=element_scalars.end();++si){
      string varname = *si ;
      data_names.push_back(varname);
      int size = elem_ids.size();
      data_size.push_back(size);    
      vector<float> valout(size, 0);  
      if(surfacePartList[i]->hasElementScalarVar(varname)) {
        vector<float> qvals, tvals, gvals ;
        surfacePartList[i]->getElementScalar(varname,qvals,tvals,gvals) ;
        int nqval = qvals.size();
        int ntval = tvals.size();
        int ngval = gvals.size();
        for (int fi=0;fi<nqval;fi++) valout[fi] = qvals[fi];
        for (int fi=0;fi<ntval;fi++) valout[fi+nqval] = tvals[fi];
        for (int fi=0;fi<ngval;fi++) valout[fi+nqval+ntval] = gvals[fi];
      }
      for (int i=0;i<size;i++) elem_data.push_back(valout[i]);
    }

    //output boundary vector
    for(si=element_vectors.begin();si!=element_vectors.end();++si) {
      string varname = *si ;
      data_names.push_back(varname);
      int size = 3 * elem_ids.size();
      data_size.push_back(size);
      
      vector<float> xvalout(elem_ids.size(),0);  
      vector<float> yvalout(elem_ids.size(),0);  
      vector<float> zvalout(elem_ids.size(),0);
      if(surfacePartList[i]->hasElementVectorVar(varname)) {
        vector<vector3d<float> > qvals, tvals, gvals ;
        surfacePartList[i]->getElementVector(varname,qvals,tvals,gvals) ;
        int nqval = qvals.size();
        int ntval = tvals.size();
        int ngval = gvals.size();
        
        for (int fi=0;fi<nqval;fi++) xvalout[fi] = qvals[fi].x;
        for (int fi=0;fi<nqval;fi++) yvalout[fi] = qvals[fi].y;
        for (int fi=0;fi<nqval;fi++) zvalout[fi] = qvals[fi].z;

        for (int fi=0;fi<ntval;fi++) xvalout[fi+nqval] = tvals[fi].x;
        for (int fi=0;fi<ntval;fi++) yvalout[fi+nqval] = tvals[fi].y;
        for (int fi=0;fi<ntval;fi++) zvalout[fi+nqval] = tvals[fi].z;

        for (int fi=0;fi<ngval;fi++) xvalout[fi+nqval+ntval] = gvals[fi].x;
        for (int fi=0;fi<ngval;fi++) yvalout[fi+nqval+ntval] = gvals[fi].y;
        for (int fi=0;fi<ngval;fi++) zvalout[fi+nqval+ntval] = gvals[fi].z;
      }
      for (int i=0;i<(int)elem_ids.size();i++) {
	elem_data.push_back(xvalout[i]);
	elem_data.push_back(yvalout[i]);
	elem_data.push_back(zvalout[i]);
      }
    }


    //create mesh positions
    {
      position.resize(3*npnts);
      vector<vector3d<float> > pos ;
      surfacePartList[i]->getPos(pos) ;
      int npt  = pos.size();
      for (int j=0;j<npt;j++) {
        position[3*j+0]   = pos[j].x; 
        position[3*j+1] = pos[j].y; 
        position[3*j+2] = pos[j].z; 
      }
    }

    ncells = elem_types.size();
    
    long long unsigned int Offset = 0;
    fprintf(fid,"    <Piece NumberOfPoints='%llu' NumberOfCells='%llu'>\n",npnts,ncells);
    fprintf(fid,"      <Points>\n");
    fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%llu'/>\n",Offset);
    fprintf(fid,"      </Points>\n");
    fprintf(fid,"      <Cells>\n");
    Offset += 3 * npnts * sizeof(float) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += (long long unsigned int)elem_conn.size() * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += ncells * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += ncells * sizeof (unsigned char) + int_size ;
    fprintf(fid,"      </Cells>\n");


    fprintf(fid,"      <PointData>\n");
    for (long long unsigned int i=0;i<(long long unsigned int)point_data_names.size();i++) { 
      int comp=-1;
      if      ( (int) point_data_size[i] ==     (int) npnts) comp = 1;
      else if ( (int) point_data_size[i] == 3 * (int) npnts) comp = 3;
      else { cout << "Wrong size" << endl; exit(1); }
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",point_data_names[i].c_str(),comp,Offset) ;
      Offset += comp * npnts * sizeof (float) + int_size; 
    }
    fprintf(fid,"      </PointData>\n");


    fprintf(fid,"      <CellData>\n");
    for (long long unsigned int i=0;i<(long long unsigned int)data_names.size();i++) { 
      int comp=-1;
      if      ( (int) data_size[i] ==     (int) ncells) comp = 1;
      else if ( (int) data_size[i] == 3 * (int) ncells) comp = 3;
      else { cout << "Wrong size" << endl; exit(1); }
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",data_names[i].c_str(),comp,Offset) ;
      Offset += comp * ncells * sizeof (float) + int_size; 
    }
    fprintf(fid,"      </CellData>\n");
    fprintf(fid,"    </Piece>\n");
    fprintf(fid,"  </UnstructuredGrid>\n");
    fprintf(fid,"  <AppendedData encoding='raw'>\n");
    fprintf(fid,"_");

    long long unsigned int Scalar = (long long unsigned int) npnts * sizeof (float), Vector = 3 * Scalar;
    long long unsigned int Cells = (long long unsigned int) elem_types.size() * sizeof(int);
    long long unsigned int CellChars = (long long unsigned int) elem_types.size() * sizeof(unsigned char);
    long long unsigned int Conn = (long long unsigned int) elem_conn.size() * sizeof(int);
    fwrite((const char *) (&Vector), int_size, 1, fid) ;
    fwrite((const char *) (&position[0]), sizeof (float), 3 * (long long unsigned int)npnts, fid) ;
    fwrite((const char *) (&Conn), int_size, 1, fid) ;
    fwrite((const char *) (&elem_conn[0]), sizeof (int), (long long unsigned int) elem_conn.size(), fid) ;
    fwrite((const char *) (&Cells), int_size, 1, fid) ;
    fwrite((const char *) (&elem_offsets[0]), sizeof (int), (long long unsigned int) elem_offsets.size(), fid) ;
    fwrite((const char *) (&CellChars), int_size, 1, fid) ;
    fwrite((const char *) (&elem_types[0]), sizeof (unsigned char), (long long unsigned int) elem_types.size(), fid) ;

    long long unsigned int curr_loc = 0;
    for (size_t j=0;j<point_data_size.size();j++) {
      long long unsigned int ndata = point_data_size[j] ;
      long long unsigned int ndata_size = ndata * sizeof (float);
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      fwrite((const char *) (&point_data[curr_loc]), sizeof (float), ndata, fid) ;
      curr_loc += ndata;
    }

    curr_loc = 0;
    for (size_t j=0;j<data_size.size();j++) {
      long long unsigned int ndata = data_size[j] ;
      long long unsigned int ndata_size = ndata * sizeof (float);
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      fwrite((const char *) (&elem_data[curr_loc]), sizeof (float), ndata, fid) ;
      curr_loc += ndata;
    }
    fprintf(fid,"  </AppendedData>\n");
    fprintf(fid,"</VTKFile>\n");

    fclose(fid);
  }
  if(particlePartList.size() > 0) {
    int npnts = 0 ;
    for(size_t i=0;i<particlePartList.size();++i)
      npnts += particlePartList[i]->getNumParticles() ;
    unsigned long long int nP = npnts ;
    set<string> particle_scalars ;
    set<string> particle_vectors ;
    for(size_t i=0;i<particlePartList.size();++i) {
      vector<string> nscalars = particlePartList[i]->getScalarVars() ;
      for(size_t j=0;j<nscalars.size();++j) {
	particle_scalars.insert(nscalars[j]) ;
      }
      vector<string> nvectors = particlePartList[i]->getVectorVars() ;
      for(size_t j=0;j<nvectors.size();++j) {
	particle_vectors.insert(nvectors[j]) ;
      }    
    }


    string filename = dirname+"/vtk_particle_"+casename+"_"+iteration+".vtu" ;
    FILE *fid = fopen(filename.c_str(),"w");
    fprintf(fid,"<?xml version='1.0'?>\n");
    fprintf(fid,"<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>\n");
    fprintf(fid,"  <UnstructuredGrid>\n");
    unsigned long long Offset = 0;
    fprintf(fid,"    <Piece NumberOfPoints='%llu' NumberOfCells='%llu'>\n",nP,nP);
    fprintf(fid,"      <Points>\n");
    fprintf(fid,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='%llu'/>\n",Offset);
    fprintf(fid,"      </Points>\n");
    fprintf(fid,"      <Cells>\n");
    int int_size = sizeof(unsigned int);
    long long unsigned int node_size = (long long unsigned int)nP;
    Offset += 3 * nP * sizeof(float) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += node_size * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += nP * sizeof (int) + int_size ;
    fprintf(fid,"        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='%llu'/>\n",Offset);
    Offset += nP * sizeof (unsigned char) + int_size ;
    fprintf(fid,"      </Cells>\n");
    fprintf(fid,"      <PointData>\n");
    set<string>::const_iterator si ;
    for(si=particle_scalars.begin();si!=particle_scalars.end();++si) {
      string name = *si ;
      int comp=1 ;
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",name.c_str(),comp,Offset) ;
      Offset += comp * nP * sizeof (float) + int_size; 
    }

    for(si=particle_vectors.begin();si!=particle_vectors.end();++si) {
      string name = *si ;
      int comp=3 ;
      fprintf(fid,"        <DataArray type='Float32' Name='%s' NumberOfComponents='%d' format='appended' offset='%llu'/>\n",name.c_str(),comp,Offset) ;
      Offset += comp * nP * sizeof (float) + int_size; 
    }
    fprintf(fid,"      </PointData>\n");
    fprintf(fid,"    </Piece>\n");
    fprintf(fid,"  </UnstructuredGrid>\n");
    fprintf(fid,"  <AppendedData encoding='raw'>\n");
    fprintf(fid,"_");
    long long unsigned int Scalar = node_size * sizeof (float), Vector = 3 * Scalar;
    long long unsigned int Cells = node_size * sizeof(int);
    long long unsigned int CellChars = node_size * sizeof(unsigned char);
    long long unsigned int Conn = node_size * sizeof(int);
    fwrite((const char *) (&Vector), int_size, 1, fid) ;
    for(size_t j=0;j<particlePartList.size();++j) {
      vector<vector3d<float> > ppos ;
      particlePartList[j]->getParticlePositions(ppos) ;
      fwrite((const char *)(&ppos[0]),sizeof(float),3*ppos.size(),fid) ;
    }
    
    //    fwrite((const char *) (&position[0]), sizeof (float), 3*nP, fid) ;
    fwrite((const char *) (&Conn), int_size, 1, fid) ;
    for(unsigned long long int i=0;i<nP;++i) 
      fwrite((const char *) (&i), sizeof (int), 1, fid) ;
    
    fwrite((const char *) (&Cells), int_size, 1, fid) ;
    for(unsigned long long int i=0;i<nP;++i) {
      int off = i+1;
      fwrite((const char *) (&off), sizeof (int), 1, fid) ;
    }
    fwrite((const char *) (&CellChars), int_size, 1, fid) ;
    unsigned char type = 1;
    
    for(unsigned long long int i=0;i<nP;++i) {
      fwrite((const char *) (&type), sizeof (unsigned char), 1, fid) ;
    }
  
    for(si=particle_scalars.begin();si!=particle_scalars.end();++si) {
      string name = *si ;
      int comp=1 ;
      long long unsigned int ndata_size = nP*sizeof(float)*comp ;
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      for(size_t i=0;i<particlePartList.size();++i) {
	int np = particlePartList[i]->getNumParticles() ;
	vector<float> val(np,0) ;
	particlePartList[i]->getParticleScalar(name,val) ;
	fwrite(&val[0],sizeof(float),np*comp,fid) ;
      }
    }
    for(si=particle_vectors.begin();si!=particle_vectors.end();++si) {
      string name = *si ;
      int comp=3 ;
      long long unsigned int ndata_size = nP*sizeof(float)*comp ;
      fwrite((const char *) (&ndata_size), int_size, 1, fid) ;
      for(size_t i=0;i<particlePartList.size();++i) {
	int np = particlePartList[i]->getNumParticles() ;
	vector<vector3d<float> > val(np,vector3d<float>(0,0,0)) ;
	particlePartList[i]->getParticleVector(name,val) ;
	fwrite(&val[0],sizeof(float),np*comp,fid) ;
      }
    }

    fprintf(fid,"  </AppendedData>\n");
    fprintf(fid,"</VTKFile>\n");
    fclose(fid);
  }
} 
