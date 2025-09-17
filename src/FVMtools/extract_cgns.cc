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

#ifdef USE_CGNS
#include "cgnslib.h"



struct surf_info {
  string name ;
  vector<cgsize_t> trias;
  vector<cgsize_t> quads;
  vector<cgsize_t> genface;
  int nquads;
  int ntrias;
  int ngenfc;
  cgsize_t start;
  cgsize_t end;
} ;
  
void read_surf_info(string casename, string iteration,
                    vector<surf_info>& surfs) {
  FATAL(Loci::MPI_processes != 1) ;

  hid_t file_id;
  
  string gridtopo = getTopoFileName(output_dir, casename, iteration) ;

  file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to get grid topo  " << gridtopo
         << endl ;
    cerr << "does file '" << gridtopo << "' exist?" << endl ;
    Loci::Abort() ;
    exit(-1) ;
  }
 
#ifdef H5_USE_16_API  
  hid_t bndg = H5Gopen(file_id,"boundaries") ;
#else
  hid_t bndg = H5Gopen(file_id,"boundaries",H5P_DEFAULT) ;
#endif
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  
  
  
  surfs.resize(num_bcs);
  
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    
    surf_info& the_surf = surfs[bc];
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;
#ifdef H5_USE_16_API
    hid_t bcg = H5Gopen(bndg,buf) ;
#else
    hid_t bcg = H5Gopen(bndg,buf,H5P_DEFAULT) ;
#endif
    
    int nquads = sizeElementType(bcg,"quads") ;
    int ntrias = sizeElementType(bcg,"triangles") ;
    int ngeneral = sizeElementType(bcg,"nside_sizes") ;
  
    if(ntrias){
      the_surf.trias.resize(ntrias*3) ;
      {
        vector<Array<int, 3> > tmp(ntrias);
        readElementType(bcg,"triangles",tmp) ;
        int count = 0; 
        for(int i = 0; i <ntrias; i++){
          for(int j = 0; j < 3; j++){
            the_surf.trias[count++] = tmp[i][j];
          }
        }
      }
    }
  
    if(nquads){
      the_surf.quads.resize(nquads*4) ;
      {
        vector<Array<int, 4> > tmp(nquads);
        readElementType(bcg,"quads", tmp) ;
        int count = 0; 
        for(int i = 0; i <nquads; i++){
          for(int j = 0; j < 4; j++){
            the_surf.quads[count++] = tmp[i][j];
          }
        }
      }
    }
    if(ngeneral){  
      vector<int> nside_sizes(ngeneral) ;
      readElementType(bcg,"nside_sizes",nside_sizes) ;
      int nside_nodes_size = sizeElementType(bcg,"nside_nodes") ;
      vector<int> nside_nodes(nside_nodes_size) ;
      readElementType(bcg,"nside_nodes",nside_nodes) ;
      the_surf.genface.resize(nside_nodes_size + ngeneral);
      
      int count1 = 0, count2 = 0; 
      for(int i = 0; i <ngeneral; i++){
        the_surf.genface[count1++] = nside_sizes[i];
        for(int j = 0; j < nside_sizes[i]; j++){
          the_surf.genface[count1++] = nside_nodes[count2++];
        }
      }
    }
    the_surf.ngenfc = ngeneral;
    the_surf.ntrias = ntrias;
    the_surf.nquads = nquads;
    the_surf.name = string(buf);
    cout<<"boundary: " << the_surf.name << endl;
  }
    
   Loci::hdf5CloseFile(file_id) ;
}
bool read_pambient(string casename, string iteration, float& pambient){

  string filename = output_dir+"/Pambient_par." + iteration +"_" + casename ;
  
  hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                               H5F_ACC_RDONLY,
                               H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to open file '" << filename << "'!" << endl ;
    return false;
  }
  fact_db facts ;
  param<float> Pambient ;
  Loci::readContainer(file_id,"Pambient",Pambient.Rep(),EMPTY,facts) ;
  Loci::hdf5CloseFile(file_id) ;
  pambient = *Pambient;
  return true;
}
void write_surf_info(int cgfile, int cgbase, int cgzone,
                     vector<surf_info>& surfs, cgsize_t& start, cgsize_t& end) {
  char name[80];
  int cgsect, cgbc;
  int num_bc = surfs.size();
  for(int bc = 0; bc < num_bc; bc++){
   
    surf_info &the_surf = surfs[bc];
    the_surf.start = start;
    
    if(the_surf.nquads > 0){
   
      end = start + the_surf.nquads - 1;
   
      memset(name, '\0', 80) ;
      sprintf(name, "QuadElements %d", bc+1);
      if (cg_section_write(cgfile, cgbase, cgzone, name,
                           CGNS_ENUMV(QUAD_4), start, end, 0,  &(the_surf.quads[0]), &cgsect))
        cg_error_exit();
      start = end + 1;
    }
    if(the_surf.ntrias > 0){
   
       end = start + the_surf.ntrias - 1;
   
      memset(name, '\0', 80) ;
      sprintf(name, "TriElements %d", bc+1);
      if (cg_section_write(cgfile, cgbase, cgzone, name,
                           CGNS_ENUMV(TRI_3), start, end, 0, &(the_surf.trias[0]), &cgsect))
        cg_error_exit();
      start = end + 1;
    }
    
    if(the_surf.ngenfc > 0){
      end = start +  the_surf.ngenfc - 1;
      memset(name, '\0', 80) ;
      sprintf(name, "GenFaces %d", bc+1);
      if (cg_section_write(cgfile, cgbase, cgzone, name,
                           CGNS_ENUMV(NGON_n), start, end, 0, &the_surf.genface[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }
    the_surf.end = end;
  }
        


    /* write BCs */
  CGNS_ENUMT(BCType_t) bctype = CGNS_ENUMV(BCTypeUserDefined);

  for (int bc = 0; bc <num_bc ; bc++) {
    cgsize_t range[2];
    surf_info& the_surf = surfs[bc];
    range[0] = the_surf.start;
    range[1] = the_surf.end;
    memset(name, '\0', 80) ;
    snprintf(name,80, "%s", the_surf.name.c_str());

    
#if CGNS_VERSION < 3100
    if (cg_boco_write(cgfile, cgbase, cgzone, name,
                      bctype, CGNS_ENUMV(ElementRange), 2, range, &cgbc))
#else
        if (cg_boco_write(cgfile, cgbase, cgzone, name,
                bctype, CGNS_ENUMV(PointRange), 2, range, &cgbc) ||
            cg_boco_gridlocation_write(cgfile, cgbase, cgzone,
                cgbc, CGNS_ENUMV(FaceCenter)))
#endif
            cg_error_exit();
  }
}


  

bool cgnsPartConverter::processesVolumeElements() const {
  return true ; 
}
bool cgnsPartConverter::processesSurfaceElements() const {
  return true ; 
}
bool cgnsPartConverter::processesParticleElements() const {
  return false; 
}


void cgnsPartConverter::exportPostProcessorFiles(string casename,
                                                 string iteration) const {
  
  string filename = casename +"_" +iteration +".cgns";

  vector<surf_info> surfs;
  read_surf_info(casename, iteration, surfs);
  cgsize_t tot_surf_elements = 0;
  for(size_t bc = 0; bc < surfs.size(); bc++){
    surf_info& s = surfs[bc];
    tot_surf_elements += s.ntrias + s.nquads + s.ngenfc; 
  }

 
  char name[80];
  
  int cgfile, cgbase, cgzone;
  int cgcoord, cgsect;
 
  cgsize_t sizes[3],  start = 0, end = 0;
  const CGNS_ENUMT(DataType_t) dtype = CGNS_ENUMV(RealSingle);
  const CGNS_ENUMT(DataType_t) ddtype = CGNS_ENUMV(RealDouble);
  
  const CGNS_ENUMT(DataType_t) datatype = CGNS_ENUMV(Integer);
  
  
  //open file for writing and create a base
  if (cg_open(filename.c_str(), CG_MODE_WRITE, &cgfile) ||
      cg_base_write(cgfile, "Base", 3, 3, &cgbase)) cg_error_exit();
  cout << " reading pambient " << endl;
  float pambient = 1.0;
  if(read_pambient(casename, iteration, pambient)){
    cout << " writing pambient " << pambient << endl;
    cgsize_t dim = 1;
    //write user defined data
    if (cg_goto (cgfile, cgbase, NULL) ||
        cg_user_data_write ("User"))cg_error_exit ();
    if (cg_goto (cgfile, cgbase, "UserDefinedData_t", 1, "end")|| 
        cg_array_write ("Pambient", CGNS_ENUMV(RealSingle), 1, &dim, &pambient))cg_error_exit();
  }
  if( cg_goto(cgfile, cgbase, "end") ||
      cg_descriptor_write("Descriptor", "Mismatched Grid") ||
      cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
      cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
                     CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))cg_error_exit();

  
  
  //each part is written into a zone
  for(size_t i=0;i<volumePartList.size();++i) {
   
    //set sizes
    size_t pnts = volumePartList[i]->getNumNodes();
    sizes[0] = pnts;
    sizes[1] = volumePartList[i]->getNumTets() + volumePartList[i]->getNumPyrm()
      + volumePartList[i]->getNumPrsm() + volumePartList[i]->getNumHexs()+volumePartList[i]->getNumGenc() ;
    sizes[1] += tot_surf_elements;
    sizes[2] = 0;
    string pName = volumePartList[i]->getPartName();
    memset(name, '\0', 80) ;
    snprintf(name,80, "%s", pName.c_str());
    //open file  
    if(cg_zone_write(cgfile, cgbase, name, sizes,
                     CGNS_ENUMV(Unstructured), &cgzone))cg_error_exit();
    cout << " Zone: " << pName <<  " ; Index  " << cgzone <<  endl;
    {
      vector<vector3d<double> > pos ;
      volumePartList[i]->getPos(pos) ;
      if(pos.size() != pnts){
        cerr << " ERROR: volume part pos size and numNodes mismatches"<< endl;
        cg_error_exit();
      }
      
      vector<double> pos_x(pnts) ;
      //write coordinates
      for(size_t k=0;k<pnts;++k) {
        pos_x[k] = pos[k].x;
      }
      if (cg_coord_write(cgfile, cgbase, cgzone, ddtype,
                         "CoordinateX", &pos_x[0], &cgcoord))
        cg_error_exit();
    
      for(size_t k=0;k<pnts;++k) {
        pos_x[k] = pos[k].y;
      }
      if( cg_coord_write(cgfile, cgbase, cgzone, ddtype,
                         "CoordinateY", &pos_x[0], &cgcoord)) cg_error_exit();
      for(size_t k=0;k<pnts;++k) {
        pos_x[k] = pos[k].z; 
      }
      if(cg_coord_write(cgfile, cgbase, cgzone, ddtype,
                        "CoordinateZ", &pos_x[0], &cgcoord))cg_error_exit();
    }
    //write elements
    start = 1;
    const int block_size=65536 ; // Size of blocking factor

    cout << " num_surf " << surfs.size() << endl; 
    if(surfs.size() > 0){
      cout<<"writing out surf info " << endl;
      write_surf_info(cgfile, cgbase, cgzone,
                      surfs, start, end);
    }
    
    //get tets
    if(volumePartList[i]->getNumTets() > 0) { 
      int tot = volumePartList[i]->getNumTets() ;
      cout<<"num_tets " << tot << endl;
      vector<cgsize_t> elems(tot*4);
      size_t elems_index = 0; 
      int bottom = 0 ;
      int top = tot + volumePartList[i]->getNumTetsIblank() ;
      while(bottom < top) {
        vector<Array<int,4> > temp ;
        volumePartList[i]->getTetBlock(temp,bottom,block_size) ;
        bottom += block_size ;
        if(temp.size() > 0) {
         
          for(size_t j=0;j<temp.size();++j) {
            for(int k = 0; k < 4; k++){
              elems[elems_index++]  = static_cast<cgsize_t>(temp[j][k]) ;
            }
          }
         
        }
      }
      end = start + tot - 1;
      memset(name, '\0', 80) ;
      sprintf(name, "TetElements");
      if (cg_section_write(cgfile, cgbase, cgzone, name,
                           CGNS_ENUMV(TETRA_4), start, end, 0, &elems[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }
    
   
    if(volumePartList[i]->getNumPyrm() > 0) { 
      int tot = volumePartList[i]->getNumPyrm() ;
      vector<cgsize_t> elems(tot*5);
      size_t elems_index = 0;
      cout<<"num_pyrm " << tot << endl;
      int bottom = 0 ;
      int top = tot + volumePartList[i]->getNumPyrmIblank() ;
      while(bottom < top) {
        vector<Array<int,5> > temp ;
        volumePartList[i]->getPyrmBlock(temp,bottom,block_size) ;
        bottom += block_size ;
        if(temp.size() > 0) {
          for(size_t j=0;j<temp.size();++j) {
            for(int k = 0; k < 5; k++){
              elems[elems_index++]  = static_cast<cgsize_t>(temp[j][k]);
            }
          }

         
          
        }
      }

      memset(name, '\0', 80) ;
      sprintf(name, "PyraElements");
      end = start + tot - 1;
      if (cg_section_write(cgfile, cgbase, cgzone, name,
                           CGNS_ENUMV(PYRA_5), start, end, 0, &elems[0], &cgsect))
        cg_error_exit();
      start = end + 1; 
    }
    
    if(volumePartList[i]->getNumPrsm() > 0) {
      int tot = volumePartList[i]->getNumPrsm() ;
      vector<cgsize_t> elems(tot*6);
      size_t elems_index = 0;
      cout<<"num_prsm " << tot << endl;
      int bottom = 0 ;
      int top = tot + volumePartList[i]->getNumPrsmIblank() ;
      while(bottom < top) {
        vector<Array<int,6> > temp ;
        volumePartList[i]->getPrsmBlock(temp,bottom,block_size) ;
        bottom += block_size ;
        if(temp.size() > 0) {
         
          for(size_t j=0;j<temp.size();++j) {
            for(int k = 0; k < 6; k++){
              elems[elems_index++]  = static_cast<cgsize_t>(temp[j][k]);
            }
          }

          
        }
      }

      memset(name, '\0', 80) ;
      sprintf(name, "PentaElements");
      end = start + tot - 1;
      if (cg_section_write(cgfile, cgbase, cgzone, name,
                           CGNS_ENUMV(PENTA_6), start, end, 0,&elems[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }

  
    if(volumePartList[i]->getNumHexs() > 0) {
      int tot = volumePartList[i]->getNumHexs() ;
      vector<cgsize_t> elems(tot*8);
      size_t elems_index = 0;
      cout<<"num_hex " << tot << endl;
      int bottom = 0 ;
      int top = tot + volumePartList[i]->getNumHexsIblank() ;
      while(bottom < top) {
        vector<Array<int,8> > temp ;
        volumePartList[i]->getHexBlock(temp,bottom,block_size) ;
        bottom += block_size ;
        if(temp.size() > 0) {
          for(size_t j=0;j<temp.size();++j) {
            for(int k = 0; k < 8; k++){
              elems[elems_index++] = static_cast<cgsize_t>(temp[j][k]) ;
            }
          }
          
        }
      }
      memset(name, '\0', 80) ;
      sprintf(name, "HexaElements");
      end = start + tot - 1;
      if (cg_section_write(cgfile, cgbase, cgzone, name,
                           CGNS_ENUMV(HEXA_8), start, end, 0, &elems[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }
    //don't know if vog format and cgns format match
    //In cgns, first general face is defined, and then general face is defined, and the face id should be negative if its normal points inward
    if(volumePartList[i]->getNumGenc() > 0) { // write out general cells
      int tot = volumePartList[i]->getNumGenc();
      cout<<"num_genc " << tot << endl;
      vector<int> genCellNfaces;//for each general cell, the num_of_face it has
      vector<int> genCellNsides;//for each general face, the num_of_node it has
      vector<int>  genCellNodes ;//face2node. 
      volumePartList[i]->getGenCell(genCellNfaces,genCellNsides,genCellNodes) ;

      size_t num_faces = genCellNsides.size();
      //first write out face2node of general faces
      vector<cgsize_t> face2node;
      int cnt = 0 ;
      for(size_t f=0;f<num_faces;++f) {
        face2node.push_back(genCellNsides[f]);
        for(int n = 0; n <genCellNsides[f]; n++){
          face2node.push_back(genCellNodes[cnt++]);
        }
      }
      end = start + num_faces - 1;
      int face_offset = start;
      if (cg_section_write(cgfile, cgbase, cgzone, "GenFaces",
                           CGNS_ENUMV(NGON_n), start, end, 0, &face2node[0], &cgsect))
        cg_error_exit();
      start = end + 1;
      //adjust the index of cell2face
      vector<cgsize_t> cell2face;
      cnt = 0;
      for(int c = 0; c < tot; ++c) {
        cell2face.push_back(genCellNfaces[c]);//num_face the cell has
        for(int f = 0; f < genCellNfaces[c]; f++){
          cell2face.push_back(cnt+face_offset);
          cnt++;
        }
      }
      //then write out cell2face
      end = start + tot - 1;
      if (cg_section_write(cgfile, cgbase, cgzone, "GenElements",
                           CGNS_ENUMV(NFACE_n), start, end, 0, &cell2face[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }

   
    vector<string> nscalars = volumePartList[i]->getNodalScalarVars() ;
    vector<string> nvectors = volumePartList[i]->getNodalVectorVars() ;
    if(nscalars.size() > 0 || nvectors.size() > 0){
      int isol, ifld;
      if (cg_sol_write(cgfile, cgbase, cgzone,
                       "NodeVariables", CGNS_ENUMV(Vertex), &isol))cg_error_exit();
      for(size_t j = 0; j < nscalars.size(); j++){//each variable is a field
        string varname = nscalars[j];
        vector<float> vals ;
        volumePartList[i]->getNodalScalar(varname,vals) ;
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", varname.c_str());
        
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals[0], &ifld))cg_error_exit();
       
      }
      
     
      for(size_t j = 0; j < nvectors.size(); j++){//each variable is a field
        string varname = nvectors[j];
        vector<vector3d<float> > vals ;
        volumePartList[i]->getNodalVector(varname,vals) ;
        vector<float> vals_c(vals.size());
        for(size_t k = 0; k < vals.size(); k++)vals_c[k] = vals[k].x;
        string nameX = varname + "X";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameX.c_str()); 
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals_c[0], &ifld))cg_error_exit();
        string nameY = varname + "Y";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameY.c_str());
        for(size_t k = 0; k < vals.size(); k++)vals_c[k] = vals[k].y;
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals_c[0], &ifld))cg_error_exit();
        string nameZ = varname + "Z";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameZ.c_str());
        for(size_t k = 0; k < vals.size(); k++)vals_c[k] = vals[k].z;
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals_c[0], &ifld))cg_error_exit(); 
      }
    }
    //write out fake ids
    if(id_required){
      int isol, ifld;
      if (cg_sol_write(cgfile, cgbase, cgzone,
                       "ElementVariables", CGNS_ENUMV(CellCenter), &isol))cg_error_exit();
      size_t num_elements = end;
      vector<int> ids(num_elements, 1);
      memset(name, '\0', 80) ;
      sprintf(name, "ElementIds");
      if (cg_field_write(cgfile, cgbase, cgzone, isol,
                         datatype, name, &ids[0], &ifld))cg_error_exit(); 
    }
  }
 
  
  //each surface part is written into a zone
  for(size_t i=0;i<surfacePartList.size();++i) {

    //set sizes
    size_t pnts = surfacePartList[i]->getNumNodes();
    sizes[0] = pnts;
    sizes[1] = surfacePartList[i]->getNumTrias() + surfacePartList[i]->getNumQuads() +surfacePartList[i]-> getNumGenfc();
    sizes[2] = 0;
    string pName = surfacePartList[i]->getPartName();
    memset(name, '\0', 80) ;
    snprintf(name,80, "%s", pName.c_str());
    
    

    //open file  
   
    if(cg_zone_write(cgfile, cgbase, name, sizes,
                     CGNS_ENUMV(Unstructured), &cgzone))cg_error_exit();
    cout << " Zone: " << pName <<  " ; Index " << cgzone <<  endl;

    {
      vector<vector3d<double> > pos ;
      surfacePartList[i]->getPos(pos) ;
     
      //write coordinates
  
      vector<double> pos_x(pnts) ;
      for(size_t j=0;j<pnts;++j) {
        pos_x[j] = pos[j].x;
      }
      if (cg_coord_write(cgfile, cgbase, cgzone, ddtype,
                         "CoordinateX", &pos_x[0], &cgcoord))
        cg_error_exit();

      for(size_t j=0;j<pnts;++j) {
        pos_x[j] = pos[j].y;
      }
    
      if( cg_coord_write(cgfile, cgbase, cgzone, ddtype,
                         "CoordinateY", &pos_x[0], &cgcoord)) cg_error_exit();
        
      for(size_t j=0;j<pnts;++j) {
        pos_x[j] = pos[j].z; 
      }
      if(cg_coord_write(cgfile, cgbase, cgzone, ddtype,
                        "CoordinateZ", &pos_x[0], &cgcoord))cg_error_exit();
    }
 
    //write elements
    start = 1;

    //get quads
   
    
    if(surfacePartList[i]->getNumQuads() > 0) { 
      int tot = surfacePartList[i]->getNumQuads() ;
      vector<Array<int,4> > elems;
      surfacePartList[i]->getQuads(elems) ;
      vector<cgsize_t> vals(elems.size()*4);
      size_t cnt = 0;
      for(size_t j = 0; j < elems.size(); j++){
        for(int k = 0; k < 4; k++){
          vals[cnt++] = elems[j][k];
        }
      }
      end = start + tot - 1;
      if (cg_section_write(cgfile, cgbase, cgzone, "QuadElements",
                           CGNS_ENUMV(QUAD_4), start, end, 0, &vals[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }

    if(surfacePartList[i]->getNumTrias() > 0) { 
      int tot = surfacePartList[i]->getNumTrias() ;
      vector<Array<int,3> > elems;
      surfacePartList[i]->getTrias(elems) ;
      vector<cgsize_t> vals(elems.size()*3);
      size_t cnt = 0;
      for(size_t j = 0; j < elems.size(); j++){
        for(int k = 0; k < 3; k++){
          vals[cnt++] = elems[j][k];
        }
      }
      end = start + tot - 1;
      if (cg_section_write(cgfile, cgbase, cgzone, "TriElements",
                           CGNS_ENUMV(TRI_3), start, end, 0, &vals[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }
   
   
    if(surfacePartList[i]->getNumGenfc() > 0) {
      int tot = surfacePartList[i]->getNumGenfc() ;
      vector<int> numGenFnodes;
      vector<int> genNodes;
      surfacePartList[i]->getGenf(numGenFnodes, genNodes);

      vector<cgsize_t> elems;
      int cnt = 0;
      for(size_t j = 0; j < numGenFnodes.size(); j++){
        elems.push_back(numGenFnodes[j]);
        for(int k = 0; k < numGenFnodes[j]; k++)elems.push_back(genNodes[cnt++]);
      }
            
      end = start + tot - 1;
      if (cg_section_write(cgfile, cgbase, cgzone, "GenFaces",
                           CGNS_ENUMV(NGON_n), start, end, 0,&elems[0], &cgsect))
        cg_error_exit();
      start = end + 1;
    }
   
    vector<string> nscalars = surfacePartList[i]->getNodalScalarVars() ;
    vector<string> nvectors = surfacePartList[i]->getNodalVectorVars() ;
    if(nscalars.size() > 0 || nvectors.size() > 0){
      int isol, ifld;
      if (cg_sol_write(cgfile, cgbase, cgzone,
                       "NodeVariables", CGNS_ENUMV(Vertex), &isol))cg_error_exit();
    
      for(size_t j = 0; j < nscalars.size(); j++){//each variable is a field
        string varname = nscalars[j];
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", varname.c_str()); 
        vector<float> vals ;
        surfacePartList[i]->getNodalScalar(varname,vals) ;
        
    
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals[0], &ifld))cg_error_exit();
      }
      
     
      for(size_t j = 0; j < nvectors.size(); j++){//each variable is a field
        string varname = nvectors[j];
        vector<vector3d<float> > vals ;
        surfacePartList[i]->getNodalVector(varname,vals) ;
        
        string nameX = varname+"X";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameX.c_str());
        vector<float> vals_c(vals.size());
        for(size_t k = 0; k < vals.size(); k++)vals_c[k] = vals[k].x;
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals_c[0], &ifld))cg_error_exit();

        string nameY = varname+"Y";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameY.c_str());
        for(size_t k = 0; k < vals.size(); k++)vals_c[k] = vals[k].y;
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals_c[0], &ifld))cg_error_exit();

        
        string nameZ = varname+"Z";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameZ.c_str());
        for(size_t k = 0; k < vals.size(); k++)vals_c[k] = vals[k].z;
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals_c[0], &ifld))cg_error_exit();
       
      }
    }

    vector<string> escalars = surfacePartList[i]->getElementScalarVars() ;
    vector<string> evectors = surfacePartList[i]->getElementVectorVars() ;
    {
      int isol, ifld;
      if (cg_sol_write(cgfile, cgbase, cgzone,
                       "ElementVariables", CGNS_ENUMV(CellCenter), &isol))cg_error_exit();
      for(size_t j = 0; j < escalars.size(); j++){//each variable is a field
        string varname = escalars[j];
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", varname.c_str());
        
        vector<float> qvals;
        vector<float> tvals; 
        vector<float> gvals;
      
        surfacePartList[i]->getElementScalar(varname,qvals, tvals, gvals) ;
      
        size_t num_elems = qvals.size()+ tvals.size()+gvals.size();
        vector<float> vals(num_elems);
        size_t cnt = 0;
        for(size_t k = 0; k < qvals.size(); k++)vals[cnt++] = qvals[k];
        for(size_t k = 0; k < tvals.size(); k++)vals[cnt++] = tvals[k];
        for(size_t k = 0; k < gvals.size(); k++)vals[cnt++] = gvals[k];
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals[0], &ifld))cg_error_exit();
      }
      //write out element id
      if(id_required) {
        vector<int> quads_ids;
        vector<int> trias_ids;
        vector<int> genface_ids;
        surfacePartList[i]->getQuadsIds(quads_ids);
        surfacePartList[i]->getTriasIds(trias_ids);
        surfacePartList[i]->getGenfIds(genface_ids);
        size_t num_elems =quads_ids.size() + trias_ids.size() + genface_ids.size();
        
        vector<int> vals(num_elems);
        size_t cnt = 0;
        for(size_t k = 0; k < quads_ids.size(); k++)vals[cnt++] = quads_ids[k];
        for(size_t k = 0; k < trias_ids.size(); k++)vals[cnt++] = trias_ids[k];
        for(size_t k = 0; k < genface_ids.size(); k++)vals[cnt++] = genface_ids[k];
        
        memset(name, '\0', 80) ;
        sprintf(name, "ElementIds");
        
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           datatype, name, &vals[0], &ifld))cg_error_exit();
      }       
      for(size_t j = 0; j < evectors.size(); j++){//each variable is a field
        string varname = evectors[j];
        vector<vector3d<float> > qvals;
        vector<vector3d<float> > tvals;
        vector<vector3d<float> > gvals;
        surfacePartList[i]->getElementVector(varname,qvals, tvals, gvals) ;

        string nameX = varname+"X";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameX.c_str());
    
        size_t num_elems = qvals.size()+ tvals.size()+gvals.size();
        vector<float> vals(num_elems) ;
        size_t cnt = 0;
        for(size_t k = 0; k < qvals.size(); k++)vals[cnt++] = qvals[k].x;
        for(size_t k = 0; k < tvals.size(); k++)vals[cnt++] = tvals[k].x;
        for(size_t k = 0; k < gvals.size(); k++)vals[cnt++] = gvals[k].x;
      
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals[0], &ifld))cg_error_exit();

        string nameY = varname+"Y";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameY.c_str());
     
        cnt = 0;
        for(size_t k = 0; k < qvals.size(); k++)vals[cnt++] = qvals[k].y;
        for(size_t k = 0; k < tvals.size(); k++)vals[cnt++] = tvals[k].y;
        for(size_t k = 0; k < gvals.size(); k++)vals[cnt++] = gvals[k].y;
      
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals[0], &ifld))cg_error_exit();

        string nameZ = varname+"Z";
        memset(name, '\0', 80) ;
        snprintf(name,80, "%s", nameZ.c_str());
    
        cnt = 0;
        for(size_t k = 0; k < qvals.size(); k++)vals[cnt++] = qvals[k].z;
        for(size_t k = 0; k < tvals.size(); k++)vals[cnt++] = tvals[k].z;
        for(size_t k = 0; k < gvals.size(); k++)vals[cnt++] = gvals[k].z;
      
        if (cg_field_write(cgfile, cgbase, cgzone, isol,
                           dtype, name, &vals[0], &ifld))cg_error_exit();

      }
    }
    
  }
  
  cg_close(cgfile); 
}

#else

bool cgnsPartConverter::processesVolumeElements() const {
  cerr << "CGNS support not compiled for Loci utilities.  Check setup and recompile to use CGNS support." << endl ;
  return false ; 
}
bool cgnsPartConverter::processesSurfaceElements() const {
  cerr << "CGNS support not compiled for Loci utilities.  Check setup and recompile to use CGNS support." << endl ;
  return false ; 
}
bool cgnsPartConverter::processesParticleElements() const {
  return false; 
}


void cgnsPartConverter::exportPostProcessorFiles(string casename,
                                                 string iteration) const {
  
}
#endif


  
 
 


// ... the end ...


