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
#include <algorithm>
#include <set>
#include <sstream> 
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
using std::set;
using std::stringstream;
#include "extract.h"

#include <sys/types.h>
#include <sys/stat.h>


/* Numeric tags (codes) for FIELDVIEW binary file format. */
 
#define MAJOR_VERSION   3
#define MINOR_VERSION   2
#define FV_MAGIC        0x00010203      /* decimal 66051 */
 
/* Content of the file (grid only, results only or combined). */
#define FV_GRIDS_FILE           1
#define FV_RESULTS_FILE         2
#define FV_COMBINED_FILE        3
 
#define FV_NODES                1001
#define FV_FACES                1002
#define FV_ELEMENTS             1003
#define FV_VARIABLES            1004
#define FV_BNDRY_VARS           1006
#define FV_ARB_POLY_FACES       1007
#define FV_ARB_POLY_ELEMENTS    1008
#define FV_ARB_POLY_BNDRY_VARS  1009
 
#define FV_TET_ELEM_ID          1
#define FV_HEX_ELEM_ID          2
#define FV_PRISM_ELEM_ID        3
#define FV_PYRA_ELEM_ID         4
#define FV_ARB_POLY_ELEM_ID     5
 
/* Values for "wall_info" array (see comments in fv_encode_elem_header). */
#ifdef __STDC__
#define A_WALL         (07u)
#define NOT_A_WALL     (0u)
#else
#define A_WALL         (07)
#define NOT_A_WALL     (0)
#endif

/* Don't change these - used by fv_encode_elem_header ! */
#define MAX_NUM_ELEM_FACES     6
#define BITS_PER_WALL  3
#define ELEM_TYPE_BIT_SHIFT    (MAX_NUM_ELEM_FACES*BITS_PER_WALL)
 
/*
** fv_encode_elem_header:  return an encoded binary element header
**
** Input:
**    elem_type:  integer element type as shown in fv_reader_tags.h
**    wall_info:  array of integer "wall" flags, one for each face of
**                the element.  The wall flags are used during streamline
**                calculation.  Currently, the only meaningful values are
**                A_WALL and NOT_A_WALL as shown in fv_reader_tags.h.
**                Streamlines are forced away from faces marked as
**                "A_WALL", by limiting velocity and position very near
**                the wall.
** Output:
**    Function return value is the encoded binary element header.
*/
 
unsigned int fv_encode_elem_header (int elem_type, int wall_info[])
{
  unsigned int header;
  int i, nfaces;
 
  switch (elem_type)
    {
    case FV_TET_ELEM_ID:
      header = (1 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 4;
      break;
    case FV_HEX_ELEM_ID:
      header = (4 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 6;
      break;
    case FV_PRISM_ELEM_ID:
      header = (3 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 5;
      break;
    case FV_PYRA_ELEM_ID:
      header = (2 << ELEM_TYPE_BIT_SHIFT);
      nfaces = 5;
      break;
    default:
      fprintf(stderr, "ERROR:  Unknown element type\n");
      return 0;
    }
 
  for (i = 0; i < nfaces; i++)
    {
      unsigned int u = wall_info[i];
      if (u > A_WALL)
        {
          fprintf(stderr, "ERROR:  Bad wall value\n");
          return 0;
        }
      header |= (u << (i*BITS_PER_WALL));
    }
  return header;
}

string convert_fv_compatible(string var) {
  if(var == "i")
    return string("i_") ;
  if(var == "j")
    return string("j_") ;
  if(var == "k")
    return string("k_") ;
  if(var == "x")
    return string("x_") ;
  if(var == "y")
    return string("y_") ;
  if(var == "z")
    return string("z_") ;
  if(var == "R")
    return string("R_") ;
  return var ;
}


bool fieldViewPartConverter::processesVolumeElements() const { return true ;}
bool fieldViewPartConverter::processesSurfaceElements() const { return true ;}
bool fieldViewPartConverter::processesParticleElements() const { return true ;}

void fieldViewPartConverter::
exportPostProcessorFiles(string casename, string iteration) const {

  string dirname ;
  string filename ;
  
  size_t npnts = 0 ;//in voluem parts and in surface parts
  size_t ntets = 0, nprsm = 0, npyrm = 0, nhexs = 0 ; //, ngen = 0 ;
  
  
  vector<vector<int> > part_nodes ;
  vector<int> elem_ids ;
  bool general_boundary = false ;
  bool first_var = true ;
  bool first_boundary  = true;
  FILE *OFP ;
  
  int ibuf[10] ;
  char fv[80] ;

  set<string> nodal_scalars ;
  set<string> nodal_vectors ;
  set<string> element_scalars ;
  set<string> element_vectors ;
 
  //open grid file and write the header
  { 
    dirname = casename + "_fv."+iteration ;
    struct stat statbuf ;
    if(stat(dirname.c_str(),&statbuf))
      mkdir(dirname.c_str(),0755) ;
    else
      if(!S_ISDIR(statbuf.st_mode)) {
        cerr << "file " << dirname
             << " should be a directory!, rename it and start again."
             << endl ;
        exit(-1) ;
      }
    
    filename = casename + "_fv_"+iteration + ".bin";
    filename = dirname + "/" + filename ;
    
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
      //      ngen += volumePartList[i]->getNumGenc();
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
      npnts += surfacePartList[i]->getNumNodes() ;
    }
   


//    double time = atof(iteration.c_str()); //? how to set up this value?
    double time = 0 ;
    int ncycle = 0 ;
    string timefile = "output/timestep_txt." + iteration + "_" + casename ;
    std::ifstream timein(timefile.c_str(),ios::in) ;
    if(!timein.fail())
      timein >> ncycle >> time ;
    
    OFP = fopen(filename.c_str(), "wb") ;
    if(OFP == NULL) {
      cerr << "can't open file " << filename << endl ;
      exit(-1) ;
    }
    
    // Write out the magic number and the fieldview version info 
    ibuf[0] = FV_MAGIC ;
    fwrite(ibuf, sizeof(int), 1, OFP) ;
    memset(fv,'\0',80) ;
    snprintf(fv,80, "FIELDVIEW") ;
    fwrite(&fv, sizeof(char), 80, OFP) ;
    ibuf[0] = 3 ;
    ibuf[1] = 0 ;
    fwrite(ibuf, sizeof(int), 2,OFP) ;

    ibuf[0] = FV_COMBINED_FILE ;
    fwrite(ibuf,sizeof(int),1,OFP) ;

    ibuf[0] = 0 ;
    fwrite(ibuf,sizeof(int),1,OFP) ;

    float TIME, FSMACH, ALPHA, RE ;
    TIME = time ;
    FSMACH = 0.0 ;
    ALPHA = 0.0 ;
    RE = 0.0 ;

    fwrite(&TIME, sizeof(float), 1, OFP) ;
    fwrite(&FSMACH, sizeof(float), 1, OFP) ;
    fwrite(&ALPHA, sizeof(float), 1, OFP) ;
    fwrite(&RE, sizeof(float), 1, OFP) ;
  
    ibuf[0] = 1 ; // Number of grids 
    fwrite(ibuf, sizeof(int), 1, OFP) ;

    vector<string> nlist ;
    vector<string> blist ;

    set<string>::const_iterator si ;
    for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {
      const string var(convert_fv_compatible(*si)) ;
      nlist.push_back(var) ;
    }
    for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
      const string var(convert_fv_compatible(*si)) ;
      string v = var ;
      nlist.push_back(v + "x ; " + v) ;
      nlist.push_back(v + "y") ;
      nlist.push_back(v + "z") ;
    }
    for(si=element_scalars.begin();si!=element_scalars.end();++si) {
      const string var(convert_fv_compatible(*si)) ;
      blist.push_back(var) ;
    }
    for(si=element_vectors.begin();si!=element_vectors.end();++si) {
      const string var(convert_fv_compatible(*si)) ;
      string v = var ;
      blist.push_back(v + "x ; " + v) ;
      blist.push_back(v + "y") ;
      blist.push_back(v + "z") ;
    }
  
  
    // Number of face types (boundary flags)
    ibuf[0] = surfacePartList.size() ;
    fwrite(ibuf, sizeof(int), 1, OFP) ;
    
    for(size_t i=0;i<surfacePartList.size();++i) {
      ibuf[0] = 0 ; // zero, no boundary variables
      if(blist.size() > 0)
        ibuf[0] = 1 ;
      ibuf[1] = 1 ; // Normals should be consistent
      fwrite(ibuf,sizeof(int),2,OFP) ;
      memset(fv,'\0',80) ;
      snprintf(fv,80, "%s",(surfacePartList[i]->getPartName()).c_str()) ;
      fwrite(&fv, sizeof(char), 80, OFP) ;
    }

    // Nodal variables
    ibuf[0] = nlist.size() ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    for(size_t i=0;i<nlist.size();++i) {
      memset(fv,'\0',80) ;
      snprintf(fv,80, "%s",nlist[i].c_str()) ;
      fwrite(&fv, sizeof(char), 80, OFP) ;
    }
    // boundary variables 
    ibuf[0] = blist.size() ;
    fwrite(ibuf,sizeof(int),1,OFP) ;
    for(size_t i=0;i<blist.size();++i) {
      memset(fv,'\0',80) ;
      snprintf(fv,80, "%s",blist[i].c_str()) ;
      fwrite(&fv, sizeof(char), 80, OFP) ;
    }
  }//finish writing header
  

  
  int nparts = volumePartList.size() + surfacePartList.size();
  vector<int> node_offset(nparts);

  cout << "npnts = " << npnts << endl ;
  //merge the node position in all volume parts and surface parts in one vector
  //and write it into file
  {
    //put the nodes from volume parts and surface parts into
    //pos_x, pos_y and pos_z
    vector<float> pos_x(npnts),pos_y(npnts),pos_z(npnts) ;
    size_t tot = 0;
    for(size_t i=0;i<volumePartList.size();++i) {
      vector<vector3d<float> > part_pos;
      volumePartList[i]->getPos(part_pos) ;
      for(size_t j = 0; j<part_pos.size(); j++){
        pos_x[j+tot] = part_pos[j].x;
        pos_y[j+tot] = part_pos[j].y;
        pos_z[j+tot] = part_pos[j].z; 
      }
      node_offset[i] = tot;
      tot += part_pos.size();
    }
    
   
  
    for(size_t i =0;i<surfacePartList.size();++i) {
      vector<vector3d<float> > part_pos;
      surfacePartList[i]->getPos(part_pos) ;
      for(size_t j = 0; j<part_pos.size(); j++){
        pos_x[j+tot] = part_pos[j].x;
        pos_y[j+tot] = part_pos[j].y;
        pos_z[j+tot] = part_pos[j].z; 
      }
      node_offset[i+volumePartList.size()] = tot;
      tot +=part_pos.size();
    }
    if(tot != npnts)cerr<<"ERROR: number of points is not consisitent" << endl;
    size_t pts = npnts ;
   
    //write pos_x, pos_y and pos_z into file
    {
     
      ibuf[0] = FV_NODES ;
      ibuf[1] = pts ;
      fwrite(ibuf,sizeof(int),2,OFP) ;
      for(size_t i=0;i<pts;++i) {
        float x = pos_x[i];
        fwrite(&x,sizeof(float),1,OFP) ;
      }
      for(size_t i=0;i<pts;++i) {
        float y = pos_y[i] ;
        fwrite(&y,sizeof(float),1,OFP) ;
      }
      for(size_t i=0;i<pts;++i) {
        float z = pos_z[i];
        fwrite(&z,sizeof(float),1,OFP) ;
      }
    }
  }//finish writing positions
  
  //write boundary elements
  vector<int> face_offset(surfacePartList.size()) ;
  int face_id = 0;
  part_nodes.resize(surfacePartList.size());
  for(size_t part_id =1;part_id<=surfacePartList.size();++part_id) {
    string name = surfacePartList[part_id-1]->getPartName();
    int part_num = volumePartList.size()+part_id -1;//index to node_offset
    size_t npt = surfacePartList[part_id-1]->getNumNodes() ;
    int local_id = 0;
   
    int num_ord_face = surfacePartList[part_id -1]->getNumQuads()
      +surfacePartList[part_id -1]->getNumTrias();

    //put trias and quads into ordinary_faces
    //in trias and quads, the index is local, and starts with 1
    //in part_nodes, the index is global, and starts with 0
    //in ordinary_faces, the index is the same as in trias and quads
    //the connectivity written into files is actually in global numbering
    vector<Array<int,4> > ordinary_faces(num_ord_face);
    part_nodes[part_id-1].resize(npt) ;//the node_set 
    for(size_t j=0;j < npt;++j){
      part_nodes[part_id-1][j] = node_offset[part_num]+j+1 ;
    }
    face_offset[part_id-1] = face_id;
        
    size_t nquads = surfacePartList[part_id -1]->getNumQuads() ;
    if(nquads > 0) {
      vector<Array<int,4> > quads ;
      surfacePartList[part_id-1]->getQuads(quads) ; 
      for(size_t i=0;i<quads.size();++i) {
        Array<int,4> a ;
        a[0] = quads[i][0] ;
        a[1] = quads[i][1] ;
        a[2] = quads[i][2] ;
        a[3] = quads[i][3] ;
        ordinary_faces[local_id++] = a ;
        elem_ids.push_back(face_id++) ;
      }
    }


    size_t ntrias = surfacePartList[part_id -1]->getNumTrias() ;
    if(ntrias > 0) {
      vector<Array<int,3> > trias ;
      surfacePartList[part_id -1]->getTrias(trias) ; 

      for(size_t i=0;i<ntrias;++i) {
        Array<int,4> a ;
        a[0] = trias[i][0] ;
        a[1] = trias[i][1] ;
        a[2] = trias[i][2] ;
        a[3] = 0 ;

        ordinary_faces[local_id++] = a ;
        elem_ids.push_back(face_id++) ;
      }
    }

    //process general face, and write all faces on this boundary into file 
    size_t ngeneral = surfacePartList[part_id -1]->getNumGenfc() ;
    vector<int> nside_sizes,nside_nodes ;
    if(ngeneral > 0) {
      surfacePartList[part_id-1]->getGenf(nside_sizes,nside_nodes) ;
    }
    for(size_t i=0;i<ngeneral;++i) {
      elem_ids.push_back(face_id++) ;
    }
     
    ibuf[1] = part_id ;
    ibuf[2] = ngeneral+ordinary_faces.size() ;
    
    if(ngeneral == 0) {
      ibuf[0] = FV_FACES ;
      fwrite(ibuf,sizeof(int),3,OFP) ;

      for(size_t i=0;i<ordinary_faces.size();++i) {
	ibuf[0]=part_nodes[part_id-1][ordinary_faces[i][0]-1] ;
	ibuf[1]=part_nodes[part_id-1][ordinary_faces[i][1]-1] ;
	ibuf[2]=part_nodes[part_id-1][ordinary_faces[i][2]-1] ;
	if(ordinary_faces[i][3] == 0) {
	  ibuf[3]= 0 ;
	} else {
	  ibuf[3]=part_nodes[part_id-1][ordinary_faces[i][3]-1] ;
	}

	fwrite(ibuf,sizeof(int),4,OFP) ;
      }
    } else {
       ibuf[0] = FV_ARB_POLY_FACES ;
       fwrite(ibuf,sizeof(int),3,OFP) ;
       general_boundary = true ;
       for(size_t i=0;i<ordinary_faces.size();++i) {
         if(ordinary_faces[i][3] == 0) {
           ibuf[0] = 3 ;
           fwrite(ibuf,sizeof(int),1,OFP) ;
            ibuf[0]= part_nodes[part_id-1][ordinary_faces[i][0]-1] ;
            ibuf[1]= part_nodes[part_id-1][ordinary_faces[i][1]-1] ;
            ibuf[2]= part_nodes[part_id-1][ordinary_faces[i][2]-1] ;
            fwrite(ibuf,sizeof(int),3,OFP) ;
         } else {
           ibuf[0] = 4 ;
            fwrite(ibuf,sizeof(int),1,OFP) ;
            ibuf[0]= part_nodes[part_id-1][ordinary_faces[i][0]-1] ;
            ibuf[1]= part_nodes[part_id-1][ordinary_faces[i][1]-1] ;
            ibuf[2]= part_nodes[part_id-1][ordinary_faces[i][2]-1] ;
            ibuf[3]= part_nodes[part_id-1][ordinary_faces[i][3]-1] ;
            fwrite(ibuf,sizeof(int),4,OFP) ;
         }
       }
       
       int cnt = 0 ;
       for(size_t i=0;i<ngeneral;++i) {
         ibuf[0] = nside_sizes[i] ;
         fwrite(ibuf,sizeof(int),1,OFP) ;
         for(int j=0;j<nside_sizes[i];++j) {
           ibuf[0] = part_nodes[part_id-1][nside_nodes[cnt++]-1] ;
           fwrite(ibuf,sizeof(int),1,OFP) ;
         }
       }
     }
  }//finish writing boundary elements
    

  //prepare to write volume elements
  {
    ibuf[0] = FV_ELEMENTS ;
    ibuf[1] = ntets ;
    ibuf[2] = nhexs ;
    ibuf[3] = nprsm ;
    ibuf[4] = npyrm ;
    fwrite(ibuf,sizeof(int),5,OFP) ;
  }
   
  //write tets
  static int tet_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                              NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
  
  unsigned int elem_header = fv_encode_elem_header(FV_TET_ELEM_ID,
                                                   tet_walls) ;

  
  const int block_size=65536 ; // Size of blocking factor
  for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumTets() > 0) { // write out tets
        int tot = volumePartList[i]->getNumTets() ;
        int start = 0 ;
        int top = tot + volumePartList[i]->getNumTetsIblank() ;
        while(start < top) {
          vector<Array<int,4> > tets ;
          volumePartList[i]->getTetBlock(tets,start,block_size) ;
          start += block_size ;
          if(tets.size() > 0) {
            for(size_t j=0;j<tets.size();++j) {
              for(int k = 0; k < 4; k++){
                tets[j][k] += node_offset[i];
              }
            }
            for(size_t j=0;j<tets.size();++j) {
              fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
              fwrite(&tets[j][0],sizeof(int),4,OFP) ;
            }
          }
        }
      }
  }
    
  //write pyrm 
    static int pyrm_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                                 NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
      
    elem_header = fv_encode_elem_header(FV_PYRA_ELEM_ID,
                                        pyrm_walls) ;
    
    for(size_t i=0;i<volumePartList.size();++i) { 
      if(volumePartList[i]->getNumPyrm() > 0) {
        int tot = volumePartList[i]->getNumPyrm() ;
        int start = 0 ;
        int top = tot + volumePartList[i]->getNumPyrmIblank() ;
        while(start < top) {
          vector<Array<int,5> > pyrm ;
          volumePartList[i]->getPyrmBlock(pyrm,start,block_size) ;
          start += block_size ;
          for(size_t j=0;j<pyrm.size();++j) {
            for(int k = 0; k < 5; k++){
              pyrm[j][k] += node_offset[i];
            }
          }
          for(size_t j=0;j<pyrm.size();++j) {
            fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
            fwrite(&pyrm[j][0],sizeof(int),5,OFP) ;
          }
        }
      }
    }

    //write prsm
    static int prsm_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                                 NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
      
    elem_header = fv_encode_elem_header(FV_PRISM_ELEM_ID,
                                        prsm_walls) ;
      
    for(size_t i=0;i<volumePartList.size();++i) { 
      if(volumePartList[i]->getNumPrsm() > 0) { // write out prsm
        int tot = volumePartList[i]->getNumPrsm() ;
        int start = 0 ;
        int top = tot + volumePartList[i]->getNumPrsmIblank() ;
        while(start < top) {
          vector<Array<int,6> > prsm ;
          volumePartList[i]->getPrsmBlock(prsm,start,block_size) ;
          start += block_size ;
          for(size_t j=0;j<prsm.size();++j) {
            for(int k = 0; k < 6; k++){
              prsm[j][k] += node_offset[i];
            }
          }
      
          for(size_t j=0;j<prsm.size();++j) {
            fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
            Array<int,6> prsml ;
            prsml[0] = prsm[j][0] ;
            prsml[1] = prsm[j][3] ;
            prsml[2] = prsm[j][4] ;
            prsml[3] = prsm[j][1] ;
            prsml[4] = prsm[j][5] ;
            prsml[5] = prsm[j][2] ;
            fwrite(&prsml[0],sizeof(int),6,OFP) ;
          }

        }
      }
    }
    //write hex
      
    static int hex_walls[6] = { NOT_A_WALL, NOT_A_WALL, NOT_A_WALL,
                                NOT_A_WALL, NOT_A_WALL, NOT_A_WALL };
      
    elem_header = fv_encode_elem_header(FV_HEX_ELEM_ID,
                                        hex_walls) ;

 
    for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumHexs() > 0) { // write out tets
        int tot = volumePartList[i]->getNumHexs() ;
        int start = 0 ;
        int top = tot + volumePartList[i]->getNumHexsIblank() ;
        while(start < top) {
          vector<Array<int,8> > hexs ;
          volumePartList[i]->getHexBlock(hexs,start,block_size) ;
          start += block_size ;
          for(size_t j=0;j<hexs.size();++j) {
            for(int k = 0; k < 8; k++){
              hexs[j][k] += node_offset[i];
            }
          }
          for(size_t j=0;j<hexs.size();++j) {
            fwrite(&elem_header,sizeof(elem_header),1,OFP) ;
            Array<int,8> hlocal ;
            hlocal[0] = hexs[j][0] ;
            hlocal[1] = hexs[j][1] ;
            hlocal[2] = hexs[j][3] ;
            hlocal[3] = hexs[j][2] ;
            hlocal[4] = hexs[j][4] ;
            hlocal[5] = hexs[j][5] ;
            hlocal[6] = hexs[j][7] ;
            hlocal[7] = hexs[j][6] ;
            fwrite(&hlocal[0],sizeof(int),8,OFP) ;
          }
        }
      }
    }
    
    //write general cell
      
    for(size_t i=0;i<volumePartList.size();++i) {
      if(volumePartList[i]->getNumGenc() > 0) { // write out general cells
        vector<int> nfaces, nsides,nodes ;
        volumePartList[i]->getGenCell(nfaces,nsides,nodes) ;
        size_t nnfaces = nfaces.size();
        size_t nnsides = nsides.size();
        size_t nnodes = nodes.size();
        for(size_t j=0;j<nodes.size();++j) {
          nodes[j] += node_offset[i];
        }
     
        
        ibuf[0] = FV_ARB_POLY_ELEMENTS ;
        ibuf[1] = nnfaces ;
        fwrite(ibuf,sizeof(int),2,OFP) ;
        int sc = 0 ;
        int nc = 0 ;
        for(size_t fi=0;fi<nnfaces;++fi) {
          int nf = nfaces[fi] ;
          int nnd = 0 ;
          for(int j=0;j<nf;++j)
            nnd+=nsides[sc+j] ;
          vector<int> tmp(nnd) ;
          for(int j=0;j<nnd;++j)
            tmp[j] = nodes[nc+j] ;
          sort(tmp.begin(),tmp.end()) ;
          tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end()) ;
    
          ibuf[0] = nf ;
          ibuf[1] = tmp.size() ;
          ibuf[2] = -1 ;
          fwrite(ibuf,sizeof(int),3,OFP) ;
          for(int j=0;j<nf;++j) {
            ibuf[0] = A_WALL ;
            fwrite(ibuf,sizeof(int),1,OFP) ;
            ibuf[0] = nsides[sc] ;
            fwrite(ibuf,sizeof(int),1,OFP) ;
            fwrite(&nodes[nc],sizeof(int),nsides[sc],OFP) ;
            nc += nsides[sc] ; // next node starting place
            sc++ ; // next face
            ibuf[0] = 0 ;
            fwrite(ibuf,sizeof(int),1,OFP) ;
          }
        }
          
        if(size_t(nc) != nnodes) {
          cerr << " something wrong!" << endl ;
        }
        if(size_t(sc) != nnsides) {
          cerr << "something worng with sides!" << endl ;
        }
      }
    }
    //finish writing volume elements

   
    set<string>::const_iterator si ;
    // write out nodal scalars values
    // vector val_out will be preset to 0, then each volume part and surfaces part
    // put its val in  
    for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {
      string varname = *si;
      vector<float> val_out(npnts, 0.0);
      if(first_var) {
       
        ibuf[0] = FV_VARIABLES ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
        first_var = false ;
      }
      
      for(size_t i =0;i<volumePartList.size();++i) {
        if(volumePartList[i]->hasNodalScalarVar(varname)) {
          vector<float> val ;
          volumePartList[i]->getNodalScalar(varname,val);
          for( size_t j = 0; j < val.size(); j++) val_out[j+node_offset[i]] = val[j] ;
        }
      }
      
      for(size_t i =0;i<surfacePartList.size();++i) {
        if(surfacePartList[i]->hasNodalScalarVar(varname)) {
          vector<float> val ;
          surfacePartList[i]->getNodalScalar(varname,val);
          for( size_t j = 0; j < val.size(); j++) val_out[j+node_offset[i+volumePartList.size()]] = val[j] ;
        }
      }
       
      fwrite(&val_out[0],sizeof(float),npnts,OFP) ;
    }
    
    //write nodal vector values
    // vector val_out will be preset to vector3d(0.0, 0.0, 0.0), then each volume part and surfaces part
    // put its val in  
   
    for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
      string varname = *si ;
      vector<float> tmp_x(npnts, 0.0) ;
      vector<float> tmp_y(npnts, 0.0) ;
      vector<float> tmp_z(npnts, 0.0) ;
      if(first_var) {
       
        ibuf[0] = FV_VARIABLES ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
        first_var = false ;
      }
        
      for(size_t i =0;i<volumePartList.size();++i) {
        if(volumePartList[i]->hasNodalVectorVar(varname)) {
          vector<vector3d<float> > val ;
          volumePartList[i]->getNodalVector(varname,val) ;
          for(size_t j=0;j<val.size();++j){
            tmp_x[j+node_offset[i]] = val[j].x ;
            tmp_y[j+node_offset[i]] = val[j].y ;
            tmp_z[j+node_offset[i]] = val[j].z ;
          }
        }
      }
      for(size_t i =0;i<surfacePartList.size();++i) {
        if(surfacePartList[i]->hasNodalVectorVar(varname)) {
          vector<vector3d<float> > val ;
          surfacePartList[i]->getNodalVector(varname,val) ;
          for(size_t j=0;j<val.size();++j){
            tmp_x[j+node_offset[i+volumePartList.size()]] = val[j].x ;
            tmp_y[j+node_offset[i+volumePartList.size()]] = val[j].y ;
            tmp_z[j+node_offset[i+volumePartList.size()]] = val[j].z ;
          }
        }
      }
              
      for(size_t i=0;i<npnts;++i) {
        float d = tmp_x[i];
        fwrite(&d,sizeof(float),1,OFP) ;
      }
      for(size_t i=0;i<npnts;++i) {
        float d = tmp_y[i] ;
        fwrite(&d,sizeof(float),1,OFP) ;
      }
      for(size_t i=0;i<npnts;++i) {
        float d = tmp_z[i] ;
        fwrite(&d,sizeof(float),1,OFP) ;
      }
    }
    
   
    // writing element scalar values
    for( si=element_scalars.begin();si!=element_scalars.end();++si){
      string varname = *si ;
      if(first_var) {
       
        ibuf[0] = FV_VARIABLES ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
        first_var = false ;

      }
      if(first_boundary) {
       
        if(general_boundary)
          ibuf[0] = FV_ARB_POLY_BNDRY_VARS ;
        else
          ibuf[0] = FV_BNDRY_VARS ;
    
        fwrite(ibuf,sizeof(int),1,OFP) ;
        first_boundary=false ;
      }

      
      vector<float> valout(elem_ids.size(), 0) ;
      for(size_t i =0;i<surfacePartList.size();++i) {
        if(surfacePartList[i]->hasElementScalarVar(varname)) {
          vector<float> qvals, tvals, gvals ;
          surfacePartList[i]->getElementScalar(varname,qvals,tvals,gvals) ;
          int nqval = qvals.size();
          int ntval = tvals.size();
          int ngval = gvals.size();
          for (int fi=0;fi<nqval;fi++) valout[fi+face_offset[i]] = qvals[fi];
          for (int fi=0;fi<ntval;fi++) valout[fi+nqval+face_offset[i]] = tvals[fi];
          for (int fi=0;fi<ngval;fi++) valout[fi+nqval+ntval+face_offset[i]] = tvals[fi];
        }
      }
     
      fwrite(&valout[0],sizeof(float),valout.size(),OFP) ;
    }  
   
    //writing element vector values
    for(si=element_vectors.begin();si!=element_vectors.end();++si) {
      string varname = *si ;
     
      if(first_var) {
       
        ibuf[0] = FV_VARIABLES ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
        first_var = false ;
      }
      if(first_boundary) {
        
        if(general_boundary)
              ibuf[0] = FV_ARB_POLY_BNDRY_VARS ;
        else
          ibuf[0] = FV_BNDRY_VARS ;
        fwrite(ibuf,sizeof(int),1,OFP) ;
        first_boundary=false ;
      } 
     
      vector<float> xvalout(elem_ids.size(),0.0);  
      vector<float> yvalout(elem_ids.size(),0.0);  
      vector<float> zvalout(elem_ids.size(),0.0);
      for(size_t i =0;i<surfacePartList.size();++i) {
        if(surfacePartList[i]->hasElementVectorVar(varname)) {
          vector<vector3d<float> > qvals, tvals, gvals ;
          surfacePartList[i]->getElementVector(varname,qvals,tvals,gvals) ;
          int nqval = qvals.size();
          int ntval = tvals.size();
          int ngval = gvals.size();
            
          for (int fi=0;fi<nqval;fi++) xvalout[fi+face_offset[i]] = qvals[fi].x;
          for (int fi=0;fi<ntval;fi++) yvalout[fi+nqval+face_offset[i]] = tvals[fi].y;
          for (int fi=0;fi<ngval;fi++) zvalout[fi+nqval+ntval+face_offset[i]] = tvals[fi].z;
          
        }
      }
      
      fwrite(&xvalout[0],sizeof(float),xvalout.size(),OFP) ;
      
      fwrite(&yvalout[0],sizeof(float),yvalout.size(),OFP) ;
     
      fwrite(&zvalout[0],sizeof(float),zvalout.size(),OFP) ;
    }
     
     
     //close grid file
    {
      if(first_boundary) {
        if(general_boundary)
          ibuf[0] = FV_ARB_POLY_BNDRY_VARS ;
        else
          ibuf[0] = FV_BNDRY_VARS ;
        
         fwrite(ibuf,sizeof(int),1,OFP) ;
         first_boundary=false ;
      }
      fclose(OFP) ;
    }
    
   
   
      
    //  output particle files
    // each particle part in a different file  
    
    if(particlePartList.size() > 0) {
      for(size_t pi=0;pi<particlePartList.size();++pi) {
        //write particle scalar values
        vector<string>  particle_scalars_names = particlePartList[pi]->getScalarVars() ;
        vector<vector<float> > particle_scalars;
        for(size_t j = 0; j <  particle_scalars_names.size(); j++){
          vector<float> val ;
          particlePartList[pi]->getParticleScalar( particle_scalars_names[j],val);
          particle_scalars.push_back(val) ;

        }
       
        //write particle vector values
        vector<string>  particle_vectors_names = particlePartList[pi]->getVectorVars() ;
        vector<vector<vector3d<float> > > particle_vectors ;
        for(size_t j = 0; j <  particle_vectors_names.size(); j++){
          vector<vector3d<float> > val;
          particlePartList[pi]->getParticleVector( particle_vectors_names[j],val);
          particle_vectors.push_back(val) ;
        }
     
        //write particle position

      
       
      
        vector<vector3d<float> > pos ;
        particlePartList[pi]->getParticlePositions(pos) ;
        int pts = pos.size();
             
        string particle_filename = dirname + "/" + casename+"_fv_"+iteration+".fvp";
        if(pi > 0 ){
          stringstream ss;//create a stringstream
          ss << pi;//add number to the stream
          string part = ss.str();
          particle_filename = dirname + "/" + casename+"_fv_"+iteration+"p_"+part+".fvp";
        }
#ifdef OLD
        // create an ASCII file handler
        ofstream of(particle_filename.c_str(), std::ios::out) ;
        if(!of.good()) {
          cerr << "can't open file '" << particle_filename
               << "' for writing particle info!" << endl ;
          return ;
        }

        // header information
        of << "FVPARTICLES 2 1" << endl ;
        of << "Tag Names" << endl ;
        of << "0" << endl ;
        of << "Variable Names" << endl ;

        int vsize = particle_scalars_names.size() +
          particle_vectors_names.size() ;

        of << vsize << endl ;

        for(size_t i=0;i<particle_scalars_names.size();++i)
          of << particle_scalars_names[i] << endl ;
        for(size_t i=0;i<particle_vectors_names.size();++i)
          of << particle_vectors_names[i] << endl ;

        // particle data

        // we will have to create "m" pathes with 1 particle per path
        
        for(int p=0;p<pts;++p) {
          of << "1" << endl ;
          // first output position data
          of << pos[p].x << " " << pos[p].y << " " << pos[p].z << " " ;
          // then output variables data
          // NOTE: variable data are already sampled, so we use "p" as index
          for(size_t j=0;j<particle_scalars.size();++j) {
            of << particle_scalars[j][p] << " " ;
          }
          // Note, it seems that fieldview particle format
          // only allows scalars, for vector variables, we
          // compute its magnitude instead
          for(size_t j=0;j<particle_vectors.size();++j) {
            float x = particle_vectors[j][p].x ;
            float y = particle_vectors[j][p].y ;
            float z = particle_vectors[j][p].z ;
            float vm = x*x + y*y + z*z ;
            vm = sqrt(vm) ;

            of << vm << " " ;
          }
          of << endl ;
        }
        of.close() ;
#endif
	FILE *FVP = 0 ;
	FVP = fopen(particle_filename.c_str(), "wb") ;
	if(FVP == NULL) {
	  cerr << "can't open file " << particle_filename << endl ;
	  exit(-1) ;
	}
	int ibuf[3] ;
	char fv[80] ;
	// Write out the magic number and the fieldview version info 
	ibuf[0] = FV_MAGIC ;
	fwrite(ibuf, sizeof(int), 1, FVP) ;
	memset(fv,'\0',80) ;
	snprintf(fv,80, "FVPARTICLES") ;
	fwrite(&fv, sizeof(char), 80, FVP) ;
	ibuf[0] = MAJOR_VERSION ;
	ibuf[1] = MINOR_VERSION ;
	ibuf[2] = 0 ; //for FV reserved use
	fwrite(ibuf, sizeof(int), 3,FVP) ;

	int vsize = particle_scalars_names.size() +
	  particle_vectors_names.size() ;
	ibuf[0] = vsize ;
	fwrite(ibuf, sizeof(int), 1, FVP) ;
	
	for(size_t i=0;i<particle_scalars_names.size();++i) {
	  memset(fv,'\0',80) ;
	  snprintf(fv,80, "%s",particle_scalars_names[i].c_str()) ;
	  fwrite(&fv, sizeof(char), 80, FVP) ;
	}

	for(size_t i=0;i<particle_vectors_names.size();++i) {
	  memset(fv,'\0',80) ;
	  snprintf(fv,80, "%s",particle_vectors_names[i].c_str()) ;
	  fwrite(&fv, sizeof(char), 80, FVP) ;
	}
	
	// particle data
	
	ibuf[0] = pts ; //number of particles
	fwrite(ibuf,sizeof(int),1,FVP) ;
	
	for(int p=0;p<pts;++p) {
	  float x = pos[p].x ;
	  float y = pos[p].y ;
	  float z = pos[p].z ;
	  fwrite(&x,sizeof(float),1,FVP) ;
	  fwrite(&y,sizeof(float),1,FVP) ;
	  fwrite(&z,sizeof(float),1,FVP) ;
	}

	// NOTE: variable data are already sampled, so we use "p" as index
	for(size_t i=0;i<particle_scalars.size();++i) {
	  for(int p=0;p<pts;++p) {
	    float var = particle_scalars[i][p] ;
	    fwrite(&var,sizeof(float),1,FVP) ;
	  }
	}

	// Note, it seems that fieldview particle format
	// only allows scalars, for vector variables, we
	// compute its magnitude instead
	for(size_t i=0;i<particle_vectors.size();++i) {
	  for(int p=0;p<pts;++p) {
	    float x = particle_vectors[i][p].x ;
	    float y = particle_vectors[i][p].y ;
	    float z = particle_vectors[i][p].z ;
	    float vm = x*x + y*y + z*z ;
	    vm = sqrt(vm) ;
	    fwrite(&vm,sizeof(float),1,FVP) ;
	  }
	}
      }
    }
       
  
}

// ... the end ...

