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
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>

#ifdef USE_NATIVE_TECPLOT
#include "TECIO.h"

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

//#define CREATE_ASCII
//#define CREATE_CLASSIC_PLT

using namespace std ;

namespace {
/**********************************************
 *            TECXXX WRAPPER FUNCTIONS        *
 **********************************************/
bool TECINI_internal(char const* problem_name,
                     char const* variables_name,
                     char const* tec_name,
                     char const* dot,
                     int*        Debug,
                     int*        v_is_double) {

bool isOk = true;
#if defined CREATE_ASCII
ostringstream oss ;
oss << tec_name ;
string filename = oss.str() ;
cout << "Initializing file " << filename << endl ;
ofstream ofile(filename.c_str(),ios::out) ;
ofile << "variables = " << variables_name << endl ;
ofile.close() ;

#else
  #if defined CREATE_CLASSIC_PLT
  INTEGER4 fileFormat = 0;
  #else
  INTEGER4 fileFormat = 1;
  #endif

  INTEGER4 fileType = 0;
  INTEGER4 debug = static_cast<INTEGER4>(*Debug);
  INTEGER4 VIsDouble = static_cast<INTEGER4>(*v_is_double);
  isOk = TECINI142(problem_name, variables_name, tec_name, dot, &fileFormat, &fileType, &debug, &VIsDouble) == 0;
#endif

  return(isOk) ;
}

bool TECZNE_internal(char const* zone_name, 
                     int*        npnts, 
                     int*        ncells,
                     char const* tec_name,
                     int const*  zone_type, 
                     double*     solution_time,
                     int*        strandid) {

  bool isOk = true;
#if defined CREATE_ASCII
  string filename = tec_name ;
  ofstream ofile(filename.c_str(),ios::app) ;
 
  ofile << "ZONE T = \"" << zone_name << '"'
        << ", DATAPACKING=BLOCK" 
        << ", n = " << *npnts
        << ", e = " << *ncells
        << ", ZONETYPE = ";

  if      (*zone_type == 1) ofile << "ORDERED\n";
  else if (*zone_type == 2) ofile << "FETRIANGLE\n";
  else if (*zone_type == 3) ofile << "FEQUADRILATERAL\n";
  else if (*zone_type == 4) ofile << "FETETRA\n";
  else if (*zone_type == 5) ofile << "FEBRICK\n";
  else if (*zone_type == 6) ofile << "FELINESEG\n";
  else if (*zone_type == 7) ofile << "FEPOLYGON\n";
  else if (*zone_type == 8) ofile << "FEPOLYHEDRON\n";
  else {
      fprintf(stderr,"Err: Invalid zone type\n");
      isOk = true;
  }

  ofile.close() ;

#else
    INTEGER4 zoneType = static_cast<INTEGER4>(*zone_type);
    INTEGER4 numberPoints = static_cast<INTEGER4>(*npnts);
    INTEGER4 numberElements = static_cast<INTEGER4>(*ncells);
    INTEGER4 nFaces = 0;
    INTEGER4 ICellMax = 0;
    INTEGER4 JCellMax = 0;
    INTEGER4 KCellMax = 0;
    INTEGER4 strandID = static_cast<INTEGER4>(*strandid);
    INTEGER4 ParentZone = 0;
    INTEGER4 IsBlock = 1;
    INTEGER4 NumFaceConnections = 0;
    INTEGER4 FaceNeighborMode = 0; // There are no face connections
    INTEGER4 TotalNumFaceNodes = 0;
    INTEGER4 NumConnectedBoundaryFaces = 0;
    INTEGER4 TotalNumBoundaryConnections = 0;
    INTEGER4 ShareConnectivityFromZone = 0;

    isOk = TECZNE142(zone_name, &zoneType, &numberPoints, &numberElements, &nFaces, &ICellMax, &JCellMax, &KCellMax, solution_time, &strandID, &ParentZone, &IsBlock, &NumFaceConnections, &FaceNeighborMode, &TotalNumFaceNodes, &NumConnectedBoundaryFaces,
                        &TotalNumBoundaryConnections, NULL, NULL, NULL, &ShareConnectivityFromZone) == 0; 
#endif

    return(isOk) ;
}

bool TECDAT_internal(int*        npnts, 
                     float*      vals, 
                     char const* tec_name) {

  bool isOk = true;
#if defined CREATE_ASCII
  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  ofstream ofile(filename.c_str(),ios::app) ;
  ofile.precision(16) ;

  for(int i = 0; i < *npnts; i++)
    ofile << vals[i] << endl ;

  ofile.close() ;

#else
    INTEGER4 N = static_cast<INTEGER4>(*npnts);
    INTEGER4 IsDouble = 0;
    isOk = TECDAT(&N, vals, &IsDouble) == 0;
#endif

  return(isOk) ;
}

bool TECNOD_internal(int*        nm, 
                     int*        ncells, 
                     int         nodesPerCell,
                     char const* tec_name) {

  bool isOk = true;
#if defined CREATE_ASCII
  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  ofstream ofile(filename.c_str(),ios::app) ;
  for(int i = 0; i < *ncells; i++) {
      int base = i*nodesPerCell;
      for (int j = 0; j < nodesPerCell; ++j)
          ofile << nm[base+j] << " ";
      ofile << endl;
  }

  ofile.close() ;

#else 
  int64_t totalNumValues = static_cast<int64_t>(*ncells)*nodesPerCell;
  int64_t curOffset = 0;
  while (isOk && curOffset < totalNumValues) {
      int64_t numValuesLeft = totalNumValues - curOffset;
      int32_t nextChunkSize = (std::numeric_limits<int32_t>::max()-1);

      if (nextChunkSize > numValuesLeft)
          nextChunkSize = static_cast<int32_t>(numValuesLeft);

      isOk = TECNODE142(&nextChunkSize,&nm[curOffset]) == 0;

      curOffset += nextChunkSize;
  }
#endif
  return(isOk) ;
}



bool TECEND_internal() {
#if defined CREATE_ASCII
  return true;
#else
  return(TECEND142() == 0);
#endif
}
}


/**********************************************
 *              ZONE WRITERS                  *
 **********************************************/
namespace {
bool writeOutSurfaceZone(surfacePartP surfacePart,
                         set<string>& scalar_varnames,
                         set<string>& vector_varnames,
                         std::string& boundary_name,
                         std::string& filename,
                         double       solutionTime,
                         int          strandid) {
  bool isOk = true;
  int npnts = surfacePart->getNumNodes();

  int ntrias = surfacePart->getNumTrias();
  int nquads = surfacePart->getNumQuads();
  int nelm = ntrias + nquads;
  int zoneType = 3;

  isOk = TECZNE_internal(boundary_name.c_str(), &npnts, &nelm, filename.c_str(), &zoneType, &solutionTime, &strandid);

  if (isOk) {
      //create mesh positions
      vector<float> pos_x(npnts),pos_y(npnts),pos_z(npnts) ;

      vector<vector3d<float> > part_pos;
      surfacePart->getPos(part_pos) ;
      for(size_t j = 0; j<part_pos.size(); j++) {
          pos_x[j] = part_pos[j].x;
          pos_y[j] = part_pos[j].y;
          pos_z[j] = part_pos[j].z; 
      }

      isOk = TECDAT_internal(&npnts,&pos_x[0],filename.c_str()) &&
             TECDAT_internal(&npnts,&pos_y[0],filename.c_str()) &&
             TECDAT_internal(&npnts,&pos_z[0],filename.c_str()) ;
   }

  if (isOk) {
      for (set<string>::const_iterator si = scalar_varnames.begin(); si != scalar_varnames.end(); ++si) {
          string varname = *si;
          vector<float> val_out(npnts, 0.0);
          if (surfacePart->hasNodalScalarVar(varname)) {
              vector<float> val;
              surfacePart->getNodalScalar(varname, val);
              isOk = TECDAT_internal(&npnts, &val[0], filename.c_str());
          }
          else {
              isOk = TECDAT_internal(&npnts, &val_out[0], filename.c_str());
          }
      }
  }
  // write nodal vectors
  if (isOk) {
      //write out nodal vector
      for (set<string>::const_iterator si = vector_varnames.begin(); isOk && si != vector_varnames.end(); ++si) {
          string varname = *si;
          vector<float> xtmp(npnts, 0.0);
          vector<float> ytmp(npnts, 0.0);
          vector<float> ztmp(npnts, 0.0);
          if (surfacePart->hasNodalVectorVar(varname)) {
              vector<vector3d<float> > val;
              surfacePart->getNodalVector(varname, val);
              for (size_t j = 0; j < val.size(); ++j) {
                  xtmp[j] = val[j].x;
                  ytmp[j] = val[j].y;
                  ztmp[j] = val[j].z;
              }
          }
          isOk = TECDAT_internal(&npnts, &xtmp[0], filename.c_str()) &&
                 TECDAT_internal(&npnts, &ytmp[0], filename.c_str()) &&
                 TECDAT_internal(&npnts, &ztmp[0], filename.c_str()) ;
      }
  }

   // And finally, write out the nodemap.

   if (isOk) {
       vector<Array<int,4> > ordinary_faces ;
       string name = surfacePart->getPartName();
       //create boundary part
       string  boundary_name = name;
       ordinary_faces.clear();

       //write out quads
       int nquads = surfacePart->getNumQuads();
       if (nquads > 0) {
           vector<Array<int, 4> > quads;
           surfacePart->getQuads(quads);
           for (int j = 0; j < nquads; ++j) {
               Array<int, 4> a;
               a[0] = quads[j][0];
               a[1] = quads[j][1];
               a[2] = quads[j][2];
               a[3] = quads[j][3];

               ordinary_faces.push_back(a);
           }
       }


       //write out trias
       int ntrias = surfacePart->getNumTrias();
       if (ntrias > 0) {
           vector<Array<int, 3> > trias;
           surfacePart->getTrias(trias);
           for (int j = 0; j < ntrias; ++j) {
               Array<int, 4> a;
               a[0] = trias[j][0];
               a[1] = trias[j][1];
               a[2] = trias[j][2];
               a[3] = trias[j][2];
               ordinary_faces.push_back(a);
           }
       }

       isOk = TECNOD_internal(&ordinary_faces[0][0],&nelm,4,filename.c_str());
   }

  if (!isOk) {
      fprintf(stdout,"Error generating tecplot datafile\n");
  }
  return isOk;
}


//TODO (byron) H 2016/07/12: Could easily convert this to write out one volume "part" at a time instead.
//                           Would actually simply things.
bool writeOutCombinedVolumeZone(vector<volumePartP> volumePartList,
                                set<string>&        scalar_varnames,
                                set<string>&        vector_varnames,
                                vector<int>&        volume_part_offset,
                                std::string&        filename,
                                double             solutionTime,
                                int                 strandid) {
    bool isOk = true;
    int npnts = 0;
    int ntets = 0, nprsm = 0, npyrm = 0, nhexs = 0;

    for(size_t i=0;i<volumePartList.size();++i) {
        npnts += volumePartList[i]->getNumNodes();
        ntets += volumePartList[i]->getNumTets();
        nprsm += volumePartList[i]->getNumPrsm();
        nhexs += volumePartList[i]->getNumHexs();
        npyrm += volumePartList[i]->getNumPyrm();
    }

    int nelm = ntets + nprsm + nhexs + npyrm;
    int zoneType = 5;

    isOk = TECZNE_internal("VOLUME_MESH", &npnts, &nelm, filename.c_str(), &zoneType, &solutionTime, &strandid);
    // write out mesh
    if (isOk) {
        vector<float> pos_x(npnts),pos_y(npnts),pos_z(npnts) ;
        {
            size_t tot = 0;
            for(size_t i=0;i<volumePartList.size();++i) {
                vector<vector3d<float> > part_pos;
                volumePartList[i]->getPos(part_pos) ;
                for(size_t j = 0; j<part_pos.size(); j++) {
                    pos_x[j+tot] = part_pos[j].x;
                    pos_y[j+tot] = part_pos[j].y;
                    pos_z[j+tot] = part_pos[j].z; 
                }
                tot += part_pos.size();
            }

            isOk = TECDAT_internal(&npnts,&pos_x[0],filename.c_str()) &&
                   TECDAT_internal(&npnts,&pos_y[0],filename.c_str()) &&
                   TECDAT_internal(&npnts,&pos_z[0],filename.c_str()) ;
        }
    }

    // write out nodal scalars
    if (isOk) {
        for(set<string>::const_iterator si=scalar_varnames.begin();si!=scalar_varnames.end();++si) {                                         
            string varname = *si;
            vector<float> val_out(npnts, 0.0);
            for(size_t i =0;isOk && i<volumePartList.size();++i) {
                if(volumePartList[i]->hasNodalScalarVar(varname)) {
                    vector<float> val ;
                    volumePartList[i]->getNodalScalar(varname,val);
                    isOk = TECDAT_internal(&npnts,&val[0],filename.c_str());
                }
                else {
                    isOk = TECDAT_internal(&npnts,&val_out[0],filename.c_str());
                }

            }
        }
    }

    if (isOk) {
        //write out nodal vector
        for(set<string>::const_iterator si=vector_varnames.begin();isOk && si!=vector_varnames.end();++si) {
            string varname = *si ;
            vector<float> xtmp(npnts, 0.0) ;
            vector<float> ytmp(npnts, 0.0) ;
            vector<float> ztmp(npnts, 0.0) ;
            for(size_t i =0;isOk && i<volumePartList.size();++i) {
                if(volumePartList[i]->hasNodalVectorVar(varname)) {
                    vector<vector3d<float> > val ;
                    volumePartList[i]->getNodalVector(varname,val) ;
                    for(size_t j=0;j<val.size();++j) {
                        xtmp[j] = val[j].x ;
                        ytmp[j] = val[j].y ;
                        ztmp[j] = val[j].z ;
                    }
                }
            }

            isOk = TECDAT_internal(&npnts,&xtmp[0],filename.c_str()) &&
                   TECDAT_internal(&npnts,&ytmp[0],filename.c_str()) &&
                   TECDAT_internal(&npnts,&ztmp[0],filename.c_str()) ;
        }
    }

    vector<Array<int, 8> > bricks ;
    // Write out the nodemap
    if (isOk) {
        // collect tets 
        const int block_size=65536 ; // Size of blocking factor
        for(size_t i=0;i<volumePartList.size();++i) {
            if(volumePartList[i]->getNumTets() > 0) { 
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
                                tets[j][k] += volume_part_offset[i];
                            }
                        }
                        for(size_t j=0;j<tets.size();++j) {
                            Array<int,8> brick ;
                            brick[0] = tets[j][0] ;
                            brick[1] = tets[j][1] ;
                            brick[2] = tets[j][2] ;
                            brick[3] = brick[2] ;
                            brick[4] = tets[j][3] ;
                            brick[5] = brick[4] ;
                            brick[6] = brick[4] ;
                            brick[7] = brick[4] ;
                            bricks.push_back(brick) ;
                        }    
                    }
                }
            }
        }

        //collect pyrm
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
                            pyrm[j][k] += volume_part_offset[i];
                        }
                    }

                    for(size_t j=0;j<pyrm.size();++j) {
                        Array<int,8> brick ;
                        brick[0] = pyrm[j][0] ;
                        brick[1] = pyrm[j][1] ;
                        brick[2] = pyrm[j][2] ;
                        brick[3] = pyrm[j][3] ;
                        brick[4] = pyrm[j][4] ;
                        brick[5] = brick[4] ;
                        brick[6] = brick[4] ;
                        brick[7] = brick[4] ;
                        bricks.push_back(brick) ;
                    }
                }
            }
        }  


        // collect prsm
        for(size_t i=0;i<volumePartList.size();++i) { 
            if(volumePartList[i]->getNumPrsm() > 0) {
                int tot = volumePartList[i]->getNumPrsm() ;
                int start = 0 ;
                int top = tot + volumePartList[i]->getNumPrsmIblank() ;
                while(start < top) {
                    vector<Array<int,6> > prsm ;
                    volumePartList[i]->getPrsmBlock(prsm,start,block_size) ;
                    start += block_size ;
                    for(size_t j=0;j<prsm.size();++j) {
                        for(int k = 0; k < 6; k++){
                            prsm[j][k] += volume_part_offset[i];
                        }
                    }

                    for(size_t j=0;j<prsm.size();++j) {
                        Array<int,8> brick ;
                        brick[0] = prsm[j][0] ;
                        brick[1] = prsm[j][1] ;
                        brick[2] = prsm[j][2] ;
                        brick[3] = prsm[j][2] ;
                        brick[4] = prsm[j][3] ;
                        brick[5] = prsm[j][4] ;
                        brick[6] = prsm[j][5] ;
                        brick[7] = prsm[j][5] ;
                        bricks.push_back(brick) ;
                    }
                }
            }
        }    

        // collect hexes
        for(size_t i=0;i<volumePartList.size();++i) {
            if(volumePartList[i]->getNumHexs() > 0) {
                int tot = volumePartList[i]->getNumHexs() ;
                int start = 0 ;
                int top = tot + volumePartList[i]->getNumHexsIblank() ;
                while(start < top) {
                    vector<Array<int,8> > hexs ;
                    volumePartList[i]->getHexBlock(hexs,start,block_size) ;
                    start += block_size ;
                    for(size_t j=0;j<hexs.size();++j) {
                        for(int k = 0; k < 8; k++){
                            hexs[j][k] += volume_part_offset[i];
                        }
                    }
                    for(size_t j=0;j<(size_t)hexs.size();++j) {
                        Array<int,8> brick ;
                        for(int k=0;k<8;++k)
                            brick[k] = hexs[j][k] ;
                        bricks.push_back(brick) ;
                    }
                }
            }
        }
    }
    //Write volume elements
    if (isOk)
    {
        int nelem = (int)bricks.size();
        isOk = TECNOD_internal(&bricks[0][0],&nelem,8,filename.c_str()) ;
    }

    return isOk;
}

}





bool tecplotPartConverter::processesVolumeElements() const {
  return true ;
}
bool tecplotPartConverter::processesSurfaceElements() const {
  return true ;
}
bool tecplotPartConverter::processesParticleElements() const {
  return false ;
}

void tecplotPartConverter::exportPostProcessorFiles(string casename, string iteration) const {
  
  string filename ;
 
  bool isOk = true;

  char *EndScan = NULL;
  double solutionTime = strtod(iteration.c_str(), &EndScan);
  if (*EndScan == '\0')
      isOk = true;
  else
  {
      isOk = false;
      cerr << "Iteration is not a number." << endl; 
  }


  if(volumePartList.size() == 0 ||
     surfacePartList.size() == 0) {
     cerr<<"ERROR: empty volume or surface part list!" << endl;
     isOk = false;
  }

  set<string> scalar_varnames ;
  set<string> vector_varnames ;
  if (isOk) {
      for(size_t i=0;i<volumePartList.size();++i) {
          vector<string> nscalars = volumePartList[i]->getNodalScalarVars() ;
          for(size_t j=0;j<nscalars.size();++j) 
              scalar_varnames.insert(nscalars[j]) ;
          vector<string> nvectors = volumePartList[i]->getNodalVectorVars() ;
          for(size_t j=0;j<nvectors.size();++j) 
              vector_varnames.insert(nvectors[j]) ;
      }
  }


  // prepare filename for file and variable string
  if (isOk) {
    filename = "tec_"+casename+"_"+iteration;
    // szplt automatically appends the .szplt extension.  Others do not.
#if defined CREATE_CLASSIC_PLT
    filename += ".plt" ;
#elif defined CREATE_ASCII
    filename += ".dat" ;
#endif
    string varstring = "x,y,z" ;

    for(set<string>::const_iterator si=scalar_varnames.begin();si!=scalar_varnames.end();++si) {                                         
      string varname = *si ;
      varstring += ","+varname;
    }
    for(set<string>::const_iterator si=vector_varnames.begin();si!=vector_varnames.end();++si) {
      string varname = *si ;
      varstring += ","+varname+".x" ;
      varstring += ","+varname+".y" ;
      varstring += ","+varname+".z" ;
    } 
  
    int Debug = 0 ;
    int VIsDouble = 0 ;
    isOk = TECINI_internal(casename.c_str(),
                           varstring.c_str(),
                           filename.c_str(),
                           ".",
                           &Debug,
                           &VIsDouble) ;
  }

  /*
   * Collect offests for volume parts.  Remove this if we elect to write out
   * separate zones for each volume part someday.
   */
  vector<int> volume_part_offset(volumePartList.size());
  volume_part_offset[0] = 0;
  for(size_t i=1;i<volumePartList.size();++i) {
      vector<vector3d<float> > part_pos;
      volumePartList[i]->getPos(part_pos) ;
      volume_part_offset[i] = volume_part_offset[i-1] + part_pos.size();
  }
  int strandID = 1;
  if (isOk)
      // use the iteration as the solution time
      isOk = writeOutCombinedVolumeZone(volumePartList,
                                        scalar_varnames,
                                        vector_varnames,
                                        volume_part_offset,
                                        filename,
                                        solutionTime,
                                        strandID);
      strandID++ ;
  if (isOk) {
      //write boundary elements
      for (size_t i = 0; isOk && i < surfacePartList.size(); ++i) {
          string boundary_name = surfacePartList[i]->getPartName();
          isOk = writeOutSurfaceZone(surfacePartList[i],
                                     scalar_varnames,
                                     vector_varnames,
                                     boundary_name,
                                     filename,
                                     solutionTime,
                                     strandID);
          strandID++;
      }

      if (isOk)
          isOk = TECEND_internal();
  }
  if (!isOk) {
      fprintf(stdout,"Error generating tecplot datafile\n");
  }
}

#else 
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
int TECINI(const char *problem_name,const char *variables_name,
           const char *tec_name,const char *dot,int *Debug,int *VIsDouble) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  cout << "Initializing file " << filename << endl ;
 
  ofstream ofile(filename.c_str(),ios::out) ;
  ofile << "variables = " << variables_name << endl ;
  ofile.close() ;

  return(0) ;
}

int TECZNE(const char *zone_name, int *npnts, int *ncells,const char *tec_name,
           const char *elem_type) {

  string filename = tec_name ;
  ofstream ofile(filename.c_str(),ios::app) ;
 
  ofile << "ZONE T = \"" << zone_name << '"'
        << ", DATAPACKING=BLOCK" 
        << ", n = " << *npnts
        << ", e = " << *ncells
        << ", ZONETYPE = " << elem_type << endl ;
  ofile.close() ;

  return(0) ;
}

int TECDAT(int *npnts, float *vals, const char *tec_name) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  ofstream ofile(filename.c_str(),ios::app) ;
  ofile.precision(16) ;

  for(int i = 0; i < *npnts; i++)
    ofile << vals[i] << endl ;

  ofile.close() ;

  return(0) ;
}

int TECNOD(int *nm, int *ncells, const char *tec_name) {

  ostringstream oss ;

  oss << tec_name ;
  string filename = oss.str() ;
  ofstream ofile(filename.c_str(),ios::app) ;
  for(int i = 0; i < *ncells; i++) {
    int j = i*8 ;
    ofile << nm[j+0] << " " << nm[j+1] << " "
          << nm[j+2] << " " << nm[j+3] << " "
          << nm[j+4] << " " << nm[j+5] << " "
          << nm[j+6] << " " << nm[j+7] << endl ;
  }

  ofile.close() ;

  return(0) ;
}

int TECEND() {
  return(0)  ;
}


bool tecplotPartConverter::processesVolumeElements() const {
  return true ;
}
bool tecplotPartConverter::processesSurfaceElements() const {
  return true ;
}
bool tecplotPartConverter::processesParticleElements() const {
  return false ;
}

void tecplotPartConverter::exportPostProcessorFiles(string casename, string iteration) const {
  
  string filename ;
  size_t npnts = 0 ;
  size_t ntets = 0, nprsm = 0, npyrm = 0, nhexs = 0 ; //, ngen = 0 ;
  int nvars = 0 ;
  vector<Array<int, 8> > bricks ;
 
  vector<Array<int,4> > ordinary_faces ;
  vector<int> node_ids ;
  vector<int> elem_ids ;
 
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
    //    ngen += volumePartList[i]->getNumGenc();
  }
  for(size_t i =0;i<surfacePartList.size();++i) {
    npnts += surfacePartList[i]->getNumNodes() ;
  }

  //open file and write header
  {
  
    nvars = 3 ;

    filename = "tec_"+casename+"_"+iteration+".dat" ;
    int ncells = ntets+nprsm+npyrm+nhexs ;
    string varstring = "\"x\", \"y\", \"z\"" ;
    set<string>::const_iterator si ;
    for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {                                         
      string varname = *si ;
      varstring += ", \"" + varname + '"' ;
      nvars++ ;
    }
    for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
      string varname = *si ;
      varstring += ", \"" + varname +".x\"" ;
      varstring += ", \"" + varname +".y\"" ;
      varstring += ", \"" + varname  +".z\"" ;
      nvars += 3;
    } 
  
    int Debug = 1 ;
    int VIsDouble = 0 ;
    TECINI(casename.c_str(),varstring.c_str(),filename.c_str(),
           ".",&Debug,&VIsDouble) ;
    int npt = npnts ;
    TECZNE("VOLUME_MESH",&npt,&ncells,filename.c_str(),"FEBRICK") ;
  
  }
  
  //create mesh positions
  int nparts = volumePartList.size() + surfacePartList.size();
  vector<int> node_offset(nparts);
  vector<float> pos_x(npnts),pos_y(npnts),pos_z(npnts) ;
  
  {
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
    int npt = npnts ;
    TECDAT(&npt,&pos_x[0],filename.c_str()) ;
    TECDAT(&npt,&pos_y[0],filename.c_str()) ;
    TECDAT(&npt,&pos_z[0],filename.c_str()) ;
  }
 
  set<string>::const_iterator si ;
  // write out nodal scalars
  for(si=nodal_scalars.begin();si!=nodal_scalars.end();++si) {                                         
    string varname = *si;
    vector<float> val_out(npnts, 0.0);
    for(size_t i =0;i<volumePartList.size();++i) {
      if(volumePartList[i]->hasNodalScalarVar(varname)) {
        vector<float> val ;
        volumePartList[i]->getNodalScalar(varname,val);
        for(long long unsigned int j = 0; j < (long long unsigned int)val.size(); j++) val_out[j+node_offset[i]] = val[j] ;
      }
    }
                                              
    int npt = npnts ;
    TECDAT(&npt,&val_out[0],filename.c_str()) ;
  }

  //write out nodal vector
  for(si=nodal_vectors.begin();si!=nodal_vectors.end();++si) {
    string varname = *si ;
    vector<float> tmp(npnts, 0.0) ;
    for(size_t i =0;i<volumePartList.size();++i) {
      if(volumePartList[i]->hasNodalVectorVar(varname)) {
        vector<vector3d<float> > val ;
        volumePartList[i]->getNodalVector(varname,val) ;
        for(size_t j=0;j<val.size();++j)
          tmp[j+node_offset[i]] = val[j].x ;
      }
    }
    int npt = npnts ;
    TECDAT(&npt,&tmp[0],filename.c_str()) ;
    
    for(size_t i =0;i<volumePartList.size();++i) {
      if(volumePartList[i]->hasNodalVectorVar(varname)) {
        vector<vector3d<float> > val ;
        volumePartList[i]->getNodalVector(varname,val) ;
        for(size_t j=0;j<val.size();++j)
          tmp[j+node_offset[i]] = val[j].y ;
      }
    }
    TECDAT(&npt,&tmp[0],filename.c_str()) ;

    for(size_t i =0;i<volumePartList.size();++i) {
      if(volumePartList[i]->hasNodalVectorVar(varname)) {
        vector<vector3d<float> > val ;
        volumePartList[i]->getNodalVector(varname,val) ;
        for(size_t j=0;j<val.size();++j)
          tmp[j+node_offset[i]] = val[j].z ;
        TECDAT(&npt,&tmp[0],filename.c_str()) ;
      }  
    }
  }
  //write out volume elements
  // write out tets 
  const int block_size=65536 ; // Size of blocking factor
  for(size_t i=0;i<volumePartList.size();++i) {
    if(volumePartList[i]->getNumTets() > 0) { 
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
            Array<int,8> brick ;
            brick[0] = tets[j][0] ;
            brick[1] = tets[j][1] ;
            brick[2] = tets[j][2] ;
            brick[3] = brick[2] ;
            brick[4] = tets[j][3] ;
            brick[5] = brick[4] ;
            brick[6] = brick[4] ;
            brick[7] = brick[4] ;
            bricks.push_back(brick) ;
          }    
        }
      }
    }
  }

  //write out pyrm
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
          Array<int,8> brick ;
          brick[0] = pyrm[j][0] ;
          brick[1] = pyrm[j][1] ;
          brick[2] = pyrm[j][2] ;
          brick[3] = pyrm[j][3] ;
          brick[4] = pyrm[j][4] ;
          brick[5] = brick[4] ;
          brick[6] = brick[4] ;
          brick[7] = brick[4] ;
          bricks.push_back(brick) ;
        }
      }
    }
  }  
  

 
  // write out prsm
  for(size_t i=0;i<volumePartList.size();++i) { 
    if(volumePartList[i]->getNumPrsm() > 0) {
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
          Array<int,8> brick ;
          brick[0] = prsm[j][0] ;
          brick[1] = prsm[j][1] ;
          brick[2] = prsm[j][2] ;
          brick[3] = prsm[j][2] ;
          brick[4] = prsm[j][3] ;
          brick[5] = prsm[j][4] ;
          brick[6] = prsm[j][5] ;
          brick[7] = prsm[j][5] ;
          bricks.push_back(brick) ;
        }
      }
    }
  }    
  
  // write out tets
  for(size_t i=0;i<volumePartList.size();++i) {
    if(volumePartList[i]->getNumHexs() > 0) {
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
        for(size_t j=0;j<nhexs;++j) {
          Array<int,8> brick ;
          for(int k=0;k<8;++k)
            brick[k] = hexs[j][k] ;
          bricks.push_back(brick) ;
        }
      }
    }
  }
  //close volume elements
  {
    int ncells = bricks.size() ;
    TECNOD(&bricks[0][0],&ncells, filename.c_str()) ;
  }
    
  //write boundary elements
  for(size_t i =0;i<surfacePartList.size();++i) {
    string name = surfacePartList[i]->getPartName(); 
    //create boundary part
    string  boundary_name = name ;
    int part_num = volumePartList.size()+i;
    int npt = surfacePartList[i]->getNumNodes() ;
    ordinary_faces.clear() ;
    node_ids.clear() ;
    for(int j=1;j<=npt;++j)
      node_ids.push_back(node_offset[part_num]+j) ;
      
  
    int face_id = 0;
    //write out quads
    int nquads = surfacePartList[i]->getNumQuads() ;
    if(nquads > 0) {
      vector<Array<int,4> > quads ;
      surfacePartList[i]->getQuads(quads) ; 
      for(int j=0;j<nquads;++j) {
        Array<int,4> a ;
        a[0] = quads[j][0];
        a[1] = quads[j][1];
        a[2] = quads[j][2];
        a[3] = quads[j][3];
      
        ordinary_faces.push_back(a) ;
        elem_ids.push_back(face_id++) ;
      }
    }
 

    //write out trias
    int ntrias = surfacePartList[i]->getNumTrias() ;
    if(ntrias > 0) {
      vector<Array<int,3> > trias ;
      surfacePartList[i]->getTrias(trias) ; 
      for(int j=0;j<ntrias;++j) {
        Array<int,4> a ;
        a[0] = trias[j][0] ;
        a[1] = trias[j][1];
        a[2] = trias[j][2];
        a[3] = trias[j][2];
        ordinary_faces.push_back(a) ;
        elem_ids.push_back(face_id++) ;
      }
    }


    //write out general faces

    ofstream ofile(filename.c_str(),ios::app) ;
       
    int nelm = ordinary_faces.size() ;
    ofile << "ZONE T = \"" << boundary_name << '"'
          << ", N = " << npnts
          << ", E = " << nelm
          << ", ZONETYPE=FEQUADRILATERAL"
          << ",VARSHARELIST = ([1" ;
    for(int j=1;j<nvars;++j)
      ofile << ',' << j+1 ;
    ofile << "]=1)" << endl ;
       
    for(int j=0;j<nelm;++j)
      ofile << node_ids[ordinary_faces[j][0]-1] << ' '
            << node_ids[ordinary_faces[j][1]-1] << ' '
            << node_ids[ordinary_faces[j][2]-1] << ' '
            << node_ids[ordinary_faces[j][3]-1] << endl ;
    
    
    //close boundary part
    ordinary_faces.clear() ;
    node_ids.clear() ;
    boundary_name = "" ;
  }
 
}

#endif
