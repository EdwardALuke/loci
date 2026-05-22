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
#include <iostream>
#include <fstream>
#include "vogtools.h"
#include "Tools/debug.h"

#include <map>

#include <sstream>
using std::istringstream ;
#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
using std::cout ;
using std::endl ;
using std::cerr ;
using std::ifstream ;
using std::ios ;
#include <ctype.h>

bool checkToken(string s1, string s2) {
  int s1l = s1.length() ;
  int s2l = s2.length() ;
  if(s1l < s2l)
    return false ;
  for(int i=0;i<s2l;++i)
    if(toupper(s1[i]) != toupper(s2[i]))
      return false ;
  return true ;
}
struct quadFace {
  Array<int,4> nodes ;
  int cell ;
  bool left ;
} ;

struct triaFace {
  Array<int,3> nodes ;
  int cell ;
  bool left ;
} ;

inline bool quadCompare(const quadFace &f1, const quadFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] == f2.nodes[2] && f1.nodes[3] < f2.nodes[3]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] == f2.nodes[2] && f1.nodes[3] == f2.nodes[3] &&
	   f1.cell > f2.cell)) ;
}


inline bool triaCompare(const triaFace &f1, const triaFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] == f2.nodes[2] && f1.cell > f2.cell)) ;
}


inline bool quadEqual(const quadFace &f1, const quadFace &f2) {
  return ((f1.nodes[0] == f2.nodes[0]) &&
          (f1.nodes[1] == f2.nodes[1]) &&
          (f1.nodes[2] == f2.nodes[2]) &&
          (f1.nodes[3] == f2.nodes[3]));
}

inline bool triaEqual(const triaFace &f1, const triaFace &f2) {
  return ((f1.nodes[0] == f2.nodes[0]) &&
          (f1.nodes[1] == f2.nodes[1]) &&
          (f1.nodes[2] == f2.nodes[2]));

}

void getnewline(std::istream &s) {
  while(!s.eof() && s.peek() != '\n')
    s.get() ;
  s.get() ;
}

void swapb(char *p, int size) {
  char *pe = p+(size-1) ;
  while(p<pe) {
    std::swap(*p,*pe) ;
    p++ ;
    pe--;
  }
}

#define READVAR(VAL)  if(binary_format) { gridfile.read((char *) &VAL,sizeof(VAL)) ; if(swapbytes) swapb((char *) &VAL,sizeof(VAL)); } else { gridfile >> VAL ; }

int main(int ac, char *av[]) {
  using namespace Loci ;
  string Lref = "NOSCALE" ;
  bool optimize = true ;
  Loci::Init(&ac,&av) ;

  if(Loci::MPI_processes != 1) {
    cerr << "msh2vog converter can only run serial at present" << endl ;
  }

  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 2 && !strcmp(av[1],"-v")) {
      cout << "Loci version: " << Loci::version() << endl ;
      if(ac == 2) {
	Loci::Finalize() ;
	exit(0) ;
      }
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-o")) {
      optimize = false ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-in")) {
      Lref = "1 inch" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-ft")) {
      Lref = "1 foot" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-cm")) {
      Lref = "1 centimeter" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-m")) {
      Lref = "1 meter" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-mm")) {
      Lref = "1 millimeter" ;
      ac-- ;
      av++ ;
   } else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      ac-- ;
      av++ ;
    }
  }

  if(Lref == "NOSCALE") {
    cerr << "Must set grid units!" << endl
         << "Use options -in, -ft, -cm, -m, or -Lref to set grid units." << endl ;
    exit(-1) ;
  }


  if(ac != 2) {
    cerr << "Usage: cobalt2vog <options> <file>" << endl
         << "Where options are listed below and <file> is the filename sans postfix" << endl
         << "flags:" << endl
         << "  -o  : disable optimization that reorders nodes and faces" << endl
         << "  -v  : display version" << endl
         << "  -in : input grid is in inches" << endl
         << "  -ft : input grid is in feet" << endl
         << "  -cm : input grid is in centimeters" << endl
         << "  -m : input grid is in meters" << endl
         << "  -Lref <units> : 1 unit in input grid is <units> long" << endl
         << endl ;
      exit(-1) ;
  }

  if(Lref == "")
    Lref = "1 meter" ;

  if(!isdigit(Lref[0])) {
    Lref = string("1") + Lref ;
  }
  Loci::UNIT_type tp ;
  istringstream iss(Lref) ;
  iss >> tp ;

  double posScale = tp.get_value_in("meter") ;

  if(Loci::MPI_rank == 0) {
    cout << "input grid file units = " << tp ;
    if(posScale != 1.0)
      cout << " = " << posScale << " meters " ;
    cout << endl ;
  }

  string input_name = av[1] ;
  string filename = VOG::ensureSuffix(input_name,".msh") ;
  string outfile = VOG::stripSuffix(input_name,".msh") + ".vog" ;

  ifstream gridfile(filename.c_str(),ios::in) ;
  if(gridfile.fail()) {
    cerr << "open failed on '" << filename << "'" << endl ;
    exit(-1) ;
  }
  bool binary_format = false ;
  bool swapbytes = false ;
  std::map<int,string> surfIds ;
  std::map<int,string> voltags ;
  std::map<int,int> physicalSurfaceTag ;
  std::map<int,int> physicalVolumeTag ;
  string text ;
  vector<vector3d<double> > vpos ;
  int ncells = 0 ;
  int cellid = 1 ;
  vector<size_t> posid ;
  vector<triaFace> tria ;
  vector<quadFace> quad ;

  size_t numNodes=0, minNodeTag=0,maxNodeTag=0 ;

  while(true) {
    getline(gridfile,text) ;
    if(gridfile.eof()) {
      break ;
    }
    if(checkToken(text,"$MeshFormat")) {
      double gmshver = 0 ;
      int bin=0, dsz=0 ;
      gridfile>>gmshver >> bin >> dsz ;
      if(gmshver < 4.0 || gmshver >=5.0) {
	cerr << "gmesh format converter only works for version 4 file format, file version=" << gmshver  << endl ;
	Loci::Abort() ;
      }
      cout << "GMesh Format: " << gmshver << endl ;
      if(bin != 0) {
	binary_format = true ;
	getnewline(gridfile) ;
	int val = 0 ;
	gridfile.read((char *)&val,sizeof(val)) ;
	if(val != 1)
	  swapbytes = true ;

	cout << "binary file " ;
	if(swapbytes) {
	  cout << " using swapbytes to correct endian form" ;
	}
	cout << endl ;

      }

      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndMeshFormat")) ;
    }
    if(checkToken(text,"$PhysicalNames")) {
      cout << "---- Begin Physical Names Section ----" << endl ;
      int numNames = 0 ;
      gridfile >> numNames ;
      for(int i=0;i<numNames;++i) {
	int dim=0,tag=0 ;
	string tagname ;
	gridfile >> dim >> tag >> tagname ;
	string dropquotes ;
	for(size_t i=0;i<tagname.size();++i) {
	  if(tagname[i] != '"')
	    dropquotes += tagname[i] ;
	}
	if(dim == 2) {
	  surfIds[tag] = dropquotes ;
	}
	if(dim == 3) {
	  voltags[tag] = dropquotes ;
	}
      }
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndPhysicalNames")) ;
      cout << "Named Surface Boundaries:" << endl ;
      for(std::map<int,string>::const_iterator mi=surfIds.begin();mi!=surfIds.end();++mi) {
	cout << mi->second << ":" <<mi->first << endl ;
      }
      if(voltags.begin() != voltags.end()) {
	cout << "Named Volumes:" << endl ;
	for(std::map<int,string>::const_iterator mi=voltags.begin();mi!=voltags.end();++mi) {
	  cout << mi->second << ":" <<mi->first << endl ;
	}
      }
      cout << "---- End Physical Names ----" << endl ;
    }
    if(checkToken(text,"$Entities")) {
      size_t npnts=0,ncurves=0,nsurf=0,nvol=0 ;
      READVAR(npnts) ;
      READVAR(ncurves) ;
      READVAR(nsurf) ;
      READVAR(nvol) ;

      for(size_t i=0;i<npnts;++i) {
	int ptag ;
	double x,y,z ;
	size_t nphystags=0 ;
	READVAR(ptag) ;
	READVAR(x) ;
	READVAR(y) ;
	READVAR(z) ;
	READVAR(nphystags) ;

	for(size_t j=0;j<nphystags;++j) {
	  int val ;
	  READVAR(val) ;
	}
      }
      for(size_t i=0;i<ncurves;++i) {
	int ctag ;
	double minX,minY,minZ,maxX,maxY,maxZ ;
	size_t numPhysicalTags ;
	READVAR(ctag) ;
	READVAR(minX) ;
	READVAR(minY) ;
	READVAR(minZ) ;
	READVAR(maxX) ;
	READVAR(maxY) ;
	READVAR(maxZ) ;
	READVAR(numPhysicalTags)

	for(size_t j=0;j<numPhysicalTags;++j) {
	  int tags ;
	  READVAR(tags) ;
	}
	size_t nbpts ;
	READVAR(nbpts) ;
	for(size_t i=0;i<nbpts;++i){
	  int val ;
	  READVAR(val) ;
	}
      }
      for(size_t i=0;i<nsurf;++i) {
	int surfTag ;
	double minX,minY,minZ,maxX,maxY,maxZ ;
	size_t numPhysicalTags ;

	READVAR(surfTag) ;
	READVAR(minX) ;
	READVAR(minY) ;
	READVAR(minZ) ;
	READVAR(maxX) ;
	READVAR(maxY) ;
	READVAR(maxZ) ;
	READVAR(numPhysicalTags)
	physicalSurfaceTag[surfTag] = 0 ;
	for(size_t j=0;j<numPhysicalTags;++j) {
	  int val ;
	  READVAR(val) ;
	  physicalSurfaceTag[surfTag] = val ;
	}
	size_t numBoundingSurfaces ;
	READVAR(numBoundingSurfaces) ;
	for(size_t j=0;j<numBoundingSurfaces;++j) {
	  int val ;
	  READVAR(val) ;
	}
      }
      for(size_t i=0;i<nvol;++i) {
	int volTag ;
	double minX,minY,minZ,maxX,maxY,maxZ ;
	size_t numPhysicalTags ;
	READVAR(volTag) ;
	READVAR(minX) ;
	READVAR(minY) ;
	READVAR(minZ) ;
	READVAR(maxX) ;
	READVAR(maxY) ;
	READVAR(maxZ) ;
	READVAR(numPhysicalTags)

	physicalVolumeTag[volTag] = 0 ;
	for(size_t j=0;j<numPhysicalTags;++j) {
	  int val ;

	  READVAR(val) ;

	  physicalVolumeTag[volTag] = val ;
	}
	size_t numBoundingSurfaces ;
	READVAR(numBoundingSurfaces) ;

	for(size_t j=0;j<numBoundingSurfaces;++j) {
	  int val ;
	  READVAR(val) ;
	}

      }
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndEntities")) ;

    }
    if(checkToken(text,"$Nodes")) {

      size_t numEntityBlocks=0 ;

      READVAR(numEntityBlocks) ;
      READVAR(numNodes) ;
      READVAR(minNodeTag) ;
      READVAR(maxNodeTag) ;

      cout << "Nodes: " << numNodes << " [" << minNodeTag << "," << maxNodeTag << "]" << endl ;
      if(maxNodeTag-minNodeTag+1 != numNodes ||
	 minNodeTag != 1) {
	cerr << "Converter expecting continuous node labels starting from 1"
	     << endl ;
	Loci::Abort() ;
      }
      vector<vector3d<double> > tmp(numNodes) ;
      vpos.swap(tmp) ;
      vector<size_t> itmp(numNodes) ;
      posid.swap(itmp) ;
      size_t cnt = 0 ;
      for(size_t i=0;i<numEntityBlocks;++i) {
	int dim=0,etag=0, parametric=0;
	size_t blksze=0 ;
	READVAR(dim) ;
	READVAR(etag) ;
	READVAR(parametric) ;
	READVAR(blksze) ;
	for(size_t j=0;j<blksze;++j){
	  size_t nid = 0 ;
	  READVAR(nid) ;
	  posid[cnt+j] = nid ;
	}
	for(size_t j=0;j<blksze;++j){
	  double x=0,y=0,z=0 ;
	  READVAR(x) ;
	  READVAR(y) ;
	  READVAR(z) ;
	  vpos[cnt+j] = vector3d<double>(x,y,z) ;
	}
	cnt += blksze ;
      }
      if(cnt != numNodes) {
	cerr << "did not read in numNodes cnt=" << cnt << " numNodes=" <<numNodes << endl ;
      }

      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndNodes")) ;
    }
    if(checkToken(text,"$Elements")) {
      size_t numEntityBlocks=0,numElements=0,minElemTag=0,maxElemTag=0 ;

      READVAR(numEntityBlocks) ;
      READVAR(numElements) ;
      READVAR(minElemTag) ;
      READVAR(maxElemTag) ;
      cout << "Elements: " << numElements << " [" << minElemTag << "," << maxElemTag << "]" << endl ;
      if(maxElemTag-minElemTag+1 != numElements ||
	 minElemTag != 1) {
	cerr << "Converter expecting continuous cell labels starting from 1"
	     << endl ;
	Loci::Abort() ;
      }
      for(size_t i=0;i<numEntityBlocks;++i) {
	int dim=0,tag=0,type=0 ;
	size_t nelem = 0 ;

	READVAR(dim) ;
	READVAR(tag) ;
	READVAR(type) ;
	READVAR(nelem) ;

	for(size_t j=0;j<nelem;++j) {
	  switch(type) {
	  case 1: // 2-node line
	    {
	      size_t id,pt1,pt2 ;
	      READVAR(id) ;
	      READVAR(pt1) ;
	      READVAR(pt2) ;
	    }
	    break ;
	  case 2: // 3-node triangle
	    {
	      size_t id, pt1,pt2,pt3 ;
	      READVAR(id) ;
	      READVAR(pt1) ;
	      READVAR(pt2) ;
	      READVAR(pt3) ;
	      int psid = physicalSurfaceTag[tag] ;
	      if(psid > 0) { // boundary surface
		triaFace tf ;
		tf.nodes[0] = pt1 ;
		tf.nodes[1] = pt2 ;
		tf.nodes[2] = pt3 ;
		tf.left = true ;
		tf.cell = -psid ;
		tria.push_back(tf) ;
	      } else {
		cerr << "id=" << id << ",psid=" << psid << endl ;
	      }

	    }
	    break ;
	  case 3: // 4-node quadrangle
	    {
	      size_t id, pt1,pt2,pt3,pt4 ;
	      READVAR(id) ;
	      READVAR(pt1) ;
	      READVAR(pt2) ;
	      READVAR(pt3) ;
	      READVAR(pt4) ;
	      int psid = physicalSurfaceTag[tag] ;
	      if(psid > 0) { // boundary surface
		quadFace qf ;
		qf.nodes[0] = pt1 ;
		qf.nodes[1] = pt2 ;
		qf.nodes[2] = pt3 ;
		qf.nodes[3] = pt4 ;
		qf.left = true ;
		qf.cell = -psid ;
		quad.push_back(qf) ;
	      }
	    }
	    break ;
	  case 4: // 4 node tetrahedron
	    {
	      size_t id, pt1,pt2,pt3,pt4 ;
	      READVAR(id) ;
	      READVAR(pt1) ;
	      READVAR(pt2) ;
	      READVAR(pt3) ;
	      READVAR(pt4) ;

	      triaFace tf ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt2 ;
	      tf.left = true ;
	      tf.cell = cellid ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt4 ;
	      tf.nodes[2] = pt3 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt2 ;
	      tf.nodes[2] = pt4 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt2 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt4 ;
	      tria.push_back(tf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 5: // 8 node hexahedron
	    {
	      size_t id, pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8 ;
	      READVAR(id) ;
	      READVAR(pt1) ;
	      READVAR(pt2) ;
	      READVAR(pt3) ;
	      READVAR(pt4) ;
	      READVAR(pt5) ;
	      READVAR(pt6) ;
	      READVAR(pt7) ;
	      READVAR(pt8) ;

	      quadFace qf ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt4 ;
	      qf.nodes[2] = pt3 ;
	      qf.nodes[3] = pt2 ;
	      qf.left = true ;
	      qf.cell = cellid ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt2 ;
	      qf.nodes[2] = pt6 ;
	      qf.nodes[3] = pt5 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt5 ;
	      qf.nodes[2] = pt8 ;
	      qf.nodes[3] = pt4 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt7 ;
	      qf.nodes[1] = pt3 ;
	      qf.nodes[2] = pt4 ;
	      qf.nodes[3] = pt8 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt7 ;
	      qf.nodes[1] = pt8 ;
	      qf.nodes[2] = pt5 ;
	      qf.nodes[3] = pt6 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt7 ;
	      qf.nodes[1] = pt6 ;
	      qf.nodes[2] = pt2 ;
	      qf.nodes[3] = pt3 ;
	      quad.push_back(qf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 6: // 6 node prism
	    {
	      size_t id, pt1,pt2,pt3,pt4,pt5,pt6 ;

	      READVAR(id) ;
	      READVAR(pt1) ;
	      READVAR(pt2) ;
	      READVAR(pt3) ;
	      READVAR(pt4) ;
	      READVAR(pt5) ;
	      READVAR(pt6) ;

	      triaFace tf ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt2 ;
	      tf.left = true ;
	      tf.cell = cellid ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt4 ;
	      tf.nodes[1] = pt5 ;
	      tf.nodes[2] = pt6 ;
	      tria.push_back(tf) ;
	      quadFace qf ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt2 ;
	      qf.nodes[2] = pt5 ;
	      qf.nodes[3] = pt4 ;
	      qf.left = true ;
	      qf.cell = cellid ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt4 ;
	      qf.nodes[2] = pt6 ;
	      qf.nodes[3] = pt3 ;
	      quad.push_back(qf) ;
	      qf.nodes[0] = pt2 ;
	      qf.nodes[1] = pt3 ;
	      qf.nodes[2] = pt6 ;
	      qf.nodes[3] = pt5 ;
	      quad.push_back(qf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 7: // 5 node pyramid
	    {
	      size_t id, pt1,pt2,pt3,pt4,pt5 ;
	      READVAR(id) ;
	      READVAR(pt1) ;
	      READVAR(pt2) ;
	      READVAR(pt3) ;
	      READVAR(pt4) ;
	      READVAR(pt5) ;

	      triaFace tf ;
	      tf.nodes[0] = pt1 ;
	      tf.nodes[1] = pt2 ;
	      tf.nodes[2] = pt5 ;
	      tf.left = true ;
	      tf.cell = cellid ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt2 ;
	      tf.nodes[1] = pt3 ;
	      tf.nodes[2] = pt5 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt3 ;
	      tf.nodes[1] = pt4 ;
	      tf.nodes[2] = pt5 ;
	      tria.push_back(tf) ;
	      tf.nodes[0] = pt4 ;
	      tf.nodes[1] = pt1 ;
	      tf.nodes[2] = pt5 ;
	      tria.push_back(tf) ;
	      quadFace qf ;
	      qf.nodes[0] = pt1 ;
	      qf.nodes[1] = pt4 ;
	      qf.nodes[2] = pt3 ;
	      qf.nodes[3] = pt2 ;
	      qf.left = true ;
	      qf.cell = cellid ;
	      quad.push_back(qf) ;
	      cellid++ ;
	      ncells++ ;
	    }
	    break ;
	  case 15: // 1 point node
	    {
	      size_t id, pt1 ;
	      READVAR(id) ;
	      READVAR(pt1) ;
	      cerr << "1 point node not supported" << endl ;
	    }
	    break ;
	  default:
	    cout << "element type= " << type << " unsupported by converter."
		 << endl ;
	    Loci::Abort() ;
	  }
	}
      }
      do {
	getline(gridfile,text) ;
      } while(!checkToken(text,"$EndElements")) ;
      cout << "Read " << ncells << " Cells" << endl ;
    }
  }

  // prepare triangle faces (sort them)
  for(size_t i=0;i<tria.size();++i) {
    // pos numbers nodes from zero
    tria[i].nodes[0] -= 1 ;
    tria[i].nodes[1] -= 1 ;
    tria[i].nodes[2] -= 1 ;

    if(tria[i].nodes[0] > tria[i].nodes[1]) {
      std::swap(tria[i].nodes[0],tria[i].nodes[1]) ;
      tria[i].left = !tria[i].left ;
    }
    if(tria[i].nodes[0] > tria[i].nodes[2]) {
      std::swap(tria[i].nodes[0],tria[i].nodes[2]) ;
      tria[i].left = !tria[i].left ;
    }
    if(tria[i].nodes[1] > tria[i].nodes[2]) {
      std::swap(tria[i].nodes[1],tria[i].nodes[2]) ;
      tria[i].left = !tria[i].left ;
    }
  }
  //sort tria
  sort(tria.begin(),tria.end(),triaCompare) ;

  // prepare quad faces (sort them, but be careful)
  for(size_t i=0;i<quad.size();++i) {
    // pos numbers nodes from zero
    quad[i].nodes[0] -=1 ;
    quad[i].nodes[1] -=1 ;
    quad[i].nodes[2] -=1 ;
    quad[i].nodes[3] -=1 ;
    // First make sure first entry is lowest number
    int tmp_face[4] ;
    int vs = quad[i].nodes[0] ;
    size_t nv = 0 ;
    for(size_t j=1;j<4;++j)
      if(vs > quad[i].nodes[j]) {
        vs = quad[i].nodes[j] ;
        nv = j ;
      }
    for(size_t j=0;j<4;++j)
      tmp_face[j] = quad[i].nodes[(j+nv)&0x3] ;
    // next make orientation so that it will match other face
    if(tmp_face[1] < tmp_face[3])
      for(int j=0;j<4;++j)
        quad[i].nodes[j] = tmp_face[j] ;
    else {
      for(size_t j=0;j<4;++j)
        quad[i].nodes[j] = tmp_face[(4 - j) &0x3 ] ;
      quad[i].left = !quad[i].left ;
    }
  }
  sort(quad.begin(),quad.end(),quadCompare) ;


  // Ok now sort out how many triangle and quad faces there are

  vector<int> tria_equals ;
  tria_equals.reserve(tria.size()/2) ;
  for(size_t i=0;i<tria.size();++i) {
    int cnt = 1 ;
    size_t j=i+1 ;
    while(j<tria.size() && triaEqual(tria[i],tria[j])) {
      cnt++ ;j++ ;
    }
    i=j-1;
    tria_equals.push_back(cnt) ;
  }

  cout << "ntrias=" << tria_equals.size() << endl ;

  vector<int> quad_equals ;
  quad_equals.reserve(quad.size()/2) ;
  for(size_t i=0;i<quad.size();++i) {
    int cnt = 1 ;
    size_t j=i+1 ;
    while(j<quad.size() && quadEqual(quad[i],quad[j])) {
      cnt++ ;j++ ;
    }
    i=j-1;
    quad_equals.push_back(cnt) ;
  }

  cout << "nquads=" << quad_equals.size() << endl ;


  int ntria = tria_equals.size() ;
  int nquad = quad_equals.size() ;
  int nfaces = ntria+nquad ;

  int facebase = std::numeric_limits<int>::min()+2048 ;

  vector<int> facesizes(MPI_processes) ;
  MPI_Allgather(&nfaces,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    facebase += facesizes[i] ;

  if(Loci::MPI_rank == 0) {
    for(int i=0;i<MPI_processes;++i)
      if(facesizes[i] == 0) {
	cerr << "Run msh2vog with fewer than " << i << " processors for this mesh!" << endl ;
	Loci::Abort() ;
	break ;
      }
  }
  entitySet faces = interval(facebase,facebase+nfaces-1) ;

  if(nfaces == 0) {
    faces = EMPTY ;
  }
  store<int> count ;
  Map cl,cr ;

  cl.allocate(faces) ;
  cr.allocate(faces) ;

  count.allocate(faces) ;
  int fc = facebase ;

  for(int i=0;i<ntria;i++)
    count[fc++] = 3 ;

  for(int i=0;i<nquad;i++)
    count[fc++] = 4 ;


  multiMap face2node ;
  face2node.allocate(count) ;
  fc = facebase ;
  int ccerror = 0 ;
  int unmatched_boundary_faces = 0 ;

  int i=0 ;
  for(size_t eq=0;eq<tria_equals.size();++eq) {
    int corner[3];
    if(tria[i].left) {
      for(int j=0;j<3;++j)
	corner[j] = tria[i].nodes[j] ;
    } else {
      for(int j=0;j<3;++j)
	corner[j] = tria[i].nodes[2-j] ;
    }
    int f1 = i ;
    int f2 = i+1 ;
    int c1 = tria[f1].cell ;
    int c2 = tria[f2].cell ;

    if(c1 < 0 ) {
      if(tria_equals[eq] == 1) {
	cerr << "isolated boundary face!" << endl
	     << "mesh doesn't seem to be tagging boundaries properly for tag " << -c1 << endl ;
	Loci::Abort() ;
      }
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      cout << "c1=" << c1 << "c2=" << c2 << " out of " << tria_equals[eq] << endl ;
      Loci::Abort() ;
    }

    if(tria_equals[eq] == 1) {
      c2 = -1024 ;
      f2 = f1 ;
      unmatched_boundary_faces++ ;
    }
    if(tria_equals[eq] > 2) {
      cerr << "There appear to be interior tagged faces.  " << endl
	   << "Interior tagged faces are not supported in the vog file format. "
	   << endl ;
    }
    cl[fc] = c1 ;
    cr[fc] = c2 ;
    if(c1 >=0 && c2 >=0 &&(tria[f1].left == tria[f2].left)) {
      cerr << "consistency error " << c1 << ' ' << c2 << endl ;
      ccerror++ ;
    }

    for(int j=0;j<3;++j) {
      face2node[fc][j] = corner[j] ;
    }
    fc++ ;
    i += tria_equals[eq] ;
  }


  i=0 ;
  for(size_t eq=0;eq<quad_equals.size();++eq) {
    int corner[4];
    if(quad[i].left) {
      for(int j=0;j<4;++j)
	corner[j] = quad[i].nodes[j] ;
    } else {
      for(int j=0;j<4;++j)
	corner[j] = quad[i].nodes[3-j] ;
    }
    int f1 = i ;
    int f2 = i+1 ;
    int c1 = quad[f1].cell ;
    int c2 = quad[f2].cell ;

    if(c1 < 0 ) {
      if(quad_equals[eq] == 1) {
	cerr << "isolated boundary face!" << endl
	     << "mesh doesn't seem to be tagging boundaries properly for tag " << -c1 << endl ;
	Loci::Abort() ;
      }
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      cout << "c1=" << c1 << "c2=" << c2 << " out of " << quad_equals[eq] << endl ;
      Loci::Abort() ;
    }

    if(quad_equals[eq] == 1) {
      c2 = -1024 ;
      f2 = f1 ;
      unmatched_boundary_faces++ ;
    }
    if(quad_equals[eq] > 2) {
      cerr << "There appear to be interior tagged faces.  " << endl
	   << "Interior tagged faces are not supported in the vog file format. "
	   << endl ;
    }
    cl[fc] = c1 ;
    cr[fc] = c2 ;
    if(c1>=0 && c2 >=0 &&(quad[f1].left == quad[f2].left)) {
      cerr << "consistency error" << endl ;
      ccerror++ ;
    }

    for(int j=0;j<4;++j) {
      face2node[fc][j] = corner[j] ;
    }
    fc++ ;
    i += quad_equals[eq] ;
  }

  if(unmatched_boundary_faces > 0) {
    cerr << "There were " << unmatched_boundary_faces << " boundary faces that were not matched!" << endl
	 << "  ** Check the mesh to make sure all boundary faces are tagged." << endl
	 << "  ** The vog file has given these faces the 1024 tag." << endl ;
  }
  if(ccerror > 0) {
    cerr << "consistency error #=" << ccerror << endl ;
    cerr << "something wrong with mesh, have you given physical tags to the surfaces?" << endl ;
    Loci::Abort() ;
  }
  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  store<vector3d<double> > pos ;
  entitySet nodes = interval(0,numNodes-1) ;
  pos.allocate(nodes) ;
  for(size_t i=0;i<numNodes;++i)
    pos[i] = vpos[i]*posScale ;

  vector<pair<int,string> > surf_ids ;
  for(std::map<int,string>::const_iterator mi=surfIds.begin();mi!=surfIds.end();++mi) {
    surf_ids.push_back(pair<int,string>(mi->first,mi->second)) ;
  }

  VOG::colorMatrix(pos,cl,cr,face2node) ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  //  if(MPI_rank == 0)
    //    cerr << "done  coloring" << endl ;
  if(optimize) {
    if(MPI_rank == 0)
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }

  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids, MPI_COMM_WORLD) ;
  Loci::Finalize() ;
  return 0 ;
}
