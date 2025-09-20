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
#include "vogtools.h"
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

void input_error() {
  cerr << "error reading file" << endl ;
  exit(-1) ;
}

int readCobaltMesh(string filename,
                   store<vector3d<double> > &pos,
                   Map &cl, Map &cr, multiMap &face2node) {
  FILE* IFP = fopen(filename.c_str(), "r") ;
  if(IFP == NULL) {
    return -1 ;
  }
  int npatch, maxppf, maxfpc ;
  int ndim, nzones ;
  int npnts, nfaces, ncells ;
  if(fscanf(IFP, "%d%d%d", &ndim, &nzones, &npatch)!=3)
    input_error() ;
  if(ndim != 3) {
    cerr << "currently only supports 3-D grid files" << endl ;
    return -1 ;
  }
  if(nzones != 1) {
    cerr << "currently only supports single zone files" << endl ;
    return -1 ;
  }
  if(fscanf(IFP, "%d%d%d%d%d", &npnts, &nfaces, &ncells, &maxppf, &maxfpc)!=5)
    input_error() ;
  cout << " ndim = " << ndim << endl ;
  cout << " npnts = " << npnts << endl ;
  cout << " nfaces = " << nfaces << endl ;
  cout << "ncells = " << ncells << endl ;

  // read in pos
  entitySet pdom = interval(0,npnts-1) ;
  pos.allocate(pdom) ;
  for(int i = 0; i < npnts; ++i) {
    if(fscanf(IFP, "%lf", &pos[i].x)!=1)
      input_error() ;
    if(fscanf(IFP, "%lf", &pos[i].y)!=1)
      input_error() ;
    if(fscanf(IFP, "%lf", &pos[i].z)!=1)
      input_error() ;
  }
  
  fpos_t after_pos ;
  fgetpos(IFP, &after_pos) ;

  int dummy_node ; 
  entitySet fdom = interval(0,nfaces-1) ;
  store<int> count ;
  count.allocate(fdom) ;
  cl.allocate(fdom) ;
  cr.allocate(fdom) ;
  
  for(int i = 0; i < nfaces; ++i) {
    if(fscanf(IFP, "%d", &count[i])!=1)
      input_error() ;
    for(int k = 0; k < count[i]; ++k)
      if(fscanf(IFP, "%d", &dummy_node)!=1)
	input_error() ; 
    if(fscanf(IFP, "%d%d", &cl[i], &cr[i])!=2)
      input_error() ; 
  }
  face2node.allocate(count) ;
  const fpos_t tmp_position = after_pos ;
  fsetpos(IFP, &tmp_position) ;
  for(int i = 0; i < nfaces; ++i) {
    int fsz = 0 ;
    if(fscanf(IFP, "%d", &fsz)!=1)
      input_error() ;
    for(int k = 0; k < fsz; ++k) 
      if(fscanf(IFP, "%d", &face2node[i][k])!=1)
	input_error() ;
    for(int k = 0; k < fsz; ++k)
      face2node[i][k] -= 1 ;
    if(fscanf(IFP, "%d%d", &cl[i], &cr[i])!=2)
      input_error() ; 
  }
  fclose(IFP) ;
  return(0) ;
}

int main(int ac, char *av[]) {
  using namespace Loci ;
  using namespace VOG ;

  bool optimize = true ;
  
  Loci::Init(&ac,&av) ;

  string Lref = "NOSCALE" ;
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

  string filename = av[1] ;
  filename += ".cog" ;

  string outfile = av[1] ;
  outfile += ".vog" ;

  store<vector3d<double> > pos;
  Map cl, cr;
  multiMap face2node;

  // Read in Mesh
  int status = readCobaltMesh(filename,pos,cl,cr,face2node) ;
  if(status != 0) {
    cerr << "unable to read file '" << filename << "'" << endl ;
  }

  if(posScale != 1.0) {
    FORALL(pos.domain(),nd) {
      pos[nd] *= posScale ;
    } ENDFORALL ;
  }

  // establish face left-right orientation
  if(MPI_rank == 0) 
    cerr << "orienting faces" << endl ;
  VOG::orientFaces(pos,cl,cr,face2node) ;
    
  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;

  if(optimize) {
    if(MPI_rank == 0) 
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }
  
  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;

  vector<BC_descriptor> bcs ;

  if(MPI_rank == 0) {
    string tagsfile = string(av[1]) + ".tags" ;
    bcs = readTags(tagsfile) ;
    if(bcs.size() != 0) {
      cout << "boundary faces:"<< endl ;
      for(size_t i=0;i<bcs.size();++i)
        cout << bcs[i].name << ' ' << bcs[i].id << endl ;
    }
  }

  vector<pair<int,string> > surf_ids ;
  for(size_t i=0;i<bcs.size();++i)
    surf_ids.push_back(pair<int,string>(bcs[i].id,bcs[i].name)) ;
  
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;

  Loci::Finalize() ;
  return 0 ;
}
  
