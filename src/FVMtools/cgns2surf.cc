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

/*cgns library can be compiled with 64-bit or 32-bit,
  which will decide cgsize_t is long int or int
  Loci library can also be compiled with 64-bit or 32-bit, which will affect gEntity
  static_cast is used in case the size of gEntity and cgsize_t are different
*/

/*
  This program take a cgns path as input, which is written as "<filename>/<basename>/<zonename>".

  CGNS file requirements:
    1. Only unstructured grid. Cell dimension is 3D, and physical dimension is 3D.
    2. Allowed element type: TRI_3, QUAD_4,
    TETRA_4,  PYRA_5,
    PENTA_6,  HEXA_8,
    3. assume boundary names are specified in ZoneBC nodes, and boundary conditions are defined on faces
*/



#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <utility>

#include <ctype.h>
#include <iostream>
#ifdef _WIN32
#include <io.h>
#define unlink _unlink
#else
#include <unistd.h>
#endif

#ifdef USE_CGNS



#include "cgnslib.h"
#include "hash.h"

#ifndef CGNS_ENUMV
# define CGNS_ENUMV(V) V
# define CGNS_ENUMT(T) T
#endif

using std::cout ;
using std::endl ;
using std::cerr ;
using std::ios ;
using std::string ;
using std::vector ;
using std::istringstream ;
using std::pair;
using std::ofstream;
using std::min;
using std::max;


typedef struct {
  cgsize_t id;
  int bcnum;
  int nnodes;
  cgsize_t nodes[4];
} FACE;

//the last nodes is bc num
typedef struct {
  int nodes[4];
} Tria;

typedef struct {
  int nodes[5];
} Quad;


static void error_exit(const char *func)
{
  printf("CGNSlib ERROR:");
  if (func != NULL && *func)
    printf("%s:", func);
  printf("%s\n", cg_get_error());
  exit(1);
}


/*--------------------------------------------------------------------*/

static int sort_points (const void *v1, const void *v2)
{
  return (int)(*((cgsize_t *)v1) - *((cgsize_t *)v2));
}

/*--------------------------------------------------------------------*/

static int find_point (cgsize_t id, cgsize_t np, cgsize_t *pts)
{
  
  cgsize_t lo = 0, hi = np - 1, mid;

  if (!np || id < pts[0]) return 0;
  if (id == pts[0]) return 1;
  if (!hi || id > pts[hi]) return 0;
  if (id == pts[hi]) return 1;

  while (lo <= hi) {
    mid = (lo + hi) >> 1;
    if (id == pts[mid]) return 1;
    if (id < pts[mid])
      hi = mid - 1;
    else
      lo = mid + 1;
  }
  return 0;
}
  




/*--------------------------------------------------------------------*/
/* set the bc_num of all faces on this boundary to nb */ 
static void unstructured_boundary (int nb,
                                   const vector<int>& fid,
                                   CGNS_ENUMT(PointSetType_t) ptype,
                                   CGNS_ENUMT(GridLocation_t) location,
                                   cgsize_t np,
                                   cgsize_t *ptset,
                                   vector<int>& fsurfid)
{
  cgsize_t  is, ie;
 
  if (ptype == CGNS_ENUMV(PointRange) ||
      ptype == CGNS_ENUMV(PointList)) {
    if (location == CGNS_ENUMV(FaceCenter) ||
        location == CGNS_ENUMV(CellCenter)) {
      ptype = (ptype == CGNS_ENUMV(PointRange) ?
               CGNS_ENUMV(ElementRange) : CGNS_ENUMV(ElementList));
    }
    else if (location == CGNS_ENUMV(Vertex)) {
      error_exit("The boundary condition is defined on Vertex"); 
    }
  }
  
  if (ptype == CGNS_ENUMV(ElementRange)) {
    if (ptset[0] < ptset[1]) {
      is = ptset[0];
      ie = ptset[1];
    }
    else {
      is = ptset[1];
      ie = ptset[0];
    }
    
    for (size_t nf = 0; nf < fid.size(); nf++) {
      int id = fid[nf];
      if(id>=is && id <=ie)fsurfid[nf] = nb;
    }
    return;
  }

  if (ptype == CGNS_ENUMV(ElementList)) {
    qsort(ptset, (size_t)np, sizeof(cgsize_t), sort_points);
    for (size_t nf = 0; nf < fid.size(); nf++) {
      int id = fid[nf];
      if (find_point(id, np, ptset)){
        fsurfid[nf] = nb;
      }
    }
    return;
  }
  cout << "ptype " << ptype << endl;
  cout << "location " << location << endl;
  cout<< "cgns2surf assume the boundary names are specified in ZoneBC node, and boundary conditions are specified on faces" << endl; 
  error_exit("The boundary condition is not defined on faces ");
}
/* find if bc_name is in boundaries  */
static bool find_bc_name(const string& bc_name, const vector<string>& boundaries){
  bool found = false;
  for(size_t i = 0 ; i < boundaries.size(); i++){
    if(bc_name == boundaries[i]){
      found = true;
      break;
    }
  }
  return found;
}

/* read in the specified boundary conditions from cgns file, and then set the bc_bum of related faces  */
static void boundary_conditions  (cgsize_t cgFile, cgsize_t cgBase, cgsize_t cgZone,
                                  const vector<string>& boundaries,//only read in the specified boundaries 
                                  const vector<int>& fid,
                                  vector<int>& fsurfid)
{
  int nb, ib, nrmlindex[3];
  cgsize_t is;
  cgsize_t np; //num_points
 
  char name[33];
  CGNS_ENUMT(BCType_t) bctype;
  CGNS_ENUMT(PointSetType_t) ptype;
  CGNS_ENUMT(GridLocation_t) location;
  CGNS_ENUMT(DataType_t) datatype;
  int dim =  1; //for unstructured mesh
  int nBocos = 0;
  if (cg_nbocos (cgFile, cgBase, cgZone, &nBocos))
    error_exit ("cg_nbocos");
  cout<< "nBcocos " << nBocos << endl;
  if (nBocos) {
    for (nb = 1; nb <= nBocos; nb++) {
      if (cg_boco_info(cgFile, cgBase, cgZone, nb, name,
                       &bctype, &ptype, &np, nrmlindex, &is,
                       &datatype, &ib))
        error_exit("cg_boco_info");
      string bc_name = string(name);
      cout<< " bc_name " << bc_name << endl;
      if(find_bc_name(bc_name, boundaries)){
        cgsize_t *ptset;//point  set
        if (cg_boco_gridlocation_read(cgFile, cgBase, cgZone,
                                      nb, &location))
          error_exit("cg_boco_gridlocation_read");
        ptset = (cgsize_t *)malloc((size_t)(np * dim) * sizeof(cgsize_t));
        if (ptset == NULL)
          error_exit( "malloc failed for boco ptset");
        if (cg_boco_read(cgFile, cgBase, cgZone, nb, ptset, 0))
          error_exit("cg_boco_read");
        cout << " name found " << endl;
        //set up the bcnum of faces
        unstructured_boundary(nb, fid, ptype, location, np, ptset, fsurfid);
      
        free(ptset);
      }
    }
   
  }
}


/*read in all the faces, and then set up bc_num of faces in boundaries
  this function assume the cgns file :
  1. is unstructured grid
  2. has a volume zone that contains all the boundary face elements and ZoneBC node
*/
static void unstructured_faces ( cgsize_t cgFile, cgsize_t cgBase, cgsize_t cgZone,
                                 const vector<string>& boundaries,
                                 vector<int>& node_set,
                                 vector<Tria>& trias,
                                 vector<Quad>& quads
                                 ){
  //read in all the faces in volume zone
  vector<int> f2size ;
  vector<int> fsurfid ;
  vector<int> f2node ;
  vector<int> fid;
 
  int nsect = 0, nn = 0, ip = 0;
  
  CGNS_ENUMT(ElementType_t) elemtype, et;
  char name[33];
  cgsize_t is, ie, ne;
  cgsize_t size, *conn;
 
  
  if (cg_nsections (cgFile, cgBase, cgZone, &nsect))
    error_exit ("cg_nsections");
  for (int ns = 1; ns <= nsect; ns++) {
    int nf = 0;
    if (cg_section_read (cgFile, cgBase, cgZone, ns,
                         name, &elemtype, &is, &ie, &nn, &ip))
      error_exit ("cg_section_read");
    if (elemtype < CGNS_ENUMV(TRI_3) ||
        (elemtype > CGNS_ENUMV(QUAD_9) &&
         elemtype != CGNS_ENUMV(MIXED))) continue; //only read in tria and quads
    ne = ie - is + 1;
    if (cg_ElementDataSize (cgFile, cgBase, cgZone, ns, &size))
      error_exit ("cg_ElementDataSize");
    conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
    if (conn == NULL)
      error_exit ("malloc failed for element connectivity");
    if (cg_elements_read (cgFile, cgBase, cgZone, ns, conn, NULL))
      error_exit ("cg_elements_read");

    et = elemtype;
    for (cgsize_t i = 0, n = 0; n < ne; n++) {
      if (elemtype == CGNS_ENUMV(MIXED))
        et = (CGNS_ENUMT(ElementType_t))conn[i++];
      switch (et) {
      case CGNS_ENUMV(TRI_3):
        nf = 3;
        break;
      case CGNS_ENUMV(QUAD_4):
        nf = 4;
        break;
      default:
        nf = 0;
        break;
      }
      if (nf) {
        f2size.push_back(nf);
        fid.push_back(n+is);
        fsurfid.push_back(-ns);
        for (nn = 0; nn < nf; nn++)
          f2node.push_back(conn[i+nn]);
      }
      cg_npe (et, &nn);
      i += nn;
    }
    free (conn);
  }
  //read boudary conditions in cgns file, and set up the surf_id for selected faces
  boundary_conditions ( cgFile, cgBase, cgZone, boundaries,fid, fsurfid);
  //remap nodes
  size_t off = 0;
  for(size_t fi = 0; fi < fsurfid.size(); fi++){
    if(fsurfid[fi]>0){
      if(f2size[fi]==3){
        Tria tria;
        for(int i = 0; i < 3; i++){
          tria.nodes[i] = f2node[i+off];
          node_set.push_back(f2node[i+off]);
        }
        tria.nodes[3] = fsurfid[fi];
        trias.push_back(tria);
      }else if(f2size[fi]==4){
        Quad quad;
        for(int i = 0; i < 4; i++){
          quad.nodes[i] = f2node[i+off];
          node_set.push_back(f2node[i+off]);
        }
        quad.nodes[4] = fsurfid[fi];
        quads.push_back(quad);
      }else{
        cerr<<" General faces not allowed " << endl;
        exit(1);
      }
    }
    off += f2size[fi];
  }

  if(node_set.size() == 0) {
    cerr << "The input boundary name should match that  in cgns file ZoneBC node !" << endl;
    cerr << "nothing to process in surf extract! Exiting."
         << endl ;
    exit(-1) ;
  }
  
  sort(node_set.begin(),node_set.end()) ;
  node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;

  int mx = node_set[0] ;
  int mn = node_set[0] ;
  for(size_t i=0;i<node_set.size();++i) {
    mx = max(mx,node_set[i]) ;
    mn = min(mn,node_set[i]) ;
  }
  
  vector<int> nmap(mx+1, -1);
  for(size_t i=0;i<node_set.size();++i) {
    nmap[node_set[i]] = i+1 ; //assume original numbering start with 1?
  }

  for(size_t i = 0; i< trias.size(); i++){
    for(int j = 0; j< 3; j++){
      trias[i].nodes[j] = nmap[trias[i].nodes[j]];
    }
  }
  
  
  for(size_t i = 0; i< quads.size(); i++){
    for(int j = 0; j< 4; j++){
      quads[i].nodes[j] = nmap[quads[i].nodes[j]];
    }
  }
}


/*read from cgns file,
  The coordinates are read in as double precision, the format in file can be in single or double precision 
*/
void readCGNS( string pathname,
               const vector<string>& boundaries,
               vector<double > &pos,
               vector<int>& node_set,
               vector<Tria>& trias,
               vector<Quad>& quads
               ) {
  
  std::size_t p1 = pathname.find('/');
  std::size_t p2 = pathname.find('/', p1+1);

  string filename = pathname.substr(0, p1);
  string basename = pathname.substr(p1+1, p2-p1-1);
  string zonename = pathname.substr(p2+1);
  cout << " pathname " << pathname << endl; 
  cout << " filename " << filename << endl;
  cout << " basename " << basename << endl;
  cout << " zonename " << zonename << endl;
 
 
  int celldim = 3, phydim = 3;
 
  //sizes: Number of vertices, cells, and boundary vertices in each (index)-dimension.
  cgsize_t sizes[3]; 
  cgsize_t start_index, end_index;
  
  int num_bases=0, num_zones = 0;
  int index_file = 0, index_base=0, index_zone=0;
  cgsize_t num_nodes=0;
  
  cgsize_t errs = 0;
  char bname[33], zname[33];
  CGNS_ENUMT(ZoneType_t) ztype;
    
  cout <<"start reading cgns file" << endl;
  
  if(cg_open (filename.c_str(), CG_MODE_READ, &index_file)) error_exit(" unable to open CGNS grid file ");
  if(cg_nbases (index_file, &num_bases))error_exit("error reading number of bases");
  for(index_base = 1; index_base<= num_bases; index_base++){
    if(cg_base_read (index_file, index_base, bname, &celldim, &phydim))error_exit("error reading base information");
    if(basename == string(bname)) break;
  }
  if(index_base > num_bases){
    cg_close(index_file);
    cerr<<"base "<< basename <<" not found"<< endl;
    exit(-1);
  }
  if(celldim != 3 || phydim != 3){
    cg_close(index_file);
    cerr << "only 3D cell and physical dimensions are allowed in CGNS file" << endl;
    exit(-1);
  }
  if(cg_nzones (index_file, index_base, &num_zones)) error_exit("error reading number of zones");
  for(index_zone = 1; index_zone <= num_zones; index_zone++){
    if(cg_zone_read(index_file, index_base, index_zone, zname, sizes))error_exit("error reading zone information");
    if(zonename == string(zname)) break;
  }
  if(index_zone > num_zones){
    cg_close(index_file);
    cerr<<"zone "<< zonename <<" not found"<< endl;
    exit(-1);
  }

  
  if(cg_zone_type (index_file, index_base, index_zone, &ztype))error_exit("error reading zone type");
  if (ztype != Unstructured) error_exit("can only handle unstructured grid in CGNS file");
  if(cg_zone_read (index_file, index_base, index_zone, zname, sizes))error_exit("error reading zone information");
    
  //3D unstructured sizes:		NVertex, NCell3D, NBoundVertex
  num_nodes = sizes[0];
  if(num_nodes < 1){
    cg_close(index_file);
    cerr << "number of nodes < 1 " << endl;
    exit(-1);
  }
  
  cgsize_t mxsize = num_nodes;
  pos.resize(mxsize*3) ;
  if(num_nodes > 0){ // Read in positions
   
    start_index = 1;
    end_index = num_nodes;
        
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateX", RealDouble,
                          &start_index, &end_index, &pos[0]);
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateY", RealDouble,
                          &start_index, &end_index, &pos[mxsize]);
    
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateZ", RealDouble,
                          &start_index, &end_index, &pos[2*mxsize]); 
    
  }
  
 

  //read in trias and quads, then hash and sort faces, get boundary information
  
  cout << "start  unstructured_faces " << endl;
  unstructured_faces (index_file, index_base, index_zone,
                      boundaries,
                      node_set,
                      trias,
                      quads
                      );

  cout << "end  unstructured_faces " << endl; 
  cout <<"finish reading cgns file" << endl;
        
}



  
void write_surf(const string& filename,
                const vector<double>& pos,
                const vector<int>& node_set,
                const vector<Tria>& trias,
                const vector<Quad>& quads
                ){
  
 
 
  ofstream ssfile(filename.c_str(),ios::out) ;
  ssfile.precision(15) ;

  int ntri = trias.size();
  int nqua = quads.size();
   
  ssfile << ntri << ' ' << nqua << ' ' 
         << node_set.size() << endl ;
  
  
  double normal_spacing = 0 ;
  int offset = pos.size()/3;
  for(size_t i=0;i<node_set.size();++i) {
    
    int nd = node_set[i]-1 ;
    ssfile << pos[nd] << ' ' << pos[nd+offset] << ' ' << pos[nd+2*offset]
           << ' '<< normal_spacing << endl ;
  }

  
  // output triangles
  for(int i = 0; i < ntri; i++){
    for(int j = 0 ; j < 4; j++){
      ssfile << trias[i].nodes[j] << ' ';
    }
    ssfile << "0 1" << endl;
  }
    
  
  // output quads
  for(int i = 0; i < nqua; i++){
    for(int j = 0 ; j < 5; j++){
      ssfile << quads[i].nodes[j] << ' ';
    }
    ssfile << "0 1" << endl;
  }
  ssfile.close();
}


int main(int ac, char* av[]) {
  
  std::string zonepath;
  vector<string> bcnames;
  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-bc")){
      bcnames.push_back(string(av[2]));
      ac -= 2 ;
      av += 2 ;
    }
  }
  if(ac == 2) {
    zonepath = av[1];
  }else{
    cerr << "Usage: cgns2surf <options> <zonepath>" << endl
         << "Where options are listed below and <zonepath> is the cgns path written as filename/basename/zonename" << endl
         << "flags:" << endl
         <<" -bc <string> : boundary name to extract, this option can repeat itself " << endl
         <<" cgns2surf assumes that boundary names are specified in ZoneBC node, and the boundary conditions are defined on faces"<< endl
         << endl ;
    exit(-1) ;
  }
  
  int loc = 0;
  loc = zonepath.find('.') ;
  std::string casename = zonepath.substr(0, loc) ;
  string outfile = string(casename) + string(".surf") ;
  if(bcnames.size() ==1)outfile = bcnames[0]+string(".surf") ;
  cout << "outfile " << outfile << endl; 
  vector<double> pos ;
  vector<int> node_set;
  vector<Tria> trias;
  vector<Quad> quads;
  readCGNS(zonepath, bcnames,
           pos,  node_set, trias, quads); 
  
  //write surf
  write_surf(outfile, pos, node_set, trias, quads);
}

#else
using namespace std ;

int main(int ac, char *av[]) {
  cerr << "Loci not compiled with CGNS support enabled! This utility cannot work!" << endl ;
  return -1 ;
}

#endif
