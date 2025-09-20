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

/*Read cgns:
  Only read in Grid specification part, i.e., grid coordinates, element connectivity, and boundary conditions    
  1. All data located in one file. The path is provided by user  as
  as "<filename>/<basename>/<zonename>".  If <zonename> is not given, read in the first zone of the base,
  If <basename is not given, read in the first base of the file.
  2. Only unstructured grid. Cell dimension is 3D, and physical dimension is 3D.
  3. The data read only contains core vertex coordinates and core elements, and boundary conditions
  4. Allowed element type: TRI_3, QUAD_4,
  TETRA_4,  PYRA_5, 
  PENTA_6,  HEXA_8
  5.assume that boundary names are specified in ZoneBC node, if there is spaces in boundary names, replace them with '_',
  and assume the boundary conditions are defined on faces
*/

/*write vog
  1. read in boundary surface and boundary names
  2. remove transparent faces
  3. post scale unit
  
*/


#include <strings.h>
#include <Loci.h>
#include <Tools/simple_partition_long.h>
#include "vogtools.h"
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

#include <stdio.h>
#include <ctype.h>
#include <string.h>
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
using Loci::Array ;
using Loci::MPI_rank ;
using Loci::MPI_processes ;

typedef struct {
  cgsize_t id;
  int bcnum;
  int nnodes;
  cgsize_t nodes[4];
} FACE;

// the faces of different CGNS elements
static int facenodes[20][5] = {
  /* tet */
  {3, 0, 2, 1, 0},
  {3, 0, 1, 3, 0},
  {3, 1, 2, 3, 0},
  {3, 2, 0, 3, 0},
  /* pyramid */
  {4, 0, 3, 2, 1},
  {3, 0, 1, 4, 0},
  {3, 1, 2, 4, 0},
  {3, 2, 3, 4, 0},
  {3, 3, 0, 4, 0},
  /* wedge */
  {4, 0, 1, 4, 3},
  {4, 1, 2, 5, 4},
  {4, 2, 0, 3, 5},
  {3, 0, 2, 1, 0},
  {3, 3, 4, 5, 0},
  /* hex */
  {4, 0, 3, 2, 1},
  {4, 0, 1, 5, 4},
  {4, 1, 2, 6, 5},
  {4, 2, 3, 7, 6},
  {4, 0, 4, 7, 3},
  {4, 4, 5, 6, 7}
};

  
struct section{
  int id;
  CGNS_ENUMT(ElementType_t) etype;
  string cname;
  cgsize_t start_index;
  cgsize_t end_index;
  int nbndry;
  int parent_flag;
  section(){}
  
  section( int i,
           ElementType_t t,
           string n,
           cgsize_t start,
           cgsize_t end,
           int nb,
           int pf):id(i),etype(t),cname(n),start_index(start),end_index(end),nbndry(nb),parent_flag(pf){}
};

struct sect{ //part of a section
  int id;
  cgsize_t start_index; //file numbering
  cgsize_t end_index;  //file numbering
  sect(int i,
       cgsize_t start,
       cgsize_t end):id(i),start_index(start),end_index(end){}
};
  
// static vector<vector<sect> > redistribute_sections(const vector<section>& sections,//sections in file numbering
//                                                    const vector<cgsize_t>& elemdist,//simple partition of etype
//                                                    CGNS_ENUMT(ElementType_t) my_type)
// {

//   if(Loci::MPI_processes ==1){
//     vector<vector<sect> > dist_section(1);
//     int num_sections = sections.size();
//     for(int i = 0; i < num_sections; i++){
//       if(sections[i].etype == my_type ||sections[i].etype == CGNS_ENUMV(MIXED)){
//         dist_section[0].push_back(sect(sections[i].id, sections[i].start_index,sections[i].end_index));
//       }
//     }
//     return dist_section;
//   }
  
    
  
//   //trasfer file_ids to gloabl_ids
//   cgsize_t global_id = 0;
//   int num_sections = sections.size();
//   vector<section> tmp_sections(sections); //copy of sections
//   for(int i = 0; i < num_sections; i++){
//     if(tmp_sections[i].etype == my_type){
//       cgsize_t start = global_id;
//       cgsize_t end = global_id + tmp_sections[i].end_index - tmp_sections[i].start_index;
//       global_id = end +1;
//       tmp_sections[i].start_index = start;
//       tmp_sections[i].end_index = end;
//     }
//   }
//   //for each section, assign it or part of it to a process
//   int current_p = 0;
//   int P = elemdist.size()-1;
//   vector<vector<sect> > dist_section(P);
//   for(int i = 0; i < num_sections; i++){
//     for(int p = current_p; p < P; p++){
//       if(tmp_sections[i].etype == my_type){
//         cgsize_t global_start = max(tmp_sections[i].start_index,elemdist[p]);
//         cgsize_t global_end = min(tmp_sections[i].end_index,elemdist[p+1]-1);
//         if(global_end >=global_start){
//           cgsize_t file_start = sections[i].start_index + global_start - tmp_sections[i].start_index;
//           cgsize_t file_end = sections[i].start_index + global_end - tmp_sections[i].start_index;
//           dist_section[p].push_back(sect(sections[i].id, file_start, file_end));
          
//         }else{
//           current_p = p+1;
//           break;
//         }
//       }
      
//     }
    
//   }
//   return dist_section;
// }
// /*for debug purpose
//  */

// static void writeASCII(string filename, const store<vector3d<double> > &pos,
//                        const vector<Array<int,5> > &qfaces, const vector<Array<int,4> > &tfaces,
//                        const  vector<Array<int,4> > &tets, const vector<Array<int,5> > &pyramids,
//                        const  vector<Array<int,6> > &prisms,const vector<Array<int,8> > &hexs,
//                        const  vector<pair<int, string> >& bc_ids) {
  
//   std::ofstream ofile(filename.c_str(),std::ios::out) ;
//   entitySet dom = pos.domain();
  
//   ofile << "num_nodes " << dom.size() << endl;
//   ofile.precision(15) ;
//   FORALL(dom, e){
  
//     ofile << pos[e].x << ' ' << pos[e].y << ' ' << pos[e].z
//           << endl ;
//   }ENDFORALL;
//   int num_elem = tfaces.size();
//   ofile << " num trias " << num_elem << endl;
 
//   for(int j=0;j<num_elem;++j){
//     for(int i = 0; i < 4; i++){
//       ofile << tfaces[j][i] << ' ';
//     }
//     ofile << endl ; 
//   }
 
//   num_elem = qfaces.size();
//   ofile << " num quads " << num_elem << endl;
 
//   for(int j=0;j<num_elem;++j){
//     for(int i = 0; i < 5; i++){
//       ofile << qfaces[j][i] << ' ';
//     }
//     ofile << endl ; 
//   }
 
//   num_elem = tets.size();
//   ofile << " num tets " << num_elem << endl;
 
//   for(int j=0;j<num_elem;++j){
//     for(int i = 0; i < 4; i++){
//       ofile << tets[j][i] << ' ';
//     }
//     ofile << endl ; 
//   }
 
//   num_elem = pyramids.size();
//   ofile << " num pyramids " << num_elem << endl;
 
//   for(int j=0;j<num_elem;++j){
//     for(int i = 0; i < 5; i++){
//       ofile << pyramids[j][i] << ' ';
//     }
//     ofile << endl ; 
//   }
 
//   num_elem = prisms.size();
//   ofile << " num prisms " << num_elem << endl;
 
//   for(int j=0;j<num_elem;++j){
//     for(int i = 0; i < 6; i++){
//       ofile << prisms[j][i] << ' ';
//     }
//     ofile << endl ; 
//   }
 
//   num_elem = hexs.size();
//   ofile << " num hexs " << num_elem << endl;
 
//   for(int j=0;j<num_elem;++j){
//     for(int i = 0; i < 8; i++){
//       ofile << hexs[j][i] << ' ';
//     }
//     ofile << endl ; 
//   }
//   ofile.close();
// }

static void error_exit(const char *func)
{
  printf("CGNSlib ERROR:");
  if (func != NULL && *func)
    printf("%s:", func);
  printf("%s\n", cg_get_error());
  exit(1);
}


/*-------------------------------------------------------------------*/

static int compare_faces (void *v1, void *v2)
{
  FACE *f1 = (FACE *)v1;
  FACE *f2 = (FACE *)v2;
  int i, k;
  cgsize_t id, nn, n1[4], n2[4];

  if (f1->nnodes != f2->nnodes)
    return (f1->nnodes - f2->nnodes);

  for (i = 0; i < f1->nnodes; i++) {
    id = f1->nodes[i];
    for (k = 0; k < i; k++) {
      if (n1[k] > id) {
        nn = n1[k];
        n1[k] = id;
        id = nn;
      }
    }
    n1[i] = id;
  }
  for (i = 0; i < f2->nnodes; i++) {
    id = f2->nodes[i];
    for (k = 0; k < i; k++) {
      if (n2[k] > id) {
        nn = n2[k];
        n2[k] = id;
        id = nn;
      }
    }
    n2[i] = id;
  }

  for (i = 0; i < f1->nnodes; i++) {
    if (n1[i] != n2[i])
      return (int)(n1[i] - n2[i]);
  }
  return 0;
}
/*-------------------------------------------------------------------*/
void ListHash(HASH hash, FACE ** Faces){
  HASH_TAB *tabp = (HASH_TAB *)hash;
  int ne = 0;
  for (size_t i = 0; i < tabp->size; i++) {
    for (BUCKET *p = tabp->table[i]; p != NULL; p = p->next){
      Faces[ne++] = (FACE *)(p->entry);
    }
  }
}

/*-------------------------------------------------------------------*/

static size_t hash_face (void *v)
{
  FACE *f = (FACE *)v;
  int n;
  size_t hash = 0;

  for (n = 0; n < f->nnodes; n++)
    hash += (size_t)f->nodes[n];
  return hash;
}



/*-------------------------------------------------------------------*/

static FACE *new_face (cgsize_t faceid)
{
  FACE *face = (FACE *)malloc(sizeof(FACE));
  if (face == NULL)
    error_exit( "malloc failed for a new face");
  face->id = faceid;
  face->bcnum = 0;
  face->nnodes = 0;
  return face;
}

/*--------------------------------------------------------------------*/

static int sort_faces (const void *v1, const void *v2)
{
  FACE **f1 = (FACE **)v1;
  FACE **f2 = (FACE **)v2;

  return (int)((*f1)->id - (*f2)->id);
}

/*--------------------------------------------------------------------*/
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
/*
  convert each element into faces and hash the faces
*/
static void unstructured_elements (HASH facehash, cgsize_t *conn, cgsize_t ne, CGNS_ENUMT(ElementType_t) elemtype)
{
 
 
  int ip = 0, nf = 0;
 
  FACE face, *pf;
  
  if (elemtype < CGNS_ENUMV(TETRA_4)) return;

  CGNS_ENUMT(ElementType_t)  et = elemtype;
  // cout << " et " << et <<  " ne " << ne << endl;
  cgsize_t i = 0;
  for ( cgsize_t n = 0; n < ne; n++) {
    if (elemtype == CGNS_ENUMV(MIXED))
      et = (CGNS_ENUMT(ElementType_t))conn[i++];
    // cout<< " n " << n << endl;
    switch (et) {
    case CGNS_ENUMV(TETRA_4):
   
      ip = 0;
      nf = 4;
      break;
    case CGNS_ENUMV(PYRA_5):
   
      ip = 4;
      nf = 5;
      break;
    case CGNS_ENUMV(PENTA_6):
   
      ip = 9;
      nf = 5;
      break;
    case CGNS_ENUMV(HEXA_8):
   
      ip = 14;
      nf = 6;
      break;
    default:
      nf = 0;
      break;
    }
    for (int j = 0; j < nf; j++) {
      face.nnodes = facenodes[ip+j][0];
      for (int nn = 0; nn < face.nnodes; nn++)
        face.nodes[nn] = conn[i+facenodes[ip+j][nn+1]];
      pf = (FACE *)HashFind(facehash, &face);
      if (NULL == pf) {
        pf = new_face(0);
        pf->nnodes = face.nnodes;
        for (int nn = 0; nn < face.nnodes; nn++)
          pf->nodes[nn] = face.nodes[nn];
        HashAdd(facehash, pf);
        // cout << " face added " << endl;
      }
      else {
        HashDelete(facehash, pf);
        free(pf);
        // cout << " face deleted " << endl;
      }
    }
    int num_nodes = 0;
    cg_npe (et, &num_nodes);
    i += num_nodes;
    // cout << " i " << i << endl;
  }
   
}
/*--------------------------------------------------------------------*/
/* set the bcnum of faces */ 
static void unstructured_boundary (int nb, cgsize_t nFaces, FACE** Faces,
                                   CGNS_ENUMT(PointSetType_t) ptype,
                                   CGNS_ENUMT(GridLocation_t) location,
                                   cgsize_t np, cgsize_t *ptset)
{
  cgsize_t nf, is, ie;
  int n;

  if (ptype == CGNS_ENUMV(PointRange) ||
      ptype == CGNS_ENUMV(PointList)) {
    if (location == CGNS_ENUMV(FaceCenter) ||
        location == CGNS_ENUMV(CellCenter)) {
      ptype = (ptype == CGNS_ENUMV(PointRange) ?
               CGNS_ENUMV(ElementRange) : CGNS_ENUMV(ElementList));
    }
    else if (location != CGNS_ENUMV(Vertex)) {
      return;
    }
  }

  if (ptype == CGNS_ENUMV(PointRange)) {
    if (ptset[0] < ptset[1]) {
      is = ptset[0];
      ie = ptset[1];
    }
    else {
      is = ptset[1];
      ie = ptset[0];
    }
    for (nf = 0; nf < nFaces; nf++) {
      for (n = 0; n < Faces[nf]->nnodes; n++) {
        if (Faces[nf]->nodes[n] < is ||
            Faces[nf]->nodes[n] > ie) break;
      }
      if (n == Faces[nf]->nnodes)
        Faces[nf]->bcnum = nb;
    }
    return;
  }

  if (ptype == CGNS_ENUMV(PointList)) {
    qsort(ptset, (size_t)np, sizeof(cgsize_t), sort_points);
    for (nf = 0; nf < nFaces; nf++) {
      for (n = 0; n < Faces[nf]->nnodes; n++) {
        if (!find_point(Faces[nf]->nodes[n], np, ptset))
          break;
      }
      if (n == Faces[nf]->nnodes)
        Faces[nf]->bcnum = nb;
    }
    return;
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
    for (nf = 0; nf < nFaces; nf++) {
      if (Faces[nf]->id >= is && Faces[nf]->id <= ie)
        Faces[nf]->bcnum = nb;
    }
    return;
  }

  if (ptype == CGNS_ENUMV(ElementList)) {
    qsort(ptset, (size_t)np, sizeof(cgsize_t), sort_points);
    for (nf = 0; nf < nFaces; nf++) {
      if (find_point(Faces[nf]->id, np, ptset))
        Faces[nf]->bcnum = nb;
    }
    return;
  }
}

/*--------------------------------------------------------------------*/
/* read in boundary conditions, set the bcbum for each face, and store the boundary names in bc_ids */
static void boundary_conditions (cgsize_t nFaces,  FACE** Faces, cgsize_t cgFile, cgsize_t cgBase, cgsize_t cgZone,
                                 vector<pair<int,string> >& bc_ids)
{
  int nb, ib, nrmlindex[3];
  cgsize_t is;
  cgsize_t np; //num_points
  cgsize_t *ptset;//point  set
  char name[33];
  CGNS_ENUMT(BCType_t) bctype;
  CGNS_ENUMT(PointSetType_t) ptype;
  CGNS_ENUMT(GridLocation_t) location;
  CGNS_ENUMT(DataType_t) datatype;
  int dim =  1; //for unstructured mesh
  int nBocos = 0;
  if (cg_nbocos (cgFile, cgBase, cgZone, &nBocos))
    error_exit ("cg_nbocos");

  if (nBocos) {
    // Bocos = (BOCO *)malloc(nBocos * sizeof(BOCO));
    // if (Bocos == NULL) error_exit( "malloc failed for bocos");
 
    for (nb = 1; nb <= nBocos; nb++) {
      if (cg_boco_info(cgFile, cgBase, cgZone, nb, name,
                       &bctype, &ptype, &np, nrmlindex, &is,
                       &datatype, &ib))
        error_exit("cg_boco_info");
      //Bocos[nb-1].type = bctype;
      //strcpy(Bocos[nb-1].name, name);
      string bc_name = string(name);
      for (unsigned int i=0; i<bc_name.length(); ++i)
        {
          if(bc_name[i] == ' ') bc_name[i] = '_';
        }
      cout << "bc_name : " << bc_name << endl;
      bc_ids.push_back(pair<int, string>(nb, bc_name));
      if (cg_boco_gridlocation_read(cgFile, cgBase, cgZone,
                                    nb, &location))
        error_exit("cg_boco_gridlocation_read");
      ptset = (cgsize_t *)malloc((size_t)(np * dim) * sizeof(cgsize_t));
      if (ptset == NULL)
        error_exit( "malloc failed for boco ptset");
      if (cg_boco_read(cgFile, cgBase, cgZone, nb, ptset, 0))
        error_exit("cg_boco_read");

      //set up the bcnum of faces
      unstructured_boundary(nb, nFaces, Faces, ptype, location, np, ptset);
      
      free(ptset);
    }
  }
  //why?
  for (np = 0; np < nFaces; np++) {
    if (Faces[np]->bcnum < 0){
      cout<< " negative bcnum" << Faces[np]->bcnum<<endl;
      Faces[np]->bcnum = nBocos - Faces[np]->bcnum;
    }
  }
}

/*read in and hash the faces, set up boundary number and boundary names*/ 
static void unstructured_faces (HASH& facehash, cgsize_t cgFile, cgsize_t cgBase, cgsize_t cgZone,
                                vector<Array<cgsize_t,5> > &qfaces, vector<Array<cgsize_t,4> > &tfaces,
                                vector<pair<int,string> > &surf_ids){

  
  
 
  int nsect = 0, nn = 0, ip = 0, nf = 0;
  CGNS_ENUMT(ElementType_t) elemtype, et;
  char name[33];
  cgsize_t i, n, is, ie, ne;
  cgsize_t size, *conn;
  FACE face, *pf;
    
  if (cg_nsections (cgFile, cgBase, cgZone, &nsect))
    error_exit ("cg_nsections");
  for (int ns = 1; ns <= nsect; ns++) {
    if (cg_section_read (cgFile, cgBase, cgZone, ns,
                         name, &elemtype, &is, &ie, &nn, &ip))
      error_exit ("cg_section_read");
    if (elemtype < CGNS_ENUMV(TRI_3) ||
        (elemtype > CGNS_ENUMV(QUAD_9) &&
         elemtype != CGNS_ENUMV(MIXED))) continue;
    ne = ie - is + 1;
    if (cg_ElementDataSize (cgFile, cgBase, cgZone, ns, &size))
      error_exit ("cg_ElementDataSize");
    conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
    if (conn == NULL)
      error_exit ("malloc failed for element connectivity");
    if (cg_elements_read (cgFile, cgBase, cgZone, ns, conn, NULL))
      error_exit ("cg_elements_read");

    et = elemtype;
    for (i = 0, n = 0; n < ne; n++) {
      if (elemtype == CGNS_ENUMV(MIXED))
        et = (CGNS_ENUMT(ElementType_t))conn[i++];
      switch (et) {
      case CGNS_ENUMV(TRI_3):
      case CGNS_ENUMV(TRI_6):
        nf = 3;
      break;
      case CGNS_ENUMV(QUAD_4):
      case CGNS_ENUMV(QUAD_8):
      case CGNS_ENUMV(QUAD_9):
        nf = 4;
      break;
      default:
        nf = 0;
        break;
      }
      if (nf) {
        face.nnodes = nf;
        for (nn = 0; nn < face.nnodes; nn++)
          face.nodes[nn] = conn[i+nn];
        pf = (FACE *)HashFind(facehash, &face);
        if (NULL == pf) {
          pf = new_face(0);
          pf->nnodes = face.nnodes;
          for (nn = 0; nn < face.nnodes; nn++)
            pf->nodes[nn] = face.nodes[nn];
          HashAdd(facehash, pf);
        }
        pf->id = is + n;
        pf->bcnum = -ns;
      }
      cg_npe (et, &nn);
      i += nn;
    }
    free (conn);
  }

  cgsize_t nFaces = HashSize(facehash);
  FACE ** Faces = (FACE **)malloc(nFaces * sizeof(FACE *));
  if (Faces == NULL)
    error_exit("malloc failed for face list");
  
  ListHash(facehash, Faces);
  HashDestroy(facehash, NULL);
    
  qsort(Faces, nFaces, sizeof(FACE *), sort_faces);

  
  boundary_conditions (nFaces, Faces, cgFile, cgBase, cgZone, surf_ids );

  qfaces.clear();
  tfaces.clear();
 
  for (cgsize_t nf = 0; nf < nFaces; nf++) {
    if (Faces[nf]->nnodes == 4){
      Array<cgsize_t,5> face;
      for(int j = 0; j<4; j++) face[j] = Faces[nf]->nodes[j];
      face[4] = Faces[nf]->bcnum;
      qfaces.push_back(face);
    }else{
      Array<cgsize_t,4> face;
      for(int j = 0; j<3; j++) face[j] = Faces[nf]->nodes[j];
      face[3] = Faces[nf]->bcnum;
      tfaces.push_back(face);
    }
  }
  //clean up
  for(cgsize_t i = 0; i < nFaces; i++){
    free(Faces[i]);
  }
  free(Faces);
}


// /*Not tested yet, 
//   The coordinates are read in as double precision, the format in file can be in single or double precision 
// */
// void readCGNS_parallel(string filename, store<vector3d<double> > &pos,
//                        vector<Array<cgsize_t,5> > &qfaces, vector<Array<cgsize_t,4> > &tfaces,
//                        vector<Array<cgsize_t,4> > &tets, vector<Array<cgsize_t,5> > &pyramids,
//                        vector<Array<cgsize_t,6> > &prisms, vector<Array<cgsize_t,8> > &hexs,
//                        vector<pair<int,string> > &surf_ids ) {

//   char  errmsg[128];


//   pos.allocate(EMPTY) ;
//   qfaces.clear() ;
//   tfaces.clear() ;
//   tets.clear() ;
//   pyramids.clear() ;
//   prisms.clear() ;
//   hexs.clear() ;
  
//   int celldim = 3, phydim = 3;
 
//   //sizes: Number of vertices, cells, and boundary vertices in each (index)-dimension.
//   //
//   // Note that for unstructured grids, the number of cells is the number of highest order elements.
//   //Thus, in three dimensions it's the number of 3-D cells, and in two dimensions it's the number of 2-D cells.

//   //Also for unstructured grids, if the nodes are sorted between internal nodes and boundary nodes,
//   //the optional parameter NBoundVertex must be set equal to the number of boundary nodes.
//   //By default, NBoundVertex equals zero, meaning that the nodes are unsorted.

//   //Note that a non-zero value for NBoundVertex only applies to unstructured grids.
//   //For structured grids, the NBoundVertex parameter always equals 0 in all directions.
  
//   cgsize_t sizes[3]; 
//   cgsize_t start_index, end_index;
//   int nbndry, parent_flag;
//   int num_bases=0, num_zones = 0, num_sections = 0;
//   int index_file = 0, index_base=0, index_zone=0, index_sect = 0;
//   cgsize_t num_nodes=0, num_sf_trias=0, num_sf_quads=0 ;
//   cgsize_t num_vol_tets=0, num_vol_pents5=0, num_vol_pents6=0, num_vol_hexs=0 ;
//   cgsize_t errs = 0;
//   char bname[33], zname[33], cname[33];
 
  
//   CGNS_ENUMT(ZoneType_t) ztype;
//   CGNS_ENUMT(ElementType_t) etype;
//   CGNS_ENUMT(ElementType_t) et;
//   vector<section> sections;

  

//   const int P = MPI_processes ;
//   const int R = MPI_rank ;
//   if(R == 0)cout <<"start reading cgns file" << endl;
//   if(R == 0) {
//     if(cg_open (filename.c_str(), CG_MODE_READ, &index_file)) error_exit(" unable to open CGNS grid file ");
//     if(cg_nbases (index_file, &num_bases))error_exit("error reading number of bases");
//     if(num_bases != 1){
//       cout<<" there are " << num_bases << " bases"<< endl;
//       cout<< "only read the first one" << endl;
//     }
//     index_base =1; //assume only one base and its index is 1
//     if(cg_base_read (index_file, index_base, bname, &celldim, &phydim))error_exit("error reading base information");
//     if(celldim != 3 || phydim != 3){
//       cg_close(index_file);
//       cerr << "only 3D cell and physical dimensions are allowed in CGNS file" << endl;
//       exit(-1);
//     }
//     if(cg_nzones (index_file, index_base, &num_zones)) error_exit("error reading number of zones");
//     if(num_zones != 1){
//       cout<<" there are " << num_zones << " zones"<< endl;
//       cout<< "only read the first one" << endl;
//     }
//     index_zone = 1;//assume only one zone and its index is 1
//     if(cg_zone_type (index_file, index_base, index_zone, &ztype))error_exit("error reading zone type");
//     if (ztype != Unstructured) error_exit("can only handle unstructured grid in CGNS file");
//     if(cg_zone_read (index_file, index_base, index_zone, zname, sizes))error_exit("error reading zone information");
    
//     //3D unstructured sizes:		NVertex, NCell3D, NBoundVertex
//     num_nodes = sizes[0];
//     if(num_nodes < 1){
//       cg_close(index_file);
//       cerr << "number of nodes < 1 " << endl;
//       exit(-1);
//     }
//     if(cg_nsections (index_file, index_base, index_zone,
//                      &num_sections))error_exit("error reading number of sections");
   
//     if(num_sections <1){
//       cg_close(index_file);
//       cerr << "number of section < 1 " << endl;
//       exit(-1);
//     }
//     sections.resize(num_sections);
//     for (index_sect = 1; index_sect <= num_sections; ++index_sect)
//       {
         
//         if(cg_section_read (index_file, index_base, index_zone,
//                             index_sect, cname, &etype,
//                             &start_index, &end_index,
//                             &nbndry, &parent_flag))error_exit("error reading section ");
         
        
//         if (parent_flag != 0)
//           {
//             // cg_close(index_file);
//             // cerr<< "parent data not allowed in CGNS grid file"<<endl;
//             //exit(-1);
//           }
       
//         sections[index_sect-1] = section(index_sect,  etype, string(cname),
//                                          start_index, end_index,
//                                          nbndry,  parent_flag);
//         cgsize_t num_elem = end_index -  start_index + 1;
//         if(etype ==CGNS_ENUMV(MIXED)) {
//           cout<< " Mixed type exists, can only run in serial" << endl;
//           cgsize_t size = 0;
//           if (cg_ElementDataSize (index_file, index_base, index_zone, index_sect, &size))error_exit ("cg_ElementDataSize");
//           cgsize_t * conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
//           if (conn == NULL)error_exit ("memory allocation failed for element connectivity");
//           if (cg_elements_read (index_file, index_base, index_zone,index_sect, conn, NULL))error_exit ("cg_elements_read");
//           int i = 0;
//           for (  int n = 0; n < num_elem; n++) {
//             et = (CGNS_ENUMT(ElementType_t))conn[i++];
            
//             switch (et) {
//             case CGNS_ENUMV(TRI_3):
//             case CGNS_ENUMV(TRI_6):
//               num_sf_trias++;
//             break;
//             case CGNS_ENUMV(QUAD_4):
//             case CGNS_ENUMV(QUAD_8):
//             case CGNS_ENUMV(QUAD_9):
//               num_sf_quads++; 
//             break;
//             case CGNS_ENUMV(TETRA_4):
//             case CGNS_ENUMV(TETRA_10):
//               num_vol_tets++; 
//             break;
//             case CGNS_ENUMV(PYRA_5):
//             case CGNS_ENUMV(PYRA_13):
//             case CGNS_ENUMV(PYRA_14):
//               num_vol_pents5++; 
//             break;
//             case CGNS_ENUMV(PENTA_6):
//             case CGNS_ENUMV(PENTA_15):
//             case CGNS_ENUMV(PENTA_18):
//               num_vol_pents6++; 
//             break;
//             case CGNS_ENUMV(HEXA_8):
//             case CGNS_ENUMV(HEXA_20):
//             case CGNS_ENUMV(HEXA_27):
//               num_vol_hexs++; 
//             break;
//             /* ignore these */
//             case CGNS_ENUMV(NODE):
//             case CGNS_ENUMV(BAR_2):
//             case CGNS_ENUMV(BAR_3):
//               break;
//             /* invalid */
//             default:
//               sprintf(errmsg,
//                       "element type %s not allowed ",
//                       cg_ElementTypeName(et));
//               error_exit( errmsg);
//               break;
//             }
//             int nn = 0;
//             if (cg_npe(et, &nn) || nn <= 0)error_exit("cg_npe");
//             i += nn;
//           }
       
//           free (conn);
//         }else{
            
//           switch(etype){
//           case CGNS_ENUMV(TRI_3):
//           case CGNS_ENUMV(TRI_6):
//             num_sf_trias += num_elem;
//           break;
//           case CGNS_ENUMV(QUAD_4):
//           case CGNS_ENUMV(QUAD_8):
//           case CGNS_ENUMV(QUAD_9):
//             num_sf_quads += num_elem;
//           break;
//           case CGNS_ENUMV(TETRA_4):
//           case CGNS_ENUMV(TETRA_10):
//             num_vol_tets += num_elem;
//           break;
//           case CGNS_ENUMV(PYRA_5):
//           case CGNS_ENUMV(PYRA_13):
//           case CGNS_ENUMV(PYRA_14):
//             num_vol_pents5 += num_elem;
//           break;
//           case CGNS_ENUMV(PENTA_6):
//           case CGNS_ENUMV(PENTA_15):
//           case CGNS_ENUMV(PENTA_18):
//             num_vol_pents6 += num_elem;
//           break;
//           case CGNS_ENUMV(HEXA_8):
//           case CGNS_ENUMV(HEXA_20):
//           case CGNS_ENUMV(HEXA_27):
//             num_vol_hexs += num_elem;
//           break;
//           /* ignore these */
//           case CGNS_ENUMV(NODE):
//           case CGNS_ENUMV(BAR_2):
//           case CGNS_ENUMV(BAR_3):
//             break;
//           /* invalid */
//           default:
//             sprintf(errmsg,
//                     "element type %s not allowed",
//                     cg_ElementTypeName(etype));
//             error_exit(errmsg);
//             break;
//           }
//         }
//       }
  
  
                
      

    
//     if( num_vol_tets+num_vol_pents5+num_vol_pents6+num_vol_hexs+num_sf_trias+ num_sf_quads != sizes[1]){
//       cerr<<" total elements: " <<num_vol_tets+num_vol_pents5+num_vol_pents6+num_vol_hexs << " sizes[1] " <<  sizes[1] << endl;
//       error_exit("number of elements and faces does not match in CGNS grid file ");
//     }
    
//     cout << "nnodes=" << num_nodes <<",ntria="<<num_sf_trias
//          << ",nquad="<<num_sf_quads<<",ntets="<<num_vol_tets
//          << ",npyrm="<<num_vol_pents5<<",nprsm="<<num_vol_pents6
//          << ",nhex="<<num_vol_hexs << endl ;

//     cout << "sections: " << endl;
//     for(unsigned int i = 0; i < sections.size(); i++){
//       cout << "section " << sections[i].id <<" : "<<sections[i].cname<< "  etype " << sections[i].etype
//            << "  start " << sections[i].start_index << " end " <<  sections[i].end_index<<endl;
//     }
//   }

//   Array<cgsize_t,7> data ;
//   data[0] = num_nodes ;
//   data[1] = num_sf_trias ;
//   data[2] = num_sf_quads ;
//   data[3] = num_vol_tets ;
//   data[4] = num_vol_pents5 ;
//   data[5] = num_vol_pents6 ;
//   data[6] = num_vol_hexs ;
//   MPI_Bcast(&data[0],7*sizeof(cgsize_t),MPI_BYTE,0,MPI_COMM_WORLD) ;
//   num_nodes      = data[0] ;
//   num_sf_trias   = data[1] ;
//   num_sf_quads   = data[2] ;
//   num_vol_tets   = data[3] ;
//   num_vol_pents5 = data[4] ;
//   num_vol_pents6 = data[5] ;
//   num_vol_hexs   = data[6] ;

//   vector<cgsize_t> node_ptns = Loci::g_simple_partition_vec<cgsize_t>(0,num_nodes-1,P) ;
//   vector<gEntitySet> local_nodes(P) ;
//   for(cgsize_t i=0;i<P;++i)
//     local_nodes[i] = gEntitySet(std::pair<gEntity, gEntity>(static_cast<gEntity>(node_ptns[i]),static_cast<gEntity>(node_ptns[i+1]-1))) ;
  
//   size_t mxsize = max(local_nodes[0].size(),
//                       local_nodes[P-1].size()) ;
  
//   pos.allocate(local_nodes[R]) ;
//   vector<double> buf(mxsize*3) ;
//   if(num_nodes > 0){ // Read in positions
//     if(R == 0) { 
   
   
    
//       start_index = local_nodes[0].Min()+1;
//       end_index = local_nodes[0].Max()+1;
    
    
//       errs = cg_coord_read (index_file, index_base, index_zone,
//                             "CoordinateX", RealDouble,
//                             &start_index, &end_index, &buf[0]);
//       errs = cg_coord_read (index_file, index_base, index_zone,
//                             "CoordinateY", RealDouble,
//                             &start_index, &end_index, &buf[mxsize]);
    
//       errs = cg_coord_read (index_file, index_base, index_zone,
//                             "CoordinateZ", RealDouble,
//                             &start_index, &end_index, &buf[2*mxsize]); 
//       cgsize_t index = 0;
//       FORALL(local_nodes[0],nd) {
//         pos[nd].x = buf[index++];
//       } ENDFORALL ;
//       index = 0;
//       FORALL(local_nodes[0],nd) {
//         pos[nd].y = buf[mxsize+index];
//         index++;
//       } ENDFORALL ;
     
//       index = 0;
//       FORALL(local_nodes[0],nd) {
//         pos[nd].z = buf[2*mxsize+index];
//         index++;
//       } ENDFORALL ;

//       //ElementRange contains the index of the first and last elements defined
//       //in ElementConnectivity. The elements are indexed with a global numbering
//       //system, starting at 1, for all element sections under a given Zone_t data
//       //structure. The global numbering insures that each element, whether it's a
//       //cell, a face, or an edge, is uniquely identified by its number. They are
//       //also listed as a continuous list of element numbers within any single
//       //element section. Therefore the number of elements in a section is:

//       //ElementSize = ElementRange.end - ElementRange.start + 1

     
    
    

//       for(cgsize_t i=1;i<P;++i) {
//         start_index = local_nodes[i].Min()+1;
//         end_index = local_nodes[i].Max()+1;
//         errs = cg_coord_read (index_file, index_base, index_zone,
//                               "CoordinateX", RealDouble,
//                               &start_index, &end_index, &buf[0]);
//         errs = cg_coord_read (index_file, index_base, index_zone,
//                               "CoordinateY", RealDouble,
//                               &start_index, &end_index, &buf[mxsize]);
      
//         errs = cg_coord_read (index_file, index_base, index_zone,
//                               "CoordinateY", RealDouble,
//                               &start_index, &end_index, &buf[2*mxsize]); 
      
//         MPI_Send(&buf[0],mxsize*3,MPI_DOUBLE,i,10,MPI_COMM_WORLD) ;
//       }
//     } else {
//       MPI_Status status ;
        
//       MPI_Recv(&buf[0],3*mxsize,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status) ;
    
//       cgsize_t index = 0;
//       FORALL(local_nodes[R],nd) {
//         pos[nd].x = buf[index++];
//       } ENDFORALL ;
    
//       index = 0;
//       FORALL(local_nodes[R],nd) {
//         pos[nd].y = buf[mxsize+index];
//         index++;
//       } ENDFORALL ;
     
//       index = 0;
//       FORALL(local_nodes[R],nd) {
//         pos[nd].z = buf[2*mxsize+index];
//         index++;
//       } ENDFORALL ;
//     }
//   }
  
 
//   cgsize_t *Parent_Data = NULL;
 
  
//   if(R == 0){

//     //create a a hash  
//     cgsize_t tot_vol_elem = num_vol_tets+num_vol_pents5+num_vol_pents6+num_vol_hexs;
//     HASH facehash;
//     facehash = HashCreate(tot_vol_elem > 1024 ? (size_t)tot_vol_elem / 3 : 127,
//                           compare_faces, hash_face);
//     if (NULL == facehash)
//       error_exit("hash table creation failed");

    
//     cout <<"start reading volume" << endl;
    
//     // Read in volume elements

//     // Read in tetrahedra
//     vector<cgsize_t> tetsdist = Loci::g_simple_partition_vec<cgsize_t>(0,num_vol_tets-1,P) ;
//     tets = vector<Array<cgsize_t,4> >(tetsdist[R+1]-tetsdist[R]) ;

//     if(R == 0) {

//       //first read in my own
//       cgsize_t tsz = max(tetsdist[1]-tetsdist[0],tetsdist[P]-tetsdist[P-1]) ;
//       vector<Array<cgsize_t,4> > send_buf(tsz) ;

//       cgsize_t local_id = 0;
//       vector<vector<sect> > sect_dist = redistribute_sections(sections,//sections in file numbering
//                                                               tetsdist,//simple partition of etype
//                                                               CGNS_ENUMV(TETRA_4));
         
         
//       for(unsigned int secti = 0; secti < sect_dist[0].size(); secti++){
        
//         if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                      sect_dist[0][secti].id, sect_dist[0][secti].start_index ,
//                                      sect_dist[0][secti].end_index,
//                                      &tets[local_id][0], Parent_Data))error_exit("error reading element information");
           
//         //hash all the faces of tets that just been read in 
     
//         unstructured_elements(facehash, &tets[local_id][0], sect_dist[0][secti].end_index -sect_dist[0][secti].start_index+1 , CGNS_ENUMV(TETRA_4));
     

      
//         local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//       }
          
//       //read for other processes
//       for(cgsize_t p=1;p<P;++p) {
//         local_id = 0;
//         for(unsigned int secti = 0; secti < sect_dist[p].size(); secti++){
//           if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                        sect_dist[p][secti].id, sect_dist[p][secti].start_index ,
//                                        sect_dist[p][secti].end_index,
//                                        &send_buf[local_id][0], Parent_Data))error_exit("error reading element information");
         
//           //hash all the faces of tets that just been read in 
        
           
//           unstructured_elements(facehash, &send_buf[local_id][0], sect_dist[p][secti].end_index -sect_dist[p][secti].start_index+1 , CGNS_ENUMV(TETRA_4));
           
        
 
       
//           local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//         }
           
//         cgsize_t ltsz = tetsdist[p+1]-tetsdist[p] ;
//         if(ltsz != 0)
//           MPI_Send(&send_buf[0][0],ltsz*4*sizeof(cgsize_t),MPI_BYTE,p,7,MPI_COMM_WORLD) ;
//       }
//     } else {
//       MPI_Status status ;
//       cgsize_t tsz = tets.size() ;
//       if(tsz != 0)
//         MPI_Recv(&tets[0][0],tsz*4*sizeof(cgsize_t),MPI_BYTE,0,7,MPI_COMM_WORLD,&status) ;
//     }
//     cout<< " end tets" << endl;
  
//     // Read in pyramids
//     vector<cgsize_t> pyrmdist = Loci::g_simple_partition_vec<cgsize_t>(0,num_vol_pents5-1,P) ;
//     pyramids = vector<Array<cgsize_t,5> >(pyrmdist[R+1]-pyrmdist[R]) ;
  
//     if(R == 0) {
//       cgsize_t tsz = max(pyrmdist[1]-pyrmdist[0],pyrmdist[P]-pyrmdist[P-1]) ;
//       vector<Array<cgsize_t,5> > send_buf(tsz) ;
    
//       cgsize_t local_id = 0;
//       vector<vector<sect> > sect_dist = redistribute_sections(sections,//sections in file numbering
//                                                               pyrmdist,//simple partition of etype
//                                                               CGNS_ENUMV(PYRA_5));
         
         
//       for(unsigned int secti = 0; secti < sect_dist[0].size(); secti++){
//         if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                      sect_dist[0][secti].id, sect_dist[0][secti].start_index ,
//                                      sect_dist[0][secti].end_index,
//                                      &pyramids[local_id][0], Parent_Data))error_exit("error reading element information");
           
           
//         //hash all the faces of pyramid that just been read in 
      
//         unstructured_elements(facehash, &pyramids[local_id][0], sect_dist[0][secti].end_index -sect_dist[0][secti].start_index+1 , CGNS_ENUMV(PYRA_5));
      
      
//         local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//       }
          
//       //read for other processes
//       for(int p=1;p<P;++p) {
//         local_id = 0;
//         for(unsigned int secti = 0; secti < sect_dist[p].size(); secti++){
//           if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                        sect_dist[p][secti].id, sect_dist[p][secti].start_index ,
//                                        sect_dist[p][secti].end_index,
//                                        &send_buf[local_id][0], Parent_Data))error_exit("error reading element information");
         
          
        
//           //hash all the faces of tets that just been read in 
         
//           unstructured_elements(facehash, &send_buf[local_id][0], sect_dist[p][secti].end_index -sect_dist[p][secti].start_index+1 , CGNS_ENUMV(PYRA_5));
         
        
//           local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//         }
    
//         cgsize_t ltsz = pyrmdist[p+1]-pyrmdist[p] ;
//         if(ltsz != 0)
//           MPI_Send(&send_buf[0][0],ltsz*5*sizeof(cgsize_t),MPI_BYTE,p,6,MPI_COMM_WORLD) ;
//       }
//     } else {
//       MPI_Status status ;
//       cgsize_t tsz = pyramids.size() ;
//       if(tsz != 0)
//         MPI_Recv(&pyramids[0][0],tsz*5*sizeof(cgsize_t),MPI_BYTE,0,6,MPI_COMM_WORLD,&status) ;
//     }

   
  
//     // Read in prisms
//     vector<cgsize_t> prsmdist = Loci::g_simple_partition_vec<cgsize_t>(0,num_vol_pents6-1,P) ;
//     prisms = vector<Array<cgsize_t,6> >(prsmdist[R+1]-prsmdist[R]) ;

//     if(R == 0) {

//       cgsize_t tsz = max(prsmdist[1]-prsmdist[0],prsmdist[P]-prsmdist[P-1]) ;
//       vector<Array<cgsize_t,6> > send_buf(tsz) ;
    
//       cgsize_t local_id = 0;
//       vector<vector<sect> > sect_dist = redistribute_sections(sections,//sections in file numbering
//                                                               prsmdist,//simple partition of etype
//                                                               CGNS_ENUMV(PENTA_6));
         
         
//       for(unsigned int secti = 0; secti < sect_dist[0].size(); secti++){
//         if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                      sect_dist[0][secti].id, sect_dist[0][secti].start_index ,
//                                      sect_dist[0][secti].end_index,
//                                      &prisms[local_id][0], Parent_Data))error_exit("error reading element information" );
           
//         //hash all the faces of pyramid that just been read in 
       
//         unstructured_elements(facehash, &prisms[local_id][0], sect_dist[0][secti].end_index -sect_dist[0][secti].start_index+1 , CGNS_ENUMV(PENTA_6));
       
     
//         local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//       }
          
//       //read for other processes
//       for(int p=1;p<P;++p) {
//         local_id = 0;
//         for(unsigned int secti = 0; secti < sect_dist[p].size(); secti++){
//           if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                        sect_dist[p][secti].id, sect_dist[p][secti].start_index ,
//                                        sect_dist[p][secti].end_index,
//                                        &send_buf[local_id][0], Parent_Data))error_exit("error reading element information");
         
          
//           //hash all the faces of tets that just been read in 
       
//           unstructured_elements(facehash, &send_buf[local_id][0], sect_dist[p][secti].end_index -sect_dist[p][secti].start_index+1 , CGNS_ENUMV(PENTA_6));
       
        
        
//           local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//         }
      
//         cgsize_t ltsz = prsmdist[p+1]-prsmdist[p] ;
//         if(ltsz != 0)
//           MPI_Send(&send_buf[0][0],ltsz*6*sizeof(cgsize_t),MPI_BYTE,p,5,MPI_COMM_WORLD) ;
//       }
//     } else {
//       MPI_Status status ;
//       cgsize_t tsz = prisms.size() ;
//       if(tsz != 0)
//         MPI_Recv(&prisms[0][0],tsz*6*sizeof(cgsize_t),MPI_BYTE,0,5,MPI_COMM_WORLD,&status) ;
//     }

//     // Read in hexahdra
//     vector<cgsize_t> hexsdist = Loci::g_simple_partition_vec<cgsize_t>(0,num_vol_hexs-1,P) ;
//     hexs = vector<Array<cgsize_t,8> >(hexsdist[R+1]-hexsdist[R]) ;

//     if(R == 0) {
//       cgsize_t tsz = max(hexsdist[1]-hexsdist[0],hexsdist[P]-hexsdist[P-1]) ;
//       vector<Array<cgsize_t,8> > send_buf(tsz) ;
    
//       cgsize_t local_id = 0;
//       vector<vector<sect> > sect_dist = redistribute_sections(sections,//sections in file numbering
//                                                               hexsdist,//simple partition of etype
//                                                               CGNS_ENUMV(HEXA_8));
         
         
//       for(unsigned int secti = 0; secti < sect_dist[0].size(); secti++){
//         if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                      sect_dist[0][secti].id, sect_dist[0][secti].start_index ,
//                                      sect_dist[0][secti].end_index,
//                                      &hexs[local_id][0], Parent_Data))error_exit("error reading element information");
           
//         //hash all the faces of pyramid that just been read in 
       
//         unstructured_elements(facehash, &hexs[local_id][0], sect_dist[0][secti].end_index -sect_dist[0][secti].start_index+1 , CGNS_ENUMV(HEXA_8));
       
     
      
//         local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//       }
          
//       //read for other processes
//       for(int p=1;p<P;++p) {
//         local_id = 0;
//         for(unsigned int secti = 0; secti < sect_dist[p].size(); secti++){
//           if(cg_elements_partial_read (index_file, index_base, index_zone,
//                                        sect_dist[p][secti].id, sect_dist[p][secti].start_index ,
//                                        sect_dist[p][secti].end_index,
//                                        &send_buf[local_id][0], Parent_Data))error_exit("error reading element information");
         
//           //hash all the faces of tets that just been read in 
       
//           unstructured_elements(facehash, &send_buf[local_id][0], sect_dist[p][secti].end_index -sect_dist[p][secti].start_index+1 , CGNS_ENUMV(HEXA_8));
       
         
       
//           local_id += sect_dist[0][secti].end_index - sect_dist[0][secti].start_index + 1 ;
//         }
//         cgsize_t ltsz = hexsdist[p+1]-hexsdist[p] ;
//         if(ltsz != 0)
//           MPI_Send(&send_buf[0][0],ltsz*8*sizeof(cgsize_t),MPI_BYTE,p,4,MPI_COMM_WORLD) ;
//       }
//     } else {
//       MPI_Status status ;
//       cgsize_t tsz = hexs.size() ;
//       if(tsz != 0)
//         MPI_Recv(&hexs[0][0],tsz*8*sizeof(cgsize_t),MPI_BYTE,0,4,MPI_COMM_WORLD,&status) ;
//     }


    
  
//     //read in trias and quads, then hash and sort faces, get boundary information
//     if(R==0){
//       cout << "start  unstructured_faces " << endl;
//       unstructured_faces (facehash, index_file, index_base, index_zone,
//                           qfaces, tfaces,
//                           surf_ids);

//       cout << "end  unstructured_faces " << endl; 
//     }
//     num_sf_trias  = tfaces.size();
//     num_sf_quads = qfaces.size();
//     {  
//       Array<cgsize_t,2> data ;
//       data[0] = num_sf_trias ;
//       data[1] = num_sf_quads ;
//       MPI_Bcast(&data[0],2*sizeof(cgsize_t),MPI_BYTE,0,MPI_COMM_WORLD) ;
//       num_sf_trias   = data[0] ;
//       num_sf_quads   = data[1] ;
//     }

//     // triangles
//     vector<cgsize_t> triadist = Loci::g_simple_partition_vec<cgsize_t>(0,num_sf_trias-1,P) ;
//     vector<cgsize_t> quaddist = Loci::g_simple_partition_vec<cgsize_t>(0,num_sf_quads-1,P) ;
//     if(R != 0){
//       qfaces = vector<Array<cgsize_t,5> >(quaddist[R+1]-quaddist[R]) ;
//       tfaces = vector<Array<cgsize_t,4> >(triadist[R+1]-triadist[R]) ;
//     }
//     if(R==0){
//       if(num_sf_trias > 0) {
      
//         for(int p=1;p<P;++p) {
//           cgsize_t ltsz = triadist[p+1]-triadist[p] ;
//           if(ltsz != 0)
//             MPI_Send(&tfaces[triadist[p]][0],ltsz*4*sizeof(cgsize_t),MPI_BYTE,p,9,MPI_COMM_WORLD) ;
//         }
//       }
      
//       if(num_sf_quads > 0){//quads
      
       
//         for(int p=1;p<P;++p) {
         
           
//           cgsize_t ltsz = quaddist[p+1]-quaddist[p] ;
//           if(ltsz != 0)
//             MPI_Send(&qfaces[quaddist[p]][0],ltsz*5*sizeof(cgsize_t),MPI_BYTE,p,8,MPI_COMM_WORLD) ;
//         }
//       }
          
//     } else {
//       MPI_Status status ;
//       cgsize_t tsz = tfaces.size() ;
//       if(tsz != 0)
//         MPI_Recv(&tfaces[0][0],tsz*4*sizeof(cgsize_t),MPI_BYTE,0,9,MPI_COMM_WORLD,&status) ;

//       cgsize_t qsz = qfaces.size();
//       if(qsz != 0)
//         MPI_Recv(&qfaces[0][0],qsz*5*sizeof(cgsize_t),MPI_BYTE,0,8,MPI_COMM_WORLD,&status) ;
//     }
 
//     if(R==0){
//       tfaces.resize(triadist[1]-triadist[0]);
//       qfaces.resize(quaddist[1]-quaddist[0]);
//     }
   
//     if(R == 0)cout <<"finish reading cgns file" << endl;
//   }
// }
/*read from cgns file,
  The coordinates are read in as double precision, the format in file can be in single or double precision 
*/
void readCGNS_serial(string zonepath,  store<vector3d<double> > &pos,
                     vector<Array<cgsize_t,5> > &qfaces, vector<Array<cgsize_t,4> > &tfaces,
                     vector<Array<cgsize_t,4> > &tets, vector<Array<cgsize_t,5> > &pyramids,
                     vector<Array<cgsize_t,6> > &prisms, vector<Array<cgsize_t,8> > &hexs,
                     vector<pair<int,string> > &surf_ids ) {

  char  errmsg[128];


  pos.allocate(EMPTY) ;
  qfaces.clear() ;
  tfaces.clear() ;
  tets.clear() ;
  pyramids.clear() ;
  prisms.clear() ;
  hexs.clear() ;
  
  int celldim = 3, phydim = 3;
 
  //sizes: Number of vertices, cells, and boundary vertices in each (index)-dimension.
  //
  // Note that for unstructured grids, the number of cells is the number of highest order elements.
  //Thus, in three dimensions it's the number of 3-D cells, and in two dimensions it's the number of 2-D cells.

  //Also for unstructured grids, if the nodes are sorted between internal nodes and boundary nodes,
  //the optional parameter NBoundVertex must be set equal to the number of boundary nodes.
  //By default, NBoundVertex equals zero, meaning that the nodes are unsorted.

  //Note that a non-zero value for NBoundVertex only applies to unstructured grids.
  //For structured grids, the NBoundVertex parameter always equals 0 in all directions.
  
  cgsize_t sizes[3]; 
  cgsize_t start_index, end_index;
  int nbndry, parent_flag;
  int num_bases=0, num_zones = 0, num_sections = 0;
  int index_file = 0, index_base=0, index_zone=0, index_sect = 0;
  cgsize_t num_nodes=0, num_sf_trias=0, num_sf_quads=0 ;
  cgsize_t num_vol_tets=0, num_vol_pents5=0, num_vol_pents6=0, num_vol_hexs=0 ;
  cgsize_t errs = 0;
  char bname[33], zname[33], cname[33];
 
  
  CGNS_ENUMT(ZoneType_t) ztype;
  CGNS_ENUMT(ElementType_t) etype;
  
 

  std::size_t p1 = zonepath.find('/');
  std::size_t p2 = zonepath.find('/', p1+1);

  if(p1 == string::npos){
    index_base = 1;
    index_zone = 1;
  }else if(p2 == string::npos){
    index_zone = 1;
  }
  string filename = zonepath.substr(0, p1);
  cout << " zonepath " << zonepath << endl; 
  cout << " filename " << filename << endl;
  cout <<"start reading cgns file" << endl;
  if(cg_open (filename.c_str(), CG_MODE_READ, &index_file)) error_exit(" unable to open CGNS grid file ");
  if(!index_base){
    string basename = zonepath.substr(p1+1, p2-p1-1);
    cout << " basename " << basename << endl;
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
  }
  if(!index_zone){
    string zonename = zonepath.substr(p2+1);
    cout << " zonename " << zonename << endl;
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
  if(cg_nsections (index_file, index_base, index_zone,
                   &num_sections))error_exit("error reading number of sections");
   
  if(num_sections <1){
    cg_close(index_file);
    cerr << "number of section < 1 " << endl;
    exit(-1);
  }
   
  for (index_sect = 1; index_sect <= num_sections; ++index_sect)
    {
         
      if(cg_section_read (index_file, index_base, index_zone,
                          index_sect, cname, &etype,
                          &start_index, &end_index,
                          &nbndry, &parent_flag))error_exit("error reading section ");
         
        
        
      cgsize_t num_elem = end_index -  start_index + 1;
      if(etype ==CGNS_ENUMV(MIXED)) {
        cout<< " Mixed type exists, can only run in serial" << endl;
        cgsize_t size = 0;
        if (cg_ElementDataSize (index_file, index_base, index_zone, index_sect, &size))error_exit ("cg_ElementDataSize");
        cgsize_t * conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
        if (conn == NULL)error_exit ("memory allocation failed for element connectivity");
        if (cg_elements_read (index_file, index_base, index_zone,index_sect, conn, NULL))error_exit ("cg_elements_read");
        int i = 0;
        for (  int n = 0; n < num_elem; n++) {
          CGNS_ENUMT(ElementType_t)  et = (CGNS_ENUMT(ElementType_t))conn[i++];
         
          switch (et) {
          case CGNS_ENUMV(TRI_3):
          case CGNS_ENUMV(TRI_6):
            num_sf_trias++;
          break;
          case CGNS_ENUMV(QUAD_4):
          case CGNS_ENUMV(QUAD_8):
          case CGNS_ENUMV(QUAD_9):
            num_sf_quads++; 
          break;
          case CGNS_ENUMV(TETRA_4):
          case CGNS_ENUMV(TETRA_10):
            num_vol_tets++; 
          break;
          case CGNS_ENUMV(PYRA_5):
          case CGNS_ENUMV(PYRA_13):
          case CGNS_ENUMV(PYRA_14):
            num_vol_pents5++; 
          break;
          case CGNS_ENUMV(PENTA_6):
          case CGNS_ENUMV(PENTA_15):
          case CGNS_ENUMV(PENTA_18):
            num_vol_pents6++; 
          break;
          case CGNS_ENUMV(HEXA_8):
          case CGNS_ENUMV(HEXA_20):
          case CGNS_ENUMV(HEXA_27):
            num_vol_hexs++; 
          break;
          /* ignore these */
          case CGNS_ENUMV(NODE):
          case CGNS_ENUMV(BAR_2):
          case CGNS_ENUMV(BAR_3):
            break;
          /* invalid */
          default:
            sprintf(errmsg,
                    "element type %s not allowed",
                    cg_ElementTypeName(et));
            error_exit( errmsg);
            break;
          }
          int nn=0;
          if (cg_npe(et, &nn) || nn <= 0)
            error_exit("cg_npe");
          i += nn;
        }
       
        free (conn);
      }else{
         
        switch(etype){
        case CGNS_ENUMV(TRI_3):
        case CGNS_ENUMV(TRI_6):
          num_sf_trias += num_elem;
        break;
        case CGNS_ENUMV(QUAD_4):
        case CGNS_ENUMV(QUAD_8):
        case CGNS_ENUMV(QUAD_9):
          num_sf_quads += num_elem;
        break;
        case CGNS_ENUMV(TETRA_4):
        case CGNS_ENUMV(TETRA_10):
          num_vol_tets += num_elem;
        break;
        case CGNS_ENUMV(PYRA_5):
        case CGNS_ENUMV(PYRA_13):
        case CGNS_ENUMV(PYRA_14):
          num_vol_pents5 += num_elem;
        break;
        case CGNS_ENUMV(PENTA_6):
        case CGNS_ENUMV(PENTA_15):
        case CGNS_ENUMV(PENTA_18):
          num_vol_pents6 += num_elem;
        break;
        case CGNS_ENUMV(HEXA_8):
        case CGNS_ENUMV(HEXA_20):
        case CGNS_ENUMV(HEXA_27):
          num_vol_hexs += num_elem;
        break;
        /* ignore these */
        case CGNS_ENUMV(NODE):
        case CGNS_ENUMV(BAR_2):
        case CGNS_ENUMV(BAR_3):
          break;
        /* invalid */
        default:
          sprintf(errmsg,
                  "element type %s not allowed",
                  cg_ElementTypeName(etype));
          error_exit(errmsg);
          break;
        }
      }
    }
  
  
                
      

    
  // if( num_vol_tets+num_vol_pents5+num_vol_pents6+num_vol_hexs + num_sf_trias + num_sf_quads != sizes[1]){
  cout<<" total elements: " <<num_vol_tets+num_vol_pents5+num_vol_pents6+num_vol_hexs << " sizes[1] " <<  sizes[1] << endl;
  //error_exit("number of elements and faces does not match in CGNS grid file ");
  //}
    
  cout << "nnodes=" << num_nodes <<",ntria="<<num_sf_trias
       << ",nquad="<<num_sf_quads<<",ntets="<<num_vol_tets
       << ",npyrm="<<num_vol_pents5<<",nprsm="<<num_vol_pents6
       << ",nhex="<<num_vol_hexs << endl ;


 

  entitySet local_nodes = interval(0, num_nodes-1);
  pos.allocate(local_nodes) ;
  cgsize_t mxsize = num_nodes;
  vector<double> buf(mxsize*3) ;
  if(num_nodes > 0){ // Read in positions
   
    start_index = local_nodes.Min()+1;
    end_index = local_nodes.Max()+1;
        
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateX", RealDouble,
                          &start_index, &end_index, &buf[0]);
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateY", RealDouble,
                          &start_index, &end_index, &buf[mxsize]);
    
    errs = cg_coord_read (index_file, index_base, index_zone,
                          "CoordinateZ", RealDouble,
                          &start_index, &end_index, &buf[2*mxsize]); 
    cgsize_t index = 0;
    FORALL(local_nodes,nd) {
      pos[nd].x = buf[index++];
    } ENDFORALL ;
    index = 0;
    FORALL(local_nodes,nd) {
      pos[nd].y = buf[mxsize+index];
      index++;
    } ENDFORALL ;
     
    index = 0;
    FORALL(local_nodes,nd) {
      pos[nd].z = buf[2*mxsize+index];
      index++;
    } ENDFORALL ;

    //ElementRange contains the index of the first and last elements defined
    //in ElementConnectivity. The elements are indexed with a global numbering
    //system, starting at 1, for all element sections under a given Zone_t data
    //structure. The global numbering insures that each element, whether it's a
    //cell, a face, or an edge, is uniquely identified by its number. They are
    //also listed as a continuous list of element numbers within any single
    //element section. Therefore the number of elements in a section is:

    //ElementSize = ElementRange.end - ElementRange.start + 1

     
    
    

     
  }
  
 

 
  
  

  //create a a hash  
  cgsize_t tot_vol_elem = num_vol_tets+num_vol_pents5+num_vol_pents6+num_vol_hexs;
  HASH facehash;
  facehash = HashCreate(tot_vol_elem > 1024 ? (size_t)tot_vol_elem / 3 : 127,
                        compare_faces, hash_face);
  if (NULL == facehash)
    error_exit("hash table creation failed");
      
      
  cout <<"start reading volume" << endl;
    
  // Read in volume elements
      
     
  tets = vector<Array<cgsize_t,4> >(num_vol_tets) ;
  cgsize_t tets_id = 0;
  hexs = vector<Array<cgsize_t,8> >(num_vol_hexs) ;
  cgsize_t hexs_id = 0;
  prisms = vector<Array<cgsize_t,6> >(num_vol_pents6) ;
  cgsize_t prisms_id = 0;
  pyramids = vector<Array<cgsize_t,5> >(num_vol_pents5) ;
  cgsize_t pyramids_id = 0;




  for (index_sect = 1; index_sect <= num_sections; ++index_sect)
    {
         
      if(cg_section_read (index_file, index_base, index_zone,
                          index_sect, cname, &etype,
                          &start_index, &end_index,
                          &nbndry, &parent_flag))error_exit("error reading section ");
      cgsize_t num_elem = end_index -  start_index + 1;
      if(etype ==CGNS_ENUMV(MIXED)) {
        cgsize_t size = 0;
        if (cg_ElementDataSize (index_file, index_base, index_zone, index_sect, &size))error_exit ("cg_ElementDataSize");
        cgsize_t * conn = (cgsize_t *) malloc ((size_t)size * sizeof(cgsize_t));
        if (conn == NULL)error_exit ("memory allocation failed for element connectivity");
            
        if (cg_elements_read (index_file, index_base, index_zone,index_sect, conn, NULL))error_exit ("cg_elements_read");
        int i = 0;
        for (  int n = 0; n < num_elem; n++) {
          CGNS_ENUMT(ElementType_t)  et = (CGNS_ENUMT(ElementType_t))conn[i++];
          switch (et) {
          case CGNS_ENUMV(TRI_3):
          case CGNS_ENUMV(TRI_6):
                
            break;
          case CGNS_ENUMV(QUAD_4):
          case CGNS_ENUMV(QUAD_8):
          case CGNS_ENUMV(QUAD_9):
             
            break;
          case CGNS_ENUMV(TETRA_4):
          case CGNS_ENUMV(TETRA_10):
            tets[tets_id][0] = conn[i+0];
          tets[tets_id][1] = conn[i+1];
          tets[tets_id][2] = conn[i+2];
          tets[tets_id][3] = conn[i+3];
          unstructured_elements(facehash, &tets[tets_id][0], 1 , CGNS_ENUMV(TETRA_4)); 
          tets_id++;
          break;
          case CGNS_ENUMV(PYRA_5):
          case CGNS_ENUMV(PYRA_13):
          case CGNS_ENUMV(PYRA_14):
            pyramids[pyramids_id][0] = conn[i+0];
          pyramids[pyramids_id][1] = conn[i+1];
          pyramids[pyramids_id][2] = conn[i+2];
          pyramids[pyramids_id][3] = conn[i+3];
          pyramids[pyramids_id][4] = conn[i+4];
          unstructured_elements(facehash, &pyramids[pyramids_id][0], 1 , CGNS_ENUMV(PYRA_5)); 
          pyramids_id++;
          break;
          case CGNS_ENUMV(PENTA_6):
          case CGNS_ENUMV(PENTA_15):
          case CGNS_ENUMV(PENTA_18):
            prisms[prisms_id][0] = conn[i+0];
          prisms[prisms_id][1] = conn[i+1];
          prisms[prisms_id][2] = conn[i+2];
          prisms[prisms_id][3] = conn[i+3];
          prisms[prisms_id][4] = conn[i+4];
          prisms[prisms_id][5] = conn[i+5];
          unstructured_elements(facehash, &prisms[prisms_id][0], 1 , CGNS_ENUMV(PENTA_6)); 
          prisms_id++;
          break;
          case CGNS_ENUMV(HEXA_8):
          case CGNS_ENUMV(HEXA_20):
          case CGNS_ENUMV(HEXA_27):
            hexs[hexs_id][0] = conn[i+0];
          hexs[hexs_id][1] = conn[i+1];
          hexs[hexs_id][2] = conn[i+2];
          hexs[hexs_id][3] = conn[i+3];
          hexs[hexs_id][4] = conn[i+4];
          hexs[hexs_id][5] = conn[i+5];
          hexs[hexs_id][6] = conn[i+6];
          hexs[hexs_id][7] = conn[i+7];
              
          unstructured_elements(facehash, &hexs[hexs_id][0], 1 , CGNS_ENUMV(HEXA_8)); 
          hexs_id++;
          break;
          /* ignore these */
          case CGNS_ENUMV(NODE):
          case CGNS_ENUMV(BAR_2):
          case CGNS_ENUMV(BAR_3):
            break;
          /* invalid */
          default:
            sprintf(errmsg,
                    "element type %s not allowed",
                    cg_ElementTypeName(et));
            error_exit( errmsg);
            break;
          }
          int nn= 0;
          cg_npe (et, &nn);
          i += nn;
        }
       
        free (conn);
      }else{
        switch (etype) {
        case CGNS_ENUMV(TRI_3):
        case CGNS_ENUMV(TRI_6):
                
          break;
        case CGNS_ENUMV(QUAD_4):
        case CGNS_ENUMV(QUAD_8):
        case CGNS_ENUMV(QUAD_9):
             
          break;
        case CGNS_ENUMV(TETRA_4):
        case CGNS_ENUMV(TETRA_10):
          if(cg_elements_read (index_file, index_base, index_zone,index_sect, &tets[tets_id][0], NULL))error_exit ("cg_elements_read");
        unstructured_elements(facehash, &tets[tets_id][0], num_elem, CGNS_ENUMV(TETRA_4));
        tets_id += num_elem;
            
             
        break;
        case CGNS_ENUMV(PYRA_5):
        case CGNS_ENUMV(PYRA_13):
        case CGNS_ENUMV(PYRA_14):
          if(cg_elements_read (index_file, index_base, index_zone,index_sect, &pyramids[pyramids_id][0], NULL))error_exit ("cg_elements_read");
        unstructured_elements(facehash, &pyramids[pyramids_id][0], num_elem , CGNS_ENUMV(PYRA_5)); 
        pyramids_id += num_elem;
        break;
        case CGNS_ENUMV(PENTA_6):
        case CGNS_ENUMV(PENTA_15):
        case CGNS_ENUMV(PENTA_18):
          if(cg_elements_read (index_file, index_base, index_zone,index_sect, &prisms[prisms_id][0], NULL))error_exit ("cg_elements_read");
        unstructured_elements(facehash, &prisms[prisms_id][0], num_elem , CGNS_ENUMV(PENTA_6)); 
        prisms_id += num_elem;
        break;
        case CGNS_ENUMV(HEXA_8):
        case CGNS_ENUMV(HEXA_20):
        case CGNS_ENUMV(HEXA_27): 
          if(cg_elements_read (index_file, index_base, index_zone,index_sect, &hexs[hexs_id][0], NULL))error_exit ("cg_elements_read");
        unstructured_elements(facehash, &hexs[hexs_id][0], num_elem , CGNS_ENUMV(HEXA_8)); 
        hexs_id += num_elem;
        break;
        /* ignore these */
        case CGNS_ENUMV(NODE):
        case CGNS_ENUMV(BAR_2):
        case CGNS_ENUMV(BAR_3):
          break;
        /* invalid */
        default:
          sprintf(errmsg,
                  "element type %s not allowed ",
                  cg_ElementTypeName(etype));
          error_exit( errmsg);
          break;
        }
              
      }
    }

            

   
  //read in trias and quads, then hash and sort faces, get boundary information
  
  cout << "start  unstructured_faces " << endl;
  unstructured_faces (facehash, index_file, index_base, index_zone,
                      qfaces, tfaces,
                      surf_ids);

  cout << "end  unstructured_faces " << endl; 
    
   
  cout <<"finish reading cgns file" << endl;
        
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

struct edge {
  pair<int,int> nodes ;
  bool triEdge ;
  bool cvEdge ;
  edge(int n1, int n2, bool te,bool cve) {
    nodes.first=min(n1,n2) ;
    nodes.second=max(n1,n2) ;
    triEdge = te ;
    cvEdge = cve ;
  }
  edge() {} 
} ;

inline bool quadCompare(const quadFace &f1, const quadFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] == f2.nodes[2] && f1.nodes[3] < f2.nodes[3])) ;
}

inline bool triaCompare(const triaFace &f1, const triaFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2])) ;
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

inline bool edgeEqual(const edge &e1, const edge &e2) {
  return ((e1.nodes.first == e2.nodes.first) &&
          (e1.nodes.second == e2.nodes.second)) ;

}

inline bool edgeCompare(const edge &e1, const edge &e2) {
  return ((e1.nodes.first < e2.nodes.first) ||
          (e1.nodes.first == e2.nodes.first && 
           e1.nodes.second < e2.nodes.second)) ;
}

inline bool pairCompare(const pair<int,int> &p1, const pair<int,int> &p2) {
  return ((p1.first < p2.first) ||
          (p1.first == p2.first && 
           p1.second < p2.second)) ;
}


/*
  Find single faces in  alist that do not have the same nodes as its previous neighbor or next neighbor,
  remove them from alist, and place them in the returned vector.
  returned: the single faces in alist
  input: before: alist contains both paired faces and single faces
  after:  alist contains only paired faces
*/
template <class T, class Cmp>
std::vector<T> get_single_face(std::vector<T> &alist, Cmp cmp) {
  vector<T> singles;
  int num_face = alist.size();
  if(num_face == 0) return singles;
  vector<T> pairs;
  int  current = 0 ;
  while(current < num_face){
    int next = current;
    next++;
    if(next == num_face){
      singles.push_back(alist[current]);
      break;
    }else if(cmp(alist[current], alist[next])){
      pairs.push_back(alist[current]);
      pairs.push_back(alist[next]);
      next++;
    }else{
      singles.push_back(alist[current]);
    }
    current = next;    
  }
  alist.swap(pairs);
  return singles;  
}

/*
  Split each quadFace in quad into 2 pairs of triaFace, push them to the end of tria.
  sort the nodes of each newly generated triaFaces in tria so that nodes[0]< nodes[1]< nodes[2]
*/
void split_quad(const vector<quadFace>& quad, vector<triaFace>& tria){
  if(quad.size()==0)return;
  int begin = tria.size();
  for(unsigned int i = 0; i < quad.size(); i++){
    triaFace face1, face2, face3, face4;
    face1.nodes[0] = quad[i].nodes[0];
    face1.nodes[1] = quad[i].nodes[1];
    face1.nodes[2] = quad[i].nodes[2];
    face1.cell = quad[i].cell;
    face1.left = quad[i].left;
    tria.push_back(face1);
    
    face2.nodes[0] = quad[i].nodes[0];
    face2.nodes[1] = quad[i].nodes[2];
    face2.nodes[2] = quad[i].nodes[3];
    face2.cell = quad[i].cell;
    face2.left = quad[i].left;
    tria.push_back(face2);

    face3.nodes[0] = quad[i].nodes[0];
    face3.nodes[1] = quad[i].nodes[1];
    face3.nodes[2] = quad[i].nodes[3];
    face3.cell = quad[i].cell;
    face3.left = quad[i].left;
    tria.push_back(face3);
    
    face4.nodes[0] = quad[i].nodes[1];
    face4.nodes[1] = quad[i].nodes[2];
    face4.nodes[2] = quad[i].nodes[3];
    face4.cell = quad[i].cell;
    face4.left = quad[i].left;
    tria.push_back(face4);
  }
  int end = tria.size();
  for(int i = begin; i < end; i++){
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
}
  
void extract_trifaces(vector<Array<int,5> > &triangles,
                      vector<Array<int,5> > &unbound,
                      int cellid,
                      const vector<Array<int,4> > &tfaces,
                      const vector<Array<int,4> > &tets, 
                      const vector<Array<int,5> > &pyramids,
                      const vector<Array<int,6> > &prisms) {
  int num_tria_faces =
    tfaces.size() + tets.size()*4 + pyramids.size()*4 + prisms.size()*2 ;

  vector<triaFace> tria(num_tria_faces) ;
  int tf = 0 ;
  for(size_t i=0;i<tfaces.size();++i) {
    for(int n=0;n<3;++n)
      tria[tf].nodes[n] = tfaces[i][n] ;
    tria[tf].left = true ;
    tria[tf++].cell = tfaces[i][3] ;
  }
  // Create faces generated by tetrahedra
  for(size_t i=0;i<tets.size();++i) {
    tria[tf].nodes[0] = tets[i][0] ;
    tria[tf].nodes[1] = tets[i][1] ;
    tria[tf].nodes[2] = tets[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = tets[i][1] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;


    tria[tf].nodes[0] = tets[i][3] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][0] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = tets[i][0] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    cellid++ ;
  }

  // create faces generated by pyramids
  for(size_t i=0;i<pyramids.size();++i) {
    tria[tf].nodes[0] = pyramids[i][4] ;
    tria[tf].nodes[1] = pyramids[i][1] ;
    tria[tf].nodes[2] = pyramids[i][2] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][4] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][3] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][0] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][0] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    cellid++ ;
  }

  // create faces generated by prisms
  for(size_t i=0;i<prisms.size();++i) {
    tria[tf].nodes[0] = prisms[i][3] ;
    tria[tf].nodes[1] = prisms[i][4] ;
    tria[tf].nodes[2] = prisms[i][5] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = prisms[i][0] ;
    tria[tf].nodes[1] = prisms[i][2] ;
    tria[tf].nodes[2] = prisms[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    cellid++ ;
  }

  if(tf != num_tria_faces) {
    cerr << "internal consistency error on triangle faces" << endl ;
    Loci::Abort() ;
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
  int tsz = tria.size() ;
  int mtsz ;
  MPI_Allreduce(&tsz,&mtsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  if(mtsz != 0) {
    Loci::parSampleSort(tria,triaCompare,MPI_COMM_WORLD) ;
  }
  vector<triaFace> single_tria = get_single_face(tria, triaEqual);
  unbound = vector<Array<int,5> >(single_tria.size()) ;
  for(size_t i=0;i<single_tria.size();++i) {
    unbound[i][0] = tria[i].nodes[0] ;
    unbound[i][1] = tria[i].nodes[2] ;
    unbound[i][2] = tria[i].nodes[3] ;
    unbound[i][3] = tria[i].cell ;
    unbound[i][4] = tria[i].left ;
  }
  
  int ntria = tria.size()/2 ;
  vector<Array<int,5> > triscratch(ntria) ;
  for(int i=0;i<ntria;++i) {
    for(int j=0;j<3;++j) {
      triscratch[i][j] = tria[i*2].nodes[j] ;
      FATAL(tria[i*2].nodes[j] != tria[i*2+1].nodes[j]) ;
    }
    int c1 = tria[i*2].cell ;
    int c2 = tria[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      Loci::Abort() ;
    }
    if(c1 < 0) {
      triscratch[i][3] = c2 ;
      triscratch[i][4] = c1 ;
      if(tria[i*2+1].left)
        for(int j=0;j<3;++j)
          triscratch[i][j] = tria[i*2+1].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          triscratch[i][j] = tria[i*2+1].nodes[2-j] ;
    } else if(c2 < 0) {
      triscratch[i][3] = c1 ;
      triscratch[i][4] = c2 ;
      if(tria[i*2].left)
        for(int j=0;j<3;++j)
          triscratch[i][j] = tria[i*2].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          triscratch[i][j] = tria[i*2].nodes[2-j] ;
    } else {
      if(!tria[i*2].left) 
        std::swap(c1,c2) ;
      triscratch[i][3] = c1 ;
      triscratch[i][4] = c2 ;
      if((tria[i*2].left && tria[i*2+1].left) ||
         (!tria[i*2].left && !tria[i*2+1].left)) {
        cerr << "consistency error" << endl ;
      }
    }
  }
  triangles.swap(triscratch) ;
}  

void convert2face(store<vector3d<double> > &pos,
                  vector<Array<int,5> > &qfaces, vector<Array<int,4> > &tfaces,
                  vector<Array<int,4> > &tets, vector<Array<int,5> > &pyramids,
                  vector<Array<int,6> > &prisms, vector<Array<int,8> > &hexs,
                  multiMap &face2node,Map &cl, Map &cr) {
  int maxid = 0 ;
  entitySet posDom = pos.domain() ;
  if(posDom != EMPTY)
    maxid = posDom.Max()+1 ;
  int cellid ;

  MPI_Allreduce(&maxid,&cellid,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  int cellbase = cellid ;
  int ncells = tets.size()+pyramids.size()+prisms.size()+hexs.size() ;
  vector<int> cellsizes(MPI_processes) ;
  MPI_Allgather(&ncells,1,MPI_INT,&cellsizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    cellid += cellsizes[i] ;

  int num_quad_faces =
    qfaces.size() + pyramids.size()+ prisms.size()*3 + hexs.size()*6 ;

  int num_tria_faces =
    tfaces.size() + tets.size()*4 + pyramids.size()*4 + prisms.size()*2 ;

  vector<triaFace> tria(num_tria_faces) ;
  vector<quadFace> quad(num_quad_faces) ;
  int tf = 0 ;
  int qf = 0 ;
  for(size_t i=0;i<tfaces.size();++i) {
    for(int n=0;n<3;++n)
      tria[tf].nodes[n] = tfaces[i][n] ;
    tria[tf].left = true ;
    tria[tf++].cell = tfaces[i][3] ;
  }
  for(size_t i=0;i<qfaces.size();++i) {
    for(int n=0;n<4;++n)
      quad[qf].nodes[n] = qfaces[i][n] ;
    quad[qf].left = true ;
    quad[qf++].cell = qfaces[i][4] ;
  }

  // Create faces generated by tetrahedra
  for(size_t i=0;i<tets.size();++i) {
    tria[tf].nodes[0] = tets[i][0] ;
    tria[tf].nodes[1] = tets[i][1] ;
    tria[tf].nodes[2] = tets[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = tets[i][1] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;


    tria[tf].nodes[0] = tets[i][3] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][0] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = tets[i][0] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    cellid++ ;
  }

  // create faces generated by pyramids
  for(size_t i=0;i<pyramids.size();++i) {
    tria[tf].nodes[0] = pyramids[i][4] ;
    tria[tf].nodes[1] = pyramids[i][1] ;
    tria[tf].nodes[2] = pyramids[i][2] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][4] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][3] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][0] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][0] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    quad[qf].nodes[0] = pyramids[i][0] ;
    quad[qf].nodes[1] = pyramids[i][1] ;
    quad[qf].nodes[2] = pyramids[i][4] ;
    quad[qf].nodes[3] = pyramids[i][3] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;
    cellid++ ;
  }

  // create faces generated by prisms
  for(size_t i=0;i<prisms.size();++i) {
    tria[tf].nodes[0] = prisms[i][3] ;
    tria[tf].nodes[1] = prisms[i][4] ;
    tria[tf].nodes[2] = prisms[i][5] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = prisms[i][0] ;
    tria[tf].nodes[1] = prisms[i][2] ;
    tria[tf].nodes[2] = prisms[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    quad[qf].nodes[0] = prisms[i][0] ;
    quad[qf].nodes[1] = prisms[i][1] ;
    quad[qf].nodes[2] = prisms[i][4] ;
    quad[qf].nodes[3] = prisms[i][3] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = prisms[i][1] ;
    quad[qf].nodes[1] = prisms[i][2] ;
    quad[qf].nodes[2] = prisms[i][5] ;
    quad[qf].nodes[3] = prisms[i][4] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = prisms[i][3] ;
    quad[qf].nodes[1] = prisms[i][5] ;
    quad[qf].nodes[2] = prisms[i][2] ;
    quad[qf].nodes[3] = prisms[i][0] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    cellid++ ;
  }

  // create faces generated by hexahedra
  for(size_t i=0;i<hexs.size();++i) {
    quad[qf].nodes[0] = hexs[i][0] ;
    quad[qf].nodes[1] = hexs[i][1] ;
    quad[qf].nodes[2] = hexs[i][5] ;
    quad[qf].nodes[3] = hexs[i][4] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][1] ;
    quad[qf].nodes[1] = hexs[i][2] ;
    quad[qf].nodes[2] = hexs[i][6] ;
    quad[qf].nodes[3] = hexs[i][5] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][2] ;
    quad[qf].nodes[1] = hexs[i][3] ;
    quad[qf].nodes[2] = hexs[i][7] ;
    quad[qf].nodes[3] = hexs[i][6] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][4] ;
    quad[qf].nodes[1] = hexs[i][7] ;
    quad[qf].nodes[2] = hexs[i][3] ;
    quad[qf].nodes[3] = hexs[i][0] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][4] ;
    quad[qf].nodes[1] = hexs[i][5] ;
    quad[qf].nodes[2] = hexs[i][6] ;
    quad[qf].nodes[3] = hexs[i][7] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][3] ;
    quad[qf].nodes[1] = hexs[i][2] ;
    quad[qf].nodes[2] = hexs[i][1] ;
    quad[qf].nodes[3] = hexs[i][0] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    cellid++ ;
  }

  if(qf != num_quad_faces) {
    cerr << "internal consistency error on quad faces" << endl ;
    Loci::Abort() ;
  }
  if(tf != num_tria_faces) {
    cerr << "internal consistency error on triangle faces" << endl ;
    Loci::Abort() ;
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

  
  //sort quad
  int qsz = quad.size() ;
  int mqsz ;
  MPI_Allreduce(&qsz,&mqsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  if(mqsz != 0) {
    Loci::parSampleSort(quad,quadCompare,MPI_COMM_WORLD) ;
  }
  //check if there is single faces in quad, if yes, remove them from quad, and split the single face into tria
  vector<quadFace> single_quad = get_single_face(quad, quadEqual);
  int num_single_quad = single_quad.size();
  int total_num_single_quad = 0;
  MPI_Allreduce(&num_single_quad,&total_num_single_quad,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
 
  if(num_single_quad != 0){
    split_quad(single_quad, tria);
  }

  //sort tria
  int tsz = tria.size() ;
  int mtsz ;
  MPI_Allreduce(&tsz,&mtsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  if(mtsz != 0) {
    Loci::parSampleSort(tria,triaCompare,MPI_COMM_WORLD) ;
  }
  
  //if there is single quad, the split process will generate extra single faces in tria
  // get_single_face will remove them
  if(total_num_single_quad != 0){
    vector<triaFace> single_tria = get_single_face(tria, triaEqual);
    int num_single_tria = single_tria.size();
    int total_num_single_tria = 0;
    MPI_Allreduce(&num_single_tria,&total_num_single_tria,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    if(total_num_single_tria != (2*total_num_single_quad)){
      cerr << "single triangle faces remain! inconsistent! "<< endl ;  
      Loci::Abort() ;
    }
  }
  
  if((tria.size() & 1) == 1) {
    cerr << "non-even number of triangle faces! inconsistent!" << endl ;
    Loci::Abort() ;
  }
  if((quad.size() & 1 ) == 1) {
    cerr << "non-even number of quad faces! inconsistent!" << endl ;
    Loci::Abort() ;
  }
  
  int ntria = tria.size()/2 ;
  int nquad = quad.size()/2 ;
  
  int nfaces = ntria+nquad ;

  ncells = 0 ;
  for(int i=0;i<MPI_processes;++i)
    ncells += cellsizes[i] ;
  int facebase = cellbase + ncells ;
  vector<int> facesizes(MPI_processes) ;
  MPI_Allgather(&nfaces,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    facebase += facesizes[i] ;
  entitySet faces = interval(facebase,facebase+nfaces-1) ;
  store<int> count ;
  cl.allocate(faces) ;
  cr.allocate(faces) ;
  count.allocate(faces) ;
  int fc = facebase ;
  for(int i=0;i<ntria;i++)
    count[fc++] = 3 ;
  for(int i=0;i<nquad;i++)
    count[fc++] = 4 ;
  face2node.allocate(count) ;
  fc = facebase ;
  for(int i=0;i<ntria;++i) {
    for(int j=0;j<3;++j) {
      face2node[fc][j] = tria[i*2].nodes[j] ;
      FATAL(tria[i*2].nodes[j] != tria[i*2+1].nodes[j]) ;
    }
    int c1 = tria[i*2].cell ;
    int c2 = tria[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      Loci::Abort() ;
    }
    if(c1 < 0) {
      cl[fc] = c2 ;
      cr[fc] = c1 ;
      if(tria[i*2+1].left)
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2+1].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2+1].nodes[2-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(tria[i*2].left)
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2].nodes[2-j] ;
    } else {
      if(!tria[i*2].left) 
        std::swap(c1,c2) ;
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if((tria[i*2].left && tria[i*2+1].left) ||
         (!tria[i*2].left && !tria[i*2+1].left)) {
        cerr << "consistency error" << endl ;
      }
    }
    fc++ ;
  }
  for(int i=0;i<nquad;++i) {
    for(int j=0;j<4;++j) {
      face2node[fc][j] = quad[i*2].nodes[j] ;
      FATAL(quad[i*2].nodes[j] != quad[i*2+1].nodes[j]) ;
    }
    int c1 = quad[i*2].cell ;
    int c2 = quad[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      Loci::Abort() ;
    }
    
    if(c1 < 0) {
      cl[fc] = c2 ;
      cr[fc] = c1 ;
      if(quad[i*2+1].left)
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2+1].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2+1].nodes[3-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(quad[i*2].left)
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2].nodes[3-j] ;
    } else {
      if(!quad[i*2].left) 
        std::swap(c1,c2) ;
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if((quad[i*2].left && quad[i*2+1].left) ||
         (!quad[i*2].left && !quad[i*2+1].left)) {
        cerr << "consistency error" << endl ;
      }
    }
    fc++ ;
  }

}

inline void find_edge(int &edge_id, int &edge_orient, 
                      const vector<edge> &edgeData,
                      int edgelist[], int nedges, int id1, int id2) {
  edge_orient = 1 ;
  for(int i=0;i<nedges;++i) {
    if(edgeData[edgelist[i]].nodes.first == id1 && 
       edgeData[edgelist[i]].nodes.second == id2) {
      edge_id = edgelist[i] ;
      return ;
    }
    if(edgeData[edgelist[i]].nodes.first == id2 && 
       edgeData[edgelist[i]].nodes.second == id1) {
      edge_id = edgelist[i] ;
      edge_orient=-1 ;
      return ;
    }
  }
  edge_id = -1 ;
  edge_orient = 0 ;
  cerr << "unable to find edge (" << id1 << "," << id2 << ")" << endl ;
  cerr << "searching: " ;
  for(int i=0;i<nedges;++i) {
    cerr << "(" 
         << edgeData[edgelist[i]].nodes.first << ","
         << edgeData[edgelist[i]].nodes.second << ") " ;
  }
  cerr << endl ;
}

// Create edge structures by looping over elements and extracting
// edges
void computeEdgeMap(vector<edge> &edgeData,
                    const vector<Array<int,4> > &tets,
                    const vector<Array<int,5> > &pyramids,
                    const vector<Array<int,6> > &prisms, 
                    const vector<Array<int,8> > &hexs) {
  vector<edge> edgeInfo(tets.size()*6+pyramids.size()*8+
                        prisms.size()*9+hexs.size()*12) ;
  // Create edges generated by tetrahedra
  int cnt= 0 ;
  for(size_t i=0;i<tets.size();++i) {
    edgeInfo[cnt++] = edge(tets[i][0],tets[i][1],true,false) ;
    edgeInfo[cnt++] = edge(tets[i][0],tets[i][2],true,false) ;
    edgeInfo[cnt++] = edge(tets[i][0],tets[i][3],true,false) ;
    edgeInfo[cnt++] = edge(tets[i][1],tets[i][2],true,false) ;
    edgeInfo[cnt++] = edge(tets[i][1],tets[i][3],true,false) ;
    edgeInfo[cnt++] = edge(tets[i][2],tets[i][3],true,false) ;
  }

  // create faces generated by pyramids
  for(size_t i=0;i<pyramids.size();++i) {
    // base
    edgeInfo[cnt++] = edge(pyramids[i][0],pyramids[i][1],true,false) ;
    edgeInfo[cnt++] = edge(pyramids[i][1],pyramids[i][4],true,false) ;
    edgeInfo[cnt++] = edge(pyramids[i][4],pyramids[i][3],true,false) ;
    edgeInfo[cnt++] = edge(pyramids[i][3],pyramids[i][0],true,false) ;
    // tip
    edgeInfo[cnt++] = edge(pyramids[i][0],pyramids[i][2],true,false) ;
    edgeInfo[cnt++] = edge(pyramids[i][1],pyramids[i][2],true,false) ;
    edgeInfo[cnt++] = edge(pyramids[i][3],pyramids[i][2],true,false) ;
    edgeInfo[cnt++] = edge(pyramids[i][4],pyramids[i][2],true,false) ;
  }

  // create faces generated by prisms
  for(size_t i=0;i<prisms.size();++i) {
    edgeInfo[cnt++] = edge(prisms[i][0],prisms[i][1],true,false) ;
    edgeInfo[cnt++] = edge(prisms[i][1],prisms[i][2],true,false) ;
    edgeInfo[cnt++] = edge(prisms[i][2],prisms[i][0],true,false) ;
    edgeInfo[cnt++] = edge(prisms[i][3],prisms[i][4],true,false) ;
    edgeInfo[cnt++] = edge(prisms[i][4],prisms[i][5],true,false) ;
    edgeInfo[cnt++] = edge(prisms[i][5],prisms[i][3],true,false) ;
     
    edgeInfo[cnt++] = edge(prisms[i][0],prisms[i][3],false,true) ;
    edgeInfo[cnt++] = edge(prisms[i][1],prisms[i][4],false,true) ;
    edgeInfo[cnt++] = edge(prisms[i][2],prisms[i][5],false,true) ;
  }

  // create faces generated by hexahedra
  for(size_t i=0;i<hexs.size();++i) {
    // bottom edges
    edgeInfo[cnt++] = edge(hexs[i][0],hexs[i][1],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][1],hexs[i][2],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][2],hexs[i][3],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][3],hexs[i][0],false,false) ;

    // top edges
    edgeInfo[cnt++] = edge(hexs[i][4],hexs[i][5],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][5],hexs[i][6],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][6],hexs[i][7],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][7],hexs[i][4],false,false) ;

    // edges connecting bottom totop
    edgeInfo[cnt++] = edge(hexs[i][0],hexs[i][4],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][1],hexs[i][5],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][2],hexs[i][6],false,false) ;
    edgeInfo[cnt++] = edge(hexs[i][3],hexs[i][7],false,false) ;
  }
  
  // remove duplicate data
  {
    Loci::parSampleSort(edgeInfo,edgeCompare,MPI_COMM_WORLD) ;
    int esz = edgeInfo.size() ;
    int ecount = 0 ;
    int split_count = 0;
    int prisme_count = 0 ;
    int conflict_count = 0;
    vector<edge> edgelist ;
    for(int i=0;i<esz;++i) {
      bool split_edge = false ;
      bool prismedge = false ;
      int j=i ;
      for(j=i;j<esz && edgeEqual(edgeInfo[i],edgeInfo[j]);++j) {
        split_edge = split_edge||edgeInfo[j].triEdge ;
        prismedge = prismedge||edgeInfo[j].cvEdge ;
      }
      edgeInfo[i].triEdge = split_edge ;
      edgeInfo[i].cvEdge = prismedge ;
      edgelist.push_back(edgeInfo[i]) ;
      i=j-1 ;
      ecount++ ;
      if(split_edge)
        split_count++ ;
      if(prismedge)
        prisme_count++ ;
      if(prismedge && split_edge)
        conflict_count++ ;
    }
    edgeData = edgelist ;
    cout << "edge count = " << ecount << ", split edges=" << split_count
         << ", prism edges = " << prisme_count 
         <<", conflicts=" << conflict_count
         << endl ;
  }
}

// compute edge2node mapping by inverting edge data references
void getEdge2Node(vector<int> &n2e, vector<int> &n2e_off,
                  const vector<edge> &edgeData) {
  // edgeData is now the edge map
  int ecount = edgeData.size() ;

  // compute edge2node map
  vector<pair<int,int> > node2edge(ecount*2) ;
  for(int i=0;i<ecount;++i) {
    node2edge[i*2].second = i ;
    node2edge[i*2].first = edgeData[i].nodes.first ;
    node2edge[i*2+1].second = i ;
    node2edge[i*2+1].first = edgeData[i].nodes.second ;
  }

  Loci::parSampleSort(node2edge,pairCompare,MPI_COMM_WORLD) ;

  if(node2edge[0].first != 1) {
    cerr << "node numbering should start with 1 in node2edge" << endl ;
  }
  int num_nodes = node2edge[(ecount-1)*2+1].first ;
  cout << "num_nodes = " << num_nodes << endl ;
  vector<int> ecounts(num_nodes,0) ;
  int edg = 0 ;
  for(int i=1;i<num_nodes+1;++i) {
    int edgo = edg ;
    while(edg < ecount*2 && node2edge[edg].first == i) {
      edg++ ;
    }
    ecounts[i-1] = edg-edgo ;
  }
  vector<int> offsets(num_nodes+1) ;
  offsets[0] = 0 ;
  for(int i=0;i<num_nodes;++i) {
    offsets[i+1] = offsets[i]+ecounts[i] ;
  }
  vector<int> n2e_tmp(node2edge.size()) ;
  for(size_t i=0;i<node2edge.size();++i)
    n2e_tmp[i] = node2edge[i].second ;

  n2e.swap(n2e_tmp) ;
  n2e_off.swap(offsets) ;
}

inline bool Array4Compare0(const Array<int,4> &p1, const Array<int,4> &p2) {
  return (p1[0] < p2[0]) ;
}
inline bool SpecialFaceCompare(const Array<int,4> &p1, const Array<int,4> &p2) {
  return (p1[0] < p2[0]) ||
    ((p1[0] == p2[0]) && (p1[3] < p2[3])) ;
}

// Currently only works in serial for all tetrahedral meshes

void convert2cellVertexface(store<vector3d<double> > &pos,
                            vector<Array<int,5> > &qfaces, vector<Array<int,4> > &tfaces,
                            vector<Array<int,4> > &tets, vector<Array<int,5> > &pyramids,
                            vector<Array<int,6> > &prisms, vector<Array<int,8> > &hexs,
                            multiMap &face2node,Map &cl, Map &cr) {

  // Compute edge data structures
  vector<edge> edgeData ;
  computeEdgeMap(edgeData,tets,pyramids,prisms,hexs) ;
  vector<int> node2edge ;
  vector<int> node2edge_off ;
  getEdge2Node(node2edge,node2edge_off,edgeData) ;
  // to access node2edge for node i access node2edge from
  // node2edge_off[i] to node2edge_off[i+1]
  // number of edges
  int ecount = edgeData.size() ;


  //-------------------------------------------------------------------------
  // Identify new control volumes
  //-------------------------------------------------------------------------

  vector<int> tet_cnts(pos.domain().size(),0) ;
  
  // count tets that neighbor node
  for(size_t i=0;i<tets.size();++i) 
    for(int j=0;j<4;++j)  
      tet_cnts[tets[i][j]-1]++ ;
  // Find all vertexes that convert to control volumes
  entitySet vertexCVs, notVertexCVs ;
  for(size_t i=0;i<tet_cnts.size();++i) {
    if(tet_cnts[i] > 2)
      vertexCVs += int(i) ;
    else
      notVertexCVs += int(i) ;
  }

  cout << "notVertexCVs = " << notVertexCVs << endl ;
 
  //-------------------------------------------------------------------------
  // Break corners off of tets
  //-------------------------------------------------------------------------
 
  // Now we collect the triangular facets that go to each node
  // the data will include the node number, the contributing cell number
  // and the numbering face nodes that form the triangular face
  vector<Array<int,5> > nodefacets ;
  
  // now loop over tets and gather and break them down
  for(size_t i=0;i<tets.size();++i) {
    // edges of tets, e0 = t0-t1, e1 = t0-t2, e2=t0-t3
    //                e3 = t1-t2, e4 = t1-t3, e5=t2-t3
    // corner  c0 = e0,e1,e2
    //         c1 = e0,e3,e4
    //         c2 = e1,e3,e5
    //         c3 = e2,e4,e5
    // edge_orientation 1 or -1
    int el[6] ; // tet edges 
    int eo[6] ; // tet edge orientation ;
    int no = node2edge_off[tets[i][0]-1] ;
    int nsearch = node2edge_off[tets[i][0]] - no ;
    find_edge(el[0],eo[0],edgeData,
              &node2edge[no],nsearch,tets[i][0],tets[i][1]) ;
    find_edge(el[1],eo[1],edgeData,
              &node2edge[no],nsearch,tets[i][0],tets[i][2]) ;
    find_edge(el[2],eo[2],edgeData,
              &node2edge[no],nsearch,tets[i][0],tets[i][3]) ;
    no = node2edge_off[tets[i][1]-1] ;
    nsearch = node2edge_off[tets[i][1]] - no ;
    find_edge(el[3],eo[3],edgeData,
              &node2edge[no],nsearch,tets[i][1],tets[i][2]) ;
    find_edge(el[4],eo[4],edgeData,
              &node2edge[no],nsearch,tets[i][1],tets[i][3]) ;
    no = node2edge_off[tets[i][2]-1] ;
    nsearch = node2edge_off[tets[i][2]] - no ;
    find_edge(el[5],eo[5],edgeData,
              &node2edge[no],nsearch,tets[i][2],tets[i][3]) ;
    // now output corner 0
    if(vertexCVs.inSet(tets[i][0]-1)) {
      Array<int,5> facet ;
      facet[0] = tets[i][0]-1 ;
      facet[1] = i ; // number tets first
      facet[2] = el[0]*2+((eo[0]>0)?0:1) ;
      facet[3] = el[1]*2+((eo[1]>0)?0:1) ;
      facet[4] = el[2]*2+((eo[2]>0)?0:1) ;
      nodefacets.push_back(facet) ;
    }
    // now output corner 1
    if(vertexCVs.inSet(tets[i][1]-1) ) {
      Array<int,5> facet ;
      facet[0] = tets[i][1]-1 ;
      facet[1] = i ; // number tets first
      facet[3] = el[0]*2+((eo[0]>0)?1:0) ;
      facet[2] = el[3]*2+((eo[3]>0)?0:1) ;
      facet[4] = el[4]*2+((eo[4]>0)?0:1) ;
      nodefacets.push_back(facet) ;
    }
    // now output corner 2
    if(vertexCVs.inSet(tets[i][2]-1) ) {
      Array<int,5> facet ;
      facet[0] = tets[i][2]-1 ;
      facet[1] = i ; // number tets first
      facet[2] = el[1]*2+((eo[1]>0)?1:0) ;
      facet[3] = el[3]*2+((eo[3]>0)?1:0) ;
      facet[4] = el[5]*2+((eo[5]>0)?0:1) ;
      nodefacets.push_back(facet) ;
    }
    // now output corner 3
    if(vertexCVs.inSet(tets[i][3]-1) ) {
      Array<int,5> facet ;
      facet[0] = tets[i][3]-1 ;
      facet[1] = i ; // number tets first
      facet[3] = el[2]*2+((eo[2]>0)?1:0) ;
      facet[2] = el[4]*2+((eo[4]>0)?1:0) ;
      facet[4] = el[5]*2+((eo[5]>0)?1:0) ;
      nodefacets.push_back(facet) ;
    }
  }
  
  vector<Loci::Array<int,4> > bfaceinfo ;
  // process boundary faces
  for(size_t i=0;i<tfaces.size();++i) {
    int n1 = tfaces[i][0] ;
    int n2 = tfaces[i][1] ;
    int n3 = tfaces[i][2] ;
    int el[3] ; // tri edges 
    int eo[3] ; // tri edge orientation ;
    int no = node2edge_off[n1-1] ;
    int nsearch = node2edge_off[n1] - no ;
    find_edge(el[0],eo[0],edgeData,
              &node2edge[no],nsearch,n1,n2) ;
    find_edge(el[2],eo[2],edgeData,
              &node2edge[no],nsearch,n3,n1) ;
    no = node2edge_off[n2-1] ;
    nsearch = node2edge_off[n2]-no ;
    find_edge(el[1],eo[1],edgeData,
              &node2edge[no],nsearch,n2,n3) ;
    if(vertexCVs.inSet(n1-1)) {
      Loci::Array<int,4> tmpface ;
      tmpface[0] = n1-1 ;
      tmpface[1] = el[2]*2+((eo[2]>0)?1:0) ;
      tmpface[2] = el[0]*2+((eo[0]>0)?0:1) ;
      tmpface[3] = tfaces[i][3] ;
      bfaceinfo.push_back(tmpface) ;
    }
    if(vertexCVs.inSet(n2-1)) {
      Loci::Array<int,4> tmpface ;
      tmpface[0] = n2-1 ;
      tmpface[1] = el[0]*2+((eo[0]>0)?1:0) ;
      tmpface[2] = el[1]*2+((eo[1]>0)?0:1) ;
      tmpface[3] = tfaces[i][3] ;
      bfaceinfo.push_back(tmpface) ;
    }
    if(vertexCVs.inSet(n3-1)) {
      Loci::Array<int,4> tmpface ;
      tmpface[0] = n3-1 ;
      tmpface[1] = el[1]*2+((eo[1]>0)?1:0) ; 
      tmpface[2] = el[2]*2+((eo[2]>0)?0:1) ;
      tmpface[3] = tfaces[i][3] ;
      bfaceinfo.push_back(tmpface) ;
    }
  }
  
  sort(bfaceinfo.begin(),bfaceinfo.end(),SpecialFaceCompare) ;

  
  
  entitySet new_bnodes ;
  vector<int> nodecounts ;
  vector<int> bface_sizes ;
  vector<int> bface_type ;
  size_t bfsz = bfaceinfo.size() ;
  for(size_t i=0;i<bfsz;) {
    int c = 0 ;
    for(size_t j=i;(j<bfsz) &&(bfaceinfo[i][0]==bfaceinfo[j][0]);j++)
      c++ ;
    nodecounts.push_back(c) ;
    bool hasnode=false ;
    for(int j=1;j<c;++j)
      if(bfaceinfo[i][3] != bfaceinfo[i+j][3])
        hasnode = true ;
    if(hasnode) {
      new_bnodes += bfaceinfo[i][0] ;
      // now divide face up  based on boundary id
      for(int j=0;j<c;) {
        int c2 = 0 ;
        for(int k=j;(k<c)&&(bfaceinfo[i+j][3]==bfaceinfo[i+k][3]);++k)
          c2++ ;
        bface_sizes.push_back(c2) ;
        bface_type.push_back(1) ;
        // Now sort this group
        // find the first entry
        for(int k=1;k<c2;++k) {
          if(bfaceinfo[i+j][1] == bfaceinfo[i+j+k][2]) {
            std::swap(bfaceinfo[i+j],bfaceinfo[i+j+k]) ;
            k=0 ;
          }
        }
        // and then sort
        for(int k=1;k<c2;++k) {
          int look = bfaceinfo[i+j+k-1][2] ;
          for(int kk=k;kk<c2;++kk)
            if(bfaceinfo[i+j+kk][1] == look)
              std::swap(bfaceinfo[i+j+k],bfaceinfo[i+j+kk]) ;
        }
        // check sort
        bool checksort = true ;
        for(int k=1;k<c2;++k) {
          if(bfaceinfo[i+j+k-1][2] != bfaceinfo[i+j+k][1]) {
            checksort = false ;
          }
        }
        if(!checksort) {
          cout << "checksort failed on boundary split face" << endl ;
          for(int k=0;k<c2;++k) {
            cout << "("<< bfaceinfo[i+j+k][1] << "," << bfaceinfo[i+j+k][2] << ")" << endl ;
          }
        }
          

        j+=c2 ;
      }
          
    } else {
      //sort this group of bfaces
      for(int j=1;j<c;++j) {
        int look = bfaceinfo[j-1+i][2] ;
        for(int k=j;k<c;++k)
          if(bfaceinfo[k+i][1] == look)
            std::swap(bfaceinfo[j+i],bfaceinfo[k+i]) ;
      }
      // check sorting
      bool checksort = true ;
      for(int j=1;j<c;++j) {
        if(bfaceinfo[j-1+i][2] != bfaceinfo[j+i][1]) {
          checksort = false ;
        }
      }
      if(!checksort) {
        cout << "unsorted face! " << bfaceinfo[i][0] << ":" << i << "," ;
        for(int j=0;j<c;++j)
          cout << "("
               << bfaceinfo[j+i][1]<<","<<bfaceinfo[j+i][2]<<")" ;
        cout << endl ;
      }
      bface_sizes.push_back(c) ;
      bface_type.push_back(0) ;
    }
    i+= c ;
  }
  

  entitySet keep_nodes = notVertexCVs + new_bnodes ;
  int numnewnodes = ecount*2 + keep_nodes.size() ;
  

  //-------------------------------------------------------------------------
  // Compute positions of nodes in new polyhedral mesh
  //-------------------------------------------------------------------------

  // compute an average radius for each node to determine
  // optimal splitting locations for edge
  vector<int> edge_cnts(pos.domain().size(),0) ;
  vector<double> edge_radius(pos.domain().size(),0) ;
  
  for(int i=0;i<ecount;++i) {
    int n1 = edgeData[i].nodes.first-1 ;
    int n2 = edgeData[i].nodes.second-1 ;
    vector3d<double> p1 = pos[n1] ;
    vector3d<double> p2 = pos[n2] ;
    double el = norm(p1-p2)/3.0 ; // individual optimal split edge into 3 parts
    edge_cnts[n1]++ ;
    edge_cnts[n2]++ ;
    edge_radius[n1] += el ;
    edge_radius[n2] += el ;
  }
  for(size_t i=0;i<edge_radius.size();++i)
    edge_radius[i] *= 1./double(edge_cnts[i]) ;

  store<vector3d<double> > newpos ;
  entitySet newdom = interval(0,numnewnodes-1) ;
  newpos.allocate(newdom) ;
  for(int i=0;i<ecount;++i) {
    int n1 = edgeData[i].nodes.first-1 ;
    int n2 = edgeData[i].nodes.second-1 ;
    vector3d<double> p1 = pos[n1] ;
    vector3d<double> p2 = pos[n2] ;
    vector3d<double> dv = p2-p1 ;
    double el = norm(dv) ;
    // bound the nodal radius so it still makes sense for this edge
    double r1 = max(min(edge_radius[n1]/el,0.45),0.1) ;
    double r2 = max(min(edge_radius[n2]/el,0.45),0.1) ;
    newpos[i*2+0] = p1+r1*dv ; //(2.*pos[n1]+pos[n2])/3.0 ;
    newpos[i*2+1] = p2-r2*dv ; //(pos[n1]+2.0*pos[n2])/3.0 ;
  }
  vector<int> nodeids(pos.domain().size(),-1) ;
  int cnt = 0 ;
  FORALL(keep_nodes,ii) {
    newpos[cnt+2*ecount] = pos[ii] ;
    nodeids[ii] = cnt+2*ecount ;
    cnt++ ;
  } ENDFORALL;

  //-------------------------------------------------------------------------
  // Establish cell numbering
  //-------------------------------------------------------------------------
  int norig_cells = tets.size()+pyramids.size()+prisms.size()+hexs.size() ;
  vector<int> vertexIds(pos.domain().size(),-1) ;
  int cellid = numnewnodes ; // id of cells base
  cnt = numnewnodes+norig_cells ; // start numbering vertexCVs after maincells
  FORALL(vertexCVs,ii) {
    vertexIds[ii] = cnt ;
    cnt++ ;
  } ENDFORALL ;
  int ncells = norig_cells+vertexCVs.size() ;

  cout << "numnodes = " << numnewnodes << ", vertexCVs=" << vertexCVs.size()
       << ", ncells=" << ncells << endl;


  //-------------------------------------------------------------------------
  // Establish faces of new mesh
  //-------------------------------------------------------------------------

  // Now we need to make the faces, first we will convert the baseline mesh to faces, and
  // the triangular faces will change to hexagons
  vector<Array<int,5> > triangles ;
  vector<Array<int,5> > tria_unbound ;

  
  extract_trifaces(triangles,tria_unbound,cellid,tfaces,tets,pyramids,prisms) ;
  if(tria_unbound.size() != 0) {
    cerr << "unable to process unbound triangle faces" << endl ;
    exit(-1) ;
  }
  
  cout << triangles.size() << " hexagonal faces generated" << endl ;

  cout << nodefacets.size() << " node faces " << endl ;

  cout << bface_sizes.size() << " boundary node faces " << endl ;
  cout << triangles.size()+nodefacets.size()+bface_sizes.size() << " total faces"<< endl ;
    

  int nfaces = triangles.size()+nodefacets.size()+bface_sizes.size() ;
  
  entitySet fdom = interval(cellid+ncells,cellid+ncells+nfaces-1) ;
  Map ncl,ncr ;
  store<int> count ;
  ncl.allocate(fdom) ;
  ncr.allocate(fdom) ;
  count.allocate(fdom) ;
  int fbase = cellid+ncells ;
  for(size_t i=0;i<triangles.size();++i) {
    count[i+fbase] = 6 ;
    for(int j=0;j<3;++j)
      if(!vertexCVs.inSet(triangles[i][j])) 
        count[i+fbase]++ ;
  }
  fbase += triangles.size() ;
  for(size_t i=0;i<nodefacets.size();++i)
    count[i+fbase] = 3 ;
  fbase += nodefacets.size() ;

  for(size_t i=0;i<bface_sizes.size();++i) {
    count[i+fbase] = bface_sizes[i] ;
    if(bface_type[i] > 0)
      count[i+fbase] = bface_sizes[i]+2 ;
  }

  multiMap nface2node ;
  nface2node.allocate(count) ;
  int fcnt =  cellid+ncells ;
  // first fill in hex faces
  for(size_t i=0;i<triangles.size();++i) {
    int n1 = triangles[i][0]+1 ;
    int n2 = triangles[i][1]+1 ;
    int n3 = triangles[i][2]+1 ;
    int el[3] ; // tet edges 
    int eo[3] ; // tet edge orientation ;
    int no = node2edge_off[n1-1] ;
    int nsearch = node2edge_off[n1] - no ;
    find_edge(el[0],eo[0],edgeData,
              &node2edge[no],nsearch,n1,n2) ;
    find_edge(el[2],eo[2],edgeData,
              &node2edge[no],nsearch,n3,n1) ;
    no = node2edge_off[n2-1] ;
    nsearch = node2edge_off[n2]-no ;
    find_edge(el[1],eo[1],edgeData,
              &node2edge[no],nsearch,n2,n3) ;
    
    int lcnt = 0 ;
    if(!vertexCVs.inSet(triangles[i][0]))
      nface2node[fcnt][lcnt++] = nodeids[triangles[i][0]] ;
    nface2node[fcnt][lcnt++] = el[0]*2+((eo[0]>0)?0:1) ;
    nface2node[fcnt][lcnt++] = el[0]*2+((eo[0]>0)?1:0) ;
    if(!vertexCVs.inSet(triangles[i][1]))
      nface2node[fcnt][lcnt++] = nodeids[triangles[i][1]] ;
    nface2node[fcnt][lcnt++] = el[1]*2+((eo[1]>0)?0:1) ;
    nface2node[fcnt][lcnt++] = el[1]*2+((eo[1]>0)?1:0) ;
    if(!vertexCVs.inSet(triangles[i][2]))
      nface2node[fcnt][lcnt++] = nodeids[triangles[i][2]] ;
    nface2node[fcnt][lcnt++] = el[2]*2+((eo[2]>0)?0:1) ;
    nface2node[fcnt][lcnt++] = el[2]*2+((eo[2]>0)?1:0) ;

    ncl[fcnt] = triangles[i][3] ;
    ncr[fcnt] = triangles[i][4] ;
    fcnt++ ;
  }
  // now fill in node cells from tetrahedra faces
  for(size_t i=0;i<nodefacets.size();++i) {
    nface2node[fcnt][0] = nodefacets[i][2] ;
    nface2node[fcnt][1] = nodefacets[i][3] ;
    nface2node[fcnt][2] = nodefacets[i][4] ;
    ncl[fcnt] = vertexIds[nodefacets[i][0]] ;
    ncr[fcnt] = nodefacets[i][1]+cellid ;
    fcnt++ ;
  }
  
  // Now fill in the boundary faces to nodal CV's
  cnt = 0 ;
  for(size_t i=0;i<bface_sizes.size();++i) {
    int bs = bface_sizes[i] ;
    for(int j=0;j<bs;++j) {
      nface2node[fcnt][j] = bfaceinfo[cnt+j][1] ;
    }
    if(bface_type[i] > 0) {
      nface2node[fcnt][bs] = bfaceinfo[cnt+bs-1][2] ;
      nface2node[fcnt][bs+1] = nodeids[bfaceinfo[cnt][0]] ;
      if(nodeids[bfaceinfo[cnt][0]] < 0)
        cerr << "nodeids not set for boundary node" << endl ;
    }
    ncl[fcnt] = vertexIds[bfaceinfo[cnt][0]] ;
    ncr[fcnt] = bfaceinfo[cnt][3] ;
    fcnt++ ;
    cnt += bs ;
  }
  
  face2node.setRep(nface2node.Rep()) ;
  cl.setRep(ncl.Rep()) ;
  cr.setRep(ncr.Rep()) ;
  pos.setRep(newpos.Rep()) ;
  //  exit(-1) ;
}

int main(int ac, char* av[]) {
  using namespace VOG ;

  bool optimize = true ;
  bool cellVertexTransform = false ;
  string zonepath;
  Loci::Init(&ac,&av) ;
  string Lref = "NOSCALE" ;
  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    
    }else if(ac >= 2 && !strcmp(av[1],"-v")) {
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

 

  if(ac == 2) {
    zonepath = av[1];
  } else {
    cerr << "Usage: cgns2vog <options> <filename>/<basename>/<zonename>" << endl
         <<"If <zonename> is not specified, the first zone of the base is read" << endl
         <<"If <basename> is not specified, the first base in the file is read" << endl 
         << "Where options are listed below and " << endl
         << "flags:" << endl
        
         << "  -o  : disable optimization that reorders nodes and faces" << endl
         << "  -v  : display version" << endl
         << "  -in : input grid is in inches" << endl
         << "  -ft : input grid is in feet" << endl
         << "  -cm : input grid is in centimeters" << endl
         << "  -m  : input grid is in meters" << endl
         << "  -Lref <units> : 1 unit in input grid is <units> long" << endl
         << endl ;
    exit(-1) ;
  }

  if(Lref == "NOSCALE") {
    cerr << "Must set grid units!" << endl
         << "Use options -in, -ft, -cm, -m, or -Lref to set grid units." << endl ;
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
  
 

  int loc = 0;
  loc = zonepath.find('.') ;
  std::string casename = zonepath.substr(0, loc) ;
  string outfile = casename + string(".vog") ;

  store<vector3d<double> > pos ;
  vector<Array<cgsize_t,5> > qfaces_long ;
  vector<Array<cgsize_t,4> > tfaces_long ;
  vector<Array<cgsize_t,4> > tets_long ;
  vector<Array<cgsize_t,5> > pyramids_long ;
  vector<Array<cgsize_t,6> > prisms_long ;
  vector<Array<cgsize_t,8> > hexs_long ;
  vector<pair<int,string> > surf_ids ;
  // if(Loci::MPI_processes > 1)readCGNS_parallel(infile, pos,qfaces_long,tfaces_long,tets_long,pyramids_long,prisms_long,hexs_long, surf_ids) ;
  //else
  readCGNS_serial(zonepath, pos,qfaces_long,tfaces_long,tets_long,pyramids_long,prisms_long,hexs_long, surf_ids) ;
  

  

   

  
  vector<Array<int,5> > qfaces ;
  vector<Array<int,4> > tfaces ;
  vector<Array<int,4> > tets ;
  vector<Array<int,5> > pyramids ;
  vector<Array<int,6> > prisms ;
  vector<Array<int,8> > hexs ;

  //static_cast all data to keep clean interface
  qfaces.resize(qfaces_long.size());
  tfaces.resize(tfaces_long.size());
  tets.resize(tets_long.size());
  pyramids.resize(pyramids_long.size());
  prisms.resize(prisms_long.size());
  hexs.resize(hexs_long.size());
  for(size_t i = 0; i < qfaces.size(); i++){
    for(int j = 0; j < 4; j++){
      qfaces[i][j] = static_cast<int>(qfaces_long[i][j]);
    }
    qfaces[i][4] = static_cast<int>(-qfaces_long[i][4]);
  }
  for(size_t i = 0; i < tfaces.size(); i++){
    for(int j = 0; j < 3; j++){
      tfaces[i][j] = static_cast<int>(tfaces_long[i][j]);
    }
    tfaces[i][3] = static_cast<int>(-tfaces_long[i][3]); 
  }
  for(size_t i = 0; i < tets.size(); i++){
    for(int j = 0; j < 4; j++){
      tets[i][j] = static_cast<int>(tets_long[i][j]);
    }
  }
  //the connectivity of pyramid is modified to ugrid connectivity.
  //So that the convert2face() of ugrid2vod can be used without modification
  for(size_t i = 0; i < pyramids.size(); i++){
    pyramids[i][0] = static_cast<int>(pyramids_long[i][0]);
    pyramids[i][1] = static_cast<int>(pyramids_long[i][3]);
    pyramids[i][2] = static_cast<int>(pyramids_long[i][4]);
    pyramids[i][3] = static_cast<int>(pyramids_long[i][1]); 
    pyramids[i][4] = static_cast<int>(pyramids_long[i][2]);   
    
  }
  for(size_t i = 0; i < prisms.size(); i++){
    for(int j = 0; j < 6; j++){
      prisms[i][j] = static_cast<int>(prisms_long[i][j]);
    }
  }
  for(size_t i = 0; i < hexs.size(); i++){
    for(int j = 0; j < 8; j++){
      hexs[i][j] = static_cast<int>(hexs_long[i][j]);
    }
  }
  vector<Array<cgsize_t,5> >().swap(qfaces_long);
  vector<Array<cgsize_t,4> >().swap(tfaces_long);
  vector<Array<cgsize_t,4> >().swap(tets_long);
  vector<Array<cgsize_t,5> >().swap(pyramids_long);
  vector<Array<cgsize_t,6> >().swap(prisms_long);
  vector<Array<cgsize_t,8> >().swap(hexs_long);

  //for debug purpose
  // writeASCII("cgns2ascii.txt", pos,qfaces,tfaces,tets,pyramids,prisms,hexs,surf_ids ) ;

  if(posScale != 1.0) {
    FORALL(pos.domain(),nd) {
      pos[nd] *= posScale ;
    } ENDFORALL ;
  }
   
  
  multiMap face2node ;
  Map cl,cr ;

  if(cellVertexTransform) {
    cout << " convert2cellVertexface " << endl;
    convert2cellVertexface(pos,qfaces,tfaces,tets,pyramids,prisms,hexs,
                           face2node,cl,cr) ;
		       
  } else {
    cout << " convert2face " << endl;
    convert2face(pos,qfaces,tfaces,tets,pyramids,prisms,hexs,
                 face2node,cl,cr) ;
  }

  
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
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;

  Loci::Finalize() ;
  
}
#else
int main(int ac, char *av[]) {
  fprintf(stderr,"Loci not compiled with CGNS support enabled! This utility cannot work!\n") ;
  return -1 ;
}
#endif
