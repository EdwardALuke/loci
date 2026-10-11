/*
 * aflr3_to_cgns.c - convert AFLR3 file to CGNS
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
using std::string;
using std::vector;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::set;
#ifdef _WIN32
#include <io.h>
#define ACCESS _access
#else
#include <unistd.h>
#define ACCESS access
#endif

#ifdef USE_CGNS

#include "binaryio.h"
#include "cgnslib.h"
#include "getargs.h"

#ifndef CGNS_ENUMV
# define CGNS_ENUMV(V) V
# define CGNS_ENUMT(T) T
#endif

int nCoords = 0;
void *xCoord, *yCoord, *zCoord;

typedef struct {
    int id, bctype;
    int start, end;
} FACESET;

typedef struct {
    int id;
    int nodes[4];
} FACE;

int nTris = 0;
FACE *Tris;
int nTriSets = 0;
FACESET *TriSets;

int nQuads = 0;
FACE *Quads;
int nQuadSets = 0;
FACESET *QuadSets;

int nTets = 0;
cgsize_t *Tets;

int nPyras = 0;
cgsize_t *Pyras;

int nPrisms = 0;
cgsize_t *Prisms;

int nHexas = 0;
cgsize_t *Hexas;

static char options[] = "48fuslb";

static const char *usgmsg[] = {
    "usage  : aflr3_to_cgns [options] AFLR3file [CGNSfile]",
    "options:",
    "   -4       = coordinates are real*4 (float)",
    "   -8       = coordinates are real*8 (double - default)",
    "   -f       = formatted (ASCII) file format",
    "   -u       = Fortran unformatted file format",
    "   -s       = Fortran stream format (binary - default)",
    "   -l       = little-endian format",
    "   -b       = big-endian format (default)",
    "default is big-endian binary (Fortran stream) AFLR3 file",
    NULL
};

int is_single = 0;

/*----------------------------------------------------------------------*/

static int sort_faces(const void *v1, const void *v2)
{
    const FACE *f1 = (const FACE *)v1;
    const FACE *f2 = (const FACE *)v2;
    return f1->id - f2->id;
}

/*----------------------------------------------------------------------*/

static void build_trisets()
{
    int n, id, ns, nf;

    for (n = 1; n < nTris; n++) {
        if (Tris[n].id < Tris[n-1].id) {
            qsort(Tris, nTris, sizeof(FACE), sort_faces);
            break;
        }
    }

    id = Tris[0].id;
    nTriSets = 1;
    for (n = 1; n < nTris; n++) {
        if (Tris[n].id != id) {
            id = Tris[n].id;
            nTriSets++;
        }
    }

    TriSets = (FACESET *)malloc(nTriSets * sizeof(FACESET));
    if (TriSets == NULL) {
        fprintf(stderr, "malloc failed for tri sets\n");
        exit(1);
    }

    id = Tris[0].id;
    ns = 0;
    nf = 0;
    for (n = 1; n < nTris; n++) {
        if (Tris[n].id != id) {
            TriSets[ns].id = id;
            TriSets[ns].bctype = -1;
            TriSets[ns].start = nf;
            TriSets[ns].end = n - 1;
            ns++;
            id = Tris[n].id;
            nf = n;
        }
    }
    TriSets[ns].id = id;
    TriSets[ns].bctype = -1;
    TriSets[ns].start = nf;
    TriSets[ns].end = nTris - 1;
}

/*----------------------------------------------------------------------*/

static void build_quadsets()
{
    int n, id, ns, nf;

    for (n = 1; n < nQuads; n++) {
        if (Quads[n].id < Quads[n-1].id) {
            qsort(Quads, nQuads, sizeof(FACE), sort_faces);
            break;
        }
    }

    id = Quads[0].id;
    nQuadSets = 1;
    for (n = 1; n < nQuads; n++) {
        if (Quads[n].id != id) {
            id = Quads[n].id;
            nQuadSets++;
        }
    }

    QuadSets = (FACESET *)malloc(nQuadSets * sizeof(FACESET));
    if (QuadSets == NULL) {
        fprintf(stderr, "malloc failed for quad sets\n");
        exit(1);
    }

    id = Quads[0].id;
    ns = 0;
    nf = 0;
    for (n = 1; n < nQuads; n++) {
        if (Quads[n].id != id) {
            QuadSets[ns].id = id;
            QuadSets[ns].bctype = -1;
            QuadSets[ns].start = nf;
            QuadSets[ns].end = n - 1;
            ns++;
            id = Quads[n].id;
            nf = n;
        }
    }
    QuadSets[ns].id = id;
    QuadSets[ns].bctype = -1;
    QuadSets[ns].start = nf;
    QuadSets[ns].end = nQuads - 1;
}

/*----------------------------------------------------------------------*/

static void read_ugrid(char *filename, int flags)
{
    int i, n, nn, counts[7], nodes[8];
    BINARYIO *bf;

    /* open file */

    bf = bf_open(filename, flags);
    if (bf == NULL) {
        fprintf(stderr, "error opening <%s> for reading\n", filename);
        exit(1);
    }

    /* get counts */

    if (7 != bf_getints(bf, 7, counts)) {
        fprintf(stderr, "error reading counts\n");
        exit(1);
    }
    nCoords = counts[0];
    nTris   = counts[1];
    nQuads  = counts[2];
    nTets   = counts[3];
    nPyras  = counts[4];
    nPrisms = counts[5];
    nHexas  = counts[6];

    /* read coordinates */

    if (is_single) {
        float coord[3];
        float *x, *y, *z;
        x = (float *)malloc(nCoords * sizeof(float));
        y = (float *)malloc(nCoords * sizeof(float));
        z = (float *)malloc(nCoords * sizeof(float));
        if (x == NULL || y == NULL || z == NULL) {
            fprintf(stderr, "malloc failed for coordinates\n");
            exit(1);
        }
        for (n = 0; n < nCoords; n++) {
            if (3 != bf_getfloats(bf, 3, coord)) {
                fprintf(stderr, "error reading coordinates\n");
                exit(1);
            }
            x[n] = coord[0];
            y[n] = coord[1];
            z[n] = coord[2];
        }
        xCoord = (void *)x;
        yCoord = (void *)y;
        zCoord = (void *)z;
    }
    else {
        double coord[3];
        double *x, *y, *z;
        x = (double *)malloc(nCoords * sizeof(double));
        y = (double *)malloc(nCoords * sizeof(double));
        z = (double *)malloc(nCoords * sizeof(double));
        if (x == NULL || y == NULL || z == NULL) {
            fprintf(stderr, "malloc failed for coordinates\n");
            exit(1);
        }
        for (n = 0; n < nCoords; n++) {
            if (3 != bf_getdoubles(bf, 3, coord)) {
                fprintf(stderr, "error reading coordinates\n");
                exit(1);
            }
            x[n] = coord[0];
            y[n] = coord[1];
            z[n] = coord[2];
        }
        xCoord = (void *)x;
        yCoord = (void *)y;
        zCoord = (void *)z;
    }

    /* read boundary tris */

    if (nTris) {
        Tris = (FACE *)malloc(nTris * sizeof(FACE));
        if (Tris == NULL) {
            fprintf(stderr, "malloc failed for boundary tris\n");
            exit(1);
        }
        for (n = 0; n < nTris; n++) {
            Tris[n].id = 0;
            if (3 != bf_getints(bf, 3, Tris[n].nodes)) {
                fprintf(stderr, "error reading tris\n");
                exit(1);
            }
        }
    }

    /* read boundary quads */

    if (nQuads) {
        Quads = (FACE *)malloc(nQuads * sizeof(FACE));
        if (Quads == NULL) {
            fprintf(stderr, "malloc failed for boundary quads\n");
            exit(1);
        }
        for (n = 0; n < nQuads; n++) {
            Quads[n].id = 0;
            if (4 != bf_getints(bf, 4, Quads[n].nodes)) {
                fprintf(stderr, "error reading quads\n");
                exit(1);
            }
        }
    }

    /* read surface IDs */

    if (nTris) {
        for (n = 0; n < nTris; n++) {
            if (1 != bf_getints(bf, 1, &Tris[n].id)) {
                fprintf(stderr, "error reading tri surface IDs\n");
                exit(1);
            }
        }
        build_trisets();
    }
    if (nQuads) {
        for (n = 0; n < nQuads; n++) {
            if (1 != bf_getints(bf, 1, &Quads[n].id)) {
                fprintf(stderr, "error reading quad surface IDs\n");
                exit(1);
            }
        }
        build_quadsets();
    }

    /* read tets */

    if (nTets) {
        Tets = (cgsize_t *)malloc(nTets * 4 * sizeof(cgsize_t));
        if (Tets == NULL) {
            fprintf(stderr, "malloc failed for tets\n");
            exit(1);
        }
        for (nn = 0, n = 0; n < nTets; n++) {
            if (4 != bf_getints(bf, 4, nodes)) {
                fprintf(stderr, "error reading tets\n");
                exit(1);
            }
            for (i = 0; i < 4; i++)
                Tets[nn++] = nodes[i];
        }
    }

    /* read pyramids */

    if (nPyras) {
        Pyras = (cgsize_t *)malloc(nPyras * 5 * sizeof(cgsize_t));
        if (Pyras == NULL) {
            fprintf(stderr, "malloc failed for pyramids\n");
            exit(1);
        }
        for (nn = 0, n = 0; n < nPyras; n++) {
            if (5 != bf_getints(bf, 5, nodes)) {
                fprintf(stderr, "error reading pyramids\n");
                exit(1);
            }
            Pyras[nn++] = nodes[0];
            Pyras[nn++] = nodes[3];
            Pyras[nn++] = nodes[4];
            Pyras[nn++] = nodes[1];
            Pyras[nn++] = nodes[2];
        }
    }

    /* read prisms */

    if (nPrisms) {
        Prisms = (cgsize_t *)malloc(nPrisms * 6 * sizeof(cgsize_t));
        if (Prisms == NULL) {
            fprintf(stderr, "malloc failed for prisms\n");
            exit(1);
        }
        for (nn = 0, n = 0; n < nPrisms; n++) {
            if (6 != bf_getints(bf, 6, nodes)) {
                fprintf(stderr, "error reading prisms\n");
                exit(1);
            }
            for (i = 0; i < 6; i++)
                Prisms[nn++] = nodes[i];
        }
    }

    /* read hexas */

    if (nHexas) {
        Hexas = (cgsize_t *)malloc(nHexas * 8 * sizeof(cgsize_t));
        if (Hexas == NULL) {
            fprintf(stderr, "malloc failed for hexas\n");
            exit(1);
        }
        for (nn = 0, n = 0; n < nHexas; n++) {
            if (8 != bf_getints(bf, 8, nodes)) {
                fprintf(stderr, "error reading hexas\n");
                exit(1);
            }
            for (i = 0; i < 8; i++)
                Hexas[nn++] = nodes[i];
        }
    }
    
    bf_close(bf);
}

/*----------------------------------------------------------------------*/

static void read_mapbc(char *filename)
{
    FILE *fp;
    int i, n, numBC, id, bc;
    char buf[256];

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "error opening <%s> for reading\n", filename);
        return;
    }
    if (NULL == fgets(buf, sizeof(buf), fp) ||
        1 != sscanf(buf, "%d", &numBC)) {
        fprintf(stderr, "error reading mapbc data\n");
        exit(1);
    }
    for (n = 0; n < numBC; n++) {
        if (NULL == fgets(buf, sizeof(buf), fp) ||
            2 != sscanf(buf, "%d %d", &id, &bc)) {
            fprintf(stderr, "error reading bc data\n");
            exit(1);
        }
        for (i = 0; i < nTriSets; i++) {
            if (TriSets[i].id == id)
                TriSets[i].bctype = bc;
        }
        for (i = 0; i < nQuadSets; i++) {
            if (QuadSets[i].id == id)
                QuadSets[i].bctype = bc;
        }
    }
    fclose(fp);
}

/*----------------------------------------------------------------------*/
 struct BC_descriptor {
    std::string name ;
    int id ;
    bool Trans ;
  } ;
 vector<BC_descriptor> readTags(string filename) {
    vector<BC_descriptor> bcs ;
    ifstream file(filename.c_str(),ios::in) ;
    if(file.fail()) {
      return bcs ;
    }
    char buf[1024] ;
    buf[1023] = '\0' ;
    file.getline(buf,1023) ;
    if(strncmp(buf,"# Generated",11) == 0) {
      cout << "Tagfile: " << buf << endl ;
      // New tag file, dump first two lines
      file.getline(buf,1023) ;
      file.getline(buf,1023) ;
      int nsp = 0 ;
      for(int i=0;i<1023;++i) {
        if(buf[i] == '\0')
          break ;
        if(isspace(buf[i])) {
          nsp++ ;
          while(isspace(buf[i]))
            i++ ;
        }
        if(strncmp(&buf[i],"Trans",5) == 0)
          break ;
      }
          
      int n2trans = nsp-1 ;
      while(file.peek() != EOF) {
        int id = -1;
        file >> id ;
        if(id == -1)
          break ;
        string name ;
        file >> name ;
        if(name == "")
          break ;
        string tmp = name ;
        size_t nsz = name.size() ;
        if(!(name[0] >= 'a' && name[0] <= 'z') &&
           !(name[0] >= 'A' && name[0] <= 'Z'))
          name[0] = '_' ;
        for(size_t i=1;i<nsz;++i) 
          if(!(name[i] >= 'a' && name[i] <= 'z') &&
             !(name[i] >= 'A' && name[i] <= 'Z') &&
             !(name[i] >= '0' && name[i] <= '9'))
            name[i] = '_' ;
        if(tmp != name) 
          cerr << "Renaming tag '" << tmp << "' to '" << name << "'!" << endl ;
        
        bool trans ;
        for(int i=0;i<n2trans;++i)
          file >> trans ;
        BC_descriptor BC ;
        BC.id = id ;
        BC.name = name ;
        BC.Trans = trans ;
        file.getline(buf,1023) ;
        bcs.push_back(BC) ;
        while(file.peek() != EOF && (isspace(file.peek())))
          file.get() ;
        while(file.peek() == '#') { // eat comments
          file.getline(buf,1023) ;
          while(file.peek() != EOF && (isspace(file.peek())))
            file.get() ;
        }
      }
    } else {
      if(strncmp(buf,"#ID:",4) != 0) {
        cerr << "unable to interpret tag file format" << endl ;
        cerr << "first line: " << buf << endl ;
        return bcs ;
      }
      int nsp = 0 ;
      for(int i=0;i<120;++i) {
        if(buf[i] == '\0')
          break ;
        if(isspace(buf[i])) {
          nsp++ ;
          while(isspace(buf[i]))
            i++ ;
        }
        if(strncmp(&buf[i],"Trans",5) == 0)
          break ;
      }
      int n2trans = nsp-1 ;
      while(file.peek() != EOF) {
        while(file.peek() != EOF&&file.peek() != '#')
          file.get() ;

        char c ;
        c = file.get() ;
        if(c != '#') {
          cerr << "expected to get a hash while reading '" << filename << "'"
               << endl ;
          bcs.clear() ;
          return bcs ;
        }
        int id = 0 ;
        file >> id ;
        c = file.get() ;
        if(c != ':') {
          cerr << "expected to get a colon while reading '" << filename << "'"
               << endl ;
          bcs.clear() ;
          return bcs ;
        }
        string name ;
        file >> name ;
        string tmp = name ;
        size_t nsz = name.size() ;
        if(!(name[0] >= 'a' && name[0] <= 'z') &&
           !(name[0] >= 'A' && name[0] <= 'Z'))
          name[0] = '_' ;
        for(size_t i=1;i<nsz;++i) 
          if(!(name[i] >= 'a' && name[i] <= 'z') &&
             !(name[i] >= 'A' && name[i] <= 'Z') &&
             !(name[i] >= '0' && name[i] <= '9'))
            name[i] = '_' ;
        if(tmp != name) 
          cerr << "Renaming tag '" << tmp << "' to '" << name << "'!" << endl ;

        bool trans ;
        for(int i=0;i<n2trans;++i)
          file >> trans ;
        BC_descriptor BC ;
        BC.id = id ;
        BC.name = name ;
        BC.Trans = trans ;
        file.getline(buf,1024) ;
        bcs.push_back(BC) ;
        while(file.peek() != EOF && (isspace(file.peek())))
          file.get() ;
      }
    }
    return bcs ;
  }
/*----------------------------------------------------------------------*/

static CGNS_ENUMV(BCType_t) get_bctype (int bctype)
{
    switch (bctype) {
        case 1:
        case 2:
        case 5000:
        case 5025:
        case 5050:
            return CGNS_ENUMV(BCFarfield);
        case 3:
        case 4:
        case 5:
        case 6:
        case 4000:
            return CGNS_ENUMV(BCWallViscous);
        case 5026:
            return CGNS_ENUMV(BCExtrapolate);
        case 6661:
        case 6662:
        case 6663:
            return CGNS_ENUMV(BCSymmetryPlane);
        case 3000:
            return CGNS_ENUMV(BCWallInviscid);
        case 5051:
        case 5052:
            return CGNS_ENUMV(BCTunnelOutflow);
        case 7011:
            return CGNS_ENUMV(BCInflowSubsonic);
        case 7012:
            return CGNS_ENUMV(BCOutflowSubsonic);
        case 7031:
            return CGNS_ENUMV(BCOutflow);
        case 7036:
        case 7100:
            return CGNS_ENUMV(BCInflow);
        default:
            if (bctype >= 7 && bctype <= 18)
                return CGNS_ENUMV(BCWallInviscid);
            break;
    }
    return CGNS_ENUMV(BCTypeUserDefined);
}

/*----------------------------------------------------------------------*/

static void write_cgns(char *filename, const vector<BC_descriptor>& bcs)
{
    char name[33];
    const char *bcname;
    int cgfile, cgbase, cgzone;
    int cgcoord, cgsect, cgbc;
    int n, i, ns, ne, nemax;
    cgsize_t sizes[3], *nodes, start, end, range[2];
    CGNS_ENUMT(DataType_t) dtype;
    CGNS_ENUMT(BCType_t) bctype = CGNS_ENUMV(BCTypeUserDefined) ;

    sizes[0] = nCoords;
    sizes[1] = nTets + nPyras + nPrisms + nHexas;
    if (sizes[1] == 0) sizes[1] = nTris + nQuads;
    sizes[2] = 0;
    
    if (cg_open(filename, CG_MODE_WRITE, &cgfile) ||
        cg_base_write(cgfile, "Base", 3, 3, &cgbase) ||
        cg_zone_write(cgfile, cgbase, "Zone", sizes,
                      CGNS_ENUMV(Unstructured), &cgzone))
      cg_error_exit();
    dtype = is_single ? CGNS_ENUMV(RealSingle) : CGNS_ENUMV(RealDouble);
    if (cg_coord_write(cgfile, cgbase, cgzone, dtype,
            "CoordinateX", xCoord, &cgcoord) ||
        cg_coord_write(cgfile, cgbase, cgzone, dtype,
            "CoordinateY", yCoord, &cgcoord) ||
        cg_coord_write(cgfile, cgbase, cgzone, dtype,
            "CoordinateZ", zCoord, &cgcoord))
        cg_error_exit();
    free(xCoord);
    free(yCoord);
    free(zCoord);

    start = 1;
    if (nTets) {
        end = start + nTets - 1;
        if (cg_section_write(cgfile, cgbase, cgzone, "TetElements",
                CGNS_ENUMV(TETRA_4), start, end, 0, Tets, &cgsect))
            cg_error_exit();
        free(Tets);
        start = end + 1;
    }
    if (nPyras) {
        end = start + nPyras - 1;
        if (cg_section_write(cgfile, cgbase, cgzone, "PyraElements",
                CGNS_ENUMV(PYRA_5), start, end, 0, Pyras, &cgsect))
            cg_error_exit();
        free(Pyras);
        start = end + 1;
    }
    if (nPrisms) {
        end = start + nPrisms - 1;
        if (cg_section_write(cgfile, cgbase, cgzone, "PentaElements",
                CGNS_ENUMV(PENTA_6), start, end, 0, Prisms, &cgsect))
            cg_error_exit();
        free(Prisms);
        start = end + 1;
    }
    if (nHexas) {
        end = start + nHexas - 1;
        if (cg_section_write(cgfile, cgbase, cgzone, "HexaElements",
                CGNS_ENUMV(HEXA_8), start, end, 0, Hexas, &cgsect))
            cg_error_exit();
        free(Hexas);
        start = end + 1;
    }

    if (nTris + nQuads == 0) {
        cg_close(cgfile);
        return;
    }

    nemax = 0;
    for (ns = 0; ns < nTriSets; ns++) {
        ne = 3 * (TriSets[ns].end - TriSets[ns].start + 1);
        if (nemax < ne) nemax = ne;
    }
    for (ns = 0; ns < nQuadSets; ns++) {
        ne = 4 * (QuadSets[ns].end - QuadSets[ns].start + 1);
        if (nemax < ne) nemax = ne;
    }
    nodes = (cgsize_t *)malloc(nemax * sizeof(cgsize_t));
    if (nodes == NULL) {
        fprintf(stderr, "malloc failed for tri/quad array\n");
        exit(1);
    }

    /* collect all the ids */
    set<int> ids;
    for (ns = 0; ns < nTriSets; ns++) ids.insert(TriSets[ns].id);
    for (ns = 0; ns < nQuadSets; ns++)ids.insert(QuadSets[ns].id);
    if(!bcs.empty() && ids.size() != bcs.size()){
      cerr<<" ERROR: number of boundaries in .tags file does not match that in .ugrid file" <<endl;
      cg_error_exit(); 
    }

    /* for each boundary, first write out its triset, then its quadset*/
   
    for(set<int>::const_iterator si=ids.begin(); si != ids.end(); si++){
      int bc_id = *si;
      /* write tris */

      for (ns = 0; ns < nTriSets; ns++) {
        if(TriSets[ns].id == bc_id){
          n = 0;
          for (ne = TriSets[ns].start; ne <= TriSets[ns].end; ne++) {
            for (i = 0; i < 3; i++)
              nodes[n++] = Tris[ne].nodes[i];
          }
          end = start + TriSets[ns].end - TriSets[ns].start;
          sprintf(name, "TriElements %d", ns+1);
          if (cg_section_write(cgfile, cgbase, cgzone, name,
                               CGNS_ENUMV(TRI_3), start, end, 0, nodes, &cgsect))
            cg_error_exit();
          TriSets[ns].start = (int)start;
          TriSets[ns].end = (int)end;
          start = end + 1;
          break;
        }
      }

      /* write quads */

      for (ns = 0; ns < nQuadSets; ns++) {
        if(QuadSets[ns].id == bc_id){
          n = 0;
          for (ne = QuadSets[ns].start; ne <= QuadSets[ns].end; ne++) {
            for (i = 0; i < 4; i++)
              nodes[n++] = Quads[ne].nodes[i];
          }
          end = start + QuadSets[ns].end - QuadSets[ns].start;
          sprintf(name, "QuadElements %d", ns+1+nTriSets);
          if (cg_section_write(cgfile, cgbase, cgzone, name,
                               CGNS_ENUMV(QUAD_4), start, end, 0, nodes, &cgsect))
            cg_error_exit();
          QuadSets[ns].start = (int)start;
          QuadSets[ns].end = (int)end;
          start = end + 1;
          break;
        }
      }
    }
    free(nodes);

    /* write BCs */
    for(set<int>::const_iterator si=ids.begin(); si != ids.end(); si++){
      int bc_id = *si;
      bool name_found = false;
      if(!bcs.empty()){
        for(unsigned int i= 0; i< bcs.size(); i++){
          if(bcs[i].id == bc_id){
            // std::strcpy (name, bcs[i].name.c_str());
            sprintf(name, "%s", bcs[i].name.c_str()); 
            name_found = true;
            break;
          }
        }
      }
      bool tri_found = false;
      for (ns = 0; ns < nTriSets; ns++) {
        if(TriSets[ns].id == bc_id){
          bctype = get_bctype(TriSets[ns].bctype);
          if(!name_found){
            bcname = cg_BCTypeName(bctype);
            if (0 == strncmp(bcname, "BC", 2)) bcname += 2;
            sprintf(name, "%s %d", bcname, ns+1);
            name_found=true;
          }
          cout<< "bc name: " << string(name) << endl;
          range[0] = TriSets[ns].start;
          range[1] = TriSets[ns].end;
          tri_found = true;
          break;
        }
      }

      bool quad_found = false;
      for (ns = 0; ns < nQuadSets; ns++) {
        if(QuadSets[ns].id == bc_id){
           bctype = get_bctype(QuadSets[ns].bctype);
           if(!name_found){
             bcname = cg_BCTypeName(bctype);
             if (0 == strncmp(bcname, "BC", 2)) bcname += 2;
             sprintf(name, "%s %d", bcname, ns+1+nTriSets);
             name_found = true;
           }
           if(tri_found){
             range[1] = QuadSets[ns].end;
           }else{
             range[0] = QuadSets[ns].start;
             range[1] = QuadSets[ns].end;
           }
           quad_found = true;
           break;
        }
      }
      if(!name_found) {
        cerr << " ERROR: boundary name not found for bc id " << bc_id << endl;
        cg_error_exit();
      }
      if(!(tri_found || quad_found)){
        cerr<<"ERROR: not face set defined for bc id " << bc_id << endl;
        cg_error_exit();
      }
          
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

    
    cg_close(cgfile);
}

/*======================================================================*/

int main (int argc, char *argv[])
{
    int n, flags = MACH_IEEE;
    char *p, *tail, *basename;

    if (argc < 2)
        print_usage(usgmsg, NULL);

    /* get options */

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case '4':
                is_single = 1;
                break;
            case '8':
                is_single = 0;
                break;
            case 'f':
                flags &= ~OPEN_FORTRAN;
                flags |= OPEN_ASCII;
                break;
            case 'u':
                flags &= ~OPEN_ASCII;
                flags |= OPEN_FORTRAN;
                break;
            case 's':
                flags &= ~(OPEN_ASCII | OPEN_FORTRAN);
                break;
            case 'l':
                flags &= ~MACH_UNKNOWN;
                flags |= MACH_BSIEEE;
                break;
            case 'b':
                flags &= ~MACH_UNKNOWN;
                flags |= MACH_IEEE;
                break;
        }
    }

    if (argind >= argc)
        print_usage(usgmsg, "AFLR3 filename not specified\n");

    basename = (char *)malloc(strlen(argv[argind]) + 7);
    if (basename == NULL) {
        fprintf(stderr, "malloc failed for filename\n");
        return 1;
    }
    strcpy(basename, argv[argind++]);
    
    /* read mesh */

    printf("reading grid from \"%s\"\n", basename);
    fflush(stdout);
    read_ugrid(basename, flags);
    printf("nCoords = %d\n", nCoords);
    if (nTris)   printf("nTris   = %d\n", nTris);
    if (nQuads)  printf("nQuads  = %d\n", nQuads);
    if (nTets)   printf("nTets   = %d\n", nTets);
    if (nPyras)  printf("nPyras  = %d\n", nPyras);
    if (nPrisms) printf("nPrisms = %d\n", nPrisms);
    if (nHexas)  printf("nHexas  = %d\n", nHexas);



    /* check for tags file */
    vector<BC_descriptor> bcs ;
    {
      string tmpname =  string(basename);
      std::size_t pos = tmpname.find('.');
      
      string tagsfile = tmpname.substr(0, pos) + ".tags" ;
      bcs = readTags(tagsfile) ;
      if(bcs.size() == 0) {
        cerr << "unable to read '" << tagsfile << "'" << endl ;
      }
      
      // for(size_t i=0;i<bcs.size();++i)
      //   if(bcs[i].Trans)
      //     transsurf.push_back(bcs[i].id) ;

      cout << "boundary faces:"<<endl; ;
      for(size_t i=0;i<bcs.size();++i)
        cout << ' ' <<  bcs[i].id << " " << bcs[i].name<< endl;
      
    }
    /* check for mapbc file */

    tail = strrchr(basename, '/');
#ifdef _WIN32
    if (tail == NULL) tail = strrchr(basename, '\\');
#endif
    if (tail == NULL)
        tail = basename;
    else
        tail++;
    if ((p = strrchr(tail, '.')) == NULL)
        p = basename + strlen(basename);
 
    strcpy(p, ".mapbc");
    n = ACCESS(basename, 0);
    if (n) {
        *p = 0;
        if ((p = strrchr(tail, '.')) != NULL) {
            strcpy(p, ".mapbc");
            n = ACCESS(basename, 0);
        }
    }
    if (0 == n) {
        printf("reading boundary conditions from \"%s\"\n", basename);
        fflush(stdout);
        read_mapbc(basename);
    }
    
    
    /* open CGNS file */

    if (argind < argc) {
        p = argv[argind];
    }
    else {
        strcpy(p, ".cgns");
        p = basename;
    }
    printf ("writing CGNS file to \"%s\"\n", p);
    write_cgns(p, bcs);

    return 0;
}

#else
int main(int ac, char *av[]) {
  fprintf(stderr,"Loci not compiled with CGNS support enabled! This utility cannot work!\n") ;
  return -1 ;
}
#endif
