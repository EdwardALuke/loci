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
//#define VERBOSE
#include <Loci>
#include <LociGridReaders.h>
#include <Tools/tools.h>
#include <map>
#include "pnn.h"
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
#include <algorithm>
using std::sort ;
using std::unique ;
#include <Tools/stream.h>

#include "dist_tools.h"
using std::cout ;

#include <iomanip>
#define vect3d vector3d<double>

namespace Loci{
  extern int getKeyDomain(entitySet dom, fact_db::distribute_infoP dist,
                          MPI_Comm comm) ;
  extern  bool useDomainKeySpaces  ;

  /// Map from local numbering to input file numbering.
  /// Assumes the union of nodes on all processors will be either all the
  /// nodes, all the faces, or all the cells. i.e., one interval in file
  /// numbering.
  storeRepP get_node_remap(fact_db &facts, entitySet nodes, int keyspace) {
    REPORTMEM() ;
    if(MPI_processes == 1) {
      int minNode = nodes.Min() ;

      Map nm ;
      nm.allocate(nodes) ;
      FORALL(nodes,nd) {
        nm[nd] = nd - minNode + 1 ;
      } ENDFORALL ;
      REPORTMEM() ;
      return nm.Rep() ;
    }

    vector<entitySet> init_ptn = facts.get_init_ptn(keyspace) ;
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;
    dMap g2f ;
    g2f = df->g2fv[keyspace].Rep() ;

    entitySet gnodes = l2g.image(nodes&l2g.domain()) ;

    Map newnum ;
    newnum.allocate(nodes) ;

    // Expand g2f to include clone regions
    entitySet out_of_dom = gnodes - init_ptn[MPI_rank] ;
    g2f.setRep(MapRepP(g2f.Rep())->expand(out_of_dom, init_ptn)) ;

    entitySet fnodes = g2f.image(gnodes) ;
    int minNode_local = fnodes.Min() ;
    int minNode = minNode_local ;
    MPI_Allreduce(&minNode_local,&minNode,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;

    FORALL(nodes,i) {
      newnum[i] = g2f[l2g[i]]-minNode+1 ;
    } ENDFORALL ;
    REPORTMEM() ;
    return newnum.Rep() ;
  }

  /// Return a map from local numbering to output file numbering.
  /// The file number is assumed to start with 1 and the nodes will be
  /// written out in the global ordering.
  ///
  /// @param[in] facts Fact database.
  /// @param[in] nodes The local store domain that needs to be output.
  storeRepP get_output_node_remap(fact_db &facts, entitySet nodes) {
    int p = MPI_processes ;
    if(p == 1) {
      int index = 1;

      Map nm ;
      nm.allocate(nodes) ;
      FORALL(nodes,nd) {
        nm[nd] = index++ ;
      } ENDFORALL ;
      return nm.Rep() ;
    }
#ifdef VERBOSE
    stopWatch s ;
    s.start() ;
#endif
    // when working in parallel we will have data scattered over
    // processors, so we will have to sort this data in a distributed
    // fashion in order to find the file numbering of the condensed
    // global numbered set
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;

    // gnodes is the global numbering for entities that this processor
    // references (because it is the result of an image, there will
    // be duplicates).  We map to the global numbering since this will
    // match the ordering of the vector write routines
    entitySet gnodes = l2g.image(nodes&l2g.domain()) ;
    if(nodes.size() != gnodes.size()){
      debugout << "ERROR: l2g.domain is smaller than dom "
               << "in get_output_node_remap " << endl ;
    }
    // extract a sorted list of global entities accessed by this processor
#ifdef VERBOSE
    debugout << "time to get gnodes = " << s.stop() << endl ;
    s.start() ;
#endif
    int gsz = gnodes.size() ;
    int color = (gsz>0)?1:0 ;
    MPI_Comm groupcomm ;
    MPI_Comm_split(MPI_COMM_WORLD,color,p,&groupcomm) ;
    if(color==0) {
      MPI_Comm_free(&groupcomm) ;
#ifdef VERBOSE
      debugout << "time to split comm and return = " << s.stop() << endl ;
      s.start() ;
#endif
      Map newnum ;
      return newnum.Rep() ;
    }
#ifdef VERBOSE
    debugout << "time to split comm = " << s.stop() << endl ;
    s.start() ;
#endif
    MPI_Comm_size(groupcomm,&p) ;
    vector<int> gnodelistg(gsz) ;
    int cnt = 0 ;

    FORALL(gnodes,ii) {
      gnodelistg[cnt] = ii ;
      cnt++ ;
    } ENDFORALL ;
    if(cnt != gsz) {
      cerr << "cnt=" << cnt << "gsz=" << gsz << endl ;
    }
    // form a partition of the global numbered entities across processors
    int lmaxnode = gnodes.Max() ;
    int lminnode = gnodes.Min() ;
    int maxnode,minnode ;
    MPI_Allreduce(&lmaxnode,&maxnode,1,MPI_INT,MPI_MAX,groupcomm) ;
    MPI_Allreduce(&lminnode,&minnode,1,MPI_INT,MPI_MIN,groupcomm) ;
    int delta = max((1+maxnode-minnode)/p,1024) ;
    // find out what processor owns each entity
    vector<int> sendsz(p,0) ;
    for(int i=0;i<gsz;++i) {
      int proc = min((gnodelistg[i]-minnode)/delta,p-1) ;
      sendsz[proc]++ ;
    }
#ifdef VERBOSE
    debugout << "time to form partition = " << s.stop() << endl ;
    s.start() ;
#endif
    // transpose sendsz to find out how much each processor partition recvs
    vector<int> recvsz(p,0) ;
    MPI_Alltoall(&sendsz[0],1,MPI_INT,&recvsz[0],1,MPI_INT,groupcomm) ;
    // setup recv buffers to communicate data
    int recvbufsz = 0;
    for(int i=0;i<p;++i) {
      recvbufsz+= recvsz[i] ;
    }
    vector<int> rdispls(p),sdispls(p) ;
    rdispls[0] = 0 ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+recvsz[i-1] ;
      sdispls[i] = sdispls[i-1]+sendsz[i-1] ;
    }
    vector<int> recvdata(recvbufsz) ;
    MPI_Alltoallv(&gnodelistg[0],&sendsz[0],&sdispls[0],MPI_INT,
		              &recvdata[0],&recvsz[0],&rdispls[0],MPI_INT,
		              groupcomm) ;
#ifdef VERBOSE
    debugout << "distributed data local size = " << recvdata.size() << endl;
    debugout << "time to distribute set data = " << s.stop() << endl ;
    s.stop() ;
#endif
    // Now sort a copy of the received data so the local data is in order
    // and duplicates are removed
    vector<int> recvdatac = recvdata ;
    sort(recvdatac.begin(),recvdatac.end()) ;
    vector<int>::iterator last = std::unique(recvdatac.begin(),
		      recvdatac.end()) ;
    recvdatac.erase(last,recvdatac.end()) ;
    // now calculate the starting index for each processor
    int ssz = recvdatac.size() ;
    vector<int> sortsz(p) ;
    int roff = 0 ;
    MPI_Scan(&ssz,&roff,1,MPI_INT,MPI_SUM,groupcomm) ;
    // adjust roff so that it is the count at the beginning of this
    // processors segment
    roff -= ssz ;
    // we are numbering the first entry with 1
    roff += 1 ;
    // Now compute a map from our global numbering to the file numbering
    // which is contiguously numbered from 1
    std::map<int,int> filemap ;
    for(int i=0;i<ssz;++i) {
      filemap[recvdatac[i]] = i+roff ;
    }
    // Now send back the file numbering
    vector<int> filenosend(recvbufsz,-1) ;
    for(int i=0;i<recvbufsz;++i) {
      filenosend[i] = filemap[recvdata[i]] ;
    }
    vector<int> gnodefnum(gnodelistg.size(),-1) ;
    // note now recvsz is sending and sendsz is recving
    MPI_Alltoallv(&filenosend[0],&recvsz[0],&rdispls[0],MPI_INT,
		              &gnodefnum[0],&sendsz[0],&sdispls[0],MPI_INT,
		              groupcomm) ;
#ifdef VERBOSE
    debugout << "time to return index = " << s.stop() << endl ;
    s.start() ;
#endif
    MPI_Comm_free(&groupcomm) ;

    // now create global 2 file map by extracting data from gnodefnum
    std::map<int,int> g2f;
    for(size_t i=0;i<gnodelistg.size();++i)  {
      g2f[gnodelistg[i]] = gnodefnum[i] ;
    }

    // join the two maps to get local 2 file mapping
    Map newnum ;
    newnum.allocate(nodes) ;

    FORALL(nodes,i) {
      newnum[i] = g2f[l2g[i]];
    } ENDFORALL ;
#ifdef VERBOSE
    debugout << "time to form map = " << s.stop() << endl ;
#endif
    return newnum.Rep() ;
  }

  /// Returns the classification of a cell based on face topology.
  ///
  /// Classification is:
  ///  0 for 4 triangles only i.e. tetrahedron.
  ///  1 if hex_test i.e. hexahedron.
  ///  2 if prism_test i.e. triangular prism.
  ///  3 for 4 triangles + 1 quad i.e. pyramid.
  ///  4 for anything else i.e. general/other polyhedron.
  ///
  /// @param[in] faces faces that make up the cell.
  /// @param[in] nfaces number of faces.
  /// @param[in] face2node mapping from face to nodes.
  int classify_cell(Entity *faces, int nfaces, const_multiMap &face2node) {
    int num_triangles = 0 ;
    int num_quads = 0 ;
    int num_others = 0 ;
    int triangle_nodes[3][2] ;
    for(int f=0;f<nfaces;++f) {
      Entity fc = faces[f] ;
      int count = face2node[fc].size() ;
      if(count == 3) {
        if(num_triangles < 2) {
          triangle_nodes[0][num_triangles] = face2node[fc][0] ;
          triangle_nodes[1][num_triangles] = face2node[fc][1] ;
          triangle_nodes[2][num_triangles] = face2node[fc][2] ;
        }
        num_triangles++ ;
      } else if(count == 4) {
        num_quads++ ;
      } else {
        num_others++ ;
      }
    }
    bool prism_test = false ;

    if((num_triangles == 2) && (num_quads == 3) && (num_others == 0)) {
      prism_test = true ;
      for(int i=0;i<3;++i) {
        for(int j=0;j<3;++j) {
          if(triangle_nodes[i][0] == triangle_nodes[j][1]) {
            prism_test = false ;
          }
        }
      }
    }

    bool hex_test = false ;
    if( (num_triangles == 0) && (num_quads == 6) && (num_others == 0)) {
      const Entity ef = faces[0] ;
      int count = 0 ;
      for(int fj = 1;fj<nfaces;++fj) {
        Entity ej = faces[fj] ;
        bool find = false ;
        for(int i=0;i<4;++i) {
          for(int j=0;j<4;++j) {
            if(face2node[ef][i] == face2node[ej][j]) {
              find = true ;
            }
          }
        }
        if(find) {
          count++ ;
        }
      }
      if(count == 4) {
        hex_test = true ;
      }
    }

    // new classification code
    if( (num_triangles == 4) && (num_quads == 0) && (num_others == 0)) {
      return 0 ;
    } else if( hex_test ) {
      return 1 ;
    } else if( prism_test ) {
      return 2 ;
    } else if( (num_triangles == 4) && (num_quads == 1) && (num_others == 0)) {
      return 3 ;
    }
    return 4 ;
  }

  /// Fills in the tet element connectivity from face connectivity.
  /// @param[out] tet The tet element connectivity to fill in.
  /// @param[in] tri_faces The triangular faces that make up the tet.
  void fillTet(Array<int,4> &tet, Array<int,3> *tri_faces) {
    tet[0] = tri_faces[0][2] ;
    tet[1] = tri_faces[0][1] ;
    tet[2] = tri_faces[0][0] ;
    for(int i=0;i<3;++i) {
      tet[3] = tri_faces[1][i] ;
      if(tet[3] != tet[0] && tet[3] != tet[1] && tet[3] != tet[2]) {
        return ;
      }
    }
    cerr << "unable to form valid tet!" << endl ;
  }

  /// Fills in the hex element connectivity from face connectivity.
  /// @param[out] hex The hex element connectivity to fill in.
  /// @param[in] quad_faces The quadrilateral faces that make up the hex.
  void fillHex(Array<int,8> &hex, Array<int,4> *quad_faces) {
    int quad_id[6] ;
    for(int i=0;i<6;++i) {
      quad_id[i] = i ;
    }
    bool degenerate = quad_faces[quad_id[0]][0] == quad_faces[quad_id[0]][3];
    for(int j=0;j<3;++j)
      if(quad_faces[quad_id[0]][j] == quad_faces[quad_id[0]][j+1])
        degenerate = true ;
    if(degenerate) {
      for(int i=1;i<6;++i) {
        degenerate = quad_faces[quad_id[i]][0] == quad_faces[quad_id[i]][3];
        for(int j=0;j<3;++j) {
          if(quad_faces[quad_id[i]][j] == quad_faces[quad_id[i]][j+1]) {
            degenerate = true ;
          }
        }
        if(!degenerate) {
          std::swap(quad_id[i],quad_id[0]) ;
          break ;
        }
      }
    }
    hex[0] = quad_faces[quad_id[0]][3] ;
    hex[1] = quad_faces[quad_id[0]][2] ;
    hex[2] = quad_faces[quad_id[0]][1] ;
    hex[3] = quad_faces[quad_id[0]][0] ;
    hex[4] = hex[0] ;
    hex[5] = hex[1] ;
    hex[6] = hex[2] ;
    hex[7] = hex[3] ;
    for(int i = 0; i < 4; i+=2) {
      int n1 = hex[i] ;
      int n2 = hex[i+1] ;

      int cnt = 0 ;
      for(int j=1; j<6; ++j) {
        for(int k=0;k<4;++k) {
          int f1 = quad_faces[quad_id[j]][k] ;
          int f2 = quad_faces[quad_id[j]][(k+1)%4] ;
          if((f1 == n1 && f2 == n2)) {
            hex[i+4] = quad_faces[quad_id[j]][(k-1+4)%4] ;
            hex[i+1+4] = quad_faces[quad_id[j]][(k+2)%4] ;
            cnt++ ;
          }
        }
      }
      if(cnt != 1) {
        cerr << "Error: Hex elem ordering screwed up " <<  endl ;
      }
    }
  }

  /// Fills in the prism element connectivity from face connectivity.
  /// @param[out] prism The prism element connectivity to fill in.
  /// @param[in] tri_faces The triangular faces that make up the prism.
  /// @param[in] quad_faces The quadrilateral faces that make up the prism.
  void fillPrism(Array<int,6> &prism, Array<int,3> *tri_faces,
                 Array<int,4> *quad_faces) {
    prism[0] = tri_faces[0][2] ;
    prism[1] = tri_faces[0][1] ;
    prism[2] = tri_faces[0][0] ;
    prism[3] = prism[0] ;
    prism[4] = prism[1] ;
    prism[5] = prism[2] ;

    int n1 = prism[0] ;
    int n2 = prism[1] ;
    int n3 = prism[2] ;
    int cnt = 0 ;
    for(int j=0;j<3;++j) {
      for(int k=0;k<4;++k) {
        int f1 = quad_faces[j][k] ;
        int f2 = quad_faces[j][(k+1)%4] ;
        if((f1 == n1 && f2 == n2)) {
          prism[3] = quad_faces[j][(k+3)%4] ;
          cnt++ ;
        }

        if((f1 == n2 && f2 == n3)) {
          prism[4] = quad_faces[j][(k+3)%4] ;
          prism[5] = quad_faces[j][(k+2)%4] ;
          cnt++ ;
        }
      }
    }
    if(cnt != 2) {
      cerr << "prism ordering screwed up" << endl ;
    }
  }

  /// Fills in the pyramid element connectivity from face connectivity.
  /// @param[out] pyramid The pyramid element connectivity to fill in.
  /// @param[in] tri_faces The triangular faces that make up the pyramid.
  /// @param[in] quad_faces The quadrilateral faces that make up the pyramid.
  void fillPyramid(Array<int,5> &pyramid, Array<int,3> *tri_faces,
                   Array<int,4> *quad_faces) {
    pyramid[0] = quad_faces[0][3] ;
    pyramid[1] = quad_faces[0][2] ;
    pyramid[2] = quad_faces[0][1] ;
    pyramid[3] = quad_faces[0][0] ;
    pyramid[4] = pyramid[0] ;
    for(int i=0;i<3;++i) {
      int nd = tri_faces[0][i] ;
      if(nd != pyramid[0] && nd != pyramid[1] &&
         nd != pyramid[2] && nd != pyramid[3]) {
        pyramid[4] = nd ;
        return ;
      }
    }
    cerr << "pyramid ordering screwed up!" << endl ;
  }

  /// Concatenates strings from all processors in the communicator.
  /// @param[in] input The input string on each processor.
  /// @param[in] comm The MPI communicator.
  /// @return The concatenated string.
  string MPIConcatStrings(string input, MPI_Comm comm) {
    int sz = input.size()+1 ;
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    int *sizes = new int[p] ;
    MPI_Allgather(&sz,1,MPI_INT,sizes,1,MPI_INT,comm) ;
    int tot = 0 ;
    for(int i=0;i<p;++i) {
      tot += sizes[i] ;
    }
    char *buf = new char[tot] ;
    int *displ = new int[p] ;
    displ[0] = 0 ;
    for(int i=1;i<p;++i) {
      displ[i] = displ[i-1]+sizes[i-1] ;
    }
    char *ibuf = new char[sz] ;
    strcpy(ibuf,input.c_str()) ;
    MPI_Allgatherv(ibuf, sz, MPI_CHAR, buf, sizes, displ, MPI_CHAR, comm) ;
    string retval ;
    for(int i=0;i<p;++i) {
      retval += string(&(buf[displ[i]])) ;
    }

    delete[] sizes ;
    delete[] buf ;
    delete[] ibuf ;
    delete[] displ ;

    return retval ;
  }

  /// Writes out the grid topology in parallel to a file.
  ///
  /// The node_id range is [1,numNodes]. The face_id range is
  /// [numNodes,numNodes+numFaces]. The cell_id range is [1,numCells]. This is
  /// because `get_node_remap()` is used to get the file number of nodes and
  /// cells, which set the start index to 1, however, g2f(l2g(ff)) is used to
  /// get the file number of faces. Since boundary variables also use the
  /// face_id, to avoid changing too many places, we keep the face_id this way.
  ///
  /// @param[in] filename The name of the output file.
  /// @param[in] upperRep The upper faces.
  /// @param[in] lowerRep The lower faces.
  /// @param[in] boundary_mapRep The boundary faces.
  /// @param[in] face2nodeRep The face to node map.
  /// @param[in] refRep The ref map.
  /// @param[in] bnamesRep The boundary names.
  /// @param[in] posRep The node positions.
  /// @param[in] localCells The local cells to write out.
  /// @param[in,out] facts The fact database.
  void parallelWriteGridTopology(const char *filename, storeRepP upperRep,
        storeRepP lowerRep, storeRepP boundary_mapRep, storeRepP face2nodeRep,
        storeRepP refRep, storeRepP bnamesRep, storeRepP posRep,
        entitySet localCells, fact_db &facts) {
    const_multiMap upper(upperRep), lower(lowerRep),
      boundary_map(boundary_mapRep),face2node(face2nodeRep) ;
    const_Map ref(refRep) ;

    store<int> elem_type ;
    elem_type.allocate(localCells) ;

    int ntets = 0 ;
    int nhexs = 0 ;
    int nprsm = 0 ;
    int npyrm = 0 ;
    //    int ngnrl = 0 ;

    // Classify Cells
    FORALL(localCells,cc) {
      int nfaces = upper[cc].size()+lower[cc].size()+boundary_map[cc].size() ;
      tmp_array<Entity> faces(nfaces) ;
      int cnt = 0 ;
      for(int i=0;i<upper[cc].size();++i) {
        faces[cnt++] = upper[cc][i] ;
      }
      for(int i=0;i<lower[cc].size();++i) {
        faces[cnt++] = lower[cc][i] ;
      }
      for(int i=0;i<boundary_map[cc].size();++i) {
        faces[cnt++] = boundary_map[cc][i] ;
      }

      elem_type[cc] = classify_cell(faces,nfaces,face2node) ;
      switch(elem_type[cc]) {
      case 0:
        ntets++; break ;
      case 1:
        nhexs++ ; break ;
      case 2:
        nprsm++ ; break ;
      case 3:
        npyrm++ ; break ;
      default:
	//        ngnrl++ ;
	break ;
      }
    } ENDFORALL ;

    // Collect 4 cell type info
    vector<Array<int,4> > tets(ntets) ;
    vector<Array<int,5> > pyrm(npyrm) ;
    vector<Array<int,6> > prsm(nprsm) ;
    vector<Array<int,8> > hexs(nhexs) ;

    //for ids
    vector<int> tets_ids(ntets) ;
    vector<int> pyrm_ids(npyrm) ;
    vector<int> prsm_ids(nprsm) ;
    vector<int> hexs_ids(nhexs) ;
    vector<int> generalCell_ids;
    int tet_no = 0 ;
    int hex_no = 0 ;
    int pyramid_no = 0 ;
    int prism_no = 0 ;

    Map node_remap ;
    node_remap = get_node_remap(facts,posRep->domain(),
                                posRep->getDomainKeySpace()) ;

    Map cell_remap;
    cell_remap = get_node_remap(facts,localCells,
                                upperRep->getDomainKeySpace()) ;

    vector<int> generalCellNfaces ;
    vector<int> generalCellNsides ;
    vector<int> generalCellNodes ;

    // Generate Cells
    FORALL(localCells,cc) {
      int nfaces = upper[cc].size()+lower[cc].size()+boundary_map[cc].size() ;
      tmp_array<int> faces(nfaces) ;
      tmp_array<int> swapface(nfaces) ;
      tmp_array<Array<int,3> > tri_faces(nfaces) ;
      tmp_array<Array<int,4> > quad_faces(nfaces) ;

      int tcnt = 0 ;
      int qcnt = 0 ;
      int nf = 0 ;
      for(int i=0;i<upper[cc].size();++i) {
        int fc = upper[cc][i] ;
        swapface[nf] = 0 ;
        faces[nf] = fc ;
        nf++ ;
        int fsz = face2node[fc].size() ;
        if(fsz == 3) {
          tri_faces[tcnt][0] = node_remap[face2node[fc][0]] ;
          tri_faces[tcnt][1] = node_remap[face2node[fc][1]] ;
          tri_faces[tcnt][2] = node_remap[face2node[fc][2]] ;
          tcnt++ ;
        }
        if(fsz == 4) {
          quad_faces[qcnt][0] = node_remap[face2node[fc][0]] ;
          quad_faces[qcnt][1] = node_remap[face2node[fc][1]] ;
          quad_faces[qcnt][2] = node_remap[face2node[fc][2]] ;
          quad_faces[qcnt][3] = node_remap[face2node[fc][3]] ;
          qcnt++ ;
        }
      }

      for(int i=0;i<lower[cc].size();++i) {
        int fc = lower[cc][i] ;
        swapface[nf] = 1 ;
        faces[nf] = fc ;
        nf++ ;
        int fsz = face2node[fc].size() ;
        if(fsz == 3) {
          tri_faces[tcnt][0] = node_remap[face2node[fc][2]] ;
          tri_faces[tcnt][1] = node_remap[face2node[fc][1]] ;
          tri_faces[tcnt][2] = node_remap[face2node[fc][0]] ;
          tcnt++ ;
        }
        if(fsz == 4) {
          quad_faces[qcnt][0] = node_remap[face2node[fc][3]] ;
          quad_faces[qcnt][1] = node_remap[face2node[fc][2]] ;
          quad_faces[qcnt][2] = node_remap[face2node[fc][1]] ;
          quad_faces[qcnt][3] = node_remap[face2node[fc][0]] ;
          qcnt++ ;
        }
      }

      for(int i=0;i<boundary_map[cc].size();++i) {
        int fc = boundary_map[cc][i] ;
        swapface[nf] = 0 ;
        faces[nf] = fc ;
        nf++ ;
        int fsz = face2node[fc].size() ;
        if(fsz == 3) {
          tri_faces[tcnt][0] = node_remap[face2node[fc][0]] ;
          tri_faces[tcnt][1] = node_remap[face2node[fc][1]] ;
          tri_faces[tcnt][2] = node_remap[face2node[fc][2]] ;
          tcnt++ ;
        }
        if(fsz == 4) {
          quad_faces[qcnt][0] = node_remap[face2node[fc][0]] ;
          quad_faces[qcnt][1] = node_remap[face2node[fc][1]] ;
          quad_faces[qcnt][2] = node_remap[face2node[fc][2]] ;
          quad_faces[qcnt][3] = node_remap[face2node[fc][3]] ;
          qcnt++ ;
        }
      }

      switch(elem_type[cc]) {
      case 0:
        tets_ids[tet_no] = cell_remap[cc];
        fillTet(tets[tet_no++],tri_faces) ;
        break ;
      case 1:
        hexs_ids[hex_no] = cell_remap[cc];
        fillHex(hexs[hex_no++],quad_faces) ;
        break ;
      case 2:
        prsm_ids[prism_no] = cell_remap[cc];
        fillPrism(prsm[prism_no++],tri_faces,quad_faces) ;
        break ;
      case 3:
        pyrm_ids[pyramid_no] = cell_remap[cc];
        fillPyramid(pyrm[pyramid_no++],tri_faces,quad_faces) ;
        break ;
      default:
        generalCellNfaces.push_back(nfaces) ;
        generalCell_ids.push_back(cell_remap[cc]);

        for(int i =0;i<nfaces;++i) {
          int fc = faces[i] ;
          int fsz = face2node[fc].size() ;
          generalCellNsides.push_back(fsz) ;
          if(swapface[i] == 1) {
            for(int j=0;j<fsz;++j) {
              generalCellNodes.push_back(node_remap[face2node[fc][fsz-j-1]]) ;
            }
          } else {
            for(int j=0;j<fsz;++j) {
              generalCellNodes.push_back(node_remap[face2node[fc][j]]) ;
            }
          }
        }
      }
    } ENDFORALL ;


    // write grid topology file
    hid_t file_id = 0, group_id = 0 ;
    file_id=writeVOGOpen(filename) ;
    if(use_parallel_io ||MPI_rank == 0 ) {
      group_id = H5Gcreate(file_id, "elements", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT) ;
    }


    writeUnorderedVector(group_id, "tetrahedra",tets) ;
    writeUnorderedVector(group_id, "tetrahedra_ids",tets_ids) ;

    writeUnorderedVector(group_id, "hexahedra",hexs) ;
    writeUnorderedVector(group_id, "hexahedra_ids",hexs_ids) ;

    writeUnorderedVector(group_id, "prism",prsm) ;
    writeUnorderedVector(group_id, "prism_ids",prsm_ids) ;

    writeUnorderedVector(group_id, "pyramid",pyrm) ;
    writeUnorderedVector(group_id, "pyramid_ids",pyrm_ids) ;

    writeUnorderedVector(group_id, "GeneralCellNfaces",generalCellNfaces) ;
    writeUnorderedVector(group_id, "GeneralCellNsides",generalCellNsides) ;
    writeUnorderedVector(group_id, "GeneralCellNodes", generalCellNodes) ;
    writeUnorderedVector(group_id, "GeneralCell_ids", generalCell_ids) ;


    if(use_parallel_io || MPI_rank == 0) {
      H5Gclose(group_id) ;
      group_id = H5Gcreate(file_id, "boundaries", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT) ;
    }

    const_store<string> boundary_names(bnamesRep) ;

    entitySet boundaries = boundary_names.domain() ;
    if(MPI_processes > 1) {
      entitySet local_boundaries ;
      int bkeyspace = bnamesRep->getDomainKeySpace() ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn(bkeyspace) ;
      Map l2g ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      l2g = df->l2g.Rep() ;
      FORALL(boundaries,bb) {
        if(init_ptn[MPI_rank].inSet(l2g[bb])) {
          local_boundaries += bb ;
        }
      } ENDFORALL ;
      boundaries = local_boundaries ;
    }
    string bnames ;
    FORALL(boundaries,bb) {
      bnames += '"' + boundary_names[bb] + '"' ;
    } ENDFORALL ;


    bnames = MPIConcatStrings(bnames,MPI_COMM_WORLD) ;

    vector<string> bnamelist ;
    for(size_t i=0;i<bnames.size();++i) {
      if(bnames[i] != '"') {
        cerr << "confused in boundary name extraction" << endl ;
        break ;
      }
      ++i ;
      string tmp ;
      while(i < bnames.size() && bnames[i] != '"') {
        tmp += bnames[i++] ;
      }
      bnamelist.push_back(tmp) ;
    }


    // Identify Boundary faces to be written by this processor
    entitySet fset = (MapRepP(boundary_mapRep)->image(localCells)+
                      MapRepP(upperRep)->image(localCells)) & ref.domain() ;

    Map l2f ;
    if(MPI_processes > 1) {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      int kd =  getKeyDomain(fset, df, MPI_COMM_WORLD) ;

      if(kd < 0) {
        cerr << "gridTopology, boundary faces not in single keyspace!"
             << endl ;
        kd = 0 ;
      }
      //debugout << "in write gridTopology, face key domain is " << kd
      //         << endl ;
      l2f = df->l2f.Rep() ;
    } else {
      l2f.allocate(fset) ;
      FORALL(fset,fc) {
	      l2f[fc] = fc ;
      } ENDFORALL ;
    }

    for(size_t i=0;i<bnamelist.size();++i) {
      hid_t bc_id = 0 ;
      string current_bc = bnamelist[i] ;
      if(use_parallel_io || MPI_rank==0) {
        bc_id = H5Gcreate(group_id, current_bc.c_str(), H5P_DEFAULT,
                          H5P_DEFAULT,H5P_DEFAULT) ;
      }

      bool found_ref = false ;
      Entity id = 0 ;
      FORALL(boundary_names.domain(),bb) {
        if(boundary_names[bb] == current_bc) {
          found_ref = true ;
          id = bb ;
        }
      } ENDFORALL ;
      entitySet bfaces ;
      if(found_ref) {
        FORALL(fset,fc) {
          if(ref[fc] == id) {
            bfaces+= fc ;
          }
        }ENDFORALL ;
      }

      int ntria=0, nquad=0, nsided =0;
      FORALL(bfaces,fc) {
        if(face2node[fc].size() == 3) {
          ntria++ ;
        } else if(face2node[fc].size() == 4) {
          nquad++ ;
        } else {
          nsided++ ;
        }
      } ENDFORALL ;
      vector<Array<int,3> > Trias(ntria) ;
      vector<Array<int,4> > Quads(nquad) ;
      vector<int> tria_ids(ntria) ;
      vector<int> quad_ids(nquad) ;
      vector<int> genc_ids(nsided) ;
      int nt = 0 ;
      int nq = 0 ;
      int ng = 0 ;

      vector<int> nsizes(nsided) ;
      vector<int> nsidenodes ;

      FORALL(bfaces,fc) {
        if(face2node[fc].size() == 3) {
          Trias[nt][0] = node_remap[face2node[fc][0]] ;
          Trias[nt][1] = node_remap[face2node[fc][1]] ;
          Trias[nt][2] = node_remap[face2node[fc][2]] ;
          tria_ids[nt] = l2f[fc] ;
          nt++ ;
        } else if(face2node[fc].size() == 4) {
          Quads[nq][0] = node_remap[face2node[fc][0]] ;
          Quads[nq][1] = node_remap[face2node[fc][1]] ;
          Quads[nq][2] = node_remap[face2node[fc][2]] ;
          Quads[nq][3] = node_remap[face2node[fc][3]] ;
          quad_ids[nq] = l2f[fc] ;
          nq++ ;
        } else {
          nsizes[ng] = face2node[fc].size() ;
          for(int i=0;i<nsizes[ng];++i) {
            nsidenodes.push_back(node_remap[face2node[fc][i]]) ;
          }
          genc_ids[ng] = l2f[fc] ;
          ng++ ;
        }
      } ENDFORALL ;


      writeUnorderedVector(bc_id,"triangles",Trias) ;
      writeUnorderedVector(bc_id,"triangles_id",tria_ids) ;

      writeUnorderedVector(bc_id,"quads",Quads) ;
      writeUnorderedVector(bc_id,"quads_id",quad_ids) ;

      writeUnorderedVector(bc_id,"nside_sizes",nsizes) ;
      writeUnorderedVector(bc_id,"nside_nodes",nsidenodes) ;
      writeUnorderedVector(bc_id,"nside_id",genc_ids) ;

      if(use_parallel_io || MPI_rank == 0) {
        H5Gclose(bc_id) ;
      }
    }

    if(use_parallel_io || MPI_rank == 0) {
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
    }

  }

  /// Returns all boundary names.
  /// @param[in] bnamesRep The boundary name store.
  /// @param[in,out] facts The fact database.
  /// @return A string vector of all boundary names.
  std::vector<string> get_boundary_names(storeRepP bnamesRep, fact_db &facts){

    const_store<string> boundary_names(bnamesRep) ;
    entitySet boundaries = boundary_names.domain() ;
    if(MPI_processes > 1) {
      entitySet local_boundaries ;
      int bkeyspace = bnamesRep->getDomainKeySpace() ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn(bkeyspace) ;
      Map l2g ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      l2g = df->l2g.Rep() ;
      FORALL(boundaries,bb) {
        if(init_ptn[MPI_rank].inSet(l2g[bb])) {
          local_boundaries += bb ;
        }
      } ENDFORALL ;
      boundaries = local_boundaries ;
    }

    string bnames ;
    FORALL(boundaries,bb) {
      bnames += '"' + boundary_names[bb] + '"' ;
    } ENDFORALL ;


    bnames = MPIConcatStrings(bnames,MPI_COMM_WORLD) ;

    vector<string> bnamelist ;
    for(size_t i=0;i<bnames.size();++i) {
      if(bnames[i] != '"') {
        cerr << "confused in boundary name extraction" << endl ;
        break ;
      }
      ++i ;
      string tmp ;
      while(i < bnames.size() && bnames[i] != '"') {
        tmp += bnames[i++] ;
      }
      bnamelist.push_back(tmp) ;
    }
    std::sort(bnamelist.begin(), bnamelist.end());
    return bnamelist;
  }

  /// Creates an directory with the name scheme: `/output/$bc_name`.
  /// @param[in] bc_name The name of the boundary.
  void get_bc_directory(string bc_name){
    string directory_name = "output/" + bc_name;
    struct stat statbuf ;
    int fid = open(directory_name.c_str(), O_RDONLY) ;
    if(fid < 0) {
      mkdir(directory_name.c_str(),0755) ;
    } else {
      fstat(fid,&statbuf) ;
      if(!S_ISDIR(statbuf.st_mode)) {
        cerr << "file " << directory_name <<" should be a directory!, rename "
             << "it and start again." << endl ;
        Loci::Abort() ;
      }
      close(fid) ;
    }
  }

  /// Opens an HDF5 file from the location: `/output/$bc_name/$file_name`.
  /// @param[in] bc_name The name of the boundary.
  /// @param[in] file_name The name of the file to open.
  /// @return An HDF5 file ID (hid_t) to the file.
  hid_t open_boundary_file(string bc_name, string file_name) {
    // open the file
    hid_t file_id = 0;
    if(MPI_rank == 0) {
      get_bc_directory(bc_name);
    }
    string dirname = "output/" + bc_name + "/";
    string filename = dirname+file_name;
    file_id=writeVOGOpen(filename.c_str()) ;
    return file_id;
  }

  /// Get boundary faces that belong to a boundary surface `current_bc`.
  /// @param[in] current_bc Boundary name.
  /// @param[in] refRep The `ref` map.
  /// @param[in] bnamesRep The boundary name store.
  /// @param[in] fset All boundary faces.
  /// @return An entitySet of all boundary faces that belong to the boundary
  ///         surface `current_bc` provided to this function.
  entitySet get_boundary_faces(string current_bc, storeRepP refRep,
                               storeRepP bnamesRep, entitySet fset) {

    const_store<string> boundary_names(bnamesRep) ;
    const_Map ref(refRep) ;

    //find its ref id
    bool found_ref = false ;
    Entity id = 0 ;
    FORALL(boundary_names.domain(),bb) {
      if(boundary_names[bb] == current_bc) {
        found_ref = true ;
        id = bb ;
      }
    } ENDFORALL ;
    entitySet bfaces ;
    if(found_ref) {
      FORALL(fset,fc) {
        if(ref[fc] == id) {
          bfaces+= fc ;
        }
      }ENDFORALL ;
    }
    return bfaces;
  }

  /// Get boundary nodes that belong to a boundary surface `current_bc`.
  /// @param[in] current_bc Boundary name.
  /// @param[in] face2nodeRep The face to node map.
  /// @param[in] refRep The `ref` map.
  /// @param[in] bnamesRep The boundary name store.
  /// @param[in] fset All boundary faces.
  /// @param[in,out] facts The fact database.
  /// @return The entitySet of all boundary nodes that belong to the boundary
  ///         surface `current_bc` provided to this function.
  entitySet get_boundary_nodes(string current_bc, storeRepP face2nodeRep,
        storeRepP refRep, storeRepP bnamesRep, entitySet fset,
        fact_db &facts ){
    entitySet bfaces = get_boundary_faces(current_bc, refRep, bnamesRep, fset);
    //get containers
    const_multiMap face2node(face2nodeRep) ;
    //get all the local nodes on the boundary
    entitySet nodes_local = MapRepP(face2node)->image(bfaces) ;
    // get the nodes that belong to this processor
    if(MPI_processes > 1) {
      entitySet dom = nodes_local ;
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      entitySet  my_entities = dist->my_entities ;
      nodes_local = my_entities & nodes_local ;
    }
    return nodes_local ;
  }

  /// Write out the boundary surface topology to an HDF5 file.
  /// @param[in] file_id The HDF5 file ID of this boundary surface.
  /// @param[in] face2nodeRep The face to node map.
  /// @param[in] bfaces The boundary faces that define this surface.
  /// @param[in,out] facts The fact database.
  void writeBoundaryTopo(hid_t file_id, storeRepP face2nodeRep,
                         entitySet bfaces, fact_db &facts ) {
#ifdef VERBOSE
    debugout << "writing out boundary surface topology" << endl ;
    stopWatch s ;
    s.start() ;
#endif

    const_multiMap face2node(face2nodeRep) ;

    //collect the local boundary nodes belong to this boundary
    entitySet nodes_local = (MapRepP(face2node)->image(bfaces)) ;

    //map each node to its file number in output file
    Map node_remap ;
    node_remap = get_output_node_remap(facts, nodes_local) ;

#ifdef VERBOSE
    debugout << "time to get_output_node_remap = " << s.stop() << endl ;
    s.start() ;
#endif
    // Compute face reordering for topo sorting
    store<int> faceorder ;
    faceorder.allocate(bfaces) ;
    int sz = bfaces.size() ;
    int off = 0 ;
    MPI_Scan(&sz,&off,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    off -= sz ;
    FORALL(bfaces,fc) {
      faceorder[fc] = off ;
      off++ ;
    } ENDFORALL ;

    //get output vectors
    int ntria=0, nquad=0, nsided =0 ;
    FORALL(bfaces,fc) {
      if(face2node[fc].size() == 3) {
        ntria++ ;
      } else if(face2node[fc].size() == 4) {
        nquad++ ;
      } else {
        nsided++ ;
      }
    } ENDFORALL ;
    vector<Array<int,3> > Trias(ntria) ;
    vector<Array<int,4> > Quads(nquad) ;
    vector<int> tria_ids(ntria) ;
    vector<int> quad_ids(nquad) ;
    vector<int> genc_ids(nsided) ;
    int nt = 0 ;
    int nq = 0 ;
    int ng = 0 ;

    vector<int> nsizes(nsided) ;
    vector<int> nsidenodes ;

    FORALL(bfaces,fc) {
      if(face2node[fc].size() == 3) {
        Trias[nt][0] = node_remap[face2node[fc][0]] ;
        Trias[nt][1] = node_remap[face2node[fc][1]] ;
        Trias[nt][2] = node_remap[face2node[fc][2]] ;
        tria_ids[nt] = faceorder[fc] ;
        nt++ ;
      } else if(face2node[fc].size() == 4) {
        Quads[nq][0] = node_remap[face2node[fc][0]] ;
        Quads[nq][1] = node_remap[face2node[fc][1]] ;
        Quads[nq][2] = node_remap[face2node[fc][2]] ;
        Quads[nq][3] = node_remap[face2node[fc][3]] ;

        quad_ids[nq] = faceorder[fc] ;
        nq++ ;
      } else {
        nsizes[ng] = face2node[fc].size() ;
        for(int i=0;i<nsizes[ng];++i) {
          nsidenodes.push_back(node_remap[face2node[fc][i]]) ;
        }

        genc_ids[ng] = faceorder[fc] ;
        ng++ ;
      }
    } ENDFORALL ;

#ifdef VERBOSE
    debugout << "time to generate topology datasets = " << s.stop() << endl ;
    s.start() ;
#endif

    //write out vectors
    writeUnorderedVector(file_id,"triangles",Trias) ;
    writeUnorderedVector(file_id,"triangles_ord",tria_ids) ;
    writeUnorderedVector(file_id,"quads",Quads) ;
    writeUnorderedVector(file_id,"quads_ord",quad_ids) ;

    writeUnorderedVector(file_id,"nside_sizes",nsizes) ;
    writeUnorderedVector(file_id,"nside_nodes",nsidenodes) ;
    writeUnorderedVector(file_id,"nside_ord",genc_ids) ;

#ifdef VERBOSE
    debugout << "time to write unordered vectors = " << s.stop() << endl ;
    debugout << "finished writing boundary topology" << endl ;
#endif
  }

  /// Finds the index of an inner edge for a given `face`, and creates one
  /// if it does not exist.
  ///
  /// An inner edge is the edge from facecenter to one of its nodes during
  /// face triangulation. If the edge already exists, the index is returned,
  /// otherwise the edge is created and registered in `inner_edges` and
  /// `edgeIndex`.
  /// @param[in] face The face ID.
  /// @param[in] node The node ID.
  /// @param[in,out] edgeIndex A map of inner edges to its index in
  ///                          `inner_edges`. The index starts with 1.
  /// @param[in,out] inner_edges The edges from facecenter and one of its
  ///                            nodes.
  /// @return The index of the edge in `inner_edges`, starts with 1.
  int create_inner_edge(int face, int node,
                       std::map<pair<int, int>, int>& edgeIndex,
                       vector<pair<int, int> >& inner_edges) {
    pair<int, int> inner_edge = make_pair(face, node) ;
    //make sure no duplicate in inner_edges
    if(edgeIndex.find(inner_edge) == edgeIndex.end()){
      inner_edges.push_back(inner_edge) ;
      edgeIndex[inner_edge] = inner_edges.size() ;
    }
    return edgeIndex[inner_edge] ;
  }

  /// If a face has more than 2 edges cut, to disambiguate the connection
  /// between the cut points, the face will be triangulated and each tri face
  /// will be registered.
  ///
  /// @param f The face entity.
  /// @param face2node The face to node map.
  /// @param face2edge The face to edge map.
  /// @param edge2node The edge to node map.
  /// @param val .
  /// @param edgesCut The edges (node to node) that are cut.
  /// @param edgeIndex A map of inner edges (face to node) to its index in
  ///                  `inner_edges`, the index starts with 1.
  /// @param intersects The pair of edges, if the second number is negative,
  ///                   it is the index to `inner_edges`.
  /// @param inner_edges The edges from facecenter and one of its nodes.
  void disambiguateFace(Entity f, const_multiMap& face2node,
        const_multiMap& face2edge, const_MapVec<2>& edge2node,
        const_store<double>& val, entitySet& edgesCut,
        std::map<pair<int, int>, int >& edgeIndex,
        vector<pair<int, int> >& intersects,
        vector<pair<int, int> >& inner_edges ) {

    int nNodes = face2node.num_elems(f) ;

    //get the value of face center
    double  newNode = 0.0 ;
    for (int i = 0; i < nNodes; ++i) {
      newNode += val[face2node[f][i]] ;
    }
    newNode /= double(nNodes) ;

    //triangulate the face
    //and register each tri face.
    int faceCut[2] ;
    for (int i = 0;i < face2edge.num_elems(f); ++i) {
      int cutsFound = 0 ;
      //check this edge
      Entity edge = face2edge[f][i];
      if (signbit(val[edge2node[edge][0]] * val[edge2node[edge][1]])) {
        edgesCut += edge ;
        faceCut[cutsFound] = edge ;
        cutsFound++ ;
      }
      //check face center to first node
      if (signbit(newNode * val[edge2node[edge][0]])) {
        int node = edge2node[edge][0] ;
        int edgeId = create_inner_edge(f, node, edgeIndex, inner_edges) ;
        faceCut[cutsFound]  = -edgeId ;
        cutsFound++ ;
      }

      //check face center to second node
      if (signbit(newNode * val[edge2node[edge][1]])) {
        if(cutsFound > 1) {
          cerr<<"ERROR: tri face has more than two cuts" << endl ;
          exit(-1) ;
        }
        int node = edge2node[edge][1] ;
        int edgeId = create_inner_edge(f, node, edgeIndex, inner_edges) ;
        faceCut[cutsFound]  = -edgeId ;
        cutsFound++ ;
      }
      if (cutsFound == 2) { intersects.push_back(make_pair(faceCut[0], faceCut[1])) ; }
      //else ????holes happen?
    }
  }

  /// Find edges cut in the face and register them into `intersects`.
  /// @param f The face entity.
  /// @param face2node The face to node map.
  /// @param face2edge The face to edge map.
  /// @param edge2node The edge to node map.
  /// @param val .
  /// @param edgesCut The edges (node 2 node) that are cut.
  /// @param edgeIndex A map of inner edges (face center 2 node) to its index
  ///                  in `inner_edges`, the index starts with 1.
  /// @param intersects The pair of edges, if the second number is negative,
  ///                   it is the index to `inner_edges`.
  /// @param inner_edges A pair<face, node>, the edges from facecenter and one
  ///                    of the face nodes.
  bool registerFace(Entity f, const_multiMap& face2node,
        const_multiMap& face2edge, const_MapVec<2>& edge2node,
        const_store<double>& val, entitySet& edgesCut,
        std::map<pair<int, int>, int >& edgeIndex,
        vector<pair<int, int> >& intersects,
        vector<pair<int, int> >& inner_edges) {

    int faceCut[2] ;
    int cutsFound = 0 ;
    //check each edge
    for(int i = 0; i < face2edge.num_elems(f); ++i) {
      Entity edge = face2edge[f][i] ;
      //if it is cut
      if (signbit(val[edge2node[edge][0]] * val[edge2node[edge][1]])) {
        if (cutsFound < 2) {
          //store it in faceCut
          faceCut[cutsFound] = edge ;
          edgesCut += edge ;
        }
        cutsFound++ ;
      }
    }
    // no ambiguation, register faceCut in inertsection otherwise, disambiguate
    // the face.
    if (cutsFound == 2) {
      intersects.push_back(make_pair(faceCut[0], faceCut[1])) ;
    } else if (cutsFound > 2) {
      disambiguateFace(f, face2node, face2edge, edge2node, val, edgesCut,
                       edgeIndex, intersects, inner_edges) ;
    }
    return (cutsFound > 0) ;
  }

  /// After all the faces of a cell are registered, check the intersects to
  /// form the faces in cutting plane.
  /// @param intersects The pair of edges, if the value number is negative, it
  ///                   is index to inner_edges.
  /// @param faceLoops The loops(faces in cutting plane) formed.
  /// @param[in] start The start index to intersects.
  /// @param[in] end The end index to intersects.
  void checkLoop(const vector<pair<int, int> >& intersects,
                 vector<vector<int > >& faceLoops, int start, int end) {
    bool loopIsGood = true ;
    vector<int > loop ;

    list<int> edges ;
    for (int i = start+1; i < end; i++) {
      edges.push_back(i) ;
    }

    int firstNode, nextNode ;
    firstNode = intersects[start].first ;
    nextNode = intersects[start].second ;
    loop.push_back(firstNode) ;
    loop.push_back(nextNode) ;

    list<int>::iterator iter ;
    while (!edges.empty() && loopIsGood) {
      for (iter = edges.begin(); iter != edges.end(); iter++) {
        if (intersects[*iter].first == nextNode) {
          nextNode = intersects[*iter].second ;
          if(firstNode != nextNode) { loop.push_back(nextNode) ; }
          edges.erase(iter) ;
          break ;
        } else if (intersects[*iter].second == nextNode) {
          nextNode = intersects[*iter].first ;
          if(firstNode != nextNode) {loop.push_back(nextNode) ; }
          edges.erase(iter) ;
          break ;
        }
      }
      //can not find the next edge
      if (iter == edges.end()) {
        loopIsGood = false ;
      }

      //a loop is formed
      if (firstNode == nextNode) {
        if(!edges.empty()) {
          faceLoops.push_back(loop) ;
          loop.clear() ;

          firstNode = intersects[edges.front()].first ;
          nextNode = intersects[edges.front()].second ;
          loop.push_back(firstNode) ;
          loop.push_back(nextNode) ;
          edges.erase(edges.begin()) ;
        } else {
          faceLoops.push_back(loop) ;
        }
      }
    }

    if (firstNode != nextNode) {
      debugout << "ERROR: ** Problem cell:  ** (failed loop test)" << endl ;
    }
  }

  /// Removes all negative values in loops after the loops are formed.
  /// `faceLoops` is the loops(faces in cutting plane) formed.
  /// @param[in] faceLoops Loops (faces in cutting plane) formed.
  /// @return Updated set of `faceLoops` with all negative values in loops
  ///         removed.
  vector<vector<int > > remove_inner_edges(
        const vector<vector<int > >& faceLoops) {
    vector<vector<int > > new_faceLoops ;
    for(unsigned int i=0; i<faceLoops.size(); i++){
      vector<int> loop ;
      for(unsigned int j=0; j<faceLoops[i].size(); j++){
        if(faceLoops[i][j] >=0) { loop.push_back(faceLoops[i][j]) ; }
      }
      if(loop.size() >0) { new_faceLoops.push_back(loop) ; }
    }
    return new_faceLoops ;
  }

  /// Returns the cutting position of an edge.
  double get_edge_weight(Entity e, const_MapVec<2>& edge2node,
                         const_store<double>& val) {

    double a = val[edge2node[e][0]] ;
    double b = val[edge2node[e][1]] ;
    return ((b)/(b - a + 1e-30)) ;
  }


  /// Generates Cutplane.
  CutPlane getCutPlane(storeRepP upperRep, storeRepP lowerRep,
                       storeRepP boundary_mapRep, storeRepP face2nodeRep,
                       storeRepP face2edgeRep, storeRepP edge2nodeRep,
                       storeRepP valRep, entitySet localCells,//all geom_cells
                       fact_db &facts) {

    const_multiMap upper(upperRep), lower(lowerRep),
      boundary_map(boundary_mapRep), face2node(face2nodeRep),
      face2edge(face2edgeRep) ;

    const_MapVec<2> edge2node(edge2nodeRep);
    const_store<double> val(valRep) ;

    // data structure:
    entitySet edgesCut ;

    // the weight for interpoplation for each edgesCut, allocated on edgesCut
    store<double> edgesWeight ;

    // the inner edges(facecenter to one of its nodes) cut, the values stored
    // are pair<local_faceid, noderank>
    vector<pair<int, int> > inner_edges ;

    // loops formed, the values stored are edge ids, which is either local edge
    // entity or index to inner_edges
    vector<vector<int > > faceLoops ;


    //extra data structure
    //the pairs of edges  or inner edges cut
    vector<pair<int, int> > intersects ;

    //map of inner edges to their indexes in inner_edges
    std::map<pair<int, int>, int > edgeIndex ;
    // for each cell
    FORALL(localCells,cc) {
      bool isCut = false ;
      int  intersectStart = intersects.size() ;

      //collect all its faces
      entitySet faces;
      for(int i=0;i<upper[cc].size();++i) {faces += upper[cc][i] ; }
      for(int i=0;i<lower[cc].size();++i) {faces += lower[cc][i] ; }
      for(int i=0;i<boundary_map[cc].size();++i) {
        faces += boundary_map[cc][i] ;
      }

      //for each face, find the edges that is cut and register it
      for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ei++) {
        if (registerFace(*ei, face2node, face2edge, edge2node, val, edgesCut,
                         edgeIndex, intersects, inner_edges)) {
          isCut = true ;
        }
      }

      //check the loops formed by this cell and register the loop
      if(isCut) {
        checkLoop(intersects, faceLoops, intersectStart, intersects.size()) ;
        intersectStart = intersects.size() ;
      }
    }ENDFORALL;
    //remove the inner edges cut from loops
    vector<vector<int > > new_faceLoops = remove_inner_edges(faceLoops) ;

    //compute the cutting postions of edges
    edgesWeight.allocate(edgesCut) ;
    FORALL(edgesCut, e){
      double t = get_edge_weight(e, edge2node, val) ;
      edgesWeight[e] = t ;
    }ENDFORALL ;

    return CutPlane(edgesWeight.Rep(), new_faceLoops) ;
  }

  /// Write out the cut plane topology to an HDF5 file.
  /// @param[in] bc_id The HDF5 file ID of this cut plane.
  /// @param[in] cp The cut plane data structure.
  /// @param[in] facts The fact database.
  void writeCutPlaneTopo(hid_t bc_id, const CutPlane& cp, fact_db &facts) {
#ifdef VERBOSE
    debugout << "write cutPlaneTopology" << endl ;
    stopWatch s ;
    s.start() ;
#endif
    //map from local numbering to output file numbering for edge entities
    Map node_remap ;
    entitySet edgesCut = (cp.edgesWeight)->domain() ;
    node_remap = get_output_node_remap(facts, edgesCut) ;
#ifdef VERBOSE
    debugout << "time to get remap = " << s.stop() << endl ;
    s.start() ;
#endif

    //write out the sizes of faceLoops
    int num_faces = cp.faceLoops.size() ;
    vector<int> nsizes(num_faces) ;
    for(int i = 0; i < num_faces; i++) {
      nsizes[i] = cp.faceLoops[i].size() ;
    }

    writeUnorderedVector(bc_id,"nside_sizes",nsizes) ;

    //write out face nodes
    vector<int> nsidenodes ;
    for(int i = 0; i < num_faces; i++) {
      for(int j = 0; j < nsizes[i]; j++) {
        int node = cp.faceLoops[i][j] ;
        if(node >=0) {
          nsidenodes.push_back(node_remap[node]) ; //edge nodes
        }
      }
    }

    writeUnorderedVector(bc_id,"nside_nodes",nsidenodes) ;

#ifdef VERBOSE
    debugout << "time to write cut plane topology=" << s.stop() << endl ;
#endif
  }

  namespace {

    /// Struct to hold boundary condition information. Used to access options
    /// that are given as suboptions to boundary conditions. For something like
    /// BC5=exampleBC(option1=value1,option2=value2), this struct holds the
    /// metadata about a BC and the options that have been provided to it.
    struct BCinfo {
      std::string name ; // Boundary condition type name
      int key ; // Boundary ID
      entitySet apply_set ; // Set of entities to which this BC applies
      options_list bc_options ; // Options for this boundary condition
      BCinfo() {}
      BCinfo(const std::string &n, int k, const entitySet &a,
             const options_list &o) :
            name(n),key(k),apply_set(a),bc_options(o) {}

    } ;
  }

  /// Builds facts related to the periodic boundary mappings as well as the
  /// `pmap` and `periodicTransform` facts.
  ///
  /// For each periodic pair, compute face centers, apply the rigid transform
  /// to the master side, match faces by nearest neighbor, verify one-to-one
  /// correspondence, then store:
  /// - `pmap`: face -> periodic partner face
  /// - `periodicTransform`: bc_id -> rigid_transform
  ///
  /// @param[in] periodic_list List of (master, slave) periodic_info pairs with
  ///                          bc_num, bset, and transform parameters.
  /// @param[in,out] facts Fact database.
  void setup_periodic_bc(list<pair<periodic_info,periodic_info> >
                         &periodic_list,fact_db &facts) {

    dMap pmap ;
    dstore<rigid_transform> periodic_transform ;


    // Compute fluid face centers
    store<vector3d<double> > pos ;
    pos = facts.get_variable("pos") ;

    // First fill in tmp_pos so that it is valid for any reference to
    // it from face2node on this processor.
    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    int fk = face2node.Rep()->getDomainKeySpace() ;
    entitySet f2n_image = MapRepP(face2node.Rep())->image(face2node.domain()) ;
    entitySet out_of_dom = f2n_image - pos.domain() ;
    dstore<vector3d<double> > tmp_pos ;
    FORALL(pos.domain(), pi) {
      tmp_pos[pi] = pos[pi] ;
    } ENDFORALL ;
    storeRepP sp = tmp_pos.Rep() ;
    int tmp_out = out_of_dom.size() ;
    std::vector<entitySet> init_ptn ;
    if(facts.is_distributed_start()) {
      int pk = pos.Rep()->getDomainKeySpace() ;
      init_ptn = facts.get_init_ptn(pk) ;
      if(GLOBAL_OR(tmp_out)) {
        fill_clone(sp, out_of_dom, init_ptn) ;
      }
    }

    list<pair<periodic_info,periodic_info> >::const_iterator ii ;

    for(ii=periodic_list.begin();ii!=periodic_list.end();++ii) {
      int bc1 = ii->first.bc_num ;
      int bc2 = ii->second.bc_num ;
      double angle = realToDouble(ii->first.angle) ;
      vector3d<double> center(realToDouble(ii->first.center.x),
			      realToDouble(ii->first.center.y),
			      realToDouble(ii->first.center.z)) ;
      vector3d<double> v(realToDouble(ii->first.v.x),
			 realToDouble(ii->first.v.y),
			 realToDouble(ii->first.v.z)) ;

      vector3d<double> trans(realToDouble(ii->first.translate.x),
			     realToDouble(ii->first.translate.y),
			     realToDouble(ii->first.translate.z)) ;


      periodic_transform[bc1] = rigid_transform(center,v,angle,trans) ;
      periodic_transform[bc2] = rigid_transform(center,v,-angle,-1.*trans) ;

      // Compute face centers for point matching
      dstore<vector3d<double> > p1center ;
      entitySet p1Set = ii->first.bset ;
      rigid_transform tran = periodic_transform[bc1] ;
      for(entitySet::const_iterator ei = p1Set.begin();ei!=p1Set.end();++ei) {
        vector3d<double> tot = vector3d<double>(0.0,0.0,0.0);
        const int sz = face2node.end(*ei)-face2node.begin(*ei) ;
        for(int i=0; i<sz; ++i) {
          tot += tmp_pos[face2node[*ei][i]] ;
        }
        tot *= double(1)/double(sz) ;
        vector3d<double> totr(tot.x,tot.y,tot.z) ;
        totr = tran.transform(totr) ;
        p1center[*ei] = realToDouble(totr) ;
      }

      dstore<vector3d<double> > p2center ;
      entitySet p2Set = ii->second.bset ;
      for(entitySet::const_iterator ei = p2Set.begin();ei!=p2Set.end();++ei) {
        vector3d<double> tot = vector3d<double>(0.0,0.0,0.0);
        const int sz = face2node.end(*ei)-face2node.begin(*ei) ;
        for(int i=0; i<sz; ++i) {
          tot += tmp_pos[face2node[*ei][i]] ;
        }
        tot *= double(1)/double(sz) ;
        p2center[*ei] = tot ;
      }

      // Find closest points
      vector<kdTree::coord3d> p1(p1center.domain().size()) ;
      vector<int> p1id(p1center.domain().size()) ;
      int cnt = 0 ;
      FORALL(p1center.domain(),fc) {
        p1[cnt][0] = realToDouble(p1center[fc].x) ;
        p1[cnt][1] = realToDouble(p1center[fc].y) ;
        p1[cnt][2] = realToDouble(p1center[fc].z) ;
        p1id[cnt] = fc ;
        cnt++ ;
      } ENDFORALL ;

      vector<kdTree::coord3d> p2(p2center.domain().size()) ;
      vector<int> p2id(p2center.domain().size()) ;
      cnt = 0 ;
      FORALL(p2center.domain(),fc) {
        p2[cnt][0] = realToDouble(p2center[fc].x) ;
        p2[cnt][1] = realToDouble(p2center[fc].y) ;
        p2[cnt][2] = realToDouble(p2center[fc].z) ;
        p2id[cnt] = fc ;
        cnt++ ;
      } ENDFORALL ;

      vector<int> p1closest(p1.size()) ;

      parallelNearestNeighbors(p2,p2id,p1,p1closest,MPI_COMM_WORLD) ;

      vector<int> p2closest(p2.size()) ;
      parallelNearestNeighbors(p1,p1id,p2,p2closest,MPI_COMM_WORLD) ;

      for(size_t i=0;i<p1.size();++i) {
        pmap[p1id[i]] = p1closest[i] ;
      }
      for(size_t i=0;i<p2.size();++i) {
        pmap[p2id[i]] = p2closest[i] ;
      }

      // Check to make sure connection is one-to-one
      dstore<int> check ;
      for(size_t i=0;i<p2.size();++i) {
        check[p2id[i]] = p2closest[i] ;
      }

      entitySet p1map = create_entitySet(p1closest.begin(),p1closest.end()) ;
      storeRepP sp = check.Rep() ;
      std::vector<entitySet> init_ptn ;
      if(facts.is_distributed_start()) {
        init_ptn = facts.get_init_ptn(fk) ;
        fill_clone(sp, p1map, init_ptn) ;
      }
      bool periodic_problem = false ;
      for(size_t i=0;i<p1id.size();++i) {
        if(check[pmap[p1id[i]]] != p1id[i]) {
          periodic_problem = true ;
        }
      }
      if(GLOBAL_OR(periodic_problem)) {
        if(MPI_rank == 0) {
          cerr << "Periodic boundary did not connect properly, is boundary "
               << "point matched?" << endl ;
          Loci::Abort() ;
        }
      }
    }

    // Add periodic datastructures to fact database
    pmap.Rep()->setDomainKeySpace(fk) ;
    MapRepP(pmap.Rep())->setRangeKeySpace(fk) ;
    facts.create_fact("pmap",pmap) ;
    facts.create_fact("periodicTransform",periodic_transform) ;

    constraint pfaces ;
    Map cl ;
    pfaces = facts.get_variable("periodicFaces") ;
  }

  /// Creates the `ci` map for boundary faces that are not periodic.
  /// @param[in,out] facts The fact database.
  void create_ci_map(fact_db &facts) {
    constraint boundary_faces ;
    boundary_faces = facts.get_variable("boundary_faces") ;
    entitySet ci_faces = *boundary_faces ;

    storeRepP pfacesP = facts.get_variable("periodicFaces") ;
    if(pfacesP != 0) {
      constraint periodicFaces ;
      periodicFaces = pfacesP ;
      debugout << "periodicFaces = " << periodicFaces << endl ;
      ci_faces -= *periodicFaces ;
    }

    Map cl,ci ;

    cl = facts.get_variable("cl") ;
    ci.allocate(ci_faces) ;

    FORALL(ci_faces,fc) {
      ci[fc] = cl[fc] ;
    } ENDFORALL ;
    ci.Rep()->setDomainKeySpace(cl.Rep()->getDomainKeySpace()) ;
    facts.create_fact("ci",ci) ;
    debugout << "boundary_faces = " << *boundary_faces << endl ;
    debugout << "ci_faces = " << ci_faces << endl ;
  }

  /// Create facts for various types of boundary conditions.
  ///
  /// Creates the `periodicFaces` fact. Creates the `<BCName>_BC` fact for
  /// each boundary, where <BCName> is the name of each boundary condition.
  /// Also creates `notPeriodicCells`, `no_symmetry_BC`, and `BC_options` facts.
  /// @param[in,out] facts The fact database.
  void setupBoundaryConditions(fact_db &facts) {
    list<BCinfo> BCinfo_list ;
    std::map<std::string,entitySet> BCsets ;

    /*Boundary Conditions*/
    entitySet periodic ;
    constraint periodic_faces;
    constraint no_symmetry_BC ;

    entitySet symmetry ;

    storeRepP tmp = facts.get_variable("boundary_names") ;
    if(tmp == 0) {
      throw(StringError("boundary_names not found in setupBoundaryConditions! Grid file read?")) ;
    }
    store<string> boundary_names ;
    store<string> boundary_tags ;
    boundary_names = tmp ;
    boundary_tags = facts.get_variable("boundary_tags") ;

    Map ref ;
    ref = facts.get_variable("ref") ;
    entitySet dom = boundary_names.domain() ;
    dom = all_collect_entitySet(dom) ;// gather the boundary surface ids

    int fk = ref.Rep()->getDomainKeySpace() ;

    param<options_list> bc_info ;
    tmp = facts.get_variable("boundary_conditions") ;
    if(tmp == 0) {
      throw(StringError("boundary_conditions not found in setupBoundaryConditions! Is vars file read?")) ;
    }
    bc_info = tmp ;

    param<double> Lref ;
    *Lref = 1.0 ;
    storeRepP p = facts.get_variable("Lref") ;
    if(p != 0) {
      Lref = p ;
    }

    vector<periodic_info> periodic_data ;

    bool fail = false ;
    // WARNING, boundaryName assignment assuming that boundary
    // names and tags have been duplicated on all processors (which
    // they have).  But could break if the grid reader changes to a
    // more standard implementation.
    for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei) {
      Entity bc = *ei ;
      entitySet bcset = interval(bc,bc) ;
      entitySet bfaces = ref.preimage(bcset).first ;

      string bname = boundary_names[bc] ;
      string tname = boundary_tags[bc] ;

      constraint bconstraint ;
      *bconstraint = bfaces ;
      bconstraint.Rep()->setDomainKeySpace(fk) ;

      facts.create_fact(tname,bconstraint) ;
      debugout << "boundary " << bname << "("<< tname << ") = "
               << *bconstraint << endl ;
      param<string> boundarySet ;
      *boundarySet = bname ;
      string factname = "boundaryName(" + bname + ")" ;
      boundarySet.set_entitySet(bfaces) ;
      boundarySet.Rep()->setDomainKeySpace(fk) ;
      facts.create_fact(factname,boundarySet) ;

      option_value_type vt = bc_info->getOptionValueType(bname) ;
      option_values ov = bc_info->getOption(bname) ;
      options_list::arg_list value_list ;
      string name ;

      switch(vt) {
      case NAME :
        ov.get_value(name) ;
        bc_info->setOption(bname,name) ;
        {
          BCinfo_list.push_back(BCinfo(name,bc,bfaces,options_list())) ;
          BCsets[name] += bfaces ;
        }
        break ;
      case FUNCTION:
        ov.get_value(name) ;
        ov.get_value(value_list) ;
        bc_info->setOption(bname,name,value_list) ;
        {
          options_list ol ;
          ol.Input(value_list) ;
          BCinfo_list.push_back(BCinfo(name,bc,bfaces,ol)) ;
          BCsets[name] += bfaces ;
        }
        break ;
      default:
        cerr << "Boundary tag '" << bname << "' exists in VOG file without "
             << "boundary_condition entry!  Please add boundary condition "
             << "specification for this boundary face!" << endl ;
        fail = true ;
      }
      if(name == "symmetry") {
        symmetry += bfaces ;
      } else if(name == "reflecting") {
        symmetry += bfaces ;
      } else if(name == "periodic") {
        periodic += bfaces ;
        periodic_info pi ;
        pi.bc_num = bc ;
        pi.bset = bfaces ;
        options_list ol ;
        ol.Input(value_list) ;
        if(ol.optionExists("rotate") || ol.optionExists("translate")) {
          pi.master = true ;
        }
        if(ol.optionExists("name")) {
          ol.getOption("name",pi.name) ;
        }
        if(ol.optionExists("center")) {
          //get_vect3dOption(ol,"center","m",pi.center,*Lref) ;
          double Lr = realToDouble(*Lref) ;
          ol.getOptionUnits("center","m",pi.center,Lr) ;
        }
        if(ol.optionExists("vector")) {
          //get_vect3d(ol,"vector",pi.v) ;
          ol.getOptionUnits("vector","",pi.v) ;
          pi.v /=  norm(pi.v) ;
        }
        if(ol.optionExists("translate")) {
          //get_vect3dOption(ol,"translate","m",pi.translate,*Lref) ;
          double Lr = realToDouble(*Lref) ;
          ol.getOptionUnits("translate","m",pi.translate,Lr) ;
        }
        if(ol.optionExists("rotate")) {
          ol.getOptionUnits("rotate","radians",pi.angle) ;
        }

        periodic_data.push_back(pi) ;
      }
    }
    if(fail) {
      Loci::Abort() ;
    }


    {
      // Create constraints for each
      std::map<std::string, entitySet>::const_iterator mi ;
      for(mi=BCsets.begin();mi!=BCsets.end();++mi) {
        if(GLOBAL_OR(mi->second.size())) {
          constraint bc_constraint ;
          bc_constraint = mi->second ;
          std::string constraint_name = mi->first + std::string("_BC") ;
	        bc_constraint.Rep()->setDomainKeySpace(fk) ;
          facts.create_fact(constraint_name,bc_constraint) ;
          if(MPI_processes == 1) {
            std::cout << constraint_name << ' ' << mi->second << endl ;
          } else if(MPI_rank == 0) {
            std::cout << "setting boundary condition " << constraint_name
                      << endl ;
          }
        }
      }
    }

    constraint cells ;
    cells = facts.get_variable("cells") ;
    store<options_list> BC_options ;
    BC_options.allocate(dom) ;

    std::map<std::string,entitySet> BC_options_args ;

    for(list<BCinfo>::iterator
          li = BCinfo_list.begin() ; li != BCinfo_list.end() ; ++li) {
      BC_options[li->key] = li->bc_options ;
      options_list::option_namelist onl=li->bc_options.getOptionNameList() ;
      options_list::option_namelist::const_iterator onli  ;

      for(onli=onl.begin();onli!=onl.end();++onli) {
        BC_options_args[*onli] += li->key ;
      }
    }

    {
      std::map<std::string, entitySet>::const_iterator mi ;
      for(mi=BC_options_args.begin();mi!=BC_options_args.end();++mi) {
        constraint bc_constraint ;
        bc_constraint = mi->second ;
        std::string constraint_name = mi->first + std::string("_BCoption") ;
        facts.create_fact(constraint_name,bc_constraint) ;
      }
    }

    if(periodic_data.size() != 0) {
      periodic_faces = periodic ;
      periodic_faces.Rep()->setDomainKeySpace(fk) ;
      facts.create_fact("periodicFaces",periodic_faces) ;

      list<pair<periodic_info,periodic_info> > periodic_list ;
      for(size_t i=0;i<periodic_data.size();++i) {
        if(!periodic_data[i].processed) {
          periodic_data[i].processed = true ;
          periodic_info p1 = periodic_data[i] ;
          periodic_info p2 ;
          p2.name = "," ;
          for(size_t j=i+1;j<periodic_data.size();++j) {
            if(periodic_data[i].name == periodic_data[j].name) {
              if(p2.name != ",") {
                cerr << "periodic name appears more than two times!" ;
                Abort() ;
              }
              p2 = periodic_data[j] ;
              periodic_data[j].processed = true ;
            }
          }
          if(p1.name != p2.name) {
            cerr << "Could not find matching periodic boundary named "
                 << p1.name << endl ;
            Abort() ;
          }
          int p1inp = p1.bset.size() ;
          int p2inp = p2.bset.size() ;
          int p1size ;
          int p2size ;
          MPI_Allreduce(&p1inp,&p1size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
          MPI_Allreduce(&p2inp,&p2size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
          if(p1size != p2size) {
            if(MPI_rank == 0) {
              cerr << "periodic boundaries " << p1.name
                   << " do not match in number of faces" << endl ;
              cerr << "master has " << p1size << " faces and slave has "
                   << p2size << " faces" << endl ;
            }
            Abort() ;
          }
          if((p1.master & p2.master) | (!p1.master & !p2.master)) {
            cerr << "only one master in periodic boundary conditons named "
                 << p1.name << endl ;
            Abort() ;
          }
          if(p1.master) {
            periodic_list.push_back(std::make_pair(p1,p2)) ;
          } else {
            periodic_list.push_back(std::make_pair(p2,p1)) ;
          }
        }
      }

      setup_periodic_bc(periodic_list,facts) ;
    } else {
      constraint notPeriodicCells ;
      *notPeriodicCells = ~EMPTY ;
      facts.create_fact("notPeriodicCells",notPeriodicCells) ;
    }


    entitySet no_symmetry ;
    constraint allfaces ;
    allfaces = facts.get_variable("faces") ;
    no_symmetry  = *allfaces - symmetry ;
    no_symmetry_BC = no_symmetry ;
    no_symmetry_BC.Rep()->setDomainKeySpace(fk) ;
    facts.create_fact("no_symmetry_BC",no_symmetry_BC) ;

    facts.create_fact("BC_options",BC_options) ;

    create_ci_map(facts) ;
  }

  //--------------------------------------------------------------------------
  //-- Create Cell Stencil Plan
  //--
  //-- Get relations from face2cell and face2node
  //-- Join to get node2cell
  //-- Join node2cell & node2cell to get cell2cell (through shared nodes)
  //-- Remove self2self references
  //-- Convert resulting structure to multiMap
  //--------------------------------------------------------------------------

  /// Computes face centers, areas, and normals for all faces.
  /// @param[in] facts The fact database.
  /// @param[out] fcenter Face centers.
  /// @param[out] area Face areas.
  /// @param[out] normal Face normals.
  void getFaceCenter(fact_db &facts, store<vector3d<double> > &fcenter,
                     store<double> &area, store<vector3d<double> > &normal) {
    // Compute face centers
    store<vector3d<double> > pos ;
    pos = facts.get_variable("pos") ;

    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    entitySet fdom = face2node.domain() ;

    // Find nodes that face2node accesses
    entitySet node_access = Loci::MapRepP(face2node.Rep())->image(fdom) ;

    std::vector<vector3d<double> > posdata ;
    int nkeyspace = pos.getDomainKeySpace() ;
    std::vector<entitySet> node_ptn = facts.get_init_ptn(nkeyspace) ;

    // gather nodal data needed for the computation.
    std::map<int,int> g2l ;
    getLocalContextMap(g2l,node_access) ;
    gatherData(posdata,pos,node_access,node_ptn) ;
    fcenter.allocate(fdom) ;
    area.allocate(fdom) ;
    normal.allocate(fdom) ;

    // check centroid type
    param<std::string> centroid ;
    centroid = facts.get_variable("centroid") ;
    bool exact = false ;
    if(*centroid == "exact") {
      exact = true ;
    }
    FORALL(fdom,fc) {
      int fsz = face2node[fc].size() ;
      vector3d<double> center = vector3d<double>(0,0,0) ;
      double wsum = 0 ;

      // loop over face edges to compute a wireframe centroid
      for(int i=1;i<fsz;++i) {
        FATAL(g2l[face2node[fc][i-1]] > int(posdata.size())) ;
        FATAL(g2l[face2node[fc][i-1]] < 0) ;
        FATAL(g2l[face2node[fc][i]] > int(posdata.size())) ;
        FATAL(g2l[face2node[fc][i]] < 0) ;
        double len = norm(posdata[g2l[face2node[fc][i-1]]]-
                          posdata[g2l[face2node[fc][i]]]) ;
        vector3d<double> eloc = 0.5*(posdata[g2l[face2node[fc][i-1]]]+
                                     posdata[g2l[face2node[fc][i]]]) ;
        center += len*eloc ;
        wsum += len ;
      }
      // final edge from last to first node
      double len = norm(posdata[g2l[face2node[fc][0]]]-
                        posdata[g2l[face2node[fc][fsz-1]]]) ;
      vector3d<double> eloc = 0.5*(posdata[g2l[face2node[fc][0]]]+
                                   posdata[g2l[face2node[fc][fsz-1]]]) ;
      center += len*eloc ;
      wsum += len ;
      center *= 1./wsum ;
      fcenter[fc] = center ;
      if(exact) {
        // Iterate to find exact centroid that is on the face for an
        // area-weighted centroid.
        vector3d<double>  tmpcenter = center ;
        const int NITER=4 ;
        for(int iter=0;iter<NITER;++iter) {
          // compute centroid using triangles formed by wireframe centroid
          vector3d<double> centroidsum(0.0,0.0,0.0) ;
          double facearea = 0 ;
          for(int i=0;i<fsz;++i) {
            int n1 = i ;
            int n2 = (i+1==fsz)?0:i+1 ;
            vector3d<double>  p1 = posdata[g2l[face2node[fc][n1]]] ;
            vector3d<double>  p2 = posdata[g2l[face2node[fc][n2]]] ;

            const vector3d<double> t_centroid = (p1 + p2 + tmpcenter)/3.0 ;
            const double t_area = 0.5*norm(cross(p1-tmpcenter,p2-tmpcenter)) ;
            centroidsum += t_area*t_centroid ;
            facearea += t_area ;
          }
          tmpcenter = centroidsum/facearea ;
        }
        fcenter[fc] = tmpcenter ;
      }
      vector3d<double> facearea = vector3d<double>(0,0,0) ;
      vector3d<double> tmpcenter = fcenter[fc] ;
      for(int i=0;i<fsz;++i) {
        int n1 = i ;
        int n2 = (i+1==fsz)?0:i+1 ;
        vector3d<double>  p1 = posdata[g2l[face2node[fc][n1]]] ;
        vector3d<double>  p2 = posdata[g2l[face2node[fc][n2]]] ;

        facearea += cross(p1-tmpcenter,p2-tmpcenter) ;
      }
      double nfacearea = norm(facearea) ;
      area[fc] = 0.5*nfacearea ;
      normal[fc] = (1./nfacearea)*facearea ;
    } ENDFORALL ;
  }

  /// Computes cell centers
  /// @param[in] facts The fact database.
  /// @param[out] ccenter Cell centers.
  /// @param[in] fcenter Face centers.
  /// @param[in] area Face areas.
  void getCellCenter(fact_db &facts, store<vector3d<double> > &ccenter,
                     store<vector3d<double> > &fcenter, store<double> &area) {

    multiMap upper,lower,boundary_map ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;
    entitySet cells = upper.domain()&lower.domain()&boundary_map.domain() ;
    entitySet faceimage  ;
    faceimage += Loci::MapRepP(upper.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(lower.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(boundary_map.Rep())->image(cells) ;
    ccenter.allocate(cells) ;

    vector<vector3d<double> > fcenterdata ;
    vector<double> areadata ;
    int fkeyspace = upper.getRangeKeySpace() ;
    std::vector<entitySet> face_ptn = facts.get_init_ptn(fkeyspace) ;
    std::map<int,int> g2l ;
    getLocalContextMap(g2l,faceimage) ;
    gatherData(fcenterdata,fcenter,faceimage,face_ptn) ;
    gatherData(areadata,area,faceimage,face_ptn) ;
    // compute wireframe centroid
    FORALL(cells,cc) {
      vector3d<double>  csum = vector3d<double> (0,0,0) ;
      double wsum = 0 ;
      int lsz = lower[cc].size() ;
      for(int i=0;i<lsz;++i) {
        int f = g2l[lower[cc][i]] ;
        double w = areadata[f] ;
        vector3d<double>  v = fcenterdata[f] ;
        csum += w*v ;
        wsum += w ;
      }
      int usz = upper[cc].size() ;
      for(int i=0;i<usz;++i) {
        int f = g2l[upper[cc][i]] ;
        double w = areadata[f] ;
        vector3d<double>  v = fcenterdata[f] ;
        csum += w*v ;
        wsum += w ;
      }
      int bsz = boundary_map[cc].size() ;
      for(int i=0;i<bsz;++i) {
        int f = g2l[boundary_map[cc][i]] ;
        double w = areadata[f] ;
        vector3d<double>  v = fcenterdata[f] ;
        csum += w*v ;
        wsum += w ;
      }
      csum *= 1./wsum ;
      ccenter[cc] = csum ;
    } ENDFORALL ;
    param<std::string> centroid ;
    centroid = facts.get_variable("centroid") ;
    if(*centroid == "exact") {
      std::cerr << "Warning: exact cell centroids not implemented yet."
                << std::endl ;
      // Here we add code to compute exact centroid.
      // face2node needs to be expanded to do this..., not trivial
    }
  }

  void get_full_cellStencil(multiMap &cellStencil,fact_db &facts) {
    using std::vector ;
    using std::pair ;
    Map cl,cr ;
    multiMap face2node ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    face2node = facts.get_variable("face2node") ;
    entitySet faces = face2node.domain() ;
    entitySet cellmask = cl.image(faces)+cr.image(faces) ;

    constraint geom_cells_c ;
    geom_cells_c = facts.get_variable("geom_cells") ;
    entitySet geom_cells = *geom_cells_c ;

    // We need to have all of the geom_cells to do the correct test in the
    // loop before, so gather with clone cells
    int gkeyspace = geom_cells_c.getDomainKeySpace() ;
    std::vector<entitySet> ptn = facts.get_init_ptn(gkeyspace) ;
    geom_cells = distribute_entitySet(geom_cells,ptn) ;
    entitySet geom_cell_expand =
      dist_expand_entitySet(geom_cells,cellmask,ptn) ;
    Loci::protoMap f2cell ;

#ifdef DEBUG
    // Check to see if there are any ORPHAN cells in geom_cells.
    entitySet accessedSet = distribute_entitySet(cellmask,ptn) ;
    WARN(GLOBAL_OR((geom_cells-accessedSet) != EMPTY))
#endif
    // Get mapping from face to geometric cells
    Loci::addToProtoMap(cl,f2cell) ;
    FORALL(faces,fc) {
      WARN(cl[fc] == cr[fc]) ;
      if(geom_cell_expand.inSet(cr[fc]))
        f2cell.push_back(pair<int,int>(fc,cr[fc])) ;
    } ENDFORALL ;

    // Get mapping from face to nodes
    Loci::protoMap f2node ;
    Loci::addToProtoMap(face2node,f2node) ;

    // Equijoin on first of pairs to get node to neighboring cell mapping
    // This will give us a mapping from nodes to neighboring cells
    Loci::protoMap n2c ;
    Loci::equiJoinFF(f2node,f2cell,n2c) ;

    // In case there are processors that have no n2c's allocated to them
    // re-balance map distribution
    Loci::balanceDistribution(n2c,MPI_COMM_WORLD) ;
    // Equijoin node2cell with itself to get cell to cell map of
    // all cells that share one or more nodes
    Loci::protoMap n2cc = n2c ;
    Loci::protoMap c2c ;
    Loci::equiJoinFF(n2c,n2cc,c2c) ;

    // Remove self references
    Loci::removeIdentity(c2c) ;

    // Create cell stencil map from protoMap
    distributed_inverseMap(cellStencil,c2c,geom_cells,geom_cells,ptn) ;

#ifdef DEBUG
    entitySet degen_cells ;
    FORALL(geom_cells,ii) {
      if(cellStencil[ii].size() == 0) {
        debugout << "geom cell " << ii << " has degenerate stencil"
                 << endl ;
        degen_cells += ii ;
      }
    } ENDFORALL ;
    entitySet tot_degen = all_collect_entitySet(degen_cells) ;
    WARN(tot_degen!=EMPTY) ;
    if(tot_degen != EMPTY) {
      debugout << "tot_degen=" << tot_degen << endl ;
      if(MPI_rank==0)
        cerr << "tot_degen = " << tot_degen << endl ;
    }
    FORALL(faces,fc) {
      if(tot_degen.inSet(cl[fc]) || tot_degen.inSet(cr[fc])) {
        debugout << "fc=" << fc << ",cl="<<cl[fc] << ",cr=" << cr[fc] << endl ;
        cerr << "fc=" << fc << ",cl="<<cl[fc] << ",cr=" << cr[fc] << endl ;
      }
    } ENDFORALL ;
#endif
  }

  void create_cell_stencil_full(fact_db &facts) {
    // Create cell stencil map from protoMap
    multiMap cellStencil ;
    get_full_cellStencil(cellStencil,facts) ;
    // Put in fact database
    facts.create_fact("cellStencil",cellStencil) ;
  }

  inline void getMaxMinDirections(int &maxid, int &minid,
        const vector3d<double> &dir, std::vector<vector3d<double> > &cdirs) {
    minid = 0 ;
    maxid = 0 ;
    double minval = dot(dir,cdirs[0]) ;
    double maxval = minval ;
    int sz = cdirs.size() ;
    for(int j=1;j<sz;++j) {
      double v = dot(dir,cdirs[j]) ;
      if(minval < v) {
        minid = j ;
        minval = v ;
      }
      if(maxval > v) {
        maxid = j ;
        maxval = v ;
      }
    }
  }

  void get_stable_cellStencil(multiMap &cellStencilFiltered, fact_db & facts) {

    // get full stencil
    multiMap cellStencil ;
    get_full_cellStencil(cellStencil,facts) ;

    entitySet cells = cellStencil.domain() ;
    multiMap upper,lower,boundary_map ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;

    entitySet faceimage  ;
    faceimage += Loci::MapRepP(upper.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(lower.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(boundary_map.Rep())->image(cells) ;
    int fkeyspace = upper.getRangeKeySpace() ;
    std::vector<entitySet> face_ptn = facts.get_init_ptn(fkeyspace) ;
    store<vector3d<double> > fcenter,normal ;
    store<double> area ;
    getFaceCenter(facts,fcenter,area,normal) ;
    store<vector3d<double> > ccenter ;
    getCellCenter(facts,ccenter,fcenter,area) ;

    // gather face data needed for the computation.
    std::map<int,int> fg2l ;
    getLocalContextMap(fg2l,faceimage) ;
    std::vector<vector3d<double> > fcenterdata ;
    std::vector<vector3d<double> > fnormaldata ;
    gatherData(fcenterdata,fcenter,faceimage,face_ptn) ;
    gatherData(fnormaldata,normal,faceimage,face_ptn) ;

    int ckeyspace = upper.getDomainKeySpace() ;
    std::vector<entitySet> cell_ptn = facts.get_init_ptn(ckeyspace) ;

    entitySet cellImage = Loci::MapRepP(cellStencil.Rep())->image(cells) ;

    std::vector<vector3d<double> > ccenterdata ;
    std::map<int,int> cg2l ;
    getLocalContextMap(cg2l,cellImage) ;
    gatherData(ccenterdata,ccenter,cellImage,cell_ptn) ;

    store<int> sizes ;
    sizes.allocate(cells) ;
    FORALL(cells,cc) {
      int csz = cellStencil[cc].size() ;
      int bsz = boundary_map[cc].size() ; ;
      vector3d<double>  ccent = ccenter[cc] ;
      std::vector<vector3d<double> > cdirs(csz+bsz) ;
      for(int i=0;i<csz;++i) {
        cdirs[i] = ccenterdata[cg2l[cellStencil[cc][i]]]-ccent ;
        cdirs[i] *= 1./norm(cdirs[i]) ;
      }
      for(int i=0;i<bsz;++i) {
        cdirs[csz+i] = fcenterdata[fg2l[boundary_map[cc][i]]]-ccent ;
        cdirs[csz+i] *= 1./norm(cdirs[csz+i]) ;
      }

      std::vector<int> flags(csz+bsz,0) ;
      int lsz = lower[cc].size() ;
      for(int i=0;i<lsz;++i) {
        vector3d<double>  dir = fcenterdata[fg2l[lower[cc][i]]] -ccent;
        dir *= 1./norm(dir) ;
        int minid = 0, maxid = 0 ;
        getMaxMinDirections(maxid,minid,dir,cdirs) ;
        flags[minid] = 1 ;
        flags[maxid] = 1 ;
        dir = fnormaldata[fg2l[lower[cc][i]]] ;
        getMaxMinDirections(maxid,minid,dir,cdirs) ;
        flags[minid] = 1 ;
        flags[maxid] = 1 ;
      }
      int usz = upper[cc].size() ;
      for(int i=0;i<usz;++i) {
        vector3d<double>  dir = fcenterdata[fg2l[upper[cc][i]]] -ccent;
        dir *= 1./norm(dir) ;
        int minid = 0, maxid = 0 ;
        getMaxMinDirections(maxid,minid,dir,cdirs) ;
        flags[minid] = 1 ;
        flags[maxid] = 1 ;

        dir = fnormaldata[fg2l[upper[cc][i]]] ;
        getMaxMinDirections(maxid,minid,dir,cdirs) ;
        flags[minid] = 1 ;
        flags[maxid] = 1 ;
      }

      for(int i=0;i<bsz;++i) {
        vector3d<double> dir = fcenterdata[fg2l[boundary_map[cc][i]]] -ccent ;
        dir *= 1./norm(dir) ;
        int minid = 0, maxid = 0 ;
        getMaxMinDirections(maxid,minid,dir,cdirs) ;
        flags[minid] = 1 ;
        flags[maxid] = 1 ;
        minid = 0 ;
        maxid = 0 ;
        dir = fnormaldata[fg2l[boundary_map[cc][i]]] ;
        getMaxMinDirections(maxid,minid,dir,cdirs) ;
        flags[minid] = 1 ;
        flags[maxid] = 1 ;
      }

      int cnt = 0 ;
      for(int i=0;i<csz;++i) {
        if(flags[i] > 0) {
                cellStencil[cc][cnt] = cellStencil[cc][i] ;
          cnt++ ;
        }
      }
      sizes[cc] = cnt ;
    } ENDFORALL ;

    cellStencilFiltered.allocate(sizes) ;
    FORALL(cells,cc) {
      for(int i=0;i<cellStencilFiltered[cc].size();++i) {
	      cellStencilFiltered[cc][i] = cellStencil[cc][i] ;
      }
    } ENDFORALL ;
  }

  void create_cell_stencil(fact_db & facts) {
    // Create cell stencil map from protoMap
    multiMap cellStencil ;
    get_stable_cellStencil(cellStencil,facts) ;
    // Put in fact database
    facts.create_fact("cellStencil",cellStencil) ;
  }

  void get_neigh_cellStencil(multiMap &cellStencil,fact_db &facts) {
    using std::vector ;
    using std::pair ;
    Map cl,cr ;
    multiMap face2node ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    face2node = facts.get_variable("face2node") ;
    entitySet faces = face2node.domain() ;
    entitySet cellmask = cl.image(faces)+cr.image(faces) ;

    constraint geom_cells_c ;
    geom_cells_c = facts.get_variable("geom_cells") ;
    entitySet geom_cells = *geom_cells_c ;

    int gkeyspace = geom_cells_c.getDomainKeySpace() ;
    std::vector<entitySet> ptn = facts.get_init_ptn(gkeyspace) ;
    geom_cells = distribute_entitySet(geom_cells,ptn) ;
    entitySet geom_cell_expand =
      dist_expand_entitySet(geom_cells,cellmask,ptn) ;
    Loci::protoMap f2cell ;

#ifdef DEBUG
    // Check to see if there are any ORPHAN cells in geom_cells.
    entitySet accessedSet = distribute_entitySet(cellmask,ptn) ;
    WARN(GLOBAL_OR((geom_cells-accessedSet) != EMPTY))
#endif

    // Get mapping from face to geometric cells
    Loci::addToProtoMap(cl,f2cell) ;
    FORALL(faces,fc) {
      WARN(cl[fc] == cr[fc]) ;
      if(geom_cell_expand.inSet(cr[fc])) {
        f2cell.push_back(pair<int,int>(fc,cr[fc])) ;
      }
    } ENDFORALL ;

    // join f2cell to get f
    Loci::protoMap f2celll = f2cell ;
    Loci::protoMap f2f ;
    Loci::equiJoinFF(f2cell,f2celll,f2f) ;
    Loci::removeIdentity(f2f) ;
    Loci::balanceDistribution(f2f,MPI_COMM_WORLD) ;

    // Create cell stencil map from protoMap
    distributed_inverseMap(cellStencil,f2f,geom_cells,geom_cells,ptn) ;

  }

  /// Builds a symmetric neighbor stencil for each cell.
  /// @param[out] cellStencilSymm
  /// @param[in] facts
  void get_symm_cellStencil(multiMap &cellStencilSymm, fact_db & facts) {
    if(Loci::MPI_rank==0) {
      cout <<"Generating symm stencil" << endl ;
    }
    using std::vector ;
    // get full stencil
    multiMap cellStencil;
    get_full_cellStencil(cellStencil,facts) ;

    entitySet cells = cellStencil.domain() ;
    multiMap upper,lower,boundary_map ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;

    entitySet faceimage  ;
    faceimage += Loci::MapRepP(upper.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(lower.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(boundary_map.Rep())->image(cells) ;
    int fkeyspace = upper.getRangeKeySpace() ;
    std::vector<entitySet> face_ptn = facts.get_init_ptn(fkeyspace) ;
    store<vector3d<double> > fcenter,normal ;
    store<double> area ;
    getFaceCenter(facts,fcenter,area,normal) ;
    store<vector3d<double> > ccenter ;
    getCellCenter(facts,ccenter,fcenter,area) ;

    // gather face data needed for the computation.
    std::map<int,int> fg2l ; // face global to local
    getLocalContextMap(fg2l,faceimage) ;
    vector<vector3d<double> > fcenterdata, fnormaldata ;
    gatherData(fcenterdata,fcenter,faceimage,face_ptn) ;
    gatherData(fnormaldata,normal,faceimage,face_ptn) ;

    int ckeyspace = upper.getDomainKeySpace() ;
    std::vector<entitySet> cell_ptn = facts.get_init_ptn(ckeyspace) ;

    entitySet cellImage = Loci::MapRepP(cellStencil.Rep())->image(cells) ;

    vector<vector3d<double> > ccenterdata ;
    std::map<int,int> cg2l ;
    getLocalContextMap(cg2l,cellImage) ;
    gatherData(ccenterdata,ccenter,cellImage,cell_ptn) ;

    multiMap neighStencil;
    get_neigh_cellStencil(neighStencil,facts);
    entitySet neighCellImage = Loci::MapRepP(neighStencil.Rep())->image(cells);
    vector<vector3d<double> > ncenterdata ;
    std::map<int,int> ncg2l ;
    getLocalContextMap(ncg2l,neighCellImage) ;
    gatherData(ncenterdata,ccenter,neighCellImage,cell_ptn) ;

    store<int> sizes ;
    sizes.allocate(cells) ;
    vector<int> cellmap ;
    double theta = 3.0*M_PI/4.0 ; // 135 degrees
    FORALL(cells,cc) {
      int csz = cellStencil[cc].size() ;
      int nsz = neighStencil[cc].size() ;
      int bsz = boundary_map[cc].size() ;

      vector3d<double>  ccent = ccenter[cc] ;
      vector<vector3d<double> > cdirs(csz+bsz) ;
      for(int i=0;i<csz;++i) {
        cdirs[i] = ccenterdata[cg2l[cellStencil[cc][i]]]-ccent ;
        cdirs[i] *= 1./norm(cdirs[i]) ;
      }
      for(int i=0;i<bsz;++i) {
        cdirs[csz+i] = fcenterdata[fg2l[boundary_map[cc][i]]]-ccent ;
        cdirs[csz+i] *= 1./norm(cdirs[csz+i]) ;
      }

      vector<int> flags(csz+bsz,0) ;
      for(int i=0;i<nsz;++i) {
        vector3d<double> target = ncenterdata[ncg2l[neighStencil[cc][i]]] ;
        vector3d<double>  ejk = target-ccent ;
        ejk *= 1./norm(ejk) ;

        int minid = 0 ;
        double dvmin = 1e13 ;
        int nid = -1 ;
        for(int j=0;j<csz;++j) {
          // exclude itself from the search with theta ...
          double dv = dot(cdirs[j],ejk) ;
          if(dv<cos(theta) && dv<dvmin) {
            minid = j ;
            dvmin = dv ;
          }
          // find neighbor id
          if(cellStencil[cc][j]==neighStencil[cc][i]) {
            nid = j ;
          }
        }
        if(nid==-1) {
            cerr <<"symm stencil error: could not find neighbor in "
                 << "cellStencil list\n\n" << endl ;
            Loci::Abort() ;
        }
        flags[nid] = 1 ;
        flags[minid] = 1 ;
        }

      for(int i=0;i<bsz;++i) {
        vector3d<double> target = fcenterdata[fg2l[boundary_map[cc][i]]] ;
        vector3d<double>  ejk = target-ccent ;
        ejk *= 1./norm(ejk) ;

        int minid = 0 ;
        double dvmin = 1e13 ;
        for(int j=0;j<csz;++j) {
          // exclude itself from the search with theta ...
          double dv = dot(cdirs[j],ejk) ;
          if(dv<cos(theta) && dv<dvmin) {
            minid = j ;
            dvmin = dv ;
          }
        }
        flags[minid] = 1 ;
      }

      int cnt = 0 ;
      for(int i=0;i<csz;++i) {
        if(flags[i] > 0) {
          cellmap.push_back(cellStencil[cc][i]) ;
          cnt++ ;
        }
      }
      sizes[cc] = cnt ;
    } ENDFORALL;

    cellStencilSymm.allocate(sizes) ;
    int cnt = 0 ;
    FORALL(cells,cc) {
      for(int i=0;i<cellStencilSymm[cc].size();++i) {
        cellStencilSymm[cc][i] = cellmap[cnt] ;
        cnt++ ;
      }
    } ENDFORALL;
  }

  void create_cell_stencil_symm(fact_db & facts) {
    // Create cell stencil map from protoMap
    multiMap cellStencil ;
    get_symm_cellStencil(cellStencil,facts) ;
    // Put in fact database
    facts.create_fact("cellStencil",cellStencil) ;
  }

  double Faug(const int nft, const vector3d<double> x1,
              const tmp_array<vector3d<double>> &xneigh) {
    const int nvar = 4 ;
    tmp_array<double> Ac(nvar*nvar) ;
    for(int i=0; i<nvar*nvar; i++) {
      Ac[i] = 0. ;
    }

    tmp_array<double> wi(nft) ;
    tmp_array<vector3d<double>> dr(nft) ;
    // wi
    int nf = 0 ;
    for(int i=0;i<nft;i++) {
      const vector3d<double> x0 = xneigh[i] ;
      const vector3d<double> ds = (x1-x0) ;
      dr[nf] = ds ;
      wi[nf++] = 1./norm(ds) ;
    }

    // A matrix
    int row = 0;
    for(int f=0; f<nft; f++) {
      row=0;
      Ac[row*nvar+0] += wi[f];
      Ac[row*nvar+1] += wi[f]*dr[f].x;
      Ac[row*nvar+2] += wi[f]*dr[f].y;
      Ac[row*nvar+3] += wi[f]*dr[f].z;

      row=1;
      Ac[row*nvar+0] += wi[f]*dr[f].x;
      Ac[row*nvar+1] += wi[f]*dr[f].x*dr[f].x;
      Ac[row*nvar+2] += wi[f]*dr[f].y*dr[f].x;
      Ac[row*nvar+3] += wi[f]*dr[f].z*dr[f].x;

      row=2;
      Ac[row*nvar+0] += wi[f]*dr[f].y;
      Ac[row*nvar+1] += wi[f]*dr[f].x*dr[f].y;
      Ac[row*nvar+2] += wi[f]*dr[f].y*dr[f].y;
      Ac[row*nvar+3] += wi[f]*dr[f].z*dr[f].y;

      row=3;
      Ac[row*nvar+0] += wi[f]*dr[f].z;
      Ac[row*nvar+1] += wi[f]*dr[f].x*dr[f].z;
      Ac[row*nvar+2] += wi[f]*dr[f].y*dr[f].z;
      Ac[row*nvar+3] += wi[f]*dr[f].z*dr[f].z;
    }

    // A^T A, Frobenius Norm of A^T*A
    double FBNorm = 0.0 ;
    for(int i=0; i<nvar; i++) {
      for(int col=0; col<nvar; col++) {
        double sum = 0.;
        for(int k=0; k<nvar; k++) {
          // AT[i*nvar+k] = Ac[k*nvar+i]
          const double AT = Ac[k*nvar+i] ;
          const double A  = Ac[k*nvar+col] ;
          sum += AT*A ;
        }
        FBNorm += sum*sum ;
      }
    }
    FBNorm = sqrt(FBNorm) ;

    double s = 0.0;
    for (int f=0; f<nft; f++) {
      const double dm = norm(dr[f]) ;
      s += wi[f]*wi[f]*dm ;
    }

    double F = s/FBNorm ;
    return F ;

  }

  void get_symmF_cellStencil(multiMap &cellStencilSymmF, fact_db & facts) {
    if(Loci::MPI_rank==0) {
      cout << "Generating symmF stencil" << endl ;
    }
    using std::vector ;
    // get full stencil
    multiMap cellStencil ;
    get_full_cellStencil(cellStencil,facts) ;

    entitySet cells = cellStencil.domain();
    multiMap upper,lower,boundary_map ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;

    entitySet faceimage  ;
    faceimage += Loci::MapRepP(upper.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(lower.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(boundary_map.Rep())->image(cells) ;
    int fkeyspace = upper.getRangeKeySpace() ;
    std::vector<entitySet> face_ptn = facts.get_init_ptn(fkeyspace) ;
    store<vector3d<double> > fcenter,normal ;
    store<double> area ;
    getFaceCenter(facts,fcenter,area,normal) ;
    store<vector3d<double> > ccenter ;
    getCellCenter(facts,ccenter,fcenter,area) ;

    // gather face data needed for the computation.
    std::map<int,int> fg2l ;
    getLocalContextMap(fg2l,faceimage) ;
    vector<vector3d<double> > fcenterdata ;
    vector<vector3d<double> > fnormaldata ;
    gatherData(fcenterdata,fcenter,faceimage,face_ptn) ;
    gatherData(fnormaldata,normal,faceimage,face_ptn) ;

    int ckeyspace = upper.getDomainKeySpace() ;
    std::vector<entitySet> cell_ptn = facts.get_init_ptn(ckeyspace) ;

    entitySet cellImage = Loci::MapRepP(cellStencil.Rep())->image(cells) ;

    vector<vector3d<double> > ccenterdata ;
    std::map<int,int> cg2l ;
    getLocalContextMap(cg2l,cellImage) ;
    gatherData(ccenterdata,ccenter,cellImage,cell_ptn) ;

    multiMap neighStencil;
    get_neigh_cellStencil(neighStencil,facts);
    entitySet neighCellImage = Loci::MapRepP(neighStencil.Rep())->image(cells);
    vector<vector3d<double> > ncenterdata ;
    std::map<int,int> ncg2l ;
    getLocalContextMap(ncg2l,neighCellImage) ;
    gatherData(ncenterdata,ccenter,neighCellImage,cell_ptn) ;

    store<int> sizes ;
    sizes.allocate(cells) ;
    vector<int> cellmap ;
    double theta = 3.*M_PI/4.0 ;
    FORALL(cells,cc) {
      int csz = cellStencil[cc].size() ;
      int nsz = neighStencil[cc].size() ;
      int bsz = boundary_map[cc].size() ;

      vector3d<double> ccent = ccenter[cc] ;
      vector<vector3d<double> > cdirs(csz+bsz) ;
      for(int i=0;i<csz;++i) {
        cdirs[i] = ccenterdata[cg2l[cellStencil[cc][i]]]-ccent ;
        cdirs[i] *= 1./norm(cdirs[i]) ;
      }
      for(int i=0;i<bsz;++i) {
        cdirs[csz+i] = fcenterdata[fg2l[boundary_map[cc][i]]]-ccent ;
        cdirs[csz+i] *= 1./norm(cdirs[csz+i]) ;
      }

      vector<int> flags(csz+bsz,0) ;
      for(int i=0;i<nsz;++i) {
        vector3d<double> target = ncenterdata[ncg2l[neighStencil[cc][i]]] ;
        vector3d<double> ejk = target-ccent ;
        ejk *= 1./norm(ejk) ;

        int minid = 0 ;
        double dvmin = 1e13 ;
        int nid = -1 ;
        for(int j=0;j<csz;++j) {
          // exclude itself from the search with theta ...
          double dv = dot(cdirs[j],ejk) ;
          if(dv<cos(theta) && dv<dvmin) {
            minid = j ;
            dvmin = dv ;
          }
          // find neighbor id
          if(cellStencil[cc][j]==neighStencil[cc][i]) {
            nid = j ;
          }
        }
        if(nid==-1) {
          cerr <<"symmF stencil error: could not find neighbor in cellStencil list\n\n" << endl ;
          Loci::Abort() ;
        }
        flags[nid] = 1 ;
        flags[minid] = 1 ;
      }

      for(int i=0;i<bsz;++i) {
        vector3d<double> target = fcenterdata[fg2l[boundary_map[cc][i]]] ;
        vector3d<double>  ejk = target-ccent ;
        ejk *= 1./norm(ejk) ;

        int minid = 0 ;
        double dvmin = 1e13 ;
        for(int j=0;j<csz;++j) {
          // exclude itself from the search with theta ...
          double dv = dot(cdirs[j],ejk) ;
          if(dv<cos(theta) && dv<dvmin) {
            minid = j ;
            dvmin = dv ;
          }
        }
        flags[minid] = 1 ;
      }

      // now do the augmentation
      tmp_array<vector3d<double> > tmpcc(csz+bsz) ;
      int ne = 0 ;
      int check_sten = 0 ;
      for(int i=0;i<csz;++i) {
        if(flags[i]) {
          tmpcc[ne++] = ccenterdata[cg2l[cellStencil[cc][i]]] ;
        }
        check_sten += (!flags[i]?1:0) ;
      }

      for(int i=0;i<bsz;++i) {
        tmpcc[ne++] = fcenterdata[fg2l[boundary_map[cc][i]]] ;
      }

      double F0 = Faug(ne,ccent,tmpcc) ;

      for(int k=0; k<check_sten;k++) {
        for(int i=0;i<csz;++i) {
          if(flags[i]) {
            continue ;
          }

          // build stencil and add this cell to the stencil
          int ne = 0 ;
          for(int j=0;j<csz;++j) {
            if(flags[j] || (i==j)) {
              tmpcc[ne++] = ccenterdata[cg2l[cellStencil[cc][j]]] ;
            }
          }

          for(int j=0;j<bsz;++j) {
            tmpcc[ne++] = fcenterdata[fg2l[boundary_map[cc][j]]] ;
          }

          double F = Faug(ne,ccent,tmpcc);
          if(F<0.85*F0) {
            flags[i] = 1 ;
            F0 = F ;
          }
        }
      }

      int cnt = 0 ;
      for(int i=0;i<csz;++i) {
        if(flags[i] > 0) {
          cellmap.push_back(cellStencil[cc][i]) ;
          cnt++ ;
        }
      }
      sizes[cc] = cnt ;
    } ENDFORALL ;

    cellStencilSymmF.allocate(sizes) ;
    int cnt = 0 ;
    FORALL(cells,cc) {
      for(int i=0;i<cellStencilSymmF[cc].size();++i) {
        cellStencilSymmF[cc][i] = cellmap[cnt] ;
        cnt++ ;
      }
    } ENDFORALL ;
  }

  void create_cell_stencil_symmF(fact_db & facts) {
    // Create cell stencil map from protoMap
    multiMap cellStencil ;
    get_symmF_cellStencil(cellStencil,facts) ;
    // Put in fact database
    facts.create_fact("cellStencil",cellStencil) ;
  }

  template<class T, class T2> void lu_4x4( T2 *A, T *x, T *b) {
    T lvar_0 = 1.0/A[0] ;
    T lvar_1 = A[12]*lvar_0 ;
    T lvar_2 = A[4]*lvar_0 ;
    T lvar_3 = -b[0]*lvar_2 + b[1] ;
    T lvar_3p5 = -A[1]*lvar_2 + A[5] ;
    T sign   = T(lvar_3p5<0.?-1.:1.) ;
    T lvar_4 = sign/max<T>(abs(lvar_3p5),T(1e-32)) ;
    T lvar_5 = lvar_4*(A[13] - A[1]*lvar_1) ;
    T lvar_6 = A[8]*lvar_0 ;
    T lvar_7 = lvar_4*(-A[1]*lvar_6 + A[9]) ;
    T lvar_8 = -b[0]*lvar_6 + b[2] - lvar_3*lvar_7 ;
    T lvar_9 = -A[2]*lvar_2 + A[6] ;
    T lvar_9p5 = A[10] - A[2]*lvar_6 - lvar_7*lvar_9 ;
    sign    = T(lvar_9p5<0.?-1.:1.) ;
    T lvar_10 = sign/max<T>(abs(lvar_9p5),T(1e-32)) ;
    T lvar_11 = lvar_10*(A[14] - A[2]*lvar_1 - lvar_5*lvar_9) ;
    T lvar_12 = -A[3]*lvar_2 + A[7] ;
    T lvar_13 = A[11] - A[3]*lvar_6 - lvar_12*lvar_7 ;
    T lvar_13p5 = (A[15] - A[3]*lvar_1 - lvar_11*lvar_13 - lvar_12*lvar_5) ;
    sign    = T(lvar_13p5<0.?-1.:1.) ;
    T lvar_14 = (-b[0]*lvar_1 + b[3] - lvar_11*lvar_8 - lvar_3*lvar_5)*sign/max<T>(abs(lvar_13p5),T(1e-32)) ;
    T lvar_15 = lvar_10*(-lvar_13*lvar_14 + lvar_8) ;
    T lvar_16 = lvar_4*(-lvar_12*lvar_14 - lvar_15*lvar_9 + lvar_3) ;
    x[0] = lvar_0*(-A[1]*lvar_16 - A[2]*lvar_15 - A[3]*lvar_14 + b[0]) ;
    x[1] = lvar_16 ;
    x[2] = lvar_15 ;
    x[3] = lvar_14 ;
  }

  template<class T> void inv_from_lu(T *A, T *Ainv) {
    int size = 4;
    T x[4];
    T b[4];
    for(int row=0; row<size; row++) {
      for(int i=0; i<size; i++) {
        b[i] = 0 ;
      }
      b[row] = 1. ;
      lu_4x4(A,&x[0],&b[0]) ;
      for(int i=0; i<size; i++) {
        Ainv[i*size+row] = x[i] ;
      }
    }
  }

  void symmC_optimize(const int nft, const vector3d<double> x1,
                      const tmp_array<vector3d<double>> &xneigh,
                      double &p1, double &p2, double &pinf) {
    const int nvar = 4 ;
    tmp_array<double> Ac(nvar*nvar) ;
    tmp_array<double> Ainv(nvar*nvar) ;
    tmp_array<double> wi(nft) ;

    for (int i=0; i<nvar*nvar; i++) {
      Ac[i] = 0. ;
      Ainv[i] = 0. ;
    }
    vector<vector3d<double>> dr(nft) ;
    // wi
    int nf = 0 ;
    for(int i=0;i<nft;i++) {
      const vector3d<double> x0 = xneigh[i] ;
      vector3d<double> ds = (x1-x0) ;
      dr[nf] = ds ;
      wi[nf++] = 1./norm(ds) ;
    }

    // A matrix
    int row = 0 ;
    for (int f=0; f<nft; f++) {
      row=0 ;
      Ac[row*nvar+0] += wi[f] ;
      Ac[row*nvar+1] += wi[f]*dr[f].x ;
      Ac[row*nvar+2] += wi[f]*dr[f].y ;
      Ac[row*nvar+3] += wi[f]*dr[f].z ;

      row=1 ;
      Ac[row*nvar+0] += wi[f]*dr[f].x;
      Ac[row*nvar+1] += wi[f]*dr[f].x*dr[f].x ;
      Ac[row*nvar+2] += wi[f]*dr[f].y*dr[f].x ;
      Ac[row*nvar+3] += wi[f]*dr[f].z*dr[f].x ;

      row=2 ;
      Ac[row*nvar+0] += wi[f]*dr[f].y ;
      Ac[row*nvar+1] += wi[f]*dr[f].x*dr[f].y ;
      Ac[row*nvar+2] += wi[f]*dr[f].y*dr[f].y ;
      Ac[row*nvar+3] += wi[f]*dr[f].z*dr[f].y ;

      row=3 ;
      Ac[row*nvar+0] += wi[f]*dr[f].z ;
      Ac[row*nvar+1] += wi[f]*dr[f].x*dr[f].z ;
      Ac[row*nvar+2] += wi[f]*dr[f].y*dr[f].z ;
      Ac[row*nvar+3] += wi[f]*dr[f].z*dr[f].z ;
    }

    // A^-1
    inv_from_lu<double>(&Ac[0], &Ainv[0]) ;

    // norms of Ac
    double normA_p1   = 0. ;
    double normA_pinf = 0. ;
    double normA_p2   = 0. ;
    for (int i=0; i<nvar; i++) {
      double tmp_p1   = 0. ;
      double tmp_pinf = 0. ;
      for (int j=0; j<nvar; j++) {
        tmp_p1   += abs(Ac[j*nvar+i]) ;
        tmp_pinf += abs(Ac[i*nvar+j]) ;
        normA_p2 += Ac[i*nvar+j]*Ac[i*nvar+j] ;
      }
      normA_p1 = max<double>(normA_p1, tmp_p1) ;
      normA_pinf = max<double>(normA_pinf, tmp_pinf) ;
    }
    normA_p2 = sqrt(normA_p2) ;

    // norms of Ainv
    double normAi_p1   = 0. ;
    double normAi_pinf = 0. ;
    double normAi_p2   = 0. ;
    for (int i=0; i<nvar; i++) {
      double tmp_p1   = 0. ;
      double tmp_pinf = 0. ;
      for (int j=0; j<nvar; j++) {
        tmp_p1 += abs(Ainv[j*nvar+i]) ;
        tmp_pinf += abs(Ainv[i*nvar+j]) ;
        normAi_p2 += Ainv[i*nvar+j]*Ainv[i*nvar+j] ;
      }
      normAi_p1 = max<double>(normAi_p1, tmp_p1) ;
      normAi_pinf = max<double>(normAi_pinf, tmp_pinf) ;
    }
    normAi_p2 = sqrt(normAi_p2);

    // condition number
    p1 = normA_p1*normAi_p1 ;
    p2 = normA_p2*normAi_p2 ;
    pinf = normA_pinf*normAi_pinf ;
  }

  void get_symmC_cellStencil(multiMap &cellStencilSymmC, fact_db & facts) {
    if(Loci::MPI_rank==0) {
      cout <<"Generating symmC stencil" << endl ;
    }
    using std::vector ;
    // get full stencil
    multiMap cellStencil ;
    get_full_cellStencil(cellStencil, facts) ;

    entitySet cells = cellStencil.domain() ;
    multiMap upper,lower,boundary_map ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;

    entitySet faceimage ;
    faceimage += Loci::MapRepP(upper.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(lower.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(boundary_map.Rep())->image(cells) ;
    int fkeyspace = upper.getRangeKeySpace() ;
    std::vector<entitySet> face_ptn = facts.get_init_ptn(fkeyspace) ;
    store<vector3d<double> > fcenter,normal ;
    store<double> area ;
    getFaceCenter(facts,fcenter,area,normal) ;
    store<vector3d<double> > ccenter ;
    getCellCenter(facts,ccenter,fcenter,area) ;

    // gather face data needed for the computation.
    std::map<int,int> fg2l ;
    getLocalContextMap(fg2l,faceimage) ;
    vector<vector3d<double> > fcenterdata ;
    vector<vector3d<double> > fnormaldata ;
    gatherData(fcenterdata,fcenter,faceimage,face_ptn) ;
    gatherData(fnormaldata,normal,faceimage,face_ptn) ;

    int ckeyspace = upper.getDomainKeySpace() ;
    std::vector<entitySet> cell_ptn = facts.get_init_ptn(ckeyspace) ;

    entitySet cellImage = Loci::MapRepP(cellStencil.Rep())->image(cells) ;

    vector<vector3d<double> > ccenterdata ;
    std::map<int,int> cg2l ;
    getLocalContextMap(cg2l,cellImage) ;
    gatherData(ccenterdata,ccenter,cellImage,cell_ptn) ;

    multiMap neighStencil;
    get_neigh_cellStencil(neighStencil,facts);
    entitySet neighCellImage = Loci::MapRepP(neighStencil.Rep())->image(cells);
    vector<vector3d<double> > ncenterdata ;
    std::map<int,int> ncg2l ;
    getLocalContextMap(ncg2l,neighCellImage) ;
    gatherData(ncenterdata,ccenter,neighCellImage,cell_ptn) ;

    store<int> sizes ;
    sizes.allocate(cells) ;
    vector<int> cellmap ;
    double theta = 3.*M_PI/4.0 ;
    FORALL(cells,cc) {
      int csz = cellStencil[cc].size() ;
      int nsz = neighStencil[cc].size() ;
      int bsz = boundary_map[cc].size() ;

      vector3d<double> ccent = ccenter[cc] ;
      vector<vector3d<double> > cdirs(csz+bsz) ;
      for(int i=0;i<csz;++i) {
        cdirs[i] = ccenterdata[cg2l[cellStencil[cc][i]]]-ccent ;
        cdirs[i] *= 1./norm(cdirs[i]) ;
      }
      for(int i=0;i<bsz;++i) {
        cdirs[csz+i] = fcenterdata[fg2l[boundary_map[cc][i]]]-ccent ;
        cdirs[csz+i] *= 1./norm(cdirs[csz+i]) ;
      }

      vector<int> flags(csz+bsz,0) ;
      for(int i=0;i<nsz;++i) {
        vector3d<double> target = ncenterdata[ncg2l[neighStencil[cc][i]]] ;
        vector3d<double>  ejk = target-ccent ;
        double nejk = norm(ejk) ;
        ejk *= 1./nejk ;

        int minid = 0 ;
        double dvmin = 1e13 ;
        int nid = -1 ;
        for(int j=0;j<csz;++j) {
          // exclude itself from the search with theta ...
          double dv = dot(cdirs[j],ejk) ;
          if(dv<cos(theta) && dv<dvmin) {
            minid = j ;
            dvmin = dv ;
          }
          // find neighbor id
          if(cellStencil[cc][j]==neighStencil[cc][i]) {
            nid = j ;
          }
        }
        flags[nid] = 1 ;
        flags[minid] = 1 ;
      }

      for(int i=0;i<bsz;++i) {
        vector3d<double> target = fcenterdata[fg2l[boundary_map[cc][i]]] ;
        vector3d<double>  ejk = target-ccent ;
        double nejk = norm(ejk) ;
        ejk *= 1./nejk ;

        int minid = 0 ;
        double dvmin = 1e13 ;
        for(int j=0;j<csz;++j) {
          // exclude itself from the search with theta ...
          double dv = dot(cdirs[j],ejk) ;
          if(dv<cos(theta) && dv<dvmin) {
            minid = j ;
            dvmin = dv ;
          }
        }
        flags[minid] = 1 ;
      }

      // now do the augmentation
      tmp_array<vector3d<double> > tmpcc(csz+bsz) ;
      int ne = 0 ;
      int check_sten = 0 ;
      for(int i=0;i<csz;++i) {
        if(flags[i]) {
          tmpcc[ne++] = ccenterdata[cg2l[cellStencil[cc][i]]] ;
        }
        check_sten += (!flags[i]?1:0) ;
      }

      for(int i=0;i<bsz;++i) {
        tmpcc[ne++] = fcenterdata[fg2l[boundary_map[cc][i]]] ;
      }

      double cn_p1 = 0. ;
      double cn_p2 = 0. ;
      double cn_pinf = 0. ;
      symmC_optimize(ne, ccent, tmpcc, cn_p1, cn_p2, cn_pinf) ;

      for(int k=0; k<check_sten;k++) {
        double cn_opt_p1 = 0. ;
        double cn_opt_p2 = 0. ;
        double cn_opt_pinf = 0. ;
        for(int i=0;i<csz;++i) {
          if(flags[i]) {
            continue ;
          }

          int ne = 0 ;
          for(int j=0;j<csz;++j) {
            if(flags[j] || (i==j)) {
              tmpcc[ne++] = ccenterdata[cg2l[cellStencil[cc][j]]] ;
            }
          }

          for(int j=0;j<bsz;++j) {
            tmpcc[ne++] = fcenterdata[fg2l[boundary_map[cc][j]]] ;
          }

          symmC_optimize(ne, ccent, tmpcc, cn_opt_p1, cn_opt_p2, cn_opt_pinf) ;
          if(cn_opt_p1<0.95*cn_p1 && cn_opt_p2<0.95*cn_p2 &&
             cn_opt_pinf<0.95*cn_pinf) {
            flags[i] = 1 ;
            cn_p1 = cn_opt_p1 ;
            cn_p2 = cn_opt_p2 ;
            cn_pinf = cn_opt_pinf ;
          }
        }
      }

      int cnt = 0 ;
      for(int i=0;i<csz;++i) {
        if(flags[i] > 0) {
          cellmap.push_back(cellStencil[cc][i]) ;
          cnt++ ;
        }
      }
      sizes[cc] = cnt ;
    } ENDFORALL ;

    cellStencilSymmC.allocate(sizes) ;
    int cnt = 0 ;
    FORALL(cells,cc) {
      for(int i=0;i<cellStencilSymmC[cc].size();++i) {
        cellStencilSymmC[cc][i] = cellmap[cnt] ;
        cnt++ ;
      }
    } ENDFORALL ;
  }

  void create_cell_stencil_symmC(fact_db & facts) {
    // Create cell stencil map from protoMap
    multiMap cellStencil ;
    get_symmC_cellStencil(cellStencil,facts) ;
    // Put in fact database
    facts.create_fact("cellStencil",cellStencil) ;
  }


  void createLowerUpper(fact_db &facts) {
    constraint geom_cells,interior_faces,boundary_faces;
    constraint faces = facts.get_variable("faces") ;
    geom_cells = facts.get_variable("geom_cells") ;
    int ckeyspace = geom_cells.getDomainKeySpace() ;
    interior_faces = facts.get_variable("interior_faces") ;
    boundary_faces = facts.get_variable("boundary_faces") ;
    entitySet bfaces = *boundary_faces ;
    entitySet ifaces = *interior_faces ;

    storeRepP pfacesP = facts.get_variable("periodicFaces") ;
    if(pfacesP != 0) {
      constraint periodicFaces ;
      periodicFaces = pfacesP ;
      bfaces -= *periodicFaces ;
      ifaces += *periodicFaces ;
    }

    int fkeyspace = faces.getDomainKeySpace() ;

    std::vector<entitySet> fptn = facts.get_init_ptn(fkeyspace) ;
    entitySet
      global_interior_faces = (distribute_entitySet(ifaces,fptn) & (*faces)) ;

    entitySet
      global_boundary_faces = (distribute_entitySet(bfaces,fptn) & (*faces)) ;

    Map cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    // Note, all_collect_entitySet has the potential to be inefficient,
    // but not in the present case.  These low level utilities need to be
    // rethought.
    entitySet global_geom_cells = all_collect_entitySet(*geom_cells) ;
    multiMap lower,upper,boundary_map ;
    distributed_inverseMap(upper, cl, global_geom_cells, global_interior_faces,
                           facts,ckeyspace) ;

    distributed_inverseMap(lower, cr, global_geom_cells, global_interior_faces,
                           facts,ckeyspace) ;
    distributed_inverseMap(boundary_map, cl, global_geom_cells,
                           global_boundary_faces, facts,ckeyspace) ;

    MapRepP clr = MapRepP(cl.Rep()) ;
    MapRepP crr = MapRepP(cr.Rep()) ;
    lower.Rep()->setDomainKeySpace(crr->getRangeKeySpace()) ;
    upper.Rep()->setDomainKeySpace(clr->getRangeKeySpace()) ;
    boundary_map.Rep()->setDomainKeySpace(clr->getRangeKeySpace()) ;
    MapRepP(lower.Rep())->setRangeKeySpace(crr->getDomainKeySpace()) ;
    MapRepP(upper.Rep())->setRangeKeySpace(clr->getDomainKeySpace()) ;
    MapRepP(boundary_map.Rep())->setRangeKeySpace(clr->getDomainKeySpace()) ;
    facts.create_fact("lower",lower) ;
    facts.create_fact("upper",upper) ;
    facts.create_fact("boundary_map",boundary_map) ;

    param<std::string> gradStencil ;
    storeRepP var = facts.get_variable("gradStencil") ;
    if(var != 0) {
      gradStencil = var ;
      if(*gradStencil == "stable")
	      create_cell_stencil(facts) ;
      if(*gradStencil == "full")
	      create_cell_stencil_full(facts) ;
      if(*gradStencil == "symm")
        create_cell_stencil_symm(facts) ;
      if(*gradStencil == "symmF")
        create_cell_stencil_symmF(facts) ;
      if(*gradStencil == "symmC")
        create_cell_stencil_symmC(facts) ;
    }
  }

  extern int factdb_allocated_base ;
  namespace { // Utilities for createEdgesPar
    
    /// Comparison that only pays attention to the first element in the array
    template<class T, size_t N> inline bool split_sort_1(const Array<T,N> &k1,
                                                         const Array<T,N> &k2) {
      return k1[0] < k2[0] ;
    }

    /// Comparison that only pays attention to the first two elements in the array
    template<class T, size_t N> inline bool split_sort_12(const Array<T,N> &k1,
                                                          const Array<T,N> &k2) {
      return k1[0] < k2[0] || (k1[0] == k2[0] && k1[1] < k2[1]) ;
    }

    /// Sort edge keys, if edge_splits is empty then generate the splits
    /// from the data, otherwise use splits provided
    template<class T, size_t N>
    void edge_sort(std::vector<Array<T,N> > &edge_keys,
                   std::vector<Array<T,N> > &edge_splits,
                   MPI_Comm comm) {
      std::sort(edge_keys.begin(),edge_keys.end(),split_sort_1<T,N>) ;
      int p = 1 ;
      MPI_Comm_size(comm,&p) ;
      // edge splits not provided, compute them using uniform algorithmic
      // partition
      if(edge_splits.size() == 0 && p > 1) {
        Array<T,N> zero ;
        for(size_t i=0;i<N;++i)
          zero[i] = 0 ;
        std::vector<Array<T,N> > tmp(p-1,zero) ;
        T emaxl = edge_keys[edge_keys.size()-1][0] ;
        T eminl = edge_keys[0][0] ;
        T emax = emaxl ;
        MPI_Allreduce(&emaxl,&emax,1,MPI_INT,MPI_MAX,comm) ;
        T emin = eminl ;
        MPI_Allreduce(&eminl,&emin,1,MPI_INT,MPI_MIN,comm) ;
        T delta = (emax/p)- (emin/p) + 1;
        edge_splits.swap(tmp) ;
        edge_splits[0][0] = emin+delta ;
        for(int i=1;i<p-1;++i) {
          edge_splits[i][0] = edge_splits[i-1][0] + delta ;
        }
      }
      parSplitSort(edge_keys,edge_splits,split_sort_1<T,N>,comm) ;
      std::sort(edge_keys.begin(),edge_keys.end(),split_sort_12<T,N>) ;
    }

    /// create edge keys from the face2node map.  The keys are stored
    /// as the file numbering of the two edge nodes first, then the
    /// global numbering second, and finally the face that generated
    /// the edge. To ensure a consistent numbering with the VOG file,
    /// edges are oriented from lowest file number node to largest
    /// file number node. Note, edges are duplicated by each face
    /// that shares the node.  Sorting will be used to identify
    /// unique edges.
    ///
    /// -- edge_key is an array of 5 values:
    /// -- edge_key[0] - file number of first node defined by edge (lowest file number)
    /// -- edge_key[1] - file number of second node defined by edge
    /// -- edge_key[2] - global number of first node defined by edge
    /// -- edge_key[3] - global number of second node defined by edge
    /// -- edge_key[4] - global number of face that generated this edge
    ///
    template<class T> void extract_edge_keys(std::vector<Array<T,5> > &edge_keys,
                                             entitySet fdom,
                                             entitySet ndom,
                                             int nkeyspace,
                                             multiMap &face2node,
                                             fact_db &facts) {
      // Find nodes that face2node accesses
      entitySet node_access = Loci::MapRepP(face2node.Rep())->image(fdom) ;

      std::vector<T> n2f ; // get node2file numbering data
      std::vector<entitySet> node_ptn = facts.get_init_ptn(nkeyspace) ;

      store<T> node_fnum;
      node_fnum.allocate(ndom) ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;

      if(df != 0) {
        dMap ng2f ;
        ng2f = df->g2fv[nkeyspace].Rep() ;
        FORALL(ndom,ii) {
          node_fnum[ii] = ng2f[ii] ;
        } ENDFORALL ;
      } else {
        FORALL(ndom,ii) {
          node_fnum[ii] = ii ;
        } ENDFORALL ;
      }
    
      // gather nodal data needed for the computation.
      std::map<int,int> g2l ;
      getLocalContextMap(g2l,node_access) ;
      gatherData(n2f,node_fnum,node_access,node_ptn) ;
      // Create list of edge keys, first count how many edge keys we will
      // generate
      size_t edge_gen = 0 ;
      FORALL(fdom,fc) {
        int fsz = face2node[fc].size() ;
        edge_gen += fsz ;
      } ENDFORALL ;
      // Allocate edge keys
      {
        vector<Array<T,5> > tmp(edge_gen) ;
        edge_keys.swap(tmp) ;
      }
      size_t edge_local = 0 ;
      // Loop over faces and create an edge key for each edge of the face
      // Record both the global and file number for the two nodes that form
      // the edge
      FORALL(fdom,fc) {
        int fsz = face2node[fc].size() ;
        T face = fc ;
        for(int i=1;i<fsz;++i) {
          T gn1 = face2node[fc][i-1] ;
          T fn1 = n2f[g2l[gn1]] ;
          T gn2 = face2node[fc][i] ;
          T fn2 = n2f[g2l[gn2]] ;
          edge_keys[edge_local][0] = fn1 ;
          edge_keys[edge_local][1] = fn2 ;
          edge_keys[edge_local][2] = gn1 ;
          edge_keys[edge_local][3] = gn2 ;
          edge_keys[edge_local][4] = face ;
          if(fn1 > fn2) {
            std::swap(edge_keys[edge_local][0],edge_keys[edge_local][1]) ;
            std::swap(edge_keys[edge_local][2],edge_keys[edge_local][3]) ;
          }
          edge_local++ ;
        }
        // Final edge
        T gn1 = face2node[fc][0] ;
        T fn1 = n2f[g2l[gn1]] ;
        T gn2 = face2node[fc][fsz-1] ;
        T fn2 = n2f[g2l[gn2]] ;
        edge_keys[edge_local][0] = fn1 ;
        edge_keys[edge_local][1] = fn2 ;
        edge_keys[edge_local][2] = gn1 ;
        edge_keys[edge_local][3] = gn2 ;
        edge_keys[edge_local][4] = face ;
        if(fn1 > fn2) {
          std::swap(edge_keys[edge_local][0],edge_keys[edge_local][1]) ;
          std::swap(edge_keys[edge_local][2],edge_keys[edge_local][3]) ;
        }
        edge_local++ ;
      } ENDFORALL ;
      FATAL(edge_local != edge_gen) ;
    }    


    /// Scan through edges sorted to processors according to the file
    /// numbering of nodes.  Use this to define a unique file ordering
    /// for edges.  As we scan through the edge array, convert the key
    /// to global nubering for following steps, and then fill in
    /// field 2 with the edge file number
    ///
    /// -- After this step, edge_key is redefined as
    /// -- edge_key[0] - global number of first node that forms edge
    /// -- edge_key[1] - global number of second node that forms edge
    /// -- edge_key[2] - number of edges in duplicate edge cluster
    /// -- edge_key[3] - edge number in file numbering
    /// -- edge_key[4] - global number of face generating this edge
    //
    template<class T>
    size_t edge_file_numbering(std::vector<Array<T,5> > &edge_keys,
                               MPI_Comm comm) {
      long long ecnt = 0 ;
      for(size_t i=0;i<edge_keys.size();) {
        // search for cluster of edges
        size_t j = i+1 ;
        for(;j<edge_keys.size();++j) {
          if((edge_keys[i][0] != edge_keys[j][0]) ||
             (edge_keys[i][1] != edge_keys[j][1]))
            break ;
        }
        size_t nedges = j-i ;
        for(size_t k=i;k<j;++k) {
          // Convert key to global numbering 
          edge_keys[k][0] = edge_keys[k][2] ;
          edge_keys[k][1] = edge_keys[k][3] ;

          FATAL(edge_keys[k][0] != edge_keys[i][0] ||
                edge_keys[k][1] != edge_keys[i][1]) ;
          // Use remaining space to store information about edges
          // Edge file number (starting from zero on each processor
          edge_keys[k][2] = ecnt ;
          // Number of edges in cluster
          edge_keys[k][3] = nedges ;
        }
        ecnt++ ;
        i=j ;
      }
      // Fix local numbering to file numbering
      long long ecnt_total = ecnt ;
      MPI_Scan(&ecnt, &ecnt_total,1,MPI_LONG_LONG,MPI_SUM,comm) ;
      long long offset = ecnt_total - ecnt + factdb_allocated_base ;
      for(size_t i=0;i<edge_keys.size();++i) 
        edge_keys[i][2] = (edge_keys[i][2]+offset) ;
      return ecnt ;
    }

    /// Extract splits consistent with the partition set given, assumes
    /// the set allocated to each processor is contiguous.
    template<class T>
    void get_splits(std::vector<Array<T,5> > &edge_splits,
                    entitySet ndom,
                    MPI_Comm comm) {
      int p=1 ;
      MPI_Comm_size(comm,&p) ;
      { // Initialize splits
        Array<T,5> zero ;
        for(size_t i=0;i<5;++i)
          zero[i] = 0 ;
        std::vector<Array<T,5> > tmp(p-1,zero) ;
        edge_splits.swap(tmp) ;
      }
      T val = ndom.Min() ;
      std::vector<T> valp(p) ;
      MPI_Allgather(&val,sizeof(T),MPI_BYTE,
                    &valp[0],sizeof(T),MPI_BYTE,
                    comm) ;
      for(int i=0;i<p-1;++i)
        edge_splits[i] = valp[i+1] ;
    }

    /// Collect edges from the sorted edge_keys
    template<class T>
    void collect_edges(vector<pair<T,T> > &edgemap,
                       vector<T> &fileno,
                       const std::vector<Array<T,5> > &edge_keys) {
#ifdef DEBUG
      // In debug mode, check consistency of edge_keys size information
      for(size_t i=0;i<edge_keys.size();) {
        // search for cluster of edges
        size_t j = i+1 ;
        for(;j<edge_keys.size();++j)
          if((edge_keys[i][0] != edge_keys[j][0]) ||
             (edge_keys[i][1] != edge_keys[j][1]))
            break ;
        int nedges = j-i ;
        FATAL(int(j-i) != edge_keys[i][3]);
        i = j ;
      }
#endif
      // Count edges
      size_t edge_cnt = 0 ;
      for(size_t i=0;i<edge_keys.size();i+=edge_keys[i][3])  edge_cnt++ ;

      // Allocate edgemap and fileno
      vector<pair<T,T> > tmp(edge_cnt) ;
      edgemap.swap(tmp) ;
      vector<T> tmp2(edge_cnt) ;
      fileno.swap(tmp2) ;
      // Extract edge information from edge_keys
      edge_cnt = 0 ;
      for(size_t i=0;i<edge_keys.size();i+=edge_keys[i][3]) {
        edgemap[edge_cnt] = std::make_pair(edge_keys[i][0],edge_keys[i][1]) ;
        fileno[edge_cnt] = edge_keys[i][2] ;
        edge_cnt++ ;
      }
    }

    /// Change the edge_keys edge number annotation to use the numbering in the
    /// supplied edges set.
    ///
    /// -- changes edge_keys[2] to the edge number in global numbering
    ///
    template<class T> void renumber_edge_keys(std::vector<Array<T,5> > &edge_keys, entitySet edges) {
      auto edge_global = edges.begin() ;
      for(size_t i=0;i<edge_keys.size();i+=edge_keys[i][3],++edge_global) {
        for(int j=0;j<edge_keys[i][3];++j) 
          edge_keys[i+j][2] = *edge_global ;
      }
      // Check to make sure that edge set is the same size as the number of edges in
      // edge_keys
      FATAL(edge_global != edges.end()) ;
    }
  }

  /// Create edge datastructure.  Note, that we need this datastructure
  /// to have a consistent numbering regardless of the partition of the
  /// entities, so we have to make reference to the original file numbering
  /// supplied by the VOG file.
  ///
  /// The general strategy is to loop over faces and create a list of edges.
  /// This list will contain duplicate edges which will need to be removed.
  /// To remove the duplicates we first represent the edges with an ordered
  /// pair of the two defining nodes (numbered in the node file numbering)
  /// sorted so that the lowest node number comes first.  Sorting this list
  /// will produce a consistent ordering of edges where duplicate edges can
  /// be removed.  This will define the file numbering of the edges, but this
  /// will not be an efficient partition of edges.  So once the file numbering
  /// is determined, the edges will need to be sorted a second time based
  /// on the partitioning of nodes to get a proper global numbering of
  /// edges.  Finally we need to reconstruct the face2edge map which will be
  /// easier if we carry the face information that generated the edges through
  /// all of these sorting processes.
  void createEdgesPar(fact_db &facts) {
    // First get node domain
    store<vector3d<double> > pos ;
    pos = facts.get_variable("pos") ;
    entitySet ndom = pos.domain() ;
    int nkeyspace = pos.getDomainKeySpace() ;
    
    // Get face2node map
    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    entitySet fdom = face2node.domain() ;

    // Extract edges from face2node
    vector<Array<int,5> > edge_keys ;
    extract_edge_keys(edge_keys, fdom, ndom, nkeyspace,
                      face2node,facts) ;

    // Sort edges based on file numbering of edge nodes
    vector<Array<int,5> > edge_splits ;
    edge_sort(edge_keys,edge_splits,MPI_COMM_WORLD) ;

    // Find file numbering from edges, change layout of
    // edge_keys to use global node numbering
    edge_file_numbering(edge_keys,MPI_COMM_WORLD) ;
    
    // Now sort based on the global node numbering to improve locality
    get_splits(edge_splits,ndom,MPI_COMM_WORLD) ;
    edge_sort(edge_keys,edge_splits,MPI_COMM_WORLD) ;

    // Collect edge information
    vector<pair<int,int> > edge_list ;
    vector<int> edge_fileno ;
    collect_edges(edge_list,edge_fileno,edge_keys) ;

    // Allocate entities for new edges
    int num_edges = edge_list.size() ;
    int ek = facts.getKeyDomain("Edges") ;
    //    if(!useDomainKeySpaces) {
    //      ek = 0 ;
    //    }
    if(Loci::MPI_processes < 2)
      ek = 0 ;
    int fk = face2node.Rep()->getDomainKeySpace() ;

    entitySet edges = facts.get_distributed_alloc(num_edges,ek).first ;

    
    //create constraint edges
    constraint edges_tag;
    *edges_tag = edges;
    edges_tag.Rep()->setDomainKeySpace(ek) ;
    facts.create_fact("edges", edges_tag);
    // Copy edge nodes into a MapVec
    MapVec<2> edge ;
    // Note, set range keyspace also (when we have a Nodes keyspace)
    edge.Rep()->setDomainKeySpace(ek) ;
    edge.allocate(edges) ;
    int cnt = 0 ;
    debugout << "Edges = " << edges << endl ;
    // Pull edge2node data out of edge_list
    for(auto ei=edges.begin();ei!=edges.end();++ei,++cnt) {
      edge[*ei][0] = edge_list[cnt].first ;
      edge[*ei][1] = edge_list[cnt].second ;
    }

    // Update file numbering map
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    if(df != 0) {
      int cnt = 0 ;
      for(auto ei=edges.begin();ei!=edges.end();++ei,++cnt) {
        df->g2fv[ek][*ei] = edge_fileno[cnt] ; // Update global to file map
      }
    }

    // Add edge2node to fact database
    facts.create_fact("edge2node",edge) ;

    // Now we are creating the face2edge map.  For this we are going to sort the
    // edge_keys based on the global face number to return the data to the generating
    // face.  Along the way we will be communicating the final global numbering of the
    // created edges.  In this way we can create the face2edge map.

    
    // First change numbering in edge_keys from file numbering of edges to global numbering
    renumber_edge_keys(edge_keys,edges) ;

    // Swap the first entry in the key with the face global number
    // and sort the data so that all of the edge data is returned to
    // the face that generated it, only this time with the edge global
    // number
    for(size_t i=0;i<edge_keys.size();++i)
      std::swap(edge_keys[i][0],edge_keys[i][4]) ;
    // Sort based on new key to move edge_keys back to processor with generating face
    get_splits(edge_splits,fdom,MPI_COMM_WORLD) ;
    edge_sort(edge_keys,edge_splits,MPI_COMM_WORLD) ;

    // Now prepare face2edge map, it should have the same number of entries as face2node
    multiMap face2edge ;
    store<int> count ;
    count.allocate(fdom) ;
    for(auto ei = fdom.begin(); ei!=fdom.end();++ei)
      count[*ei] = face2node[*ei].size() ;
    face2edge.allocate(count) ;

    // Now we have to search over the edges generated by the face in edge_keys and find the
    // specific edge to be consistent with face2node
    auto fptr = fdom.begin() ;
    for(size_t i=0;i<edge_keys.size();) {
      // search for cluster edges tied to each face
      size_t j = i+1 ;
      for(;j<edge_keys.size();++j)
        if((edge_keys[i][0] != edge_keys[j][0]))
          break ;
      int fsz = j-i ;
      // We should have found a number of edges that matches the face2node size
      FATAL(fsz != face2node[*fptr].size());
      // Also our face should match the face given in the edge_keys
      FATAL(*fptr != edge_keys[i][0]) ;

      // Now scan over the candidate edges to find the matching edge
      for(int f=0;f<fsz;++f) {
        Entity n1 = face2node[*fptr][f] ;
        Entity n2 = face2node[*fptr][(f+1==fsz)?0:f+1] ;
        // e1 and e2 are the consistent edge node numbers
        Entity e1 = min(n1,n2) ; 
        Entity e2 = max(n1,n2) ;
        int k=0 ;
        for(k=0;k<fsz;++k) {
          // m1 and m2 are the consistent edge node numbers from the search set.
          int m1 = min(edge_keys[i+k][1],edge_keys[i+k][4]) ;
          int m2 = max(edge_keys[i+k][1],edge_keys[i+k][4]) ;
          // If we find it, grab the edge global number from the edge keys
          if(e1 == m1 && e2 == m2) {
            face2edge[*fptr][f] = edge_keys[i+k][2] ;
            break ;
          }
        }
        // Search failed to find corresponding edge, something wrong with data
        FATAL(k==fsz) ;
      }
      // advance to the next face
      fptr++ ;
      i = j ;
    }

    // Add face2edge to the fact database
    face2edge.Rep()->setDomainKeySpace(fk) ;
    MapRepP(face2edge.Rep())->setRangeKeySpace(ek) ;
    facts.create_fact("face2edge",face2edge) ;
  } // end of createEdgesPar

  void setupOverset(fact_db &facts) {
    storeRepP sp = facts.get_variable("componentGeometry") ;
    if(sp == 0) {
      return ;
    }
    sp = facts.get_variable("ci") ;
    if(sp == 0) {
      return ;
    }
    entitySet bfaces = sp->domain() ;
    sp = facts.get_variable("interface_BC") ;
    if(sp != 0) {
      constraint tmp ;
      tmp = sp ;
      bfaces -= *tmp ;
    }
    sp = facts.get_variable("symmetry_BC") ;
    if(sp != 0) {
      constraint tmp ;
      tmp = sp ;
      bfaces -= *tmp ;
    }
    sp = facts.get_variable("face2node") ;

    MapRepP mp = MapRepP(sp->getRep()) ;
    entitySet nodeSet = mp->image(bfaces) ;

    store<vector3d<double> > pos ;

    pos = facts.get_variable("pos") ;
    entitySet dom = pos.domain() ;
    vector<entitySet> posptn = all_collect_vectors(dom,MPI_COMM_WORLD) ;
    entitySet surfNodes = dist_collect_entitySet(nodeSet,posptn) ;

    // Now get volume tags
    variableSet vars = facts.get_extensional_facts() ;
    std::map<string,entitySet> volMap ;

    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      if(variable(*vi).get_arg_list().size() > 0 &&
         variable(*vi).get_info().name == "volumeTag") {
        param<string> vname(facts.get_variable(*vi)) ;

	      std::ostringstream vn ;
        vn << *vi ;
        string name = vn.str() ;
        volMap[name] = vname.domain() ;
      }
    }
    Map cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    entitySet domf = cl.domain()+cr.domain() ;
    entitySet domc = cl.image(domf)+cr.image(domf) ;

    // If no volume tags (Weird), then make default tag.
    if(volMap.begin() == volMap.end()) {
      volMap[string("Main")] = ~EMPTY ;
    }
    vector<entitySet> volSets ;
    std::map<string,entitySet>::const_iterator mi ;
    for(mi=volMap.begin();mi!=volMap.end();++mi) {
      int ckeyspace = cl.getRangeKeySpace() ;
      std::vector<entitySet> cptn = facts.get_init_ptn(ckeyspace) ;
      entitySet volgather = distribute_entitySet(mi->second,cptn) ;
      volgather = dist_expand_entitySet(volgather,domc,cptn) ;
      volSets.push_back(volgather) ;
    }

    // Now get face associations with volumes
    vector<entitySet> facesets ;
    int sz = volSets.size() ;
    for(int i=0;i<sz;++i) {
      entitySet faces = (cr.preimage(volSets[i]).first +
                         cl.preimage(volSets[i]).first) ;
      facesets.push_back(faces) ;
    }

    // Now get node associations with volumes

    vector<entitySet> nodesets ;
    for(int i=0;i<sz;++i) {
      entitySet nodes = mp->image(facesets[i]) ;
      int nkeyspace = pos.getDomainKeySpace() ;
      std::vector<entitySet> nptn = facts.get_init_ptn(nkeyspace) ;

      nodes = distribute_entitySet(nodes,nptn) ;
      nodesets.push_back(nodes) ;
    }


    Map min_node2surf_loc ;
    min_node2surf_loc.allocate(pos.domain()) ;

    entitySet excludeSet ;

    for(int i=0;i<sz;++i) {
      entitySet nodeSet = nodesets[i] & pos.domain() ;
      entitySet nodeSetsurf = nodesets[i] & surfNodes ;
      if(!GLOBAL_OR(nodeSetsurf.size()!=0)) {
        excludeSet += nodeSet ;
        continue ;
      }
      vector<Loci::kdTree::coord3d> bcnodes_pts(nodeSetsurf.size()) ;
      vector<int> bcnodes_ids(nodeSetsurf.size()) ;

      int cnt = 0 ;
      FORALL(nodeSetsurf,nd) {
        bcnodes_pts[cnt][0] = pos[nd].x ;
        bcnodes_pts[cnt][1] = pos[nd].y ;
        bcnodes_pts[cnt][2] = pos[nd].z ;
        bcnodes_ids[cnt] = nd ;
        cnt++ ;
      } ENDFORALL ;

      vector<Loci::kdTree::coord3d> node_pts(nodeSet.size()) ;
      vector<int> closest(nodeSet.size(),-1) ;
      cnt = 0 ;
      FORALL(nodeSet,nd) {
        node_pts[cnt] = pos[nd] ;
        cnt++ ;
      } ENDFORALL ;

      Loci::parallelNearestNeighbors(bcnodes_pts,bcnodes_ids,node_pts,closest,
                                     MPI_COMM_WORLD) ;

      cnt = 0 ;
      FORALL(nodeSet,nd) {
        min_node2surf_loc[nd] = closest[cnt] ;
        cnt++ ;
      } ENDFORALL ;
    }

    Map min_node2surf ;

    if(excludeSet == EMPTY) {
      min_node2surf.setRep(min_node2surf_loc.Rep()) ;
    } else {
      entitySet dom = pos.domain()-excludeSet ;
      min_node2surf.allocate(dom) ;
      FORALL(dom,nd) {
        min_node2surf[nd] = min_node2surf_loc[nd] ;
      } ENDFORALL ;
    }

    facts.create_fact("node2surf",min_node2surf) ;
  }

}
