//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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

#define vect3d vector3d<double>

namespace Loci{
  extern int getKeyDomain(entitySet dom, fact_db::distribute_infoP dist, MPI_Comm comm) ;
  extern  bool useDomainKeySpaces  ;
  // Convert container from local numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm);
  // Convert container from local numbering to output file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // fact_db pointer  (facts)
  // MPI Communicator
  storeRepP Local2FileOrder_output(storeRepP sp, entitySet dom,
                                   fact_db& facts, MPI_Comm comm);
  
  
  //map from local numbering to input file numbering
  //assume the union of nodes on all processors will be either all the nodes,
  //all the faces, or all the cells. i.e., one interval in file numbering
  storeRepP get_node_remap(fact_db &facts,entitySet nodes,
                           int keyspace) {
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

  //return a map from local numbering to output file numbering
  //the file number starts with 1
  //input: nodes : the local store domain need to be output,
  //assume nodes will be written out in the global ordering.  

  storeRepP get_output_node_remap(fact_db &facts,entitySet nodes) {
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
      debugout<<"ERROR: l2g.domain is smaller than  dom in get_output_node_remap " << endl;
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
    for(int i=0;i<p;++i)
      recvbufsz+= recvsz[i] ;
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
    // Now sort a copy of the recieved data so the local data is in order
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
    for(int i=0;i<ssz;++i)
      filemap[recvdatac[i]] = i+roff ;
    // Now send back the file numbering
    vector<int> filenosend(recvbufsz,-1) ;
    for(int i=0;i<recvbufsz;++i) 
      filenosend[i] = filemap[recvdata[i]] ;
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
  
  int classify_cell(Entity *faces,int nfaces,const_multiMap &face2node) {
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
      } else if(count == 4)
        num_quads++ ;
      else
        num_others++ ;
    }
    bool prism_test = false ;

    if((num_triangles == 2) && (num_quads == 3) && (num_others == 0)) {
      prism_test = true ;
      for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
          if(triangle_nodes[i][0] == triangle_nodes[j][1])
            prism_test = false ;
    }

    bool hex_test = false ;
    if( (num_triangles == 0) && (num_quads == 6) && (num_others == 0)) {
      const Entity ef = faces[0] ;
      int count = 0 ;
      for(int fj = 1;fj<nfaces;++fj) {
        Entity ej = faces[fj] ;
        bool find = false ;
        for(int i=0;i<4;++i)
          for(int j=0;j<4;++j)
            if(face2node[ef][i] == face2node[ej][j])
              find = true ;
        if(find)
          count++ ;
      }
      if(count == 4)
        hex_test = true ;
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

  void fillTet(Array<int,4> &tet, Array<int,3> *tri_faces) {
    tet[0] = tri_faces[0][2] ;
    tet[1] = tri_faces[0][1] ;
    tet[2] = tri_faces[0][0] ;
    for(int i=0;i<3;++i) {
      tet[3] = tri_faces[1][i] ;
      if(tet[3] != tet[0] && tet[3] != tet[1] && tet[3] != tet[2])
        return ;
    }
    cerr << "unable to form valid tet!" << endl ;
  }

  void fillHex(Array<int,8> &hex, Array<int,4> *quad_faces) {
    int quad_id[6] ;
    for(int i=0;i<6;++i) 
      quad_id[i] = i ;
    bool degenerate = quad_faces[quad_id[0]][0] == quad_faces[quad_id[0]][3];
    for(int j=0;j<3;++j) 
      if(quad_faces[quad_id[0]][j] == quad_faces[quad_id[0]][j+1])
        degenerate = true ;
    if(degenerate) {
      for(int i=1;i<6;++i) {
        degenerate = quad_faces[quad_id[i]][0] == quad_faces[quad_id[i]][3];
        for(int j=0;j<3;++j) 
          if(quad_faces[quad_id[i]][j] == quad_faces[quad_id[i]][j+1])
            degenerate = true ;
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


  string MPIConcatStrings(string input, MPI_Comm comm) {
    int sz = input.size()+1 ;
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    int *sizes = new int[p] ;
    MPI_Allgather(&sz,1,MPI_INT,sizes,1,MPI_INT,comm) ;
    int tot = 0 ;
    for(int i=0;i<p;++i)
      tot += sizes[i] ;
    char *buf = new char[tot] ;
    int *displ = new int[p] ;
    displ[0] = 0 ;
    for(int i=1;i<p;++i)
      displ[i] = displ[i-1]+sizes[i-1] ;
    char *ibuf = new char[sz] ;
    strcpy(ibuf,input.c_str()) ;
    MPI_Allgatherv(ibuf,sz,MPI_CHAR,
                   buf, sizes, displ,MPI_CHAR,comm) ;
    string retval ;
    for(int i=0;i<p;++i)
      retval += string(&(buf[displ[i]])) ;

    delete[] sizes ;
    delete[] buf ;
    delete[] ibuf ;
    delete[] displ ;
    
    return retval ;
  }
  //node_id range is 1~numNodes
  //face_id range is numNodes~numNodes+numFaces
  //cell_id range is 1~numCells
  //Reason: in FVMstuff.cc, parallelWriteGridTopo(),
  //get_node_remap is used to get the file number of nodes and cells, which set the start index to 1
  //however, g2f(l2g(ff)) is used to get the file number of faces.
  //Since boundary variables also use the face_id, to avoid changing too many places, we keep the face_id this way.
  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells,
                                 fact_db &facts) {
    const_multiMap upper(upperRep),lower(lowerRep),
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
            for(int j=0;j<fsz;++j)
              generalCellNodes.push_back(node_remap[face2node[fc][fsz-j-1]]) ;
          } else { 
            for(int j=0;j<fsz;++j)
              generalCellNodes.push_back(node_remap[face2node[fc][j]]) ;
          }
        }
      }
    } ENDFORALL ;


    // write grid topology file
    hid_t file_id = 0, group_id = 0 ;
    file_id=writeVOGOpen(filename) ;
    if(use_parallel_io ||MPI_rank == 0 ) {
      group_id = H5Gcreate(file_id,"elements",
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

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
      group_id = H5Gcreate(file_id,"boundaries",
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
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
        if(init_ptn[MPI_rank].inSet(l2g[bb]))
          local_boundaries += bb ;
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
      while(i < bnames.size() && bnames[i] != '"')
        tmp += bnames[i++] ;
      bnamelist.push_back(tmp) ;
    }


    // Identify Boundary faces to be written by this processor
    entitySet  fset = (MapRepP(boundary_mapRep)->image(localCells)+
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
      //      debugout << "in write gridTopology, face key domain is " << kd
      //	       << endl ;
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
        bc_id = H5Gcreate(group_id,current_bc.c_str(),
			  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
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
          if(ref[fc] == id)
            bfaces+= fc ;
        }ENDFORALL ;
      }
      
      int ntria=0, nquad=0, nsided =0;
      FORALL(bfaces,fc) {
        if(face2node[fc].size() == 3)
          ntria++ ;
        else if(face2node[fc].size() == 4)
          nquad++ ;
        else
          nsided++ ;
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
          for(int i=0;i<nsizes[ng];++i)
            nsidenodes.push_back(node_remap[face2node[fc][i]]) ;
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
  
  //collect all boundary names
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
        if(init_ptn[MPI_rank].inSet(l2g[bb]))
          local_boundaries += bb ;
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
      while(i < bnames.size() && bnames[i] != '"')
        tmp += bnames[i++] ;
      bnamelist.push_back(tmp) ;
    }
    std::sort(bnamelist.begin(), bnamelist.end()); 
    return bnamelist;
  }
  
  //mkdir /output/$bc_name
  void get_bc_directory(string bc_name){
    string directory_name = "output/"+bc_name;
    struct stat statbuf ;
    int fid = open(directory_name.c_str(), O_RDONLY) ;
    if(fid < 0) {
      mkdir(directory_name.c_str(),0755) ;
    } else {
      fstat(fid,&statbuf) ;
      if(!S_ISDIR(statbuf.st_mode)) {
        cerr << "file " << directory_name <<" should be a directory!, rename it and start again."
             << endl ;
        Loci::Abort() ;
      }
      close(fid) ;
    }
  }
  //open /output/$bc_name/$file_name
  hid_t open_boundary_file(string bc_name,
                           string file_name
                           ){
    // open the file
    hid_t file_id = 0;
    if(MPI_rank == 0) {
      get_bc_directory(bc_name);
    }
    string dirname = "output/"+bc_name+"/";
    string filename = dirname+file_name;
    file_id=writeVOGOpen(filename.c_str()) ;
    return file_id;
  }
  
 
  
  //get boundary faces that belong to a boundary surface current_bc
  entitySet get_boundary_faces(string current_bc,//boundary name
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset //all boundary faces 
                               ){
    
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
        if(ref[fc] == id)
          bfaces+= fc ;
      }ENDFORALL ;
    }
    return bfaces;
  }

  //get boundary nodes that belong to a boundary surface current_bc
  entitySet get_boundary_nodes(string current_bc,//boundary name
                               storeRepP face2nodeRep,
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset, //all boundary faces 
                               fact_db &facts ){
    entitySet bfaces = get_boundary_faces(current_bc, refRep, bnamesRep, fset);
    //get containers
    const_multiMap face2node(face2nodeRep) ;
    //get all the local nodes on the boundary 
    entitySet nodes_local = MapRepP(face2node)->image(bfaces);
    // get the nodes that belong to this processor
    if(MPI_processes > 1) {
      entitySet dom = nodes_local;
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      entitySet  my_entities = dist->my_entities ;
      nodes_local = my_entities & nodes_local ;
    }
    return nodes_local;
  }
       
  void writeBoundaryTopo(hid_t file_id, //file_id of this boudnary surface
                         storeRepP face2nodeRep,
                         entitySet bfaces,//boundary faces define this surface 
                         fact_db &facts ){ 
#ifdef VERBOSE
    debugout << "writing out boundary surface topology" << endl ;
    stopWatch s ;
    s.start() ;
#endif

    const_multiMap face2node(face2nodeRep) ;

    //collect the local boundary nodes belong to this boundary
    entitySet nodes_local = (MapRepP(face2node)->image(bfaces));
          
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
    int sz = bfaces.size();
    int off = 0 ;
    MPI_Scan(&sz,&off,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    off -= sz ;
    FORALL(bfaces,fc) {
      faceorder[fc] = off ;
      off++ ;
    } ENDFORALL ;
    //get output vectors
    int ntria=0, nquad=0, nsided =0;
    FORALL(bfaces,fc) {
      if(face2node[fc].size() == 3)
        ntria++ ;
      else if(face2node[fc].size() == 4)
        nquad++ ;
      else
        nsided++ ;
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
        for(int i=0;i<nsizes[ng];++i)
          nsidenodes.push_back(node_remap[face2node[fc][i]]) ;
        
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
      
  //this function find the index of an inner edge.
  //an inner edge is the edge from facecenter to one of its nodes during face triangulation
  //if the edge already exists, the index is returned
  //otherwise, the edge is created and registered in inner_edges and edgeIndex
  //the return value is the index of the edge in inner_edges, starts with 1
  int create_inner_edge(int face, //face id
                        int node, //node id
                        std::map<pair<int, int>, int>& edgeIndex, //map of inner edges to its index in inner_edges, the index starts with 1 
                        vector<pair<int, int> >& inner_edges) //the edges from facecenter and one of its node
  {
    pair<int, int> inner_edge = make_pair(face, node);
    //make sure no duplicate in inner_edges
    if(edgeIndex.find(inner_edge) == edgeIndex.end()){
      inner_edges.push_back(inner_edge);
      edgeIndex[inner_edge] = inner_edges.size();
    }
    return edgeIndex[inner_edge];
  }
  
  //If a face has more than 2 edges cut, to disambiguate the connection between the cut points,
  //the face will be triangulated and each tri face will be registered
  void disambiguateFace(Entity f, 
                        const_multiMap& face2node,
                        const_multiMap& face2edge,
                        const_MapVec<2>& edge2node,
                        const_store<double>& val,
                        entitySet& edgesCut, //the edges(node 2 node) cut
                        std::map<pair<int, int>, int >& edgeIndex,//map of inner edges(face 2 node) to its index in inner_edges, the index starts with 1 
                        vector<pair<int, int> >& intersects, //pair of edges, if the second number is negative, it is index to inner_edges
                        vector<pair<int, int> >& inner_edges //the edges from facecenter and one of its node 
                        ){
    
    int nNodes = face2node.num_elems(f);
    
    //get the value of face center
    double  newNode = 0.0 ;
    for (int i = 0; i < nNodes; ++i) {
      newNode += val[face2node[f][i]];
    }
    newNode /= double(nNodes);
    
    //triangulate the face
    //and register each tri face.
    int faceCut[2];
    for (int i = 0;i < face2edge.num_elems(f); ++i) {
      int cutsFound = 0;
      //check this edge
      Entity edge = face2edge[f][i];
      if (signbit(val[edge2node[edge][0]] * val[edge2node[edge][1]])) {
        edgesCut += edge;
        faceCut[cutsFound] = edge;
        cutsFound++; 
      }
      //check face center to first node
      if (signbit(newNode * val[edge2node[edge][0]])) {
        int node = edge2node[edge][0];
        int edgeId = create_inner_edge(f, node, edgeIndex,inner_edges);
        faceCut[cutsFound]  = -edgeId;
        cutsFound++;                     
      }
  
      //check face center to second node
      if (signbit(newNode * val[edge2node[edge][1]])) {
        if(cutsFound >1){
          cerr<<"ERROR: tri face has more than two cuts" << endl;
          exit(-1);
        }
        int node = edge2node[edge][1];
        int edgeId = create_inner_edge(f, node, edgeIndex,inner_edges);
        faceCut[cutsFound]  = -edgeId;
        cutsFound++;                     
      }
      if (cutsFound == 2) intersects.push_back(make_pair(faceCut[0], faceCut[1]));
      //else ????holes happen?
    }
  }

  // Find edges cut in the face and register them into intersects
  bool registerFace(Entity f,
                    const_multiMap& face2node,
                    const_multiMap& face2edge,
                    const_MapVec<2>& edge2node,
                    const_store<double>& val,
                    entitySet& edgesCut,//the edges(node 2 node) cut
                    std::map<pair<int, int>, int >& edgeIndex,//map of inner edges(face center 2 node) to its index in inner_edges, the index starts with 1
                    vector<pair<int, int> >& intersects, //pair of edges, if the value is negative, it is index to inner_edges
                    vector<pair<int, int> >& inner_edges//pair<face, node>, the edges from facecenter and one of the face node 
                    ){
    
    int faceCut[2];
    int cutsFound = 0;
    //check each edge
    for(int i = 0; i < face2edge.num_elems(f); ++i) {
      Entity edge = face2edge[f][i];
      //if it is cut
      if (signbit(val[edge2node[edge][0]] * val[edge2node[edge][1]])) {
        if (cutsFound < 2) {
          //store it in faceCut
          faceCut[cutsFound] = edge;
          edgesCut += edge;
        }
        cutsFound++; 
      }
    }
    //no ambiguation, register faceCut in inertsection
    //otherwise, disambiguate the face
    if (cutsFound == 2)
      intersects.push_back(make_pair(faceCut[0], faceCut[1]));
    else if (cutsFound > 2)
      disambiguateFace( f,
                        face2node,
                        face2edge,
                        edge2node,
                        val,
                        edgesCut,
                        edgeIndex,
                        intersects,
                        inner_edges
                        );
    return (cutsFound > 0);
  }
    
  //after all the faces of a cell is registered,
  //check the intersects to form the faces in cutting plane
  void checkLoop( const vector<pair<int, int> >& intersects,//pair of edges, if the value number is negative, it is index to inner_edges
                  vector<vector<int > >& faceLoops,//the loops(faces in cutting plane) formed 
                  int start,//the start index to intersects 
                  int end) { //the end index to intersects
    bool loopIsGood = true;
    vector<int > loop;
  
    list<int> edges;
    for (int i = start+1; i < end; i++)
      edges.push_back(i);

    int firstNode, nextNode;
    firstNode = intersects[start].first;
    nextNode = intersects[start].second;
    loop.push_back(firstNode);
    loop.push_back(nextNode);

    list<int>::iterator iter;
    while (!edges.empty() && loopIsGood) {
      for (iter = edges.begin(); iter != edges.end(); iter++) {
        if (intersects[*iter].first == nextNode) {
          nextNode = intersects[*iter].second;
          if(firstNode != nextNode)loop.push_back(nextNode);
          edges.erase(iter);
          break;
        } else if (intersects[*iter].second == nextNode) {
          nextNode = intersects[*iter].first;
          if(firstNode != nextNode) loop.push_back(nextNode);
          edges.erase(iter);
          break;
        }
      }
      //can not find the next edge
      if (iter == edges.end())
        loopIsGood = false;

      //a loop is formed
      if (firstNode == nextNode){
        if(!edges.empty()) {
          faceLoops.push_back(loop);
          loop.clear();
        
          firstNode = intersects[edges.front()].first;
          nextNode = intersects[edges.front()].second;
          loop.push_back(firstNode);
          loop.push_back(nextNode);
          edges.erase(edges.begin());
        }else{
          faceLoops.push_back(loop);
        }
      }
    }
      
    if (firstNode != nextNode){
      debugout << "ERROR: ** Problem cell:  ** (failed loop test)" << endl;
    }
  }
  //after the loops are formed, this function will remove all negative value in loops.
  vector<vector<int > > remove_inner_edges(const vector<vector<int > >& faceLoops)//the loops(faces in cutting plane) formed 
  {
    vector<vector<int > >  new_faceLoops;
    for(unsigned int i = 0; i < faceLoops.size(); i++){
      vector<int> loop;
      for(unsigned int j = 0; j < faceLoops[i].size(); j++){
        if(faceLoops[i][j] >=0) loop.push_back(faceLoops[i][j]);
      }
      if(loop.size() >0)new_faceLoops.push_back(loop);
    }
    return new_faceLoops;
  }
    
  //return the cutting position of an edge 
  double get_edge_weight(Entity e,
                         const_MapVec<2>& edge2node,
                         const_store<double>& val){
    
    double a = val[edge2node[e][0]];
    double b = val[edge2node[e][1]];
    return ((b)/(b - a));
  }
 

  //generate Cutplane
  CutPlane getCutPlane(storeRepP upperRep,
                       storeRepP lowerRep,
                       storeRepP boundary_mapRep,
                       storeRepP face2nodeRep,
                       storeRepP face2edgeRep,
                       storeRepP edge2nodeRep,
                       storeRepP valRep,
                       entitySet localCells,//all geom_cells
                       fact_db &facts) {
    
    const_multiMap upper(upperRep),lower(lowerRep),
      boundary_map(boundary_mapRep),face2node(face2nodeRep), face2edge(face2edgeRep) ;
    
    const_MapVec<2> edge2node(edge2nodeRep);
    const_store<double> val(valRep);

    // data structure:
    entitySet edgesCut;
    store<double> edgesWeight; //the weight for interpoplation for each edgesCut, allocated on edgesCut
    vector<pair<int, int> > inner_edges; //the inner edges(facecenter to one of its nodes) cut, the values stored are pair<local_faceid, noderank>  
    vector<vector<int > > faceLoops;  //loops formed, the values stored are edge ids, which is either local edge entity or index to inner_edges
   
   
    //extra data structure
    vector<pair<int, int> > intersects; //the pairs of edges  or inner edges cut
    std::map<pair<int, int>, int > edgeIndex; //map of inner edges to their indexes in inner_edges
    // for each cell
    FORALL(localCells,cc) {
      bool  isCut = false;
      int  intersectStart = intersects.size();
      
      //collect all its faces
      entitySet faces;
      for(int i=0;i<upper[cc].size();++i)faces += upper[cc][i];
      for(int i=0;i<lower[cc].size();++i)faces += lower[cc][i];
      for(int i=0;i<boundary_map[cc].size();++i)faces += boundary_map[cc][i];
      
      //for each face, find the edges that is cut and register it
      for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ei++){
        if (registerFace(*ei,
                         face2node,
                         face2edge,
                         edge2node,
                         val,
                         edgesCut,
                         edgeIndex,
                         intersects, 
                         inner_edges 
                         )) isCut = true;
      }

      //check the loops formed by this cell and register the loop
      if(isCut){
        checkLoop(intersects, faceLoops, intersectStart, intersects.size());
        intersectStart = intersects.size();
      }
    }ENDFORALL;
    //remove the inner edges cut from loops
    vector<vector<int > > new_faceLoops = remove_inner_edges(faceLoops);
    
    //compute the cutting postions of edges 
    edgesWeight.allocate(edgesCut);
    FORALL(edgesCut, e){
      double t = get_edge_weight(e, edge2node, val);
      edgesWeight[e] = t;
    }ENDFORALL;
    
   
    return CutPlane(edgesWeight.Rep(),
                    new_faceLoops);
                     
   
  }

  
  void writeCutPlaneTopo(hid_t bc_id,
                         const CutPlane& cp,
                         fact_db &facts){ 
#ifdef VERBOSE
    debugout << "write cutPlaneTopology" << endl ;
    stopWatch s ;
    s.start() ;
#endif
    Map node_remap;//map from local numbering to output file numbering for edge entities
    entitySet edgesCut = (cp.edgesWeight)->domain();
    node_remap = get_output_node_remap(facts, edgesCut);
#ifdef VERBOSE
    debugout << "time to get remap = " << s.stop() << endl ;
    s.start() ;
#endif
   
    //write out the sizes of faceLoops  
    int num_faces = cp.faceLoops.size();
    vector<int> nsizes(num_faces);
    for(int i = 0; i < num_faces; i++){
      nsizes[i] = cp.faceLoops[i].size();
    }
    
    
    writeUnorderedVector(bc_id,"nside_sizes",nsizes) ;
    
    //write out face nodes
    vector<int> nsidenodes ;
    for(int i = 0; i < num_faces; i++){
      for(int j = 0; j < nsizes[i]; j++){
        int node = cp.faceLoops[i][j];
        if(node >=0){
          nsidenodes.push_back(node_remap[node]); //edge nodes
        }
      }
    }
    
   
    writeUnorderedVector(bc_id,"nside_nodes",nsidenodes) ;

#ifdef VERBOSE
    debugout << "time to write cut plane topology=" << s.stop() << endl ;
#endif
  }

  namespace {
    
    struct BCinfo {
      std::string name ;
      int key ;
      entitySet apply_set ;
      options_list bc_options ;
      BCinfo() {}
      BCinfo(const std::string &n, int k,
             const entitySet &a, const options_list &o) :
        name(n),key(k),apply_set(a),bc_options(o) {}
      
    } ; 
  }


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
      if(GLOBAL_OR(tmp_out)) 
        fill_clone(sp, out_of_dom, init_ptn) ;
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
      
      for(size_t i=0;i<p1.size();++i)
        pmap[p1id[i]] = p1closest[i] ;
      for(size_t i=0;i<p2.size();++i)
        pmap[p2id[i]] = p2closest[i] ;

      // Check to make sure connection is one-to-one
      dstore<int> check ;
      for(size_t i=0;i<p2.size();++i)
        check[p2id[i]] = p2closest[i] ;

      entitySet p1map = create_entitySet(p1closest.begin(),p1closest.end()) ;
      storeRepP sp = check.Rep() ;
      std::vector<entitySet> init_ptn ;
      if(facts.is_distributed_start()) {
        init_ptn = facts.get_init_ptn(fk) ; 
        fill_clone(sp, p1map, init_ptn) ;
      }
      bool periodic_problem = false ;
      for(size_t i=0;i<p1id.size();++i) 
        if(check[pmap[p1id[i]]] != p1id[i])
          periodic_problem = true ;
      if(GLOBAL_OR(periodic_problem)) {
        if(MPI_rank == 0) {
          cerr << "Periodic boundary did not connect properly, is boundary point matched?" << endl ;
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

  void setupBoundaryConditions(fact_db &facts) {
    list<BCinfo> BCinfo_list ;
    std::map<std::string,entitySet> BCsets ;
    
    /*Boundary Conditions*/
    entitySet periodic ;
    constraint periodic_faces;
    constraint no_symmetry_BC ;

    entitySet symmetry ;

    storeRepP tmp = facts.get_variable("boundary_names") ;
    if(tmp == 0) 
      throw(StringError("boundary_names not found in setupBoundaryConditions! Grid file read?")) ;
      
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
    if(tmp == 0)
      throw(StringError("boundary_conditions not found in setupBoundaryConditions! Is vars file read?")) ;
    bc_info = tmp ;
    
    param<double> Lref ;
    *Lref = 1.0 ;
    storeRepP p = facts.get_variable("Lref") ;
    if(p != 0)
      Lref = p ;
    
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
      
      option_value_type vt =
        bc_info->getOptionValueType(bname);
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
        cerr << "Boundary tag '" << bname << "' exists in VOG file without boundary_condition entry!  Please add boundary condition specification for this boundary face!" << endl ;
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
	  //          get_vect3dOption(ol,"center","m",pi.center,*Lref) ;
	  double Lr = realToDouble(*Lref) ;
	  ol.getOptionUnits("center","m",pi.center,Lr) ;
        }
        if(ol.optionExists("vector")) {
	  //          get_vect3d(ol,"vector",pi.v) ;
	  ol.getOptionUnits("vector","",pi.v) ;
          pi.v /=  norm(pi.v) ;
        }
        if(ol.optionExists("translate")) {
          // get_vect3dOption(ol,"translate","m",pi.translate,*Lref) ;
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
      std::map<std::string, entitySet>::const_iterator mi ;
      for(mi=BCsets.begin();mi!=BCsets.end();++mi) {
        if(GLOBAL_OR(mi->second.size())) {
          constraint bc_constraint ;
          bc_constraint = mi->second ;
          std::string constraint_name = mi->first + std::string("_BC") ;
	  bc_constraint.Rep()->setDomainKeySpace(fk) ;
          facts.create_fact(constraint_name,bc_constraint) ;
          if(MPI_processes == 1)
            std::cout << constraint_name << ' ' << mi->second << endl ;
          else if(MPI_rank == 0)
            std::cout << "setting boundary condition " << constraint_name << endl ;
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
          for(size_t j=i+1;j<periodic_data.size();++j)
            if(periodic_data[i].name == periodic_data[j].name) {
              if(p2.name != ",") {
                cerr << "periodic name appears more than two times!" ;
                Abort() ;
              }
              p2 = periodic_data[j] ;
              periodic_data[j].processed = true ;
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
          if(p1.master)
            periodic_list.push_back(std::make_pair(p1,p2)) ;
          else
            periodic_list.push_back(std::make_pair(p2,p1)) ;
        
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

  // Create Cell Stencil Plan:

  // Get relations from face2cell and face2node
  // Join to get node2cell
  // join node2cell & node2cell to get cell2cell (through shared nodes)
  // remove self2self references
  // Convert resulting structure to multiMap

  void getFaceCenter(fact_db &facts, dstore<vector3d<double> > &fcenter, dstore<double> &area, dstore<vector3d<double> > &normal) {
    // Compute face centers
    store<vector3d<double> > pos ;
    pos = facts.get_variable("pos") ;

    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    entitySet fdom = face2node.domain() ;
    entitySet f2n_image = Loci::MapRepP(face2node.Rep())->image(fdom) ;

    entitySet out_of_dom = f2n_image - pos.domain() ;

    dstore<vector3d<double> > tmp_pos ;
    FORALL(pos.domain(), pi) {
      tmp_pos[pi] = pos[pi] ;
    } ENDFORALL ;
    Loci::storeRepP sp = tmp_pos.Rep() ;
    int tmp_out = out_of_dom.size() ;
    if(facts.is_distributed_start()) {
      int nkeyspace = pos.getDomainKeySpace() ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn(nkeyspace) ;
      if(GLOBAL_OR(tmp_out)) 
	fill_clone(sp, out_of_dom, init_ptn) ;
    }

    param<std::string> centroid ;
    centroid = facts.get_variable("centroid") ;
    bool exact = false ;
    if(*centroid == "exact") 
      exact = true ;

    FORALL(fdom,fc) {
      int fsz = face2node[fc].size() ;
      vector3d<double>  center = vector3d<double> (0,0,0) ;
      double wsum = 0 ;
	
      for(int i=1;i<fsz;++i) {
	double len = norm(tmp_pos[face2node[fc][i-1]]-tmp_pos[face2node[fc][i]]) ;
	vector3d<double>  eloc = 0.5*(tmp_pos[face2node[fc][i-1]]+tmp_pos[face2node[fc][i]]) ;
	center += len*eloc ;
	wsum += len ;
      }
      double lenr = norm(tmp_pos[face2node[fc][0]]-tmp_pos[face2node[fc][fsz-1]]) ;
      vector3d<double>  elocr = 0.5*(tmp_pos[face2node[fc][0]]+tmp_pos[face2node[fc][fsz-1]]) ;
      center += lenr*elocr ;
      wsum += lenr ;
      center *= 1./wsum ;
      fcenter[fc] = center ;
      if(exact) {	
	// Iterate to find exact centroid that is on the face
	vector3d<double>  tmpcenter = center ;
	const int NITER=4 ;
	for(int iter=0;iter<NITER;++iter) {
    
	  // compute centroid using triangles formed by wireframe centroid
	  vector3d<double>  centroidsum(0.0,0.0,0.0) ;
	  double facearea = 0 ;
	  for(int i=0;i<fsz;++i) {
	    int n1 = i ;
	    int n2 = (i+1==fsz)?0:i+1 ;
	    vector3d<double>  p1 = tmp_pos[face2node[fc][n1]] ;
	    vector3d<double>  p2 = tmp_pos[face2node[fc][n2]] ;
	    
	    const vector3d<double>  t_centroid = (p1 + p2 + tmpcenter)/3.0 ;
	    const double t_area = 0.5*norm(cross(p1-tmpcenter,p2-tmpcenter)) ;
	    centroidsum += t_area*t_centroid ;
	    facearea += t_area ;
	  }
	  tmpcenter = centroidsum/facearea ;
	}
	fcenter[fc] = tmpcenter ;
      }
      vector3d<double> facearea = vector3d<double>(0,0,0) ;
      vector3d<double>  tmpcenter = fcenter[fc] ;
      for(int i=0;i<fsz;++i) {
	int n1 = i ;
	int n2 = (i+1==fsz)?0:i+1 ;
	vector3d<double>  p1 = tmp_pos[face2node[fc][n1]] ;
	vector3d<double>  p2 = tmp_pos[face2node[fc][n2]] ;
	
	facearea += cross(p1-tmpcenter,p2-tmpcenter) ;
      }
      double nfacearea = norm(facearea) ;
      area[fc] = 0.5*nfacearea ;
      normal[fc] = (1./nfacearea)*facearea ;
    } ENDFORALL ;
  }

  void getCellCenter(fact_db &facts, dstore<vector3d<double> > &fcenter,
		     dstore<vector3d<double> > &fnormal,
		     dstore<vector3d<double> > &ccenter) {

    dstore<double> area ;
    getFaceCenter(facts,fcenter,area,fnormal) ;
    multiMap upper,lower,boundary_map ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;
    entitySet cells = upper.domain()&lower.domain()&boundary_map.domain() ;
    entitySet faceimage  ;
    faceimage += Loci::MapRepP(upper.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(lower.Rep())->image(cells) ;
    faceimage += Loci::MapRepP(boundary_map.Rep())->image(cells) ;
    entitySet out_of_dom = faceimage - fcenter.domain() ;
    int tmp_out = out_of_dom.size() ;
    if(facts.is_distributed_start()) {
      
      int fkeyspace = upper.getRangeKeySpace() ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn(fkeyspace) ;
      if(GLOBAL_OR(tmp_out)) {
	Loci::storeRepP sp = fcenter.Rep() ;
	fill_clone(sp, out_of_dom, init_ptn) ;
	sp = area.Rep() ;
	fill_clone(sp, out_of_dom, init_ptn) ;
	sp = fnormal.Rep() ;
	fill_clone(sp, out_of_dom, init_ptn) ;
      }
    }
    dstore<vector3d<double> > wccenter ;
    // compute wireframe centroid
    FORALL(cells,cc) {
      vector3d<double>  csum = vector3d<double> (0,0,0) ;
      double wsum = 0 ;
      int lsz = lower[cc].size() ;
      for(int i=0;i<lsz;++i) {
	int f = lower[cc][i] ;
	double w = area[f] ;
	vector3d<double>  v = fcenter[f] ;
	csum += w*v ;
	wsum += w ;
      }
      int usz = upper[cc].size() ;
      for(int i=0;i<usz;++i) {
	int f = upper[cc][i] ;
	double w = area[f] ;
	vector3d<double>  v = fcenter[f] ;
	csum += w*v ;
	wsum += w ;
      }
      int bsz = boundary_map[cc].size() ;
      for(int i=0;i<bsz;++i) {
	int f = boundary_map[cc][i] ;
	double w = area[f] ;
	vector3d<double>  v = fcenter[f] ;
	csum += w*v ;
	wsum += w ;
      }
      csum *= 1./wsum ;
      wccenter[cc] = csum ;
      ccenter[cc] = csum ;
    } ENDFORALL ;
    param<std::string> centroid ;
    centroid = facts.get_variable("centroid") ;
    if(*centroid == "exact")  {
      // Here we add code to compute exact centroid.
      // face2node needs to be expanded to do this...
    }
    
  }

  void create_cell_stencil_full(fact_db &facts) {
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
    dist_expand_entitySet(geom_cells,cellmask,ptn) ;
    Loci::protoMap f2cell ;

    // Get mapping from face to geometric cells
    Loci::addToProtoMap(cl,f2cell) ;
    FORALL(faces,fc) {
      if(geom_cells.inSet(cr[fc]))
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

    //
    // Create cell stencil map from protoMap
    multiMap cellStencil ;
    distributed_inverseMap(cellStencil,c2c,geom_cells,geom_cells,ptn) ;
    // Put in fact database
    facts.create_fact("cellStencil",cellStencil) ;
  }
  
  void create_cell_stencil(fact_db & facts) {
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
    // loop before, so gather clone cells
    int gkeyspace = geom_cells_c.getDomainKeySpace() ;
    std::vector<entitySet> ptn = facts.get_init_ptn(gkeyspace) ;
    geom_cells = distribute_entitySet(geom_cells,ptn) ;
    dist_expand_entitySet(geom_cells,cellmask,ptn) ;
    
    Loci::protoMap f2cell ;

    // Get mapping from face to geometric cells
    Loci::addToProtoMap(cl,f2cell) ;
    FORALL(faces,fc) {
      if(geom_cells.inSet(cr[fc]))
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

    //
    // Create cell stencil map from protoMap
    multiMap cellStencil ;
    //    std::vector<entitySet> ptn = facts.get_init_ptn() ;
    distributed_inverseMap(cellStencil,c2c,geom_cells,geom_cells,ptn) ;

    // Now downselect cells
    dstore<vector3d<double> > fcenter ;
    dstore<vector3d<double> > fnormal ;
    dstore<vector3d<double> > ccenter ;
    getCellCenter(facts, fcenter, fnormal, ccenter) ;

    entitySet cells = cellStencil.domain() ;
    multiMap upper,lower,boundary_map ;
    upper = facts.get_variable("upper") ;
    lower = facts.get_variable("lower") ;
    boundary_map = facts.get_variable("boundary_map") ;

    entitySet cellImage = Loci::MapRepP(cellStencil.Rep())->image(cells) ;
    cellImage -= cells ;
    if(facts.is_distributed_start()) {
      int gkeyspace = geom_cells_c.getDomainKeySpace() ;
      std::vector<entitySet> ptn = facts.get_init_ptn(gkeyspace) ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn(gkeyspace) ;
      int tmp_out = cellImage.size() ;
      if(GLOBAL_OR(tmp_out)) {
	Loci::storeRepP sp = ccenter.Rep() ;
	fill_clone(sp, cellImage, init_ptn) ;
      }
    }

    store<int> sizes ;
    sizes.allocate(cells) ;
    vector<int> cellmap ;
    FORALL(cells,cc) {
      int csz = cellStencil[cc].size() ;
      int bsz = boundary_map[cc].size() ; ;
      vector3d<double>  ccent = ccenter[cc] ;
      vector<vector3d<double> > cdirs(csz+bsz) ;
      for(int i=0;i<csz;++i) {
	cdirs[i] = ccenter[cellStencil[cc][i]]-ccent ;
	cdirs[i] *= 1./norm(cdirs[i]) ;
      }
      for(int i=0;i<bsz;++i) {
	cdirs[csz+i] = fcenter[boundary_map[cc][i]]-ccent ;
	cdirs[csz+i] *= 1./norm(cdirs[csz+i]) ;
      }
	
      vector<int> flags(csz+bsz,0) ;
      int lsz = lower[cc].size() ;
      for(int i=0;i<lsz;++i) {
	vector3d<double>  dir = fcenter[lower[cc][i]] -ccent;
	dir *= 1./norm(dir) ;
	int minid = 0 ;
	int maxid = 0 ;
	double minval = dot(dir,cdirs[0]) ;
	double maxval = minval ;
	for(int j=1;j<csz+bsz;++j) {
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
	flags[minid] = 1 ;
	flags[maxid] = 1 ;
	minid = 0 ;
	maxid = 0 ;
	dir = fnormal[lower[cc][i]] ;
	minval = dot(dir,cdirs[0]) ;
	maxval = minval ;
	for(int j=1;j<csz+bsz;++j) {
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
	flags[minid] = 1 ;
	flags[maxid] = 1 ;
      }
      int usz = upper[cc].size() ;
      for(int i=0;i<usz;++i) {
	vector3d<double>  dir = fcenter[upper[cc][i]] -ccent;
	dir *= 1./norm(dir) ;
	int minid = 0 ;
	int maxid = 0 ;
	double minval = dot(dir,cdirs[0]) ;
	double maxval = minval ;
	for(int j=1;j<csz+bsz;++j) {
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
	flags[minid] = 1 ;
	flags[maxid] = 1 ;
	minid = 0 ;
	maxid = 0 ;
	dir = fnormal[upper[cc][i]] ;
	minval = dot(dir,cdirs[0]) ;
	maxval = minval ;
	for(int j=1;j<csz+bsz;++j) {
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
	flags[minid] = 1 ;
	flags[maxid] = 1 ;
      }

      for(int i=0;i<bsz;++i) {
	vector3d<double>  dir = fcenter[boundary_map[cc][i]] -ccent;
	dir *= 1./norm(dir) ;
	int minid = 0 ;
	int maxid = 0 ;
	double minval = dot(dir,cdirs[0]) ;
	double maxval = minval ;
	for(int j=1;j<csz+bsz;++j) {
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
	flags[minid] = 1 ;
	flags[maxid] = 1 ;
	minid = 0 ;
	maxid = 0 ;
	dir = fnormal[boundary_map[cc][i]] ;
	minval = dot(dir,cdirs[0]) ;
	maxval = minval ;
	for(int j=1;j<csz+bsz;++j) {
	  double v = dot(dir,cdirs[i]) ;
	  if(minval < v) {
	    minid = j ;
	    minval = v ;
	  }
	  if(maxval > v) {
	    maxid = j ;
	    maxval = v ;
	  }
	}
	flags[minid] = 1 ;
	flags[maxid] = 1 ;
      }

      int cnt = 0 ;
      for(int i=0;i<csz;++i)
	if(flags[i] > 0) {
	  cellmap.push_back(cellStencil[cc][i]) ;
	  cnt++ ;
	}
      sizes[cc] = cnt ;
    } ENDFORALL ;

    multiMap cellStencilFiltered ;
    cellStencilFiltered.allocate(sizes) ;
    int cnt = 0 ;
    FORALL(cells,cc) {
      for(int i=0;i<cellStencilFiltered[cc].size();++i) {
	cellStencilFiltered[cc][i] = cellmap[cnt] ;
	cnt++ ;
      }
    } ENDFORALL ;
    // Put in fact database
    facts.create_fact("cellStencil",cellStencilFiltered) ;
  }

  
  void createLowerUpper(fact_db &facts) {
    constraint geom_cells,interior_faces,boundary_faces ;
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
    }

  }

  // this is a general routine that balances the pair vector
  // on each process, it redistributes the pair vector
  // among the processes such that each process holds
  // roughly the same number of elements, it maintains the
  // original global elements ordering in the redistribution
  void
  parallel_balance_pair_vector(vector<pair<int,int> >& vp,
                               MPI_Comm comm) {
    int num_procs = 0 ;
    MPI_Comm_size(comm,&num_procs) ;

    // we still use an all-to-all personalized communication
    // algorithm to balance the element numbers on processes.
    // we pick (p-1) equally spaced element as the splitters
    // and then re-split the global vector sequence to balance
    // the number of elements on processes.

    int vp_size = vp.size() ;
    int global_vp_size = 0 ;
    MPI_Allreduce(&vp_size, &global_vp_size,
                  1, MPI_INT, MPI_SUM, comm) ;

    int space = global_vp_size / num_procs ;
    // compute a global range for the elements on each process
    int global_end = 0 ;
    MPI_Scan(&vp_size, &global_end, 1, MPI_INT, MPI_SUM, comm) ;
    int global_start = global_end - vp_size ;

    vector<int> splitters(num_procs) ;
    // splitters are just global index number
    splitters[0] = space ;
    for(int i=1;i<num_procs-1;++i)
      splitters[i] = splitters[i-1] + space ;
    splitters[num_procs-1] = global_vp_size ;

    // split and communicate the vector of particles
    vector<int> send_counts(num_procs, 0) ;
    int part_start = global_start ;
    for(int idx=0;idx<num_procs;++idx) {
      if(part_start == global_end)
        break ;
      if(splitters[idx] > part_start) {
        int part_end ;
        if(splitters[idx] < global_end)
          part_end = splitters[idx] ;
        else
          part_end = global_end ;
        send_counts[idx] = part_end - part_start ;
        part_start = part_end ;
      }
    }
      
    for(size_t i=0;i<send_counts.size();++i)
      send_counts[i] *= 2 ;

    vector<int> send_displs(num_procs) ;
    send_displs[0] = 0 ;
    for(int i=1;i<num_procs;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    vector<int> recv_counts(num_procs) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;

    vector<int> recv_displs(num_procs) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<num_procs;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    int total_recv_size = recv_displs[num_procs-1] +
      recv_counts[num_procs-1] ;

    // prepare send and recv buffer
    vector<int> send_buf(vp_size*2) ;
    int count = 0 ;
    for(int i=0;i<vp_size;++i) {
      send_buf[count++] = vp[i].first ;
      send_buf[count++] = vp[i].second ;
    }
    // release vp buffer to save some memory because we no longer need it
    vector<pair<int,int> >().swap(vp) ;
    // prepare recv buffer
    vector<int> recv_buf(total_recv_size) ;

    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_INT,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_INT, comm) ;
    // finally extract the data to fill the pair vector
    // release send_buf first to save some memory
    vector<int>().swap(send_buf) ;
    vp.resize(total_recv_size/2) ;
    count = 0 ;
    for(int i=0;i<total_recv_size;i+=2,count++)
      vp[count] = pair<int,int>(recv_buf[i], recv_buf[i+1]) ;
  }
  
  // a parallel sample sort for vector<pair<int, int> >
  // the passed in vector is the local SORTED data. NOTE:
  // the precondition to this routine is that the passed
  // in vector is sorted!!!
  // after sorting, this function puts the new sorted pairs
  // that are local to a processor in the data argument.
  void
  par_sort(vector<pair<int,int> >& data, MPI_Comm comm) {
    // first get the processor id and total number of processors
    int my_id, num_procs ;
    MPI_Comm_size(comm, &num_procs) ;
    MPI_Comm_rank(comm, &my_id) ;
    if(num_procs <= 1)
      return ;                  // single process, no need to proceed
    // get the number of local elements
    int local_size = data.size() ;
    // then select num_procs-1 equally spaced elements as splitters
    int* splitters = new int[num_procs] ;
    int even_space = local_size / (num_procs-1) ;
    int start_idx = even_space / 2 ;
    int space_idx = start_idx ;
    for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
      splitters[i] = data[space_idx].first ;
    // gather the splitters to all processors as samples
    int sample_size = num_procs * (num_procs-1) ;
    int* samples = new int[sample_size] ;
    MPI_Allgather(splitters, num_procs-1, MPI_INT,
                  samples, num_procs-1, MPI_INT, comm) ;
    // now we've obtained all the samples, first we sort them
    sort(samples, samples+sample_size) ;
    // select new splitters in the sorted samples
    even_space = sample_size / (num_procs-1) ;
    start_idx = even_space / 2 ;
    space_idx = start_idx ;
    for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
      splitters[i] = samples[space_idx] ;
    // the last one set as maximum possible integer
    splitters[num_procs-1] = std::numeric_limits<int>::max() ;

    // now we can assign local elements to buckets (processors)
    // according to the new splitters. first we will compute
    // the size of each bucket and communicate them first
    int* scounts = new int[num_procs] ;
    for(int i=0;i<num_procs;++i)
      scounts[i] = 0 ;
    { // using a block just to make the definition of "i" and "j" local
      int i, j ;
      for(j=i=0;i<local_size;++i) {
        if(data[i].first < splitters[j])
          scounts[j]++ ;
        else {
          ++j ;
          while(data[i].first >= splitters[j]) {
            scounts[j] = 0 ;
            ++j ;
          }
          scounts[j]++ ;
        }
      }
    }
    // but since one local element contains two integers (a pair of int),
    // we will need to double the size
    for(int i=0;i<num_procs;++i)
      scounts[i] *= 2 ;
    // now we compute the sending displacement for each bucket
    int* sdispls = new int[num_procs] ;
    sdispls[0] = 0 ;
    for(int i=1;i<num_procs;++i)
      sdispls[i] = sdispls[i-1] + scounts[i-1] ;
    // communicate this information to all processors so that each will
    // know how many elements are expected from every other processor
    int* rcounts = new int[num_procs] ;
    MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm) ;
    // then based on the received info. we will need to compute the
    // receive displacement
    int* rdispls = new int[num_procs] ;
    rdispls[0] = 0 ;
    for(int i=1;i<num_procs;++i)
      rdispls[i] = rdispls[i-1] + rcounts[i-1] ;
    // then we will need to pack the elements in local into
    // a buffer and communicate them
    int* local_pairs = new int[local_size*2] ;
    int count = 0 ;
    for(int i=0;i<local_size;++i) {
      local_pairs[count++] = data[i].first ;
      local_pairs[count++] = data[i].second ;
    }
    // then we allocate buffer for new local elements
    int new_local_size = rdispls[num_procs-1] + rcounts[num_procs-1] ;
    int* sorted_pairs = new int[new_local_size] ;
    // finally we communicate local_pairs to each processor
    MPI_Alltoallv(local_pairs, scounts, sdispls, MPI_INT,
                  sorted_pairs, rcounts, rdispls, MPI_INT, comm) ;
    // release buffers
    delete[] splitters ;
    delete[] samples ;
    delete[] scounts ;
    delete[] sdispls ;
    delete[] rcounts ;
    delete[] rdispls ;
    delete[] local_pairs ;
    // finally we unpack the buffer into a vector of pairs
    data.resize(new_local_size/2) ;
    int data_idx = 0 ;
    for(int i=0;i<new_local_size;i+=2,data_idx++)
      data[data_idx] = pair<int,int>(sorted_pairs[i],sorted_pairs[i+1]) ;
    // release the final buffer
    delete[] sorted_pairs ;
    // finally we sort the new local vector
    sort(data.begin(), data.end()) ;
  }

  namespace {
    // a utility that returns the global sum
    int
    global_sum(int l) {
      int g ;
      MPI_Allreduce(&l, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      return g ;
    }
    // a utility function that takes an entitySet from a processor
    // and returns a vector of entitySet gathered from all processors
    vector<entitySet>
    gather_all_entitySet(const entitySet& eset) {
      int local_size = eset.size() ;
      int global_size = global_sum(local_size) ;
      // compute receive counts from all processors
      int* recv_counts = new int[MPI_processes] ;
      MPI_Allgather(&local_size, 1, MPI_INT,
                    recv_counts, 1, MPI_INT, MPI_COMM_WORLD) ;
      // then compute receive displacement
      int* recv_displs = new int[MPI_processes] ;
      recv_displs[0] = 0 ;
      for(int i=1;i<MPI_processes;++i)
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
      // pack the local eset into an array
      int* local_eset = new int[local_size] ;
      int count = 0 ;
      for(entitySet::const_iterator ei=eset.begin();
          ei!=eset.end();++ei,++count)
        local_eset[count] = *ei ;
      // allocate the entire array for all data from all processors
      int* global_eset = new int[global_size] ;
      // communicate to obtain all esets from every processors
      MPI_Allgatherv(local_eset, local_size, MPI_INT,
                     global_eset, recv_counts, recv_displs,
                     MPI_INT, MPI_COMM_WORLD) ;
      delete[] local_eset ;
      delete[] recv_counts ;
      // unpack the raw buffer into a vector<entitySet>
      vector<entitySet> ret(MPI_processes) ;
      int k = 0 ;
      for(int i=0;i<MPI_processes;++i) {
        int limit ;
        if(i == MPI_processes-1)
          limit = global_size ;
        else
          limit = recv_displs[i+1] ;
        for(;k<limit;++k)
          ret[i] += global_eset[k] ;
      }
      delete[] recv_displs ;
      delete[] global_eset ;

      return ret ;
    }
    void
    par_sort2(vector<pair<pair<int,int>, int> >& data, MPI_Comm comm){
      // first get the processor id and total number of processors
      int my_id, num_procs ;
      MPI_Comm_size(comm, &num_procs) ;
      MPI_Comm_rank(comm, &my_id) ;
      if(num_procs <= 1)
        return ;                  // single process, no need to proceed
      // get the number of local elements
      int local_size = data.size() ;
      // then select num_procs-1 equally spaced elements as splitters
      pair<int, int> *splitters = new pair<int, int>[num_procs] ;
      int even_space = local_size / (num_procs-1) ;
      int start_idx = even_space / 2 ;
      int space_idx = start_idx ;
      for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
        splitters[i] = data[space_idx].first ;
      // gather the splitters to all processors as samples
      int sample_size = num_procs * (num_procs-1) ;
      pair<int, int>* samples = new pair<int, int>[sample_size] ;
      MPI_Allgather(splitters, (num_procs-1)*2, MPI_INT,
                    samples, (num_procs-1)*2, MPI_INT, comm) ;
      // now we've obtained all the samples, first we sort them
      sort(samples, samples+sample_size) ;
      // select new splitters in the sorted samples
      even_space = sample_size / (num_procs-1) ;
      start_idx = even_space / 2 ;
      space_idx = start_idx ;
      for(int i=0;i<num_procs-1;++i,space_idx+=even_space)
        splitters[i] = samples[space_idx] ;
      // the last one set as maximum possible integer
      int maxnumber = std::numeric_limits<int>::max();
      splitters[num_procs-1] =pair<int, int>(maxnumber, maxnumber);
                               
  
      // now we can assign local elements to buckets (processors)
      // according to the new splitters. first we will compute
      // the size of each bucket and communicate them first
      int *scounts = new int[num_procs] ;
      for(int i=0;i<num_procs;++i)
        scounts[i] = 0;
      { // using a block just to make the definition of "i" and "j" local
        int i, j ;
        for(j=i=0;i<local_size;++i) {
          if(data[i].first < splitters[j])
            scounts[j]++ ;
          else {
            ++j ;
            while(data[i].first >= splitters[j]) {
              scounts[j] = 0 ;
              ++j ;
            }
            scounts[j]++ ;
          }
        }
      }
      // but since one local element contains two integers (a pair of int),
      // we will need to double the size
      for(int i=0;i<num_procs;++i)
        scounts[i] *= 3 ;
      // now we compute the sending displacement for each bucket
      int* sdispls = new int[num_procs] ;
      sdispls[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        sdispls[i] = sdispls[i-1] + scounts[i-1] ;
      // communicate this information to all processors so that each will
      // know how many elements are expected from every other processor
      int* rcounts = new int[num_procs] ;
      MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm) ;
      // then based on the received info. we will need to compute the
      // receive displacement
      int* rdispls = new int[num_procs] ;
      rdispls[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        rdispls[i] = rdispls[i-1] + rcounts[i-1] ;
      // then we will need to pack the elements in local into
      // a buffer and communicate them
      int* local_pairs = new int[local_size*3] ;
      int count = 0 ;
      for(int i=0;i<local_size;++i) {
        local_pairs[count++] = data[i].first.first ;
        local_pairs[count++] = data[i].first.second;
        local_pairs[count++] = data[i].second;
      }
      // then we allocate buffer for new local elements
      int new_local_size = rdispls[num_procs-1] + rcounts[num_procs-1] ;
      int* sorted_pairs = new int[new_local_size] ;
      // finally we communicate local_pairs to each processor
      MPI_Alltoallv(local_pairs, scounts, sdispls, MPI_INT,
                    sorted_pairs, rcounts, rdispls, MPI_INT, comm) ;
      // release buffers
      delete[] splitters ;
      delete[] samples ;
      delete[] scounts ;
      delete[] sdispls ;
      delete[] rcounts ;
      delete[] rdispls ;
      delete[] local_pairs ;
      // finally we unpack the buffer into a vector of pairs
      data.resize(new_local_size/3) ;
      int data_idx = 0 ;
      for(int i=0;i<new_local_size;i+=3,data_idx++)
        data[data_idx] = pair<pair<int,int>, int>(pair<int, int>(sorted_pairs[i],sorted_pairs[i+1]), sorted_pairs[i+2]) ;
      // release the final buffer
      delete[] sorted_pairs ;
      // finally we sort the new local vector
      sort(data.begin(), data.end()) ;
    }


    void
    parallel_balance_pair2_vector(vector<pair<pair<int,int>, int> >& vp,
                                  MPI_Comm comm) {
      int num_procs = 0 ;
      MPI_Comm_size(comm,&num_procs) ;
  
      // we still use an all-to-all personalized communication
      // algorithm to balance the element numbers on processes.
      // we pick (p-1) equally spaced element as the splitters
      // and then re-split the global vector sequence to balance
      // the number of elements on processes.
  
      int vp_size = vp.size() ;
      int global_vp_size = 0 ;
      MPI_Allreduce(&vp_size, &global_vp_size,
                    1, MPI_INT, MPI_SUM, comm) ;
  
      int space = global_vp_size / num_procs ;
      // compute a global range for the elements on each process
      int global_end = 0 ;
      MPI_Scan(&vp_size, &global_end, 1, MPI_INT, MPI_SUM, comm) ;
      int global_start = global_end - vp_size ;
  
      vector<int> splitters(num_procs) ;
      // splitters are just global index number
      splitters[0] = space ;
      for(int i=1;i<num_procs-1;++i)
        splitters[i] = splitters[i-1] + space ;
      splitters[num_procs-1] = global_vp_size ;
  
      // split and communicate the vector of particles
      vector<int> send_counts(num_procs, 0) ;
      int part_start = global_start ;
      for(int idx=0;idx<num_procs;++idx) {
        if(part_start == global_end)
          break ;
        if(splitters[idx] > part_start) {
          int part_end ;
          if(splitters[idx] < global_end)
            part_end = splitters[idx] ;
          else
            part_end = global_end ;
          send_counts[idx] = part_end - part_start ;
          part_start = part_end ;
        }
      }
  
      for(size_t i=0;i<send_counts.size();++i)
        send_counts[i] *= 3 ;
  
      vector<int> send_displs(num_procs) ;
      send_displs[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
  
      vector<int> recv_counts(num_procs) ;
      MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                   &recv_counts[0], 1, MPI_INT, comm) ;
  
      vector<int> recv_displs(num_procs) ;
      recv_displs[0] = 0 ;
      for(int i=1;i<num_procs;++i)
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
  
      int total_recv_size = recv_displs[num_procs-1] +
        recv_counts[num_procs-1] ;
  
      // prepare send and recv buffer
      vector<int> send_buf(vp_size*3) ;
      int count = 0 ;
      for(int i=0;i<vp_size;++i) {
        send_buf[count++] = vp[i].first.first ;
        send_buf[count++] = vp[i].first.second ;
        send_buf[count++] = vp[i].second;
      }
      // release vp buffer to save some memory because we no longer need it
      vector<pair<pair<int,int>, int>  >().swap(vp) ;
      // prepare recv buffer
      vector<int> recv_buf(total_recv_size) ;
  
      MPI_Alltoallv(&send_buf[0], &send_counts[0],
                    &send_displs[0], MPI_INT,
                    &recv_buf[0], &recv_counts[0],
                    &recv_displs[0], MPI_INT, comm) ;
      // finally extract the data to fill the pair vector
      // release send_buf first to save some memory
      vector<int>().swap(send_buf) ;
      vp.resize(total_recv_size/3) ;
      count = 0 ;
      for(int i=0;i<total_recv_size;i+=3,count++)
        vp[count] = pair<pair<int,int>, int>(pair<int, int>(recv_buf[i], recv_buf[i+1]), recv_buf[i+2]) ;
    }

 
  }
  
 
  void
  createEdgesPar(fact_db &facts) {
    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    entitySet faces = face2node.domain() ;
    
    // Loop over faces and create list of edges (with duplicates)
    vector<pair<Entity,Entity> > emap ;
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      for(int i=0;i<sz-1;++i) {
        Entity e1 = face2node[*ei][i] ;
        Entity e2 = face2node[*ei][i+1] ;
        emap.push_back(pair<Entity,Entity>(min(e1,e2),max(e1,e2))) ;
      }
      Entity e1 = face2node[*ei][0] ;
      Entity e2 = face2node[*ei][sz-1] ;
      emap.push_back(pair<Entity,Entity>(min(e1,e2),max(e1,e2))) ;
    }

    // before we do the parallel sorting, we perform a check
    // to see if every process at least has one data element in
    // the "emap", if not, then the parallel sample sort would fail
    // and we pre-balance the "emap" on every process before the
    // sorting
    if(GLOBAL_OR(emap.empty())) {
      parallel_balance_pair_vector(emap, MPI_COMM_WORLD) ;
    }
    // Sort edges and remove duplicates
    sort(emap.begin(),emap.end()) ;
    vector<pair<Entity,Entity> >::iterator uend ;
    uend = unique(emap.begin(), emap.end()) ;
    emap.erase(uend, emap.end()) ;
    // then sort emap in parallel
    // but we check again to see if every process has at least one
    // element, if not, that means that the total element number is
    // less than the total number of processes, we split the communicator
    // so that only those do have elements would participate in the
    // parallel sample sorting
    if(GLOBAL_OR(emap.empty())) {
      MPI_Comm sub_comm ;
      int color = emap.empty() ;
      MPI_Comm_split(MPI_COMM_WORLD, color, MPI_rank, &sub_comm) ;
      if(!emap.empty())
        par_sort(emap, sub_comm) ;
      MPI_Comm_free(&sub_comm) ;
    } else {
      par_sort(emap, MPI_COMM_WORLD) ;
    }
    // remove duplicates again in the new sorted vector
    uend = unique(emap.begin(), emap.end()) ;
    emap.erase(uend, emap.end()) ;
#ifdef BOUNDARY_DUPLICATE_DETECT
    if(MPI_processes > 1) {
      // then we will need to remove duplicates along the boundaries
      // we send the first element in the vector to the left neighbor
      // processor (my_id - 1) and each processor compares its last
      // element with the received element. if they are the same,
      // then the processor will remove its last element
      
      // HOWEVER if the parallel sort was done using the sample sort
      // algorithm, then this step is not necessary. Because in the
      // sample sort, elements are partitioned to processors according
      // to sample splitters, it is therefore guaranteed that no
      // duplicates will be crossing the processor boundaries.
      int sendbuf[2] ;
      int recvbuf[2] ;
      if(!emap.empty()) {
        sendbuf[0] = emap[0].first ;
        sendbuf[1] = emap[0].second ;
      } else {
        // if there is no local data, we set the send buffer
        // to be the maximum integer so that we don't have problems
        // in the later comparing stage
        sendbuf[0] = std::numeric_limits<int>::max() ;
        sendbuf[1] = std::numeric_limits<int>::max() ;
      }
      MPI_Status status ;
      if(MPI_rank == 0) {
        // rank 0 only receives from 1, no sending needed
        MPI_Recv(recvbuf, 2, MPI_INT,
                 1/*source*/, 0/*msg tag*/,
                 MPI_COMM_WORLD, &status) ;
      } else if(MPI_rank == MPI_processes-1) {
        // the last processes only sends to the second last processes,
        // no receiving is needed
        MPI_Send(sendbuf, 2, MPI_INT,
                 MPI_rank-1/*dest*/, 0/*msg tag*/, MPI_COMM_WORLD) ;
      } else {
        // others will send to MPI_rank-1 and receive from MPI_rank+1
        MPI_Sendrecv(sendbuf, 2, MPI_INT, MPI_rank-1/*dest*/,0/*msg tag*/,
                     recvbuf, 2, MPI_INT, MPI_rank+1/*source*/,0/*tag*/,
                     MPI_COMM_WORLD, &status) ;
      }
      // then compare the results with last element in local emap
      if( (MPI_rank != MPI_processes-1) && (!emap.empty())){
        const pair<Entity,Entity>& last = emap.back() ;
        if( (recvbuf[0] == last.first) &&
            (recvbuf[1] == last.second)) {
          emap.pop_back() ;
        }
      }
    } // end if(MPI_Processes > 1)
#endif
    
    // Allocate entities for new edges
    int num_edges = emap.size() ;
    int ek = facts.getKeyDomain("Edges") ;
    if(!useDomainKeySpaces)
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
    edge.Rep()->setDomainKeySpace(ek) ;
    edge.allocate(edges) ;
    vector<pair<Entity,Entity> >::iterator pi = emap.begin() ;
    for(entitySet::const_iterator ei=edges.begin();
        ei!=edges.end();++ei,++pi) {
      edge[*ei][0] = pi->first ;
      edge[*ei][1] = pi->second ;
    }

    // Add edge2node data structure to fact databse
    // facts.create_fact("edge2node",edge) ;

    // Now create face2edge data-structure
    // We need to create a lower node to edge mapping to facilitate the
    // searches.  First get map from edge to lower node
    Map el ; // Lower edge map
    el.allocate(edges) ;
    for(entitySet::const_iterator ei=edges.begin();
        ei!=edges.end();++ei,++pi) {
      el[*ei] = edge[*ei][0] ;
    }
    
    // Now invert this map to get nodes-> edges that have this as a first entry
    multiMap n2e ;
    // Get nodes
    // Get mapping from nodes to edges from lower numbered node
    
    // note inorder to use the distributed_inverseMap, we need
    // to provide a vector of entitySet partitions. for this 
    // case, it is NOT the node (pos.domain()) distribution,
    // instead it is the el Map image distribution
    entitySet el_image = el.image(el.domain()) ;
    vector<entitySet> el_image_partitions =
      gather_all_entitySet(el_image) ;
    distributed_inverseMap(n2e, el, el_image, edges, el_image_partitions) ;

    // Now create face2edge map with same size as face2node
    multiMap face2edge ;
    store<int> count ;
    count.allocate(faces) ;
    for(entitySet::const_iterator ei = faces.begin();
        ei!=faces.end();++ei) 
      count[*ei] = face2node[*ei].size() ;
    face2edge.allocate(count) ;

    // before computing the face2edge map, we will need to gather
    // necessary info among all processors since the edge map is
    // distributed across all the processors. we need to retrieve
    // those that are needed from other processors.

    // we will first need to figure out the set of edges we need
    // but are not on the local processor

    // but we need to access the n2e map in the counting and it
    // is possible that the local n2e map does not have enough
    // data we are looking for, therefore we need to expand it
    // first to include possible clone regions
    entitySet nodes_accessed ;
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      for(int i=0;i<sz-1;++i) {
        Entity t1 = face2node[*ei][i] ;
        Entity t2 = face2node[*ei][i+1] ;
        Entity e1 = min(t1,t2) ;
        nodes_accessed += e1 ;
      }
      // Work on closing edge
      Entity t1 = face2node[*ei][0] ;
      Entity t2 = face2node[*ei][sz-1] ;
      Entity e1 = min(t1,t2) ;
      nodes_accessed += e1 ;
    }
    // we then expand the n2e map
    entitySet nodes_out_domain = nodes_accessed - n2e.domain() ;    
    n2e.setRep(MapRepP(n2e.Rep())->expand(nodes_out_domain,
                                          el_image_partitions)) ;
    // okay, then we are going to expand the edge map
    // first count all the edges we need
    entitySet edges_accessed ;
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      for(int i=0;i<sz-1;++i) {
        Entity t1 = face2node[*ei][i] ;
        Entity t2 = face2node[*ei][i+1] ;
        Entity e1 = min(t1,t2) ;
        for(int j=0;j<n2e[e1].size();++j) {
          int e = n2e[e1][j] ;
          edges_accessed += e ;
        }
      }
      // Work on closing edge
      Entity t1 = face2node[*ei][0] ;
      Entity t2 = face2node[*ei][sz-1] ;
      Entity e1 = min(t1,t2) ;
      for(int j=0;j<n2e[e1].size();++j) {
        int e = n2e[e1][j] ;
        edges_accessed += e ;
      }
    }
    vector<entitySet> edge_partitions = gather_all_entitySet(edge.domain()) ;
    entitySet edges_out_domain = edges_accessed - edge.domain() ;
    // but since there is no expand method implemented for
    // MapVec at this time, we will just do a hack to convert
    // the MapVec to a multiMap to reuse the expand code.
    multiMap edge2 ;
    store<int> edge2_count ;
    entitySet edge_dom = edge.domain() ;
    edge2_count.allocate(edge_dom) ;
    for(entitySet::const_iterator ei=edge_dom.begin();
        ei!=edge_dom.end();++ei)
      edge2_count[*ei] = 2 ;
    edge2.allocate(edge2_count) ;
    for(entitySet::const_iterator ei=edge_dom.begin();
        ei!=edge_dom.end();++ei) {
      edge2[*ei][0] = edge[*ei][0] ;
      edge2[*ei][1] = edge[*ei][1] ;
    }
    edge2.setRep(MapRepP(edge2.Rep())->expand(edges_out_domain,
                                              edge_partitions)) ;
    // we are now ready for the face2edge map

    // Now loop over faces, for each face search for matching edge and
    // store in the new face2edge structure
    for(entitySet::const_iterator ei=faces.begin();
        ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      // Loop over edges of the face
      for(int i=0;i<sz-1;++i) {
        Entity t1 = face2node[*ei][i] ;
        Entity t2 = face2node[*ei][i+1] ;
        Entity e1 = min(t1,t2) ;
        Entity e2 = max(t1,t2) ;
        face2edge[*ei][i] = -1 ;
        // search for matching edge
        for(int j=0;j<n2e[e1].size();++j) {
          int e = n2e[e1][j] ;
          if(edge2[e][0] == e1 && edge2[e][1] == e2) {
            face2edge[*ei][i] = e ;
            break ;
          }
        }
        if(face2edge[*ei][i] == -1)
          cerr << "ERROR: not able to find edge for face " << *ei << endl ;
      }
      // Work on closing edge
      Entity t1 = face2node[*ei][0] ;
      Entity t2 = face2node[*ei][sz-1] ;
      Entity e1 = min(t1,t2) ;
      Entity e2 = max(t1,t2) ;
      face2edge[*ei][sz-1] = -1 ;
      for(int j=0;j<n2e[e1].size();++j) {
        int e = n2e[e1][j] ;
        if(edge2[e][0] == e1 && edge2[e][1] == e2) {
          face2edge[*ei][sz-1] = e ;
          break ;
        }
      }
      if(face2edge[*ei][sz-1] == -1)
        cerr << "ERROR: not able to find edge for face " << *ei << endl ;
      
    }
    // Add face2edge to the fact database
    face2edge.Rep()->setDomainKeySpace(fk) ;
    MapRepP(face2edge.Rep())->setRangeKeySpace(ek) ;
    facts.create_fact("face2edge",face2edge) ;


    //sort edge2node according to fileNumbering
    if(MPI_processes > 1){    
      //create Map node_l2f
      entitySet nodes;
      Map node_l2f;  
      FORALL(edges, e){
        nodes += edge[e][0];
        nodes += edge[e][1];
      }ENDFORALL;
    
      
      storeRepP pos = facts.get_variable("pos");
      int nkeyspace = pos->getDomainKeySpace() ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn(nkeyspace) ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      dMap g2f ;
      g2f = df->g2fv[nkeyspace].Rep() ; 
      //don't use nodes & init_ptn to define local nodes,
      //because nodes may not cover all nodes in init_ptn
      entitySet localNodes = pos->domain()&init_ptn[MPI_rank] ;
      node_l2f.allocate(localNodes);
      FORALL(localNodes, d){
        node_l2f[d] = g2f[d];
      }ENDFORALL;
       
      entitySet out_of_dom = nodes - localNodes;
      // vector<entitySet> tmp_ptn = gather_all_entitySet(localNodes);
      node_l2f.setRep(MapRepP(node_l2f.Rep())->expand(out_of_dom, init_ptn)) ;
     
      
      //end of create Map
       
      FORALL(edge.domain(), e){
        if(node_l2f[edge[e][0] ]> node_l2f[edge[e][1]]){
          std:: swap(edge[e][0], edge[e][1]);
        }
               
      }ENDFORALL;
    
      
  
  
    
           
  

      //then update fact_db so that the file number of edges is consistent with the file number of nodes

      //give each edge a file number
      vector<pair<pair<Entity, Entity> , Entity> > edge2global(num_edges);
      int eindex = 0;
      FORALL(edges, ei){
        edge2global[eindex++] = pair<pair<Entity, Entity>, Entity>(pair<Entity, Entity>(node_l2f[edge[ei][0]], node_l2f[edge[ei][1]]), ei); 
      }ENDFORALL;


      if(GLOBAL_OR(edge2global.empty())) {
        parallel_balance_pair2_vector(edge2global, MPI_COMM_WORLD) ;
      }
      // Sort edges and remove duplicates
      sort(edge2global.begin(),edge2global.end()) ;
      vector<pair<pair<Entity,Entity>, Entity> >::iterator uend2 ;
      uend2 = unique(edge2global.begin(), edge2global.end()) ;
      edge2global.erase(uend2, edge2global.end()) ;
      // then sort emap in parallel
      // but we check again to see if every process has at least one
      // element, if not, that means that the total element number is
      // less than the total number of processes, we split the communicator
      // so that only those do have elements would participate in the
      // parallel sample sorting
      if(GLOBAL_OR(edge2global.empty())) {
        MPI_Comm sub_comm ;
        int color = edge2global.empty() ;
        MPI_Comm_split(MPI_COMM_WORLD, color, MPI_rank, &sub_comm) ;
        if(!edge2global.empty())
          par_sort2(edge2global, sub_comm) ;
        MPI_Comm_free(&sub_comm) ;
      } else {
        par_sort2(edge2global, MPI_COMM_WORLD) ;
      }
      // remove duplicates again in the new sorted vector
      uend2 = unique(edge2global.begin(), edge2global.end()) ;
      edge2global.erase(uend2, edge2global.end()) ;
#ifdef BOUNDARY_DUPLICATE_DETECT
      if(MPI_processes > 1) {
        // then we will need to remove duplicates along the boundaries
        // we send the first element in the vector to the left neighbor
        // processor (my_id - 1) and each processor compares its last
        // element with the received element. if they are the same,
        // then the processor will remove its last element

        // HOWEVER if the parallel sort was done using the sample sort
        // algorithm, then this step is not necessary. Because in the
        // sample sort, elements are partitioned to processors according
        // to sample splitters, it is therefore guaranteed that no
        // duplicates will be crossing the processor boundaries.
        int sendbuf[3] ;
        int recvbuf[3] ;
        if(!edge2global.empty()) {
          sendbuf[0] = edge2global[0].first.first ;
          sendbuf[1] = edge2global[0].first.second ;
          sendbuf[2] = edge2global[0].second;
        } else {
          // if there is no local data, we set the send buffer
          // to be the maximum integer so that we don't have problems
          // in the later comparing stage
          sendbuf[0] = std::numeric_limits<int>::max() ;
          sendbuf[1] = std::numeric_limits<int>::max() ;
          sendbuf[2] = std::numeric_limits<int>::max() ;
        }
        MPI_Status status ;
        if(MPI_rank == 0) {
          // rank 0 only receives from 1, no sending needed
          MPI_Recv(recvbuf, 3, MPI_INT,
                   1/*source*/, 0/*msg tag*/,
                   MPI_COMM_WORLD, &status) ;
        } else if(MPI_rank == MPI_processes-1) {
          // the last processes only sends to the second last processes,
          // no receiving is needed
          MPI_Send(sendbuf, 3, MPI_INT,
                   MPI_rank-1/*dest*/, 0/*msg tag*/, MPI_COMM_WORLD) ;
        } else {
          // others will send to MPI_rank-1 and receive from MPI_rank+1
          MPI_Sendrecv(sendbuf, 3, MPI_INT, MPI_rank-1/*dest*/,0/*msg tag*/,
                       recvbuf, 3, MPI_INT, MPI_rank+1/*source*/,0/*tag*/,
                       MPI_COMM_WORLD, &status) ;
        }
        // then compare the results with last element in local emap
        if( (MPI_rank != MPI_processes-1) && (!edge2global.empty())){
          const pair<pair<Entity,Entity>, Entity>& last = edge2global.back() ;
          if( (recvbuf[0] == last.first.first) &&
              (recvbuf[1] == last.first.second)&&
              (recvbuf[2] == last.secon)) {
            edge2global.pop_back() ;
          }
        }
      } // end if(MPI_Processes > 1)
#endif

    
      int local_num_edge = edge2global.size();
      vector<int> edge_sizes(MPI_processes);
      MPI_Allgather(&local_num_edge,1,MPI_INT,&edge_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    
      int file_num_offset = 0;
      for(int i = 0; i < MPI_rank; i++){
        file_num_offset += edge_sizes[i];
      }
    
      vector<pair<Entity, Entity> > file2global(edge2global.size());
      int index = file_num_offset;
    
      entitySet input_image = edges;
      entitySet input_preimage;
      for(int i = 0; i < local_num_edge; i++){
        input_preimage += index;
        file2global[i] = pair<Entity, Entity>(index, edge2global[i].second);
        index++;
      }
      
      multiMap global2file;
    
      //the input_image is not really the image, it should be the global2file's domain
      //if it doesn't include all corresponding entities in init_ptn, error will occur
      //also input_preimage is never used
      std::vector<entitySet> init_ptne = facts.get_init_ptn(ek) ;

      Loci::distributed_inverseMap(global2file,
                                   file2global,
                                   input_image,
                                   input_preimage,
                                   init_ptne);
    
    
      if(global2file.domain() != edges){
        cerr<<"the inversed map doesn't match edges" << endl;
        cerr <<"domain: " << global2file.domain() << " edge:   " << edges << endl;
        Loci::Abort();
      }
    
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;

      FORALL(global2file.domain(), ei){
        dist->g2fv[0][ei] = global2file[ei][0] ;
      }ENDFORALL;
    
    }
    //before put edge2node to fact_db, make sure each edge point from lower
    //file number node to higher file number node
    MapVec<2> edge3 ;
    edge3.allocate(edges) ;
    FORALL(edges, ei) {
      edge3[ei][0] = edge[ei][0];
      edge3[ei][1] = edge[ei][1];
    }ENDFORALL;
    
    // Add edge3node data structure to fact databse
    edge3.Rep()->setDomainKeySpace(ek) ;
    facts.create_fact("edge2node",edge3) ;  
    
  } // end of createEdgesPar


  void setupOverset(fact_db &facts) {
    using namespace Loci ;
    using std::map ;
    storeRepP sp = facts.get_variable("componentGeometry") ;
    if(sp == 0)
      return ;
    sp = facts.get_variable("ci") ;
    if(sp == 0)
      return ;
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
    map<string,entitySet> volMap ;

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
    map<string,entitySet>::const_iterator mi ;
    for(mi=volMap.begin();mi!=volMap.end();++mi) {
      int ckeyspace = cl.getRangeKeySpace() ;
      std::vector<entitySet> cptn = facts.get_init_ptn(ckeyspace) ;
      entitySet volgather = distribute_entitySet(mi->second,cptn) ;
      dist_expand_entitySet(volgather,domc,cptn) ;
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
