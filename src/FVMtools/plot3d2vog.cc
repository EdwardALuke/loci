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
#include "vogtools.h"
#include <Tools/stream.h>
//#include <Tools/tools.h>
#include <algorithm>
using std::sort ;
using std::swap ;
#include <map>
using std::map ;
#include <vector>
using std::vector ;
using std::string ;
#include <set>
using std::set ;

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

typedef double real ;
typedef Loci::vector3d<real> vect3d ;
#include <Tools/digraph.h>
using Loci::digraph ;
using Loci::component_sort ;

static char bcnames[6][4] = {"IJ1","IJN","IK1","IKN","JK1","JKN"} ;

struct block_topo {
  int ni,nj,nk ;
  int nodes_base ;
  int faces_base ;
  int cells_base ;
  block_topo() { ni=0,nj=0,nk=0,nodes_base=-1 ;}
  block_topo(int ni,int nj,int nk) ;
  ~block_topo() {};
  int num_nodes() { return ni*nj*nk ; }
  int num_faces() { return ni*(nj-1)*(nk-1)+(ni-1)*nj*(nk-1)+
		      (ni-1)*(nj-1)*(nk) ; }
  int num_cells() { return (ni-1)*(nj-1)*(nk-1) ; }

  int naddr(int i,int j,int k) { return i+j*ni+k*ni*nj + nodes_base; }
  int caddr(int i,int j,int k)
  { return (i)+(j)*(ni-1)+(k)*(ni-1)*(nj-1)+cells_base ; }

  int num_interior_faces() {
    return ( (ni-2)*(nj-1)*(nk-1) +
	     (ni-1)*(nj-2)*(nk-1) +
	     (ni-1)*(nj-1)*(nk-2) ) ; }

  void fill_interior_face(int f, Loci::Array<int,4> &nds, int &cl, int &cr) {
    if(f< (ni-2)*(nj-1)*(nk-1)) {
      int i=1+f/((nj-1)*(nk-1)) ;
      int r = f - (i-1)*((nj-1)*(nk-1)) ;
      int j = r/(nk-1) ;
      int k = r-j*(nk-1) ;
      nds[0] = naddr(i,j+0,k+0) ;
      nds[1] = naddr(i,j+1,k+0) ;
      nds[2] = naddr(i,j+1,k+1) ;
      nds[3] = naddr(i,j+0,k+1) ;
      cl = caddr(i-1,j,k) ;
      cr = caddr(i,j,k) ;
    } else {
      f-= (ni-2)*(nj-1)*(nk-1) ;
      if(f < (ni-1)*(nj-2)*(nk-1) ) {
	int j = 1+f/((ni-1)*(nk-1)) ;
	int r = f - (j-1)*((ni-1)*(nk-1)) ;
	int i = r/(nk-1) ;
	int k = r-i*(nk-1) ;
	nds[3] = naddr(i+0,j,k+0) ;
	nds[2] = naddr(i+1,j,k+0) ;
	nds[1] = naddr(i+1,j,k+1) ;
	nds[0] = naddr(i+0,j,k+1) ;
	cl = caddr(i,j-1,k);
	cr = caddr(i,j,k) ;
      } else {
	f -= (ni-1)*(nj-2)*(nk-1) ;
	int k = 1+f/((ni-1)*(nj-1)) ;
	int r = f - (k-1)*((ni-1)*(nj-1)) ;
	int i = r/(nj-1) ;
	int j = r - i*(nj-1) ;
	nds[0] = naddr(i+0,j+0,k) ;
	nds[1] = naddr(i+1,j+0,k) ;
	nds[2] = naddr(i+1,j+1,k) ;
	nds[3] = naddr(i+0,j+1,k) ;
	cl = caddr(i,j,k-1);
	cr = caddr(i,j,k) ;
      }
    }
  }


} ;

block_topo::block_topo(int NI, int NJ, int NK) {
  block_topo::ni = NI ;
  block_topo::nj = NJ ;
  block_topo::nk = NK ;

  cout << "block {" << ni << "," << nj << "," << nk << "}" << endl ;
  nodes_base = -1 ;
}

// Compute a normalized irrational direction vector (one which
// it is unlikely that grid lines will align with)
const double Axx = sqrt(3.0) ;
const double Ayy = 1./sqrt(3.1415927) ;
const double Azz = exp(1.0);

const double Ax = Axx/sqrt(Axx*Axx+Ayy*Ayy+Azz*Azz) ;
const double Ay = Ayy/sqrt(Axx*Axx+Ayy*Ayy+Azz*Azz) ;
const double Az = Azz/sqrt(Axx*Axx+Ayy*Ayy+Azz*Azz) ;

struct block_face_info {
  int is,ie,js,je, tag ;
} ;

struct block_boundary_info {
  vector<block_face_info> block_face[6] ;
} ;


void usage() {
  cout << "Parallel plot3d file to VOG file grid converter." << endl
       << endl
       << "Input file is ascii unformatted multiblock plot3d file with postfix '.grd'" << endl
       << "or binary unformatted multiblock plot3d file with postfix '.grd.b8'" << endl
       << endl << endl ;
  cout << "Usage:" << endl
       << endl
       << "plot3d2vog <options> <basename>" << endl
       << endl
       << "where <basename> is the root name of the file (sans the '.grd' postfix)." << endl
       << " and options must include one of the following grid units specifications:" << endl
       << "   -m  - grid is in meters" << endl
       << "   -cm - grid is in centimeters" << endl
       << "   -mm - grid is in millimeters" << endl
       << "   -in - grid is in inches" << endl
       << "   -ft - grid is in feet" << endl
       << "   -Lref \"1.5 ft\" - grid has specified reference length" << endl
       << endl<<endl ;
  cout << "Other useful options:" << endl
       << endl
       << "   -tol <float> :gluing tolerance (as percentage of min edge length)" << endl
       << "   -maxAspect <float> :maximum face aspect ratio on boundaries for gluing" << endl
       << "   -lefthanded :set this if grid not in right handed coordinates" << endl
       << "   -creaseAngle <float> :used to combine neighboring bc faces" << endl
       << "                         set to -1 to disable feature" << endl
       << "   -bc <filename> :Use this file to specify boundary surfaces" << endl
       << "   -o :disable mesh optimization that reorders nodes in the vog file" << endl ;
  cout << endl ;

  cout << "If you want to control how tags are assigned to boundary" <<endl
       << " faces, then you can do so by providing a supplemental" << endl
       << " boundary specification file with the flag -bc" << endl << endl
       << "Example:" << endl
       << " plot3d2vog -bc grid.bc grid" << endl << endl ;

  cout << "The boundary condition specification file is contains the" <<endl
       << " following information:" << endl << endl
       << " <Number of boundary groups> " << endl
       << " For each group: " << endl
       << " <boundary tag> <number of segments> <boundary name> " << endl
       << " For each segment:" << endl
       << " <block number> <faceid> <index1 start> <index1 end> <index2 start> <index2 end>" << endl
       << "  Where:" << endl
       << " <faceid> is one of six strings: [IJ1,IJN,JK1,JKN,IK1,IKN]" << endl
       << " indices are given in the order indicated by the faceid string."
       << endl <<endl;

}

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#endif

int fscanread(double *entries, int nentries, FILE *fp) {
  static double rval = 0.0;
  static int nvals = 0 ;
  int nr = 0 ;
  for(int i=0;i<nentries;++i) {
    int cnt = 1 ;
    if(nvals  > 0) {
      entries[i] = rval ;
      nvals-- ;
    } else {
      cnt = fscanf(fp,"%*[, \n\r\t]%lf",&entries[i]) ;
      char buf[512] ;
      if(cnt == 1 && fscanf(fp,"%[*]",buf)) {
	nvals = int(entries[i]) ;
	cnt = fscanf(fp,"%lf",&rval) ;
	if(cnt != 1)
	  return nr ;
	entries[i] = rval ;
	nvals-- ;
      }
      if(cnt != 1)
	return nr ;
    }

    if(cnt == 1)
      nr++ ;
    else {

      char buf[512] ;
      fscanf(fp,"%s",buf) ;
      cerr << "failure reading grid near " << buf << endl ;
      return nr ;


    }
  }
  return nr ;
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

bool readP3DGrid(string filename,vector<block_topo> &blockInfo, vector<vect3d> &positions, vector<int>&pos_sizes, int read_type, MPI_Comm Comm) {

  bool binary_mode = false ;
  if(1 == read_type)
    binary_mode = true ;
  int rank = 0 ;
  MPI_Comm_rank(Comm,&rank) ;

  FILE *fp = 0 ;
  int num_blocks = -1 ;
  int num_nodes = -1 ;
  bool fail = 0 ;
  if( 0 == rank) {
    fp = fopen(filename.c_str(),"r") ;
    if(NULL == fp) {
      cerr << "unable to open file '" << filename << "'" << endl ;
      fail = 1 ;
    } else {
      if(binary_mode) {
	int nr = fread(&num_blocks,sizeof(int),1, fp) ;
	if(nr != 1) {
	  cerr << "failed to read number of blocks" << endl ;
	  fail = 1 ;
	}
      } else {
	int nr = fscanf(fp,"%d",&num_blocks) ;
	if(nr != 1) {
	  cerr << "failed to read number of blocks" << endl ;
	  fail = 1 ;
	}
      }
      if(num_blocks < 1 || num_blocks > 10000)  {
	cerr << "failure reading sensible numblocks, check file '" << filename
	     << "'" << endl ;
	fail = 1 ;
      } else {
	blockInfo = vector<block_topo>(num_blocks) ;
	for(int b=0;b<num_blocks;++b) {
	  int ni=-1,nj=-1,nk=-1 ;
	  if(binary_mode) {
	    int nr = fread(&ni,sizeof(int),1,fp) ;
	    nr += fread(&nj,sizeof(int),1,fp) ;
	    nr += fread(&nk,sizeof(int),1,fp) ;
	    if(nr != 3) {
	      cerr << "failed in reading block sizes" << endl ;
	      fail = 1 ;
	    }
	  } else {
	    int nr = fscanf(fp,"%*[, \n\r]%d%*[, \n\r]%d%*[, \n\r]%d",&ni,&nj,&nk) ;
	    if(nr != 3) {
	      cerr << "failed in reading block sizes" << endl ;
	      fail = 1 ;
	    }
	  }
	  blockInfo[b] = block_topo(ni,nj,nk) ;
	  if(ni<0 || nj < 0 || nk < 0) {
	    cerr << "error reading in block sizes from file '" << filename
		 << "'" << endl ;
	    fail = 1 ;
	    break ;
	  }
	}
	int totsz = 0 ;
	for(int b=0;b<num_blocks;++b) {
	  blockInfo[b].nodes_base = totsz ;
	  totsz += blockInfo[b].num_nodes() ;
	}
	num_nodes = totsz ;
	int base = 0 ;
	for(int b=0;b<num_blocks;++b) {
	  blockInfo[b].cells_base = base ;
	  base += blockInfo[b].num_cells() ;
	}
	base = 0 ;
	for(int b=0;b<num_blocks;++b) {
	  blockInfo[b].faces_base = base ;
	  base += blockInfo[b].num_faces() ;
	}
      }
    }
  }


  MPI_Bcast(&fail, 1, MPI_INT,0,Comm) ;
  if(fail != 0)
    return false ;
  MPI_Bcast(&num_blocks,1, MPI_INT,0,Comm) ;
  if(rank != 0)
    blockInfo = vector<block_topo>(num_blocks) ;
  MPI_Bcast(&blockInfo[0],sizeof(block_topo)*num_blocks, MPI_BYTE,0,Comm) ;
  MPI_Bcast(&num_nodes,1,MPI_INT,0,Comm) ;

  // Now allocate the nodes to processors
  // if more than 1 processors don't allocate any to processor zero to give space
  // for gluing (which will be done on processor 0).
  int nP = -1 ;
  MPI_Comm_size(Comm,&nP) ;
  vector<int> node_allocation(nP,0) ;
  if(1 == nP) {
    node_allocation[0] = num_nodes ;
  } else if (2 == nP) {
    node_allocation[1] = num_nodes ;
  } else {
    const int Pdist = nP -1 ;
    const int npp = num_nodes/Pdist ;
    const int rem = num_nodes - npp*Pdist ;
    int check = 0 ;
    for(int p=1;p<nP;++p) {
      node_allocation[p] = npp ;
      if(p <= rem)
	node_allocation[p] += 1 ;
      check+=node_allocation[p] ;
    }
    if(check != num_nodes) {
      cerr << "inconsistency in node allocation!" << endl;
      cerr << "num_nodes = " << num_nodes << endl ;
      cerr << "allocation = [" ;
      for(int i=0;i<nP;++i)
        cerr << " " << node_allocation[i]  ;
      cerr << " ]" << endl ;
    }
  }
  positions = vector<vect3d>(node_allocation[rank]) ;
  pos_sizes = node_allocation ;
  // Read in positions
  if(1 == nP) {
    vector<double> scratch(positions.size()) ;
    for(int b=0;b<num_blocks;++b) {
      int sz = blockInfo[b].num_nodes() ;
      int o = blockInfo[b].nodes_base ;
      //read in x coordinates from block b
      if(binary_mode) {
	int nr = fread(&scratch[0],sizeof(double),sz,fp) ;
	if(nr != sz) {
	  cerr << "read failed reading mesh positions" << endl ;
	  fail = 1 ;
	}
      } else {
	int nr = fscanread(&scratch[0],sz,fp) ;
	if(nr != sz) {
	  cerr << "read failed reading mesh positions" << endl ;
	  fail = 1 ;
	}
      }
      for(int i=0;i<sz;++i)
	positions[o+i].x = scratch[i] ;

      //read in y coordinates from block b
      if(binary_mode) {
	int nr = fread(&scratch[0],sizeof(double),sz,fp) ;
	if(nr != sz) {
	  cerr << "read failed reading mesh positions" << endl ;
	  fail = 1 ;
	}
      } else {
	int nr = fscanread(&scratch[0],sz,fp) ;
	if(nr != sz) {
	  cerr << "read failed reading mesh positions" << endl ;
	  fail = 1 ;
	}
      }
      for(int i=0;i<sz;++i)
	positions[o+i].y = scratch[i] ;

      // read in z coordinates from block b
      if(binary_mode) {
	int nr = fread(&scratch[0],sizeof(double),sz,fp) ;
	if(nr != sz) {
	  cerr << "read failed reading mesh positions" << endl ;
	  fail = 1 ;
	}
      } else {
	int nr = fscanread(&scratch[0],sz,fp) ;
	if(nr != sz) {
	  cerr << "read failed reading mesh positions" << endl ;
	  fail = 1 ;
	}
      }
      for(int i=0;i<sz;++i)
	positions[o+i].z = scratch[i] ;
    }
  } else {
    // multiprocessor reading, not zero nodes are allocated to processor 0
    using Loci::entitySet ;
    using Loci::EMPTY ;
    vector<int> limits(nP,0) ;
    limits[0] = 0 ;

    for(int i=1;i<nP;++i)
      limits[i] = limits[i-1]+node_allocation[i] ;
    vector<entitySet> distribution(nP) ;
    distribution[0] = EMPTY ;
    for(int i=1;i<nP;++i)
      if(0 == node_allocation[i])
	distribution[i] = EMPTY ;
      else
	distribution[i] = entitySet(interval(limits[i-1],limits[i]-1)) ;

    vector<double> scratch(node_allocation[1]) ;

    if(0 == rank) {
      // processor zero reads
      for(int b=0;b<num_blocks;++b) {
	int start = blockInfo[b].nodes_base ;
	int end = start+blockInfo[b].num_nodes() ;
	entitySet blockSet = interval(start,end-1) ;

	// read in x coordinates
	for(int i=1;i<nP;++i) {
	  entitySet sendSet = blockSet & distribution[i] ;
	  if(sendSet != EMPTY) {
	    int sendsz = sendSet.size() ;
	    if(binary_mode) {
	      int nr = fread(&scratch[0],sizeof(double),sendsz,fp) ;
	      if(nr != sendsz) {
		cerr << "read failed reading in mesh positions" << endl ;
		fail = 1 ;
	      }
	    } else {
	      int nr = fscanread(&scratch[0],sendsz,fp) ;
	      if(nr != sendsz) {
		cerr << "read failed reading in mesh positions" << endl ;
		fail = 1 ;
	      }
	    }
	    MPI_Send(&scratch[0],sendsz,MPI_DOUBLE,i,0,Comm) ;
	  }
	}

	// read in y coordinates
	for(int i=1;i<nP;++i) {
	  entitySet sendSet = blockSet & distribution[i] ;
	  if(sendSet != EMPTY) {
	    int sendsz = sendSet.size() ;
	    if(binary_mode) {
	      int nr = fread(&scratch[0],sizeof(double),sendsz,fp) ;
	      if(nr != sendsz) {
		cerr << "read failed reading in mesh positions" << endl ;
		fail = 1 ;
	      }
	    } else {
	      int nr = fscanread(&scratch[0],sendsz,fp) ;
	      if(nr != sendsz) {
		cerr << "read failed reading in mesh positions" << endl ;
		fail = 1 ;
	      }
	    }
	    MPI_Send(&scratch[0],sendsz,MPI_DOUBLE,i,0,Comm) ;
	  }
	}

	// read in z coordinates
	for(int i=1;i<nP;++i) {
	  entitySet sendSet = blockSet & distribution[i] ;
	  if(sendSet != EMPTY) {
	    int sendsz = sendSet.size() ;
	    if(binary_mode) {
	      int nr = fread(&scratch[0],sizeof(double),sendsz,fp) ;
	      if(nr != sendsz) {
		cerr << "read failed reading in mesh positions" << endl ;
		fail = 1 ;
	      }
	    } else {
	      int nr = fscanread(&scratch[0],sendsz,fp) ;
	      if(nr != sendsz) {
		cerr << "read failed reading in mesh positions" << endl ;
		fail = 1 ;
	      }
	    }
	    MPI_Send(&scratch[0],sendsz,MPI_DOUBLE,i,0,Comm) ;
	  }
	}
      }

    } else {
      //other processors receive
      for(int b=0;b<num_blocks;++b) {
	int start = blockInfo[b].nodes_base ;
	int end = start+blockInfo[b].num_nodes() ;
	entitySet blockSet = interval(start,end-1) ;
	entitySet recvSet = blockSet & distribution[rank] ;
	if(recvSet != EMPTY) {
	  // this processor receiving this block
	  int recvsz = recvSet.size() ;
	  MPI_Status stat ;

	  // receive x coordinates
	  MPI_Recv(&scratch[0],
		   recvsz,MPI_DOUBLE,0,0,Comm,&stat) ;
	  FORALL(recvSet,ii) {
	    int loc = ii-distribution[rank].Min() ;
	    int sloc = ii-recvSet.Min() ;
	    positions[loc].x = scratch[sloc] ;
	  } ENDFORALL ;
	  // receive y coordinates
	  MPI_Recv(&scratch[0],recvsz,MPI_DOUBLE,0,0,Comm,&stat) ;
	  FORALL(recvSet,ii) {
	    int loc = ii-distribution[rank].Min() ;
	    int sloc = ii-recvSet.Min() ;
	    positions[loc].y = scratch[sloc] ;
	  } ENDFORALL ;
	  // receive z coordinates
	  MPI_Recv(&scratch[0],recvsz,MPI_DOUBLE,0,0,Comm,&stat) ;
	  FORALL(recvSet,ii) {
	    int loc = ii-distribution[rank].Min() ;
	    int sloc = ii-recvSet.Min() ;
	    positions[loc].z = scratch[sloc] ;
	  } ENDFORALL ;
	}
      }
    }
  }

  MPI_Bcast(&fail, 1, MPI_INT,0,Comm) ;
  if(fail == 1)
    return false ;
  return true ;
}

void fixupGluedNodesFace(vector<Loci::Array<int,4> > &fd,
                         const map<int,int>  &gmap) {

  map<int,int>::const_iterator mi ;
  int fdsz = fd.size() ;
  int dgen_face = 0 ;
  //  int remap_cnt = 0 ;
  for(int i=0;i<fdsz;++i) {
    bool remap = false ;
    for(int j=0;j<4;++j) {
      if((mi=gmap.find(fd[i][j])) != gmap.end()) {
        fd[i][j] = mi->second ;
        remap=true ;
      }
    }
    if(remap) {
      //      remap_cnt++ ;
      // Check to see if the face became a triangle, or degenerated
      int k=-1 ;
      int cnt = 0 ;
      if(fd[i][0] == fd[i][2] || fd[i][1] == fd[i][3]) {
	//	cout << "diagonal fold" << endl ;
	cnt += 2 ;
      }
      for(int j=0;j<3;++j) {
        if(fd[i][j] == fd[i][j+1]) {
          k = j ;
          cnt++ ;
        }
      }
      if(fd[i][0] == fd[i][3]) {
        if(k< 0)
          k=4 ;
        cnt++ ;
      }
      if(k != -1) {
        for(int j=k+1;j<3;++j)
          fd[i][j] = fd[i][j+1] ;
        fd[i][3] = -1 ;
      }
      if(cnt > 1) {
        // degenerated face!
        dgen_face++ ;
        for(int j=0;j<4;++j)
          fd[i][j] = -1 ;
      }


    }
  }
  if(dgen_face > 0) {
    cout << "mesh has " << dgen_face << " degenerate faces." << endl ;
  }
}

inline int findOffset(const vector<pair<int,int> > &glueSet,
                      int i ) {
  if(i==-1)
    return 0 ;

  // prevent out of bounds access if an empty glueSet is passed.
  if(glueSet.size() < 2)
    return 0 ;
  static int hint = 0 ;
  if(i <= glueSet[0].first)
    return glueSet[0].second ;
  if(i > glueSet[hint].first && i <= glueSet[hint+1].first)
    return glueSet[hint+1].second ;
  int split = hint ;
  int low = 0 ;
  int hi = split ;
  if(i > glueSet[split].first) {
    low = split ;
    hi = glueSet.size()-2 ;
  }
  for(;;) {
    if(i > glueSet[split].first && i <= glueSet[split+1].first) {
      hint =  split ;
      return glueSet[split+1].second ;
    }
    if(i> glueSet[split].first) {
      low = split ;
      split = (low+hi+1)/2 ;
    } else {
      hi = split ;
      split = (low+hi)/2 ;
    }
  }
}

void compressGluedNodesFaces(vector<Loci::Array<int,4> > &fd,
                        const vector<pair<int,int> > &glueSet) {
  int fdsz = fd.size() ;
  for(int i=0;i<fdsz;++i) {
    for(int j=0;j<4;++j) {
      fd[i][j] += findOffset(glueSet,fd[i][j]) ;
    }
  }
}


void  compressPositions(vector<vect3d> &positions, vector<int> &pos_sizes,
                        const vector<pair<int,int> > &glueSet,
                        MPI_Comm comm) {

  // Compress positions vector to remove elements in glueSet
  int nP ;
  MPI_Comm_size(comm,&nP) ;
  int rank ;
  MPI_Comm_rank(comm,&rank) ;
  if(nP == 1) { // Single processor case
    int gcnt = 0 ;
    int npos = positions.size() ;
    for(int i=0;i<npos;++i) {
      if(i == glueSet[gcnt].first) {
        gcnt++ ;
      } else {
        if(gcnt != 0)
          positions[i-gcnt] = positions[i] ;
      }
    }
    if(gcnt != 0) {
      positions.resize(npos-gcnt) ;
    }
    pos_sizes[0] = npos-gcnt ;
  } else {
    int nstart = 0 ;
    for(int i=0;i<rank;++i)
      nstart += pos_sizes[i] ;
    int gcnt = 0 ;
    int npos = positions.size() ;
    int goffset = 0 ;
    for(goffset=0;(nstart>glueSet[goffset].first);++goffset)
    /*NULL Statement*/ ;

    for(int i=0;i<npos;++i) {
      if(nstart+i == glueSet[gcnt+goffset].first) {
        gcnt++ ;
      } else {
        if(gcnt != 0)
          positions[i-gcnt] = positions[i] ;
      }
    }
    if(gcnt != 0) {
      positions.resize(npos-gcnt) ;
    }

    npos = npos-gcnt ;

    if(rank == 0) {
      int rsize = 0 ;
      MPI_Status stat ;
      MPI_Recv(&rsize,1, MPI_INT,1,0,MPI_COMM_WORLD,&stat) ;
      positions = vector<vect3d>(rsize) ;
      MPI_Recv(&positions[0],3*rsize,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&stat) ;
      npos = rsize ;
    }

    if(rank == 1) {
      int rsize = npos/2 ;
      MPI_Send(&rsize,1,MPI_INT,0,0,MPI_COMM_WORLD) ;
      MPI_Send(&positions[0],3*rsize,MPI_DOUBLE,0,0,MPI_COMM_WORLD) ;
      int rem = npos -rsize ;
      vector<vect3d> tmp(rem);
      for(int i=0;i<rem;++i)
        tmp[i] = positions[i+rsize] ;
      positions.swap(tmp) ;
      npos = rem ;
    }

    MPI_Allgather(&npos,1,MPI_INT,&pos_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;

  }
}

bool glueCompare(const Loci::Array<int,5> &a1,
               const Loci::Array<int,5> &a2) {
  for(int i=0;i<5;++i) {
    if(a1[i] < a2[i]) // first one less, then true
      return true ;
    if(a1[i] > a2[i]) // first one greater, then false
      return false ;
  }
  return false ; // two strings equal
}

void setupGlueFaces(vector<Loci::Array<int,4> > &glueFaces,
                    vector<int> &glueBC,
                    vector<int> &glueCell,
                    vector<Loci::Array<int,4> > &noglueFaces,
                    vector<int> &noglueBC,
                    vector<int> &noglueCell)
{
  // First remove any degenerated faces from the list
  int gsz = glueFaces.size() ;

  int dgen_faces = 0 ;
  while((dgen_faces<gsz) && (glueFaces[gsz-1-dgen_faces][0]==-1))
    dgen_faces++ ;
  for(int i=0;i<gsz;++i)
    if((i < gsz-1-dgen_faces) && (glueFaces[i][0] == -1)) {
      swap(glueFaces[i],glueFaces[gsz-1-dgen_faces]) ;
      swap(glueBC[i],glueBC[gsz-1-dgen_faces]) ;
      swap(glueCell[i],glueCell[gsz-1-dgen_faces]) ;
      while((dgen_faces<gsz) && (glueFaces[gsz-1-dgen_faces][0]==-1))
	dgen_faces++ ;
    }

  gsz -= dgen_faces ;
  glueFaces.resize(gsz) ;
  glueBC.resize(gsz) ;
  glueCell.resize(gsz) ;

  int ngsz = noglueFaces.size() ;
  dgen_faces = 0 ;
  while((dgen_faces<ngsz) && (noglueFaces[ngsz-1-dgen_faces][0]==-1))
    dgen_faces++ ;
  for(int i=0;i<ngsz;++i)
    if((i < ngsz-
        dgen_faces) && (noglueFaces[i][0] < -1)) {
      swap(noglueFaces[i],noglueFaces[ngsz-1-dgen_faces]) ;
      swap(noglueBC[i],noglueBC[ngsz-1-dgen_faces]) ;
      swap(noglueCell[i],noglueCell[ngsz-1-dgen_faces]) ;
      dgen_faces++ ;
    }
  ngsz -= dgen_faces ;
  noglueFaces.resize(ngsz) ;
  noglueBC.resize(ngsz) ;
  noglueCell.resize(ngsz) ;

  vector<Loci::Array<int,5> > glueTmp(gsz) ;

  for(int i=0;i<gsz;++i) {
    for(int j=0;j<4;++j)
      glueTmp[i][j] = glueFaces[i][j] ;
    glueTmp[i][4] = i ;
    // Sort the face (using simple sort)
    for(int j=0;j<4;++j)
      for(int k=j+1;k<4;++k)
        if(glueTmp[i][j] > glueTmp[i][k])
          swap(glueTmp[i][j],glueTmp[i][k]) ;
  }
  sort(glueTmp.begin(),glueTmp.end(),glueCompare) ;
  int gluecnt = 0 ;
  int glueerror = 0 ;
  Loci::entitySet errorSet ;
  vector<Loci::Array<int,4> > fset ;
  vector<int> cl,cr ;
  for(int i=0;i<gsz;++i) {
    if(i < gsz-1 &&
       glueTmp[i][0] == glueTmp[i+1][0] &&
       glueTmp[i][1] == glueTmp[i+1][1] &&
       glueTmp[i][2] == glueTmp[i+1][2] &&
       glueTmp[i][3] == glueTmp[i+1][3]) {

      int f1 = glueTmp[i][4] ;
      int f2 = glueTmp[i+1][4] ;

      if(glueCell[f1] < glueCell[f2]) {
	swap(f1,f2) ;
      }
      fset.push_back(glueFaces[f1]) ;
      cl.push_back(glueCell[f1]) ;
      cr.push_back(glueCell[f2]) ;

      ++i ;
      ++gluecnt ;
    } else {
      int fc = glueTmp[i][4] ;
      if(glueBC[fc] < 0) {
	errorSet += glueBC[fc] ;
        glueerror++ ;
      } else {
        noglueFaces.push_back(glueFaces[fc]) ;
        noglueCell.push_back(glueCell[fc]) ;
        noglueBC.push_back(glueBC[fc]) ;
      }

    }
  }
  if(glueerror > 0) {
    cerr << "Faces not identified in the bc specification file were not glued"
         << endl ;
    cerr << "Problem Blocks:" << endl ;
    for(Loci::entitySet::const_iterator ei=errorSet.begin();
	ei != errorSet.end();++ei) {
      int b = -(*ei+1)/6 ;
      int f = -(*ei+1)%6 ;
      cerr << "block = " << b+1 << ", face = " << bcnames[f] << endl ;
    }
    Loci::Abort() ;
  }
  cout << "glued " << gluecnt << " faces." << endl ;
  glueFaces = fset ;
  glueCell = cl ;
  glueBC = cr ;
}

///  creaseGroup(noglueFaces,noglueBC,positions,pos_sizes,creaseThreshold) ;

void creaseGroup(const vector<Loci::Array<int,4> > &glueFaces,
		 vector<int> &glueBC,
		 const vector<vect3d> &positions,
		 const vector<int> &pos_sizes,
		 double creaseThreshold)  {
  // compute face edges
  vector<pair<int,int> > facepairs ; // faces that share an edge
  map<pair<int,int>,double> eweight ; // weights between boundary creases
  vector<int> nodeSet ; // nodes needed to compute edge normals
  if(Loci::MPI_rank == 0) {
    vector<pair<pair<int,int>,int> > edge_list ;
    int nfaces = glueFaces.size() ;
    for(int i=0;i<nfaces;++i) {
      if(glueFaces[i][3] == -1) {
	pair<int,int> e1(min(glueFaces[i][0],glueFaces[i][1]),
			 max(glueFaces[i][0],glueFaces[i][1])) ;
	pair<int,int> e2(min(glueFaces[i][1],glueFaces[i][2]),
			 max(glueFaces[i][1],glueFaces[i][2])) ;
	pair<int,int> e3(min(glueFaces[i][2],glueFaces[i][0]),
			 max(glueFaces[i][2],glueFaces[i][0])) ;
	edge_list.push_back(pair<pair<int,int>,int>(e1,i)) ;
	edge_list.push_back(pair<pair<int,int>,int>(e2,i)) ;
	edge_list.push_back(pair<pair<int,int>,int>(e3,i)) ;
      } else {
	pair<int,int> e1(min(glueFaces[i][0],glueFaces[i][1]),
			 max(glueFaces[i][0],glueFaces[i][1])) ;
	pair<int,int> e2(min(glueFaces[i][1],glueFaces[i][2]),
			 max(glueFaces[i][1],glueFaces[i][2])) ;
	pair<int,int> e3(min(glueFaces[i][2],glueFaces[i][3]),
			 max(glueFaces[i][2],glueFaces[i][3])) ;
	pair<int,int> e4(min(glueFaces[i][3],glueFaces[i][0]),
			 max(glueFaces[i][3],glueFaces[i][0])) ;
	edge_list.push_back(pair<pair<int,int>,int>(e1,i)) ;
	edge_list.push_back(pair<pair<int,int>,int>(e2,i)) ;
	edge_list.push_back(pair<pair<int,int>,int>(e3,i)) ;
	edge_list.push_back(pair<pair<int,int>,int>(e4,i)) ;
      }
    }
    sort(edge_list.begin(),edge_list.end()) ;
    int esz = edge_list.size() ;
    for(int i=0;i<esz-1;++i) {
      if(edge_list[i].first == edge_list[i+1].first) {
	int f1 = edge_list[i].second ;
	int f2 = edge_list[i+1].second ;
	int bc1 = glueBC[f1] ;
	int bc2 = glueBC[f2] ;
	if(bc1 != bc2) {
	  facepairs.push_back(pair<int,int>(f1,f2)) ;
	  pair<int,int> ebc(min(bc1,bc2),max(bc1,bc2)) ;
	  eweight[ebc] = 1.0 ;
	}
	++i ;
      }
    }

    for(size_t i=0;i<facepairs.size();++i) {
      int f1 = facepairs[i].first ;
      int f2 = facepairs[i].second ;

      nodeSet.push_back(glueFaces[f1][0]) ;
      nodeSet.push_back(glueFaces[f1][1]) ;
      nodeSet.push_back(glueFaces[f1][2]) ;
      if(glueFaces[f1][3] >= 0)
	nodeSet.push_back(glueFaces[f1][3]) ;
      nodeSet.push_back(glueFaces[f2][0]) ;
      nodeSet.push_back(glueFaces[f2][1]) ;
      nodeSet.push_back(glueFaces[f2][2]) ;
      if(glueFaces[f2][3] >= 0)
	nodeSet.push_back(glueFaces[f2][3]) ;
    }
    sort(nodeSet.begin(),nodeSet.end()) ;
    vector<int>::iterator it ;
    it = unique(nodeSet.begin(),nodeSet.end()) ;
    nodeSet.resize(it-nodeSet.begin()) ;
  }

  // Now send glue node positions to processor 0
  int gsize = nodeSet.size() ;
  vector<vect3d> gluePos(gsize) ;
  MPI_Bcast(&gsize,1,MPI_INT,0,MPI_COMM_WORLD) ;
  if(Loci::MPI_rank != 0) {
    nodeSet = vector<int>(gsize) ;
  }
  MPI_Bcast(&nodeSet[0],gsize,MPI_INT,0,MPI_COMM_WORLD) ;
  if(Loci::MPI_processes == 1) {
    for(int i=0;i<gsize;++i)
      gluePos[i] = positions[nodeSet[i]] ;
  } else  {
    int possz = positions.size() ;
    vector<int> dist(Loci::MPI_processes,0) ;
    MPI_Allgather(&possz,1,MPI_INT,&dist[0],1,MPI_INT,MPI_COMM_WORLD) ;
    if(Loci::MPI_rank == 0) {
      int cnt = 0 ;
      int psz = positions.size() ;
      for(size_t i=0;i<nodeSet.size() && nodeSet[i] < psz;++i) {
	gluePos[cnt] = positions[nodeSet[i]] ;
	cnt++ ;
      }
      for(int i=1;i<Loci::MPI_processes;++i) {
	MPI_Send(&cnt,1,MPI_INT,i,0,MPI_COMM_WORLD) ;
	int sz = 0 ;
	MPI_Status stat ;
	MPI_Recv(&sz,1,MPI_INT,i,0,MPI_COMM_WORLD,&stat) ;
	MPI_Recv(&gluePos[cnt],sz*3,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&stat) ;
	cnt += sz ;
      }
      if(cnt != int(gluePos.size())) {
	cerr << "problem gathering positions for glue" << endl ;
      }
    } else {
      MPI_Status stat ;
      int nstrt = 0 ;
      for(int i=0;i<Loci::MPI_rank;++i)
	nstrt += dist[i] ;
      int nend = nstrt+dist[Loci::MPI_rank] ;
      int fn = -1,ln = -1 ;
      for(int i=0;i<gsize;++i) {
	if(nodeSet[i] < nstrt)
	  continue ;
	if(nodeSet[i] >= nend) {
	  ln = i ;
	  break ;
	}
	if(fn == -1)
	  fn = i ;
      }
      if(ln == -1)
	ln = nodeSet.size() ;
      int size = 0 ;
      if(fn >= 0)
	size = ln-fn ;
      vector<vect3d> sendPos(size) ;
      for(int i=0;i<size;++i)
	sendPos[i] = positions[nodeSet[fn+i]-nstrt] ;

      int flag = 0;
      MPI_Recv(&flag,1,MPI_INT,0,0,MPI_COMM_WORLD,&stat) ;
      MPI_Send(&size,1,MPI_INT,0,0,MPI_COMM_WORLD) ;
      MPI_Send(&sendPos[0],size*3,MPI_DOUBLE,0,0,MPI_COMM_WORLD) ;
    }
  }

  if(Loci::MPI_rank == 0) {
    map<int,int> g2l ;
    for(int i=0;i<gsize;++i)
      g2l[nodeSet[i]] = i ;

    int maxbc = 0 ;
    for(size_t i=0;i<facepairs.size();++i) {
      int f1 = facepairs[i].first ;
      int f2 = facepairs[i].second ;
      vect3d p10 = gluePos[g2l[glueFaces[f1][0]]] ;
      vect3d p11 = gluePos[g2l[glueFaces[f1][1]]] ;
      vect3d p12 = gluePos[g2l[glueFaces[f1][2]]] ;
      bool t1 = (glueFaces[f1][3] == -1) ;
      vect3d p13 = p10 ;
      if(!t1)
	p13 = gluePos[g2l[glueFaces[f1][3]]] ;
      vect3d n1 = cross(p12-p10,p11-p13) ;
      n1 *= 1./max(norm(n1),1e-30) ;

      vect3d p20 = gluePos[g2l[glueFaces[f2][0]]] ;
      vect3d p21 = gluePos[g2l[glueFaces[f2][1]]] ;
      vect3d p22 = gluePos[g2l[glueFaces[f2][2]]] ;
      bool t2 = (glueFaces[f2][3] == -1) ;
      vect3d p23 = p10 ;
      if(!t2)
	p23 = gluePos[g2l[glueFaces[f2][3]]] ;
      vect3d n2 = cross(p22-p20,p21-p23) ;
      n2 *= 1./max(norm(n2),1e-30) ;
      int bc1 = glueBC[f1] ;
      int bc2 = glueBC[f2] ;
      maxbc = max(maxbc,max(bc1,bc2)) ;

      pair<int,int> ebc(min(bc1,bc2),max(bc1,bc2)) ;
      eweight[ebc] = min(eweight[ebc],dot(n1,n2)) ;
    }

    double th = cos(creaseThreshold*0.0174532927778) ;
    vector<entitySet> equal(maxbc+1) ;

    map<pair<int,int>,double>::const_iterator mi ;

    for(mi=eweight.begin();mi!=eweight.end();++mi) {
      if(mi->second > th) {
	equal[mi->first.first] += mi->first.first ;
	equal[mi->first.first] += mi->first.second ;
	equal[mi->first.second] += mi->first.first ;
	equal[mi->first.second] += mi->first.second ;
      }
    }
    // now perform transitive closure of equivalence sets
    bool finished = true ;
    do {
      finished = true ;
      for(int i=0;i<maxbc+1;++i) {
	if(equal[i] != EMPTY) {
	  entitySet tmp = equal[i] ;
	  tmp += i ; // Add self
	  entitySet all = tmp ; // Gather all equivalences
	  FORALL(tmp,j) {
	    all += equal[j] ;
	  } ENDFORALL ;
	  FORALL(all,j) { // Check to see if we converged
	    if(equal[j] != all)
	      finished = false ;
	    equal[j] += all ;
	  } ENDFORALL ;
	}
      }
    } while(!finished) ;


    vector<int> map_bc(maxbc+1) ;
    for(int i=0;i<maxbc+1;++i) {
      if(equal[i] == EMPTY)
	map_bc[i] = i ;
      else {
	map_bc[i] = equal[i].Min() ;
      }
    }
    for(size_t i=0;i<glueBC.size();++i)
      glueBC[i] = map_bc[glueBC[i]] ;
  }


}

int main(int ac, char* av[]) {
  using namespace Loci ;
  using namespace VOG ;

  bool optimize = true ;
  bool lefthanded = false ;
  double creaseThreshold = 5 ;
  Loci::Init(&ac,&av) ;

  if(ac == 1) {
    usage() ;
    exit(1) ;
  }

  bool boundary_file = false ;
  string boundary_filename ;

  double tol = 1e-3 ;
  string Lref = "NOSCALE" ;
  double max_aspect = 1e6 ;
  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 2 && !strcmp(av[1],"-lefthanded")) {
      lefthanded = true ;
      ac-- ;
      av++ ;
    } else if(ac >= 3 && !strcmp(av[1],"-tol")) {
      tol = atof(av[2]) ;
      if(tol >= 1.0) {
	cerr << "tolerance cannot be greater than 1" << endl ;
	tol = 0.9;
      }
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 3 && !strcmp(av[1],"-maxAspect")) {
      max_aspect = atof(av[2]) ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 3 && !strcmp(av[1],"-creaseAngle")) {
      creaseThreshold = atof(av[2]) ;
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
    } else if (ac >= 2 && !strcmp(av[1],"-bc")) {
      boundary_file = true ;
      boundary_filename = string(av[2]) ;
      ac -= 2 ;
      av += 2 ;
    }else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      ac-- ;
      av++ ;
    }
  }

  if(Lref == "NOSCALE") {
    if(Loci::MPI_rank == 0)
      cerr << "Must set grid units!" << endl
	   << "Use options -in, -ft, -cm, -m, or -Lref to set grid units." << endl ;
    Loci::Abort() ;
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

  MPI_Bcast(&posScale,1,MPI_DOUBLE,0,MPI_COMM_WORLD) ;
  int read_type = 0 ;
  bool found_file = false ;
  string input_name = av[1] ;
  bool explicit_ascii = VOG::hasSuffix(input_name,".grd") ;
  bool explicit_binary = VOG::hasSuffix(input_name,".grd.b8") ;
  string case_name = input_name ;
  if(explicit_binary)
    case_name = VOG::stripSuffix(input_name,".grd.b8") ;
  else if(explicit_ascii)
    case_name = VOG::stripSuffix(input_name,".grd") ;
  string file = input_name ;
  if(Loci::MPI_rank == 0) {
    struct stat finfo ;
    if(explicit_ascii || explicit_binary) {
      if(stat(file.c_str(),&finfo)==0 && (S_ISREG(finfo.st_mode)))
        found_file = true ;
      if(explicit_binary)
        read_type = 1 ;
      if(!found_file) {
        cerr << "unable to find file '"<< file << "'" << endl ;
        Loci::Abort() ;
      }
    } else {
      bool found_ascii=false ;
      bool found_binary = false ;
      string ascii_file = input_name + ".grd" ;
      string binary_file = input_name + ".grd.b8" ;
      if(stat(ascii_file.c_str(),&finfo)==0 && (S_ISREG(finfo.st_mode))) {
        found_ascii = true ;
        found_file = true ;
      }
      if(stat(binary_file.c_str(),&finfo)==0 && (S_ISREG(finfo.st_mode))) {
        found_binary = true ;
        found_file = true ;
      }
      if(found_binary && found_ascii) {
        cerr << "both ascii file, '" << ascii_file << "' and binary file, '"
             << binary_file << "', found.  Remove one file to disambiguate."
             << endl ;
        Loci::Abort() ;
      }
      if(found_binary) {
        file = binary_file ;
        read_type = 1 ;
      } else {
        file = ascii_file ;
      }
      if(!found_file) {
        cerr << "unable to find file '"<< ascii_file << "' or '" << binary_file
             << "'" << endl ;
        Loci::Abort() ;
      }
    }
  }

  vector<block_topo> blockInfo ;
  vector<vect3d> positions ;
  vector<int> pos_sizes ;
  if(!readP3DGrid(file,blockInfo,positions,pos_sizes,read_type,MPI_COMM_WORLD)) {
    cerr << "error reading grid '" << file << "'" << endl ;
    Loci::Abort() ;
  }

  if(blockInfo.size() == 0) {
    cerr << "error reading grid '" << file << "'" << endl ;
    Loci::Abort() ;
  }
  // Check ni,nj,nk
  int minni=blockInfo[0].ni ;
  int minnj=blockInfo[0].nj ;
  int minnk=blockInfo[0].nk ;

  for(size_t b = 1 ; b!= blockInfo.size();++b) {
    minni = min(minni,blockInfo[b].ni) ;
    minnj = min(minni,blockInfo[b].nj) ;
    minnk = min(minni,blockInfo[b].nk) ;
  }
  if(minni == 1 || minnj == 1)  {
    cerr << "unable to process plot3d file with only ni==1 or nj == 1" << endl ;
  }
  if(minnk == 1) {
    cout << "NOTE:  Extruding grid in z direction to have at least two planes" << endl ;
    if(Loci::MPI_processes > 1) {
      cerr << "extruding only supported in serial mode!" << endl ;
      Loci::Abort() ;
    }
    vector<vect3d> pos2(positions.size()*2) ;
    pos_sizes[0] *= 2 ;
    int noffset = 0 ;
    int noffset2 = 0 ;
    for(size_t b = 0 ; b!= blockInfo.size();++b) {
      int nnodes = blockInfo[b].num_nodes() ;
      for(int i=0;i<nnodes;++i)
	pos2[noffset2+i] = positions[noffset+i] ;
      vect3d delta(0.,0.,0.01) ;
      for(int i=0;i<nnodes;++i)
	pos2[noffset2+i+nnodes] = positions[noffset+i]+delta ;
      noffset2 += nnodes*2 ;
      noffset += nnodes ;
      if(blockInfo[b].nk != 1) {
	cerr << "all blocks must have nk=1 for the extrusion to work!" << endl ;
	Loci::Abort() ;
      }
      blockInfo[b].nk = 2 ;
    }
    positions.swap(pos2) ;
    // Now update block offsets
    int ntot = 0 ;
    for(size_t b = 0 ; b!= blockInfo.size();++b) {
      blockInfo[b].nodes_base = ntot ;
      ntot += blockInfo[b].num_nodes() ;
    }
    int ctot = 0 ;
    for(size_t b = 0 ; b!= blockInfo.size();++b) {
      blockInfo[b].cells_base = ctot ;
      ctot += blockInfo[b].num_cells() ;
    }
    int ftot = 0 ;
    for(size_t b = 0 ; b!= blockInfo.size();++b) {
      blockInfo[b].faces_base = ftot ;
      ftot += blockInfo[b].num_faces() ;
    }
  }


  // Scale the grid
  for(size_t i=0;i<positions.size();++i)
    positions[i] *= posScale ;

  // Read in the bc file
  vector<block_boundary_info> bface_info (blockInfo.size()) ;
  for(size_t b = 0 ; b!= blockInfo.size();++b) {
    for(int i=0;i<6;++i) {
      block_face_info tmp;
      tmp.tag = b*6+i+1 ;
      tmp.is = -2 ;
      tmp.ie = -2 ;
      tmp.js = -2 ;
      tmp.je = -2 ;
      bface_info[b].block_face[i].push_back(tmp) ;
    }
  }

  map<int, string> bcnamelist ;
  vector<int> glueNodes ;
  vector<Loci::Array<int,4> > glueFaces ;
  vector<int> glueCell,glueBC ;
  vector<Loci::Array<int,4> > noglueFaces ;
  vector<int> noglueCell,noglueBC ;

  if(Loci::MPI_rank == 0) {
    if(boundary_file) {
      for(size_t b = 0 ; b!= blockInfo.size();++b) {
        for(int i=0;i<6;++i) {
	  bface_info[b].block_face[i].clear() ;
        }
      }
      ifstream bfile(boundary_filename.c_str(),ios::in) ;
      if(bfile.fail() || bfile.eof()) {
        cerr << "error opening boundary file " << boundary_filename << endl ;
        Loci::Abort() ;
      }
      int numbcs = 0;
      bfile >> numbcs ;
      Loci::parse::kill_white_space(bfile) ;

      for(int i=0;i<numbcs;++i) {
        if(bfile.fail() || bfile.eof()) {
          cerr << "error reading boundary file " << boundary_filename << endl ;
          exit(-1) ;
        }
        int boundary_flag = 0 ;
        bfile >> boundary_flag ;
        int boundary_segments = 0 ;
        bfile >> boundary_segments ;

        string boundary_name = Loci::parse::get_name(bfile) ;
        bcnamelist[boundary_flag] = boundary_name ;

        for(int j=0;j<boundary_segments;++j) {
          int block, face, Is,Ie,Js,Je ;
          bfile >> block ;
          block-- ;
          Loci::parse::kill_white_space(bfile) ;
          string fname = Loci::parse::get_name(bfile) ;
          if(block < 0 || block > int(blockInfo.size())) {
            cerr << "block id out of bounds for face name '" << fname << "'"
                 << endl ;
          }
          face = -1 ;
          for(int f=0;f<6;++f)
            if(fname == bcnames[f])
              face = f ;
          if(face == -1) {
            cerr << "unknown face type " << fname
               << " in file " << boundary_filename
                 << endl ;
            exit(-1) ;
          }
          Loci::parse::kill_white_space(bfile) ;
          bfile >> Is >> Ie >> Js >> Je ;

          if(Is == -1)
            Ie = 1000000000 ;
          if(Js == -1)
            Je = 1000000000 ;
	  block_face_info tmp;
	  tmp.tag = boundary_flag ;
	  tmp.is = Is ;
	  tmp.ie = Ie ;
	  tmp.js = Js ;
	  tmp.je = Je ;

          bface_info[block].block_face[face].push_back(tmp);
        }
      }
      for(size_t b = 0 ; b!= blockInfo.size();++b) {
	for(int i=0;i<6;++i) {
	  if(bface_info[b].block_face[i].empty()) {
	    block_face_info tmp;
	    tmp.tag = -(b*6+i+1) ;
	    tmp.is = -2 ;
	    tmp.ie = -2 ;
	    tmp.js = -2 ;
	    tmp.je = -2 ;
	    bface_info[b].block_face[i].push_back(tmp) ;
	  }
	}
      }
    }

    // Now find nodes to glue
    for(size_t b = 0 ; b!= blockInfo.size();++b) {
      const int ni = blockInfo[b].ni ;
      const int nj = blockInfo[b].nj ;
      const int nk = blockInfo[b].nk ;
      // IJ1 face and IJN face
      for(int i=0;i<ni-1;++i)
        for(int j=0;j<nj-1;++j) {
          bool noglue0 = (i+1 >= bface_info[b].block_face[0][0].is &&
                          i+1 <  bface_info[b].block_face[0][0].ie &&
                          j+1 >= bface_info[b].block_face[0][0].js &&
                          j+1 <  bface_info[b].block_face[0][0].je) ;
	  int bc0 = bface_info[b].block_face[0][0].tag ;
	  for(size_t l=1;l<bface_info[b].block_face[0].size();++l) {
	    bool ng = (i+1 >= bface_info[b].block_face[0][l].is &&
		       i+1 <  bface_info[b].block_face[0][l].ie &&
		       j+1 >= bface_info[b].block_face[0][l].js &&
		       j+1 <  bface_info[b].block_face[0][l].je) ;
	    if(ng) {
	      noglue0 = ng ;
	      bc0 = bface_info[b].block_face[0][l].tag ;
	    }
	  }
          bool noglue1 = (i+1 >= bface_info[b].block_face[1][0].is &&
                          i+1 <  bface_info[b].block_face[1][0].ie &&
                          j+1 >= bface_info[b].block_face[1][0].js &&
                          j+1 <  bface_info[b].block_face[1][0].je) ;
	  int bc1 = bface_info[b].block_face[1][0].tag ;
	  for(size_t l=1;l<bface_info[b].block_face[1].size();++l) {
	    bool ng = (i+1 >= bface_info[b].block_face[1][l].is &&
		       i+1 <  bface_info[b].block_face[1][l].ie &&
		       j+1 >= bface_info[b].block_face[1][l].js &&
		       j+1 <  bface_info[b].block_face[1][l].je) ;
	    if(ng) {
	      noglue1 = ng ;
	      bc1 = bface_info[b].block_face[1][l].tag ;
	    }
	  }

	  Loci::Array<int,4> face ;
	  face[3] = blockInfo[b].naddr(i  ,j  ,0) ;
	  face[2] = blockInfo[b].naddr(i+1,j  ,0) ;
	  face[1] = blockInfo[b].naddr(i+1,j+1,0) ;
	  face[0] = blockInfo[b].naddr(i  ,j+1,0) ;
	  int cell = blockInfo[b].caddr(i,j,0) ;
          if(noglue0) {
	    noglueFaces.push_back(face) ;
	    noglueCell.push_back(cell) ;
	    noglueBC.push_back(bc0) ;
	  } else {
	    glueFaces.push_back(face) ;
	    glueCell.push_back(cell) ;
	    glueBC.push_back(bc0) ;
	    for(int i=0;i<4;++i)
	      glueNodes.push_back(face[i]) ;
	  }

	  face[0] = blockInfo[b].naddr(i  ,j  ,nk-1) ;
	  face[1] = blockInfo[b].naddr(i+1,j  ,nk-1) ;
	  face[2] = blockInfo[b].naddr(i+1,j+1,nk-1) ;
	  face[3] = blockInfo[b].naddr(i  ,j+1,nk-1) ;
	  cell = blockInfo[b].caddr(i,j,nk-2) ;
          if(noglue1) {
	    noglueFaces.push_back(face) ;
	    noglueCell.push_back(cell) ;
	    noglueBC.push_back(bc1) ;
	  } else {
	    glueFaces.push_back(face) ;
	    glueCell.push_back(cell) ;
	    glueBC.push_back(bc1) ;
	    for(int i=0;i<4;++i)
	      glueNodes.push_back(face[i]) ;
	  }
        }
      // IK1 face and IKN face
      for(int i=0;i<ni-1;++i)
        for(int k=0;k<nk-1;++k) {
          bool noglue0 = (i+1 >= bface_info[b].block_face[2][0].is &&
                          i+1 <  bface_info[b].block_face[2][0].ie &&
                          k+1 >= bface_info[b].block_face[2][0].js &&
                          k+1 <  bface_info[b].block_face[2][0].je) ;
	  int bc0 = bface_info[b].block_face[2][0].tag ;
	  for(size_t l=1;l<bface_info[b].block_face[2].size();++l) {
	    int ng = (i+1 >= bface_info[b].block_face[2][l].is &&
		      i+1 <  bface_info[b].block_face[2][l].ie &&
		      k+1 >= bface_info[b].block_face[2][l].js &&
		      k+1 <  bface_info[b].block_face[2][l].je) ;
	    if(ng) {
	      noglue0 = ng ;
	      bc0 = bface_info[b].block_face[2][l].tag ;
	    }
	  }

          bool noglue1 = (i+1 >= bface_info[b].block_face[3][0].is &&
                          i+1 <  bface_info[b].block_face[3][0].ie &&
                          k+1 >= bface_info[b].block_face[3][0].js &&
                          k+1 <  bface_info[b].block_face[3][0].je) ;
	  int bc1 = bface_info[b].block_face[3][0].tag ;
	  for(size_t l=1;l<bface_info[b].block_face[3].size();++l) {
	    bool ng = (i+1 >= bface_info[b].block_face[3][l].is &&
		       i+1 <  bface_info[b].block_face[3][l].ie &&
		       k+1 >= bface_info[b].block_face[3][l].js &&
		       k+1 <  bface_info[b].block_face[3][l].je) ;
	    if(ng) {
	      noglue1 = ng ;
	      bc1 = bface_info[b].block_face[3][l].tag ;
	    }
	  }
	  Loci::Array<int,4> face ;
	  face[0] = blockInfo[b].naddr(i  ,0 ,k) ;
	  face[1] = blockInfo[b].naddr(i+1,0 ,k) ;
	  face[2] = blockInfo[b].naddr(i+1,0 ,k+1) ;
	  face[3] = blockInfo[b].naddr(i  ,0 ,k+1) ;
	  int cell = blockInfo[b].caddr(i,0,k) ;
          if(noglue0) {
	    noglueFaces.push_back(face) ;
	    noglueCell.push_back(cell) ;
	    noglueBC.push_back(bc0) ;
	  } else {
	    glueFaces.push_back(face) ;
	    glueCell.push_back(cell) ;
	    glueBC.push_back(bc0) ;
	    for(int i=0;i<4;++i)
	      glueNodes.push_back(face[i]) ;
	  }

	  face[3] = blockInfo[b].naddr(i  ,nj-1 ,k) ;
	  face[2] = blockInfo[b].naddr(i+1,nj-1 ,k) ;
	  face[1] = blockInfo[b].naddr(i+1,nj-1 ,k+1) ;
	  face[0] = blockInfo[b].naddr(i  ,nj-1 ,k+1) ;
	  cell = blockInfo[b].caddr(i,nj-2,k) ;
          if(noglue1) {
	    noglueFaces.push_back(face) ;
	    noglueCell.push_back(cell) ;
	    noglueBC.push_back(bc1) ;
	  } else {
	    glueFaces.push_back(face) ;
	    glueCell.push_back(cell) ;
	    glueBC.push_back(bc1) ;
	    for(int i=0;i<4;++i)
	      glueNodes.push_back(face[i]) ;
	  }
        }
      // JK1 face and JKN face
      for(int j=0;j<nj-1;++j)
        for(int k=0;k<nk-1;++k) {
          bool noglue0 = (j+1 >= bface_info[b].block_face[4][0].is &&
                          j+1 <  bface_info[b].block_face[4][0].ie &&
                          k+1 >= bface_info[b].block_face[4][0].js &&
                          k+1 <  bface_info[b].block_face[4][0].je) ;
	  int bc0 = bface_info[b].block_face[4][0].tag ;
	  for(size_t l=1;l<bface_info[b].block_face[4].size();++l) {
	    bool ng = (j+1 >= bface_info[b].block_face[4][l].is &&
		       j+1 <  bface_info[b].block_face[4][l].ie &&
		       k+1 >= bface_info[b].block_face[4][l].js &&
		       k+1 <  bface_info[b].block_face[4][l].je) ;
	    if(ng) {
	      noglue0 = ng ;
	      bc0 = bface_info[b].block_face[4][l].tag ;
	    }
	  }

          bool noglue1 = (j+1 >= bface_info[b].block_face[5][0].is &&
                          j+1 <  bface_info[b].block_face[5][0].ie &&
                          k+1 >= bface_info[b].block_face[5][0].js &&
                          k+1 <  bface_info[b].block_face[5][0].je) ;
	  int bc1 = bface_info[b].block_face[5][0].tag ;
	  for(size_t l=1;l<bface_info[b].block_face[5].size();++l) {
	    bool ng = (j+1 >= bface_info[b].block_face[5][l].is &&
		       j+1 <  bface_info[b].block_face[5][l].ie &&
		       k+1 >= bface_info[b].block_face[5][l].js &&
		       k+1 <  bface_info[b].block_face[5][l].je) ;
	    if(ng) {
	      noglue1 = ng ;
	      bc1 = bface_info[b].block_face[5][l].tag ;
	    }
	  }
	  Loci::Array<int,4> face ;
	  face[3] = blockInfo[b].naddr(0 ,j  ,k) ;
	  face[2] = blockInfo[b].naddr(0 ,j+1,k) ;
	  face[1] = blockInfo[b].naddr(0 ,j+1,k+1) ;
	  face[0] = blockInfo[b].naddr(0 ,j  ,k+1) ;
	  int cell = blockInfo[b].caddr(0,j,k) ;
          if(noglue0) {
	    noglueFaces.push_back(face) ;
	    noglueCell.push_back(cell) ;
	    noglueBC.push_back(bc0) ;
	  } else {
	    glueFaces.push_back(face) ;
	    glueCell.push_back(cell) ;
	    glueBC.push_back(bc0) ;
	    for(int i=0;i<4;++i)
	      glueNodes.push_back(face[i]) ;
	  }

	  face[0] = blockInfo[b].naddr(ni-1 ,j  ,k) ;
	  face[1] = blockInfo[b].naddr(ni-1 ,j+1,k) ;
	  face[2] = blockInfo[b].naddr(ni-1 ,j+1,k+1) ;
	  face[3] = blockInfo[b].naddr(ni-1 ,j  ,k+1) ;
	  cell = blockInfo[b].caddr(ni-2,j,k) ;
          if(noglue1) {
	    noglueFaces.push_back(face) ;
	    noglueCell.push_back(cell) ;
	    noglueBC.push_back(bc1) ;
	  } else {
	    glueFaces.push_back(face) ;
	    glueCell.push_back(cell) ;
	    glueBC.push_back(bc1) ;
	    for(int i=0;i<4;++i)
	      glueNodes.push_back(face[i]) ;
	  }
        }
      sort(glueNodes.begin(),glueNodes.end()) ;
      vector<int>::iterator it ;
      it = unique(glueNodes.begin(),glueNodes.end()) ;
      glueNodes.resize(it-glueNodes.begin()) ;
    }
  }
  // Now send glue node positions to processor 0
  int gsize = glueNodes.size() ;
  vector<vect3d> gluePos(gsize) ;
  MPI_Bcast(&gsize,1,MPI_INT,0,MPI_COMM_WORLD) ;
  if(Loci::MPI_rank != 0) {
    glueNodes = vector<int>(gsize) ;
  }
  MPI_Bcast(&glueNodes[0],gsize,MPI_INT,0,MPI_COMM_WORLD) ;
  if(Loci::MPI_processes == 1) {
    for(int i=0;i<gsize;++i)
      gluePos[i] = positions[glueNodes[i]] ;
  } else  {
    int possz = positions.size() ;
    vector<int> dist(Loci::MPI_processes,0) ;
    MPI_Allgather(&possz,1,MPI_INT,&dist[0],1,MPI_INT,MPI_COMM_WORLD) ;
    if(Loci::MPI_rank == 0) {
      int cnt = 0 ;
      for(int i=1;i<Loci::MPI_processes;++i) {
	MPI_Send(&cnt,1,MPI_INT,i,0,MPI_COMM_WORLD) ;
	int sz = 0 ;
	MPI_Status stat ;
	MPI_Recv(&sz,1,MPI_INT,i,0,MPI_COMM_WORLD,&stat) ;
	MPI_Recv(&gluePos[cnt],sz*3,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&stat) ;
	cnt += sz ;
      }
      if(cnt != int(gluePos.size())) {
	cerr << "problem gathering positions for glue" << endl ;
      }
    } else {
      MPI_Status stat ;
      int nstrt = 0 ;
      for(int i=1;i<Loci::MPI_rank;++i)
	nstrt += dist[i] ;
      int nend = nstrt+dist[Loci::MPI_rank] ;
      int fn = -1,ln = -1 ;
      for(int i=0;i<gsize;++i) {
	if(glueNodes[i] < nstrt)
	  continue ;
	if(glueNodes[i] >= nend) {
	  ln = i ;
	  break ;
	}
	if(fn == -1)
	  fn = i ;
      }
      if(ln == -1)
	ln = glueNodes.size() ;
      int size = 0 ;
      if(fn >= 0)
	size = ln-fn ;
      vector<vect3d> sendPos(size) ;
      for(int i=0;i<size;++i)
	sendPos[i] = positions[glueNodes[fn+i]-nstrt] ;

      int flag = 0;
      MPI_Recv(&flag,1,MPI_INT,0,0,MPI_COMM_WORLD,&stat) ;
      MPI_Send(&size,1,MPI_INT,0,0,MPI_COMM_WORLD) ;
      MPI_Send(&sendPos[0],size*3,MPI_DOUBLE,0,0,MPI_COMM_WORLD) ;
    }
  }

  map<int,int> glue2local ;
  for(int i=0;i<gsize;++i)
    glue2local[glueNodes[i]] = i ;

  vector<pair<int,int> > glueSet ;
  if(Loci::MPI_rank == 0) {
    int naspect = 0 ;
    // Compute glue length
    vector<double> gluelen(gsize,1e30) ;
    int gfsz = glueFaces.size() ;
    for(int i=0;i<gfsz;++i) {
      int n1 = glue2local[glueFaces[i][0]] ;
      int n2 = glue2local[glueFaces[i][1]] ;
      int n3 = glue2local[glueFaces[i][2]] ;
      int n4 = glue2local[glueFaces[i][3]] ;
      vect3d p1 = gluePos[n1] ;
      vect3d p2 = gluePos[n2] ;
      vect3d p3 = gluePos[n3] ;
      vect3d p4 = gluePos[n4] ;
      vect3d e1 = 0.5*(p1+p2) ;
      vect3d e2 = 0.5*(p3+p4) ;
      vect3d e3 = 0.5*(p1+p4) ;
      vect3d e4 = 0.5*(p2+p3) ;
      double l1 = norm(e1-e2) ;
      double l2 = norm(e3-e4) ;

      double lenmx = max(l1,l2) ;
      double lenmn = min(l1,l2) ;
      // Aspect ratio control, no aspect ratio face over max_aspect
      double len = (lenmx> max_aspect*lenmn)?lenmx:lenmn ;
      if(lenmx > max_aspect*lenmn && lenmn != 0) {
	naspect++ ;

      }

      gluelen[n1] = min(gluelen[n1],len) ;
      gluelen[n2] = min(gluelen[n2],len) ;
      gluelen[n3] = min(gluelen[n3],len) ;
      gluelen[n4] = min(gluelen[n4],len) ;
    }
    if(naspect > 0) {
      cerr << "high-aspect ratio faces glued: " << naspect << endl ;
    }
    // Limit tolerance to reasonable limits for double precision arithmetic
    // and use this to compute the distance needed before a node is glued to
    // another node
    const double EPS = 1e-11 ;
    tol = max(tol,EPS*EPS) ;

    for(int i=0;i<gsize;++i) {
      double eps = EPS*max(fabs(gluePos[i].x),
				 max(fabs(gluePos[i].y),fabs(gluePos[i].z))) ;
      gluelen[i] = max(tol*gluelen[i],eps) ;
    }

    vector<pair<double,int> > nodeSort(gsize) ;
    vect3d dir(Ax,Ay,Az) ;
    for(int i=0;i<gsize;++i) {
      nodeSort[i].first = dot(dir,gluePos[i]) ;
      nodeSort[i].second = i ;
    }
    sort(nodeSort.begin(),nodeSort.end()) ;

    vector<entitySet> equal(gluePos.size()) ;
    for(int i=0;i<gsize;++i) {
      int ti = nodeSort[i].second ;
      double del = gluelen[ti] ;
      vect3d p1 = gluePos[ti] ;
      for(int j=i+1;j<gsize;++j) {
	if((nodeSort[j].first-nodeSort[i].first) > del)
	  break ;
	double mindel = min(del,gluelen[nodeSort[j].second]) ;
	vect3d dv = p1-gluePos[nodeSort[j].second] ;
	if(dot(dv,dv) < mindel*mindel)
	  equal[ti] += nodeSort[j].second ;
      }
    }

    // now perform transitive closure of equivalence sets
    bool finished = true ;
    do {
      finished = true ;
      for(int i=0;i<gsize;++i) {
	if(equal[i] != EMPTY) {
	  entitySet tmp = equal[i] ;
	  tmp += i ; // Add self
	  entitySet all = tmp ; // Gather all equivalences
	  FORALL(tmp,j) {
	    all += equal[j] ;
	  } ENDFORALL ;
	  FORALL(all,j) { // Check to see if we converged
	    if(equal[j] != all)
	      finished = false ;
	    equal[j] += all ;
	  } ENDFORALL ;
	}
      }
    } while(!finished) ;

    for(int i=0;i<gsize;++i) {
      if(equal[i] != EMPTY) {
	int first = -1 ;
	entitySet set = equal[i] ;
	FORALL(set,j) {
	  equal[j] = EMPTY ;
	  if(first == -1)
	    first = glueNodes[j] ;
	  else
	    glueSet.push_back(pair<int,int>(glueNodes[j],first)) ;
	} ENDFORALL ;

      }
    }
    sort(glueSet.begin(),glueSet.end()) ;
    cout << "nodes to be removed due to gluing=" << glueSet.size() << endl ;
  }
  int gsz = glueSet.size() ;
  MPI_Bcast(&gsz,1,MPI_INT,0,MPI_COMM_WORLD) ;
  if(Loci::MPI_rank != 0)
    glueSet = vector<pair<int,int> >(gsz) ;
  MPI_Bcast(&glueSet[0],gsz*2,MPI_INT,0,MPI_COMM_WORLD) ;

  // Now lets make the rest of the faces

  // Compute total interior faces
  int tot_interior_faces = 0 ;
  for(size_t b=0;b<blockInfo.size();++b) {
    tot_interior_faces += blockInfo[b].num_interior_faces() ;
  }

  // Compute face distribution
  vector<int> iface_dist(Loci::MPI_processes,0) ;
  if(Loci::MPI_processes == 1)
    iface_dist[0] = tot_interior_faces ;
  else {
    int tot = 0 ;
    for(int i=1;i<Loci::MPI_processes-1;++i) {
      iface_dist[i] = tot_interior_faces/(Loci::MPI_processes-1) ;
      tot+= iface_dist[i] ;
    }
    iface_dist[Loci::MPI_processes-1] = tot_interior_faces-tot ;
  }


  int mydist = iface_dist[Loci::MPI_rank] ;
  vector<Loci::Array<int,4> > interior_face(mydist) ;
  vector<int> interior_cl(mydist), interior_cr(mydist) ;

  int iface_start = 0 ;
  for(int i=0;i<Loci::MPI_rank;++i)
    iface_start += iface_dist[i] ;

  int iface_end = iface_start+iface_dist[Loci::MPI_rank] ;

  // Fill in the faces
  int loc = 0;
  int gf = 0;
  for(size_t b=0;b<blockInfo.size();++b) {
    int nextloc = loc + blockInfo[b].num_interior_faces() ;
    if(nextloc >= iface_start && loc < iface_end) {
      int rstart = max(loc,iface_start)-loc  ;
      int rend = min(iface_end,nextloc)-loc ;

      for(int f=rstart;f<rend;++f) {
        blockInfo[b].fill_interior_face(f,interior_face[gf],
                                        interior_cl[gf],
                                        interior_cr[gf]) ;
        gf++ ;
      }
    }
    loc = nextloc ;
  }

  // remove duplicate nodes from face data-structures
  map<int,int> gmap ;
  for(size_t i=0;i<glueSet.size();++i) {
    gmap[glueSet[i].first] = glueSet[i].second ;
    glueSet[i].second = -int(i) ;
  }
  glueSet.push_back(pair<int,int>(2147483647,-int(gsz))) ;

  fixupGluedNodesFace(interior_face,gmap) ;
  fixupGluedNodesFace(glueFaces,gmap) ;
  fixupGluedNodesFace(noglueFaces,gmap) ;

  // Now compress the nodes removing the unused nodes
  compressGluedNodesFaces(interior_face,glueSet) ;
  compressGluedNodesFaces(glueFaces,glueSet) ;
  compressGluedNodesFaces(noglueFaces,glueSet) ;
  compressPositions(positions,pos_sizes,glueSet,MPI_COMM_WORLD) ;

  if(Loci::MPI_rank == 0) {
    setupGlueFaces(glueFaces,glueBC,glueCell,noglueFaces,noglueBC,noglueCell) ;
  }

  MPI_Bcast(&creaseThreshold,1,MPI_DOUBLE,0,MPI_COMM_WORLD) ;
  if(!boundary_file && creaseThreshold > 0.0)
    creaseGroup(noglueFaces,noglueBC,positions,pos_sizes,creaseThreshold) ;

  // Now create pos store
  int r = Loci::MPI_rank ;
  int p = Loci::MPI_processes ;
  int nnodes = 0 ;
  int snodes = 0 ;
  for(int i=0;i<p;++i) {
    nnodes += pos_sizes[i] ;
    if(i < r)
      snodes += pos_sizes[i] ;
  }
  store<vect3d> pos ;
  entitySet pdom ;
  if(pos_sizes[r] > 0)
    pdom = interval(snodes,snodes+pos_sizes[r]-1) ;
  pos.allocate(pdom) ;
  for(int i=0;i<pos_sizes[r];++i)
    pos[i+snodes] = positions[i] ;

  { vector<vect3d> tmp ; positions.swap(tmp) ;  } // clear out positions

  int ncells = 0 ;
  for(size_t b=0;b<blockInfo.size();++b) {
    ncells+=blockInfo[b].num_cells() ;
  }

  int nloc_faces = glueFaces.size()+noglueFaces.size()+interior_face.size() ;
  vector<int> fsizes(p) ;
  MPI_Allgather(&nloc_faces,1,MPI_INT,&fsizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  int fstart = nnodes+ncells ;
  for(int i=0;i<r;++i)
    fstart += fsizes[i] ;
  entitySet fdom ;
  if(nloc_faces > 0)
    fdom = interval(fstart,fstart+nloc_faces-1) ;

  store<int> count ;
  count.allocate(fdom) ;
  Map cl,cr ;
  cl.allocate(fdom) ;
  cr.allocate(fdom) ;
  int fcnt = fstart ;
  for(size_t i=0;i<interior_face.size();++i) {
    cl[fcnt] = interior_cl[i]+nnodes ;
    cr[fcnt] = interior_cr[i]+nnodes ;
    count[fcnt] = 4 ;
    if(interior_face[i][3] == -1)
      count[fcnt] = 3 ;
    fcnt++ ;
  }
  for(size_t i=0;i<glueFaces.size();++i) {
    cl[fcnt] = glueCell[i]+nnodes ;
    cr[fcnt] = glueBC[i]+nnodes ;
    count[fcnt] = 4 ;
    if(glueFaces[i][3] == -1)
      count[fcnt] = 3 ;
    fcnt++ ;
  }
  entitySet bcset ;
  for(size_t i=0;i<noglueFaces.size();++i) {
    cl[fcnt] = noglueCell[i]+nnodes ;
    bcset += noglueBC[i] ;
    cr[fcnt] = -noglueBC[i] ;
    count[fcnt] = 4 ;
    if(noglueFaces[i][3] == -1)
      count[fcnt] = 3 ;
    fcnt++ ;
  }

  multiMap face2node ;
  face2node.allocate(count) ;
  fcnt = fstart ;
  for(size_t i=0;i<interior_face.size();++i) {
    for(int j=0;j<count[fcnt];++j)
      face2node[fcnt][j] = interior_face[i][j] ;
    fcnt++ ;
  }
  for(size_t i=0;i<glueFaces.size();++i) {
    for(int j=0;j<count[fcnt];++j)
      face2node[fcnt][j] = glueFaces[i][j] ;
    fcnt++ ;
  }
  for(size_t i=0;i<noglueFaces.size();++i) {
    for(int j=0;j<count[fcnt];++j)
      face2node[fcnt][j] = noglueFaces[i][j] ;
    fcnt++ ;
  }

  // release memory
  count.allocate(EMPTY) ;
  {vector<Loci::Array<int,4> > tmp; interior_face.swap(tmp); }
  {vector<Loci::Array<int,4> > tmp; glueFaces.swap(tmp) ; }
  {vector<Loci::Array<int,4> > tmp; noglueFaces.swap(tmp) ; }
  {vector<int> tmp; interior_cl.swap(tmp) ; }
  {vector<int> tmp; interior_cr.swap(tmp) ; }
  {vector<int> tmp; glueCell.swap(tmp) ; }
  {vector<int> tmp; glueBC.swap(tmp) ;}
  {vector<int> tmp; noglueCell.swap(tmp) ; }
  {vector<int> tmp; noglueBC.swap(tmp) ;}


  if(lefthanded) { // Left handed coordinate system
    // establish face left-right orientation
    if(MPI_rank == 0)
      cout << "orienting faces" << endl ;
    VOG::orientFaces(pos,cl,cr,face2node) ;

  }

  // Establish face orientation to be consistent with matrix coloring
  // color matrix according to the numbering of the cells
  FORALL(fdom,fc) {
    if(cl[fc] > 0 && cr[fc] > 0 && cl[fc] > cr[fc]) {
      // change face orientation to match matrix coloring
      std::swap(cl[fc],cr[fc]) ;
      int i = 0 ;
      int j = face2node[fc].size() - 1;
      while(i < j) {
        std::swap(face2node[fc][i],face2node[fc][j]) ;
        i++ ;
        j-- ;
      }
    }
  }ENDFORALL ;

  if(optimize) {
    if(MPI_rank == 0)
      cout << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }

  string outfile = case_name + ".vog" ;
  if(MPI_rank == 0)
    cout << "writing VOG file" << endl ;

  vector<pair<int,string> > surf_ids ;
  set<string> iset ;

  FORALL(bcset,bc) {
    map<int,string>::const_iterator mi ;
    if((mi = bcnamelist.find(bc)) != bcnamelist.end()) {
      surf_ids.push_back(pair<int,string>(bc,mi->second)) ;
    } else {
      char buf[512] ;
      bzero(buf,512) ;
      snprintf(buf,511,"BC_%d",bc) ;
      surf_ids.push_back(pair<int,string>(bc, string(buf))) ;
    }
  } ENDFORALL ;

  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;

  Loci::Finalize() ;
  return 0 ;
}
