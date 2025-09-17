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
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

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

bool reverse_byteorder = false ;

void check_order() {
  static int test = 15 ;
  char *p = (char *)&test ;
  if(int(*p) == test) {
    reverse_byteorder = true ;
  }
}


void ug_io_reverse_byte_order
(void * Data,
 size_t Size,
 int Number)

{

  /*
   * Set file format and host to big or little endian byte ordering.
   *
   */

  char *Data_Byte_Ptr;
  char Temp_Data_Byte;

  int Byte_Index, Index, Number_of_Bytes, Reverse_Byte_Index;

  Number_of_Bytes = int(Size);

  Data_Byte_Ptr = (char *) Data;

  for (Index = 0; Index < Number; ++Index)
    {
      Reverse_Byte_Index = Number_of_Bytes;

      for (Byte_Index = 0; Byte_Index < Number_of_Bytes/2; ++Byte_Index)
	{
	  --Reverse_Byte_Index;

	  Temp_Data_Byte = Data_Byte_Ptr[Byte_Index];

	  Data_Byte_Ptr[Byte_Index] = Data_Byte_Ptr[Reverse_Byte_Index];

	  Data_Byte_Ptr[Reverse_Byte_Index] = Temp_Data_Byte;
	}

      Data_Byte_Ptr += Number_of_Bytes;
    }

  return;
}

void input_error() {
  cerr << "error reading file" << endl ;
  exit(-1) ;
}

void cfread(void *ptr, size_t size, size_t nmemb, FILE *stream) {
  size_t nread = fread(ptr,size,nmemb,stream) ;
  if(nread != nmemb)
    input_error() ;
}

void cfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream) {
  size_t nwrite = fwrite(ptr,size,nmemb,stream) ;
  if(nwrite != nmemb) {
    cerr << "fwrite failed" << endl ;
    exit(-1) ;
  }
}
void readUGRID(string filename,bool binary, store<vector3d<double> > &pos,
               vector<Array<int,5> > &qfaces, vector<Array<int,4> > &tfaces,
               vector<Array<int,4> > &tets, vector<Array<int,5> > &pyramids,
               vector<Array<int,6> > &prisms, vector<Array<int,8> > &hexs) {

  pos.allocate(EMPTY) ;
  qfaces.clear() ;
  tfaces.clear() ;
  tets.clear() ;
  pyramids.clear() ;
  prisms.clear() ;
  hexs.clear() ;

  int num_nodes=0, num_sf_trias=0, num_sf_quads=0 ;
  int num_vol_tets=0, num_vol_pents5=0, num_vol_pents6=0, num_vol_hexs=0 ;

  FILE* IFP = NULL ;

  const int P = MPI_processes ;
  const int R = MPI_rank ;
  if(R == 0) {
    if(!binary)
      IFP = fopen(filename.c_str(), "r") ;
    else
      IFP = fopen(filename.c_str(), "rb") ;
    if(IFP == NULL) {
      cerr << "can't open '" << filename << "'" << endl ;
      exit(-1) ;
    }


    if(!binary) {
      if(fscanf(IFP, "%d%d%d", &num_nodes, & num_sf_trias, & num_sf_quads)!=3)
	input_error() ;
      if(fscanf(IFP, "%d%d%d%d", &num_vol_tets, &num_vol_pents5, &num_vol_pents6, &num_vol_hexs)!=4)
	input_error() ;
    } else {
      cfread(&num_nodes, sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&num_nodes,sizeof(int),1) ;
      cfread(&num_sf_trias, sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&num_sf_trias,sizeof(int),1) ;
      cfread(&num_sf_quads, sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&num_sf_quads,sizeof(int),1) ;
      cfread(&num_vol_tets, sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&num_vol_tets,sizeof(int),1) ;
      cfread(&num_vol_pents5, sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&num_vol_pents5,sizeof(int),1) ;
      cfread(&num_vol_pents6, sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&num_vol_pents6,sizeof(int),1) ;
      cfread(&num_vol_hexs, sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&num_vol_hexs,sizeof(int),1) ;
    }

    cout << "nnodes=" << num_nodes <<",ntria="<<num_sf_trias
         << ",nquad="<<num_sf_quads<<",ntets="<<num_vol_tets
         << ",npyrm="<<num_vol_pents5<<",nprsm="<<num_vol_pents6
         << ",nhex="<<num_vol_hexs << endl ;
  }

  Array<int,7> data ;
  data[0] = num_nodes ;
  data[1] = num_sf_trias ;
  data[2] = num_sf_quads ;
  data[3] = num_vol_tets ;
  data[4] = num_vol_pents5 ;
  data[5] = num_vol_pents6 ;
  data[6] = num_vol_hexs ;
  MPI_Bcast(&data[0],7,MPI_INT,0,MPI_COMM_WORLD) ;
  num_nodes      = data[0] ;
  num_sf_trias   = data[1] ;
  num_sf_quads   = data[2] ;
  num_vol_tets   = data[3] ;
  num_vol_pents5 = data[4] ;
  num_vol_pents6 = data[5] ;
  num_vol_hexs   = data[6] ;

  vector<int> node_ptns = VOG::simplePartitionVec(0,num_nodes-1,P) ;
  vector<entitySet> local_nodes(P) ;
  for(int i=0;i<P;++i)
    local_nodes[i] = interval(node_ptns[i],node_ptns[i+1]-1) ;

  pos.allocate(local_nodes[R]) ;

  if(R == 0) { // Read in positions

    FORALL(local_nodes[0],nd) {
      double ptmp[3] ;
      if(!binary) {
        for(int i = 0; i < 3; ++i) {
          if(fscanf(IFP, "%lf", &ptmp[i])!=1)
	    input_error() ;
        }
      } else {
        cfread(ptmp,sizeof(double),3,IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(ptmp,sizeof(double),3) ;
      }
      pos[nd] = vector3d<double>(ptmp[0],ptmp[1],ptmp[2]) ;
    } ENDFORALL ;

    int mxsize = max(local_nodes[0].size(),
                     local_nodes[P-1].size()) ;
    vector<double> buf(mxsize*3) ;

    for(int i=1;i<P;++i) {
      int cnt = 0 ;
      FORALL(local_nodes[i],nd) {
        double ptmp[3] ;
        if(!binary) {
          for(int i = 0; i < 3; ++i) {
            if(fscanf(IFP, "%lf", &ptmp[i])!=1)
	      input_error() ;
          }
        } else {
          cfread(ptmp,sizeof(double),3,IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(ptmp,sizeof(double),3) ;
        }
        buf[cnt++] = ptmp[0] ;
        buf[cnt++] = ptmp[1] ;
        buf[cnt++] = ptmp[2] ;
      } ENDFORALL ;
      MPI_Send(&buf[0],local_nodes[i].size()*3,MPI_DOUBLE,i,9,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    entitySet nodeSet = local_nodes[R] ;
    int recv_count = nodeSet.size()*3 ;
    vector<double> tmp_pos(recv_count) ;

    MPI_Recv(&tmp_pos[0],recv_count,MPI_DOUBLE,0,9,MPI_COMM_WORLD,&status) ;

    int tmp = 0 ;
    FORALL(nodeSet,nd) {
      vector3d<double> t(tmp_pos[tmp], tmp_pos[tmp+1], tmp_pos[tmp+2]) ;
      tmp += 3 ;
      pos[nd] = t ;
    } ENDFORALL ;

  }

  vector<int> triadist = VOG::simplePartitionVec(0,num_sf_trias-1,P) ;
  vector<int> quaddist = VOG::simplePartitionVec(0,num_sf_quads-1,P) ;

  qfaces = vector<Array<int,5> >(quaddist[R+1]-quaddist[R]) ;
  tfaces = vector<Array<int,4> >(triadist[R+1]-triadist[R]) ;

  // Read in boundary information
  if(R == 0) {
    { // triangles
      for(int i=triadist[R];i<triadist[R+1];++i) {
        tfaces[i][3] = 0 ;
        if(!binary) {
          if(fscanf(IFP, "%d%d%d", 
		    &tfaces[i][0], &tfaces[i][1], &tfaces[i][2])!=3)
	    input_error() ;
        } else {
          cfread(&tfaces[i][0], sizeof(int), 3, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&tfaces[i][0],sizeof(int),3) ;
        }
      }

      int tsz = max(triadist[1]-triadist[0],triadist[P]-triadist[P-1]) ;
      vector<Array<int,4> > ttmp(tsz) ;
      for(int p=1;p<P;++p) {
        int ltsz = triadist[p+1]-triadist[p] ;
        for(int i=0;i<ltsz;++i) {
          ttmp[i][3] = 0 ;
          if(!binary) {
            if(fscanf(IFP, "%d%d%d", &ttmp[i][0], &ttmp[i][1], &ttmp[i][2])!= 3)
	      input_error() ;
	  } else {
            cfread(&ttmp[i][0], sizeof(int), 3, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),3) ;
          }
        }
        if(ltsz != 0)
          MPI_Send(&ttmp[0][0],ltsz*4,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

    { // Quads
      for(int i=quaddist[R];i<quaddist[R+1];++i) {
        qfaces[i][4] = 0 ;
        if(!binary) {
          if(fscanf(IFP, "%d%d%d%d", 
		    &qfaces[i][0], &qfaces[i][1], &qfaces[i][2],
		    &qfaces[i][3])!= 4)
	    input_error() ;
        } else {
          cfread(&qfaces[i][0], sizeof(int), 4, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&qfaces[i][0],sizeof(int),4) ;
        }
      }

      int qsz = max(quaddist[1]-quaddist[0],quaddist[P]-quaddist[P-1]) ;
      vector<Array<int,5> > qtmp(qsz) ;
      for(int p=1;p<P;++p) {
        int lqsz = quaddist[p+1]-quaddist[p] ;
        for(int i=0;i<lqsz;++i) {
          qtmp[i][4] = 0 ;
          if(!binary) {
            if(fscanf(IFP, "%d%d%d%d",
		      &qtmp[i][0], &qtmp[i][1], &qtmp[i][2], &qtmp[i][3])!=4)
	      input_error() ;
          } else {
            cfread(&qtmp[i][0], sizeof(int), 4, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&qtmp[i][0],sizeof(int),4) ;
          }
        }
        if(lqsz !=0)
          MPI_Send(&qtmp[0][0],lqsz*5,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

    // Boundary tags

    { // triangles
      for(int i=triadist[R];i<triadist[R+1];++i) {
        if(!binary) {
          if(fscanf(IFP, "%d", &tfaces[i][3])!=1)
	    input_error() ;
        } else {
          cfread(&tfaces[i][3], sizeof(int), 1, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&tfaces[i][3],sizeof(int),1) ;
        }
        tfaces[i][3] = -tfaces[i][3] ; // We use negative bc tags
      }

      int tsz = max(triadist[1]-triadist[0],triadist[P]-triadist[P-1]) ;
      vector<int> ttmp(tsz) ;
      for(int p=1;p<P;++p) {
        int ltsz = triadist[p+1]-triadist[p] ;
        for(int i=0;i<ltsz;++i) {
          if(!binary) {
            if(fscanf(IFP, "%d", &ttmp[i])!=1)
	      input_error() ;
          } else {
            cfread(&ttmp[i], sizeof(int), 1, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&ttmp[i],sizeof(int),1) ;
          }
        }
        if(ltsz != 0)
          MPI_Send(&ttmp[0],ltsz,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

    { // Quads
      for(int i=quaddist[R];i<quaddist[R+1];++i) {
        if(!binary) {
          if(fscanf(IFP, "%d", &qfaces[i][4])!= 1)
	    input_error() ;
        } else {
          cfread(&qfaces[i][4], sizeof(int), 1, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&qfaces[i][4],sizeof(int),1) ;
        }
        qfaces[i][4] = -qfaces[i][4] ;
      }

      int qsz = max(quaddist[1]-quaddist[0],quaddist[P]-quaddist[P-1]) ;
      vector<int> qtmp(qsz) ;
      for(int p=1;p<P;++p) {
        int lqsz = quaddist[p+1]-quaddist[p] ;
        for(int i=0;i<lqsz;++i) {
          if(!binary) {
            if(fscanf(IFP, "%d", &qtmp[i])!=1)
	      input_error() ;
          } else {
            cfread(&qtmp[i], sizeof(int), 1, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&qtmp[i],sizeof(int),1) ;
          }
        }
        if(lqsz != 0)
          MPI_Send(&qtmp[0],lqsz,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

  } else {
    MPI_Status status ;
    int tsz = tfaces.size() ;
    if(tsz != 0)
      MPI_Recv(&tfaces[0][0],tsz*4,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;

    int qsz = qfaces.size();
    if(qsz != 0)
      MPI_Recv(&qfaces[0][0],qsz*5,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;

    vector<int> tbuf(tsz) ;
    if(tsz !=0)
      MPI_Recv(&tbuf[0],tsz,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;
    for(int i=0;i<tsz;++i)
      tfaces[i][3] = -tbuf[i] ;

    vector<int> qbuf(qsz) ;
    if(qsz != 0)
      MPI_Recv(&qbuf[0],qsz,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;
    for(int i=0;i<qsz;++i)
      qfaces[i][4] = -qbuf[i] ;
  }

  // Read in volume elements

  // Read in tetrahedra
  vector<int> tetsdist = VOG::simplePartitionVec(0,num_vol_tets-1,P) ;
  tets = vector<Array<int,4> >(tetsdist[R+1]-tetsdist[R]) ;

  if(R == 0) {
    for(int i=tetsdist[R];i<tetsdist[R+1];++i) {
      if(!binary) {
        if(fscanf(IFP, "%d%d%d%d", 
		  &tets[i][0], &tets[i][1],
		  &tets[i][2], &tets[i][3])!=4)
	  input_error() ;
      } else {
        cfread(&tets[i][0], sizeof(int), 4, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tets[i][0],sizeof(int),4) ;
      }
    }

    int tsz = max(tetsdist[1]-tetsdist[0],tetsdist[P]-tetsdist[P-1]) ;
    vector<Array<int,4> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = tetsdist[p+1]-tetsdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary) {
          if(fscanf(IFP, "%d%d%d%d", 
		    &ttmp[i][0], &ttmp[i][1], 
		    &ttmp[i][2], &ttmp[i][3])!=4)
	    input_error() ;
        } else {
          cfread(&ttmp[i][0], sizeof(int), 4, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),4) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*4,MPI_INT,p,7,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = tets.size() ;
    if(tsz != 0)
      MPI_Recv(&tets[0][0],tsz*4,MPI_INT,0,7,MPI_COMM_WORLD,&status) ;
  }

  // Read in pyramids
  vector<int> pyrmdist = VOG::simplePartitionVec(0,num_vol_pents5-1,P) ;
  pyramids = vector<Array<int,5> >(pyrmdist[R+1]-pyrmdist[R]) ;

  if(R == 0) {
    for(int i=pyrmdist[R];i<pyrmdist[R+1];++i) {
      if(!binary) {
        if(fscanf(IFP, "%d%d%d%d%d", 
		  &pyramids[i][0], &pyramids[i][1],
		  &pyramids[i][2], &pyramids[i][3], &pyramids[i][4])!=5)
	  input_error() ;
      } else {
        cfread(&pyramids[i][0], sizeof(int), 5, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&pyramids[i][0],sizeof(int),5) ;
      }
    }

    int tsz = max(pyrmdist[1]-pyrmdist[0],pyrmdist[P]-pyrmdist[P-1]) ;
    vector<Array<int,5> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = pyrmdist[p+1]-pyrmdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary) {
          if(fscanf(IFP, "%d%d%d%d%d", 
		    &ttmp[i][0], &ttmp[i][1], 
		    &ttmp[i][2], &ttmp[i][3], &ttmp[i][4])!=5)
	    input_error() ;
        } else {
          cfread(&ttmp[i][0], sizeof(int), 5, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),5) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*5,MPI_INT,p,6,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = pyramids.size() ;
    if(tsz != 0)
      MPI_Recv(&pyramids[0][0],tsz*5,MPI_INT,0,6,MPI_COMM_WORLD,&status) ;
  }

  // Read in prisms
  vector<int> prsmdist = VOG::simplePartitionVec(0,num_vol_pents6-1,P) ;
  prisms = vector<Array<int,6> >(prsmdist[R+1]-prsmdist[R]) ;

  if(R == 0) {
    for(int i=prsmdist[R];i<prsmdist[R+1];++i) {
      if(!binary) {
        if(fscanf(IFP, "%d%d%d%d%d%d", 
		  &prisms[i][0], &prisms[i][1], &prisms[i][2], 
		  &prisms[i][3], &prisms[i][4], &prisms[i][5]) != 6)
	  input_error() ;
      } else {
        cfread(&prisms[i][0], sizeof(int), 6, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&prisms[i][0],sizeof(int),6) ;
      }
    }

    int tsz = max(prsmdist[1]-prsmdist[0],prsmdist[P]-prsmdist[P-1]) ;
    vector<Array<int,6> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = prsmdist[p+1]-prsmdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary) {
          if(fscanf(IFP, "%d%d%d%d%d%d", 
		    &ttmp[i][0], &ttmp[i][1], &ttmp[i][2],
		    &ttmp[i][3], &ttmp[i][4], &ttmp[i][5])!=6)
	    input_error() ;
	} else {
          cfread(&ttmp[i][0], sizeof(int), 6, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),6) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*6,MPI_INT,p,5,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = prisms.size() ;
    if(tsz != 0)
      MPI_Recv(&prisms[0][0],tsz*6,MPI_INT,0,5,MPI_COMM_WORLD,&status) ;
  }

  // Read in hexahdra
  vector<int> hexsdist = VOG::simplePartitionVec(0,num_vol_hexs-1,P) ;
  hexs = vector<Array<int,8> >(hexsdist[R+1]-hexsdist[R]) ;

  if(R == 0) {
    for(int i=hexsdist[R];i<hexsdist[R+1];++i) {
      if(!binary) {
        if(fscanf(IFP, "%d%d%d%d%d%d%d%d", 
		  &hexs[i][0], &hexs[i][1], &hexs[i][2], &hexs[i][3],
		  &hexs[i][4], &hexs[i][5], &hexs[i][6], &hexs[i][7])!=8)
	  input_error() ;
      } else {
        cfread(&hexs[i][0], sizeof(int), 8, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&hexs[i][0],sizeof(int),8) ;
      }
    }

    int tsz = max(hexsdist[1]-hexsdist[0],hexsdist[P]-hexsdist[P-1]) ;
    vector<Array<int,8> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = hexsdist[p+1]-hexsdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary) {
          if(fscanf(IFP, "%d%d%d%d%d%d%d%d", 
		    &ttmp[i][0], &ttmp[i][1], &ttmp[i][2], &ttmp[i][3], 
		    &ttmp[i][4], &ttmp[i][5], &ttmp[i][6], &ttmp[i][7])!=8)
	    input_error() ;
        } else {
          cfread(&ttmp[i][0], sizeof(int), 8, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),8) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*8,MPI_INT,p,4,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = hexs.size() ;
    if(tsz != 0)
      MPI_Recv(&hexs[0][0],tsz*8,MPI_INT,0,4,MPI_COMM_WORLD,&status) ;
  }
}

void writeUGRID(string filename, vector<vector3d<double> > &pos,
               vector<Array<int,4> > &tfaces,
		vector<Array<int,4> > &tets) {

  int num_nodes=pos.size(), num_sf_trias=tfaces.size(), num_sf_quads=0 ;
  int num_vol_tets=tets.size(), num_vol_pents5=0, num_vol_pents6=0, num_vol_hexs=0 ;

  FILE* IFP = NULL ;

  IFP = fopen(filename.c_str(), "wb") ;
  if(IFP == NULL) {
    cerr << "can't open '" << filename << "'" << endl ;
    exit(-1) ;
  }
  
  int tmp = num_nodes ;
  if(reverse_byteorder)
    ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
  cfwrite(&tmp, sizeof(int), 1, IFP) ;
  
  tmp = num_sf_trias ;
  if(reverse_byteorder)
    ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
  cfwrite(&tmp, sizeof(int), 1, IFP) ;
  
  tmp = num_sf_quads ;
  if(reverse_byteorder)
    ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
  cfwrite(&tmp, sizeof(int), 1, IFP) ;
 
  tmp = num_vol_tets ;
  if(reverse_byteorder)
    ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
  cfwrite(&tmp, sizeof(int), 1, IFP) ;
 
  tmp = num_vol_pents5 ;
  if(reverse_byteorder)
    ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
  cfwrite(&tmp, sizeof(int), 1, IFP) ;

  tmp = num_vol_pents6 ;
  if(reverse_byteorder)
    ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
  cfwrite(&tmp, sizeof(int), 1, IFP) ;

  tmp = num_vol_hexs ;
  if(reverse_byteorder)
    ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
  cfwrite(&tmp, sizeof(int), 1, IFP) ;

  cout << "writing file" << endl ;
  cout << "nnodes=" << num_nodes <<",ntria="<<num_sf_trias
         << ",nquad="<<num_sf_quads<<",ntets="<<num_vol_tets
         << ",npyrm="<<num_vol_pents5<<",nprsm="<<num_vol_pents6
         << ",nhex="<<num_vol_hexs << endl ;

  for(int i=0;i<num_nodes;++i) {
    double ptmp[3] ;
    ptmp[0] = pos[i].x ;
    ptmp[1] = pos[i].y ;
    ptmp[2] = pos[i].z ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(ptmp,sizeof(double),3) ;
    cfwrite(ptmp,sizeof(double),3,IFP) ;
  }

  for(int i = 0;i<num_sf_trias;++i) {
    for(int j=0;j<3;++j) {
      tmp = tfaces[i][j] ;
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
      cfwrite(&tmp, sizeof(int), 1, IFP) ;
    }
  }


    // Boundary tags

  for(int i = 0;i<num_sf_trias;++i) {
    tmp = -tfaces[i][3] ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
    cfwrite(&tmp, sizeof(int), 1, IFP) ;
  }

  for(int i=0;i<num_vol_tets;++i) {
    for(int j=0;j<4;++j) {
      tmp = tets[i][j] ;
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&tmp,sizeof(int),1) ;
      cfwrite(&tmp, sizeof(int), 1, IFP) ;
    }
  }
  fclose(IFP) ;
}


int main(int ac, char* av[]) {
  using namespace VOG ;

  Loci::Init(&ac,&av) ;
  const char *filename ;
  std::string tmp_str ;
  bool binary = 0;
  bool mirrorx = true ;
  bool mirrory = false ;
  bool mirrorz = false ;
  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 2 && !strcmp(av[1],"-v")) {
      cout << "Loci version: " << Loci::version() << endl ;
      if(ac == 2) {
        Loci::Finalize() ;
        exit(0) ;
      }
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-b")) {
      binary = 1 ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-mirrorx")) {
	mirrorx=true ;
	mirrory=false ;
	mirrorz=false ;
	ac-- ;
	av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-mirrory")) {
	mirrorx=false ;
	mirrory=true ;
	mirrorz=false ;
	ac-- ;
	av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-mirrorz")) {
	mirrorx=false ;
	mirrory=false ;
	mirrorz=true ;
	ac-- ;
	av++ ;
    } else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      ac-- ;
      av++ ;
    }
  }
      
  if(ac == 2) {
    tmp_str.append(av[1]) ;
  } else {
    cerr << "Usage: ugridmirror <options> <file>" << endl
         << "Where options are listed below and <file> is the filename sans postfix" << endl
         << "flags:" << endl
         << "  -b  : assume input file is binary" << endl
	 << "  -mirrorx : mirror in x coordinate" << endl 
	 << "  -mirrory : mirror in y coordinate" << endl 
	 << "  -mirrorz : mirror in z coordinate" << endl 
         << endl ;
    exit(-1) ;
  }
      

  // Check machines internal byte order
  check_order() ;

  int loc = 0;
  loc = tmp_str.find('.') ;
  std::string new_str = tmp_str.substr(0, loc) ;
  filename = new_str.c_str() ;
  char buf[512] ;
  bzero(buf,512) ;
  if(!binary) {
    struct stat fstat ;
    snprintf(buf,511,"%s.ugrid",filename) ;
    if(stat(buf,&fstat)<0) {
      binary = true ;
    }
  }

  if(!binary)
    snprintf(buf,511,"%s.ugrid",filename) ;
  else
    snprintf(buf,511,"%s.b8.ugrid",filename) ;
  string infile = buf ;

  string outfile = string(filename) ;
  if(mirrorx)
    outfile += "x" ;
  if(mirrory)
    outfile += "y" ;
  if(mirrorz)
    outfile += "z" ;
  outfile += string(".b8.ugrid") ;

  store<vector3d<double> > pos ;
  vector<Array<int,5> > qfaces ;
  vector<Array<int,4> > tfaces ;
  vector<Array<int,4> > tets ;
  vector<Array<int,5> > pyramids ;
  vector<Array<int,6> > prisms ;
  vector<Array<int,8> > hexs ;

  readUGRID(infile, binary, pos,qfaces,tfaces,tets,pyramids,prisms,hexs) ;

  if(qfaces.size()!=0 || pyramids.size()!=0 || prisms.size()!=0 || hexs.size() != 0) {
    cerr << "only supprting all tet meshes" << endl ;
  }
  
  int npos = pos.domain().size() ;
  vector<vector3d<double> > posm(npos);
  vector<int> nposid(npos*2) ;

  for(int i=0;i<npos;++i) {
    posm[i]=pos[i] ;
    nposid[i] = i+1 ;
  }
  int cnt = npos ;
  for(int i=0;i<npos;++i)  {
    vector3d<double> nvpos = posm[i] ;
    nvpos.x = mirrorx?(-nvpos.x):(nvpos.x) ;
    nvpos.y = mirrory?(-nvpos.y):(nvpos.y) ;
    nvpos.z = mirrorz?(-nvpos.z):(nvpos.z) ;
    double dist = norm(nvpos-posm[i]) ;
    if(dist > 1e-7) { // unique point
      nposid[npos+i] = cnt+1 ;
      posm.push_back(nvpos) ;
      cnt++ ;
    } else { // cloned point is removed
      nposid[npos+i] = i+1 ;
    }
  }
  cout << "glued nodes = " <<  npos-(cnt-npos) << endl ;
  // Now glue faces
  vector<Array<int,4> > ntfaces ;
  for(size_t i=0;i<tfaces.size();++i) {
    Array<int,4> face1 = tfaces[i] ;
    Array<int,4> face2 = face1 ;
    for(int j=0;j<3;++j) 
      face2[j] = nposid[face1[j]+npos-1] ;
    if((face1[0] != face2[0]) || 
       (face1[1] != face2[1]) || 
       (face1[2] != face2[2])) {
      // unique face, but before adding lets check to see the normal aligns
      // with the mirror direction
      vector3d<double> dv = posm[face2[0]-1]-posm[face1[0]-1] ;
      double dist = (norm(posm[face2[1]-1]-posm[face1[1]-1]-dv) +
		     norm(posm[face2[2]-1]-posm[face1[2]-1]-dv)) ;
      if(dist < 1e-8) { // parallel face, so lets add offset to 
	// facid so it is a different face
	face2[3] += -16 ;
      }
      std::swap(face2[0],face2[2]) ;
      ntfaces.push_back(face1) ;
      ntfaces.push_back(face2) ;
    }
  }
  int num_tets = tets.size() ;
  vector<Array<int,4> > ntets(num_tets*2) ;
  for(int i=0;i<num_tets;++i) {
    ntets[i] = tets[i] ;
    for(int j=0;j<4;++j) {
      ntets[i+num_tets][j] = nposid[tets[i][j]+npos-1] ;
    }
    std::swap(ntets[i+num_tets][0],ntets[i+num_tets][2]) ; // mirror image swap order
  }
  // Now write out new mesh

  writeUGRID(outfile,posm,ntfaces,ntets) ;
  

  Loci::Finalize() ;
}
