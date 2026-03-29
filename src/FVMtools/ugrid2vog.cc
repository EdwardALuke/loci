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
bool split_file_exist = false;


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
inline bool pairEqual(const pair<int,int> &p1, const pair<int,int> &p2) {
  return ((p1.first == p2.first) &&
          (p1.second == p2.second)) ;
}
//as shown  in simcenter/system/release/doc/ug_io/3d_input_output_grids.html SPLIT FACE HEX ELEMENT CONNECTIVITY
struct quadSplit {
  Array<int,5> nodes ;
  Array<int,4> corners ;
  int cell ;
  int face;
} ;

struct edgeSplit {
  int n1;
  int n2;
  int mp; //middle point

  edgeSplit(int nn1, int nn2, int mpp) {
    n1=min(nn1,nn2) ;
    n2=max(nn1,nn2) ;
    mp = mpp;
  }
  edgeSplit() {}
} ;

inline bool edgeSplitEqual(const edgeSplit &e1, const edgeSplit &e2) {
  return ((e1.n1 == e2.n1) &&
          (e1.n2 == e2.n2)) ;

}

inline bool edgeSplitCompare(const edgeSplit &e1, const edgeSplit &e2) {
  return ((e1.n1 < e2.n1) ||
          (e1.n1 == e2.n1 &&
	   e1.n2 < e2.n2)) ;
}

inline bool quadSplitCompare(const pair<int, int>& e1, const pair<int, int> &e2) {
  return ((e1.first < e2.first) ||
          (e1.first == e2.first &&
	   e1.second < e2.second)) ;
}


//collect all edges that are split, each process has all the edges that are split
void collect_edge_split( vector<edgeSplit>& es, vector<quadSplit>& qs){
  es.clear();

  size_t num_face = qs.size();
  if(num_face > 0) es.reserve(num_face*4);
  for(size_t i = 0; i < num_face; i++){
    for(int j = 0; j < 4; j++){
      if(qs[i].nodes[j]>0)es.push_back(edgeSplit((qs[i].corners[j])-1,(qs[i].corners[(j+1)%4])-1, (qs[i].nodes[j])-1));
    }
  }
  std::sort(es.begin(), es.end(), edgeSplitCompare) ;
  vector<edgeSplit>::iterator itr =  std::unique(es.begin(), es.end(), edgeSplitEqual);
  es.resize(itr-es.begin());

  const int P = Loci::MPI_processes;
  //int my_size = es.size()*3;
  // if(P>1){
  //   vector<int> sizes(P);
  //   MPI_Allgather(&my_size,1, MPI_INT, &sizes[0], 1, MPI_INT, MPI_COMM_WORLD) ;
  //   int total_size = 0;
  //   for(int i = 0; i < P; i++)total_size += sizes[i];
  //   vector<edgeSplit> recv(total_size/3);
  //   vector<int> recv_disp(P);
  //   recv_disp[0] = 0 ;
  //   for(int i=1;i<P;++i)
  //     recv_disp[i] = recv_disp[i-1]+sizes[i-1];
  //   MPI_Allgatherv(&es[0], my_size, MPI_INT, &recv[0], &sizes[0], &recv_disp[0], MPI_INT, MPI_COMM_WORLD);
  //   std::sort(recv.begin(), recv.end(), edgeSplitCompare) ;
  //   vector<edgeSplit>::iterator itr =  std::unique(recv.begin(), recv.end(), edgeSplitEqual);
  //   recv.resize(itr-recv.begin());
  //   es.swap(recv);
  // }

  if(P > 1){
    Loci::parSampleSort(es,edgeSplitCompare,MPI_COMM_WORLD) ;
    vector<edgeSplit>::iterator itr =  std::unique(es.begin(), es.end(), edgeSplitEqual);
    es.resize(itr-es.begin());
  }
}



//collect the edges each processor owns, then sort and unique them locally
void collect_edge_from_face(const vector<quadFace>& quad,
                            const vector<triaFace>& tria,
                            vector<pair<int, int> >& es){
  es.clear();

  size_t num_quad = quad.size();
  size_t num_tria = tria.size();

  try{
    es.reserve(num_quad*4+num_tria*3);
  }catch(const std::exception& e) { // caught by reference to base
    std::cerr << " a standard exception was caught, with message '"
              << e.what() << endl;
  }

  for(size_t i = 0; i < num_quad; i++){
    for(size_t j = 0; j < 4; j++){
      int n1 = quad[i].nodes[j];
      int n2 = quad[i].nodes[(j+1)%4];
      es.push_back(pair<int, int>(min(n1, n2),max(n1, n2)));
    }
  }
  for(size_t i = 0; i < num_tria; i++){
    for(int j = 0; j < 3; j++){
      int n1 = tria[i].nodes[j];
      int n2 = tria[i].nodes[(j+1)%3];
      es.push_back(pair<int, int>(min(n1, n2), max(n1, n2)));
    }
  }
  std::sort(es.begin(), es.end(), pairCompare) ;
  vector<pair<int, int> >::iterator itr =  std::unique(es.begin(), es.end(), pairEqual);
  es.resize(itr-es.begin());
}
//using binary search to check if an edge is split
int find_edge_split( vector<edgeSplit>& edge_splits, int node1, int node2){
  if(edge_splits.empty())return -1;
  int first = 0;
  int last = edge_splits.size()-1;
  int middle = (first+last)/2;
  edgeSplit e(node1, node2, 0);
  while(first <=last){
    if(edgeSplitEqual(edge_splits[middle], e)) return edge_splits[middle].mp;
    else if(edgeSplitCompare(edge_splits[middle], e)) first = middle +1;
    else last = middle -1;
    middle = (first+last)/2;
  }
  return -1;
}
//redistribute edge_splits according to the edges each processor owns
//es: all edges each processor owns, sorted
//edge_splits:
//     before: sorted, distributed according to global order
//     after: sorted, distributed according to the edges each processor owns
void redistribute_edge_split( const vector<pair<int, int> >& es,  vector<edgeSplit>& edge_splits){
  int P = MPI_processes;
  size_t num_edges = es.size();
  size_t num_splits = edge_splits.size();
  if(MPI_rank == 0) cerr << "      gather  the first splits of each processor " << endl;

  //first all gather the first splits of each processor
  int dummy_number = std::numeric_limits<int>::max();
  vector<int> first_split(2, dummy_number);
  if(num_splits> 0){
    first_split[0] = (edge_splits[0]).n1;
    first_split[1] = (edge_splits[0]).n2;
  }
  vector<pair<int, int> > first_splits(P+1);
  MPI_Allgather(&first_split[0],2, MPI_INT, &first_splits[0], 2, MPI_INT, MPI_COMM_WORLD) ;
  first_splits[P] = pair<int, int>(dummy_number, dummy_number);

  if(MPI_rank == 0) cerr << "     for each edge in es, find which processor owns its splits   " << endl;
  vector<int> owner_p(num_edges, -1);
  vector<int> send_count(P, 0);
  vector<int> recv_count(P, 0);
  int current_p = 0;

  for(size_t i = 0; i < num_edges; i++){
    if(pairCompare(es[i], first_splits[current_p+1])){
      send_count[current_p] += 2;
      owner_p[i] = current_p;
    }else{
      while(current_p < P && (!pairCompare(es[i], first_splits[current_p+1]))) current_p++;
      if(pairCompare(es[i], first_splits[current_p+1])){
        send_count[current_p] += 2;
        owner_p[i] = current_p;
      }
    }
  }
  //send the edges to the processor that owns the splits
  MPI_Alltoall(&send_count[0],1,MPI_INT, &recv_count[0], 1, MPI_INT,MPI_COMM_WORLD) ;

  vector<int> send_displacement(P) ;
  vector<int> recv_displacement(P);
  send_displacement[0] = 0 ;
  recv_displacement[0] = 0 ;
  for(int i=1;i<P;++i) {
    send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
    recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
  }

  int mp = P-1 ;
  size_t send_sizes = send_displacement[mp]+send_count[mp] ;
  size_t recv_sizes = recv_displacement[mp]+recv_count[mp] ;
  if(send_sizes != num_edges*2){
    cerr<<" ERROR: Not all edges sent in redistribute_edge_split()" << endl;
    Loci::Abort();
  }

  int *send_set_buf = new int[send_sizes] ;
  int *recv_set_buf = new int[recv_sizes] ;
  size_t ind = 0;

  try{
    for(int p = 0; p < P; p++){
      for(int j = 0; j < send_count[p]/2; j++){
        send_set_buf[send_displacement[p]+2*j] = (es.at(ind)).first ;
        send_set_buf[send_displacement[p]+2*j+1] = (es.at(ind)).second ;
        ind++;
      }
    }
  }catch (const std::exception& e) { // caught by reference to base
    std::cerr << " a standard exception was caught, with message '"
              << e.what() << endl;
    exit(1);
  }

  if(ind != num_edges){
    cerr<<"ERROR: not all edges are packed" << endl;
    Loci::Abort();
  }

  MPI_Alltoallv(send_set_buf, &send_count[0], &send_displacement[0] , MPI_INT,
                recv_set_buf, &recv_count[0], &recv_displacement[0], MPI_INT,
                MPI_COMM_WORLD) ;

  //allocate memory for middle points
  vector<int> send_count2(P, 0);
  vector<int> recv_count2(P, 0);
  for(int i = 0; i < P; i++){
    send_count2[i] = recv_count[i]/2;
    recv_count2[i] = send_count[i]/2;
  }
  int send_sizes2 = recv_sizes/2 ;
  int recv_sizes2 = send_sizes/2 ;
  int* send_buf2 = new int[send_sizes2] ;
  int* send_displacement2 = new int[P] ;
  int* recv_displacement2 = new int[P] ;

  send_displacement2[0] = 0 ;
  recv_displacement2[0] = 0 ;
  for(int i=1;i<P;++i) {
    send_displacement2[i] = send_displacement2[i-1]+send_count2[i-1] ;
    recv_displacement2[i] = recv_displacement2[i-1]+recv_count2[i-1] ;
  }

  //find the middle point of each received edge
  for(int i=0;i<P;++i) {
    for(int j=0;j<recv_count[i]/2;++j) {
      int i1 = recv_set_buf[recv_displacement[i]+j*2  ] ;
      int i2 = recv_set_buf[recv_displacement[i]+j*2+1] ;
      int middle = find_edge_split(edge_splits, i1, i2);
      send_buf2[send_displacement2[i]+j] = middle;
    }
  }

  delete[] recv_set_buf ;
  delete[] send_set_buf ;

  int* recv_buf2 = new int[recv_sizes2] ;
  //send the middle points back
  MPI_Alltoallv(send_buf2, &send_count2[0], send_displacement2 , MPI_INT,
                recv_buf2, &recv_count2[0], recv_displacement2, MPI_INT,
                MPI_COMM_WORLD) ;
  delete[] send_buf2 ;
  delete[] send_displacement2;
  delete[] recv_displacement2;


  //get resulted edge_splits
  vector<edgeSplit> result_edge_splits;
  result_edge_splits.reserve(num_edges);
  for(size_t i = 0 ; i < num_edges; i++){
    if(recv_buf2[i] >= 0){
      result_edge_splits.push_back(edgeSplit(es[i].first, es[i].second, recv_buf2[i]));
    }
  }
  edge_splits.swap(result_edge_splits);
  delete[] recv_buf2 ;

}







//parallel read split file
void readSplit(string filename, const vector<int>& hex_min, vector<quadSplit>& splits){

  FILE* IFP = NULL ;
  const int P = MPI_processes ;
  const int R = MPI_rank ;
  int num_splits = 0;
  if(R == 0) {
    IFP = fopen(filename.c_str(), "r") ;
    if(IFP == NULL) {
      cerr << "can't open '" << filename << "'" << endl ;
      exit(-1) ;
    }

    if(fscanf(IFP, "%d", &num_splits)!=1)input_error() ;
    cout<< "num_splits " << num_splits << endl;
    quadSplit quad_face;
    int current_processor = 0;
    vector<quadSplit> buf;
    for(int i = 0; i < num_splits; i++){
      if(fscanf(IFP, "%d%d%d%d%d%d%d%d%d%d%d", &(quad_face.cell), &(quad_face.face), &(quad_face.nodes[0]),
                &(quad_face.nodes[1]), &(quad_face.nodes[2]), &(quad_face.nodes[3]), &(quad_face.nodes[4]),
                &(quad_face.corners[0]), &(quad_face.corners[1]), &(quad_face.corners[2]), &(quad_face.corners[3]) ) != 11)input_error() ;

      if(quad_face.cell < hex_min[current_processor+1]){//belongs to current_processor
        buf.push_back(quad_face);//store it
      }else{//not belongs to current_processor,
        //first finish with current processor, either copy or send the data in buf, and then clear the buf
        if(current_processor == 0) { //if currrent_processor is root
          splits = vector<quadSplit>(buf); //copy
          buf.clear(); // and clear buf

        }else{ //current_processor is not root
          int size  = buf.size();
          MPI_Send(&size,1,MPI_INT,current_processor,1,MPI_COMM_WORLD) ; //send the size
          if(size > 0)MPI_Send(&buf[0], size*11,MPI_INT,current_processor,2,MPI_COMM_WORLD) ; // and buf
          buf.clear();//then clean buf
        }

        //Next which processor it belongs
        while(current_processor < P){
          current_processor++; //next?
          if(quad_face.cell >= hex_min[current_processor+1]){ //not next
            int size = 0;
            MPI_Send(&size,1,MPI_INT,current_processor,1,MPI_COMM_WORLD) ; //send the size
          }else{ //yes, found the processor
            break;
          }
        }
        if(current_processor == P){
          cerr<<" Proc " << R << " : can not fount processor for hex number " << quad_face.cell << endl;
          exit(-1);
        }
        buf.push_back(quad_face);
      }//end the case not belongs to current_processor


      if(i == (num_splits-1)){//finsh reading, send buf
        if(current_processor == 0) { //if currrent_processor is root
          splits = vector<quadSplit>(buf); //copy
          buf.clear(); // and clear buf
          //tell others their size is 0
          int size = 0;
          for(int other_processor = 1; other_processor < P; other_processor++){
            MPI_Send(&size,1,MPI_INT,other_processor,1,MPI_COMM_WORLD) ; //send the size
          }
        }else{ //current_processor is not root
          int size  = buf.size();
          MPI_Send(&size,1,MPI_INT,current_processor,1,MPI_COMM_WORLD) ; //send the size
          if(size > 0)MPI_Send(&buf[0], size*11,MPI_INT,current_processor,2,MPI_COMM_WORLD) ; // and buf
          buf.clear();//then clean buf
          //tell others their size is 0
          if((current_processor+1) < P){
            int size = 0;
            for(int other_processor = current_processor+1; other_processor < P; other_processor++){
              MPI_Send(&size,1,MPI_INT,other_processor,1,MPI_COMM_WORLD) ; //send the size
            }

          }
        }
      }
    }
    fclose(IFP);

  }else{//non-root processors
    int size = 0;
    int recv_count = 1;
    MPI_Status status ;
    MPI_Recv(&size,recv_count,MPI_INT,0,1,MPI_COMM_WORLD,&status) ;
    if(size > 0){
      recv_count = size*11;
      splits.resize(size);
      MPI_Recv(&splits[0],recv_count,MPI_INT,0,2,MPI_COMM_WORLD,&status) ;
    }
  }
  //sanity check
  int total_num_splits = 0;
  int my_num_splits = splits.size();
  MPI_Reduce(&my_num_splits,&total_num_splits,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if(R==0 && total_num_splits != num_splits){
    cerr<<"ERROR: not all splits are read in" << endl;
    exit(-1);
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
      cerr << "c1=" << c1 << "c2=" << c2 << endl ;
      cerr << "nodes=" ;
      for(int i=0;i<3;++i)
        cerr << " " << tria[i*2].nodes[i] ;
      cerr << endl ;
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
//get the 4 nodes of face face_id in a hex, the face normal always points outward
void get_corners(int face_id, const Array<int,8> &hex, Array<int, 4> &corners){
  switch(face_id){
  case 3:
    corners[0] = hex[0] ;
    corners[1] = hex[1] ;
    corners[2] = hex[5] ;
    corners[3] = hex[4] ;
    break;
  case 2:

    corners[0] = hex[1] ;
    corners[1] = hex[2] ;
    corners[2] = hex[6] ;
    corners[3] = hex[5] ;
    break;

  case 4:
    corners[0] = hex[2] ;
    corners[1] = hex[3] ;
    corners[2] = hex[7] ;
    corners[3] = hex[6] ;
    break;


  case 1:
    corners[0] = hex[4] ;
    corners[1] = hex[7] ;
    corners[2] = hex[3] ;
    corners[3] = hex[0] ;
    break;

  case 6:
    corners[0] = hex[4] ;
    corners[1] = hex[5] ;
    corners[2] = hex[6] ;
    corners[3] = hex[7] ;
    break;
  case 5:
    corners[0] = hex[3] ;
    corners[1] = hex[2] ;
    corners[2] = hex[1] ;
    corners[3] = hex[0] ;
    break;
  default:
    cerr<<"ERROR: illegal face_id value in get_corners(): " << face_id << endl;
    exit(1);
  }
}
//using binary search to check if a face is split
int find_quad_split( vector<quadSplit>& splits,  int hex_id, int face_id, int start){
  if(splits.empty())return -1;
  int first = start;
  int last = splits.size()-1;
  int middle = (first+last)/2;
  while(first <=last){
    if(splits[middle].cell==hex_id && splits[middle].face==face_id) return middle;
    else if(quadSplitCompare(pair<int, int>(splits[middle].cell,splits[middle].face), pair<int, int>(hex_id, face_id))) first = middle +1;
    else last = middle -1;
    middle = (first+last)/2;
  }
  return -1;
}

bool split(const quadSplit& s,
           vector<quadFace> &quad, int& qf,
           vector<triaFace> &tria, int& tf,
           int cellid){
  // 5 mid points, 4 quads
  if(s.nodes[0] > 0 && s.nodes[1] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 && s.nodes[4] > 0){
    //lower right
    quad[qf].nodes[0] = s.nodes[0] ;
    quad[qf].nodes[1] = s.corners[1] ;
    quad[qf].nodes[2] = s.nodes[1] ;
    quad[qf].nodes[3] = s.nodes[4] ;
    quad[qf].left = false ;
    quad[qf++].cell = cellid ;

    //upper right
    quad[qf].nodes[0] = s.nodes[4] ;
    quad[qf].nodes[1] = s.nodes[1] ;
    quad[qf].nodes[2] = s.corners[2] ;
    quad[qf].nodes[3] = s.nodes[2] ;
    quad[qf].left = false ;
    quad[qf++].cell = cellid ;

    //upper left
    quad[qf].nodes[0] = s.nodes[3] ;
    quad[qf].nodes[1] = s.nodes[4] ;
    quad[qf].nodes[2] = s.nodes[2] ;
    quad[qf].nodes[3] = s.corners[3];
    quad[qf].left = false ;
    quad[qf++].cell = cellid ;

    //lower left
    quad[qf].nodes[0] = s.corners[0] ;
    quad[qf].nodes[1] = s.nodes[0] ;
    quad[qf].nodes[2] = s.nodes[4] ;
    quad[qf].nodes[3] = s.nodes[3];
    quad[qf].left = false ;
    quad[qf++].cell = cellid ;
    return true;

  }
  if(s.nodes[4] == 0){
    // 4 mid points, 1 quad, 4 trias
    if(s.nodes[0] > 0  && s.nodes[1] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      quad[qf].nodes[0] = s.nodes[0] ;
      quad[qf].nodes[1] = s.nodes[1] ;
      quad[qf].nodes[2] = s.nodes[2] ;
      quad[qf].nodes[3] = s.nodes[3];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    // 3 mid points, 1 quad, 3 trias
    if(s.nodes[0] == 0  && s.nodes[1] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      quad[qf].nodes[0] = s.corners[0] ;
      quad[qf].nodes[1] = s.corners[1] ;
      quad[qf].nodes[2] = s.nodes[1] ;
      quad[qf].nodes[3] = s.nodes[3];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[1] ;
      tria[tf].nodes[1] = s.nodes[2] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[0] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      quad[qf].nodes[0] = s.corners[1] ;
      quad[qf].nodes[1] = s.corners[2] ;
      quad[qf].nodes[2] = s.nodes[2] ;
      quad[qf].nodes[3] = s.nodes[0];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[0] ;
      tria[tf].nodes[1] = s.nodes[2] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[2] == 0  && s.nodes[0] > 0 && s.nodes[1] > 0 && s.nodes[3] > 0 ){
      quad[qf].nodes[0] = s.corners[2] ;
      quad[qf].nodes[1] = s.corners[3] ;
      quad[qf].nodes[2] = s.nodes[3] ;
      quad[qf].nodes[3] = s.nodes[1];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[1] ;
      tria[tf].nodes[1] = s.nodes[3] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[3] == 0  && s.nodes[0] > 0 && s.nodes[1] > 0 && s.nodes[2] > 0 ){
      quad[qf].nodes[0] = s.corners[3] ;
      quad[qf].nodes[1] = s.corners[0] ;
      quad[qf].nodes[2] = s.nodes[0] ;
      quad[qf].nodes[3] = s.nodes[2];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[0] ;
      tria[tf].nodes[1] = s.nodes[1] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.nodes[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    //2 mid points, 4 trias
    if(s.nodes[0] == 0  && s.nodes[1] == 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      tria[tf].nodes[0] = s.nodes[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[1] ;
      tria[tf].nodes[1] = s.nodes[2] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[1] ;
      tria[tf].nodes[1] = s.nodes[3] ;
      tria[tf].nodes[2] = s.corners[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[2] == 0 && s.nodes[3] > 0 && s.nodes[0] > 0 ){
      tria[tf].nodes[0] = s.nodes[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[2] ;
      tria[tf].nodes[1] = s.nodes[3] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[2] ;
      tria[tf].nodes[1] = s.nodes[0] ;
      tria[tf].nodes[2] = s.corners[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[2] == 0  && s.nodes[3] == 0 && s.nodes[0] > 0 && s.nodes[1] > 0 ){
      tria[tf].nodes[0] = s.nodes[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[3] ;
      tria[tf].nodes[1] = s.nodes[0] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[3] ;
      tria[tf].nodes[1] = s.nodes[1] ;
      tria[tf].nodes[2] = s.corners[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[3] == 0  && s.nodes[0] == 0 && s.nodes[1] > 0 && s.nodes[2] > 0 ){
      tria[tf].nodes[0] = s.nodes[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[0] ;
      tria[tf].nodes[1] = s.nodes[1] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[0] ;
      tria[tf].nodes[1] = s.nodes[2] ;
      tria[tf].nodes[2] = s.corners[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;

      tria[tf].nodes[0] = s.corners[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    //2 mid points, 2 quads
    if(s.nodes[0] == 0  && s.nodes[2] == 0 && s.nodes[1] > 0 && s.nodes[3] > 0 ){

      quad[qf].nodes[0] = s.corners[0] ;
      quad[qf].nodes[1] = s.corners[1] ;
      quad[qf].nodes[2] = s.nodes[1] ;
      quad[qf].nodes[3] = s.nodes[3];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;

      quad[qf].nodes[0] = s.corners[2] ;
      quad[qf].nodes[1] = s.corners[3] ;
      quad[qf].nodes[2] = s.nodes[3] ;
      quad[qf].nodes[3] = s.nodes[1];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[3] == 0 && s.nodes[0] > 0 && s.nodes[2] > 0 ){
      quad[qf].nodes[0] = s.corners[1] ;
      quad[qf].nodes[1] = s.corners[2] ;
      quad[qf].nodes[2] = s.nodes[2] ;
      quad[qf].nodes[3] = s.nodes[0];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;

      quad[qf].nodes[0] = s.corners[3] ;
      quad[qf].nodes[1] = s.corners[0] ;
      quad[qf].nodes[2] = s.nodes[0] ;
      quad[qf].nodes[3] = s.nodes[2];
      quad[qf].left = false ;
      quad[qf++].cell = cellid ;
      return true;
    }

    //1 mid point, 3 trias
    if(s.nodes[0] == 0  && s.nodes[1] == 0 && s.nodes[2] == 0 && s.nodes[3] > 0 ){
      tria[tf].nodes[0] = s.corners[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[3] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[2] == 0 && s.nodes[3] == 0 && s.nodes[0] > 0 ){
      tria[tf].nodes[0] = s.corners[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[0] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[2] == 0  && s.nodes[3] == 0 && s.nodes[0] == 0 && s.nodes[1] > 0 ){
      tria[tf].nodes[0] = s.corners[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[2] ;
      tria[tf].nodes[1] = s.corners[3] ;
      tria[tf].nodes[2] = s.nodes[1] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }

    if(s.nodes[3] == 0  && s.nodes[0] == 0 && s.nodes[1] == 0 && s.nodes[2] > 0 ){
      tria[tf].nodes[0] = s.corners[0] ;
      tria[tf].nodes[1] = s.corners[1] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[1] ;
      tria[tf].nodes[1] = s.corners[2] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;


      tria[tf].nodes[0] = s.corners[3] ;
      tria[tf].nodes[1] = s.corners[0] ;
      tria[tf].nodes[2] = s.nodes[2] ;
      tria[tf].left = false ;
      tria[tf++].cell = cellid ;
      return true;
    }
  }
  cerr<<" split() has cases that can not handle" << s.nodes[0] << " "<< s.nodes[1] << " "<<s.nodes[2] << " "<<s.nodes[3] << " "<<s.nodes[4] <<endl;
  return false;
}


bool count_split(const quadSplit& s,
                 int& qf,
                 int& tf){
  // 5 mid points, 4 quads
  if(s.nodes[0] > 0 && s.nodes[1] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 && s.nodes[4] > 0){
    qf += 4;
    return true;
  }
  if(s.nodes[4] == 0){
    // 4 mid points, 1 quad, 4 trias
    if(s.nodes[0] > 0  && s.nodes[1] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      qf += 1;
      tf += 4;
      return true;
    }

    // 3 mid points, 1 quad, 3 trias
    if(s.nodes[0] == 0  && s.nodes[1] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      qf += 1;
      tf += 3;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[0] > 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      qf += 1;
      tf += 3;
      return true;
    }

    if(s.nodes[2] == 0  && s.nodes[0] > 0 && s.nodes[1] > 0 && s.nodes[3] > 0 ){
      qf += 1;
      tf += 3;
      return true;
    }

    if(s.nodes[3] == 0  && s.nodes[0] > 0 && s.nodes[1] > 0 && s.nodes[2] > 0 ){
      qf += 1;
      tf += 3;
      return true;
    }

    //2 mid points, 4 trias
    if(s.nodes[0] == 0  && s.nodes[1] == 0 && s.nodes[2] > 0 && s.nodes[3] > 0 ){
      tf += 4;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[2] == 0 && s.nodes[3] > 0 && s.nodes[0] > 0 ){
      tf += 4;
      return true;
    }

    if(s.nodes[2] == 0  && s.nodes[3] == 0 && s.nodes[0] > 0 && s.nodes[1] > 0 ){
      tf += 4;
      return true;
    }

    if(s.nodes[3] == 0  && s.nodes[0] == 0 && s.nodes[1] > 0 && s.nodes[2] > 0 ){
      tf += 4;
      return true;
    }

    //2 mid points, 2 quads
    if(s.nodes[0] == 0  && s.nodes[2] == 0 && s.nodes[1] > 0 && s.nodes[3] > 0 ){

      qf += 2;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[3] == 0 && s.nodes[0] > 0 && s.nodes[2] > 0 ){
      qf += 2;
      return true;
    }

    //1 mid point, 3 trias
    if(s.nodes[0] == 0  && s.nodes[1] == 0 && s.nodes[2] == 0 && s.nodes[3] > 0 ){
      tf += 3;
      return true;
    }

    if(s.nodes[1] == 0  && s.nodes[2] == 0 && s.nodes[3] == 0 && s.nodes[0] > 0 ){
      tf += 3;
      return true;
    }

    if(s.nodes[2] == 0  && s.nodes[3] == 0 && s.nodes[0] == 0 && s.nodes[1] > 0 ){
      tf += 3;
      return true;
    }

    if(s.nodes[3] == 0  && s.nodes[0] == 0 && s.nodes[1] == 0 && s.nodes[2] > 0 ){
      tf += 3;
      return true;
    }
  }
  cerr<<" split() has cases that can not handle" << s.nodes[0] << " "<< s.nodes[1] << " "<<s.nodes[2] << " "<<s.nodes[3] << " "<<s.nodes[4] <<endl;
  return false;
}


void convert2face(store<vector3d<double> > &pos,
                  vector<Array<int,5> > &qfaces, vector<Array<int,4> > &tfaces,
                  vector<Array<int,4> > &tets, vector<Array<int,5> > &pyramids,
                  vector<Array<int,6> > &prisms, vector<Array<int,8> > &hexs,
                  const vector<int>& hex_min, vector<quadSplit>& splits,
                  multiMap &face2node,Map &cl, Map &cr) {

  if(MPI_rank==0)cerr<<"start convert2face" << endl;

  //compute the start cellid of each process
  //  int maxid = 0 ; //std::numeric_limits<int>::min()+2048 ;
  entitySet posDom = pos.domain() ;
  //  if(posDom != EMPTY)
  //    maxid = posDom.Max()+1 ;
  int cellid  = 0;
  //  MPI_Allreduce(&maxid,&cellid,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  //  int cellbase = cellid ;
  int ncells = tets.size()+pyramids.size()+prisms.size()+hexs.size() ;
  vector<int> cellsizes(MPI_processes) ;
  MPI_Allgather(&ncells,1,MPI_INT,&cellsizes[0],1,MPI_INT,MPI_COMM_WORLD) ;

  for(int i=0;i<MPI_rank;++i)
    cellid += cellsizes[i] ;


  //to save memory cost, count num_faces generated by splits
  int num_quads_from_splits = 0;
  int num_tria_from_splits = 0;
  int num_splits =  splits.size();
  for(int i = 0 ; i < num_splits; i++)count_split(splits[i],num_quads_from_splits, num_tria_from_splits);



  int num_quad_faces =
    qfaces.size() + pyramids.size()+ prisms.size()*3 + hexs.size()*6 - num_splits + num_quads_from_splits;
  int num_tria_faces =
    tfaces.size() + tets.size()*4 + pyramids.size()*4 + prisms.size()*2  + num_tria_from_splits;

  if(MPI_rank==0)cerr<<" allocate memory for faces, num_tria_faces: " <<num_tria_faces << " num_quad_faces: " << num_quad_faces<<  endl;

  vector<triaFace> tria(num_tria_faces) ;
  vector<quadFace> quad(num_quad_faces) ;

  if(MPI_rank==0)cerr<<" collecting faces" << endl;

  int tf = 0 ;
  int qf = 0 ;

  //put boundary faces in
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
  // /usr/local/simsys/doc/ug_io/3d_input_output_grids.html provides the SPLIT FACE HEX ELEMENT CONNECTIVITY

  //cellid here and hex_id in .split file are different, hex_id start with 1

  if(split_file_exist) {
    int split_found = 0;
    int start_ind = 0;
    for(size_t i=0;i<hexs.size();++i) {
      int hex_id = hex_min[Loci::MPI_rank] + i;
      for(int face_id = 1; face_id <=6; face_id++){
	//find face
	int ind = find_quad_split(splits, hex_id, face_id, start_ind);
	if(ind < 0){
	  Array<int, 4> corners;
	  get_corners(face_id, hexs[i], corners);
	  quad[qf].nodes[0] = corners[0] ;
	  quad[qf].nodes[1] = corners[1] ;
	  quad[qf].nodes[2] = corners[2] ;
	  quad[qf].nodes[3] = corners[3] ;
	  quad[qf].left = true ;
	  quad[qf++].cell = cellid ;
	}else{
	  split_found++;
	  start_ind = ind;
	  split(splits[ind], quad, qf, tria, tf, cellid);
	}
      }
      cellid++ ;
    }
    //sanity check again
    if((unsigned int)split_found != splits.size()){
      cerr<<" ERROR: not all splits are processed " << endl;
      exit(-1);
    }
  } else {
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
  }
  if(num_quad_faces != qf || num_tria_faces != tf){
    cerr<<"ERROR: face counting not correct" << endl;
    exit(-1);
  }


  //release memory for cells
  vector<Array<int,5> >().swap(qfaces) ;
  vector<Array<int,4> >().swap(tfaces) ;
  vector<Array<int,4> >().swap(tets) ;
  vector<Array<int,5> >().swap(pyramids) ;
  vector<Array<int,6> >().swap(prisms) ;
  vector<Array<int,8> >().swap(hexs) ;

  if(MPI_rank==0)cerr<<" preparing tria faces" << endl;

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


  if(MPI_rank==0)cerr<<" preparing quad faces" << endl;
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

  if(MPI_rank==0)cerr<<" sorting quad faces" << endl;
  //sort quad
  int qsz = quad.size() ;
  int mqsz = 0;
  MPI_Allreduce(&qsz,&mqsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;

  if(mqsz != 0) {
    Loci::parSampleSort(quad,quadCompare,MPI_COMM_WORLD) ;
  }

  if(MPI_rank==0)cerr<<" getting single quad faces" << endl;
  //check if there is single faces in quad, if yes, remove them from quad, and split the single face into tria
  vector<quadFace> single_quad = get_single_face(quad, quadEqual);
  int num_single_quad = single_quad.size();
  //vector<quadFace> single_quad_copy(single_quad);
  int total_num_single_quad = 0;
  MPI_Allreduce(&num_single_quad,&total_num_single_quad,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;


  if(MPI_rank==0)cerr<<" splitting single quad faces" << endl;
  if(num_single_quad != 0){
    split_quad(single_quad, tria);
  }
  if(MPI_rank==0)cerr<<" sorting tria faces" << endl;
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
      cerr << "single triangle faces remain! inconsistent! "<< " single tria " << total_num_single_tria << " single quad " <<  total_num_single_quad << endl ;
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

  if(split_file_exist && MPI_rank==0)cerr<<" collecting edge splits" << endl;

  //collect edge splits before quad splits are modified
  vector<edgeSplit> edge_splits;
  if(split_file_exist)collect_edge_split(edge_splits, splits);
  //release the memory for quad splits
  vector<quadSplit>().swap(splits);

  if(split_file_exist && MPI_processes > 1) {
    vector<pair<int, int> > edges;

    if(MPI_rank == 0)
      cerr << " collecting edges from faces" << endl ;

    collect_edge_from_face(quad,tria, edges );

    if(MPI_rank == 0)
      cerr << " redistribute edge splits" << endl ;

    redistribute_edge_split( edges, edge_splits);

  }


  //due to edge split, trias and quads can turn into general faces
  int ntria = tria.size()/2 ;
  int nquad = quad.size()/2 ;
  int nfaces = ntria+nquad ;
  if(MPI_rank==0)cerr<<" creating face2node, cl, cr" << endl;
  ncells = 0 ;
  for(int i=0;i<MPI_processes;++i)
    ncells += cellsizes[i] ;
  int facebase = std::numeric_limits<int>::min()+2048 ;

  vector<int> facesizes(MPI_processes) ;
  MPI_Allgather(&nfaces,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    facebase += facesizes[i] ;

  if(Loci::MPI_rank == 0) {
    for(int i=0;i<MPI_processes;++i) 
      if(facesizes[i] == 0) {
	cerr << "Run ugrid2vog with fewer than " << i << " processors for this mesh!" << endl ;
	Loci::Abort() ;
	break ;
      }
  }
  entitySet faces = interval(facebase,facebase+nfaces-1) ;

  if(nfaces == 0) {
    faces = EMPTY ;
  }
  store<int> count ;
  cl.allocate(faces) ;
  cr.allocate(faces) ;

  count.allocate(faces) ;
  int fc = facebase ;

  if(split_file_exist) {
    for(int i=0;i<ntria;i++){
      int nedge = 3;
      for(int j =0; j < 3; j++){
	int mp = find_edge_split(edge_splits, tria[i*2].nodes[j],tria[i*2].nodes[(j+1)%3]);
	if(mp >= 0)nedge++;
      }
      count[fc++] = nedge ;
    }

    for(int i=0;i<nquad;i++){
      int nedge = 4;
      for(int j =0; j < 4; j++){
	int mp = find_edge_split(edge_splits, quad[i*2].nodes[j],quad[i*2].nodes[(j+1)%4]);
	if(mp >= 0) nedge++;
      }
      count[fc++] = nedge ;
    }
  } else {
    for(int i=0;i<ntria;i++)
      count[fc++] = 3 ;

    for(int i=0;i<nquad;i++)
      count[fc++] = 4 ;
  }

  face2node.allocate(count) ;
  fc = facebase ;

  for(int i=0;i<ntria;++i) {
    int corner[3];
    for(int j=0;j<3;++j) {
      corner[j] = tria[i*2].nodes[j] ;
      FATAL(tria[i*2].nodes[j] != tria[i*2+1].nodes[j]) ;
    }
    int c1 = tria[i*2].cell ;
    int c2 = tria[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a transparent surface is causing the problem!"<< endl ;
      cerr << "c1=" << c1 << "c2=" << c2 << endl ;
      cerr << "nodes=" ;
      for(int i=0;i<3;++i)
        cerr << " " << tria[i*2].nodes[i] ;
      cerr << endl ;
      Loci::Abort() ;
    }
    if(c1 < 0) {
      cl[fc] = c2 ;
      cr[fc] = c1 ;
      if(tria[i*2+1].left)
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2+1].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2+1].nodes[2-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(tria[i*2].left)
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          corner[j] = tria[i*2].nodes[2-j] ;
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

    if(count[fc] ==3){
      for(int j=0;j<3;++j) {
        face2node[fc][j] = corner[j] ;
      }
    }else{

      int nedge = 0;
      for(int j=0;j<3;++j) {
        face2node[fc][nedge++] = corner[j] ;
        int mp = find_edge_split(edge_splits, corner[j],corner[(j+1)%3]);
        if(mp >=0) face2node[fc][nedge++] = mp;
      }
      if(nedge != count[fc]){
        cerr<<"ERROR: face " << fc << " has count " << count[fc] << " num_node " << nedge << endl;
      }
    }
    fc++ ;
  }

  for(int i=0;i<nquad;++i) {
    int corner[4];
    for(int j=0;j<4;++j) {
      corner[j] = quad[i*2].nodes[j] ;
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
          corner[j] = quad[i*2+1].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          corner[j] = quad[i*2+1].nodes[3-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(quad[i*2].left)
        for(int j=0;j<4;++j)
          corner[j] = quad[i*2].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          corner[j] = quad[i*2].nodes[3-j] ;
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

    if(count[fc] ==4){
      for(int j=0;j<4;++j) {
        face2node[fc][j] = corner[j] ;
      }
    }else{
      int nedge = 0;
      for(int j=0;j<4;++j) {
        face2node[fc][nedge++] = corner[j] ;
        int mp = find_edge_split(edge_splits, corner[j],corner[(j+1)%4]);
        if(mp >=0) face2node[fc][nedge++] = mp;
      }
      if(nedge != count[fc]){
        cerr<<"ERROR: face " << fc << " has count " << count[fc] << " num_node " << nedge << endl;
      }

    }
    fc++ ;
  }
  MPI_Barrier(MPI_COMM_WORLD) ;
  if(MPI_rank==0)cerr<<"done with convert2face" << endl;
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


  size_t nfaces = triangles.size()+nodefacets.size()+bface_sizes.size() ;

  //  entitySet fdom = interval(cellid+ncells,cellid+ncells+nfaces-1) ;
  entitySet fdom = interval(0,nfaces-1) ;
  Map ncl,ncr ;
  store<int> count ;
  ncl.allocate(fdom) ;
  ncr.allocate(fdom) ;
  count.allocate(fdom) ;
  int fbase = fdom.Min() ;
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
  int fcnt = fdom.Min() ;
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
  Loci::Init(&ac,&av) ;
  std::string case_name ;
  bool binary = 0;
  string Lref = "NOSCALE" ;
  bool swapbyte = false ;
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
    } else if(ac >= 2 && !strcmp(av[1],"-b")) {
      binary = 1 ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-swapb")) {
      swapbyte = true ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-o")) {
      optimize = false ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && (!strcmp(av[1],"-cvtransform") ||
                          !strcmp(av[1],"-cvt"))) {
      cellVertexTransform = true ;
      ac-- ;
      av++ ;
      cout << "cellVertexTransformMode active" << endl ;
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

  if(ac == 2) {
    case_name = av[1] ;
  } else {
    cerr << "Usage: ugrid2vog <options> <file>" << endl
         << "Where options are listed below and <file> is the filename sans postfix" << endl
         << "flags:" << endl
         << "  -b  : assume input file is binary" << endl
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

  // Check machines internal byte order
  check_order() ;

  if(swapbyte)
    reverse_byteorder = !reverse_byteorder ;
  if(VOG::hasSuffix(case_name,".b8.ugrid"))
    case_name = VOG::stripSuffix(case_name,".b8.ugrid") ;
  else if(VOG::hasSuffix(case_name,".ugrid"))
    case_name = VOG::stripSuffix(case_name,".ugrid") ;

  if(!binary) {
    struct stat fstat ;
    string ascii_file = case_name + ".ugrid" ;
    if(stat(ascii_file.c_str(),&fstat)<0) {
      binary = true ;
    }
  }

  string infile ;
  if(!binary)
    infile = case_name + ".ugrid" ;
  else
    infile = case_name + ".b8.ugrid" ;

  string outfile = case_name + string(".vog") ;


  store<vector3d<double> > pos ;
  vector<Array<int,5> > qfaces ;
  vector<Array<int,4> > tfaces ;
  vector<Array<int,4> > tets ;
  vector<Array<int,5> > pyramids ;
  vector<Array<int,6> > prisms ;
  vector<Array<int,8> > hexs ;
  readUGRID(infile, binary, pos,qfaces,tfaces,tets,pyramids,prisms,hexs) ;

  string splitfile = case_name + string(".split") ;
  vector<quadSplit> splits;
  int num_hex = hexs.size();
  int P = Loci::MPI_processes;
  vector<int> hex_min(P+1) ;
  MPI_Allgather(&num_hex,1,MPI_INT,&hex_min[1],1,MPI_INT,MPI_COMM_WORLD) ;
  hex_min[0] = 1; //the index starts with 1
  for(int i=1;i<=P;++i)
    hex_min[i] = hex_min[i]+hex_min[i-1] ;


  struct stat buffer;
  int stat_out = 0 ;
  if(Loci::MPI_rank == 0)
    stat_out = stat (splitfile.c_str(), &buffer) ;
  MPI_Bcast(&stat_out,1,MPI_INT,0,MPI_COMM_WORLD) ;
  split_file_exist = (stat_out == 0) ;

  if(split_file_exist)readSplit(splitfile,hex_min, splits);

  if(posScale != 1.0) {
    FORALL(pos.domain(),nd) {
      pos[nd] *= posScale ;
    } ENDFORALL ;
  }

  vector<BC_descriptor> bcs ;
  vector<int> transsurf ;

  if(MPI_rank == 0) {
    string tagsfile = case_name + ".tags" ;
    bcs = readTags(tagsfile) ;
    if(bcs.size() == 0) {
      cerr << "unable to read '" << tagsfile << "'" << endl ;
    }

    for(size_t i=0;i<bcs.size();++i)
      if(bcs[i].Trans)
        transsurf.push_back(bcs[i].id) ;

    cout << "boundary faces:" ;
    for(size_t i=0;i<bcs.size();++i)
      cout << ' ' << bcs[i].name ;
    cout << endl ;

  }

  int trans_size = transsurf.size() ;
  MPI_Bcast(&trans_size,1,MPI_INT,0,MPI_COMM_WORLD) ;
  if(trans_size != 0) {
    if(MPI_rank != 0)
      transsurf = vector<int>(trans_size) ;
    MPI_Bcast(&transsurf[0],trans_size,MPI_INT,0,MPI_COMM_WORLD) ;

    if(MPI_rank == 0)
      cout << "removing transparent surfaces" << endl ;

    vector<Array<int,5> > qtfaces ;
    for(size_t i=0;i<qfaces.size();++i) {
      bool trans = false ;
      for(int j=0;j<trans_size;++j)
        if(qfaces[i][4] == -transsurf[j])
          trans = true ;
      if(!trans)
        qtfaces.push_back(qfaces[i]) ;
    }
    qfaces.swap(qtfaces) ;

    vector<Array<int,4> > ttfaces ;
    for(size_t i=0;i<tfaces.size();++i) {
      bool trans = false ;
      for(int j=0;j<trans_size;++j)
        if(tfaces[i][3] == -transsurf[j])
          trans = true ;
      if(!trans)
        ttfaces.push_back(tfaces[i]) ;
    }
    tfaces.swap(ttfaces) ;
  }

  vector<pair<int,string> > surf_ids ;
  for(size_t i=0;i<bcs.size();++i)
    surf_ids.push_back(pair<int,string>(bcs[i].id,bcs[i].name)) ;




  multiMap face2node ;
  Map cl,cr ;
  if(cellVertexTransform) {
    convert2cellVertexface(pos,qfaces,tfaces,tets,pyramids,prisms,hexs,
                           face2node,cl,cr) ;
  } else {
    convert2face(pos,qfaces,tfaces,tets,pyramids,prisms,hexs, hex_min,splits,
                 face2node,cl,cr) ;
  }


  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  if(MPI_rank == 0)
    cerr << "done  coloring" << endl ;
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
