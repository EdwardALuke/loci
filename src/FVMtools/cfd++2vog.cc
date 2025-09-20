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
#include <fstream>
#include <sstream>
//#include <rpc/rpc.h>
//#include <rpc/xdr.h>
#include <algorithm> 
#include <vector>

using std::vector ;

using std::cout ;
using std::endl ;
using std::cerr ;
using std::ios ;
using std::istringstream ;
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
using std::ifstream ;

void eatSpace(ifstream &file) {
    while(file.peek() != EOF && (isspace(file.peek())))
      file.get() ;
}

void eatComments(ifstream &file) {
  char buf[1024] ;
  buf[1023] = '\0' ;
  eatSpace(file) ;
  while(file.peek() == '#') { // eat comments
    file.getline(buf,1023) ;
    eatSpace(file) ;
  }
}

vector<VOG::BC_descriptor> readCFDppTags(string filename) {
  vector<VOG::BC_descriptor> bcs ;
  ifstream file(filename.c_str(),ios::in) ;
  if(file.fail()) {
    return bcs ;
  }
  string firstline;
  getline(file,firstline);
  if(firstline.find("ANSA") != string::npos) {
    //ANSA generated bc file 
    file.seekg(0); //rewind
    int nlines=0;
    string tmpstr;
    while(getline(file,tmpstr)){
       nlines++;
    }
    int nbcs=(nlines-2)/2;
    cout<<"ANSA File with # BCs = "<<nbcs<<std::endl;
    file.clear();
    file.seekg(0,std::ios::beg);
    getline(file,tmpstr);
    getline(file,tmpstr);
    string name;
    for(int i=0;i<nbcs;++i){
      getline(file,name);
      cout<<name<<endl;
      int id, field2, field3, field4 ;
      file >> id >> field2 >> field3 >> field4;
      getline(file,tmpstr);//read cr lf
      name.erase(0,1);//name has # in front - remove
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
        
      VOG::BC_descriptor BC ;
      BC.id = id ;
      BC.name = name ;
      BC.Trans = false ;
      bcs.push_back(BC) ;
    }
  }
  else {
    eatComments(file) ;
    string command ;
    int nbcs = 0 ;
    file >> command >> nbcs ;
    if(command != "mbcons") {
      cerr << "Boundary condition file format error!" << endl
	     << "Expecting 'mbcons' got '"<<command << "'" << endl ;
      return bcs ;
    }

    for(int i=0;i<nbcs;++i) {
      eatComments(file) ;
      eatSpace(file) ;

      int id, field2, field3, field4 ;
      string name ;
      file >> id >> field2 >> field3 >> field4 >> name ;

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
        
      VOG::BC_descriptor BC ;
      BC.id = id ;
      BC.name = name ;
      BC.Trans = false ;
      bcs.push_back(BC) ;
    }
  }
  return bcs ;
}

void Usage() {
    cerr << "Usage: cfd++2vog <options> <file>" << endl
         << "Where options are listed below and <file> is the filename sans postfix" << endl
         << "flags:" << endl
	 << "  -name <output filename> : set output grid file name" << endl 
	 << "  -help: Give usage information" <<endl
	 << "  -o   : disable grid reordering optimization" << endl 
         << "  -v   : display version" << endl
         << "  -in  : input grid is in inches" << endl
         << "  -ft  : input grid is in feet" << endl
         << "  -cm  : input grid is in centimeters" << endl
         << "  -m   : input grid is in meters" << endl
         << "  -Lref <units> : 1 unit in input grid is <units> long" << endl
         << endl ;

}

int main(int ac, char* av[]) {
  using namespace VOG ;
  Loci::Init(&ac,&av) ;

  bool optimize = true ;
  string gridName = "grid.vog" ;
  string Lref = "m" ; // default to meters for cfd++  
  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 3 && !strcmp(av[1],"-name")) {
      gridName = av[2] ;
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
    } else if(ac >= 2 && !strcmp(av[1],"-help")) {
      Usage() ;
    } else {      
      cerr << "argument " << av[1] << " is not understood." << endl ;
      Usage() ;
      ac-- ;
      av++ ;
    }
  }


  check_order() ;

  double posScale = 1.0 ;

  if(Lref!="m") {
    Loci::UNIT_type tp ;
    istringstream iss(Lref) ;
    iss >> tp ;
    tp.get_value_in("meter") ;
  }


  FILE *NFP=0,*CFP=0,*BFP=0 ;
  int mnodes = 0 ;
  const int R = MPI_rank ;
  const int P = MPI_processes ;
  if(R==0) {
    if((NFP=fopen("nodesin.bin","rb")) == NULL) {
      cerr << "unable to open 'nodesin.bin'" << endl ;
      exit(-1) ;
    }

    // Read in nodes file
    int version = 0 ;
    int rsz = fread(&version, sizeof(int), 1, NFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading nodes" << endl ;
      exit(-1) ;
    }
    rsz = fread(&mnodes, sizeof(int), 1, NFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading nodes" << endl ;
      exit(-1) ;
    }
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&mnodes,sizeof(int),1) ;
    cout << "reading " << mnodes << " nodes" << endl ;
    int nodvar_l = 0 ;
    rsz = fread(&nodvar_l, sizeof(int), 1, NFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading nodes" << endl ;
      exit(-1) ;
    }
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&nodvar_l,sizeof(int),1) ;
  }
    //  cout << "nodvar_l="<<nodvar_l << endl ;
  store<vector3d<double> > pos ;
  MPI_Bcast(&mnodes,1,MPI_INT,0,MPI_COMM_WORLD) ;
  vector<int> node_ptns = VOG::simplePartitionVec(0,mnodes-1,P) ;
  vector<entitySet> local_nodes(P) ;
  for(int i=0;i<P;++i)
    local_nodes[i] = interval(node_ptns[i],node_ptns[i+1]-1) ;

  pos.allocate(local_nodes[R]) ;

  if(R == 0) {
    FORALL(local_nodes[0],nd) {
      int nodtype = 0 ;
      int rsz = fread(&nodtype, sizeof(int), 1, NFP) ;
      if(rsz != 1) {
	cerr << "fread failed reading nodes" << endl ;
	exit(-1) ;
      }    
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&nodtype,sizeof(int),1) ;
      Array<double,3> postmp ;
      rsz = fread(&(postmp[0]),sizeof(double),3,NFP) ;
      if(rsz != 3) {
	cerr << "fread failed reading nodes" << endl ;
	exit(-1) ;
      }
      if(nodtype != 0) {
	cerr << "converter only supports type 0 nodes" << endl ;
	exit(-1) ;
      }

      if(reverse_byteorder) 
	ug_io_reverse_byte_order(&postmp[0],sizeof(double),3) ;
      pos[nd].x = postmp[0] ;
      pos[nd].y = postmp[1] ;
      pos[nd].z = postmp[2] ;
    } ENDFORALL ;
    int mxsize = max(local_nodes[0].size(),
                     local_nodes[P-1].size()) ;
    vector<double> buf(mxsize*3) ;
    for(int i=1;i<P;++i) {
      int cnt = 0 ;
      FORALL(local_nodes[i],nd) {
	int nodtype = 0 ;
	int rsz = fread(&nodtype, sizeof(int), 1, NFP) ;
	if(rsz != 1) {
	  cerr << "fread failed reading nodes" << endl ;
	  exit(-1) ;
	}    
	if(reverse_byteorder)
	  ug_io_reverse_byte_order(&nodtype,sizeof(int),1) ;
	Array<double,3> postmp ;
	rsz = fread(&(postmp[0]),sizeof(double),3,NFP) ;
	if(rsz != 3) {
	  cerr << "fread failed reading nodes" << endl ;
	  exit(-1) ;
	}
	if(nodtype != 0) {
	  cerr << "converter only supports type 0 nodes" << endl ;
	  exit(-1) ;
	}
	
	if(reverse_byteorder) 
	  ug_io_reverse_byte_order(&postmp[0],sizeof(double),3) ;
	buf[cnt++] = postmp[0] ;
	buf[cnt++] = postmp[1] ;
	buf[cnt++] = postmp[2] ;
      } ENDFORALL ;

      MPI_Send(&buf[0],local_nodes[i].size()*3,MPI_DOUBLE,i,9,MPI_COMM_WORLD) ;
    }
    fclose(NFP) ;
    
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
  
  if(posScale != 1.0) {
    FORALL(pos.domain(),nd) {
      pos[nd] *= posScale ;
    } ENDFORALL ;
  }

  if((CFP=fopen("cellsin.bin","rb")) == NULL) {
    cerr << "unable to open 'cellsin.bin'" << endl ;
    exit(-1) ;
  }

  // Face Information
  vector<triaFace> tria ;
  vector<quadFace> quad ;

  int mcells = 0 ;
  if(R==0) {
    // Read in cells file
    int version = 0 ;
    int rsz = fread(&version, sizeof(int), 1, CFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading cells" << endl ;
      exit(-1) ;
    }
    rsz = fread(&mcells, sizeof(int), 1, CFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading cells" << endl ;
      exit(-1) ;
    }
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&mcells,sizeof(int),1) ;
    cout << "reading " << mcells << " cells" << endl ;
    int info_length = 0 ;
    rsz = fread(&info_length, sizeof(int), 1, CFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading cells" << endl ;
      exit(-1) ;
    }
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&info_length,sizeof(int),1) ;
    //  cout << "info_length = " << info_length << endl ;
  }
  MPI_Bcast(&mcells,1,MPI_INT,0,MPI_COMM_WORLD) ;
  

  const int hexmap[6][4] = {{0,4,6,2},{1,3,7,5},{0,1,5,4},
                            {2,6,7,3},{0,2,3,1},{4,5,7,6}} ;
  const int prsmmap[5][4] = {{0,3,5,2},{1,2,5,4},{0,1,4,3},
                             {0,2,1,-1},{3,4,5,-1}} ;
  const int tetmap[4][3] = {{0,2,1},{0,1,3},{1,2,3},{0,3,2}} ;
  const int pyrmap[4][3] = {{0,4,2},{1,3,4},{0,1,4},{2,4,3}} ;
  
  MPI_Bcast(&mcells,1,MPI_INT,0,MPI_COMM_WORLD) ;
  vector<int> cell_ptns = VOG::simplePartitionVec(0,mcells-1,P) ;

  if(R==0) {
    for(int i=cell_ptns[0];i<cell_ptns[1];++i) {
      int celtype = 0 ;
      int rsz = fread(&celtype, sizeof(int), 1, CFP) ;
      if(rsz != 1) {
	cerr << "fread failed reading cells" << endl ;
	exit(-1) ;
      }    
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&celtype,sizeof(int),1) ;
      switch(celtype) {
      case 0: // Hexahedron
	{
	  Array<int,8> hex ;
	  rsz = fread(&hex[0], sizeof(int), 8, CFP) ;
	  if(rsz != 8) {
	    cerr << "fread failed reading hex cell" << endl ;
	    exit(-1) ;
	  }        
	  if(reverse_byteorder)
	    ug_io_reverse_byte_order(&hex[0],sizeof(int),8) ;
	  quadFace qFace ;
	  qFace.cell = i+1 ;
	  qFace.left = true ;
	  
	  for(int f=0;f<6;++f) {
	    qFace.nodes[0] = hex[hexmap[f][0]] ;
	    qFace.nodes[1] = hex[hexmap[f][1]] ;
	    qFace.nodes[2] = hex[hexmap[f][2]] ;
	    qFace.nodes[3] = hex[hexmap[f][3]] ;
	    quad.push_back(qFace) ;
	  }
	}
	break ;
      case 1: // Prism
	{
	  Array<int,6> prsm ;
	  rsz = fread(&prsm[0], sizeof(int), 6, CFP) ;
	  if(rsz != 6) {
	    cerr << "fread failed reading prsm cell" << endl ;
	    exit(-1) ;
	  }
	  if(reverse_byteorder)
	    ug_io_reverse_byte_order(&prsm[0],sizeof(int),6) ;
	  quadFace qFace ;
	  qFace.cell = i+1 ;
	  qFace.left = true ;
	  for(int f=0;f<3;++f) {
	    qFace.nodes[0] = prsm[prsmmap[f][0]] ;
	    qFace.nodes[1] = prsm[prsmmap[f][1]] ;
	    qFace.nodes[2] = prsm[prsmmap[f][2]] ;
	    qFace.nodes[3] = prsm[prsmmap[f][3]] ;
	    quad.push_back(qFace) ;
	  }
	  triaFace tFace ;
	  tFace.left = true ;
	  tFace.cell = i+1 ;
	  
	  for(int f=3;f<5;++f) {
	    tFace.nodes[0] = prsm[prsmmap[f][0]] ;
	    tFace.nodes[1] = prsm[prsmmap[f][1]] ;
	    tFace.nodes[2] = prsm[prsmmap[f][2]] ;
	    tria.push_back(tFace) ;
	  }
	}
	break ;
      case 2: // Tetrahedron
	{
	  Array<int,4> tet ;
	  rsz = fread(&tet[0], sizeof(int), 4, CFP) ;
	  if(rsz != 4) {
	    cerr << "fread failed reading tet cell" << endl ;
	    exit(-1) ;
	  }        
	  if(reverse_byteorder)
	    ug_io_reverse_byte_order(&tet[0],sizeof(int),4) ;
	  triaFace tFace ;
	  tFace.left = true ;
	  tFace.cell = i+1 ;
	  for(int f=0;f<4;++f) {
	    tFace.nodes[0] = tet[tetmap[f][0]] ;
	    tFace.nodes[1] = tet[tetmap[f][1]] ;
	    tFace.nodes[2] = tet[tetmap[f][2]] ;
	    tria.push_back(tFace) ;
	  }
	}
	break ;
      case 6: // Pyramid
	{
	  Array<int,5> pyrm ;
	  rsz = fread(&pyrm[0], sizeof(int), 5, CFP) ;
	  if(rsz != 5) {
	    cerr << "fread failed reading pyrm cell" << endl ;
	    exit(-1) ;
	  }        
	  if(reverse_byteorder)
	    ug_io_reverse_byte_order(&pyrm[0],sizeof(int),5) ;
	  triaFace tFace ;
	  tFace.left = true ;
	  tFace.cell = i+1 ;
	  for(int f=0;f<4;++f) {
	    tFace.nodes[0] = pyrm[pyrmap[f][0]] ;
	    tFace.nodes[1] = pyrm[pyrmap[f][1]] ;
	    tFace.nodes[2] = pyrm[pyrmap[f][2]] ;
	    tria.push_back(tFace) ;
	  }
	  quadFace qFace ;
	  qFace.left = true ;
	  qFace.cell = i+1 ;
	  
	  qFace.nodes[0] = pyrm[0] ;
	  qFace.nodes[1] = pyrm[2] ;
	  qFace.nodes[2] = pyrm[3] ;
	  qFace.nodes[3] = pyrm[1] ;
	  quad.push_back(qFace) ;
	}
	break ;
      default:
	cerr << "unknown cell type " << celtype << endl ;
	break ;
      }
    }
    int mxbufsz =max(cell_ptns[1]-cell_ptns[0],cell_ptns[P]-cell_ptns[P-1]) ;
    vector<int> buf(9*mxbufsz,0) ;
    for(int p=1;p<P;++p) {
      for(int i=cell_ptns[p];i<cell_ptns[p+1];++i) {
	int offset = 9*(i-cell_ptns[p]) ;
	int rsz = fread(&buf[offset], sizeof(int), 1, CFP) ;
	if(rsz != 1) {
	  cerr << "fread failed reading cells" << endl ;
	  exit(-1) ;
	}
	if(reverse_byteorder)
	  ug_io_reverse_byte_order(&buf[offset],sizeof(int),1) ;

	int celtype = buf[offset] ;
	int npnts = 0 ;
	switch(celtype) {
	case 0: // Hexahedron
	  npnts = 8 ;
	  break ;
	break ;
	case 1: // Prism
	  npnts = 6 ;
	  break ;
	break ;
	case 2: // Tetrahedron
	  npnts = 4 ;
	  break ;
	case 6: // Pyramid
	  npnts = 5 ;
	break ;
	default:
	  cerr << "unknown cell type " << celtype << endl ;
	  break ;
	}
      
	rsz = fread(&buf[offset+1], sizeof(int), npnts, CFP) ;
	if(rsz != npnts) {
	  cerr << "fread failed reading cell with " << npnts << " nodes" << endl ;
	  exit(-1) ;
	}        
	if(reverse_byteorder) 
	  ug_io_reverse_byte_order(&buf[offset+1],sizeof(int),npnts) ;
      }
      int lsz = cell_ptns[p+1]-cell_ptns[p] ;
     MPI_Send(&buf[0],lsz*9,MPI_INT,p,10,MPI_COMM_WORLD) ;
    }
    fclose(CFP) ;
  } else {
    MPI_Status status ;
    int ncells = cell_ptns[R+1]-cell_ptns[R] ;
    int recv_count = ncells*9 ;
    vector<int> tmp_cell(recv_count) ;

    MPI_Recv(&tmp_cell[0],recv_count,MPI_INT,0,10,MPI_COMM_WORLD,&status) ;

    int cnt = 0 ;
    for(int i=cell_ptns[R];i<cell_ptns[R+1];++i,++cnt) {
      int offset = cnt*9 ;
      int celtype = tmp_cell[offset] ;
      switch(celtype) {
      case 0: // Hexahedron
	{
	  Array<int,8> hex ;
	  for(int j=0;j<8;++j)
	    hex[j] = tmp_cell[offset+1+j] ;
	  quadFace qFace ;
	  qFace.cell = i+1 ;
	  qFace.left = true ;
	  
	  for(int f=0;f<6;++f) {
	    qFace.nodes[0] = hex[hexmap[f][0]] ;
	    qFace.nodes[1] = hex[hexmap[f][1]] ;
	    qFace.nodes[2] = hex[hexmap[f][2]] ;
	    qFace.nodes[3] = hex[hexmap[f][3]] ;
	    quad.push_back(qFace) ;
	  }
	}
	break ;
      case 1: // Prism
	{
	  Array<int,6> prsm ;
	  for(int j=0;j<6;++j)
	    prsm[j] = tmp_cell[offset+1+j] ;

	  quadFace qFace ;
	  qFace.cell = i+1 ;
	  qFace.left = true ;
	  for(int f=0;f<3;++f) {
	    qFace.nodes[0] = prsm[prsmmap[f][0]] ;
	    qFace.nodes[1] = prsm[prsmmap[f][1]] ;
	    qFace.nodes[2] = prsm[prsmmap[f][2]] ;
	    qFace.nodes[3] = prsm[prsmmap[f][3]] ;
	    quad.push_back(qFace) ;
	  }
	  triaFace tFace ;
	  tFace.left = true ;
	  tFace.cell = i+1 ;
	  
	  for(int f=3;f<5;++f) {
	    tFace.nodes[0] = prsm[prsmmap[f][0]] ;
	    tFace.nodes[1] = prsm[prsmmap[f][1]] ;
	    tFace.nodes[2] = prsm[prsmmap[f][2]] ;
	    tria.push_back(tFace) ;
	  }
	}
	break ;
      case 2: // Tetrahedron
	{
	  Array<int,4> tet ;
	  for(int j=0;j<4;++j)
	    tet[j] = tmp_cell[offset+1+j] ;

	  triaFace tFace ;
	  tFace.left = true ;
	  tFace.cell = i+1 ;
	  for(int f=0;f<4;++f) {
	    tFace.nodes[0] = tet[tetmap[f][0]] ;
	    tFace.nodes[1] = tet[tetmap[f][1]] ;
	    tFace.nodes[2] = tet[tetmap[f][2]] ;
	    tria.push_back(tFace) ;
	  }
	}
	break ;
      case 6: // Pyramid
	{
	  Array<int,5> pyrm ;
	  for(int j=0;j<5;++j)
	    pyrm[j] = tmp_cell[offset+1+j] ;
	  triaFace tFace ;
	  tFace.left = true ;
	  tFace.cell = i+1 ;
	  for(int f=0;f<4;++f) {
	    tFace.nodes[0] = pyrm[pyrmap[f][0]] ;
	    tFace.nodes[1] = pyrm[pyrmap[f][1]] ;
	    tFace.nodes[2] = pyrm[pyrmap[f][2]] ;
	    tria.push_back(tFace) ;
	  }
	  quadFace qFace ;
	  qFace.left = true ;
	  qFace.cell = i+1 ;
	  
	  qFace.nodes[0] = pyrm[0] ;
	  qFace.nodes[1] = pyrm[2] ;
	  qFace.nodes[2] = pyrm[3] ;
	  qFace.nodes[3] = pyrm[1] ;
	  quad.push_back(qFace) ;
	}
	break ;
      default:
	break ;
      }
    }
  }

  if(R==0) {
    if((BFP=fopen("exbcsin.bin","rb")) == NULL) {
      cerr << "unable to open 'exbcsin.bin'" << endl ;
      exit(-1) ;
    }  
    int version = 0 ;
    int rsz = fread(&version, sizeof(int), 1, BFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading boundary faces" << endl ;
      exit(-1) ;
    }
    int mbcs = 0 ;
    rsz = fread(&mbcs, sizeof(int), 1, BFP) ;
    if(rsz != 1) {
      cerr << "fread failed reading boundary faces" << endl ;
      exit(-1) ;
    }
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&mbcs,sizeof(int),1) ;
    cout << "number of boundary faces = " << mbcs << endl ;
    
    for(int i=0;i<mbcs;++i) {
      int bc = 0 ;
      rsz = fread(&bc, sizeof(int), 1, BFP) ;
      if(rsz != 1) {
	cerr << "fread failed reading boundary faces" << endl ;
	exit(-1) ;
      }    
      rsz = fread(&bc, sizeof(int), 1, BFP) ;
      if(rsz != 1) {
	cerr << "fread failed reading boundary faces" << endl ;
	exit(-1) ;
      }
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&bc,sizeof(int),1) ;
      int sz = 0 ;
      rsz = fread(&sz, sizeof(int), 1, BFP) ;
      if(rsz != 1) {
	cerr << "fread failed reading boundary faces" << endl ;
	exit(-1) ;
      }
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&sz,sizeof(int),1) ;
      if(sz == 3) {
	triaFace tFace ;
	tFace.left = true ;
	tFace.cell = -bc ;

	rsz = fread(&tFace.nodes[0], sizeof(int), 3, BFP) ;
	if(rsz != 3) {
	  cerr << "fread failed reading boundary faces" << endl ;
	  exit(-1) ;
	}
	if(reverse_byteorder)
	  ug_io_reverse_byte_order(&tFace.nodes[0],sizeof(int),3) ;
	tria.push_back(tFace) ;
      } else if(sz == 4) {
	quadFace qFace ;
	qFace.left = true ;
	qFace.cell = -bc ;
	rsz = fread(&qFace.nodes[0], sizeof(int), 4, BFP) ;
	if(rsz != 4) {
	  cerr << "fread failed reading boundary faces" << endl ;
	  exit(-1) ;
	}      
	if(reverse_byteorder)
	  ug_io_reverse_byte_order(&qFace.nodes[0],sizeof(int),4) ;
	quad.push_back(qFace) ;
      } else {
	cerr << "face size " << sz << " not supported" << endl ;
      }
    }
    fclose(BFP) ;
  }
  
  vector<BC_descriptor> bcs ;

  if(MPI_rank == 0) {
    string tagsfile = "mcfd.bc" ;
    bcs = readCFDppTags(tagsfile) ;
    if(bcs.size() == 0) {
      cerr << "unable to read '" << tagsfile << "'" << endl ;
    }

    cout << "boundary faces:" ;
    for(size_t i=0;i<bcs.size();++i)
      cout << " '" << bcs[i].name << "'" ;
    cout << endl ;

  }


  // prepare triangle faces (sort them)
  for(size_t i=0;i<tria.size();++i) {
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


  if(MPI_rank==0)cout<<" preparing quad faces" << endl;
  // prepare quad faces (sort them, but be careful)
  for(size_t i=0;i<quad.size();++i) {
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

  if(MPI_rank==0)cout<<" sorting quad faces" << endl;
  //sort quad
  int qsz = quad.size() ;
  int mqsz = 0;
  MPI_Allreduce(&qsz,&mqsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;

  if(mqsz != 0) {
    Loci::parSampleSort(quad,quadCompare,MPI_COMM_WORLD) ;
  }

  if(MPI_rank==0)cout<<" sorting tria faces" << endl;
  //sort tria
  int tsz = tria.size() ;
  int mtsz ;
  MPI_Allreduce(&tsz,&mtsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;

  if(mtsz != 0) {
    Loci::parSampleSort(tria,triaCompare,MPI_COMM_WORLD) ;
  }

  if((tria.size() & 1) == 1) {
    cerr << "non-even number of triangle faces! inconsistent!" << endl ;
    Loci::Abort() ;
  }
  if((quad.size() & 1 ) == 1) {
    cerr << "non-even number of quad faces! inconsistent!" << endl ;
    Loci::Abort() ;
  }

  //due to edge split, trias and quads can turn into general faces
  int ntria = tria.size()/2 ;
  int nquad = quad.size()/2 ;
  int nfaces = ntria+nquad ;
  if(MPI_rank==0)cout<<" creating face2node, cl, cr" << endl;

  //  int ncells = 0 ;
  //  for(int i=0;i<MPI_processes;++i)
  //    ncells += cellsizes[i] ;
  int facebase = std::numeric_limits<int>::min()+2048 ;

  vector<int> facesizes(MPI_processes) ;
  MPI_Allgather(&nfaces,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    facebase += facesizes[i] ;

  if(Loci::MPI_rank == 0) {
    for(int i=0;i<MPI_processes;++i) 
      if(facesizes[i] == 0) {
	cerr << "Run cfd++2vog with fewer than " << i << " processors for this mesh!" << endl ;
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

    for(int j=0;j<3;++j) {
      face2node[fc][j] = corner[j] ;
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

    for(int j=0;j<4;++j) {
      face2node[fc][j] = corner[j] ;
    }
    fc++ ;
  }
  MPI_Barrier(MPI_COMM_WORLD) ;
  if(MPI_rank==0)cout<<"done with convert2face" << endl;
  if(MPI_rank == 0)
    cout << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;
  MPI_Barrier(MPI_COMM_WORLD) ;
  if(MPI_rank == 0)
    cout << "done  coloring" << endl ;
  if(optimize) {
    if(MPI_rank == 0)
      cout << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }

  if(MPI_rank == 0)
    cout << "writing VOG file to " << gridName << endl ;

  vector<pair<int,string> > surf_ids ;
  for(size_t i=0;i<bcs.size();++i)
    surf_ids.push_back(pair<int,string>(bcs[i].id,bcs[i].name)) ;

  string outfile = gridName ;
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;

  Loci::Finalize() ;
  return 0 ; 
}
