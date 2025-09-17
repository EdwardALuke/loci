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
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std ;
using namespace Loci ;
namespace Loci {
  inline bool operator<(const Array<int,3> &a1,
                        const Array<int,3> &a2) {
    return a1[0]<a2[0] || (a1[0] == a2[0] && a1[1] < a2[1]) ||
      (a1[0] == a2[0] && a1[1] == a2[1] && a1[2] < a2[2]) ;
  }
  inline bool operator<(const Array<int,2> &a1,
                        const Array<int,2> &a2) {
    return a1[0]<a2[0] || (a1[0] == a2[0] && a1[1] < a2[1]) ;
  }
}

void kill_comment(istream &s) {
  while(s.peek() == ' ' || s.peek() == '\t' ||
        s.peek() == '\n' || s.peek() == '\r')
    s.get() ;
  
  while(s.peek() == '#') {
    char buf[1024] ;
    s.getline(buf,1023) ;
    while(s.peek() == ' ' || s.peek() == '\t' ||
          s.peek() == '\n' || s.peek() == '\r')
      s.get() ;
  }
}

struct bc_info {
  string name ;
  vector<pair<int,int> > edge_list ;
} ;

void Usage() {
  cerr << "Usage: extruder <options> surface.gsurf" << endl ;
  cerr << endl ;
  cerr << "This program can be used to create a cobalt grid by extruding a "
       << endl
       << "generalized surface mesh (.gsurf) file that can be extracted from"
       << endl
       << "a case run using extract -surf -bc <surfid> <casename> <iteration>"
       << endl ;
  cerr << endl ;
  cerr << "options to extruder are as follows:" << endl
       << endl
       << "  -delta <float value> " << endl
       << "     specifies the extusion direction delta (default 0.1)" <<endl
       << "  -growth <float value> " << endl
       << "     specifies the growth rate of delta (defaut 1.0)" << endl
       << "  -deltamax <float value> " << endl
       << "     maximum amount delta can grow (default 1000.0)" << endl
       << "  -nplanes <integer value> " << endl
       << "     number of planes to extrude (default 2)" << endl
       << "  -nx <float value>" << endl
       << "     extrusion normal vector x component (default 0.0)" << endl
       << "  -ny <float value>" << endl
       << "     extrusion normal vector y component (default 0.0)" << endl
       << "  -nz <float value>" << endl
       << "     extrusion normal vector z component (default 1.0)" << endl
       << "  -xr"<< endl
       << "     project cylindrical coordinates onto xy surface" << endl 
       << "  -no_hanging_nodes" << endl
       << "     mesh does not contain hanging nodes from mesh adaptation so don't merge edges" << endl ;
  cerr << endl ;
  exit(-1) ;
}

int main(int ac, char *av[]) {
  enum {NONE, XY_CUT,YZ_CUT,XZ_CUT,XR_CUT} cut_type =NONE;
  double nx=0,ny=0,nz=1.0 ;
  double delta = 0.1 ;
  double deltamax = 1000 ;
  double growth = 1.0 ;
  int nplanes = 2 ;
  bool use_delta_file = false ;
  bool merge_hanging_nodes = true ;
  bool radial = false ;
  string delta_file = "" ;
  while(ac > 2 ) {
    if(!strcmp(av[1],"-xy")) {
      cut_type = XY_CUT ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-yz")) {
      cut_type = YZ_CUT ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-xz")) {
      cut_type = XZ_CUT ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-xr")) {
      cut_type = XR_CUT ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-nx") && ac > 3) {
      ac-- ;
      av++ ;
      nx = atof(av[1]) ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-ny") && ac > 3) {
      ac-- ;
      av++ ;
      ny = atof(av[1]) ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-nz") && ac > 3) {
      ac-- ;
      av++ ;
      nz = atof(av[1]) ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-radial") && ac > 2) {
      radial = true ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-delta_file") && ac > 3) {
      ac-- ;
      av++ ;
      delta_file = av[1] ;
      use_delta_file = true ;
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-delta") && ac > 3) {
      ac-- ;
      av++ ;
      delta = atof(av[1]) ;
      if(delta <= 0.0) {
        cerr << "delta must be positive!" << endl ;
        Usage() ;
      }
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-deltamax") && ac > 3) {
      ac-- ;
      av++ ;
      deltamax = atof(av[1]) ;
      if(deltamax <=0.0) {
        cerr << "deltamax must be positive!" << endl ;
        Usage() ;
      }
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-growth") && ac > 3) {
      ac-- ;
      av++ ;
      growth = atof(av[1]) ;
      if(growth <=0.0) {
        cerr << "growth must be positive!" << endl ;
        Usage() ;
      }
      
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-nplanes") && ac > 3) {
      ac-- ;
      av++ ;
      nplanes = atoi(av[1]) ;
      if(nplanes < 2) {
        cerr << "nplanes must be two or greater!" << endl ;
        Usage() ;
      }
      ac-- ;
      av++ ;
    } else if(!strcmp(av[1],"-no_hanging_nodes") && ac > 2) {
      ac-- ;
      av++ ;
      merge_hanging_nodes = false ;
    } else {
      cerr << "unknown argument '" << av[1] << endl ;
      Usage() ;
      break ;
    }
  } 
  vector<double> delta_list ;
  if(use_delta_file) {
    ifstream dfile(delta_file.c_str(),ios::in); 
    int ndeltas = 0 ;
    dfile >> ndeltas ;
    if(ndeltas <= 0) {
      cerr << "invalid delta_file provided" << endl ;
      exit(-1) ;
    }
    vector<double> tmp(ndeltas) ;
    for(int i=0;i<ndeltas;++i) {
      dfile >> tmp[i] ;
    }
    delta_list = tmp ;
    nplanes = ndeltas+1 ;
  }
    
      
  double nnorm = sqrt(nx*nx+ny*ny+nz*nz) ;
  nx *= 1./nnorm ;
  ny *= 1./nnorm ;
  nz *= 1./nnorm ;
  if(ac != 2) {
    cerr << "takes surface file as input" << endl ;
    Usage() ;
    return -1 ;
  }

  ifstream isurf(av[1],ios::in) ;

  if(isurf.fail()) {
    cerr << "unable to open file " << av[1] << endl ;
    return -1 ;
  }

  kill_comment(isurf) ;
  int npnts = 0 ;
  isurf >> npnts ;
  kill_comment(isurf) ;
  vector<double> x(npnts),y(npnts),z(npnts) ;
  for(int i=0;i<npnts;++i)
    isurf >> x[i] >> y[i] >> z[i] ;
  kill_comment(isurf) ;
  int nfaces = 0 ;
  isurf >> nfaces ;
  vector<int> f2size(nfaces) ;
  vector<int> f2node ;
  kill_comment(isurf) ;
  for(int i=0;i<nfaces;++i) {
    isurf >> f2size[i] ;
    for(int j=0;j<f2size[i];++j) {
      int nd = 0 ;
      isurf >> nd ;
      f2node.push_back(nd) ;
    }
  }

  kill_comment(isurf) ;
  int nbcs = 0;
  isurf >> nbcs ;
  vector<bc_info> bcs(nbcs) ;


  // Get grid edges
  vector<Array<int,3> > edge_map ;
  
  for(int i=0;i<nbcs;++i) {
    kill_comment(isurf) ;
    isurf >> bcs[i].name ;
    cout << "BC: " << bcs[i].name << " - " << (i+1) << endl ;
    kill_comment(isurf) ;
    int nedges =0 ;
    isurf >> nedges ;
    kill_comment(isurf) ;
    for(int j=0;j<nedges;++j) {
      int n1, n2 ;
      isurf >> n1 >> n2 ;
      Array<int,3> edge_info ;
      edge_info[0] = min(n1,n2) ;
      edge_info[1] = max(n1,n2) ;
      edge_info[2] = -(i+1) ;
      edge_map.push_back(edge_info) ;
      bcs[i].edge_list.push_back(pair<int,int>(n1,n2)) ;
    }
  }

  int off = 0 ;
  for(int i=0;i<nfaces;++i) {
    int fsz = f2size[i] ;
    Array<int,3> edge_info ;
    int n1 = f2node[off] ;
    int n2 = f2node[off+fsz-1] ;
    edge_info[0] = min(n1,n2) ;
    edge_info[1] = max(n1,n2) ;
    edge_info[2] = i+1 ; // Cell number
    edge_map.push_back(edge_info) ;
    for(int j=1;j<fsz;++j) {
      n1 = f2node[off+j-1] ;
      n2 = f2node[off+j] ;
      edge_info[0] = min(n1,n2) ;
      edge_info[1] = max(n1,n2) ;
      edge_info[2] = i+1 ;
      edge_map.push_back(edge_info) ;
    }
    off += fsz ;
  }

  if((edge_map.size() & 1) != 0) {
    cerr << "unmatched edges in surface?  Are all boundaries set?" << endl ;
    return -1 ;
  }
  sort(edge_map.begin(),edge_map.end()) ;

  vector<Array<int,2> > boundary_edges ;


  int emsz = edge_map.size() ;
  for(int i=0;i<emsz;i+=2) {
    if((edge_map[i][0] != edge_map[i+1][0] ||
	edge_map[i][1] != edge_map[i+1][1])) {
      cerr << "unmatched edge in surface?  Are all boundaries set?" << endl ;
    }
    if(edge_map[i][2] == edge_map[i+1][2]) {
      cerr << "same cell on both sides of face?" << endl ;
    }
    if(edge_map[i][2] < 0 && edge_map[i+1][2] < 0) {
      cerr << "bc setting on both sides of face?" << endl ;
    }
    if(edge_map[i][2] < 0) {
      Array<int,2> tmp ;
      tmp[0] = edge_map[i][0] ;
      tmp[1] = i ;
      boundary_edges.push_back(tmp) ;
      tmp[0] = edge_map[i][1] ;
      boundary_edges.push_back(tmp) ;
    }
    if(edge_map[i][0] == 0 || edge_map[i+1][0] == 0 ||
       edge_map[i][1] == 0 || edge_map[i+1][1] == 0 ||
       edge_map[i][2] == 0 || edge_map[i+1][2] == 0) {
      cerr << "zero node or cell number" << endl ;
    }
  }
  
  sort(boundary_edges.begin(),boundary_edges.end()) ;
  vector<int> edgepairs ;
  vector<int> edge_map_bpairs(edge_map.size(),0) ;
  if(merge_hanging_nodes) {
    for(size_t i=0;i<boundary_edges.size(); i+=2) {
      if(boundary_edges[i][0] != boundary_edges[i+1][0]) {
	cerr << "dangling end of boundary edges!" << endl ;
      }
      int b1 = boundary_edges[i][1] ;
      int b2 = boundary_edges[i+1][1] ;
      if(edge_map_bpairs[b1] ==0 && edge_map_bpairs[b2] == 0 &&
	 edge_map[b1][2] == edge_map[b2][2] &&
	 edge_map[b1+1][2] == edge_map[b2+1][2]) {
	edge_map_bpairs[b1] = 1 ;
	edge_map_bpairs[b1+1] = 1 ;
	edge_map_bpairs[b2] = 1 ;
	edge_map_bpairs[b2+1] = 1 ;
	edgepairs.push_back(i) ;
      }
    }
  }

  cout << "num boundary edges = " << boundary_edges.size()/2 << endl ;
  if(edgepairs.size() > 0)
    cout << "paired boundary edges= " << edgepairs.size() << endl ;
  map<int,pair<int,int> > node2bedge ;
  ofstream ofile("grid.cog",ios::out) ;
  if(ofile.fail()) {
    cerr << "can't open grid.cog" << endl ;
    return  -1 ;
  }

  int splitEdgeCnt = edgepairs.size() ;
  ofile.precision(16) ;
  int npatch = nbcs+2 ;
  ofile << "3 1 " << npatch << endl ;
  int vol_nfaces = nfaces*nplanes + (emsz/2-splitEdgeCnt)*(nplanes-1) ;
  int vol_ncells = nfaces*(nplanes-1) ;
  int vol_npnts = npnts*nplanes ;
  int mxppf = 4 ;
  int mxfpc = 6 ;
  ofile << vol_npnts<< " " << vol_nfaces << " " << vol_ncells << " " << mxppf << " " << mxfpc << endl ;

  // output positions
  for(int i=0;i<npnts;++i) {
    double X,Y,Z ;
    switch(cut_type) {
    case XY_CUT:
      X = x[i] ;
      Y = y[i] ;
      Z = -0.5 ;
      break ;
    case YZ_CUT:
      X = y[i] ;
      Y = z[i] ;
      Z = -0.5 ;
      break ;
    case XZ_CUT:
      X = x[i] ;
      Y = z[i] ;
      Z = -0.5 ;
      break ;
    case XR_CUT:
      X = x[i] ;
      Y = sqrt(y[i]*y[i]+z[i]*z[i]) ;
      Z = -0.5 ;
      break ;
    case NONE:
    default:
      X = x[i] ;
      Y = y[i] ;
      Z = z[i] ;
      break ;      
    }
    x[i] = X ;
    y[i] = Y ;
    z[i] = Z ;
    ofile << x[i] << ' ' << y[i] << ' ' << z[i] << endl ;
  }
  double dn = delta ;
  double ddn = dn ;
  if(use_delta_file)
    dn = delta_list[0] ;
  for(int j=1;j<nplanes;++j) {
    double dx = nx*dn ;
    double dy = ny*dn ;
    double dz = nz*dn ;
    double dnsave = dn ;
    ddn *= growth ;
    ddn = min(ddn,deltamax) ;
    if(use_delta_file) {
      if(j+1 != nplanes) 
	dn += delta_list[j] ;
    } else {
      dn += ddn ;
    }
    if(radial) {
      for(int i=0;i<npnts;++i) {
	double rlen = 1./sqrt(y[i]*y[i]+z[i]*z[i]) ;
	double dx = 0 ;
	double dy = dnsave*rlen*y[i] ;
	double dz = dnsave*rlen*z[i] ;
	ofile << x[i]+dx << ' ' << y[i]+dy << ' ' << z[i]+dz << endl ;
      }	
    } else {
      for(int i=0;i<npnts;++i) {
	ofile << x[i]+dx << ' ' << y[i]+dy << ' ' << z[i]+dz << endl ;
      }
    }
  }
  // now output mesh faces, starting with edge extrusion
  
  for(int j=1;j<nplanes;++j) {
    int of1 = (j-1)*npnts ;
    int of2 = j*npnts ;
    int cof = (j-1)*nfaces ;
    //    int cnt = 0;
    for(int i=0;i<emsz;i+=2) {
      if(edge_map_bpairs[i] == 0) {
	// all extruded faces are quads
	ofile << 4 << ' ' ;
	ofile << edge_map[i][0]+of1 << ' ' << edge_map[i][1]+of1 << ' '
	      << edge_map[i][1]+of2 << ' ' << edge_map[i][0]+of2 << ' ' ;
	int c1 = edge_map[i+1][2] ; 
	int c2 = edge_map[i][2] ;
	c1 += (c1<0)?0:cof ;
	c2 += (c2<0)?0:cof ;
	ofile << c1 << ' ' << c2 << endl ;
      } //else
	//	cnt++ ;
    }

    for(size_t i=0;i<edgepairs.size();++i) {
      int ent = edgepairs[i] ;
      int ncent = boundary_edges[ent][0] ;
      int b1 = boundary_edges[ent][1] ;
      int b2 = boundary_edges[ent+1][1] ;
      int n1 = edge_map[b1][0]==ncent?edge_map[b1][1]:edge_map[b1][0] ;
      int n2 = edge_map[b2][0]==ncent?edge_map[b2][1]:edge_map[b2][0] ;
      // all extruded  split edges have quads with 6 edges (2 split)
	ofile << 6 << ' ' ;
	ofile << n1+of1 << ' ' << ncent+of1 << ' '
	      << n2+of1 << ' ' << n2+of2 << ' ' 
	      << ncent+of2 << ' ' << n1+of2 << ' ' ;
	int c1 = edge_map[b1+1][2] ; 
	int c2 = edge_map[b1][2] ;
	c1 += (c1<0)?0:cof ;
	c2 += (c2<0)?0:cof ;
	ofile << c1 << ' ' << c2 << endl ;
    }
  }
  entitySet nodes1 ;
  // Now build symmetry planes
  int mnfsz = 100000000 ;
  int mxfsz = 0 ;
  int bc = -nbcs-1 ;
  off = 0 ;
  int cplane = 0 ;
  for(int i=0;i<nfaces;++i) {
    int fsz = f2size[i] ;
    mxfsz = max(mxfsz,fsz) ;
    mnfsz = min(mnfsz,fsz) ;
    ofile << fsz  ;
    for(int j=0;j<fsz;++j)
      ofile << ' ' << f2node[off+j] ;
    ofile << ' ' << i+1+cplane << ' ' << bc << endl ;
    for(int j=0;j<fsz;++j)
      nodes1 += f2node[off+j] ;
    off += fsz ;
  }
  for(int j=1;j<nplanes-1;++j) {
    cplane++ ;
    off = 0 ;
    for(int i=0;i<nfaces;++i) {
      int fsz = f2size[i] ;
      mxfsz = max(mxfsz,fsz) ;
      mnfsz = min(mnfsz,fsz) ;
      ofile << fsz  ;
      for(int j=0;j<fsz;++j)
        ofile << ' ' << f2node[off+j] +npnts*cplane;
      ofile << ' ' << i+1+(cplane-1)*nfaces
            << ' ' << i+1+cplane*nfaces << endl ;
      for(int j=0;j<fsz;++j)
        nodes1 += f2node[off+j]+npnts*cplane ;
      off += fsz ;
    }
  }    
  cplane++ ;
  bc -= 1 ;
  off = 0 ;
  for(int i=0;i<nfaces;++i) {
    int fsz = f2size[i] ;
    ofile << fsz  ;
    for(int j=0;j<fsz;++j)
      ofile << ' ' << f2node[off+j]+npnts*cplane ;
    ofile << ' ' << i+1+(cplane-1)*nfaces << ' ' << bc << endl ;
    for(int j=0;j<fsz;++j)
      nodes1 += f2node[off+j]+npnts*cplane ;
    off += fsz ;
  }
  cout << "nodes1= " << nodes1 << endl ;
  cout << "mxfsz = " << mxfsz << ", mnfsz = " << mnfsz << endl ;
  
  ofstream otfile("grid.tags",ios::out) ;
  if(otfile.fail()) {
    cerr << "can't open grid.tags" << endl ;
    return  -1 ;
  }
  otfile << "#ID:    Group                   BC      Visc    Recon   Source  Trans   Rebuild" << endl ;
  for(int i=0;i<nbcs;++i) {
    otfile << "#" << i+1 << ":     " << bcs[i].name << "      "
           << " 1 0 0 0 0 0" << endl ;
  }
  otfile << "#" << nbcs+1 << ":     " << "plane1" << "      "
           << " 1 0 0 0 0 0" << endl ;
  otfile << "#" << nbcs+2 << ":     " << "plane2" << "      "
           << " 1 0 0 0 0 0" << endl ;
  otfile.close() ;
    
  
  return 0;
}
