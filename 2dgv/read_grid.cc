#include "grid.h"
#include "parse.h"
#include "debug.h"
#include "tools.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <strings.h>
#include <string.h>
#include "pens.h"

using std::ifstream ;
using std::ios ;
using std::cerr ;
using std::cin ;
using std::endl ;

using std::cout ;
using std::cerr ;
using std::endl ;
using std::istream ;
using std::vector ;
using std::list ;

using std::max ;
using std::min ;
using std::string ;
//using std::pow ;
//using std::log10 ;

#define OFFSET_FACTOR 1e-5

list<grid> grids ;

const double epsilon = 1e-30 ;

void grid::optimize_edge_list() {
  const size_t num_nodes = pos.size() ;
  const size_t num_edges = edge_list.size() ;

  vector<int> node2edgec(num_nodes) ;
  for(size_t i=0;i<num_nodes;++i)
    node2edgec[i] = 0 ;
  
  for(size_t i=0;i<num_edges;++i) {
    node2edgec[edge_list[i].l]++ ;
    node2edgec[edge_list[i].r]++ ;
  }

  vector<int> node2edgeo(num_nodes+1) ;

    int count = 0 ;
  for(size_t i=0;i<num_nodes;++i) {
    node2edgeo[i] = count ;
    count += node2edgec[i] ;
  }
  node2edgeo[num_nodes] = count ;
  vector<size_t> node2edge(count) ;

  for(size_t i=0;i<num_edges;++i) {
    size_t l = edge_list[i].l ;
    size_t r = edge_list[i].r ;
    node2edge[node2edgeo[l]+(--node2edgec[l])] = i ;
    node2edge[node2edgeo[r]+(--node2edgec[r])] = i ;
  }
  
  vector<bool> edge_visited(num_edges) ;
  vector<edges> interior_edges,exterior_edges ;

  for(size_t i=0;i<num_nodes;++i) 
    for(int j=node2edgeo[i];j<node2edgeo[i+1];++j) {
      size_t je = node2edge[j] ;
      FATAL(je>=edge_visited.size()) ;
      if(!edge_visited[je]) {
        edge_visited[je] = true ;
        bool unique = true ;
        for(int k=j+1;k<node2edgeo[i+1];++k) 
          if(edge_list[je] == edge_list[node2edge[k]]) {
            interior_edges.push_back(edge_list[je]) ;
            edge_visited[node2edge[k]] = true ;
            unique = false ;
            break ;
          }
        if(unique)
          exterior_edges.push_back(edge_list[je]) ;
      }
    }

    
  interior = interior_edges.size() ;
  if(!pos.empty()) {
    minpos = pos[0] ;
    maxpos = minpos ;
  } else {
    minpos = positions(0,0) ;
    maxpos = positions(0,0) ;
  }
  for(size_t i=0;i<exterior_edges.size();++i) {
    interior_edges.push_back(exterior_edges[i]) ;
    maxpos.x = max(maxpos.x,pos[exterior_edges[i].l].x) ;
    maxpos.x = max(maxpos.x,pos[exterior_edges[i].r].x) ;
    minpos.x = min(minpos.x,pos[exterior_edges[i].l].x) ;
    minpos.x = min(minpos.x,pos[exterior_edges[i].r].x) ;
    maxpos.y = max(maxpos.y,pos[exterior_edges[i].l].y) ;
    maxpos.y = max(maxpos.y,pos[exterior_edges[i].r].y) ;
    minpos.y = min(minpos.y,pos[exterior_edges[i].l].y) ;
    minpos.y = min(minpos.y,pos[exterior_edges[i].r].y) ; 
  }
  minview = minpos ;
  maxview = maxpos ;

  edge_list = interior_edges ;
}

void grid::input_2dgv(istream &in,bool read_values) {
  bool multi = false ;
  valid = true ;
  has_values = read_values ;
  if(parse::is_name(in)) {
    string key = parse::get_name(in) ;
    if(key != "multi") {
      cerr << "improper key in 2dgv file '" << key << "'" << endl ;
      valid = false ;
      return ;
    }
    multi = true ;
  }
  do {

    int ni=-1,nj=-1 ;
    in >> ni ;
    if(in.eof()) break ;
    if(in.peek() == ',')
      in.get() ;
    in >> nj ;
    if(ni<0 || nj<0 || ni>10000 || nj> 10000) {
      cerr << "unreasonable block size (" << ni <<","<<nj<<")"<< endl;
      ni=0;
      nj=0;
      valid = false ;
      return ;
    }
    int pbase = pos.size() ;
    
    for(int i=0;i<ni*nj;++i) {
      double x=0,y=0;
      in >> x ;
      if(in.peek() == ',')
        in.get() ;
      in >> y ;
      pos.push_back(positions(x,y)) ;
      is_node.push_back(true) ;
    }

    if(read_values) {
      for(int i=0;i<ni*nj;++i) {
        double v = 0 ;
        in >> v ;
        val.push_back(v) ;
      }
      
      WARN(vbase != pbase) ;
    }
    

    for(int i=0;i<ni-1;++i)
      for(int j=0;j<nj-1;++j) {
        int n1 = pbase + i + j*ni ;
        int n2 = pbase + (i+1) + j*ni ;
        int n3 = pbase + (i+1) + (j+1)*ni ;
        int n4 = pbase + i + (j+1)*ni ;
        positions center = (pos[n1]+pos[n2]+pos[n3]+pos[n4])*0.25 ;
        double vc = (val[n1]+val[n2]+val[n3]+val[n4])*0.25 ;
        int nc = pos.size() ;
        pos.push_back(center) ;
        is_node.push_back(false) ;
        if(read_values) {
          val.push_back(vc) ;
          triangle_list.push_back(triangles(n1,n2,nc)) ;
          triangle_list.push_back(triangles(n2,n3,nc)) ;
          triangle_list.push_back(triangles(n3,n4,nc)) ;
          triangle_list.push_back(triangles(n4,n1,nc)) ;
        }
        edge_list.push_back(edges(n1,n2)) ;
        edge_list.push_back(edges(n2,n3)) ;
        edge_list.push_back(edges(n3,n4)) ;
        edge_list.push_back(edges(n4,n1)) ;
      }
  } while(multi) ;
}

void grid::input_generalized(std::istream &in,bool read_values) {
  valid = true ;
  has_values = read_values ;
  if(parse::is_name(in)) {
    string key = parse::get_name(in) ;
    if(key != "general") {
      cerr << "improper key in 2dgv file '" << key << "'" << endl ;
      valid = false ;
      return ;
    }
  } else {
    cerr << "confused in input_generalized method" << endl ;
    return ;
  }
  size_t nnodes, nodes_start ;
  in >> nnodes >> nodes_start ;
  pos.reserve(nnodes) ;
  for(size_t i=0;i<nnodes;++i) {
    double x = 0, y = 0 ;
    in >> x >> y ;
    pos.push_back(positions(x,y)) ;
    is_node.push_back(true) ;
  }
  int nfaces, faces_start, ncells, cells_start ;
  in >> nfaces >> faces_start >> ncells >> cells_start ;
  pos.reserve(nnodes+ncells) ;

  edge_list.reserve(nfaces) ;
  vector<edges> bedges ;
  vector<edges> bf2c ;
  vector<edges> f2c ;
  f2c.reserve(nfaces) ;
  bf2c.reserve(nfaces) ;
  bedges.reserve(nfaces) ;
  for(int i=0;i<nfaces;++i) {
    int n1,n2,c1,c2 ;
    in >> n1 >> n2 >> c1 >> c2 ;
    if(c1<0)
      std::swap(c1,c2) ;
    c1 = c1-cells_start+nnodes ;
    c2 = c2<0?c2:c2-cells_start+nnodes ;
    n1 -= nodes_start ;
    n2 -= nodes_start ;
    if(c2 > 0) {
      f2c.push_back(edges(c1,c2)) ;
      edge_list.push_back(edges(n1,n2)) ;
    } else {
      bf2c.push_back(edges(c1,c2)) ;
      bedges.push_back(edges(n1,n2)) ;
    }
  }
  interior = edge_list.size() ;

  if(!pos.empty()) {
    minpos = pos[0] ;
    maxpos = minpos ;
  } else {
    minpos = positions(0,0) ;
    maxpos = positions(0,0) ;
  }
  
  for(size_t i=0;i<bedges.size();++i) {
    maxpos.x = max(maxpos.x,pos[bedges[i].l].x) ;
    maxpos.x = max(maxpos.x,pos[bedges[i].r].x) ;
    minpos.x = min(minpos.x,pos[bedges[i].l].x) ;
    minpos.x = min(minpos.x,pos[bedges[i].r].x) ;
    maxpos.y = max(maxpos.y,pos[bedges[i].l].y) ;
    maxpos.y = max(maxpos.y,pos[bedges[i].r].y) ;
    minpos.y = min(minpos.y,pos[bedges[i].l].y) ;
    minpos.y = min(minpos.y,pos[bedges[i].r].y) ; 

    edge_list.push_back(bedges[i]) ;
    f2c.push_back(bf2c[i]) ;
  }
  minview = minpos ;
  maxview = maxpos ;
    
  val.reserve(nnodes+ncells) ;
  if(parse::is_name(in)) {
    string key = parse::get_name(in) ;
    if(key != "data") {
      cerr << "improper data key in 2dgv file '" << key << "'" << endl ;
      has_values = false ;
      return ;
    }
    int nvec=1, vecnum =1;
    if(parse::is_int(in))
      nvec = parse::get_int(in) ;
    if(parse::is_int(in))
      vecnum = parse::get_int(in) ;
    if(!parse::is_string(in)) {
      cerr << "file name not specified on data statement" << endl ;
      has_values = false ;
      return ;
    }
    string dname = parse::get_string(in) ;
    ifstream data(dname.c_str(),ios::in) ;
    if(data.fail()) {
      cerr << "unable to open data file `" << dname << "'" << endl ;
      has_values = false;
      return ;
    }
    while(data.peek() == '#') {
      int ch ;
      while((ch = data.get()) != EOF && ch != '\n')
        /*NULL STATEMENT */ ;
    }
    vector<double> v(nvec) ;
    for(size_t i=0;i<nnodes;++i) {
      if(data.fail() || data.eof()) {
        has_values = false ;
        cerr << "ran out of values" << endl ;
        break ;
      }
      for(int j=0;j<nvec;++j)
        data >> v[j] ;
      val.push_back(v[vecnum]) ;
    }
  } else
    for(size_t i=0;i<nnodes;++i) {
      if(in.fail() || in.eof()) {
        has_values = false ;
        break ;
      }
      double v ;
      in >> v ;
      val.push_back(v) ;
    }
  
  for(int i=0;i<ncells;++i) {
    pos.push_back(positions(0,0)) ;
    is_node.push_back(false) ;
    if(has_values)
      val.push_back(double(0)) ;
  }
  vector<double> weight(pos.size()) ;
  triangle_list.reserve(nfaces*2) ;
  for(int i=0;i<nfaces;++i) {
    int n1 = edge_list[i].l, n2 = edge_list[i].r ;
    
    positions fc = .5*(pos[n1]+pos[n2]) ;
    double vc = .5*(val[n1]+val[n2]) ;
    double len = sqrt(pos[n1].x*pos[n1].x + pos[n2].y*pos[n2].y) ;
    if(f2c[i].l >= nnodes && f2c[i].l< pos.size()) {
      pos[f2c[i].l] = pos[f2c[i].l] + len*fc ;
      if(has_values)
        val[f2c[i].l] += len*vc ;
      weight[f2c[i].l] += len ;
      triangle_list.push_back(triangles(n1,n2,f2c[i].l)) ;
    }
    if(f2c[i].r >= nnodes && f2c[i].r< pos.size()) {
      pos[f2c[i].r] = pos[f2c[i].r] + len*fc ;
      weight[f2c[i].r] += len ;
      if(has_values)
        val[f2c[i].r] += len*vc ;
      triangle_list.push_back(triangles(n1,n2,f2c[i].r)) ;
    }
  }
  for(size_t i=nnodes;i<pos.size();++i) {
    double rweight = 1./weight[i] ;
    pos[i] = pos[i]*rweight ;
    val[i] *= rweight ;
  }
}

void grid::input_cobalt(std::string filename,bool read_values) {
  valid = true ;
  has_values = read_values ;

  ifstream in(filename.c_str(),ios::in) ;
  if(in.fail()) {
    cerr << "Open failed for " << filename << endl ;
    has_values = false ;
    return ;
  }
  int ndm, nzones, npatch ;
  int npnts, nfaces, ncells, maxppf, maxfpc ;
  in >> ndm >> nzones >> npatch;
  in >> npnts >> nfaces >> ncells >> maxppf >> maxfpc ;

  if(ndm != 2) {
    cerr << "can only read two dimensional cobalt grids!" ;
    has_values = false ;
    return ;
  }
  if(nzones != 1) {
    cerr << "can only read one zone cobalt grids!" ;
    has_values = false ;
    return ;
  }
  
  size_t nnodes, nodes_start ;
  nnodes = npnts ;
  nodes_start = 1 ;
  pos.reserve(nnodes+ncells) ;

  int faces_start, cells_start ;
  faces_start = 1 ;
  cells_start = 1 ;
  
  for(size_t i=0;i<nnodes;++i) {
    double x = 0, y = 0 ;
    in >> x >> y ;
    pos.push_back(positions(x,y)) ;
    is_node.push_back(true) ;
  }

  edge_list.reserve(nfaces) ;
  vector<edges> bedges ;
  vector<edges> bf2c ;
  vector<edges> f2c ;
  f2c.reserve(nfaces) ;
  bf2c.reserve(nfaces) ;
  bedges.reserve(nfaces) ;
  for(int i=0;i<nfaces;++i) {
    int fsz, n1,n2,c1,c2 ;
    in >> fsz ;
    if(fsz != 2) {
      cerr << "two dimensional grids can only have faces with 2 nodes!"
           << endl ;
      return ;
    }
    in >> n1 >> n2 >> c1 >> c2 ;
    if(c1<0)
      std::swap(c1,c2) ;
    c1 = c1-cells_start+nnodes ;
    c2 = c2<0?c2:c2-cells_start+nnodes ;
    n1 -= nodes_start ;
    n2 -= nodes_start ;
    if(c2 > 0) {
      f2c.push_back(edges(c1,c2)) ;
      edge_list.push_back(edges(n1,n2)) ;
    } else {
      bf2c.push_back(edges(c1,c2)) ;
      bedges.push_back(edges(n1,n2)) ;
    }
  }
  interior = edge_list.size() ;

  if(!pos.empty()) {
    minpos = pos[0] ;
    maxpos = minpos ;
  } else {
    minpos = positions(0,0) ;
    maxpos = positions(0,0) ;
  }
  
  for(size_t i=0;i<bedges.size();++i) {
    maxpos.x = max(maxpos.x,pos[bedges[i].l].x) ;
    maxpos.x = max(maxpos.x,pos[bedges[i].r].x) ;
    minpos.x = min(minpos.x,pos[bedges[i].l].x) ;
    minpos.x = min(minpos.x,pos[bedges[i].r].x) ;
    maxpos.y = max(maxpos.y,pos[bedges[i].l].y) ;
    maxpos.y = max(maxpos.y,pos[bedges[i].r].y) ;
    minpos.y = min(minpos.y,pos[bedges[i].l].y) ;
    minpos.y = min(minpos.y,pos[bedges[i].r].y) ; 

    edge_list.push_back(bedges[i]) ;
    f2c.push_back(bf2c[i]) ;
  }
  minview = minpos ;
  maxview = maxpos ;
    
  val.reserve(nnodes+ncells) ;

  string pixfilename ;
  for(int i=filename.size()-1;i>=0;--i) {
    if(filename[i] == '.') {
      for(int j=0;j<i;++j)
        pixfilename += filename[j] ;
      i=0 ;
    }
  }
  if(pixfilename == "") {
    cerr << "unable to get pixfile name!" << endl ;
    has_values = false ;
    return ;
  }
  pixfilename += ".pix" ;
  
  ifstream data(pixfilename.c_str(),ios::in) ;
  if(data.fail()) {
    cout << "unable to open pix file `" << pixfilename << "'" << endl ;
    has_values = false;
    return ;
  }
  while(data.peek() == '#') {
    int ch ;
    while((ch = data.get()) != EOF && ch != '\n')
      /*NULL STATEMENT */ ;
  }
  int pix_nodes,v1,v2,v3,v4 ;
  int nvec = 1 ;
  data >> pix_nodes >> nvec >> v1 ;
  data >> v2 >> v3 >> v4 ;

  for(size_t i=0;i<nnodes;++i) {
    if(data.fail() || data.eof()) {
      has_values = false ;
      cerr << "ran out of values" << endl ;
      break ;
    }
    double vin ;
    data >> vin ;
    val.push_back(vin) ;
  }

  
  for(int i=0;i<ncells;++i) {
    pos.push_back(positions(0,0)) ;
    is_node.push_back(false) ;
    if(has_values)
      val.push_back(double(0)) ;
  }

  vector<double> weight(pos.size()) ;
  triangle_list.reserve(nfaces*2) ;
  for(int i=0;i<nfaces;++i) {
    int n1 = edge_list[i].l, n2 = edge_list[i].r ;
    
    positions fc = (pos[n1]+pos[n2]) ;
    double vc = val[n1]+val[n2] ;
    if(f2c[i].l >= nnodes && f2c[i].l< pos.size()) {
      pos[f2c[i].l] = pos[f2c[i].l] + fc ;
      if(has_values)
        val[f2c[i].l] += vc ;
      weight[f2c[i].l] += 2.0 ;
      triangle_list.push_back(triangles(n1,n2,f2c[i].l)) ;
    }
    if(f2c[i].r >= nnodes && f2c[i].r< pos.size()) {
      pos[f2c[i].r] = pos[f2c[i].r] + fc ;
      weight[f2c[i].r] += 2.0 ;
      if(has_values)
        val[f2c[i].r] += vc ;
      triangle_list.push_back(triangles(n1,n2,f2c[i].r)) ;
    }
  }
  for(size_t i=nnodes;i<pos.size();++i) {
    double rweight = 1./weight[i] ;
    pos[i] = pos[i]*rweight ;
    val[i] *= rweight ;
  }
}

void grid::set_value_range()
{
  min_val = (val.size() == 0)?0:val[0] ;
  max_val = min_val ;
  for(size_t i=0;i<val.size();++i) {
    min_val = min(min_val,val[i]) ;
    max_val = max(max_val,val[i]) ;
  }
}

void grid::input(std::istream &in,bool read_values ) {
  valid = true ;
  int ichar = in.get() ;
  if(in.eof()) {
    cerr << "premature EOF while reading grid" << endl ;
    valid = false ;
  }
  in.putback(ichar) ;
  if(ichar == 'g') {
    input_generalized(in,read_values) ;
    set_value_range() ;
  } else {
    input_2dgv(in,read_values) ;
    set_value_range() ;
    optimize_edge_list() ;
  }
  if(has_values) {
    const double numcontours = 10.0 ;
    double range = (max_val-min_val)/numcontours ;
    double base = range==0?1:pow(10.0,double(floor(log10(range)))) ;
    double contour_lim =  max((fabs(max_val)+fabs(min_val))*1e-8,1e-10) ;
    generate_contour_curves(max(floor(range/base)*base,contour_lim)) ;
    //    generate_shading() ;
  }
  show_contours = true ;
  show_grid = false ;
  show_shading = false ;
}

void grid::input(const char *file,bool read_values) {
  filename = file ;
  char *newfilename = tilde(file) ;
  filename = newfilename ;
  char *rp = rindex(newfilename,'.') ;
  if(rp!=0 && !strcmp(rp,".cfg")) {
    input_cobalt(string(filename),read_values) ;
    set_value_range() ;
    if(has_values) {
      const double numcontours = 10.0 ;
      double range = (max_val-min_val)/numcontours ;
      double base = pow(10.0,double(floor(log10(range)))) ;
      double contour_lim =  max((fabs(max_val)+fabs(min_val))*1e-8,1e-10) ;
      generate_contour_curves(max(floor(range/base)*base,contour_lim)) ;
    }      
    show_contours = true ;
    show_grid = false ;
    show_shading = false ;
  } else {
    ifstream fio(newfilename,ios::in) ;
    if(fio.fail()) {
      cerr << "Open failed for " << newfilename << endl ;
    } else {
      free(newfilename) ;
      input(fio,read_values) ;
    }
  }
}

struct vertigo {
  positions p ;
  double v ;
  vertigo(positions pi, double vi) : p(pi),v(vi) {}
} ;

double sqr(double v) { return v*v ; }

void grid::generate_contour_curves(double cs) {
  contour_spacing = cs ;
  if(!has_values)
    return ;
  int num_contours = int(ceil((max_val-min_val)/contour_spacing)) ;
  int contour_base = int(ceil(min_val/contour_spacing)) ;

  if(num_contours <=0)
    num_contours = 1 ;
  num_contours++ ;

  if(num_contours > 500) {
    double contour_lim =  max((fabs(max_val)+fabs(min_val))*1e-8,1e-10) ;
    contour_spacing =max((max_val-min_val)/500.0,contour_lim) ;
    num_contours = int(ceil((max_val-min_val)/contour_spacing)) ;
    contour_base = int(ceil(min_val/contour_spacing)) ;
    num_contours++ ;
  }





    
  contour_curves.clear() ;
  contour_values.clear() ;
  for(int i=0;i<num_contours;++i) {
    double currval = double(i+contour_base)*contour_spacing
      +contour_spacing*OFFSET_FACTOR ;
    contour_values.push_back(contour_info(currval)) ;
  }
    
  double bignormal = sqrt(sqr(maxview.x-minview.x)+sqr(maxview.y-minview.y)) ;
  
  for(size_t i=0;i<triangle_list.size();++i) {
    const triangles &tp = triangle_list[i] ;
    vertigo v1(pos[tp.t1],val[tp.t1]),v2(pos[tp.t2],val[tp.t2]),
      v3(pos[tp.t3],val[tp.t3]) ;
    if(v1.v > v2.v)
      std::swap(v1,v2) ;
    if(v2.v > v3.v)
      std::swap(v2,v3) ;
    if(v1.v > v2.v)
      std::swap(v2,v1) ;
    if(v1.v == v3.v)
      continue ; // Skip degenerate case of constant value function
    int contour_index = int(ceil(max(v1.v,min_val)/contour_spacing)) - contour_base ;
    double currval = double(contour_index+contour_base)*contour_spacing
      +contour_spacing*OFFSET_FACTOR ;
    const double rdv31 = 1./(v3.v-v1.v+epsilon) ;
    const double rdv32 = 1./(v3.v-v2.v+epsilon) ;
    const double rdv12 = 1./(v1.v-v2.v+epsilon) ;
    
    while(currval <= v3.v) {
      if(currval > max_val)
	break ;
      positions p1 = v1.p + (v3.p-v1.p)*(currval-v1.v)*rdv31 ;
      
      positions p2 ;
      if(currval > v2.v) 
	p2 = v2.p + (v3.p-v2.p)*(currval-v2.v)*rdv32 ;
      else 
	p2 = v2.p + (v1.p-v2.p)*(currval-v2.v)*rdv12 ;
        
      contour_curves.push_back(segments(p1,p2)) ;
      double gradient = (v3.v-v1.v)*bignormal/
	(contour_spacing*sqrt(sqr(v1.p.x-v3.p.x)+sqr(v1.p.y-v3.p.y))) ;
      if(contour_index < int(contour_values.size())) 
	if(gradient < contour_values[contour_index].mingrad) {
	  contour_values[contour_index].mingrad = gradient ;
	  contour_values[contour_index].mingradpos = 0.5*(p1+p2) ;
	}
      contour_index++ ;
      currval += contour_spacing ;
    }      
  }
}

inline void filltriangle(vector<positions> pnts[], vector<int> pntindex[],
                    int cc,
                    const positions &p1, const positions &p2,
                    const positions &p3)
{
  cc = max(min(cc,MAXPENS-1),0) ;
  int indx = pntindex[cc].back() + 3 ;
  pntindex[cc].push_back(indx) ;
  pnts[cc].push_back(p1) ;
  pnts[cc].push_back(p2) ;
  pnts[cc].push_back(p3) ;
}

inline void fillquad(vector<positions> pnts[], vector<int> pntindex[],int cc,
                     const positions &p1, const positions &p2,
                     const positions &p3, const positions &p4)
{
  cc = max(min(cc,MAXPENS-1),0) ;
  int indx = pntindex[cc].back() + 4 ;
  pntindex[cc].push_back(indx) ;
  pnts[cc].push_back(p1) ;
  pnts[cc].push_back(p2) ;
  pnts[cc].push_back(p3) ;
  pnts[cc].push_back(p4) ;
}

inline void fillpent(vector<positions> pnts[], vector<int> pntindex[],int cc,
                     const positions &p1, const positions &p2,
                     const positions &p3, const positions &p4,
                     const positions &p5)
{
  cc = max(min(cc,MAXPENS-1),0) ;
  int indx = pntindex[cc].back() + 5 ;
  pntindex[cc].push_back(indx) ;
  pnts[cc].push_back(p1) ;
  pnts[cc].push_back(p2) ;
  pnts[cc].push_back(p3) ;
  pnts[cc].push_back(p4) ;
  pnts[cc].push_back(p5) ;
}
  
  
                     
void grid::generate_shading() {
  double contour_lim =  max((fabs(max_val)+fabs(min_val))*1e-8,1e-10) ;
  double cspace = max((max_val-min_val)/double(MAXPENS-1),contour_lim) ;
  double contour_base = min_val ;

  for(int i=0;i<MAXPENS+1;++i) {
    pnts[i].clear() ;
    pntindex[i].clear() ;
    pntindex[i].push_back(0) ;
  }
 
  for(size_t i=0;i<triangle_list.size();++i) {
    const triangles &tp = triangle_list[i] ;
    vertigo v1(pos[tp.t1],val[tp.t1]),v2(pos[tp.t2],val[tp.t2]),
      v3(pos[tp.t3],val[tp.t3]) ;
    v1.v = max(min(v1.v,max_val),min_val) ;
    v2.v = max(min(v2.v,max_val),min_val) ;
    v3.v = max(min(v3.v,max_val),min_val) ;
    if(v1.v > v2.v)
      std::swap(v1,v2) ;
    if(v2.v > v3.v)
      std::swap(v2,v3) ;
    if(v1.v > v2.v)
      std::swap(v2,v1) ;
    int curr_contour = int(ceil((v1.v-contour_base+epsilon)/cspace)) ;
    curr_contour = max(min(curr_contour,MAXPENS-1),0) ;
    //    FATAL(curr_contour > MAXPENS) ;
    //    FATAL(curr_contour < 0) ;
    double currval = contour_base + cspace *double(curr_contour)
      +cspace*OFFSET_FACTOR ;
    if(currval >= v3.v) {
      filltriangle(pnts,pntindex,curr_contour,v1.p,v2.p,v3.p) ;
      continue ;
    }

    const double rdv31 = 1./(v3.v-v1.v+epsilon) ;
    const double rdv32 = 1./(v3.v-v2.v+epsilon) ;
    const double rdv12 = 1./(v1.v-v2.v+epsilon) ;

    positions p1 = v1.p + (v3.p-v1.p)*(currval-v1.v)*rdv31 ;
    positions p2 ;
    bool dopent = true ;
    if(currval <= v2.v) {
      p2 = v2.p + (v1.p-v2.p)*(currval-v2.v)*rdv12 ;
      filltriangle(pnts,pntindex,curr_contour,v1.p,p1,p2) ;
    } else {
      p2 = v2.p + (v3.p-v2.p)*(currval-v2.v)*rdv32 ;
      fillquad(pnts,pntindex,curr_contour,v1.p,p1,p2,v2.p) ;
      dopent = false ;
    }
    currval += cspace ;
    curr_contour++ ;
    while(currval <= v2.v && curr_contour <= MAXPENS) {
      positions p3 = v1.p + (v3.p-v1.p)*(currval-v1.v)*rdv31 ;
      positions p4 = v2.p + (v1.p-v2.p)*(currval-v2.v)*rdv12 ;
      fillquad(pnts,pntindex,curr_contour,p1,p2,p4,p3) ;
      p1 = p3 ;
      p2 = p4 ;
      currval += cspace ;
      curr_contour++ ;  
    }
    if(dopent) {
      positions p3 = v1.p + (v3.p-v1.p)*(currval-v1.v)*rdv31 ;
      if(currval >= v3.v) {
        fillquad(pnts,pntindex,curr_contour,p1,p2,v2.p,v3.p) ;
        continue ;
      }
      positions p4 = v2.p + (v3.p-v2.p)*(currval-v2.v)*rdv32 ;
      fillpent(pnts,pntindex,curr_contour,p1,p2,v2.p,p4,p3) ;
      p1 = p3 ;
      p2 = p4 ;
      currval += cspace ;
      curr_contour++ ;
    }
    while(currval <= v3.v && curr_contour <= MAXPENS) {
      positions p3 = v1.p + (v3.p-v1.p)*(currval-v1.v)*rdv31 ;
      positions p4 = v2.p + (v3.p-v2.p)*(currval-v2.v)*rdv32 ;
      fillquad(pnts,pntindex,curr_contour,p1,p2,p4,p3) ;
      p1 = p3 ;
      p2 = p4 ;
      currval += cspace ;
      curr_contour++ ;
    }
    filltriangle(pnts,pntindex,curr_contour,p1,p2,v3.p) ;
  }  
}

void grid::zoom(const grid &gin, const positions &pmn, const positions &pmx) {

  minview = pmn ;
  maxview = pmx ;
  positions pc = 0.5*(pmn+pmx) ;
  double d = sqrt(sqr(pmx.x-pmn.x)+sqr(pmx.y-pmn.y)) ;
  const positions bmx = pmx + 0.1*d*positions(1.,1.) ;
  const positions bmn = pmn - 0.1*d*positions(1.,1.) ;
  std::vector<int> pos_alloc(gin.pos.size()) ;
  for(size_t i=0;i<pos_alloc.size();++i)
    pos_alloc[i] = -1 ;
  interior = ~(unsigned int)(0) ;
  contour_spacing = gin.contour_spacing ;
  minpos = bmn ;
  maxpos = bmx ;
  valid = gin.valid ;
  show_contours = gin.show_contours ;
  show_grid = gin.show_grid ;
  show_shading = gin.show_shading ;
  shading_computed = false ;
  filename = gin.filename ;
  rulers = gin.rulers ;
  has_probe = gin.has_probe ;
  probe_pos = gin.probe_pos ;
  probe_has_val = gin.probe_has_val ;
  probe_val = gin.probe_val ;
  /* Filter edges first */
  for(size_t i=0;i<gin.edge_list.size();++i) {
    if(i==gin.interior)
      interior = edge_list.size() ;
    const edges &e = gin.edge_list[i] ;
    if((gin.pos[e.l].x <= bmx.x && gin.pos[e.l].y <= bmx.y &&
        gin.pos[e.l].x >= bmn.x && gin.pos[e.l].y >= bmn.y) ||
       (gin.pos[e.r].x <= bmx.x && gin.pos[e.r].y <= bmx.y &&
        gin.pos[e.r].x >= bmn.x && gin.pos[e.r].y >= bmn.y)) {
      if(pos_alloc[e.l] == -1) {
        pos_alloc[e.l] = pos.size() ;
        pos.push_back(gin.pos[e.l]) ;
        is_node.push_back(gin.is_node[e.l]) ;
        if(gin.has_values)
          val.push_back(gin.val[e.l]) ;
      }
      if(pos_alloc[e.r] == -1) {
        pos_alloc[e.r] = pos.size() ;
        pos.push_back(gin.pos[e.r]) ;
        is_node.push_back(gin.is_node[e.r]) ;
        if(gin.has_values)
          val.push_back(gin.val[e.r]) ;
      }
      edge_list.push_back(edges(pos_alloc[e.l],pos_alloc[e.r])) ;
    }
  }
  if(interior == ~(unsigned int)(0))
    interior = edge_list.size() ;
  /* Filter triangles */
  for(size_t i=0;i<gin.triangle_list.size();++i) {
    const triangles &t = gin.triangle_list[i] ;
    if((gin.pos[t.t1].x <= bmx.x && gin.pos[t.t1].y <= bmx.y &&
        gin.pos[t.t1].x >= bmn.x && gin.pos[t.t1].y >= bmn.y) ||
       (gin.pos[t.t2].x <= bmx.x && gin.pos[t.t2].y <= bmx.y &&
        gin.pos[t.t2].x >= bmn.x && gin.pos[t.t2].y >= bmn.y) ||
       (gin.pos[t.t3].x <= bmx.x && gin.pos[t.t3].y <= bmx.y &&
        gin.pos[t.t3].x >= bmn.x && gin.pos[t.t3].y >= bmn.y)) {
      if(pos_alloc[t.t1] == -1) {
        pos_alloc[t.t1] = pos.size() ;
        pos.push_back(gin.pos[t.t1]) ;
        is_node.push_back(gin.is_node[t.t1]) ;
        if(gin.has_values)
          val.push_back(gin.val[t.t1]) ;
      }
      if(pos_alloc[t.t2] == -1) {
        pos_alloc[t.t2] = pos.size() ;
        pos.push_back(gin.pos[t.t2]) ;
        is_node.push_back(gin.is_node[t.t2]) ;
        if(gin.has_values)
          val.push_back(gin.val[t.t2]) ;
      }
      if(pos_alloc[t.t3] == -1) {
        pos_alloc[t.t3] = pos.size() ;
        pos.push_back(gin.pos[t.t3]) ;
        is_node.push_back(gin.is_node[t.t3]) ;
        if(gin.has_values)
          val.push_back(gin.val[t.t3]) ;
      }
      triangle_list.push_back(triangles(pos_alloc[t.t1],pos_alloc[t.t2],
                                         pos_alloc[t.t3])) ;
    }
  }
  has_values = (val.size() > 0) ;

  if(has_values) {
    max_val = val[0] ;
    min_val = val[0] ;
    for(size_t i=0;i<val.size();++i) {
      max_val = std::max(val[i],max_val) ;
      min_val = std::min(val[i],min_val) ;
    }
    const double numcontours = 10.0 ;
    double range = (max_val-min_val)/numcontours ;
    double base = range==0?1:pow(10.0,double(floor(log10(range)))) ;
    double contour_lim =  max((fabs(max_val)+fabs(min_val))*1e-8,1e-10) ;
    generate_contour_curves(max(floor(range/base)*base,contour_lim)) ;
    min_val = gin.min_val ;
    max_val = gin.max_val ;
    if(show_shading)
      generate_shading() ;
  }
}

