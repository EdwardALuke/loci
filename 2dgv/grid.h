#ifndef GRID_H
#define GRID_H

#include "pens.h"

#include <vector>
#include <iostream>
#include <list>
#include <string>
#include <math.h>

struct edges {
  size_t l,r ;
  edges(size_t li,size_t ri) : l(li),r(ri) {}
  edges() {}
  bool operator==(const edges &e)
    { return (l == e.l && r == e.r) || (l == e.r && r == e.l) ;}
} ;

struct triangles {
  int t1,t2,t3 ;
  triangles(int ti1,int ti2, int ti3) :
    t1(ti1),t2(ti2),t3(ti3) {}
  triangles() {}
} ;

struct positions {
  double x,y ;
  positions(double xi,double yi): x(xi),y(yi) {}
  positions() {}
} ;

inline positions operator+(const positions &p1, const positions &p2)
{ return positions(p1.x+p2.x,p1.y+p2.y) ; }

inline positions operator-(const positions &p1, const positions &p2)
{ return positions(p1.x-p2.x,p1.y-p2.y) ; }

inline positions operator*(const positions &p, double s)
{ return positions(p.x*s,p.y*s) ; }

inline positions operator*(double s,const positions &p)
{ return positions(p.x*s,p.y*s) ; }

inline positions operator/(const positions &p, double s)
{ return positions(p.x/s,p.y/s) ; }

inline positions operator/(double s,const positions &p)
{ return positions(p.x/s,p.y/s) ; }

struct segments {
  positions p1,p2 ;
  segments(positions p1i,positions p2i) : p1(p1i),p2(p2i) {}
} ;

struct contour_info {
  positions mingradpos ;
  double mingrad ;
  double value ;
  contour_info() {mingrad = 1e30; mingradpos = positions(0,0) ;}
  contour_info(double val)
    {mingrad = 1e30; mingradpos = positions(0,0) ; value = val ;}
} ;

struct rule {
  positions start ;
  double length ;
  bool xdir ;
  rule(const positions &p1,const positions &p2) {
    if(fabs(p1.x-p2.x) > fabs(p1.y-p2.y)) {
      xdir = true ;
      double y = 0.5*(p1.y+p2.y) ;
      start = positions(std::min(p1.x,p2.x),y) ;
      length = fabs(p1.x-p2.x) ;
    } else {
      xdir = false ;
      double x = 0.5*(p1.x+p2.x) ;
      start = positions(x,std::min(p1.y,p2.y)) ;
      length = fabs(p1.y-p2.y) ;
    }
  }
} ;

struct grid {
  std::vector<positions> pos ;
  positions minpos,maxpos ;
  positions minview,maxview ;

  bool has_values ;
  double min_val,max_val ;
  std::vector<double> val ;
  std::vector<bool> is_node ;
  
  std::vector<edges> edge_list ;
  unsigned int interior ;

  std::vector<triangles> triangle_list ;


  std::string filename ;

  double contour_spacing ;
  std::vector<segments> contour_curves ;
  std::vector<contour_info> contour_values ;


  bool shading_computed ;
  std::vector<positions>  pnts[MAXPENS+1] ;
  std::vector<int>        pntindex[MAXPENS+1] ;
  
  bool show_contours, show_grid, show_shading, valid ;

  std::vector<rule> rulers ;

  // Right-click node probe: the nearest mesh node to the last probe
  // click, its scalar value (if the grid has one), and its position.
  // Follows the same pattern as rulers -- lives on the grid so it's
  // drawn on every redraw and carried across zoom (see grid::zoom()),
  // not just the instant after the click.
  bool has_probe ;
  positions probe_pos ;
  bool probe_has_val ;
  double probe_val ;

  grid() {
    show_contours = false, show_grid = false, show_shading = false;
    shading_computed = false ; has_values = false ; filename = "none";
    valid = false ; has_probe = false ;  }
  void optimize_edge_list() ;
  void generate_contour_curves(double cs) ;
  void generate_shading() ;

  void set_value_range() ;
  void reset_value_range(double v1,double v2) {
    min_val = std::min(v1,v2) ;
    max_val = std::max(v1,v2) ;
  }

  void zoom(const grid &gin, const positions &pmn, const positions &pmx) ;

  void input_2dgv(std::istream &in,bool read_values) ;
  void input_generalized(std::istream &in,bool read_values) ;
  void input_cobalt(std::string filename, bool read_values) ;
  void input(std::istream &in,bool read_vadrlues) ;
  void input(const char *filename, bool read_values = true) ;
} ;

extern std::list<grid> grids ;

#endif 
