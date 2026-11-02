#ifndef DRAW_H
#define DRAW_H

#include "device.h"
#include "grid.h"
#include "transform.h"
#include "debug.h"

struct drawfigure {
  const grid &figure ;
  grDevice *currdev ;
  transform2d figtransform ;
  positions titlepos ;
  positions legendpos ;

  drawfigure(const grid &drawgrid) : figure(drawgrid) {currdev = 0;}
  void selectdev(grDevice &dev)  { currdev = &dev ; }
  void scalefigure() ;

  void SelectPen(int p) { currdev->SelectPen(p) ; }
  void DrawLine(positions p1, positions p2)
    { currdev->DrawLine(figtransform.apply_transform(p1),
                        figtransform.apply_transform(p2)); }
  void DrawLabel(const std::string &text, positions location,
                 int pos, int size, double degrees) 
    { currdev->DrawLabel(text,figtransform.apply_transform(location),
                         pos,size,degrees) ; }
  void FillPoly(vector<positions>::const_iterator pi, int npnts) {
    positions pnts[10] ;
    FATAL(npnts>=10) ;
    for(int i=0;i<npnts;++i,++pi) 
      pnts[i] = figtransform.apply_transform(*pi) ;
    currdev->FillPoly(pnts,npnts) ;
  }
    

  void drawtriangles() ;

  void drawborder() ;
  void drawgrid() ;
  void drawcontours() ;
  void drawshading() ;
  void drawtitle() ;
  void drawrulers() ;
  void drawprobe() ;
  void computeViewport(double &xscale, double &yscale,
                        double &viewMinX, double &viewMaxX,
                        double &viewMinY, double &viewMaxY) ;
  void drawlegend() ;
  void draw() ;
} ;

    
void draw(const grid &drawgrid, grDevice &drawdev) ;


#endif
