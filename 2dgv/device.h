#ifndef DEVICE_H
#define DEVICE_H

#include <string>
#include <vector>

using std::vector ;
using std::string ;

#include "grid.h"
#include "pens.h"

struct DevBox {
  double xleft,xright,ytop,ybottom ;
} ;
  

class grDevice {
protected:
  DevBox drawbox ;  
public:
  grDevice() {}
  virtual ~grDevice() {}
  const DevBox &getDrawBox() { return drawbox ;}
  virtual double getTickSize()  = 0 ;
  virtual bool isLabelDev() = 0 ;
  virtual void OpenDevice() = 0 ;
  virtual void CloseDevice() = 0 ;
  virtual void HardCloseDevice() = 0 ;
  virtual void SelectPen(int Pen) = 0 ;
  virtual void DrawLine(const positions &p1,const positions &p2) = 0 ;
  virtual void DrawPoly(const positions *pnts, int num_pnts);
  virtual void DrawLabel(const std::string &text, positions location,
                         int pos, int size, double degrees) = 0 ;
  virtual void FillPoly(const positions *pnts, int num_pnts) = 0 ;
  // Fills a polygon in a fixed color for UI markers (currently black),
  // independent of the pen/color table (which is a 64-entry contour
  // colormap gradient -- see the probe node marker in draw.cc).
  virtual void FillPolyMarker(const positions *pnts, int num_pnts) = 0 ;
  virtual void FillShade(int shade) = 0 ;
} ;


#endif
