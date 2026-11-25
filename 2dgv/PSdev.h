#ifndef PSDEV_H
#define PSDEV_H

#include "device.h"
#include "pens.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

class PSDevCore : public grDevice {
protected:
  int currentPen ;
  int currentShade ;
  double r[MAXPENS+1] ;// colors
  double g[MAXPENS+1] ;
  double b[MAXPENS+1] ;
  FILE *outf ;
  //  std::ofstream out ;
  void def_funcs() ;
  void setdrawbox(const DevBox &b) { drawbox = b ; }
public:
  PSDevCore() ;

  virtual ~PSDevCore() {}
  virtual double getTickSize() ;
  virtual void SelectPen(int Pen) ;
  virtual void DrawLine(const positions &p1, const positions &p2) ;
  virtual void DrawLabel(const std::string &text, positions location,
                         int pos, int size, double degrees) ;
  virtual void FillShade(int shade) ;
  virtual void FillPoly(const positions *pnts, int num_pnts) ;
  virtual void FillPolyMarker(const positions *pnts, int num_pnts) ;
} ;


class PSDevEPS : public PSDevCore {
  std::string file ;
public:
  
  PSDevEPS(const DevBox &box,const std::string &filename) ;
  virtual bool isLabelDev()  ;
  virtual void OpenDevice() ;
  virtual void CloseDevice() ;
  virtual void HardCloseDevice() ;
} ;

class PSDevPSFile : public PSDevCore {
  std::string file ;
public:
  PSDevPSFile(const std::string &filename) ;
  virtual bool isLabelDev() ;
  virtual void OpenDevice() ;
  virtual void CloseDevice() ;
  virtual void HardCloseDevice() ;

} ;

class PSDevPipe : public PSDevCore {
  std::string pipe_command ;
public:
  PSDevPipe(const std::string &pcom) ;
  virtual bool isLabelDev() ;
  virtual void OpenDevice() ;
  virtual void CloseDevice() ;
  virtual void HardCloseDevice() ;
} ;
  

#endif
