#ifndef XDEV_H
#define XDEV_H
#include <X11/Xlib.h>
#include <X11/Xos.h>

#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/StringDefs.h>


#include <X11/Intrinsic.h>
#include <X11/IntrinsicP.h>
#include <X11/CoreP.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Box.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Shell.h>
#include <X11/Xaw/List.h>
#include <X11/Xaw/Label.h>
#include <X11/Xaw/Paned.h>

#include "device.h"

extern int X_R[],X_G[],X_B[] ;

class XDevice : public grDevice {
  friend void mouseButtonPress(Widget w, XtPointer *, XEvent *, Boolean *) ;
  friend void Zoomin(Widget w, XtPointer p1, XtPointer p2) ;
  friend void ToggleRuler(Widget w, XtPointer p1, XtPointer p2) ;
  int currentPen ;
  int currentShade ;
  int drawwidth,drawheight ;
  Pixmap drawpixmap ;
  bool pixmap_allocated ;
  Window drawwindow ;
  int drawdepth ;
  Display *dpy ;
  GC linec ;
  GC forec ;
  GC Penc[MAXPENS] ;
  GC markerGC ; // dedicated GC for the probe node marker (see realize())
  bool drawnbox ;
  int bx1,by1,bx2,by2 ;
  double ticksize ;
public:
  void realize(Widget w) ;
  void Xdrawbox() {
    if(drawnbox) {
      XDrawLine(dpy,drawwindow,Penc[0],bx1,by1,bx2,by1) ;
      XDrawLine(dpy,drawwindow,Penc[0],bx2,by1,bx2,by2) ;
      XDrawLine(dpy,drawwindow,Penc[0],bx2,by2,bx1,by2) ;
      XDrawLine(dpy,drawwindow,Penc[0],bx1,by2,bx1,by1) ;
    }
  }
  void Xclearbox() {
    if(drawnbox) {
      int x1 = std::min(bx1,bx2) ;
      int x2 = std::max(bx1,bx2) ;
      int y1 = std::min(by1,by2) ;
      int y2 = std::max(by1,by2) ;
      XCopyArea(dpy,drawpixmap,drawwindow,linec,x1,y1,1,y2-y1+1,x1,y1);
      XCopyArea(dpy,drawpixmap,drawwindow,linec,x1,y2,x2-x1+1,1,x1,y2);
      XCopyArea(dpy,drawpixmap,drawwindow,linec,x2,y1,1,y2-y1+1,x2,y1);
      XCopyArea(dpy,drawpixmap,drawwindow,linec,x1,y1,x2-x1+1,1,x1,y1);
    }
  }

  void clear()
    { XFillRectangle(dpy,drawpixmap,forec,0,0,drawwidth,drawheight) ; }
  void refresh() {
    XCopyArea(dpy,drawpixmap,drawwindow,linec,0,0,drawwidth,drawheight,0,0);
    Xdrawbox() ;
  }
  void refresh(int x, int y, int width, int height) {
    XCopyArea(dpy,drawpixmap,drawwindow,linec,x,y,width,height,x,y);
    Xdrawbox() ;
  }
  void resize(int width,int height) ;
    
public:

  XDevice(Widget drawwidget) ;
  
  virtual ~XDevice() { XCloseDisplay(dpy) ; }
  virtual double getTickSize() { return ticksize ; }
  virtual bool isLabelDev() { return true ; }
  virtual void OpenDevice();
  virtual void CloseDevice() ;
  virtual void HardCloseDevice() ;
  virtual void SelectPen(int Pen) ;
  virtual void DrawLine(const positions &p1,const positions &p2) ;
  virtual void DrawLabel(const std::string &text, positions location,
                         int pos, int size, double degrees) ;
  virtual void FillShade(int shade) ;
  virtual void FillPoly(const positions *pnts, int num_pnts) ;
  virtual void FillPolyMarker(const positions *pnts, int num_pnts) ;
} ;

void xinteract(int argc, char *argv[]) ;

#endif
