#include "PSdev.h"
#include "Xdev.h"
#include <stdlib.h>
using std::max ;
using std::min ;
using std::endl ;

namespace {
#ifdef NOT_IMPL
  typedef int RGBENTRY[3] ;
  RGBENTRY rgb_table[] = {
    0x00, 0x00, 0x00, // Black
    0x00, 0x00, 0xff, // blue
    0x00, 0x22, 0xdd,
    0x00, 0x44, 0xbb,
    0x00, 0x66, 0x99,
    0x00, 0x88, 0x77,
    0x00, 0xaa, 0x55,
    0x00, 0xcc, 0x22,
    0x00, 0xff, 0x00, // Green Pen 9
    0x44, 0xff, 0x00,
    0x66, 0xdd, 0x00,
    0x77, 0xcc, 0x00,
    0x88, 0x88, 0x00,
    0xaa, 0x66, 0x00,
    0xcc, 0x44, 0x00,
    0xdd, 0x22, 0x00,
    0xff, 0x00, 0x00,
    0xff, 0xff, 0xff
  } ;
#endif
  double font_size_factor = 1.0 ;
}
PSDevCore::PSDevCore() {
  r[0] = 0 ;
  g[0] = 0 ;
  b[0] = 0 ;
  for(int i=0;i<MAXPENS-1;++i) {
    r[i+1] = double(X_R[i])/255.0 ;
    g[i+1] = double(X_G[i])/255.0 ;
    b[i+1] = double(X_B[i])/255.0 ;
  }
  r[MAXPENS] = 1.0 ;
  g[MAXPENS] = 1.0 ;
  b[MAXPENS] = 1.0 ;
}

void PSDevCore::def_funcs() {
  if(outf==NULL)
    return ;
  fprintf(outf,"/L { newpath moveto lineto stroke } def\n") ;
  fprintf(outf,"/P { newpath moveto {count 0 eq {exit} if lineto} loop fill } def\n") ;
  fprintf(outf,"/W {setlinewidth} def\n") ;
  fprintf(outf,"/S {setgray} def\n") ;
  fprintf(outf,"/C {setrgbcolor} def\n") ;
  fprintf(outf,"\n") ;
  fprintf(outf,"/T\n") ;
  fprintf(outf,"{\t/Courier-Bold findfont exch scalefont setfont gsave \n") ;
  fprintf(outf,"\ttranslate rotate translate moveto true charpath fill stroke grestore\n") ;
  fprintf(outf,"} def\n") ;
  fprintf(outf,"0 W\n") ;
  fprintf(outf,"\n") ;
  currentPen = -1 ;
}
  
double PSDevCore::getTickSize() {
  return 10.0 ;
}

void PSDevCore::SelectPen(int Pen)
{
  Pen = Pen % MAXPENS ;
  if(currentPen != Pen) {
    currentPen = Pen ;
    if(outf)
      fprintf(outf,"%f %f %f C\n",r[currentPen],g[currentPen],b[currentPen]) ;
  }
}
    
void PSDevCore::DrawLine(const positions &p1, const positions &p2)
{
  if(outf)
    fprintf(outf,"%f %f %f %f L\n",p1.x,p1.y,p2.x,p2.y) ;
}

void PSDevCore::DrawLabel(const std::string &text, positions location,
                          int pos, int size, double degrees)
{
  
  int len = text.size();
  if (len == 0)
    return;
  double sh =  size *font_size_factor;
  double sw = .7 * sh ;

  double x0=0,y0=0, x1=0,y1=0, xs,ys, xt,yt ;
  double fx = location.x ;
  double fy = location.y ;
  
  switch( pos )
    {   case 0:     /*  center at cursor */
      x0 = fx - sw * ( len / 2 );
        x1 = fx + sw * ( len / 2 + 1 );
        y0 = fy;
        y1 = fy;
        break ;
    case 1:     /*  center above cursor */
      x0 = fx - sw * len / 2;
      x1 = fx + sw * ( len / 2 + 1 );
      y0 = fy + sh;
      y1 = fy + sh;
      break ;
    case 2:     /*  left justify above cursor */
      x0 = fx;
      x1 = fx + sw * ( len + 1 );
      y0 = fy + sh;
      y1 = fy + sh;
      break ;
    case 3:     /*  left justify to the right of cursor */
      x0 = fx + sw;
      x1 = fx + sw * ( len + 1 );
      y0 = fy;
      y1 = fy;
      break ;
    case 4:     /*  left justify below cursor */
      x0 = fx + sw;
      x1 = fx + sw * ( len + 1 );
      y0 = fy - sh;
      y1 = fy - sh;
      break ;
    case 5:     /*  center below cursor */
      x0 = fx - sw * len / 2;
      x1 = fx + sw * ( len / 2 + 1 );
      y0 = fy - sh;
      y1 = fy - sh;
      break;
    case 6:     /*  right justify below cursor */
      x0 = fx - sw * len;
      x1 = fx + sw;
      y0 = fy - sh;
      y1 = fy - sh;
      break;
    case 7:     /*  right justify to the left of cursor */
      x0 = fx - sw * len;
      x1 = fx + sw;
      y0 = fy;
      y1 = fy;
      break ;
    case 8:     /*  right justify above cursor */
      x0 = fx - sw * len;
      x1 = fx + sw;
      y0 = fy + sh;
      y1 = fy + sh;
      break ;
    case 9:     /*  center two lines below cursor */
      x0 = fx - sw * len;
      x1 = fx + sw * ( len + 1 );
      y0 = fy - 2 * sh;
      y1 = fy - 2 * sh;
      break;
    case 10:    /*  center two lines above cursor */
      x0 = fx - sw * len;
      x1 = fx + sw * ( len + 1 );
      y0 = fy + 2 * sh;
      y1 = fy + 2 * sh;
      break;
    default :
      break ;
    }
  /***************************************
   * now rotate (x0,y0) and (x1,y1)       *
   * about (x,y) to produce new (x0,y0)   *
   * and (x1,y1).                         *
   ***************************************/
  double theta = degrees * 3.14159 / 180;

  x0 -= fx;
  y0 -= fy;
  x1 -= fx;
  y1 -= fy;
    
  xt = x0;
  yt = y0;

  x0 = xt * cos(theta) + yt * sin( theta);
  y0 = yt * cos(theta) - xt * sin( theta);
        
  xt = x1;
  yt = y1;
  x1 = xt * cos(theta) + yt * sin( theta);
  y1 = yt * cos(theta) - xt * sin( theta);
    
  x0 += fx;
  y0 += fy;
  x1 += fx;
  y1 += fy; 
  for(int n = 0; n < len; ++n) {
    xs = x0 + ( ( x1 - x0 ) / len ) * n;
    ys = y0 + ( ( y1 - y0 ) / len ) * n;
    if(outf)
      fprintf(outf,"(%c) %f %f %f %f %f %f %f %d T\n",text[n],
              xs,ys,-xs,-ys,degrees,xs,ys,int(sh)) ;
  }
}

void PSDevCore::FillShade(int Pen)
{
  if(currentPen != Pen) {
    currentPen = Pen ;

    if(outf)
      fprintf(outf,"%f %f %f C\n",r[currentPen],g[currentPen],b[currentPen]) ;
  }
}

void PSDevCore::FillPoly(const positions *pnts, int num_pnts)
{
  if(num_pnts < 2)
    return ;
  for(int i=0;i<num_pnts;++i)
    if(outf)
      fprintf(outf,"%f %f ",pnts[i].x, pnts[i].y) ;
  if(outf)
    fprintf(outf,"P\n") ;
}

void PSDevCore::FillPolyMarker(const positions *pnts, int num_pnts)
{
  if(num_pnts < 2)
    return ;
  // Fixed marker color (currently black), set directly rather than
  // through SelectPen/FillShade's cached pen index -- currentPen is
  // invalidated afterward so the next normal draw call still emits its
  // own color correctly.
  if(outf)
    fprintf(outf,"0 0 0 C\n") ;
  for(int i=0;i<num_pnts;++i)
    if(outf)
      fprintf(outf,"%f %f ",pnts[i].x, pnts[i].y) ;
  if(outf)
    fprintf(outf,"P\n") ;
  currentPen = -1 ;
}

PSDevEPS::PSDevEPS(const DevBox &b, const std::string &filename) {
  file = filename ;
  DevBox b2 ;
  b2.xleft = min(b.xleft,b.xright) ;
  b2.xright = max(b.xleft,b.xright) ;
  b2.ybottom = min(b.ybottom,b.ytop) ;
  b2.ytop = max(b.ybottom,b.ytop) ;
  setdrawbox(b2) ;
}

bool PSDevEPS::isLabelDev() {
  return false ;
}

void PSDevEPS::OpenDevice() {
  if((outf = fopen(file.c_str(),"w")) == NULL) {
    std::cerr << "can't open " << file << endl ;
    return ;
  }
  fprintf(outf,"%%!PS-Adobe-2.0 EPSF-1.2\n") ;
  DevBox b = getDrawBox() ;
  fprintf(outf,"%%%%BoundingBox: %f %f %f %f\n",
          b.xleft,b.ybottom,b.xright,b.ytop) ;
  fprintf(outf, "%% Generated from 2dgv\n") ;
  def_funcs() ;
  double bxmin = b.xleft ;
  double bymin = b.ybottom ;
  double bxmax = b.xright ;
  double bymax = b.ytop ;
  double maxscale =  max(fabs(bxmin-bxmax),fabs(bymin-bymax)) ;
  font_size_factor = maxscale/1000.0 ;
  fprintf(outf,
          "newpath %f %f moveto "
          "%f %f lineto %f %f lineto %f %f lineto %f %f lineto "
          "closepath clip\n",
          bxmin, bymin, bxmax, bymin, bxmax, bymax,
          bxmin, bymax, bxmin, bymin) ;
}

void PSDevEPS::CloseDevice() {
  if(outf)
    fclose(outf) ;
}

void PSDevEPS::HardCloseDevice() {
  if(outf) {
    fprintf(outf,"%% Improperly terminated\n") ;
    fclose(outf) ;
  }
}

PSDevPSFile::PSDevPSFile(const std::string &filename) {
  file = filename ;
  DevBox b2 ;
  b2.xleft = 340 ;
  b2.xright = 2375 ;
  b2.ybottom = 175 ;
  b2.ytop = 3170 ;
  font_size_factor = 2000.0/1000.0 ;
  setdrawbox(b2) ;
}

bool PSDevPSFile::isLabelDev() {
  return true ;
}

void PSDevPSFile::OpenDevice() {
  if((outf = fopen(file.c_str(),"w")) == NULL) {
    std::cerr << "can't open " << file << endl ;
    return ;
  }
  fprintf(outf,"%%!PS 2dgv 2.0 Postscript Output File\n")  ;
  fprintf(outf,"initgraphics grestoreall 72.0 300.0 div 72.0 300.0 div scale\n") ;
  DevBox b = getDrawBox() ;
  def_funcs() ;
  double bxmin = b.xleft ;
  double bymin = b.ybottom ;
  double bxmax = b.xright ;
  double bymax = b.ytop ;
  fprintf(outf,
          "newpath %f %f moveto "
          "%f %f lineto %f %f lineto %f %f lineto %f %f lineto "
          "closepath clip\n",
          bxmin, bymin, bxmax, bymin, bxmax, bymax,
          bxmin, bymax, bxmin, bymin) ;
}

void PSDevPSFile::CloseDevice() {
  if(outf) {
    fprintf(outf,"showpage\n") ;
    fclose(outf) ;
  }
}

void PSDevPSFile::HardCloseDevice() {
  if(outf) {
    fprintf(outf,"%% Improperly terminated\n") ;
    fclose(outf) ;
  }
}

PSDevPipe::PSDevPipe(const std::string &pcom) {
  pipe_command = pcom ;
  DevBox b2 ;
  b2.xleft = 340 ;
  b2.xright = 2375 ;
  b2.ybottom = 175 ;
  b2.ytop = 3170 ;
  font_size_factor = 2000.0/1000.0 ;
  setdrawbox(b2) ;
}

bool PSDevPipe::isLabelDev() {
  return true ;
}

void PSDevPipe::OpenDevice() {
  //  if((outf = fopen(file.c_str(),"w")) == NULL) {
  //    std::cerr << "can't open " << file << endl ;
  //    return ;
  //  }
  if((outf = popen(pipe_command.c_str(),"w")) == NULL) {
    std::cerr << "popen(" << pipe_command << ") failed" << endl ;
    perror("2dgv!:") ;
    if((outf = fopen("/tmp/crash.file","w")) == NULL) {
      perror("ARGH! /tmp/crash.file") ;
      exit(-1) ;
    }
  }
  fprintf(outf,"%%!PS 2dgv 2.0 Postscript Output File\n")  ;
  fprintf(outf,"initgraphics grestoreall 72.0 300.0 div 72.0 300.0 div scale\n") ;
  DevBox b = getDrawBox() ;
  def_funcs() ;
  double bxmin = b.xleft ;
  double bymin = b.ybottom ;
  double bxmax = b.xright ;
  double bymax = b.ytop ;
  fprintf(outf,
          "newpath %f %f moveto "
          "%f %f lineto %f %f lineto %f %f lineto %f %f lineto "
          "closepath clip\n",
          bxmin, bymin, bxmax, bymin, bxmax, bymax,
          bxmin, bymax, bxmin, bymin) ;
}

void PSDevPipe::CloseDevice() {
  if(outf) {
    fprintf(outf,"showpage\n") ;
    fclose(outf) ;
  }
}

void PSDevPipe::HardCloseDevice() {
  if(outf) {
    fprintf(outf,"%% Improperly terminated\n") ;
    fclose(outf) ;
  }
}
