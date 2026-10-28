#include "Xdev.h"
#include "PSdev.h"
#include "grid.h"
#include "draw.h"
#include "debug.h"

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"

#include <iostream>
#include <fstream>

using std::ifstream ;
using std::ios ;
using std::cerr ;
using std::cout ;
using std::cin ;
using std::endl ;

#include <utility>
using std::min ;
using std::max ;

static const int XDRAWSIZE = 600 ;

static XDevice *Xdev = 0 ;

static Widget toplevel ;



static char mdebugString[] = "-debug" ;
static char sdebugString[] = "*debug" ;

static XrmOptionDescRec opt_table[] = {
  {mdebugString, sdebugString,	XrmoptionNoArg,	NULL},
};

namespace STR{
  char defaultfore[] = "defaultfore" ;
  char DefaultFore[] = "DefaultFore" ;
  char defaultback[] = "defaultback" ;
  char DefaultBack[] = "DefaultBack" ;
  char black[] = "black" ;
  char gray95[] = "gray95" ;
  char pen1[] = "pen1" ;
  char Pen1[] = "Pen1" ;

  char xpen[64][7] =
    {
      "xpen02" ,
      "xpen03" ,
      "xpen04" ,
      "xpen05" ,
      "xpen06" ,
      "xpen07" ,
      "xpen08" ,
      "xpen09" ,
      "xpen10" ,
      "xpen11" ,
      "xpen12" ,
      "xpen13" ,
      "xpen14" ,
      "xpen15" ,
      "xpen16" ,
      "xpen17" ,
      "xpen18" ,
      "xpen19" ,
      "xpen20" ,
      "xpen21" ,
      "xpen22" ,
      "xpen23" ,
      "xpen24" ,
      "xpen25" ,
      "xpen26" ,
      "xpen27" ,
      "xpen28" ,
      "xpen29" ,
      "xpen30" ,
      "xpen31" ,
      "xpen32" ,
      "xpen33" ,
      "xpen02" ,
      "xpen03" ,
      "xpen04" ,
      "xpen05" ,
      "xpen06" ,
      "xpen07" ,
      "xpen08" ,
      "xpen09" ,
      "xpen10" ,
      "xpen11" ,
      "xpen12" ,
      "xpen13" ,
      "xpen14" ,
      "xpen15" ,
      "xpen16" ,
      "xpen17" ,
      "xpen18" ,
      "xpen19" ,
      "xpen20" ,
      "xpen21" ,
      "xpen22" ,
      "xpen23" ,
      "xpen24" ,
      "xpen25" ,
      "xpen26" ,
      "xpen27" ,
      "xpen28" ,
      "xpen29" ,
      "xpen30" ,
      "xpen31" ,
      "xpen32" ,
      "xpen33" ,
    } ;
  char col[64][8] =
    {
      "#0000FF",
      "#0022DD" ,
      "#0044BB" ,
      "#006699" ,
      "#008877" ,
      "#00AA55" ,
      "#00CC22" ,
      "#00FF00" ,
      "#44FF00" ,
      "#66DD00" ,
      "#77CC00" ,
      "#888800" ,
      "#AA6600" ,
      "#CC4400" ,
      "#DD2200" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#0000FF" ,
      "#0022DD" ,
      "#0044BB" ,
      "#006699" ,
      "#008877" ,
      "#00AA55" ,
      "#00CC22" ,
      "#00FF00" ,
      "#44FF00" ,
      "#66DD00" ,
      "#77CC00" ,
      "#888800" ,
      "#AA6600" ,
      "#CC4400" ,
      "#DD2200" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
      "#FF0000" ,
    } ;
}

int X_R[256], X_G[256], X_B[256] ;

struct AppResourceRec {
  Cardinal suggest_width, suggest_height ;
  Pixel defaultfore;	/* Default foreground */
  Pixel defaultback;	/* Default background */
  Pixel Pens[MAXPENS] ;
}  ;

static AppResourceRec appRes ;

static XtResource appResourceSpec[] = {
  { XtNwidth, XtCWidth, XtRInt, sizeof(Cardinal),
    XtOffsetOf(AppResourceRec,suggest_width),XtRImmediate,
    (XtPointer) XDRAWSIZE},
  { XtNheight, XtCHeight, XtRInt, sizeof(Cardinal),
    XtOffsetOf(AppResourceRec,suggest_height),XtRImmediate,
    (XtPointer) XDRAWSIZE},
  { &STR::defaultfore[0], &STR::DefaultFore[0], XtRPixel, sizeof(Pixel),
    XtOffsetOf(AppResourceRec,defaultfore), XtRString, &STR::black[0]},
  { &STR::defaultback[0], &STR::DefaultBack[0], XtRPixel, sizeof(Pixel),
    XtOffsetOf(AppResourceRec,defaultback), XtRString, &STR::gray95[0]},
  {&STR::pen1[0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[0]), XtRString, &STR::black[0]},
  {&STR::xpen[0][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[1]), XtRString, &STR::col[0][0]},
  {&STR::xpen[1][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[2]), XtRString, &STR::col[1][0]},
  {&STR::xpen[2][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[3]), XtRString, &STR::col[2][0]},
  {&STR::xpen[3][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[4]), XtRString, &STR::col[3][0]},
  {&STR::xpen[4][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[5]), XtRString, &STR::col[4][0]},
  {&STR::xpen[5][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[6]), XtRString, &STR::col[5][0]},
  {&STR::xpen[6][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[7]), XtRString, &STR::col[6][0]},
  {&STR::xpen[7][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[8]), XtRString, &STR::col[7][0]},
  {&STR::xpen[8][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[9]), XtRString, &STR::col[8][0]},
  {&STR::xpen[9][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[10]), XtRString, &STR::col[9][0]},
  {&STR::xpen[10][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[11]), XtRString, &STR::col[10][0]},
  {&STR::xpen[11][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[12]), XtRString, &STR::col[11][0]},
  {&STR::xpen[12][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[13]), XtRString, &STR::col[12][0]},
  {&STR::xpen[13][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[14]), XtRString, &STR::col[13][0]},
  {&STR::xpen[14][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[15]), XtRString, &STR::col[14][0]},
  {&STR::xpen[15][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[16]), XtRString, &STR::col[15][0]},
  {&STR::xpen[16][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[17]), XtRString, &STR::col[16][0]},
  {&STR::xpen[17][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[18]), XtRString, &STR::col[17][0]},
  {&STR::xpen[18][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[19]), XtRString, &STR::col[18][0]},
  {&STR::xpen[19][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[20]), XtRString, &STR::col[19][0]},
  {&STR::xpen[20][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[21]), XtRString, &STR::col[20][0]},
  {&STR::xpen[21][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[22]), XtRString, &STR::col[21][0]},
  {&STR::xpen[22][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[23]), XtRString, &STR::col[22][0]},
  {&STR::xpen[23][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[24]), XtRString, &STR::col[23][0]},
  {&STR::xpen[24][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[25]), XtRString, &STR::col[24][0]},
  {&STR::xpen[25][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[26]), XtRString, &STR::col[25][0]},
  {&STR::xpen[26][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[27]), XtRString, &STR::col[26][0]},
  {&STR::xpen[27][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[28]), XtRString, &STR::col[27][0]},
  {&STR::xpen[28][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[29]), XtRString, &STR::col[28][0]},
  {&STR::xpen[29][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[30]), XtRString, &STR::col[29][0]},
  {&STR::xpen[30][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[31]), XtRString, &STR::col[30][0]},
  {&STR::xpen[31][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[32]), XtRString, &STR::col[31][0]},
  {&STR::xpen[32][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[33]), XtRString, &STR::col[32][0]},
  {&STR::xpen[33][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[34]), XtRString, &STR::col[33][0]},
  {&STR::xpen[34][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[35]), XtRString, &STR::col[34][0]},
  {&STR::xpen[35][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[36]), XtRString, &STR::col[35][0]},
  {&STR::xpen[36][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[37]), XtRString, &STR::col[36][0]},
  {&STR::xpen[37][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[38]), XtRString, &STR::col[37][0]},
  {&STR::xpen[38][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[39]), XtRString, &STR::col[38][0]},
  {&STR::xpen[39][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[40]), XtRString, &STR::col[39][0]},
  {&STR::xpen[40][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[41]), XtRString, &STR::col[40][0]},
  {&STR::xpen[41][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[42]), XtRString, &STR::col[41][0]},
  {&STR::xpen[42][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[43]), XtRString, &STR::col[42][0]},
  {&STR::xpen[43][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[44]), XtRString, &STR::col[43][0]},
  {&STR::xpen[44][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[45]), XtRString, &STR::col[44][0]},
  {&STR::xpen[45][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[46]), XtRString, &STR::col[45][0]},
  {&STR::xpen[46][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[47]), XtRString, &STR::col[46][0]},
  {&STR::xpen[47][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[48]), XtRString, &STR::col[47][0]},
  {&STR::xpen[48][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[49]), XtRString, &STR::col[48][0]},
  {&STR::xpen[49][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[50]), XtRString, &STR::col[49][0]},
  {&STR::xpen[50][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[51]), XtRString, &STR::col[50][0]},
  {&STR::xpen[51][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[52]), XtRString, &STR::col[51][0]},
  {&STR::xpen[52][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[53]), XtRString, &STR::col[52][0]},
  {&STR::xpen[53][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[54]), XtRString, &STR::col[53][0]},
  {&STR::xpen[54][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[55]), XtRString, &STR::col[54][0]},
  {&STR::xpen[55][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[56]), XtRString, &STR::col[55][0]},
  {&STR::xpen[56][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[57]), XtRString, &STR::col[56][0]},
  {&STR::xpen[57][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[58]), XtRString, &STR::col[57][0]},
  {&STR::xpen[58][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[59]), XtRString, &STR::col[58][0]},
  {&STR::xpen[59][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[60]), XtRString, &STR::col[59][0]},
  {&STR::xpen[60][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[61]), XtRString, &STR::col[60][0]},
  {&STR::xpen[61][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[62]), XtRString, &STR::col[61][0]},
  {&STR::xpen[62][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[63]), XtRString, &STR::col[62][0]},
  {&STR::xpen[63][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
   XtOffsetOf(AppResourceRec,Pens[64]), XtRString, &STR::col[63][0]},

} ;



/*
  static XtResource resources[] = {
  {XtNwidth, XtCWidth, XtRInt, sizeof(int),
  suggest_width, XtRImmediate, (caddr_t) XDRAWSIZE},
  {XtNheight, XtCHeight, XtRInt, sizeof(int),
  suggest_height, XtRImmediate, (caddr_t) XDRAWSIZE},
  {&STR::defaultfore[0], &STR::DefaultFore[0], XtRPixel, sizeof(Pixel),
  defaultfore, XtRString, &STR::black[0]},
  {&STR::defaultback[0], &STR::DefaultBack[0], XtRPixel, sizeof(Pixel),
  defaultback, XtRString, &STR::gray95[0]},
  {&STR::pen1[0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[0], XtRString, &STR::black[0]},
  {&STR::xpen[0][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[1], XtRString, &STR::col[0][0]},
  {&STR::xpen[1][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[2], XtRString, &STR::col[1][0]},
  {&STR::xpen[2][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[3], XtRString, &STR::col[2][0]},
  {&STR::xpen[3][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[4], XtRString, &STR::col[3][0]},
  {&STR::xpen[4][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[5], XtRString, &STR::col[4][0]},
  {&STR::xpen[5][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[6], XtRString, &STR::col[5][0]},
  {&STR::xpen[6][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[7], XtRString, &STR::col[6][0]},
  {&STR::xpen[7][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[8], XtRString, &STR::col[7][0]},
  {&STR::xpen[8][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[9], XtRString, &STR::col[8][0]},
  {&STR::xpen[9][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[10], XtRString, &STR::col[9][0]},
  {&STR::xpen[10][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[11], XtRString, &STR::col[10][0]},
  {&STR::xpen[11][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[12], XtRString, &STR::col[11][0]},
  {&STR::xpen[12][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[13], XtRString, &STR::col[12][0]},
  {&STR::xpen[13][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[14], XtRString, &STR::col[13][0]},
  {&STR::xpen[14][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[15], XtRString, &STR::col[14][0]},
  {&STR::xpen[15][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[16], XtRString, &STR::col[15][0]},
  {&STR::xpen[16][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[17], XtRString, &STR::col[16][0]},
  {&STR::xpen[17][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[18], XtRString, &STR::col[17][0]},
  {&STR::xpen[18][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[19], XtRString, &STR::col[18][0]},
  {&STR::xpen[19][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[20], XtRString, &STR::col[19][0]},
  {&STR::xpen[20][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[21], XtRString, &STR::col[20][0]},
  {&STR::xpen[21][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[22], XtRString, &STR::col[21][0]},
  {&STR::xpen[22][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[23], XtRString, &STR::col[22][0]},
  {&STR::xpen[23][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[24], XtRString, &STR::col[23][0]},
  {&STR::xpen[24][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[25], XtRString, &STR::col[24][0]},
  {&STR::xpen[25][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[26], XtRString, &STR::col[25][0]},
  {&STR::xpen[26][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[27], XtRString, &STR::col[26][0]},
  {&STR::xpen[27][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[28], XtRString, &STR::col[27][0]},
  {&STR::xpen[28][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[29], XtRString, &STR::col[28][0]},
  {&STR::xpen[29][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[30], XtRString, &STR::col[29][0]},
  {&STR::xpen[30][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[31], XtRString, &STR::col[30][0]},
  {&STR::xpen[31][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[32], XtRString, &STR::col[31][0]},
  {&STR::xpen[32][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[33], XtRString, &STR::col[32][0]},
  {&STR::xpen[33][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[34], XtRString, &STR::col[33][0]},
  {&STR::xpen[34][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[35], XtRString, &STR::col[34][0]},
  {&STR::xpen[35][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[36], XtRString, &STR::col[35][0]},
  {&STR::xpen[36][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[37], XtRString, &STR::col[36][0]},
  {&STR::xpen[37][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[38], XtRString, &STR::col[37][0]},
  {&STR::xpen[38][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[39], XtRString, &STR::col[38][0]},
  {&STR::xpen[39][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[40], XtRString, &STR::col[39][0]},
  {&STR::xpen[40][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[41], XtRString, &STR::col[40][0]},
  {&STR::xpen[41][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[42], XtRString, &STR::col[41][0]},
  {&STR::xpen[42][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[43], XtRString, &STR::col[42][0]},
  {&STR::xpen[43][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[44], XtRString, &STR::col[43][0]},
  {&STR::xpen[44][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[45], XtRString, &STR::col[44][0]},
  {&STR::xpen[45][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[46], XtRString, &STR::col[45][0]},
  {&STR::xpen[46][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[47], XtRString, &STR::col[46][0]},
  {&STR::xpen[47][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[48], XtRString, &STR::col[47][0]},
  {&STR::xpen[48][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[49], XtRString, &STR::col[48][0]},
  {&STR::xpen[49][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[50], XtRString, &STR::col[49][0]},
  {&STR::xpen[50][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[51], XtRString, &STR::col[50][0]},
  {&STR::xpen[51][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[52], XtRString, &STR::col[51][0]},
  {&STR::xpen[52][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[53], XtRString, &STR::col[52][0]},
  {&STR::xpen[53][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[54], XtRString, &STR::col[53][0]},
  {&STR::xpen[54][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[55], XtRString, &STR::col[54][0]},
  {&STR::xpen[55][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[56], XtRString, &STR::col[55][0]},
  {&STR::xpen[56][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[57], XtRString, &STR::col[56][0]},
  {&STR::xpen[57][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[58], XtRString, &STR::col[57][0]},
  {&STR::xpen[58][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[59], XtRString, &STR::col[58][0]},
  {&STR::xpen[59][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[60], XtRString, &STR::col[59][0]},
  {&STR::xpen[60][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[61], XtRString, &STR::col[60][0]},
  {&STR::xpen[61][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[62], XtRString, &STR::col[61][0]},
  {&STR::xpen[62][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[63], XtRString, &STR::col[62][0]},
  {&STR::xpen[63][0], &STR::Pen1[0], XtRPixel, sizeof(Pixel),
  Pens[64], XtRString, &STR::col[63][0]},
  };
*/
const int  NUMBUTTONARGS = 1 ;

void SetWMGetInput(Widget w) {
  XWMHints wmhints;
  wmhints.flags = InputHint | StateHint ;
  wmhints.input = True;
  wmhints.initial_state = NormalState;
  Window win = XtWindowOfObject(w) ;
  Display *dpy = XtDisplay(w) ;
  XSetWMHints(dpy, win, &wmhints);
}

Widget
Button(const char *name,Widget parent,XtCallbackProc action)
{
  XtCallbackRec *cbp = (XtCallbackRec *) malloc(sizeof(XtCallbackRec)*3);
  Arg           *alp = (Arg *)      malloc(sizeof(Arg)*NUMBUTTONARGS) ;

  cbp[0].callback = action ;
  cbp[0].closure = NULL ;
  cbp[1].callback = NULL ;
  cbp[2].closure = NULL ;
  XtSetArg(alp[0],XtNcallback,(XtArgVal) cbp) ;
  return( XtCreateManagedWidget(name,commandWidgetClass,parent,alp,NUMBUTTONARGS)) ;
}

void XDevice::resize(int width,int height) {
  if(width != drawwidth || height != drawheight) {
    drawwidth = width ;
    drawheight = height ;
    drawbox.xleft = 0 ;
    drawbox.xright = drawwidth ;
    drawbox.ytop = 0 ;
    drawbox.ybottom = drawheight ;

    if(pixmap_allocated)
      XFreePixmap(dpy,drawpixmap) ;
    drawpixmap = XCreatePixmap(dpy,drawwindow,drawwidth,drawheight, drawdepth) ;
    clear() ;
    if(grids.size() != 0)
      draw(grids.front(),*Xdev) ;
    refresh() ;
  }
}


static void ClassInitialize() {}

static void Initialize(Widget request, Widget w, ArgList args, Cardinal *cp)
{

}


void ResizeEvent(Widget w,XtPointer data,XEvent *event,
                 Boolean *continue_dispatch)
{
  FATAL(Xdev == 0) ;
  Xdev->resize(w->core.width,w->core.height) ;
}

// Forward declarations: both are defined later in this file but are
// called from mouseButtonPress below.
void Zoomin(Widget w, XtPointer p1, XtPointer p2) ;
void ProbeNode(int screenX, int screenY) ;

void mouseButtonPress(Widget w,XtPointer *data, XEvent *event,
                      Boolean *continue_dispatch)
{
  switch(event->type) {
  case ButtonPress:
    if(event->xbutton.button == Button2) {
      // Middle click is a shortcut for pressing the Zoom button: zooms
      // into the drawn box if one exists, or steps back one zoom level
      // if not (Zoomin's own else-branch). It doesn't start a box drag
      // of its own.
      Zoomin(w,NULL,NULL) ;
      break ;
    }
    if(event->xbutton.button == Button3) {
      // Right click probes the nearest mesh node to the click. It
      // doesn't start a box drag of its own.
      ProbeNode(event->xbutton.x,event->xbutton.y) ;
      break ;
    }
    if(event->xbutton.button != Button1)
      break ;
    Xdev->Xclearbox() ;
    Xdev->bx1 = event->xbutton.x ;
    Xdev->by1 = event->xbutton.y ;
    Xdev->drawnbox = false ;
    break ;
  case ButtonRelease:
    // Only left button's release completes a box drag; middle/right
    // click are handled entirely on press above and shouldn't affect
    // box state on release.
    if(event->xbutton.button != Button1)
      break ;
    Xdev->Xclearbox() ;
    Xdev->bx2 = event->xbutton.x ;
    Xdev->by2 = event->xbutton.y ;
    if(Xdev->bx1 == Xdev->bx2 && Xdev->by1 == Xdev->by2)
      Xdev->drawnbox = false ;
    else
      Xdev->drawnbox = true ;
    Xdev->Xdrawbox() ;
    break ;
  case MotionNotify:
    // Previously folded into the same fallthrough as ButtonRelease,
    // which updated the box preview on any motion regardless of which
    // button was actually held -- so dragging with the middle or right
    // button (even though pressing them doesn't start a box) could
    // still leave a stale box behind via motion events alone. This now
    // only updates the box preview while the left button specifically
    // is held.
    if(!(event->xmotion.state & Button1Mask))
      break ;
    Xdev->Xclearbox() ;
    Xdev->bx2 = event->xmotion.x ;
    Xdev->by2 = event->xmotion.y ;
    if(Xdev->bx1 == Xdev->bx2 && Xdev->by1 == Xdev->by2)
      Xdev->drawnbox = false ;
    else
      Xdev->drawnbox = true ;
    Xdev->Xdrawbox() ;
    break ;
  default:
    break ;
  }
  *continue_dispatch = True ;
  return ;
}

void XDevice::realize(Widget w) {
  XGCValues values ;
  values.foreground = appRes.defaultfore ;
  values.background = w->core.background_pixel ;

  linec = XtGetGC(w, (XtGCMask) GCForeground | GCBackground, &values) ;
  values.foreground =  appRes.defaultback ;
  forec = XtGetGC(w, (XtGCMask) GCForeground | GCBackground, &values) ;
  XFontStruct *font_info = XLoadQueryFont(dpy,"9x15") ;
  XSetFont(dpy,linec,font_info->fid) ;
  for(int i=0; i < MAXPENS; i++) {
    values.foreground = appRes.Pens[i] ;
    Penc[i] = XtGetGC(w, (XtGCMask) GCForeground | GCBackground, &values) ;
    XSetFont(dpy,Penc[i],font_info->fid) ;
  }

  // Dedicated GC for the probe node marker (see FillPolyMarker below):
  // the Pens[] table above is a fixed 64-entry contour colormap
  // gradient, with no way to get a specific marker color out of it, so
  // this is allocated directly instead.
  XColor markerColor ;
  Colormap cmap = DefaultColormap(dpy,DefaultScreen(dpy)) ;
  XParseColor(dpy,cmap,"black",&markerColor) ;
  XAllocColor(dpy,cmap,&markerColor) ;
  values.foreground = markerColor.pixel ;
  markerGC = XtGetGC(w, (XtGCMask) GCForeground | GCBackground, &values) ;

  drawwidth = w->core.width ;
  drawheight = w->core.height ;
  drawwindow = w->core.window ;
  drawdepth = w->core.depth ;

  drawbox.xleft = 0 ;
  drawbox.xright = drawwidth ;
  drawbox.ytop = 0 ;
  drawbox.ybottom = drawheight ;

  drawpixmap = XCreatePixmap(dpy, drawwindow,drawwidth,drawheight,w->core.depth) ;
  pixmap_allocated = true ; ;
  clear() ;
  if(grids.size() != 0)
    draw(grids.front(),*Xdev) ;
}

static void Realize(Widget w, Mask *valueMask,
                    XSetWindowAttributes *attributes)
{
  attributes->save_under = False;
  attributes->backing_store = NotUseful ;

  XtCreateWindow(w, (unsigned int) InputOutput, (Visual *) CopyFromParent,
                 *valueMask | CWSaveUnder,
                 attributes) ;

  FATAL(Xdev == 0) ;
  Xdev->realize(w) ;
}


static void HandleExpose(Widget w, XEvent *event, Region r)
{
  int x,y,width,height;
  if(event->xexpose.count) return ;
  x = event->xexpose.x ;
  y = event->xexpose.y ;
  width = event->xexpose.width ;
  height = event->xexpose.height ;
  FATAL(Xdev == 0) ;
  Xdev->refresh(x,y,width,height) ;
}

static Widget lpop,linput  ;
static Widget pop ;

void popquit(Widget w=0,XtPointer p1=0,XtPointer p2=0)
{
  XtRemoveGrab(pop) ;
  XtDestroyWidget(pop) ;
}


void lpopquit(Widget w=0,XtPointer p1=0,XtPointer p2=0)
{
  XtRemoveGrab(lpop) ;
  XtDestroyWidget(lpop) ;

}



void Quit(Widget w, XtPointer p1, XtPointer p2)
{
  XtDestroyWidget(toplevel) ;
  delete Xdev ;
  printf("Bye\n") ;
  exit(0) ;
}

void Zoomin(Widget w, XtPointer p1, XtPointer p2)
{
  if(grids.size() == 0) {
    Xdev->Xclearbox() ;
    Xdev->drawnbox = false ;
    return ;
  }
  if(Xdev->drawnbox) {
    drawfigure fig(grids.front()) ;
    fig.selectdev(*Xdev) ;
    fig.scalefigure() ;
    transform2d it = fig.figtransform.inverse_transform() ;
    positions xp1(Xdev->bx1,Xdev->by1), xp2(Xdev->bx2,Xdev->by2) ;
    positions pmx = it.apply_transform(xp1) ;
    positions pmn = it.apply_transform(xp2) ;
    if(pmx.x < pmn.x)
      std::swap(pmx.x,pmn.x) ;
    if(pmx.y < pmn.y)
      std::swap(pmx.y,pmn.y) ;
    positions pc = 0.5*( pmx + pmn) ;
    double rzoom_box = (pmx.x-pmn.x)/(pmx.y-pmn.y) ;
    DevBox b = Xdev->getDrawBox() ;
    double rdisp_box = fabs((b.xleft-b.xright)/(b.ytop-b.ybottom)) ;

    if(rzoom_box > rdisp_box) {
      double d = (pmx.x-pmn.x)/rdisp_box ;
      pmx.y = pc.y + 0.5*d ;
      pmn.y = pc.y - 0.5*d ;
    } else {
      double d = (pmx.y-pmn.y)*rdisp_box ;
      pmx.x = pc.x + 0.5*d ;
      pmn.x = pc.x - 0.5*d ;
    }

    const grid &current_grid = grids.front() ;
    grids.push_front(grid()) ;
    grids.front().zoom(current_grid,pmn,pmx) ;
    Xdev->drawnbox = false ;
    Xdev->clear() ;
    draw(grids.front(),*Xdev) ;
    Xdev->refresh() ;
  } else {
    if(grids.size() > 1) {
      grids.pop_front() ;
      Xdev->clear() ;
      draw(grids.front(),*Xdev) ;
      Xdev->refresh() ;
    }

  }

}


void PSFile(Widget w, XtPointer p1,XtPointer p2)
{
  if(grids.size()!= 0) {
    string filename = grids.front().filename ;

    filename += ".ps" ;
    std::cout << "psfile = " << filename << endl ;

    PSDevPSFile PSdev(filename) ;
    draw(grids.front(),PSdev) ;
  }
  popquit();
}
void EPSFile(Widget w, XtPointer p1,XtPointer p2)
{
  if(grids.size() != 0) {
    string filename = grids.front().filename ;
    filename += ".eps" ;
    std::cout << "epsfile = " << filename << endl ;

    PSDevEPS EPSdev(Xdev->getDrawBox(),filename) ;
    draw(grids.front(),EPSdev) ;
  }
  popquit();
}


void
warnok(Widget w, XtPointer userdata, XtPointer userdata2)
{
  XtRemoveGrab(pop) ;
  XtDestroyWidget(pop) ;
}

Widget
Label(char *string, Widget parent)
{
  static Arg labarg[] = {
    {XtNlabel, (XtArgVal) 0}
  } ;

  labarg[0].value = (XtArgVal) string ;
  return XtCreateManagedWidget("label",labelWidgetClass,parent,labarg,XtNumber(labarg)) ;
}



void
Warn(char *warnstring)
{
  static Arg pargs[] = {
    {XtNx,     (XtArgVal) 500},
    {XtNy,     (XtArgVal) 500},
  } ;
  Widget box ;

  pop = XtCreatePopupShell("Popup",overrideShellWidgetClass, toplevel,
                           pargs, XtNumber(pargs)) ;

  box = XtCreateManagedWidget("Box",boxWidgetClass,pop,NULL,0) ;
  Label(warnstring,box) ;
  Button("  OK  ", box, warnok) ;
  XtAddGrab(pop, True, True) ;
  XtManageChild(pop) ;
  XtRealizeWidget(pop) ;
  SetWMGetInput(pop) ;
}

void doload()
{
  char *s = XawDialogGetValueString(linput) ;

  if(index(s, '\n'))
    *index(s, '\n') = '\0' ;
  lpopquit() ;

  grids.clear() ;
  grids.push_front(grid()) ;
  grids.front().input(s) ;
  if(!grids.front().valid) {
    grids.clear() ;
    Xdev->clear() ;
    Xdev->refresh() ;
  } else {
    draw(grids.front(),*Xdev) ;
    Xdev->refresh() ;
  }
}

void doload(Widget w, XEvent *e, String *ss, Cardinal *c) {
  doload() ;
}

void doload(Widget w, XtPointer p1, XtPointer p2) {
  doload() ;
}


void dosetcontourspacing()
{
  char *s = XawDialogGetValueString(linput) ;

  while(*s == '\n')
    s++ ;
  double cs = atof(s) ;
  lpopquit() ;
  if(grids.size() == 0)
    return ;
  double contour_lim =  max((fabs(grids.front().max_val)+fabs(grids.front().min_val))*1e-8,1e-10) ;
  grids.front().generate_contour_curves(max(cs,contour_lim)) ;
  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}

void
dosetcontourspacing(Widget w, XEvent *e, String *ss, Cardinal *c)
{
  dosetcontourspacing() ;
}

void dosetcontourspacing(Widget w, XtPointer p1, XtPointer p2)
{
  dosetcontourspacing() ;
}

void dosetmaxcontour()
{
  char *s = XawDialogGetValueString(linput) ;

  while(*s == '\n')
    s++ ;
  double cs = atof(s) ;
  lpopquit() ;
  if(grids.size() == 0)
    return ;
  if(cs < grids.front().min_val) {
    grids.front().max_val = grids.front().min_val ;
    grids.front().min_val = cs ;
  } else
    grids.front().max_val = cs ;
  const double numcontours = 10.0 ;
  double range = (grids.front().max_val-grids.front().min_val)/numcontours ;
  double base = pow(10.0,double(floor(log10(range)))) ;
  double contour_lim =  max((fabs(grids.front().max_val)+fabs(grids.front().min_val))*1e-8,1e-10) ;
  grids.front().generate_contour_curves(max(floor(range/base)*base,contour_lim)) ;
  grids.front().generate_shading() ;
  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}

void
dosetmaxcontour(Widget w, XEvent *e, String *ss, Cardinal *c)
{
  dosetmaxcontour() ;
}

void dosetmaxcontour(Widget w, XtPointer p1, XtPointer p2)
{
  dosetmaxcontour() ;
}

void dosetmincontour()
{
  char *s = XawDialogGetValueString(linput) ;

  while(*s == '\n')
    s++ ;
  double cs = atof(s) ;
  lpopquit() ;
  if(grids.size() == 0)
    return ;
  if(cs > grids.front().max_val) {
    grids.front().min_val = grids.front().max_val ;
    grids.front().max_val = cs ;
  } else 
    grids.front().min_val = cs ;
  const double numcontours = 10.0 ;
  double range = (grids.front().max_val-grids.front().min_val)/numcontours ;
  double base = pow(10.0,double(floor(log10(range)))) ;
  double contour_lim =  max((fabs(grids.front().max_val)+fabs(grids.front().min_val))*1e-8,1e-10) ;
  grids.front().generate_contour_curves(max(floor(range/base)*base,contour_lim)) ;
  grids.front().generate_shading() ;
  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}

void
dosetmincontour(Widget w, XEvent *e, String *ss, Cardinal *c)
{
  dosetmincontour() ;
}

void dosetmincontour(Widget w, XtPointer p1, XtPointer p2)
{
  dosetmincontour() ;
}

void ContourSpacing(Widget w, XtPointer p1, XtPointer p2)
{
  static Arg pargs[] = {
    {XtNx,     (XtArgVal) 400},
    {XtNy,     (XtArgVal) 400},
    {0,0}
  } ;

  static char ContourEnter[] = "ContourEnter" ;
  static XtActionsRec actionTable[] = {
    {ContourEnter,   dosetcontourspacing },
  } ;
  static char texttrans[] = "\
<Key>Return:        ContourEnter() \n\
Ctrl<Key>M:         ContourEnter() \n\
Ctrl<Key>J:         ContourEnter() \n\
<Key>Linefeed:      ContourEnter() \n\
" ;


  static Arg dargs[] = {
    {XtNlabel, (XtArgVal) "Enter Contour Spacing:"},
    {XtNvalue, (XtArgVal) "1.0"}
  } ;

  if(grids.size() == 0)
    return ;
  XtSetArg(pargs[2],XtNtranslations,XtParseTranslationTable(texttrans)) ;
  static char defbuf[512] ;
  bzero(defbuf,512) ;
  snprintf(defbuf,511,"%s",fourSigDigs(grids.front().contour_spacing)) ;
  XtSetArg(dargs[1],XtNvalue,defbuf) ;
  lpop = XtCreatePopupShell("Popup",transientShellWidgetClass, toplevel,
                            pargs, XtNumber(pargs)) ;
  linput = XtCreateManagedWidget("file", dialogWidgetClass, lpop,
                                 dargs, XtNumber(dargs)) ;
  Button("  OK  ", linput, dosetcontourspacing) ;
  Button("CANCEL", linput, lpopquit) ;
  XtAppAddActions(XtWidgetToApplicationContext(linput),
                  actionTable, XtNumber(actionTable)) ;
  XtAddGrab(lpop, True, True) ;
  XtManageChild(lpop) ;
  XtRealizeWidget(lpop) ;
  SetWMGetInput(lpop) ;
}


void setMaxContour(Widget w, XtPointer p1, XtPointer p2)
{
  static Arg pargs[] = {
    {XtNx,     (XtArgVal) 400},
    {XtNy,     (XtArgVal) 400},
    {0,0}
  } ;

  static char MaxContourEnter[] = "MaxContourEnter" ;
  static XtActionsRec actionTable[] = {
    {MaxContourEnter,   dosetmaxcontour },
  } ;
  static char texttrans[] = "\
<Key>Return:        MaxContourEnter() \n\
Ctrl<Key>M:         MaxContourEnter() \n\
Ctrl<Key>J:         MaxContourEnter() \n\
<Key>Linefeed:      MaxContourEnter() \n\
" ;


  static Arg dargs[] = {
    {XtNlabel, (XtArgVal) "Enter Max Contour"},
    {XtNvalue, (XtArgVal) "1.0"}
  } ;

  if(grids.size() == 0)
    return ;
  XtSetArg(pargs[2],XtNtranslations,XtParseTranslationTable(texttrans)) ;
  static char defbuf[512] ;
  bzero(defbuf,512) ;
  snprintf(defbuf,511,"%s",fourSigDigs(grids.front().max_val)) ;
  XtSetArg(dargs[1],XtNvalue,defbuf) ;
  lpop = XtCreatePopupShell("Popup",transientShellWidgetClass, toplevel,
                            pargs, XtNumber(pargs)) ;
  linput = XtCreateManagedWidget("file", dialogWidgetClass, lpop,
                                 dargs, XtNumber(dargs)) ;
  Button("  OK  ", linput, dosetmaxcontour) ;
  Button("CANCEL", linput, lpopquit) ;
  XtAppAddActions(XtWidgetToApplicationContext(linput),
                  actionTable, XtNumber(actionTable)) ;
  XtAddGrab(lpop, True, True) ;
  XtManageChild(lpop) ;
  XtRealizeWidget(lpop) ;
  SetWMGetInput(lpop) ;
}

void setMinContour(Widget w, XtPointer p1, XtPointer p2)
{
  static Arg pargs[] = {
    {XtNx,     (XtArgVal) 400},
    {XtNy,     (XtArgVal) 400},
    {0,0}
  } ;

  static char MinContourEnter[] = "MinContourEnter" ;
  static XtActionsRec actionTable[] = {
    {MinContourEnter,   dosetmincontour },
  } ;
  static char texttrans[] = "\
<Key>Return:        MinContourEnter() \n\
Ctrl<Key>M:         MinContourEnter() \n\
Ctrl<Key>J:         MinContourEnter() \n\
<Key>Linefeed:      MinContourEnter() \n\
" ;


  static Arg dargs[] = {
    {XtNlabel, (XtArgVal) "Enter Min Contour"},
    {XtNvalue, (XtArgVal) "1.0"}
  } ;

  if(grids.size() == 0)
    return ;
  XtSetArg(pargs[2],XtNtranslations,XtParseTranslationTable(texttrans)) ;
  static char defbuf[512] ;
  bzero(defbuf,512) ;
  snprintf(defbuf,511,"%s",fourSigDigs(grids.front().min_val)) ;
  XtSetArg(dargs[1],XtNvalue,defbuf) ;
  lpop = XtCreatePopupShell("Popup",transientShellWidgetClass, toplevel,
                            pargs, XtNumber(pargs)) ;
  linput = XtCreateManagedWidget("file", dialogWidgetClass, lpop,
                                 dargs, XtNumber(dargs)) ;
  Button("  OK  ", linput, dosetmincontour) ;
  Button("CANCEL", linput, lpopquit) ;
  XtAppAddActions(XtWidgetToApplicationContext(linput),
                  actionTable, XtNumber(actionTable)) ;
  XtAddGrab(lpop, True, True) ;
  XtManageChild(lpop) ;
  XtRealizeWidget(lpop) ;
  SetWMGetInput(lpop) ;
}


void setoutputprinter()
{
  char *s = XawDialogGetValueString(linput) ;

  cerr << "setting pipecommand" << endl ;
  while(isspace(*s))
    s++ ;
  char pipecommand[512] ;

  for(int i=0;*s != '\0' && i<500;++i) {
    pipecommand[i] = *s ;
    s++ ;
    pipecommand[i+1] = '\0' ;
  }

  lpopquit();
  cerr << "printing " << endl ;
  if(grids.size()!= 0) {
    string filename = pipecommand ;

    PSDevPipe PSdev(filename) ;
    cerr << "printing to "<< filename << endl ;
    draw(grids.front(),PSdev) ;
  }
}

void setoutputprinter(Widget w, XEvent *e, String *ss, Cardinal *c) {
  setoutputprinter() ;
}


void setoutputprinter(Widget w, XtPointer p1, XtPointer p2) {
  setoutputprinter() ;
}


void
printerPick(Widget w, XtPointer p1, XtPointer p2)
{

  const char *pipedev;

  static Arg pargs[] = {
    {XtNx,     (XtArgVal) 400},
    {XtNy,     (XtArgVal) 400},
    {0,0}
  } ;

  static char PrintEnter[] = "PrintEnter" ;
  static XtActionsRec actionTable[] = {
    {PrintEnter,   setoutputprinter },
  } ;
  static char texttrans[] = "\
<Key>Return:        PrintEnter() \n\
Ctrl<Key>M:         PrintEnter() \n\
Ctrl<Key>J:         PrintEnter() \n\
<Key>Linefeed:      PrintEnter() \n\
" ;


  static Arg dargs[] = {
    {XtNlabel, (XtArgVal) "Enter Printer:"},
    {XtNvalue, (XtArgVal) "lpr -P"}
  } ;

  static char defbuf[512] ;

  if((pipedev = getenv("LWOUTPUT")) == NULL)
    pipedev = "lpr -P" ;

  //    if (outPrinterSet) strcpy (pipedev,outPrinter);

  bzero(defbuf,512) ;
  snprintf(defbuf,511,"%s",pipedev) ;
  XtSetArg(pargs[2],XtNtranslations,XtParseTranslationTable(texttrans)) ;
  XtSetArg(dargs[1],XtNvalue,defbuf) ;
  lpop = XtCreatePopupShell("Popup",transientShellWidgetClass, toplevel,
                            pargs, XtNumber(pargs)) ;
  linput = XtCreateManagedWidget("printer", dialogWidgetClass, lpop,
                                 dargs, XtNumber(dargs)) ;
  Button("  OK  ", linput, setoutputprinter) ;
  Button("CANCEL", linput, lpopquit) ;
  XtAppAddActions(XtWidgetToApplicationContext(linput),
                  actionTable, XtNumber(actionTable)) ;
  XtAddGrab(lpop, True, True) ;
  XtManageChild(lpop) ;
  XtRealizeWidget(lpop) ;
  SetWMGetInput(lpop) ;
}

void filestuff(Widget w, XtPointer p1, XtPointer p2)
{
  static Arg pargs[] = {
    {XtNx,     (XtArgVal) 500},
    {XtNy,     (XtArgVal) 500}
  } ;

  Widget box ;

  if(grids.size() == 0)
    return ;

  pop = XtCreatePopupShell("Popup",transientShellWidgetClass, toplevel,
                           pargs, XtNumber(pargs)) ;

  box = XtCreateManagedWidget("Box",boxWidgetClass,pop,NULL,0) ;
  Button("Printer",box, printerPick) ;
  Button("PostScript File", box, PSFile) ;
  Button("Encapsulated Postscript",box,EPSFile) ;
  Button("  CANCEL  ", box, popquit);
  XtAddGrab(pop, True, True) ;
  XtManageChild(pop) ;
  XtRealizeWidget(pop) ;
  SetWMGetInput(pop) ;
}



void Load(Widget w, XtPointer p1, XtPointer p2)
{
  static Arg pargs[] = {
    {XtNx,     (XtArgVal) 400},
    {XtNy,     (XtArgVal) 400},
    {0, 0}
  } ;
  static char LoadEnter[] = "LoadEnter" ;

  static XtActionsRec actionTable[] = {
    {LoadEnter,   doload },
  } ;
  static char texttrans[] = "\
<Key>Return:        LoadEnter() \n\
Ctrl<Key>M:         LoadEnter() \n\
Ctrl<Key>J:         LoadEnter() \n\
<Key>Linefeed:      LoadEnter() \n\
" ;

  static Arg dargs[] = {
    {XtNlabel, (XtArgVal) "       Enter filename:       "},
    {XtNvalue, (XtArgVal) ""}
  } ;



  XtSetArg(pargs[2],XtNtranslations,XtParseTranslationTable(texttrans)) ;
  lpop = XtCreatePopupShell("LoadFile",transientShellWidgetClass,
                            toplevel, pargs, XtNumber(pargs)) ;
  linput = XtCreateManagedWidget("file", dialogWidgetClass, lpop,
                                 dargs, XtNumber(dargs)) ;
  XtAppAddActions(XtWidgetToApplicationContext(linput),
                  actionTable, XtNumber(actionTable)) ;
  Button("  OK  ", linput, doload) ;
  Button("CANCEL", linput, lpopquit) ;
  XtAddGrab(linput, True, True) ;
  XtManageChild(lpop) ;
  XtRealizeWidget(lpop) ;
  SetWMGetInput(lpop) ;
}



static void Destroy(Widget w) {}
static void Resize(Widget w) {}
static Boolean SetValues(Widget o, Widget r, Widget n,
                         ArgList a, Cardinal *cp) {return 0;}
static Boolean TakeFocus(Widget w, Time *tp) {return 0; }


typedef struct _DrawClassRec {
  CoreClassPart   core_class ;
} DrawClassRec ;

typedef struct _DrawRec {
  CorePart core ;
} DrawRec, *DrawWidget ;

char JdrawString[] = "Jdraw" ;
DrawClassRec drawClassRec = {
  {
    /* core_class fields */
    /* superclass      */         (WidgetClass) &widgetClassRec,
    /* class_name      */         JdrawString,
    /* widget_size     */         sizeof(DrawRec),
    /* class_initialize*/         ClassInitialize,
    /* class_part_init */         NULL,
    /* class_inited    */         FALSE,
    /* initialize      */         Initialize,
    /* initialize_hook */         NULL,
    /* realize         */         Realize,
    /* actions         */         NULL,
    /* num_actions     */         (Cardinal) 0,
    /* resources       */         NULL,
    /* num_resources   */         (Cardinal) 0,
    /* xrm_class       */         NULLQUARK,
    /* compress_motion */         TRUE,
    /* compress_exposur*/         TRUE,
    /* compress_enterle*/         TRUE,
    /* visible_interest*/         FALSE,
    /* destroy         */         Destroy,
    /* resize          */         Resize,
    /* expose          */         HandleExpose,
    /* set_values      */         SetValues,
    /* set_values_hook */         NULL,
    /* set_values_almos*/         NULL,
    /* get_values_hook */         NULL,
    /* accept_focus    */         TakeFocus,
    /* version         */         XtVersion,
    /* callback_private*/         NULL,
    /* tm_table        */         NULL, /*  defaultTranslation, */
    /* query_geometry  */         NULL,
    /* display_accelera*/         NULL,
    /* externsion`     */         NULL
  }
} ;

WidgetClass drawWidgetClass = (WidgetClass) & drawClassRec ;

void ToggleShading(Widget w, XtPointer p1, XtPointer p2)
{
  if(grids.size() == 0)
    return ;
  grids.front().show_shading = !grids.front().show_shading ;
  if(grids.front().show_shading && !grids.front().shading_computed)
    grids.front().generate_shading() ;
  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}


void ToggleGrid(Widget w, XtPointer p1, XtPointer p2)
{
  if(grids.size() == 0)
    return ;
  grids.front().show_grid = !grids.front().show_grid ;
  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}

void ToggleRuler(Widget w, XtPointer p1, XtPointer p2)
{
  if(grids.size() == 0) {
    Xdev->Xclearbox() ;
    Xdev->drawnbox = false ;
    return ;
  }
  if(Xdev->drawnbox) {
    Xdev->Xclearbox() ;
    Xdev->drawnbox = false ;

    drawfigure fig(grids.front()) ;
    fig.selectdev(*Xdev) ;
    fig.scalefigure() ;
    transform2d it = fig.figtransform.inverse_transform() ;
    positions xp1(Xdev->bx1,Xdev->by1), xp2(Xdev->bx2,Xdev->by2) ;
    positions pmx = it.apply_transform(xp1) ;
    positions pmn = it.apply_transform(xp2) ;
    if(pmx.x < pmn.x)
      std::swap(pmx.x,pmn.x) ;
    if(pmx.y < pmn.y)
      std::swap(pmx.y,pmn.y) ;

    grids.front().rulers.push_back(rule(pmx,pmn)) ;
    fig.drawrulers() ;
    Xdev->refresh() ;
    return ;
  }

  if(grids.front().rulers.size() > 0)
    grids.front().rulers.clear() ;
  else {
    positions left = grids.back().minpos ;
    positions right = positions(grids.back().maxpos.x,grids.back().minpos.y) ;
    grids.front().rulers.push_back(rule(left,right)) ;
  }
  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}

// Right-click node probe. A second right-click while a probe is already
// active clears it rather than moving it -- a third click then starts
// fresh, finding whatever node is nearest that new position.
void ProbeNode(int screenX, int screenY)
{
  if(grids.size() == 0)
    return ;
  grid &g = grids.front() ;

  if(g.has_probe) {
    g.has_probe = false ;
  } else {
    drawfigure fig(g) ;
    fig.selectdev(*Xdev) ;
    fig.scalefigure() ;
    transform2d it = fig.figtransform.inverse_transform() ;
    positions clickPos = it.apply_transform(positions(double(screenX),double(screenY))) ;

    // Plain linear nearest-neighbor search over pos, skipping synthetic
    // cell-center points (is_node[i]==false -- see grid::input_2dgv/
    // input_generalized/input_cobalt) so this only ever snaps to a real
    // mesh vertex. zoom() already keeps pos trimmed to nodes near the
    // current view, so this stays fast without needing a spatial index.
    int best = -1 ;
    double bestDist2 = 0.0 ;
    for(size_t i=0;i<g.pos.size();++i) {
      if(!g.is_node[i])
        continue ;
      double dx = g.pos[i].x-clickPos.x ;
      double dy = g.pos[i].y-clickPos.y ;
      double d2 = dx*dx+dy*dy ;
      if(best < 0 || d2 < bestDist2) {
        best = int(i) ;
        bestDist2 = d2 ;
      }
    }
    if(best >= 0) {
      g.has_probe = true ;
      g.probe_pos = g.pos[best] ;
      g.probe_has_val = g.has_values ;
      g.probe_val = g.has_values ? g.val[best] : 0.0 ;
    }
  }

  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}

void ToggleContour(Widget w, XtPointer p1, XtPointer p2)
{
  if(grids.size() == 0)
    return ;
  grids.front().show_contours = !grids.front().show_contours ;
  Xdev->clear() ;
  draw(grids.front(),*Xdev) ;
  Xdev->refresh() ;
}

void hexbit(char *pt,int val){
  int v1 = val/16 ;
  int v2 = val%16 ;
  pt[0] = '0'+v1 ;
  if(v1>9)
    pt[0] = 'A'+v1-10 ;
  pt[1] = '0'+v2 ;
  if(v2>9)
    pt[1] = 'A'+v2-10 ;
}


void xinteract(int argc, char *argv[])
{
  static Arg args[] = {{XtNwidth, 0},{XtNheight, 0},};

  toplevel = XtInitialize("2dgv", "2dgv", opt_table, XtNumber(opt_table),
			  &argc, argv);

  for(int i=0;i<MAXPENS-1;++i) {
    char c1 = '0' + i/10 ;
    char c2 = '0' + i%10 ;
    STR::xpen[i][4] = c1 ;
    STR::xpen[i][5] = c2 ;
    double s = (double)(i)/(MAXPENS-2) ;

    //    int R = min(int(128*(1.+cos(thr))),255) ;
    int R = min(int(800*s),255) ;
    //    int G = min(int(128*(1.+cos(thg))),255) ;
    int G = min(int(256*s*s),255) ;
    int B = max(min(50+int(512*(0.5-s)),255),
                min(int(512*(s-0.5)),255)) ;
    //    int B = min(int(128*(1.+cos(thb))),255) ;
    X_R[i] = R ;
    X_G[i] = G ;
    X_B[i] = B ;
    hexbit(&STR::col[i][1],R) ;
    hexbit(&STR::col[i][3],G) ;
    hexbit(&STR::col[i][5],B) ;
  }

  XtGetApplicationResources(toplevel, (XtPointer)&appRes,
                            appResourceSpec, XtNumber(appResourceSpec),
                            NULL, 0);

  args[0].value = (XtArgVal) appRes.suggest_width;
  args[1].value = (XtArgVal) appRes.suggest_height+20 ;


  Widget pane = XtCreateManagedWidget("pane", panedWidgetClass, toplevel,
                                      NULL,0) ;
  Widget box  = XtCreateManagedWidget("commandbox", boxWidgetClass, pane,
                                      NULL,0) ;
  Button("Load", box, Load) ;
  Button("Output",box, filestuff) ;
  Button("Contour",box,ToggleContour) ;
  Button("C Spacing",box,ContourSpacing) ;
  Button("C Max",box,setMaxContour) ;
  Button("C Min",box,setMinContour) ;
  Button("Grid",box,ToggleGrid) ;
  Button("Ruler",box,ToggleRuler) ;
  Button("Shading",box,ToggleShading) ;
  Button("Zoom",box,Zoomin) ;

  Widget drawwidget = XtCreateWidget("field", drawWidgetClass, pane,
                                     args, XtNumber(args));
  Xdev = new XDevice(drawwidget) ;

  Button("Quit", box, Quit) ;

  XtRealizeWidget(toplevel);

  Display *mainDisplay=XtDisplay(toplevel) ;
  int mainScreen = DefaultScreen(mainDisplay) ;
  Window rootWindow = RootWindow(mainDisplay,mainScreen) ;
  XSelectInput(mainDisplay, rootWindow, KeyPressMask | StructureNotifyMask
	       | VisibilityChangeMask | FocusChangeMask);
  XMapWindow(mainDisplay, rootWindow);
  XWMHints wmhints;

  wmhints.flags = InputHint | StateHint ;
  wmhints.input = True;
  wmhints.icon_pixmap = None;
  wmhints.initial_state = NormalState;
  Window topwin = XtWindowOfObject(toplevel) ;
  XSetWMHints(mainDisplay, topwin, &wmhints);

  XtMainLoop();
}


XDevice::XDevice(Widget dw) {
  drawwindow = XtWindow(dw) ;
  currentPen = 0 ;
  ticksize = 5. ;
  dpy = XtDisplay(dw) ;
  pixmap_allocated = false ;
  drawnbox = false ;

  drawdepth = dw->core.depth ;

  XtAddEventHandler(dw,ButtonPressMask|ButtonReleaseMask|ButtonMotionMask,
                    False,(XtEventHandler)mouseButtonPress,NULL) ;

  XtAddEventHandler(dw,StructureNotifyMask, False, ResizeEvent, NULL) ;

  XGCValues values ;
  values.foreground = appRes.defaultfore ;
  values.background = dw->core.background_pixel ;

  linec = XtGetGC(dw, (XtGCMask) GCForeground | GCBackground, &values) ;
  values.foreground =  appRes.defaultback ;
  forec = XtGetGC(dw, (XtGCMask) GCForeground | GCBackground, &values) ;
  //  XFontStruct *font_info = XLoadQueryFont(dpy,"6x13") ;
  //  XSetFont(dpy,linec,font_info->fid) ;
  for(int i=0; i < MAXPENS; i++) {
    values.foreground = appRes.Pens[i] ;
    Penc[i] = XtGetGC(dw, (XtGCMask) GCForeground | GCBackground, &values) ;
  }


  XtManageChild(dw);
}

void XDevice::OpenDevice()
{
  clear() ;
}

void XDevice::CloseDevice()
{
}

void XDevice::HardCloseDevice()
{
}

void XDevice::DrawLine(const positions &p1,const positions &p2)
{
  int x1 = (int(p1.x+.5)) ;
  int y1 = (int(p1.y+.5)) ;
  int x2 = (int(p2.x+.5)) ;
  int y2 = (int(p2.y+.5)) ;

  XDrawLine(dpy,drawpixmap,Penc[currentPen],x1,y1,x2,y2) ;
}


void XDevice::SelectPen(int color)
{
  if(color < 0 || color >= MAXPENS) {
    cerr <<" invalid pen # = " << color << endl ;
    return ;
  }
  currentPen = color ;
}

void XDevice::FillShade(int Shade)
{
  SelectPen(Shade);
}

void XDevice::DrawLabel(const std::string &text, positions location,
                        int pos, int size, double degrees)
{
  int x = int(location.x+.5) ;
  int y = int(location.y+.5) ;
  XDrawString(dpy, drawpixmap, Penc[currentPen], x,y, text.c_str(),text.size());
}

void XDevice::FillPoly(const positions *pnts, int num_pnts)
{
  XPoint arry[100] ;

  if (currentPen == (MAXPENS - 1)) currentPen--;
  if (currentPen == 0) currentPen = (MAXPENS - 1);


  for(int i=0;i<num_pnts;i++) {
    arry[i].x =  int(pnts[i].x + .5) ;
    arry[i].y =  int(pnts[i].y + .5) ;
  }

  XFillPolygon(dpy, drawpixmap, Penc[currentPen], arry, num_pnts,
               Convex, CoordModeOrigin);
}

// Same as FillPoly above but uses the dedicated markerGC (see realize())
// instead of the cycling Penc[] pen table, and doesn't need FillPoly's
// currentPen adjustment since it's not drawing with a pen index at all.
void XDevice::FillPolyMarker(const positions *pnts, int num_pnts)
{
  XPoint arry[100] ;

  for(int i=0;i<num_pnts;i++) {
    arry[i].x =  int(pnts[i].x + .5) ;
    arry[i].y =  int(pnts[i].y + .5) ;
  }

  XFillPolygon(dpy, drawpixmap, markerGC, arry, num_pnts,
               Convex, CoordModeOrigin);
}
