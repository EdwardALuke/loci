#include "draw.h"
#include "debug.h"
#include "tools.h"
#include "pens.h"

#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>

using std::min ;
using std::max ;
//using std::abs ;
using std::cerr ;
using std::cout ;
using std::endl ;

inline double sgn(double x) { return x>0?1.0:-1.0; }

const double EPSILON = 1e-20 ;

void drawfigure::scalefigure()
{
  const double minx = figure.minview.x ;
  const double miny = figure.minview.y ;
  const double maxx = figure.maxview.x ;
  const double maxy = figure.maxview.y ;
  double ticksize = currdev->getTickSize() ;
  double sgnx = sgn(currdev->getDrawBox().xright-currdev->getDrawBox().xleft) ;
  double sgny = sgn(currdev->getDrawBox().ytop-currdev->getDrawBox().ybottom) ;
  const double xleft   = currdev->getDrawBox().xleft   + sgnx*ticksize ;
  const double xright  = currdev->getDrawBox().xright  - sgnx*ticksize ;
  const double ytop    = currdev->getDrawBox().ytop    - sgny*ticksize ;
  const double adj = (currdev->isLabelDev())?8.0:1.0 ;
  const double ybottom = currdev->getDrawBox().ybottom + sgny*adj*ticksize ;
  titlepos = positions(xleft,ybottom-sgny*7.0*ticksize) ;
  legendpos = positions(xright-sgnx*30.0*ticksize,(ytop+ybottom)/2.0+sgny*double(MAXPENS/2)*ticksize*2.0) ;
  
  const double tscalex = (xright-xleft)/(maxx-minx+EPSILON) ;
  const double tscaley = (ytop-ybottom)/(maxy-miny+EPSILON) ;
  const double scale = min(fabs(tscalex),fabs(tscaley)) ;
  const double scalex = scale *sgn(tscalex) ;
  const double scaley = scale *sgn(tscaley) ;

  // Set scaleing factor
  figtransform.identity() ;
  figtransform.scale_transform(scalex,scaley) ;

  double basex = minx ;
  double basey = miny ;
  double offsetx = xleft ;
  double offsety = ybottom ;
  if(fabs(tscalex) < fabs(tscaley)) 
    offsety += (ytop-ybottom-fabs(maxy-miny)*scale*sgny)/2.0 ;
  else 
    offsetx += (xright-xleft-fabs(maxx-minx)*scale*sgnx)/2.0 ;

  // Pan transform to center of drawing space 
  figtransform.pan_transform(positions(offsetx-basex*scalex,
                                       offsety-basey*scaley)) ;
}

void drawfigure::drawgrid() {
  int nedges = figure.edge_list.size() ;

  SelectPen(0) ;
  for(int e=0;e<nedges;++e) {
    const edges &ed = figure.edge_list[e] ;
    FATAL(ed.l >= npos || ed.r >= npos) ;
    DrawLine(figure.pos[ed.l],figure.pos[ed.r]) ;
  }
}

void drawfigure::drawborder() {
  int nedges = figure.edge_list.size() ;

  SelectPen(0) ;
  for(int e=figure.interior;e<nedges;++e) {
    const edges &ed = figure.edge_list[e] ;
    FATAL(ed.l >= npos || ed.r >= npos) ;
    DrawLine(figure.pos[ed.l],figure.pos[ed.r]) ;
  }
}

void drawfigure::drawtriangles() {
  int ntri = figure.triangle_list.size() ;

  SelectPen(0) ;
  for(int t=0;t<ntri;++t) {
    const triangles &td = figure.triangle_list[t] ;
    FATAL(td.t1 >= npos || td.t2 >= npos  || td.t3>=npos) ;
    DrawLine(figure.pos[td.t1],figure.pos[td.t2]) ;
    DrawLine(figure.pos[td.t2],figure.pos[td.t3]) ;
    DrawLine(figure.pos[td.t3],figure.pos[td.t1]) ;
  }
}

void drawfigure::drawcontours()
{
  SelectPen(0) ;
  for(size_t i=0;i<figure.contour_curves.size();++i) {
    const segments &si=figure.contour_curves[i] ;
    DrawLine(si.p1,si.p2) ;
  }
  
  for(size_t i=0;i<figure.contour_values.size();++i) {
    if(figure.contour_values[i].mingrad < 100.0) {
      char buf[256] ;
      bzero(buf,256) ;
      snprintf(buf,255,"%s",fourSigDigs(figure.contour_values[i].value)) ;
      DrawLabel(buf,figure.contour_values[i].mingradpos,0,10,0) ;
    }
  }
}

void drawfigure::drawshading()
{
  for(int i=0;i<MAXPENS;++i) {
    if(figure.pnts[i].size() > 0) {
      SelectPen(i) ;
      vector<positions>::const_iterator pi = figure.pnts[i].begin();
      for(vector<int>::const_iterator vi=figure.pntindex[i].begin()+1;
          vi!=figure.pntindex[i].end();++vi) {
        int num_pnts = *vi - *(vi-1)  ;
        FillPoly(pi,num_pnts) ;
        pi = pi + num_pnts ;
      }
    }
  }    
}

void drawfigure::drawtitle() {
  char buf[512] ;
  char val[256] ;
  char date[21] ;
  time_t tloc ;

  time(&tloc) ;
  strncpy(buf,ctime(&tloc),4) ;
  strncpy(date,buf+4,12) ;
  strncpy(date+12,buf+19,5) ;
  date[17]='\0' ;

  if(!figure.valid) {
    snprintf(buf,512,"No grid available") ;
  } else if(figure.has_values) {
    snprintf(val,256,"%s",fourSigDigs(figure.contour_spacing)) ;
    snprintf(buf,512,"%s contour spacing = %s [Function: max=%4g min=%4g] Date: %s",
            figure.filename.c_str(),val,figure.max_val,figure.min_val,date) ;
  } else {
    snprintf(buf,512,"Grid filename: %s        Date: %s",figure.filename.c_str(),date) ;
  }
  currdev->DrawLabel(buf,titlepos,3,10,0) ;
}

void drawfigure::drawlegend()
{
  double ticksize = currdev->getTickSize()*2.0 ;
  positions rect[4] ;
  for(int i=1;i<MAXPENS;++i) {
    SelectPen(i) ;
    for(int j=0;j<4;++j) 
      rect[j] = positions(legendpos.x,legendpos.y+ticksize*double(i-1)) ;
    
    rect[1].x += ticksize*2.0 ;
    rect[2].x += ticksize*2.0 ;
    rect[2].y += ticksize ;
    rect[3].y += ticksize ;
    currdev->FillPoly(rect,4) ;
  }
  double maxv = figure.max_val ;
  double minv = figure.min_val ;
  for(size_t i=0;i<figure.contour_values.size();++i) {
    double v = figure.contour_values[i].value ;
    double s = (v-minv)/(maxv-minv) ;
    if(s<0.0 || s>1.0)
      continue ;
    positions p1(legendpos.x,legendpos.y+s*ticksize*double(MAXPENS-1.01)) ;
    positions p2(p1.x+ticksize*2.0,p1.y) ;
    positions p3(p2.x+ticksize,p1.y) ;
    SelectPen(int(1+s*double(MAXPENS-1.01))) ;
    currdev->DrawLine(p1,p2) ;
    char buf[256] ;
    snprintf(buf,256,"%s",fourSigDigs(figure.contour_values[i].value)) ;
    currdev->DrawLabel(buf,p3,0,10,0) ;
  }
  
}

// Pick a "nice" (floor-to-power-of-ten) tick increment sized so that
// ticks land roughly every target_px screen pixels along a ruler whose
// on-screen length is pixel_length (and whose length in data units is
// data_length). Returns 0.0 for a degenerate (zero/negative/non-finite)
// input so callers can skip tick drawing instead of propagating NaN.
static double tickIncrement(double data_length, double pixel_length,
                             double target_px)
{
  if(data_length <= 0.0 || !isfinite(data_length) ||
     pixel_length <= 0.0 || !isfinite(pixel_length))
    return 0.0 ;
  double num_graduations = pixel_length/target_px ;
  if(num_graduations < 1.0)
    num_graduations = 1.0 ; // always allow at least one tick
  double step = data_length/num_graduations ;
  if(step <= 0.0 || !isfinite(step))
    return 0.0 ;
  double base = pow(10.0,floor(log10(step))) ;
  double inc = floor(step/base)*base ;
  if(inc <= 0.0 || !isfinite(inc))
    return 0.0 ;
  return inc ;
}

// "9x15" is a fixed-width X11 core font, and the name is exactly its
// cell size in pixels, loaded once at startup by XDevice::realize().
// XDevice::DrawLabel ignores the `size` argument entirely and always
// uses this font regardless of what's requested (see the comment on
// DrawLabel's `pos` argument in drawrulers() below for the related
// justification quirk), so these are the font's real, exact dimensions,
// not estimates. Shared by drawrulers() and drawprobe().
static const double CHAR_WIDTH_PX  = 9. ;
static const double CHAR_HEIGHT_PX = 15. ;

// Computes the screen scale (pixels per data unit, x and y) and the true
// visible data-space viewport for the current transform. The viewport
// can be LARGER than figure.minview/maxview -- scalefigure() preserves
// the mesh's aspect ratio by applying the same scale factor to both axes
// (the smaller of the two fits), so whichever axis doesn't need the full
// drawable area gets centered with a margin on both sides. That margin
// is real, currently visible screen space, but figure.minview/maxview
// only ever holds the "requested" view bounds -- initially just the
// mesh's own bounding box -- not this actual rendered extent. Inverse-
// transforming the device's actual draw box gives the true bounds
// instead (the same approach ToggleRuler/Zoomin use to convert a dragged
// screen box into data space). Shared by drawrulers() and drawprobe().
void drawfigure::computeViewport(double &xscale, double &yscale,
                                  double &viewMinX, double &viewMaxX,
                                  double &viewMinY, double &viewMaxY)
{
  positions t1(0.,0.),t2(1.,0.),p1,p2 ;
  p1 = figtransform.apply_transform(t1) ;
  p2 = figtransform.apply_transform(t2) ;
  positions d = p1-p2 ;
  xscale = sqrt(d.x*d.x+d.y*d.y) ;
  t2 = positions(0.,1.) ;
  p2 = figtransform.apply_transform(t2) ;
  d = p1-p2 ;
  yscale = sqrt(d.x*d.x+d.y*d.y) ;

  const DevBox &box = currdev->getDrawBox() ;
  transform2d itrans = figtransform.inverse_transform() ;
  positions vc1 = itrans.apply_transform(positions(box.xleft,box.ybottom)) ;
  positions vc2 = itrans.apply_transform(positions(box.xright,box.ytop)) ;
  viewMinX = std::min(vc1.x,vc2.x) ;
  viewMaxX = std::max(vc1.x,vc2.x) ;
  viewMinY = std::min(vc1.y,vc2.y) ;
  viewMaxY = std::max(vc1.y,vc2.y) ;
}

void drawfigure::drawrulers()
{
  double ticksize = currdev->getTickSize() ;
  positions p1,p2 ;
  double xscale,yscale,viewMinX,viewMaxX,viewMinY,viewMaxY ;
  computeViewport(xscale,yscale,viewMinX,viewMaxX,viewMinY,viewMaxY) ;

  // Cap on how many screen pixels a ruler's tick spacing can be computed
  // over: the device's actual drawable size. Without this, an old ruler
  // that's much longer than the current (zoomed-in) view can compute an
  // enormous "pixel length" from its full data-space extent, producing
  // an absurdly fine tick spacing and a huge number of ticks.
  const DevBox &box = currdev->getDrawBox() ;
  const double maxPxX = fabs(box.xright-box.xleft) ;
  const double maxPxY = fabs(box.ytop-box.ybottom) ;

  // DrawLabel's `pos` argument (a text-justification code -- see the
  // PSDevCore::DrawLabel switch statement for the full set) is only
  // honored by the PostScript backend. XDevice::DrawLabel, used for the
  // interactive window, ignores both `pos` and `size` -- it always draws
  // left-justified (a plain XDrawString call) using one fixed font (see
  // CHAR_WIDTH_PX/CHAR_HEIGHT_PX above).

  // Target on-screen spacing between ticks. Vertical ruler labels are
  // always exactly one line tall regardless of how many characters they
  // have, so a spacing derived directly from the known font height
  // guarantees adjacent labels can never overlap vertically. Horizontal
  // labels are a different story -- their width grows with digit count
  // (more decimals, a sign, scientific notation) -- so target_px_x below
  // is only a starting point; it gets widened per ruler once the actual
  // widest label at that ruler's precision is known (see "guaranteed"
  // comment further down).
  const double target_px_x = 130. ;
  const double target_px_y = CHAR_HEIGHT_PX + 2.*CHAR_WIDTH_PX ;

  // A double has about 15-17 reliable significant decimal digits.
  // Several floating point operations happen upstream of any tick value
  // used below (the inverse view transform above, subtracting nearly-
  // equal view bounds at deep zoom, etc.), each able to erode a bit of
  // that precision, so a conservative floor well short of the
  // theoretical machine epsilon (~2.22e-16 relative) is used: an
  // increment finer than a ruler's own magnitude times this floor
  // doesn't represent real, distinguishable ticks, no matter how many
  // decimal places are printed. Shared by both branches below.
  const double MIN_RELATIVE_INCREMENT = 1e-12 ;

  // Hard cap on interior ticks per ruler, shared by both branches: a
  // defensive backstop against floating point precision loss at extreme
  // zoom (nearly-identical view bounds can make the tick math drift),
  // well above what any reasonable target_px should ever generate.
  const int MAX_TICKS = 200 ;

  char buf[512] ;
  std::vector<rule>::const_iterator i ;
  for(i=figure.rulers.begin();i!=figure.rulers.end();++i) {
    p1 = i->start ;
    double len = i->length ;
    if(i->xdir) {
      p2 = p1 + len*positions(1.,0.) ;

      // Clip to the ruler's overlap with the current view. This keeps
      // the tick increment based on what's actually on screen (so it
      // rescales correctly as you zoom, and can't compute a spacing far
      // A ruler is only actually on screen if BOTH its own extent
      // overlaps the view along its own axis AND its fixed cross-axis
      // position is within the view along the other axis. Only the
      // former was checked before -- a horizontal ruler whose X-extent
      // happens to still overlap the view, but whose fixed Y-position
      // (p1.y) has scrolled outside the view's Y-range, would still be
      // drawn: its baseline and labels get positioned at that far
      // out-of-range Y, which can transform to an extreme, effectively
      // arbitrary screen coordinate (this is what caused a vertical
      // ruler's labels to show up detached, elsewhere on screen, once
      // zoomed away from its own position -- the same failure mode is
      // possible here, just swapping which axis is "along the ruler"
      // and which is "fixed").
      const double visStartX = std::max(i->start.x,viewMinX) ;
      const double visEndX   = std::min(p2.x,viewMaxX) ;
      const double visLenX   = std::max(0.0,visEndX-visStartX) ;
      if(visLenX <= 0.0)
        continue ; // ruler doesn't overlap the current view at all
      if(p1.y < viewMinY || p1.y > viewMaxY)
        continue ; // ruler's fixed Y-position has scrolled off screen

      DrawLine(positions(visStartX,p1.y),positions(visEndX,p1.y)) ;

      double targetX = target_px_x ;
      const double pxLenX = std::min(visLenX*xscale,maxPxX) ;
      double incx = tickIncrement(visLenX,pxLenX,targetX) ;

      // Clamping up front (rather than trying to hide the problem after
      // the fact) means the ruler still shows as much genuine resolution
      // as it can, and gracefully shows just its endpoints once zoomed
      // in past what's actually resolvable, instead of either
      // fabricating precision-free ticks or blanking out entirely.
      const double repMagX = std::max(fabs(visStartX),fabs(visEndX)) ;
      const double minIncrementX = repMagX*MIN_RELATIVE_INCREMENT ;
      if(incx > 0.0 && incx < minIncrementX)
        incx = minIncrementX ;

      // Guaranteed no-overlap: checking only the two endpoints isn't
      // enough, because rulerDigits() decides scientific-vs-fixed
      // notation independently for EACH tick from its own value. A
      // ruler that straddles zero can have an INTERIOR tick land in the
      // near-zero band (0 < |value| < 0.009) and get formatted in
      // scientific notation -- a much longer string -- even when neither
      // endpoint is anywhere close to zero. Only checking the endpoints
      // misses that worst case entirely (this only bites when the
      // increment happens to be fine enough for a tick to actually land
      // in that band, which is why it only showed up at some zoom
      // levels). So every candidate tick that would actually be drawn at
      // this increment is checked directly instead of trying to guess
      // which ones bound the worst case.
      //
      // tickIncrement always rounds DOWN to a "nice" value, so the
      // actual resulting on-screen spacing (incx*xscale) is not the same
      // as the target_px it was given -- in the worst case it can be as
      // little as half of it. A single "bump target_px and recompute"
      // pass isn't a real guarantee, since the recomputed incx can still
      // fall short by the same mechanism. This verifies the ACTUAL
      // spacing directly and keeps widening the target until it's
      // genuinely satisfied (bounded to a handful of iterations; each
      // step can only make incx equal or coarser, so this converges).
      // neededPx is kept for minGap below, since it's the real, verified
      // minimum -- targetX is just the input knob used to get there.
      double neededPx = targetX ;
      if(incx > 0.0) {
        for(int guard=0; guard<8; ++guard) {
          int worstLen = 0 ;
          const double sn = ceil(visStartX/incx) ;
          for(int n=0; n<200; ++n) {
            double xval = (sn+n)*incx ;
            if(xval > visEndX + incx*1e-9)
              break ;
            int l = (int)strlen(rulerDigits(xval,incx)) ;
            if(l > worstLen)
              worstLen = l ;
          }
          // Endpoints explicitly too, in case they don't land exactly on
          // a "nice" grid point covered by the loop above.
          int lenS = (int)strlen(rulerDigits(visStartX,incx)) ;
          int lenE = (int)strlen(rulerDigits(visEndX,incx)) ;
          if(lenS > worstLen) worstLen = lenS ;
          if(lenE > worstLen) worstLen = lenE ;

          neededPx = CHAR_WIDTH_PX*worstLen + 2.*CHAR_WIDTH_PX ;
          const double actualPx = incx*xscale ;
          if(actualPx >= neededPx)
            break ;
          // tickIncrement's output only changes when step/base crosses
          // an integer boundary -- a flat multiplier can land just short
          // of that boundary and produce the exact same incx again,
          // which the "no improvement" check below would then wrongly
          // treat as "nothing more can be done" while still under-
          // spaced. Scaling by the actual measured shortfall (with a 2x
          // safety margin, since the boundary itself isn't known
          // in advance) reliably crosses it instead.
          const double bump = std::max(neededPx/std::max(actualPx,1e-300)*2.0,1.5) ;
          targetX = std::max(targetX,neededPx)*bump ;
          double newIncx = tickIncrement(visLenX,pxLenX,targetX) ;
          if(newIncx <= incx)
            break ; // no further improvement possible (e.g. already one tick)
          incx = newIncx ;
        }
        // No further fallback needed beyond this point: the minimum-
        // increment clamp above already guarantees incx reflects real,
        // distinguishable precision, and if the visible range is
        // narrower than even that minimum, the interior tick loop below
        // naturally produces zero ticks (its first candidate already
        // exceeds visEndX) -- gracefully leaving just the two endpoint
        // ticks, which are drawn unconditionally regardless.
      }

      // Draws one tick + label at data-space x value `xval`. Endpoint
      // ticks cross the ruler line on both sides with their label above
      // it; interior ticks are single-sided with their label below, so
      // the two never collide.
      auto drawTickX = [&](double xval, double inc, bool endpoint){
        if(endpoint) {
          positions ta = positions(xval,p1.y-ticksize/yscale) ;
          positions tb = positions(xval,p1.y+ticksize/yscale) ;
          DrawLine(ta,tb) ;
          snprintf(buf,256,"%s",rulerDigits(xval,inc)) ;
          DrawLabel(buf,positions(xval,p1.y+4.*ticksize/yscale),5,10,0) ;
        } else {
          positions t(xval,p1.y) ;
          positions t3 = positions(xval,p1.y-ticksize/yscale) ;
          DrawLine(t,t3) ;
          snprintf(buf,256,"%s",rulerDigits(xval,inc)) ;
          DrawLabel(buf,positions(xval,p1.y-4.*ticksize/yscale),5,10,0) ;
        }
      } ;

      // Always tick+label the ruler's exact start/end when they're on
      // screen, so the labeled range visibly matches the drawn line
      // rather than relying on a "nice" interior tick happening to land
      // near the edge. Interior ticks fill in between.
      if(i->start.x >= viewMinX && i->start.x <= viewMaxX)
        drawTickX(i->start.x, incx>0.0?incx:len, true) ;
      if(len > 0.0 && p2.x >= viewMinX && p2.x <= viewMaxX)
        drawTickX(p2.x, incx>0.0?incx:len, true) ;

      if(incx > 0.0) {
        // Minimum gap an interior tick must keep from the visible edge
        // so it doesn't collide with the endpoint tick's label.
        //
        // This used to be half of targetX. That assumed labels are
        // roughly centered on their tick, but XDevice::DrawLabel ignores
        // the `pos` justification argument entirely and always draws
        // left-justified -- so a label actually extends its FULL width
        // away from its own tick, not half of it. An endpoint's label,
        // drawn at the visible edge, reaches its full rendered width
        // into the interior; halving the clearance let an interior tick
        // (and its own label) land inside that reach, most visibly for
        // labels one character wider than expected -- a leading "-", or
        // scientific notation's extra characters. neededPx above is the
        // actual verified minimum spacing (already includes a gap), so
        // it's used directly here rather than re-deriving a fraction of
        // it.
        const double minGap = neededPx/xscale ;
        const double startn = ceil(visStartX/incx) ;
        for(int n=0; n<MAX_TICKS; ++n) {
          double xval = (startn+n)*incx ;
          if(xval > visEndX - minGap + incx*1e-9)
            break ;
          if(fabs(xval) < incx*1e-3)
            xval = 0. ; // snap negligible floating point noise near zero
          if(xval > visStartX + minGap)
            drawTickX(xval,incx,false) ;
        }
      }
    } else {
      p2 = p1 + len*positions(0.,1.) ;

      // Same as the x-direction branch above, mirrored for y -- and same
      // cross-axis check: the ruler's fixed X-position (p1.x) must also
      // be within the view, not just its Y-extent overlapping it.
      const double visStartY = std::max(i->start.y,viewMinY) ;
      const double visEndY   = std::min(p2.y,viewMaxY) ;
      const double visLenY   = std::max(0.0,visEndY-visStartY) ;
      if(visLenY <= 0.0)
        continue ;
      if(p1.x < viewMinX || p1.x > viewMaxX)
        continue ; // ruler's fixed X-position has scrolled off screen

      DrawLine(positions(p1.x,visStartY),positions(p1.x,visEndY)) ;

      double targetY = target_px_y ;
      const double pxLenY = std::min(visLenY*yscale,maxPxY) ;
      double incy = tickIncrement(visLenY,pxLenY,targetY) ;

      // Same double-precision floor as the x-direction branch (constant
      // defined once, above).
      const double repMagY = std::max(fabs(visStartY),fabs(visEndY)) ;
      const double minIncrementY = repMagY*MIN_RELATIVE_INCREMENT ;
      if(incy > 0.0 && incy < minIncrementY)
        incy = minIncrementY ;

      // Same reasoning as the x-direction branch: a vertical label's
      // height is fixed (one font line, unlike width it doesn't grow
      // with digit count), but tickIncrement's own floor-rounding can
      // still undershoot target_px_y by up to half, so the actual
      // spacing is verified directly rather than assumed.
      double neededPxY = targetY ;
      if(incy > 0.0) {
        for(int guard=0; guard<8; ++guard) {
          neededPxY = target_px_y ;
          const double actualPxY = incy*yscale ;
          if(actualPxY >= neededPxY)
            break ;
          const double bumpY = std::max(neededPxY/std::max(actualPxY,1e-300)*2.0,1.5) ;
          targetY = std::max(targetY,neededPxY)*bumpY ;
          double newIncy = tickIncrement(visLenY,pxLenY,targetY) ;
          if(newIncy <= incy)
            break ;
          incy = newIncy ;
        }
        // As in the x-direction branch, no further fallback is needed:
        // the minimum-increment clamp above already guarantees incy
        // reflects real, distinguishable precision.
      }

      // Endpoint ticks/labels mirror drawTickX, with the label on the
      // left instead of above. Interior labels sit to the right of their
      // tick and grow further rightward (harmless, since that's away
      // from the tick), but an endpoint label needs to grow leftward,
      // away from the tick, instead -- otherwise it runs back over the
      // tick mark it belongs to. Since neither backend gives us the
      // label's actual rendered width up front, the right edge is
      // pinned at a fixed gap left of the tick by shifting the anchor
      // left by the label's own known width (CHAR_WIDTH_PX above),
      // so this holds for any number of digits or notation.
      auto drawTickY = [&](double yval, double inc, bool endpoint){
        if(endpoint) {
          positions ta = positions(p1.x-ticksize/xscale,yval) ;
          positions tb = positions(p1.x+ticksize/xscale,yval) ;
          DrawLine(ta,tb) ;
          snprintf(buf,256,"%s",rulerDigits(yval,inc)) ;
          double textWidth = CHAR_WIDTH_PX*(double)strlen(buf)/xscale ;
          DrawLabel(buf,positions(p1.x-2.*ticksize/xscale-textWidth,yval),5,10,0) ;
        } else {
          positions t(p1.x,yval) ;
          positions t3 = positions(p1.x+ticksize/xscale,yval) ;
          DrawLine(t,t3) ;
          snprintf(buf,256,"%s",rulerDigits(yval,inc)) ;
          DrawLabel(buf,positions(p1.x+2.*ticksize/xscale,yval),5,10,0) ;
        }
      } ;

      if(i->start.y >= viewMinY && i->start.y <= viewMaxY)
        drawTickY(i->start.y, incy>0.0?incy:len, true) ;
      if(len > 0.0 && p2.y >= viewMinY && p2.y <= viewMaxY)
        drawTickY(p2.y, incy>0.0?incy:len, true) ;

      if(incy > 0.0) {
        const double minGap = neededPxY/yscale ;
        const double startn = ceil(visStartY/incy) ;
        for(int n=0; n<MAX_TICKS; ++n) {
          double yval = (startn+n)*incy ;
          if(yval > visEndY - minGap + incy*1e-9)
            break ;
          if(fabs(yval) < incy*1e-3)
            yval = 0. ;
          if(yval > visStartY + minGap)
            drawTickY(yval,incy,false) ;
        }
      }
    }
  }
}

// Draws the right-click node probe: a small marker circle at the
// snapped node, and a val/posx/posy info list anchored on whichever
// side of the node keeps it clear of the node and (as much as possible)
// inside the visible plot area. The anchor side is chosen per-axis: the
// list goes below a node in the top half of the plot (and above one in
// the bottom half), and right of a node in the left half (and left of
// one in the right half) -- so it opens away from the nearest edge
// rather than off the edge of the window.
void drawfigure::drawprobe()
{
  if(!figure.has_probe)
    return ;

  positions nodeScreen = figtransform.apply_transform(figure.probe_pos) ;

  // Small circle marker, built directly in screen space (not through
  // the data-space FillPoly wrapper) so it stays a fixed, visible size
  // on screen regardless of zoom level -- an octagon reads as a circle
  // at this size.
  const double radiusPx = 4. ;
  const int NPTS = 8 ;
  positions circle[NPTS] ;
  for(int k=0;k<NPTS;++k) {
    double ang = 6.283185307179586*double(k)/double(NPTS) ;
    circle[k] = positions(nodeScreen.x+radiusPx*cos(ang),
                           nodeScreen.y+radiusPx*sin(ang)) ;
  }
  currdev->FillPolyMarker(circle,NPTS) ;

  // Python-style shortest-round-trip display, not fourSigDigs()'s
  // 4-significant-figure display used elsewhere (contour labels, the
  // title bar's min/max) -- this is for inspecting exact values, so
  // nothing should be rounded off, but it also shouldn't pad out clean
  // values with meaningless extra digits (see shortestReprDouble() in
  // tools.cc). probe_pos is a double throughout; probe_val is read into
  // a double but its true source precision is commonly a 32-bit float
  // from the upstream CFD solver, so shortestReprFloat() rounds through
  // an actual float first. Both return a shared static buffer, so each
  // call is fully consumed by snprintf before the next one.
  char lines[3][256] ;
  int nlines = 0 ;
  if(figure.probe_has_val) {
    snprintf(lines[nlines],256,"val: %s",shortestReprFloat(figure.probe_val)) ;
    ++nlines ;
  }
  snprintf(lines[nlines],256,"posx: %s",shortestReprDouble(figure.probe_pos.x)) ;
  ++nlines ;
  snprintf(lines[nlines],256,"posy: %s",shortestReprDouble(figure.probe_pos.y)) ;
  ++nlines ;

  // Anchor side is decided by which named edge of the drawable area the
  // node is closer to (rather than comparing against a numeric
  // midpoint), so it doesn't need to know whether this device's Y or X
  // axis numerically increases or decreases across the window.
  const DevBox &box = currdev->getDrawBox() ;
  const bool anchorBelow = fabs(nodeScreen.y-box.ytop) < fabs(nodeScreen.y-box.ybottom) ;
  const bool anchorRight = fabs(nodeScreen.x-box.xleft) < fabs(nodeScreen.x-box.xright) ;
  // +1/-1 multiplier for "toward the bottom edge" / "toward the right
  // edge" in this device's own coordinate direction, whichever way its
  // numbers happen to run.
  const double dirBottom = (box.ybottom>=box.ytop) ? 1.0 : -1.0 ;
  const double dirRight  = (box.xright>=box.xleft) ? 1.0 : -1.0 ;

  const double blockGapY = CHAR_HEIGHT_PX ;
  const double lineHeight = CHAR_HEIGHT_PX + 2. ;
  const double blockGapX = CHAR_WIDTH_PX*2. ;
  const double anchorX = nodeScreen.x + (anchorRight?1.0:-1.0)*dirRight*blockGapX ;

  SelectPen(0) ;
  for(int k=0;k<nlines;++k) {
    // When anchored below, line 0 sits nearest the node and later lines
    // move further away; when anchored above, that order is reversed so
    // the block still reads top-to-bottom (val/posx/posy) while its far
    // side, not its near side, touches the node.
    int rank = anchorBelow ? k : (nlines-1-k) ;
    double y = nodeScreen.y +
      (anchorBelow?1.0:-1.0)*dirBottom*(blockGapY+double(rank)*lineHeight) ;
    // XDevice::DrawLabel always draws left-justified regardless of the
    // `pos` argument (see the note in drawrulers()), so anchoring left
    // of the node needs the same manual right-justify shift used for
    // vertical ruler endpoint labels: move the anchor left by this
    // line's own width so its right edge, not its left edge, lands next
    // to the node.
    double lineX = anchorX ;
    if(!anchorRight)
      lineX -= CHAR_WIDTH_PX*double(strlen(lines[k])) ;
    currdev->DrawLabel(lines[k],positions(lineX,y),5,10,0) ;
  }
}

void drawfigure::draw()
{
  if(figure.show_shading)
    drawshading() ;
  drawborder() ;
  if(figure.show_grid)
    drawgrid() ;
  if(figure.show_contours)
    drawcontours() ;
  // Label the plot
  if(currdev->isLabelDev())
    drawtitle() ;
  if(figure.rulers.size() > 0)
    drawrulers() ;
  drawprobe() ;

  //  drawlegend() ;
}

void draw(const grid &drawgrid, grDevice &drawdev)
{
  drawfigure fig(drawgrid) ;
  drawdev.OpenDevice() ;
  fig.selectdev(drawdev) ;
  if(drawgrid.valid) {
    fig.scalefigure() ;
    fig.draw() ;
  } else {
    fig.drawtitle() ;
  }
  drawdev.CloseDevice() ;
}
