#include "tools.h"
#include <stdlib.h>
#include <stdio.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>

char *tilde(const char *file)
{
    struct passwd *pass ;
    char buf[512] ;
    const char *remains ;
    char *newfile ;
    
    if(*file != '~') {
        newfile = (char *) malloc(strlen(file)+1) ;
        strcpy(newfile,file) ;
        return(newfile) ;
    }    
    file++ ;
    if(*file == '/') {
        pass = getpwuid(getuid()) ;
        if(pass == NULL) {
            fprintf(stderr,
                    "warning: Can't get password entry for your userid!\n") ;
            --file ;
            newfile = (char *) malloc(strlen(file)+1) ;
            strcpy(newfile,file) ;
            return newfile ;
        }
    } else {
        char *t ;
        
        remains = file - 1 ;
        t = buf ;
        while(*file != '/' && *file != '\0')
          *t++ = *file++ ;
        *t = '\0' ;
        pass = getpwnam(buf) ;
        if(pass == NULL) {
            fprintf(stderr,
                  "warning: Can't get password entry for userid '%s'\n",buf) ;
            newfile = (char *) malloc(strlen(remains)+1) ;
            strcpy(newfile,remains) ;
            return newfile ;
        }
    }
    newfile = (char *)malloc(strlen(file)+strlen(pass->pw_dir)+1) ;
    strcpy(newfile,pass->pw_dir) ;
    strcat(newfile,file) ;
    return(newfile) ;
}

// Prints num the way Python's str()/repr() prints a float: the shortest
// decimal string that reads back to the exact same IEEE 754 double --
// e.g. 0.1 prints as "0.1", not the much longer digit string %.17g gives
// for the same bit pattern, and 2.0 prints as "2.0", not "2". Finding
// the minimal round-tripping digit count and choosing fixed vs.
// scientific notation are done as two separate steps: doing it in one
// step by trying "%.*g" with increasing precision does not work, because
// %g's own notation-switching threshold is tied to that same precision
// value, so a low precision (found while searching for the shortest
// digit count) can trigger scientific notation even for ordinary
// numbers like 100 ("1e+02" instead of "100"). Using "%.*e" for the
// round-trip search avoids that coupling, since %e is always scientific
// regardless of precision; the exponent it reports is then compared
// against the same threshold Python's repr() uses (fixed notation for
// exponents in [-4,16), scientific outside that range, verified
// empirically against Python) to decide how to render the digits found.
// This is a slower way to reach shortest-round-trip formatting than the
// dedicated algorithms built for it (Grisu/Ryu, or C++17's
// std::to_chars), but that's fine here since it only runs once per probe
// click, not in a hot loop.
char * shortestReprDouble(double num)
{
    static char buf[64] ;
    if(num == 0.0) {
      snprintf(buf,64,signbit(num)?"-0.0":"0.0") ;
      return buf ;
    }

    char ebuf[64] ;
    int sig = 0 ;
    for(sig=0; sig<=16; ++sig) {
      snprintf(ebuf,64,"%.*e",sig,num) ;
      if(strtod(ebuf,NULL) == num)
        break ;
    }

    char *epos = strchr(ebuf,'e') ;
    int exponent = epos ? atoi(epos+1) : 0 ;

    if(exponent >= -4 && exponent < 16) {
      int decimals = sig-exponent ;
      if(decimals < 1)
        decimals = 1 ; // Python always shows at least one decimal digit
      snprintf(buf,64,"%.*f",decimals,num) ;
    } else {
      snprintf(buf,64,"%s",ebuf) ;
    }
    return buf ;
}

// Same idea as shortestReprDouble(), but for a value whose true source
// precision is a 32-bit float even though it's stored here as a double
// (grid::val is read from the data file into a double, but solution
// fields from the upstream CFD solver are commonly single precision).
// Rounding through an actual float first, then finding the shortest
// decimal that round-trips to THAT float, shows exactly what the
// original float held -- nothing invented by widening it to double, and
// nothing longer than necessary.
char * shortestReprFloat(double num)
{
    static char buf[64] ;
    float f = (float)num ;
    if(f == 0.0f) {
      snprintf(buf,64,signbit(f)?"-0.0":"0.0") ;
      return buf ;
    }

    char ebuf[64] ;
    int sig = 0 ;
    for(sig=0; sig<=8; ++sig) {
      snprintf(ebuf,64,"%.*e",sig,(double)f) ;
      if((float)strtod(ebuf,NULL) == f)
        break ;
    }

    char *epos = strchr(ebuf,'e') ;
    int exponent = epos ? atoi(epos+1) : 0 ;

    if(exponent >= -4 && exponent < 16) {
      int decimals = sig-exponent ;
      if(decimals < 1)
        decimals = 1 ;
      snprintf(buf,64,"%.*f",decimals,(double)f) ;
    } else {
      snprintf(buf,64,"%s",ebuf) ;
    }
    return buf ;
}

char * fourSigDigs(double num)
{
    static char buf[256] ;
    
    if(num == 0.0)
      snprintf(buf,256,"0.0") ;
    else if(fabs(num) < 1.0|| fabs(num)> 99999.9) {
        if(fabs(num) > 0.009 && fabs(num)<1.0)
          snprintf(buf,256,"%3.3f",num) ;
        else
          snprintf(buf,256,"%3.2e",num) ;
    } else if(fabs(num) > 999.9)
      snprintf(buf,256,"%g",num) ;
    else
      snprintf(buf,256,"%3.1f",num) ;
    return buf ;

}

// Format a ruler tick label with enough decimal places to resolve the
// current tick spacing (increment). Uses scientific notation for very
// large or very small magnitudes (same thresholds as fourSigDigs(), but
// decided per value so a ruler isn't forced into one notation across its
// whole range). EXTRA_DIGITS adds precision beyond what's strictly
// needed to distinguish adjacent ticks.
char * rulerDigits(double num, double increment)
{
    static char buf[256] ;
    const int EXTRA_DIGITS = 1 ; // bump this up for even more precision

    if(increment == 0.0 || !isfinite(increment))
      increment = 1.0 ;
    increment = fabs(increment) ;

    int decimals = -(int)floor(log10(increment)) + EXTRA_DIGITS ;
    if(decimals < 0)
      decimals = 0 ;
    if(decimals > 15)
      decimals = 15 ; // a double has about 15-17 reliable significant
                       // decimal digits; beyond this, requested precision
                       // isn't real -- draw.cc's minimum-increment clamp
                       // (see drawrulers()) keeps the increment from ever
                       // asking for more than this in practice

    if(num == 0.0) {
      snprintf(buf,256,"0") ;
    } else if(fabs(num) > 99999.9 || fabs(num) < 0.009) {
      // Scientific notation mantissa uses a small fixed precision --
      // the exponent already encodes magnitude, so it shouldn't grow
      // just because the increment is very small or very large.
      int sig = 2 + EXTRA_DIGITS ;
      snprintf(buf,256,"%.*e",sig,num) ;
    } else {
      snprintf(buf,256,"%.*f",decimals,num) ;
    }

    return buf ;
}
