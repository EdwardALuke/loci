// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_GEOMETRY_BOUNDS_H
#define LOCI_FVMADAPT_GEOMETRY_BOUNDS_H
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace Loci {
namespace FVMAdapt {
namespace detail {
  /// Outward-rounded bounds; input doubles denote exact values.
  struct Bounds {
    double low = 0, high = 0 ;
    Bounds() = default ;
    Bounds(double value) : low(value), high(value) {}
    Bounds(double low, double high)
      : low(std::nextafter(low,-std::numeric_limits<double>::infinity())),
        high(std::nextafter(high,std::numeric_limits<double>::infinity())) {}
  } ;
  inline Bounds operator+(Bounds a, Bounds b) { return {a.low+b.low,a.high+b.high} ; }
  inline Bounds operator-(Bounds a, Bounds b) { return {a.low-b.high,a.high-b.low} ; }
  inline Bounds operator*(Bounds a, Bounds b) {
    const std::array<double,4> products{{a.low*b.low,a.low*b.high,a.high*b.low,a.high*b.high}} ;
    for(double value: products) if(!std::isfinite(value))
      return {-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()} ;
    const auto range = std::minmax_element(products.begin(),products.end()) ; return {*range.first,*range.second} ;
  }
  inline Bounds operator/(Bounds a, double positive) { return {a.low/positive,a.high/positive} ; }
}
}
}
#endif
