// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_XML_REGION_H
#define LOCI_FVMADAPT_XML_REGION_H

#include "parametric/reference_element.h"
#include <functional>
#include <memory>
#include <string>

struct _xmlNode ;

namespace Loci {
namespace FVMAdapt {
  /// Immutable XML region evaluated in physical coordinates.
  class XMLRegion {
    std::function<bool(const Point &)> contains_ ;
  public:
    explicit XMLRegion(const std::string &filename) ;
    explicit XMLRegion(_xmlNode *root) ;
    bool contains(const Point &point) const { return contains_(point) ; }
  } ;
  using XMLRegionP=std::shared_ptr<const XMLRegion> ;
}
}
#endif
