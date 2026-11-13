// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_MESH_QUALITY_H
#define LOCI_FVMADAPT_MESH_QUALITY_H

#include "geometry_map.h"
#include <string>

namespace Loci {
namespace FVMAdapt {
  /// Optional limits on the emitted Cartesian mesh. A limit of -1 is disabled.
  struct MeshQualityOptions {
    std::string mode = "off" ; ///< off, report, or reject.
    std::string centroid = "wireframe" ; ///< FVM metric convention: wireframe or exact.
    double maxNonOrthogonality = -1 ; ///< Degrees, including physical boundary faces.
    double maxFaceWarpage = -1 ; ///< Maximum triangle-normal deviation, in degrees.
    double maxChildVolumeRatio = -1 ; ///< Largest/smallest immediate child, summing active descendants.
    void validate() const ;
  } ;

  /// Polygon measurements and tetrahedral first moments; no mesh ownership or connectivity.
  struct FaceMeasure {
    Point center{} ; ///< Face center using the selected FVM convention.
    Point area{} ; ///< Oriented area vector following the vertex loop.
    Point wireCenter{} ; ///< Edge-length-weighted mean of edge midpoints.
    /// moment[3*d+j] sums (a[d]+b[d])*cross(a,b)[j]/24; a,b are relative to wireCenter.
    double moment[9]{} ;
    double warpage = 0 ; ///< Largest triangle-normal deviation from area, in degrees.
  } ;
  /// Finite-volume cell center and signed boundary volume.
  struct CellMeasure {
    Point center{} ;
    double volume = 0 ;
  } ;
  /// Oriented use of a measured face; +1 is outward, -1 is inward.
  struct MeasuredFace { const FaceMeasure *face ; int sign ; } ;
  /// Measure the output polygon using the solver's wireframe or exact center convention.
  /// These centers do not define split points or integrate a retained curved surface.
  FaceMeasure measurePolygon(const std::vector<Point> &, const std::string &centroid) ;
  /// Measure a closed cell from its oriented faces using the same center convention.
  CellMeasure measurePolyhedron(const std::vector<MeasuredFace> &, const std::string &centroid) ;
  /// Angle between an oriented area vector and the left-to-right center displacement.
  double faceNonOrthogonality(const FaceMeasure &, const Point &left, const Point &right) ;

  /// Original file cell ID and current root-local leaf index; face is an original face ID, or -1 internally.
  struct QualityValue {
    double value = 0 ;
    int root = -1, leaf = -1, face = -1 ;
  } ;
  /// Bounded summary; invalid geometry is distinct from a user-limit rejection.
  struct MeshQualityReport {
    QualityValue nonOrthogonality, faceWarpage, childVolumeRatio ;
    int invalidRoot = -1, invalidLeaf = -1, invalidFace = -1 ;
    bool invalid() const { return invalidRoot >= 0 ; }
    bool exceeds(const MeshQualityOptions &) const ;
    std::string description(const MeshQualityOptions &) const ;
    void include(const MeshQualityReport &) ;
  } ;
  void includeQualityValue(QualityValue &, double value, int root, int leaf, int face = -1) ;
  /// Compare immediate active child volumes at every split. Each leaf path is its child-ordinal history.
  QualityValue childVolumeRatio(const std::vector<std::vector<unsigned>> &paths, const std::vector<CellMeasure> &, int root) ;
}
}
#endif
