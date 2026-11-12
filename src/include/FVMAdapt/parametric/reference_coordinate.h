//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#ifndef LOCI_FVMADAPT_REFERENCE_COORDINATE_H
#define LOCI_FVMADAPT_REFERENCE_COORDINATE_H

#include "reference_element.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace Loci {
namespace FVMAdapt {

  /// Arbitrary-precision reduced fraction with a positive denominator.
  class Fraction {
  public:
    using Integer = boost::multiprecision::cpp_int ;
    using Rational = boost::multiprecision::cpp_rational ;
  private:
    Rational value_ ;
  public:
    /// A zero denominator throws invalid_argument.
    Fraction(std::int64_t numerator = 0, std::int64_t denominator = 1) ;
    Fraction(const Integer &numerator, const Integer &denominator) ;
    explicit Fraction(const Rational &value) ;
    Integer numerator() const { return boost::multiprecision::numerator(value_) ; }
    Integer denominator() const { return boost::multiprecision::denominator(value_) ; }
    const Rational &rational() const { return value_ ; }
    double value() const { return value_.convert_to<double>() ; }
  } ;

  /// Convert an exact integer without truncation; out-of-range values throw overflow_error.
  template<class Integer> Integer checkedExactInteger(const Fraction::Integer &value) {
    static_assert(std::is_integral<Integer>::value,"an integral destination is required") ;
    if(value < std::numeric_limits<Integer>::min() || value > std::numeric_limits<Integer>::max())
      throw std::overflow_error("exact integer exceeds the destination range") ;
    return value.convert_to<Integer>() ;
  }

  /// Caller-selected bounds for exact constructions and retained records.
  struct ExactLimits {
    std::size_t maxBits = 4096 ;
    std::size_t maxOperations = 1000000 ;
    std::size_t maxBytes = 64*1024*1024 ;
  } ;

  /// Reject numerator or denominator magnitudes exceeding maxBits with length_error.
  void checkExact(const Fraction &value, const ExactLimits &limits) ;

  /// Local work and storage accounting; no arithmetic state is shared between callers.
  class ExactBudget {
    ExactLimits limits_ ;
    std::size_t operations_ = 0, bytes_ = 0 ;
  public:
    explicit ExactBudget(ExactLimits limits = {}) : limits_(limits) {}
    const ExactLimits &limits() const { return limits_ ; }
    void check(const Fraction &value) const { checkExact(value,limits_) ; }
    /// Charge work or retained bytes; exhaustion throws length_error without changing the counter.
    void operations(std::size_t count = 1) ;
    void retain(std::size_t bytes) ;
  } ;

  Fraction operator+(const Fraction &a, const Fraction &b) ;
  Fraction operator-(const Fraction &a, const Fraction &b) ;
  Fraction operator-(const Fraction &a) ;
  Fraction operator*(const Fraction &a, const Fraction &b) ;
  /// Division by zero throws domain_error.
  Fraction operator/(const Fraction &a, const Fraction &b) ;
  bool operator==(const Fraction &a, const Fraction &b) ;
  bool operator!=(const Fraction &a, const Fraction &b) ;
  bool operator<(const Fraction &a, const Fraction &b) ;
  bool operator>(const Fraction &a, const Fraction &b) ;
  bool operator<=(const Fraction &a, const Fraction &b) ;
  bool operator>=(const Fraction &a, const Fraction &b) ;

  /// Exact logical coordinates; physical map evaluation uses approximate(point).
  typedef std::array<Fraction,3> ReferencePoint ;

  Point approximate(const ReferencePoint &point) ;
  /// Return a native vertex; invalid indices throw out_of_range.
  ReferencePoint referenceVertex(Shape shape, unsigned vertex) ;
  /// Exact domain membership, including zero unused components. Unknown shapes throw invalid_argument.
  bool containsExact(Shape shape, const ReferencePoint &point) ;
  /// Test boundary membership without constructing boundary-local coordinates.
  /// Outside points return false; invalid shapes or boundary indices throw.
  bool onBoundary(Shape shape, unsigned boundary, const ReferencePoint &point) ;

  /// Embed a point from a boundary into its parent reference coordinates.
  /// Points outside the boundary throw invalid_argument; invalid indices throw out_of_range.
  ReferencePoint boundaryReferencePoint(Shape shape, unsigned boundary,
                                       const ReferencePoint &point) ;
  /// Convert a point on the chosen boundary to that boundary's local coordinates.
  /// Points outside the parent or off the boundary throw invalid_argument.
  ReferencePoint boundaryCoordinates(Shape shape, unsigned boundary,
                                     const ReferencePoint &point) ;
  /// Apply the same vertex permutation as orientPoint; points outside the shape throw invalid_argument.
  ReferencePoint orientReferencePoint(Shape shape, unsigned orientation,
                                     const ReferencePoint &point) ;

}
}

#endif
