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

#ifndef LOCI_FVMADAPT_RECORD_IO_H
#define LOCI_FVMADAPT_RECORD_IO_H

#include <FVMAdapt/parametric/subdivision.h>

#include <cmath>
#include <cstring>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace Loci {
namespace FVMAdapt {
namespace detail {

  class RecordWriter ;
  class RecordReader ;
  /// Shared request encoding for histories and module input records.
  void writeSplitOperation(RecordWriter &out, const SplitOperation &operation) ;
  SplitOperation readSplitOperation(RecordReader &in, std::size_t limit) ;

  /// Write canonical records while enforcing byte and exact-arithmetic limits.
  class RecordWriter {
    std::vector<unsigned char> bytes_ ;
    std::size_t maxBytes_ ;
    ExactLimits exact_ ;
    ExactBudget budget_ ;
  public:
    RecordWriter(std::size_t maxBytes, ExactLimits exact) : maxBytes_(maxBytes), exact_(exact), budget_(exact) {}

    /// Append one byte within the record size limit.
    void byte(unsigned char value) {
      if(bytes_.size() >= maxBytes_) throw std::length_error("snapshot byte limit exceeded") ;
      bytes_.push_back(value) ;
    }

    /// Write an unsigned 32-bit value in little-endian order.
    void u32(std::uint32_t value) {
      for(unsigned i=0; i<4; ++i) { byte(value & 255) ; value >>= 8 ; }
    }

    /// Write an unsigned 64-bit value in little-endian order.
    void u64(std::uint64_t value) {
      for(unsigned i=0; i<8; ++i) { byte(value & 255) ; value >>= 8 ; }
    }

    /// Append encoded bytes within the record size limit.
    void append(const std::vector<unsigned char> &bytes) {
      if(bytes.size() > maxBytes_-bytes_.size()) throw std::length_error("snapshot byte limit exceeded") ;
      bytes_.insert(bytes_.end(),bytes.begin(),bytes.end()) ;
    }

    /// Write a finite binary64 value, normalizing negative zero.
    void real(double value) {
      static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
                    "snapshots require IEEE binary64") ;
      if(!std::isfinite(value)) throw std::invalid_argument("nonfinite snapshot coordinate") ;
      if(value==0) value=0 ;
      std::uint64_t bits ;
      std::memcpy(&bits,&value,8) ;
      u64(bits) ;
    }

    /// Write an integer sign and its little-endian magnitude.
    void integer(const Fraction::Integer &value) {
      const Fraction::Integer magnitude=value<0 ? -value : value ;
      const std::size_t bits=magnitude==0 ? 0 : boost::multiprecision::msb(magnitude)+1 ;
      if(bits>exact_.maxBits) throw std::length_error("snapshot integer bit limit exceeded") ;
      const std::size_t count=bits/8+(bits%8!=0) ;
      byte(value==0 ? 0 : value<0 ? 2 : 1) ;
      u64(count) ;
      if(count>maxBytes_-bytes_.size()) throw std::length_error("snapshot byte limit exceeded") ;
      if(count) boost::multiprecision::export_bits(magnitude,std::back_inserter(bytes_),8,false) ;
    }

    /// Charge construction work against the shared exact-arithmetic budget.
    void operations(std::size_t count) { budget_.operations(count) ; }

    /// Write a reduced fraction as its numerator and denominator.
    void fraction(const Fraction &value) {
      budget_.operations() ;
      checkExact(value,exact_) ;
      integer(value.numerator()) ; integer(value.denominator()) ;
    }

    const std::vector<unsigned char> &bytes() const { return bytes_ ; }

    /// Transfer ownership of the completed byte buffer to the caller.
    std::vector<unsigned char> take() { return std::move(bytes_) ; }
  } ;

  /// Read bounded, canonical records from a byte buffer.
  class RecordReader {
    const std::vector<unsigned char> &bytes_ ;
    std::size_t position_=0 ;
    ExactLimits exact_ ;
    ExactBudget budget_ ;
  public:
    RecordReader(const std::vector<unsigned char> &bytes, std::size_t maxBytes,
                 ExactLimits exact) : bytes_(bytes), exact_(exact), budget_(exact) {
      if(bytes.size()>maxBytes) throw std::length_error("snapshot byte limit exceeded") ;
    }

    std::size_t remaining() const { return bytes_.size()-position_ ; }

    /// Read one byte, rejecting truncated input.
    unsigned char byte() {
      if(!remaining()) throw std::invalid_argument("truncated snapshot") ;
      return bytes_[position_++] ;
    }

    /// Read an unsigned 32-bit value in little-endian order.
    std::uint32_t u32() {
      std::uint32_t result=0 ;
      for(unsigned i=0; i<4; ++i) result |= std::uint32_t(byte())<<(8*i) ;
      return result ;
    }

    /// Read an unsigned 64-bit value in little-endian order.
    std::uint64_t u64() {
      std::uint64_t result=0 ;
      for(unsigned i=0; i<8; ++i) result |= std::uint64_t(byte())<<(8*i) ;
      return result ;
    }

    /// Read a bounded count and check its minimum encoded payload size.
    std::size_t count(std::size_t limit, std::size_t minimumBytes=1) {
      const auto value=u64() ;
      if(value>limit || value>std::numeric_limits<std::size_t>::max())
        throw std::length_error("snapshot count limit exceeded") ;
      if(minimumBytes && value>remaining()/minimumBytes)
        throw std::invalid_argument("truncated snapshot array") ;
      return static_cast<std::size_t>(value) ;
    }

    /// Read a finite binary64 value and reject noncanonical negative zero.
    double real() {
      const auto bits=u64() ;
      double value ;
      std::memcpy(&value,&bits,8) ;
      if(!std::isfinite(value) || (value==0 && std::signbit(value)))
        throw std::invalid_argument("noncanonical snapshot coordinate") ;
      return value ;
    }

    /// Read a signed magnitude, rejecting redundant bytes and invalid signs.
    Fraction::Integer integer() {
      const auto sign=byte() ;
      if(sign>2) throw std::invalid_argument("unknown snapshot integer sign") ;
      const std::size_t maxBytes=exact_.maxBits/8+(exact_.maxBits%8!=0) ;
      const auto size=count(maxBytes) ;
      if((size==0)!=(sign==0) || (size && bytes_[position_+size-1]==0))
        throw std::invalid_argument("noncanonical snapshot integer") ;
      Fraction::Integer result=0 ;
      if(size) {
        boost::multiprecision::import_bits(result,bytes_.begin()+position_,bytes_.begin()+position_+size,8,false) ;
        position_+=size ;
        if(boost::multiprecision::msb(result)+1>exact_.maxBits)
          throw std::length_error("snapshot integer bit limit exceeded") ;
      }
      return sign==2 ? -result : result ;
    }

    /// Read a reduced fraction with a positive denominator.
    Fraction fraction() {
      budget_.operations() ;
      const auto numerator=integer(), denominator=integer() ;
      if(denominator<=0) throw std::invalid_argument("nonpositive snapshot denominator") ;
      const Fraction result(numerator,denominator) ;
      if(result.numerator()!=numerator || result.denominator()!=denominator)
        throw std::invalid_argument("unreduced snapshot fraction") ;
      return result ;
    }

    /// Require that the complete record has been consumed.
    void finish() const {
      if(remaining()) throw std::invalid_argument("trailing snapshot bytes") ;
    }
  } ;

}
}
}

#endif
