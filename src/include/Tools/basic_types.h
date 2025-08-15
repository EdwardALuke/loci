//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
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
#ifndef LOCI_BASIC_TYPES_H
#define LOCI_BASIC_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#ifndef NO_AUTODIFF
#include <Tools/autodiff.h>
#endif

#include <Tools/debug.h>

namespace Loci {

#ifdef USE_AUTODIFF
#ifdef AUTODIFF2ND
  typedef Loci::FAD2d real_t ;
#else
#ifdef MULTIFAD
  typedef Loci::MFADd real_t ;
#else
  typedef Loci::FADd real_t ;
#endif
#endif
#else
  // No autodiff
  typedef double real_t ;
#endif

  //---------------------------------------------------------------------------
  // Expression Template Helper Classes


  // copy version of unary operator
  template<typename T, typename Op >
  struct ArrayUnaryC{
    const Op op1;
    ArrayUnaryC(const Op& a): op1(a) {}
    T operator[](const size_t i) const { return -op1[i] ;  }
  } ;
  // base version of unary operator
  template <typename T, typename Op>
  struct ArrayUnaryB{
    const Op &op1;
    ArrayUnaryB(const Op& a): op1(a) {}
    T operator[](const size_t i) const { return -op1[i] ;  }
  } ;

  template <typename T, typename Op>
  struct ArrayOperator {
    const Op op ;
    ArrayOperator(const Op &t):op(t) {}
    T operator[](const size_t i) const { return op[i]; }
    const Op &data() const { return op; }
    ArrayOperator<T,ArrayUnaryC<T,Op> > operator-() const { return ArrayOperator<T,ArrayUnaryC<T,Op> >(op) ; }
  } ;

  //---------------------------------------------------------------------------
  // Addition Operators
  //
  // Case when both left and right are concrete (not anonymous)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayAdd{
    const Op1& op1;
    const Op2& op2;
    ArrayAdd(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { return op1[i] + op2[i] ;  }
  } ;
  // Case where left argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayAddL{
    const Op1 op1;
    const Op2& op2;
    ArrayAddL(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] + op2[i] ;  }
  } ;

  // case where right argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayAddR{
    const Op1 &op1;
    const Op2 op2;
    ArrayAddR(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] + op2[i] ;  }
  } ;

  template<typename T, typename Op1 , typename Op2>
  struct ArrayAddLR{
    const Op1 op1;
    const Op2 op2;
    ArrayAddLR(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { return op1[i] + op2[i] ;  }
  } ;
  
  // Case where left side is a double constant
  template<typename T, typename Op2>
  struct ArrayAddLC {
    const double op1;
    const Op2& op2;
    ArrayAddLC(const double& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1 + op2[i] ;  }
  } ;
  // Case where right side is a double constant
  template<typename T, typename Op1 >
  struct ArrayAddRC{
    const Op1 &op1;
    const double op2;
    ArrayAddRC(const Op1& a, const double &b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] + op2 ;  }
  } ;

  // General operators for packaged expressions
  // function template for the + operator
  template<typename T, typename R1, typename R2>
  inline ArrayOperator<T, ArrayAddLR<T, R1, R2> >
  operator+ (const ArrayOperator<T, R1>& a, const ArrayOperator<T, R2>& b){
    return ArrayOperator<T, ArrayAddLR<T, R1, R2> >(ArrayAddLR<T, R1, R2 >(a.data(), b.data()));
  }
  // case where expression is added to a constant on the right
  // function template for the + operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayAddRC<T, R1> >
  operator+ (const ArrayOperator<T, R1>& a, const double& b){
    return ArrayOperator<T, ArrayAddRC<T, R1> >(ArrayAddRC<T, R1 >(a.data(), b));
  }
  // case where constant on the left is added to an expression
  // function template for the + operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayAddLC<T, R1> >
  operator+ (const double &a, const ArrayOperator<T, R1>& b){
    return ArrayOperator<T, ArrayAddLC<T, R1> >(ArrayAddLC<T, R1 >(a, b.data()));
  }


  //---------------------------------------------------------------------------
  // Subtraction Operators
  //
  // Case when both left and right are concrete (not anonymous)
  template<typename T, typename Op1 , typename Op2>
  struct ArraySub{
    const Op1& op1;
    const Op2& op2;
    ArraySub(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { return op1[i] - op2[i] ;  }
  } ;
  // Case where left argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArraySubL{
    const Op1 op1;
    const Op2& op2;
    ArraySubL(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] - op2[i] ;  }
  } ;

  // case where right argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArraySubR{
    const Op1 &op1;
    const Op2 op2;
    ArraySubR(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] - op2[i] ;  }
  } ;

  template<typename T, typename Op1 , typename Op2>
  struct ArraySubLR{
    const Op1 op1;
    const Op2 op2;
    ArraySubLR(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { return op1[i] - op2[i] ;  }
  } ;
  // Case where left side is a double constant
  template<typename T, typename Op2>
  struct ArraySubLC {
    const double op1;
    const Op2& op2;
    ArraySubLC(const double& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1 - op2[i] ;  }
  } ;
  // Case where right side is a double constant
  template<typename T, typename Op1 >
  struct ArraySubRC{
    const Op1 &op1;
    const double op2;
    ArraySubRC(const Op1& a, const double &b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] - op2 ;  }
  } ;

  // General operators for packaged expressions
  // function template for the - operator
  template<typename T, typename R1, typename R2>
  inline ArrayOperator<T, ArraySubLR<T, R1, R2> >
  operator- (const ArrayOperator<T, R1>& a, const ArrayOperator<T, R2>& b){
    return ArrayOperator<T, ArraySubLR<T, R1, R2> >(ArraySubLR<T, R1, R2 >(a.data(), b.data()));
  }
  // case where expression is added to a constant on the right
  // function template for the - operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArraySubRC<T, R1> >
  operator- (const ArrayOperator<T, R1>& a, const double& b){
    return ArrayOperator<T, ArraySubRC<T, R1> >(ArraySubRC<T, R1 >(a.data(), b));
  }
  // case where constant on the left is added to an expression
  // function template for the - operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArraySubLC<T, R1> >
  operator- (const double &a, const ArrayOperator<T, R1>& b){
    return ArrayOperator<T, ArraySubLC<T, R1> >(ArraySubLC<T, R1 >(a, b.data()));
  }

  //---------------------------------------------------------------------------
  // Multiplication Operators
  //
  // Case when both left and right are concrete (not anonymous)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayMul{
    const Op1& op1;
    const Op2& op2;
    ArrayMul(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { return op1[i] * op2[i] ;  }
  } ;
  // Case where left argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayMulL{
    const Op1 op1;
    const Op2& op2;
    ArrayMulL(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] * op2[i] ;  }
  } ;

  // case where right argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayMulR{
    const Op1 &op1;
    const Op2 op2;
    ArrayMulR(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] * op2[i] ;  }
  } ;

  // case where right argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayMulLR{
    const Op1 op1;
    const Op2 op2;
    ArrayMulLR(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] * op2[i] ;  }
  } ;

  // Case where left side is a double constant
  template<typename T, typename Op2>
  struct ArrayMulLC {
    const double op1;
    const Op2& op2;
    ArrayMulLC(const double& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1 * op2[i] ;  }
  } ;
  // Case where right side is a double constant
  template<typename T, typename Op1 >
  struct ArrayMulRC{
    const Op1 &op1;
    const double op2;
    ArrayMulRC(const Op1& a, const double &b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] * op2 ;  }
  } ;

  // General operators for packaged expressions
  // function template for the * operator
  template<typename T, typename R1, typename R2>
  inline ArrayOperator<T, ArrayMulLR<T, R1, R2> >
  operator* (const ArrayOperator<T, R1>& a, const ArrayOperator<T, R2>& b){
    return ArrayOperator<T, ArrayMulLR<T, R1, R2> >(ArrayMulLR<T, R1, R2 >(a.data(), b.data()));
  }
  // case where expression is added to a constant on the right
  // function template for the * operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayMulRC<T, R1> >
  operator* (const ArrayOperator<T, R1>& a, const double& b){
    return ArrayOperator<T, ArrayMulRC<T, R1> >(ArrayMulRC<T, R1 >(a.data(), b));
  }
  // case where constant on the left is added to an expression
  // function template for the * operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayMulLC<T, R1> >
  operator* (const double &a, const ArrayOperator<T, R1>& b){
    return ArrayOperator<T, ArrayMulLC<T, R1> >(ArrayMulLC<T, R1 >(a, b.data()));
  }

  //---------------------------------------------------------------------------
  // Division Operators
  //
  // Case when both left and right are concrete (not anonymous)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayDiv{
    const Op1& op1;
    const Op2& op2;
    ArrayDiv(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { return op1[i] / op2[i] ;  }
  } ;
  // Case where left argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayDivL{
    const Op1 op1;
    const Op2& op2;
    ArrayDivL(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] / op2[i] ;  }
  } ;

  // case where right argument of add is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayDivR{
    const Op1 &op1;
    const Op2 op2;
    ArrayDivR(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] / op2[i] ;  }
  } ;
  template<typename T, typename Op1 , typename Op2>
  struct ArrayDivLR{
    const Op1 op1;
    const Op2 op2;
    ArrayDivLR(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { return op1[i] / op2[i] ;  }
  } ;

  // Case where left side is a double constant
  template<typename T, typename Op2>
  struct ArrayDivLC {
    const double op1;
    const Op2& op2;
    ArrayDivLC(const double& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1 / op2[i] ;  }
  } ;
  // Case where right side is a double constant
  template<typename T, typename Op1 >
  struct ArrayDivRC{
    const Op1 &op1;
    const double op2;
    ArrayDivRC(const Op1& a, const double &b): op1(a), op2(b){}
    T operator[](const size_t i) const { return op1[i] / op2 ;  }
  } ;

  // General operators for packaged expressions
  // function template for the / operator
  template<typename T, typename R1, typename R2>
  inline ArrayOperator<T, ArrayDivLR<T, R1, R2> >
  operator/ (const ArrayOperator<T, R1>& a, const ArrayOperator<T, R2>& b){
    return ArrayOperator<T, ArrayDivLR<T, R1, R2> >(ArrayDivLR<T, R1, R2 >(a.data(), b.data()));
  }
  // case where expression is added to a constant on the right
  // function template for the / operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayDivRC<T, R1> >
  operator/ (const ArrayOperator<T, R1>& a, const double& b){
    return ArrayOperator<T, ArrayDivRC<T, R1> >(ArrayDivRC<T, R1 >(a.data(), b));
  }
  // case where constant on the left is added to an expression
  // function template for the / operator
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayDivLC<T, R1> >
  operator/ (const double &a, const ArrayOperator<T, R1>& b){
    return ArrayOperator<T, ArrayDivLC<T, R1> >(ArrayDivLC<T, R1 >(a, b.data()));
  }

  //---------------------------------------------------------------------------
  // Trigonometric functions
  //-----
  // cos
  template<typename T, typename Op1 >
  struct ArraycosB{
    const Op1 &op1;
    ArraycosB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::cos ; return cos(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraycosC{
    const Op1 op1;
    ArraycosC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::cos ; return cos(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraycosC<T, R1> >
  cos (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraycosC<T, R1> >(ArraycosC<T, R1 >(a.data()));
  }
  //-----
  // sin
  template<typename T, typename Op1 >
  struct ArraysinB{
    const Op1 &op1;
    ArraysinB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::sin; return sin(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraysinC{
    const Op1 op1;
    ArraysinC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::sin ; return sin(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraysinC<T, R1> >
  sin (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraysinC<T, R1> >(ArraysinC<T, R1 >(a.data()));
  }
  //-----
  // tan
  template<typename T, typename Op1 >
  struct ArraytanB{
    const Op1 &op1;
    ArraytanB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::tan ; return tan(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraytanC{
    const Op1 op1;
    ArraytanC(const Op1& a): op1(a){}
    T operator[](const size_t i) const { return tan(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraytanC<T, R1> >
  tan (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraytanC<T, R1> >(ArraytanC<T, R1 >(a.data()));
  }
  //-----
  // acos
  template<typename T, typename Op1 >
  struct ArrayacosB{
    const Op1 &op1;
    ArrayacosB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::acos ; return acos(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayacosC{
    const Op1 op1;
    ArrayacosC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::acos; return acos(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayacosC<T, R1> >
  acos (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayacosC<T, R1> >(ArrayacosC<T, R1 >(a.data()));
  }
  //-----
  // asin
  template<typename T, typename Op1 >
  struct ArrayasinB{
    const Op1 &op1;
    ArrayasinB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::asin; return asin(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayasinC{
    const Op1 op1;
    ArrayasinC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::asin ; return asin(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayasinC<T, R1> >
  asin (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayasinC<T, R1> >(ArrayasinC<T, R1 >(a.data()));
  }
  //-----
  // atan
  template<typename T, typename Op1 >
  struct ArrayatanB{
    const Op1 &op1;
    ArrayatanB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::atan ; return atan(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayatanC{
    const Op1 op1;
    ArrayatanC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::atan; return atan(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayatanC<T, R1> >
  atan (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayatanC<T, R1> >(ArrayatanC<T, R1 >(a.data()));
  }

  //---------------------------------------------------------------------------
  // Hyperbolic functions
  //-----
  // cosh
  template<typename T, typename Op1 >
  struct ArraycoshB{
    const Op1 &op1;
    ArraycoshB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::cosh ; return cosh(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraycoshC{
    const Op1 op1;
    ArraycoshC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::cosh; return cosh(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraycoshC<T, R1> >
  cosh (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraycoshC<T, R1> >(ArraycoshC<T, R1 >(a.data()));
  }
  //-----
  // sinh
  template<typename T, typename Op1 >
  struct ArraysinhB{
    const Op1 &op1;
    ArraysinhB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::sinh; return sinh(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraysinhC{
    const Op1 op1;
    ArraysinhC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::sinh; return sinh(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraysinhC<T, R1> >
  sinh (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraysinhC<T, R1> >(ArraysinhC<T, R1 >(a.data()));
  }
  //-----
  // tanh
  template<typename T, typename Op1 >
  struct ArraytanhB{
    const Op1 &op1;
    ArraytanhB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::tanh ; return tanh(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraytanhC{
    const Op1 op1;
    ArraytanhC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::tanh; return tanh(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraytanhC<T, R1> >
  tanh (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraytanhC<T, R1> >(ArraytanhC<T, R1 >(a.data()));
  }
  //-----
  // acosh
  template<typename T, typename Op1 >
  struct ArrayacoshB{
    const Op1 &op1;
    ArrayacoshB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::acosh ; return acosh(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayacoshC{
    const Op1 op1;
    ArrayacoshC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::acosh ; return acosh(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayacoshC<T, R1> >
  acosh (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayacoshC<T, R1> >(ArrayacoshC<T, R1 >(a.data()));
  }
  //-----
  // asinh
  template<typename T, typename Op1 >
  struct ArrayasinhB{
    const Op1 &op1;
    ArrayasinhB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::asinh ; return asinh(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayasinhC{
    const Op1 op1;
    ArrayasinhC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::asinh ; return asinh(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayasinhC<T, R1> >
  asinh (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayasinhC<T, R1> >(ArrayasinhC<T, R1 >(a.data()));
  }
  //-----
  // atanh
  template<typename T, typename Op1 >
  struct ArrayatanhB{
    const Op1 &op1;
    ArrayatanhB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::atanh ;return atanh(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayatanhC{
    const Op1 op1;
    ArrayatanhC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::atanh ; return atanh(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayatanhC<T, R1> >
  atanh (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayatanhC<T, R1> >(ArrayatanhC<T, R1 >(a.data()));
  }
  
  //---------------------------------------------------------------------------
  // exponential and logrithmic functions
  //
  //----
  // exp
  template<typename T, typename Op1 >
  struct ArrayexpB{
    const Op1 &op1;
    ArrayexpB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::exp ; return exp(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayexpC{
    const Op1 op1;
    ArrayexpC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::exp;  return exp(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayexpC<T, R1> >
  exp (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayexpC<T, R1> >(ArrayexpC<T, R1 >(a.data()));
  }
  //----
  // log
  template<typename T, typename Op1 >
  struct ArraylogB{
    const Op1 &op1;
    ArraylogB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::log ; return log(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraylogC{
    const Op1 op1;
    ArraylogC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::log ; return log(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraylogC<T, R1> >
  log (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraylogC<T, R1> >(ArraylogC<T, R1 >(a.data()));
  }
  //----
  // log10
  template<typename T, typename Op1 >
  struct Arraylog10B{
    const Op1 &op1;
    Arraylog10B(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::log10; return log10(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct Arraylog10C{
    const Op1 op1;
    Arraylog10C(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::log10 ; return log10(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, Arraylog10C<T, R1> >
  log10 (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, Arraylog10C<T, R1> >(Arraylog10C<T, R1 >(a.data()));
  }

  //---------------------------------------------------------------------------
  // error and gamma functions
  //
  //-----
  // erf
  template<typename T, typename Op1 >
  struct ArrayerfB{
    const Op1 &op1;
    ArrayerfB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::erf ; return erf(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayerfC{
    const Op1 op1;
    ArrayerfC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::erf ; return erf(op1[i]) ;  }
  } ;
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayerfC<T, R1> >
  erf (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayerfC<T, R1> >(ArrayerfC<T, R1 >(a.data()));
  }
  //-----
  // erfc
  template<typename T, typename Op1 >
  struct ArrayerfcB{
    const Op1 &op1;
    ArrayerfcB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::erfc ; return erfc(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayerfcC{
    const Op1 op1;
    ArrayerfcC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::erfc ; return erfc(op1[i]) ;  }
  } ;
  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayerfcC<T, R1> >
  erfc (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayerfcC<T, R1> >(ArrayerfcC<T, R1 >(a.data()));
  }
  //-----
  // tgamma
  template<typename T, typename Op1 >
  struct ArraytgammaB{
    const Op1 &op1;
    ArraytgammaB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::tgamma ; return tgamma(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraytgammaC{
    const Op1 op1;
    ArraytgammaC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::tgamma ; return tgamma(op1[i]) ;  }
  } ;
  template<typename T, typename R1>
  inline ArrayOperator<T, ArraytgammaC<T, R1> >
  tgamma (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraytgammaC<T, R1> >(ArraytgammaC<T, R1 >(a.data()));
  }
  //-----
  // lgamma
  template<typename T, typename Op1 >
  struct ArraylgammaB{
    const Op1 &op1;
    ArraylgammaB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::lgamma ; return lgamma(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraylgammaC{
    const Op1 op1;
    ArraylgammaC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::lgamma ; return lgamma(op1[i]) ;  }
  } ;
  template<typename T, typename R1>
  inline ArrayOperator<T, ArraylgammaC<T, R1> >
  lgamma (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraylgammaC<T, R1> >(ArraylgammaC<T, R1 >(a.data()));
  }

  //---------------------------------------------------------------------------
  // misc functions, ceil, floor,trunc,round,fabs,abs
  //
  //----
  // ceil
  template<typename T, typename Op1 >
  struct ArrayceilB{
    const Op1 &op1;
    ArrayceilB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::ceil ; return ceil(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayceilC{
    const Op1 op1;
    ArrayceilC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::ceil ; return ceil(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayceilC<T, R1> >
  ceil (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayceilC<T, R1> >(ArrayceilC<T, R1 >(a.data()));
  }
  //----
  // floor
  template<typename T, typename Op1 >
  struct ArrayfloorB{
    const Op1 &op1;
    ArrayfloorB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::floor ; return floor(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayfloorC{
    const Op1 op1;
    ArrayfloorC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::floor ; return floor(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayfloorC<T, R1> >
  floor (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayfloorC<T, R1> >(ArrayfloorC<T, R1 >(a.data()));
  }
  //----
  // trunc
  template<typename T, typename Op1 >
  struct ArraytruncB{
    const Op1 &op1;
    ArraytruncB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { return trunc(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraytruncC{
    const Op1 op1;
    ArraytruncC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { return trunc(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraytruncC<T, R1> >
  trunc (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraytruncC<T, R1> >(ArraytruncC<T, R1 >(a.data()));
  }
  //----
  // round
  template<typename T, typename Op1 >
  struct ArrayroundB{
    const Op1 &op1;
    ArrayroundB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::round ; return round(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayroundC{
    const Op1 op1;
    ArrayroundC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::round ; return round(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayroundC<T, R1> >
  round (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayroundC<T, R1> >(ArrayroundC<T, R1 >(a.data()));
  }
  //----
  // fabs
  template<typename T, typename Op1 >
  struct ArrayfabsB{
    const Op1 &op1;
    ArrayfabsB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::fabs ; return fabs(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayfabsC{
    const Op1 op1;
    ArrayfabsC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::fabs ; return fabs(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayfabsC<T, R1> >
  fabs (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayfabsC<T, R1> >(ArrayfabsC<T, R1 >(a.data()));
  }
  //----
  // abs
  template<typename T, typename Op1 >
  struct ArrayabsB{
    const Op1 &op1;
    ArrayabsB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::abs ; return abs(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArrayabsC{
    const Op1 op1;
    ArrayabsC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::abs ; return abs(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArrayabsC<T, R1> >
  abs (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArrayabsC<T, R1> >(ArrayabsC<T, R1 >(a.data()));
  }

  //---------------------------------------------------------------------------
  // sqrt operator
  // 
  template<typename T, typename Op1 >
  struct ArraysqrtB{
    const Op1 &op1;
    ArraysqrtB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::sqrt ; return sqrt(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraysqrtC{
    const Op1 op1;
    ArraysqrtC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::sqrt ; return sqrt(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraysqrtC<T, R1> >
  sqrt (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraysqrtC<T, R1> >(ArraysqrtC<T, R1 >(a.data()));
  }
  //---------------------------------------------------------------------------
  // cbrt operator
  // 
  template<typename T, typename Op1 >
  struct ArraycbrtB{
    const Op1 &op1;
    ArraycbrtB(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::cbrt ; return cbrt(op1[i]) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraycbrtC{
    const Op1 op1;
    ArraycbrtC(const Op1& a): op1(a){}
    T operator[](const size_t i) const
    { using std::cbrt ; return cbrt(op1[i]) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraycbrtC<T, R1> >
  cbrt (const ArrayOperator<T, R1>& a){
    return ArrayOperator<T, ArraycbrtC<T, R1> >(ArraycbrtC<T, R1 >(a.data()));
  }

  //---------------------------------------------------------------------------
  // pow operator with integer exponent
  // 
  template<typename T, typename Op1 >
  struct ArraypowiB{
    const Op1 &op1;
    const int e ;
    ArraypowiB(const Op1& a,const int &b): op1(a),e(b) {}
    T operator[](const size_t i) const
    { using std::pow ; return pow(op1[i],e) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraypowiC{
    const Op1 op1;
    const int e ;
    ArraypowiC(const Op1& a, const int &b): op1(a),e(b){}
    T operator[](const size_t i) const
    { using std::pow ; return pow(op1[i],e) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraypowiC<T, R1> >
  pow (const ArrayOperator<T, R1>& a,const int &e){
    return ArrayOperator<T, ArraypowiC<T, R1> >(ArraypowiC<T, R1 >(a.data(),e));
  }
  //---------------------------------------------------------------------------
  // pow operator with double exponent
  // 
  template<typename T, typename Op1 >
  struct ArraypowdB{
    const Op1 &op1;
    const double e ;
    ArraypowdB(const Op1& a,const double &b): op1(a),e(b) {}
    T operator[](const size_t i) const
    { using std::pow ; return pow(op1[i],e) ;  }
  } ;
  template<typename T, typename Op1 >
  struct ArraypowdC{
    const Op1 op1;
    const double e ;
    ArraypowdC(const Op1& a, const double &b): op1(a),e(b){}
    T operator[](const size_t i) const
    { using std::pow ; return pow(op1[i],e) ;  }
  } ;

  template<typename T, typename R1>
  inline ArrayOperator<T, ArraypowdC<T, R1> >
  pow (const ArrayOperator<T, R1>& a,const double &e){
    return ArrayOperator<T, ArraypowdC<T, R1> >(ArraypowdC<T, R1 >(a.data(),e));
  }

  //---------------------------------------------------------------------------
  // General Pow Operator
  //
  // Case when both left and right are concrete (not anonymous)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayPow{
    const Op1& op1;
    const Op2& op2;
    ArrayPow(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const { using std::pow ; return pow(op1[i] , op2[i]) ;  }
  } ;
  // Case where left argument of pow is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayPowL{
    const Op1 op1;
    const Op2& op2;
    ArrayPowL(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const
    { using std::pow ; return pow(op1[i], op2[i]) ;  }
  } ;

  // case where right argument of pow is a temporary (thus we need to copy)
  template<typename T, typename Op1 , typename Op2>
  struct ArrayPowR{
    const Op1 &op1;
    const Op2 op2;
    ArrayPowR(const Op1& a, const Op2& b): op1(a), op2(b){}
    T operator[](const size_t i) const
    { using std::pow ; return pow(op1[i],op2[i]) ;  }
  } ;

  template<typename T, typename Op1 , typename Op2>
  struct ArrayPowLR{
    const Op1 op1;
    const Op2 op2;
    ArrayPowLR(const Op1& a, const Op2& b): op1(a), op2(b) {}
    T operator[](const size_t i) const
    { using std::pow ; return pow(op1[i] , op2[i]) ;  }
  } ;
  // General operators for packaged expressions
  // function template for the * operator
  template<typename T, typename R1, typename R2>
  inline ArrayOperator<T, ArrayPowLR<T, R1, R2> >
  pow (const ArrayOperator<T, R1>& a, const ArrayOperator<T, R2>& b){
    return ArrayOperator<T, ArrayPowLR<T, R1, R2> >(ArrayPowLR<T, R1, R2 >(a.data(), b.data()));
  }

  //---------------------Array----------------------//
  template <class T,size_t n> class Array {
    T x[n] ;
  public:
    typedef T * iterator ;
    typedef const T * const_iterator ;
    
    //    Array() {} ;
    //    Array(const Array<T,n> &v)
    //    { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; } 
    //    Array<T,n> &operator=(const Array<T,n> &v)
    //    { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; return *this ; } 

    Array<T,n> &operator=(const T &v) {
      for(size_t i=0;i<n;++i)
	x[i] = v ;
      return *this ;
    }
    // Expression template interface for unary operator
    ArrayOperator<T,ArrayUnaryB<T,Array<T,n> > > operator-() {
      return ArrayOperator<T,ArrayUnaryB<T,Array<T,n> > >(ArrayUnaryB<T,Array<T,n> >(*this)) ;
    }
						       
    // Expression template loop hoisting for assignment operators
    template <typename R> Array<T,n> &operator=(const ArrayOperator<T,R> &v) {
      for(size_t i=0;i<n;++i) {
	x[i] = v[i];
      }
      return *this ;
    }
    template <typename R> Array<T,n> &operator+=(const ArrayOperator<T,R> &v) {
      for(size_t i=0;i<n;++i)
	x[i] += v[i];
      return *this ;
    }
    template <typename R> Array<T,n> &operator-=(const ArrayOperator<T,R> &v) {
      for(size_t i=0;i<n;++i)
	x[i] -= v[i];
      return *this ;
    }
    template <typename R> Array<T,n> &operator*=(const ArrayOperator<T,R> &v) {
      for(size_t i=0;i<n;++i)
	x[i] *= v[i];
      return *this ;
    }
    template <typename R> Array<T,n> &operator/=(const ArrayOperator<T,R> &v) {
      for(size_t i=0;i<n;++i)
	x[i] /= v[i];
      return *this ;
    }

    Array<T,n> &operator +=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] += v.x[i] ; return *this ; }
    Array<T,n> &operator -=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; }
    Array<T,n> &operator *=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; }
    Array<T,n> &operator /=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] /= v.x[i] ; return *this ; }

    T &operator[](size_t indx) { return x[indx]; }
    const T &operator[](size_t indx) const { return x[indx] ; }
    T &operator[](int indx) { return x[indx]; }
    const T &operator[](int indx) const { return x[indx] ; }
    T &operator[](unsigned char indx) { return x[indx]; }
    const T &operator[](unsigned char indx) const { return x[indx] ; }
    T &operator[](unsigned int indx) { return x[indx]; }
    const T &operator[](unsigned int indx) const { return x[indx] ; }

    iterator begin() { return &x[0] ; }
    iterator end() { return begin()+n ; }
    const_iterator begin() const { return &x[0] ; }
    const_iterator end() const { return begin()+n ; }

    size_t size() const  { return n ; }
    
  } ;

  //---------------------------------------------------------------------------
  // Gateway operators for expression template facility
  // Note that we have several operators here to control pass by reference
  // and pass by value to prevent orphaned data from the destruction of
  // intermediate operators
  //
  //---------------------------------------------------------------------------
  // Summation
  // Operator for the sum of two arrays
  template<typename T, size_t n>
  inline ArrayOperator<T, ArrayAdd<T,Array<T,n>,Array<T,n> > >
  operator+ (const Array<T, n>& a, const Array<T, n>& b){
    return ArrayOperator<T, ArrayAdd<T, Array<T,n>, Array<T,n> > >(ArrayAdd<T, Array<T,n>, Array<T,n> >(a, b));
  }
  // Adding an expression to an array
  template<typename T, size_t n, typename R>
  inline ArrayOperator<T,ArrayAddL<T,R,Array<T,n> > >
  operator+(const R &a, const Array<T,n> &b) {
    return ArrayOperator<T,ArrayAddL<T,R,Array<T,n> > >(ArrayAddL<T,R,Array<T,n> >(a.data(),b)) ;
  }
  // Adding an array to an expression
  template<typename T, size_t n,typename R >
  inline ArrayOperator<T,ArrayAddR<T,Array<T,n>,R > >
  operator+(const Array<T,n> &a, const  R &b) {
    return ArrayOperator<T,ArrayAddR<T,Array<T,n>,R > >(ArrayAddR<T,Array<T,n>,R>(a,b.data())) ;
  }
  // Operator for constant + array
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayAddLC<T,Array<T,n> > >
  operator+(const double a, const Array<T,n> &b) {
    return ArrayOperator<T,ArrayAddLC<T,Array<T,n> > >(ArrayAddLC<T,Array<T,n> >(a,b)) ;
  }

  // Operator for array + constant
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayAddRC<T,Array<T,n> > >
  operator+(const Array<T,n> &a, const  double b) {
    return ArrayOperator<T,ArrayAddRC<T,Array<T,n> > >(ArrayAddRC<T,Array<T,n> >(a,b)) ;
  }
  //---------------------------------------------------------------------------
  // Subtraction
  // Operator for the sum of two arrays
  template<typename T, size_t n>
  inline ArrayOperator<T, ArraySub<T,Array<T,n>,Array<T,n> > >
  operator- (const Array<T, n>& a, const Array<T, n>& b){
    return ArrayOperator<T, ArraySub<T, Array<T,n>, Array<T,n> > >(ArraySub<T, Array<T,n>, Array<T,n> >(a, b));
  }
  // Subing an expression to an array
  template<typename T, size_t n, typename R>
  inline ArrayOperator<T,ArraySubL<T,R,Array<T,n> > >
  operator-(const R &a, const Array<T,n> &b) {
    return ArrayOperator<T,ArraySubL<T,R,Array<T,n> > >(ArraySubL<T,R,Array<T,n> >(a.data(),b)) ;
  }
  // Subing an array to an expression
  template<typename T, size_t n,typename R >
  inline ArrayOperator<T,ArraySubR<T,Array<T,n>,R > >
  operator-(const Array<T,n> &a, const  R &b) {
    return ArrayOperator<T,ArraySubR<T,Array<T,n>,R > >(ArraySubR<T,Array<T,n>,R>(a,b.data())) ;
  }
  // Operator for constant - array
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraySubLC<T,Array<T,n> > >
  operator-(const double a, const Array<T,n> &b) {
    return ArrayOperator<T,ArraySubLC<T,Array<T,n> > >(ArraySubLC<T,Array<T,n> >(a,b)) ;
  }

  // Operator for array - constant
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraySubRC<T,Array<T,n> > >
  operator-(const Array<T,n> &a, const  double b) {
    return ArrayOperator<T,ArraySubRC<T,Array<T,n> > >(ArraySubRC<T,Array<T,n> >(a,b)) ;
  }
  //---------------------------------------------------------------------------
  // Multiplication
  // Operator for the sum of two arrays
  template<typename T, size_t n>
  inline ArrayOperator<T, ArrayMul<T,Array<T,n>,Array<T,n> > >
  operator* (const Array<T, n>& a, const Array<T, n>& b){
    return ArrayOperator<T, ArrayMul<T, Array<T,n>, Array<T,n> > >(ArrayMul<T, Array<T,n>, Array<T,n> >(a, b));
  }
  // Muling an expression to an array
  template<typename T, size_t n, typename R>
  inline ArrayOperator<T,ArrayMulL<T,R,Array<T,n> > >
  operator*(const R &a, const Array<T,n> &b) {
    return ArrayOperator<T,ArrayMulL<T,R,Array<T,n> > >(ArrayMulL<T,R,Array<T,n> >(a.data(),b)) ;
  }
  // Muling an array to an expression
  template<typename T, size_t n,typename R >
  inline ArrayOperator<T,ArrayMulR<T,Array<T,n>,R > >
  operator*(const Array<T,n> &a, const  R &b) {
    return ArrayOperator<T,ArrayMulR<T,Array<T,n>,R > >(ArrayMulR<T,Array<T,n>,R>(a,b.data())) ;
  }
  // Operator for constant * array
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayMulLC<T,Array<T,n> > >
  operator*(const double a, const Array<T,n> &b) {
    return ArrayOperator<T,ArrayMulLC<T,Array<T,n> > >(ArrayMulLC<T,Array<T,n> >(a,b)) ;
  }

  // Operator for array * constant
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayMulRC<T,Array<T,n> > >
  operator*(const Array<T,n> &a, const  double b) {
    return ArrayOperator<T,ArrayMulRC<T,Array<T,n> > >(ArrayMulRC<T,Array<T,n> >(a,b)) ;
  }

  //----------------------------------------------------------------------------
  // Division
  // Operator for the sum of two arrays
  template<typename T, size_t n>
  inline ArrayOperator<T, ArrayDiv<T,Array<T,n>,Array<T,n> > >
  operator/ (const Array<T, n>& a, const Array<T, n>& b){
    return ArrayOperator<T, ArrayDiv<T, Array<T,n>, Array<T,n> > >(ArrayDiv<T, Array<T,n>, Array<T,n> >(a, b));
  }
  // Diving an expression to an array
  template<typename T, size_t n, typename R>
  inline ArrayOperator<T,ArrayDivL<T,R,Array<T,n> > >
  operator/(const R &a, const Array<T,n> &b) {
    return ArrayOperator<T,ArrayDivL<T,R,Array<T,n> > >(ArrayDivL<T,R,Array<T,n> >(a.data(),b)) ;
  }
  // Diving an array to an expression
  template<typename T, size_t n,typename R >
  inline ArrayOperator<T,ArrayDivR<T,Array<T,n>,R > >
  operator/(const Array<T,n> &a, const  R &b) {
    return ArrayOperator<T,ArrayDivR<T,Array<T,n>,R > >(ArrayDivR<T,Array<T,n>,R>(a,b.data())) ;
  }
  // Operator for constant / array
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayDivLC<T,Array<T,n> > >
  operator/(const double a, const Array<T,n> &b) {
    return ArrayOperator<T,ArrayDivLC<T,Array<T,n> > >(ArrayDivLC<T,Array<T,n> >(a,b)) ;
  }

  // Operator for array / constant
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayDivRC<T,Array<T,n> > >
  operator/(const Array<T,n> &a, const  double b) {
    return ArrayOperator<T,ArrayDivRC<T,Array<T,n> > >(ArrayDivRC<T,Array<T,n> >(a,b)) ;
  }

  //---------------------------------------------------------------------------
  // trigonometric 
  //----
  // cos
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraycosB<T,Array<T,n> > > cos(const Array<T,n> &a) {
    return ArrayOperator<T,ArraycosB<T,Array<T,n> > >(ArraycosB<T,Array<T,n> >(a)) ;
  }
  //----
  // sin
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraysinB<T,Array<T,n> > > sin(const Array<T,n> &a) {
    return ArrayOperator<T,ArraysinB<T,Array<T,n> > >(ArraysinB<T,Array<T,n> >(a)) ;
  }
  //----
  // tan
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraytanB<T,Array<T,n> > > tan(const Array<T,n> &a) {
    return ArrayOperator<T,ArraytanB<T,Array<T,n> > >(ArraytanB<T,Array<T,n> >(a)) ;
  }
  //----
  // acos
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayacosB<T,Array<T,n> > > acos(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayacosB<T,Array<T,n> > >(ArrayacosB<T,Array<T,n> >(a)) ;
  }
  //----
  // asin
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayasinB<T,Array<T,n> > > asin(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayasinB<T,Array<T,n> > >(ArrayasinB<T,Array<T,n> >(a)) ;
  }
  //----
  // atan
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayatanB<T,Array<T,n> > > atan(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayatanB<T,Array<T,n> > >(ArrayatanB<T,Array<T,n> >(a)) ;
  }
  //---------------------------------------------------------------------------
  // hyperbolic 
  //----
  // cosh
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraycoshB<T,Array<T,n> > > cosh(const Array<T,n> &a) {
    return ArrayOperator<T,ArraycoshB<T,Array<T,n> > >(ArraycoshB<T,Array<T,n> >(a)) ;
  }
  //----
  // sinh
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraysinhB<T,Array<T,n> > > sinh(const Array<T,n> &a) {
    return ArrayOperator<T,ArraysinhB<T,Array<T,n> > >(ArraysinhB<T,Array<T,n> >(a)) ;
  }
  //----
  // tanh
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraytanhB<T,Array<T,n> > > tanh(const Array<T,n> &a) {
    return ArrayOperator<T,ArraytanhB<T,Array<T,n> > >(ArraytanhB<T,Array<T,n> >(a)) ;
  }
  //----
  // acosh
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayacoshB<T,Array<T,n> > > acosh(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayacoshB<T,Array<T,n> > >(ArrayacoshB<T,Array<T,n> >(a)) ;
  }
  //----
  // asinh
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayasinhB<T,Array<T,n> > > asinh(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayasinhB<T,Array<T,n> > >(ArrayasinhB<T,Array<T,n> >(a)) ;
  }
  //----
  // atanh
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayatanhB<T,Array<T,n> > > atanh(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayatanhB<T,Array<T,n> > >(ArrayatanhB<T,Array<T,n> >(a)) ;
  }

  //---------------------------------------------------------------------------
  // exponential operators
  //
  //----
  // exp
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayexpB<T,Array<T,n> > > exp(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayexpB<T,Array<T,n> > >(ArrayexpB<T,Array<T,n> >(a)) ;
  }
  //----
  // log
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraylogB<T,Array<T,n> > > log(const Array<T,n> &a) {
    return ArrayOperator<T,ArraylogB<T,Array<T,n> > >(ArraylogB<T,Array<T,n> >(a)) ;
  }
  //----
  // log10
  template<typename T, size_t n >
  inline ArrayOperator<T,Arraylog10B<T,Array<T,n> > > log10(const Array<T,n> &a) {
    return ArrayOperator<T,Arraylog10B<T,Array<T,n> > >(Arraylog10B<T,Array<T,n> >(a)) ;
  }
  
  //---------------------------------------------------------------------------
  // Error and gamma functions
  //----
  // erf
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayerfB<T,Array<T,n> > > erf(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayerfB<T,Array<T,n> > >(ArrayerfB<T,Array<T,n> >(a)) ;
  }
  //----
  // erfc
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayerfcB<T,Array<T,n> > > erfc(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayerfcB<T,Array<T,n> > >(ArrayerfcB<T,Array<T,n> >(a)) ;
  }
  //----
  // tgamma
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraytgammaB<T,Array<T,n> > > tgamma(const Array<T,n> &a) {
    return ArrayOperator<T,ArraytgammaB<T,Array<T,n> > >(ArraytgammaB<T,Array<T,n> >(a)) ;
  }
  //----
  // lgamma
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraylgammaB<T,Array<T,n> > > lgamma(const Array<T,n> &a) {
    return ArrayOperator<T,ArraylgammaB<T,Array<T,n> > >(ArraylgammaB<T,Array<T,n> >(a)) ;
  }

  //---------------------------------------------------------------------------
  // misc functions ceil, floor,trunc,round,fabs,abs
  //
  //----
  // ceil
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayceilB<T,Array<T,n> > > ceil(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayceilB<T,Array<T,n> > >(ArrayceilB<T,Array<T,n> >(a)) ;
  }
  //----
  // floor
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayfloorB<T,Array<T,n> > > floor(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayfloorB<T,Array<T,n> > >(ArrayfloorB<T,Array<T,n> >(a)) ;
  }
  //----
  // trunc
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraytruncB<T,Array<T,n> > > trunc(const Array<T,n> &a) {
    return ArrayOperator<T,ArraytruncB<T,Array<T,n> > >(ArraytruncB<T,Array<T,n> >(a)) ;
  }
  //----
  // round
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayroundB<T,Array<T,n> > > round(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayroundB<T,Array<T,n> > >(ArrayroundB<T,Array<T,n> >(a)) ;
  }
  //----
  // fabs
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayfabsB<T,Array<T,n> > > fabs(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayfabsB<T,Array<T,n> > >(ArrayfabsB<T,Array<T,n> >(a)) ;
  }
  //----
  // abs
  template<typename T, size_t n >
  inline ArrayOperator<T,ArrayabsB<T,Array<T,n> > > abs(const Array<T,n> &a) {
    return ArrayOperator<T,ArrayabsB<T,Array<T,n> > >(ArrayabsB<T,Array<T,n> >(a)) ;
  }

  //---------------------------------------------------------------------------
  // sqrt
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraysqrtB<T,Array<T,n> > > sqrt(const Array<T,n> &a) {
    return ArrayOperator<T,ArraysqrtB<T,Array<T,n> > >(ArraysqrtB<T,Array<T,n> >(a)) ;
  }
  //---------------------------------------------------------------------------
  // cbrt
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraycbrtB<T,Array<T,n> > > cbrt(const Array<T,n> &a) {
    return ArrayOperator<T,ArraycbrtB<T,Array<T,n> > >(ArraycbrtB<T,Array<T,n> >(a)) ;
  }
  //---------------------------------------------------------------------------
  // pow with integer exponent
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraypowiB<T,Array<T,n> > > pow(const Array<T,n> &a, const int b) {
    return ArrayOperator<T,ArraypowiB<T,Array<T,n> > >(ArraypowiB<T,Array<T,n> >(a,b)) ;
  }
  //---------------------------------------------------------------------------
  // pow with double exponent
  template<typename T, size_t n >
  inline ArrayOperator<T,ArraypowdB<T,Array<T,n> > > pow(const Array<T,n> &a, const double b) {
    return ArrayOperator<T,ArraypowdB<T,Array<T,n> > >(ArraypowdB<T,Array<T,n> >(a,b)) ;
  }
  //---------------------------------------------------------------------------
  // general pow (arrays for both arguments)
  //---------------------------------------------------------------------------
  template<typename T, size_t n>
  inline ArrayOperator<T, ArrayPow<T,Array<T,n>,Array<T,n> > >
  pow(const Array<T, n>& a, const Array<T, n>& b){
    return ArrayOperator<T, ArrayPow<T, Array<T,n>, Array<T,n> > >(ArrayPow<T, Array<T,n>, Array<T,n> >(a, b));
  }
  template<typename T, size_t n, typename R>
  inline ArrayOperator<T,ArrayPowL<T,R,Array<T,n> > >
  pow(const R &a, const Array<T,n> &b) {
    return ArrayOperator<T,ArrayPowL<T,R,Array<T,n> > >(ArrayPowL<T,R,Array<T,n> >(a.data(),b)) ;
  }
  template<typename T, size_t n,typename R >
  inline ArrayOperator<T,ArrayPowR<T,Array<T,n>,R > >
  pow(const Array<T,n> &a, const  R &b) {
    return ArrayOperator<T,ArrayPowR<T,Array<T,n>,R > >(ArrayPowR<T,Array<T,n>,R>(a,b.data())) ;
  }

  template<typename T, typename F, typename Op>
  struct ArrayFunctorB {
    const F &f ;
    const Op &op ;
    ArrayFunctorB(const F &fin,const Op o):f(fin),op(o) {}
    T operator[](const size_t i) const
    { return f(op[i]) ; }
  } ;
  template<typename T, typename F, typename Op>
  struct ArrayFunctorC {
    const F &f ;
    const Op op ;
    ArrayFunctorC(const F &fin,const Op o):f(fin),op(o) {}
    T operator[](const size_t i) const
    { return f(op[i]) ; }
  } ;
  
  template<typename T, typename F,typename R>
  inline  ArrayOperator<T,ArrayFunctorC<T,F,R> > mapArray(const F &f, const ArrayOperator<T,R> &a) {
    return ArrayOperator<T,ArrayFunctorC<T,F,R> >(ArrayFunctorC<T,F,R>(f,a.data())) ;
  }
  template<typename T, typename F, size_t n>
  inline ArrayOperator<T,ArrayFunctorB<T,F,Array<T,n> > > mapArray(const F &f, const Array<T,n> &a) {
    return ArrayOperator<T,ArrayFunctorB<T,F,Array<T,n> > >(ArrayFunctorB<T,F,Array<T,n> >(f,a)) ;
  }

  //---------------------vector3d------------------//
  template <class T> 
  struct vector3d {
    T x,y,z ;
    vector3d() {} 
    vector3d(const vector3d &v): x(v.x),y(v.y),z(v.z) {}
    vector3d(T xx,T yy, T zz) : x(xx),y(yy),z(zz) {}
    template <class S> vector3d(const vector3d<S> &v): x(v.x),y(v.y),z(v.z) {}
    //    template <class S> vector3d(const Array<S,3> &a) {x=a[0];y=a[1];z=a[2];}
    // template <class S> operator Array<S,3>() {
    //   Array<S,3> a ;
    //   a[0] = x ;
    //   a[1] = y ;
    //   a[2] = z ;
    //   return a;
    // }
    T &operator[](size_t i) { FATAL(i>2) ; return (&x)[i] ; }
    const T &operator[](size_t i) const { FATAL(i>2) ; return (&x)[i] ; }
  } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector3d<T> &v)
  {
    s << v.x << ' ' << v.y << ' ' << v.z << ' ' ;
    return s ;
  }

  template <class T> inline std::istream &operator>>(std::istream &s, vector3d<T> &v)
  {
    s >> v.x >> v.y >> v.z ;
    return s ;
  }

  template <class T> inline T dot(const vector3d<T> &v1, const vector3d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z ;
  }

  template <class T> inline T norm(const vector3d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.y*v2.z-v1.z*v2.y,
                       v1.z*v2.x-v1.x*v2.z,
                       v1.x*v2.y-v1.y*v2.x) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const T ra2[]) {
    return vector3d<T>(v1.y*ra2[2]-v1.z*ra2[1],
                       v1.z*ra2[0]-v1.x*ra2[2],
                       v1.x*ra2[1]-v1.y*ra2[0]) ;
  }
  template<class T> inline vector3d<T> cross(const T ra1[], const vector3d<T> &v2) {
    return vector3d<T>(ra1[1]*v2.z-ra1[2]*v2.y,
                       ra1[2]*v2.x-ra1[0]*v2.z,
                       ra1[0]*v2.y-ra1[1]*v2.x) ;
  }

  template<class T,class S> inline vector3d<T> &operator*=(vector3d<T> &target, S val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T,class S> inline vector3d<T> &operator/=(vector3d<T> &target, S val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }
  
  template<class T> inline vector3d<T> operator+=(vector3d<T> &target, const vector3d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    target.z += val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator-=(vector3d<T> &target, const vector3d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    target.z -= val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z) ;
  }

  template<class T> inline vector3d<T> operator-(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z) ;
  }

  template<class T,class S> inline vector3d<T> operator*(const vector3d<T> &v1, S r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }
  template<class T,class S> inline vector3d<T> operator*(S r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T,class S> inline vector3d<T> operator/(const vector3d<T> &v1, S r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T>  struct tensor3d : public vector3d<vector3d< T > > {
    tensor3d() {}
    tensor3d(vector3d<T> xx,vector3d<T> yy, vector3d<T> zz)
      : vector3d<vector3d< T> > (xx,yy,zz) {}
    tensor3d(const tensor3d &v) : vector3d<vector3d< T> >(v) {}
  } ;

  template<class T> inline vector3d<T> dot(const tensor3d<T> &t,
                                           const vector3d<T> &v) {
    return vector3d<T>(dot(t.x,v),dot(t.y,v),dot(t.z,v)) ;
  }

  template<class T> inline tensor3d<T> product(const tensor3d<T> &t1,
                                               const tensor3d<T> &t2) {
    tensor3d<T> temp ;
    temp.x.x = t1.x.x*t2.x.x+t1.x.y*t2.y.x+t1.x.z*t2.z.x ;
    temp.y.x = t1.y.x*t2.x.x+t1.y.y*t2.y.x+t1.y.z*t2.z.x ;
    temp.z.x = t1.z.x*t2.x.x+t1.z.y*t2.y.x+t1.z.z*t2.z.x ;

    temp.x.y = t1.x.x*t2.x.y+t1.x.y*t2.y.y+t1.x.z*t2.z.y ;
    temp.y.y = t1.y.x*t2.x.y+t1.y.y*t2.y.y+t1.y.z*t2.z.y ;
    temp.z.y = t1.z.x*t2.x.y+t1.z.y*t2.y.y+t1.z.z*t2.z.y ;

    temp.x.z = t1.x.x*t2.x.z+t1.x.y*t2.y.z+t1.x.z*t2.z.z ;
    temp.y.z = t1.y.x*t2.x.z+t1.y.y*t2.y.z+t1.y.z*t2.z.z ;
    temp.z.z = t1.z.x*t2.x.z+t1.z.y*t2.y.z+t1.z.z*t2.z.z ;

    return temp ;
  }

  inline vector3d<float> realToFloat(vector3d<double> v) { return vector3d<float>(float(v.x),float(v.y),float(v.z)); }
  inline vector3d<double> realToDouble(vector3d<double> v) { return v ; }


#ifndef NO_AUTODIFF
  inline vector3d<float> realToFloat(vector3d<MFADd> v) { return vector3d<float>(realToFloat(v.x),realToFloat(v.y),realToFloat(v.z)); }
  inline vector3d<double> realToDouble(vector3d<MFADd> v) { return vector3d<double>(realToDouble(v.x),realToDouble(v.y),realToDouble(v.z)); }

  inline vector3d<float> realToFloat(vector3d<FADd> v) { return vector3d<float>(realToFloat(v.x),realToFloat(v.y),realToFloat(v.z)); }
  inline vector3d<double> realToDouble(vector3d<FADd> v) { return vector3d<double>(realToDouble(v.x),realToDouble(v.y),realToDouble(v.z)); }

  inline vector3d<float> realToFloat(vector3d<FAD2d> v) { return vector3d<float>(realToFloat(v.x),realToFloat(v.y),realToFloat(v.z)); }
  inline vector3d<double> realToDouble(vector3d<FAD2d> v) { return vector3d<double>(realToDouble(v.x),realToDouble(v.y),realToDouble(v.z)); }
#endif
  
  //---------------------vector2d------------------//
  template <class T> 
  struct vector2d {
    T x,y ;
    vector2d() {} 
    vector2d(T xx,T yy) : x(xx),y(yy) {}
    vector2d(const vector2d &v) {x=v.x;y=v.y;}
    T &operator[](int i) {
      FATAL(i>1) ;
      switch(i) {
      case 0:
        return x ;
      case 1:
        return y ;
      default:
        return y ;
      }
    }
    const T &operator[](int i) const {
      FATAL(i>1) ;
      switch(i) {
      case 0:
        return x ;
      case 1:
        return y ;
      default:
        return y ;
      }
    }
  } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector2d<T> &v)
  {
    s << v.x << ' ' << v.y << ' ' ;
    return s ;
  }

  template <class T> inline std::istream &operator>>(std::istream &s, vector2d<T> &v)
  {
    s >> v.x >> v.y ;
    return s ;
  }

  template <class T> inline T dot(const vector2d<T> &v1, const vector2d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y ;
  }

  template <class T> inline T dot(const vector2d<T> &v1, const T ra2[]) {
    return v1.x*ra2[0] + v1.y*ra2[1] ;
  }

  template <class T> inline T norm(const vector2d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y) ;
  }

  template<class T> inline T cross(const vector2d<T> &v1, const vector2d<T> &v2) {
    return v1.x*v2.y-v1.y*v2.x ;
  }

  template<class T> inline T cross(const vector2d<T> &v1, const T ra2[]) {
    return v1.x*ra2[1]-v1.y*ra2[0] ;
  }


  template<class T,class S> inline vector2d<T> &operator*=(vector2d<T> &target, S val) {
    target.x *= val ;
    target.y *= val ;
    return target ;
  }

  template<class T,class S> inline vector2d<T> &operator/=(vector2d<T> &target, S val) {
    target.x /= val ;
    target.y /= val ;
    return target ;
  }

  template<class T> inline vector2d<T> operator+=(vector2d<T> &target, const vector2d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    return target ;
  }

  template<class T> inline vector2d<T> operator-=(vector2d<T> &target, const vector2d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    return target ;
  }

  template<class T> inline vector2d<T> operator+(const vector2d<T> &v1, const vector2d<T> &v2) {
    return vector2d<T>(v1.x+v2.x,v1.y+v2.y) ;
  }

  template<class T> inline vector2d<T> operator-(const vector2d<T> &v1, const vector2d<T> &v2) {
    return vector2d<T>(v1.x-v2.x,v1.y-v2.y) ;
  }

  template<class T,class S> inline vector2d<T> operator*(const vector2d<T> &v1, S r2) {
    return vector2d<T>(v1.x*r2,v1.y*r2) ;
  }

  template<class T,class S> inline vector2d<T> operator*(S r1, const vector2d<T> &v2) {
    return vector2d<T>(v2.x*r1,v2.y*r1) ;
  }

  template<class T,class S> inline vector2d<T> operator/(const vector2d<T> &v1, S r2) {
    return vector2d<T>(v1.x/r2,v1.y/r2) ;
  }

  // For allocating temporary arrays of small size
  const int tmp_array_internal_SIZE=25 ;
  template <class T> class tmp_array {
    int sz ;
    T data[tmp_array_internal_SIZE] ;
    T * p ;
    void alloc(int size) {
      sz = size ;
      p = data ;
      if(sz > tmp_array_internal_SIZE)
        p = new T[sz] ;
    }
    void free() {
      if(sz > tmp_array_internal_SIZE)
        delete[] p ;
    }
    tmp_array() { alloc(0) ; }
  public:
    tmp_array(int size) {
      alloc(size) ;
    }
    tmp_array(const tmp_array &ta) {
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
    }
    tmp_array &operator=(const tmp_array &ta) {
      free() ;
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
      return *this ;
    }
    ~tmp_array() { free(); }
    T & operator[](int i) { FATAL(i<0 || i>=sz) ; return p[i] ; }
    T & operator[](int i) const { FATAL(i<0 || i>=sz) ; return p[i] ; }
    T & operator[](size_t i) { FATAL(i>=sz) ; return p[i] ; }
    T & operator[](size_t i) const { FATAL(i>=sz) ; return p[i] ; }
    T & operator[](unsigned int i) { FATAL(i>=sz) ; return p[i] ; }
    T & operator[](unsigned int i) const { FATAL(i>=sz) ; return p[i] ; }
    T & operator[](unsigned char i) { FATAL(i>=sz) ; return p[i] ; }
    T & operator[](unsigned char i) const { FATAL(i>=sz) ; return p[i] ; }
    T & operator[](unsigned short i) { FATAL(i>=sz) ; return p[i] ; }
    T & operator[](unsigned short i) const { FATAL(i>=sz) ; return p[i] ; }
    operator T *() { return p ; }
    operator const T *() const { return p ; }
  } ;
}

#endif
