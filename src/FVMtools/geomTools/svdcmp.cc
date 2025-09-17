//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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
/*******************************************************************************
Singular value decomposition program, svdcmp, from "Numerical Recipes in C"
(Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
and B.P. Flannery
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include "dmatrix.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
inline double DMAX(double a, double b) { return a>b?a:b; }
inline int IMIN(int a,int b) { return a<b?a:b ; }

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

//c version
void reorder(dmatrix &u, int m , int n, dvector &w, dmatrix &v) { 
  /*Given the output of decompose, 
    this routine sorts the singular values, 
    and corresponding columns of u and v, by decreasing magnitude.
    Also, signs of corresponding columns are flipped 
    so as to maximize the number of positive elements. */

  int i,j,k,s,inc=1; 
  double sw; 
  dvector su = dvector(1,m); 
  dvector sv = dvector(1,n); 
  do { inc *= 3; inc++; } while (inc <= n); 
  /*Sort. The method is Shell's sort.*/ 
  do { 
    /*(The work is negligible as compared to that already done in decompose)*/
    inc /= 3; 
    for (i=inc;i<n;i++){ 
      sw = w[i+1]; 
      for (k=0;k<m;k++) su[k+1] = u[k+1][i+1]; 
      for (k=0;k<n;k++) sv[k+1] = v[k+1][i+1]; 
      j = i; 
      while (w[j-inc+1] < sw){ 
        w[j+1] = w[j-inc+1]; 
        for (k=0;k<m;k++) u[k+1][j+1] = u[k+1][j-inc+1]; 
        for (k=0;k<n;k++) v[k+1][j+1] = v[k+1][j-inc+1]; 
        j -= inc; 
        if (j < inc) break; 
      } 
      w[j+1] = sw; 
      for (k=0;k<m;k++) u[k+1][j+1] = su[k+1]; 
      for (k=0;k<n;k++) v[k+1][j+1] = sv[k+1];      
    } 
  } while (inc > 1); 
  for (k=0;k<n;k++){ 
    /*Flip signs.*/ 
    s=0; 
    for (i=0;i<m;i++) if (u[i+1][k+1] < 0.) s++; 
    for (j=0;j<n;j++) if (v[j+1][k+1] < 0.) s++; 
    if (s > (m+n)/2) { 
      for (i=0;i<m;i++) u[i+1][k+1] = -u[i+1][k+1]; 
      for (j=0;j<n;j++) v[j+1][k+1] = -v[j+1][k+1]; 
    } 
  }
} 
void print(dmatrix &a, int m, int n){
  printf("\n");
  for(int i = 1; i <= m ; i++){
    for(int j = 1; j <=n; j++){
      printf(" %f; ", a[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

/******************************************************************************/
void svdcmp(dmatrix &a, int m, int n, dvector &w, dmatrix &v)
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
{
  // printf("start svdcmp \n");
  //print(a,m,n);
  dmatrix a_copy = dmatrix(1,m, 1,n);
  for(int i = 1; i <= m; i++)
    for(int j = 1; j <=n; j++) a_copy[i][j] = a[i][j];
  
  
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
	
  dvector rv1=dvector(1,n);
  g=scale=anorm=0.0; /* Householder reduction to bidiagonal form */
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm = DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) { /* Accumulation of right-hand transformations. */
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++) /* Double division to avoid possible underflow. */
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) { /* Accumulation of left-hand transformations. */
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) { /* Diagonalization of the bidiagonal form. */
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) { /* Test for splitting. */
	nm=l-1; /* Note that rv1[1] is always zero. */
	if ((double)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0; /* Cancellation of rv1[l], if l > 1. */
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) { /* Convergence. */
	if (z < 0.0) { /* Singular value is made nonnegative. */
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30){
	printf("no convergence in 30 svdcmp iterations");
	print(a_copy, m, n);
      }
      x=w[l]; /* Shift from bottom 2-by-2 minor. */
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0; /* Next QR transformation: */
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z; /* Rotation can be arbitrary if z = 0. */
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }

  reorder(a,  m , n, w,v);  
}

  
