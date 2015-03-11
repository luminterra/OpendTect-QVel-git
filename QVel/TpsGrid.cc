/* 
 * Additions Copyright (C) 2012, JB West, LuminTerra, LLC.
 * Derived from code with the following Copyright:
 */
/*
 *  Thin Plate Spline demo/example in C++
 *
 *  - a simple TPS editor, using the Boost uBlas library for large
 *    matrix operations and OpenGL + GLUT for 2D function visualization
 *    (curved plane) and user interface
 *
 *  Copyright (C) 2003,2005 by Jarno Elonen
 *
 *  TPSDemo is Free Software / Open Source with a very permissive
 *  license:
 *
 *  Permission to use, copy, modify, distribute and sell this software
 *  and its documentation for any purpose is hereby granted without fee,
 *  provided that the above copyright notice appear in all copies and
 *  that both that copyright notice and this permission notice appear
 *  in supporting documentation.  The authors make no representations
 *  about the suitability of this software for any purpose.
 *  It is provided "as is" without express or implied warranty.
 *
 *  TODO:
 *    - implement TPS approximation 3 as suggested in paper
 *      Gianluca Donato and Serge Belongie, 2002: "Approximation
 *      Methods for Thin Plate Spline Mappings and Principal Warps"
 */

#include <commondefs.h>
#include <boost/numeric/ublas/matrix.hpp>

#include "linalg3d.h"
#include "ludecomposition.h"

#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace boost::numeric::ublas;

#include "TpsGrid.h"

// ========= BEGIN INTERESTING STUFF  =========

namespace QVel
{

TpsGrid::TpsGrid() :
	 nrows_(0)
	,ncols_(0)
	,grid_(NULL)
	,regularization_(0)
{
  
}

TpsGrid::TpsGrid(int nrow, int ncol)
{
  setSize(nrow, ncol);
}

TpsGrid::~TpsGrid()
{
  for (int i =0; i < nrows_; i++)
    delete [] grid_[i];

  delete[] grid_;
  
}
void TpsGrid::setSize(int nrow, int ncol)
{
  nrows_ = nrow;
  ncols_ = ncol;
  grid_   = new float* [nrows_];
  for (int i =0; i < nrows_; i++)
    grid_[i] = new float[ncols_];
}

void TpsGrid::setControlPoints(std::vector< Vec > cp)
{
   controlPoints_ = cp;
}

double TpsGrid::tps_base_func(double r)
{
  if ( r == 0.0 )
    return 0.0;
  else
    return r*r * log(r);
}

/*
 *  Calculate Thin Plate Spline (TPS) weights from
 *  control points and build a new height grid_ by
 *  interpolating with them.
 */
int TpsGrid::calcTps()
{
  // You We need at least 3 points to define a plane
  if ( controlPoints_.size() < 3 )
    return -1;

  unsigned p = controlPoints_.size();

  // Allocate the matrix and vector
  matrix<double> mtx_l(p+3, p+3);
  matrix<double> mtx_v(p+3, 1);
  matrix<double> mtx_orig_k(p, p);

  // Fill K (p x p, upper left of L) and calculate
  // mean edge length from control points
  //
  // K is symmetrical so we really have to
  // calculate only about half of the coefficients.
  double a = 0.0;
  for ( unsigned i=0; i<p; ++i )
  {
    for ( unsigned j=i+1; j<p; ++j )
    {
      Vec pt_i = controlPoints_[i];
      Vec pt_j = controlPoints_[j];
      pt_i.y = pt_j.y = 0;
      double elen = (pt_i - pt_j).len();
      mtx_l(i,j) = mtx_l(j,i) =
        mtx_orig_k(i,j) = mtx_orig_k(j,i) =
          tps_base_func(elen);
      a += elen * 2; // same for upper & lower tri
    }
  }
  a /= (double)(p*p);

  // Fill the rest of L
  for ( unsigned i=0; i<p; ++i )
  {
    // diagonal: reqularization parameters (lambda * a^2)
    mtx_l(i,i) = mtx_orig_k(i,i) =
      regularization_ * (a*a);

    // P (p x 3, upper right)
    mtx_l(i, p+0) = 1.0;
    mtx_l(i, p+1) = controlPoints_[i].x;
    mtx_l(i, p+2) = controlPoints_[i].z;

    // P transposed (3 x p, bottom left)
    mtx_l(p+0, i) = 1.0;
    mtx_l(p+1, i) = controlPoints_[i].x;
    mtx_l(p+2, i) = controlPoints_[i].z;
  }
  // O (3 x 3, lower right)
  for ( unsigned i=p; i<p+3; ++i )
    for ( unsigned j=p; j<p+3; ++j )
      mtx_l(i,j) = 0.0;


  // Fill the right hand vector V
  for ( unsigned i=0; i<p; ++i )
    mtx_v(i,0) = controlPoints_[i].y;
  mtx_v(p+0, 0) = mtx_v(p+1, 0) = mtx_v(p+2, 0) = 0.0;

  // Solve the linear system "inplace"
  if (0 != LU_Solve(mtx_l, mtx_v))
  {
    puts( "Singular matrix! Aborting." );
    return -1;
  }

  // Interpolate grid_ heights
  for ( int x=-nrows_/2; x<nrows_/2; ++x )
  {
    for ( int z=-ncols_/2; z<ncols_/2; ++z )
    {
      double h = mtx_v(p+0, 0) + mtx_v(p+1, 0)*x + mtx_v(p+2, 0)*z;
      Vec pt_i, pt_cur(x,0,z);
      for ( unsigned i=0; i<p; ++i )
      {
        pt_i = controlPoints_[i];
        pt_i.y = 0;
        h += mtx_v(i,0) * tps_base_func( ( pt_i - pt_cur ).len());
      }
      grid_[x+nrows_/2][z+ncols_/2] = h;
    }
  }

  // Calc bending energy
  matrix<double> w( p, 1 );
  for ( int i=0; i<p; ++i )
    w(i,0) = mtx_v(i,0);
  matrix<double> be = prod( prod<matrix<double> >( trans(w), mtx_orig_k ), w );
  bendingenergy_ = be(0,0);
  return 0;
}

} // namespace

