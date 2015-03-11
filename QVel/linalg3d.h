#ifndef __GEOM_H
#define __GEOM_H

/*
 *  One file long C++ library of linear algebra primitives for
 *  simple 3D programs
 *
 *  Copyright (C) 2001-2003 by Jarno Elonen
 *
 *  Permission to use, copy, modify, distribute and sell this software
 *  and its documentation for any purpose is hereby granted without fee,
 *  provided that the above copyright notice appear in all copies and
 *  that both that copyright notice and this permission notice appear
 *  in supporting documentation.  The authors make no representations
 *  about the suitability of this software for any purpose.
 *  It is provided "as is" without express or implied warranty.
 */

#include <cmath>

#define EPSILON 0.00001f
#define PI 3.1415926
#define Deg2Rad(Ang) ((float)( Ang * PI / 180.0 ))
#define Rad2Deg(Ang) ((float)( Ang * 180.0 / PI ))

// =========================================
// 3-vector
// =========================================
class Vec
{
public:

  // Position
  float x, y, z;

  // Default constructor
  Vec()
  : x( 0 ), y( 0 ), z( 0 ) {}

  // Element constructor
  Vec( float x1, float y1, float z1 )
  : x( x1 ), y( y1 ), z( z1 ) {}

  // Copy constructor
  Vec( const Vec& a )
  : x( a.x ), y( a.y ), z( a.z ) {}

  // Norm (len^2)
  inline float norm() const { return x*x + y*y + z*z; }

  // Length of the vector
  inline float len() const { return (float)sqrt(norm()); }

  Vec &operator += ( const Vec &src ) { x += src.x; y += src.y; z += src.z; return *this; }
  Vec operator + ( const Vec &src ) const { Vec tmp( *this ); return ( tmp += src ); }
  Vec &operator -= ( const Vec &src ) { x -= src.x; y -= src.y; z -= src.z; return *this; }
  Vec operator - ( const Vec &src ) const { Vec tmp( *this ); return ( tmp -= src ); }

  Vec operator - () const { return Vec(-x,-y,-z); }

  Vec &operator *= ( const float src ) { x *= src; y *= src; z *= src;  return *this; }
  Vec operator * ( const float src ) const { Vec tmp( *this ); return ( tmp *= src ); }
  Vec &operator /= ( const float src ) { x /= src; y /= src; z /= src; return *this; }
  Vec operator / ( const float src ) const { Vec tmp( *this ); return ( tmp /= src ); }

  bool operator == ( const Vec& b) const { return ((*this)-b).norm() < EPSILON; }
  //bool operator == ( const Vec& b) const { return x==b.x && y==b.y && z==b.z; }
};

  // Left hand float multplication
  inline Vec operator * ( const float src, const Vec& v ) { Vec tmp( v ); return ( tmp *= src ); }

  // Dot product
  inline float dot( const Vec& a, const Vec& b )
  { return a.x*b.x + a.y*b.y + a.z*b.z; }

  // Cross product
  inline Vec cross( const Vec &a, const Vec &b )
  { return Vec( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x ); }

#endif
