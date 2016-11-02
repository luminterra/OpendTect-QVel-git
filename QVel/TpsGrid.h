/*
 Copyright (C) 2012 LuminTerra, LLC All rights reserved.

 This file is part of OpendTect and may be used under the terms of:

 The GNU General Public License version 3 or higher, as published by
 the Free Software Foundation.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

 Ver 1.1 JB West 4/2012

*/

#ifndef tpsgrid_h
#define tpsgrid_h

#include "commondefs.h"
#include "linalg3d.h"


namespace QVel
{
class TpsGrid
{
public:
  TpsGrid();
  ~TpsGrid();
  TpsGrid(int nrow, int ncol);
  void setSize(int nrows, int ncols);
  void setControlPoints(std::vector< Vec > cp);
  int  calcTps();
  float            **grid_;
  void setRegularization( float r ) { regularization_ = r; }

private:
  double tps_base_func(double r);
  std::vector< Vec > controlPoints_;
  int                nrows_;
  int                ncols_;
  double             regularization_;
  double			 bendingenergy_;
};

}
#endif // tpsgrid_h
