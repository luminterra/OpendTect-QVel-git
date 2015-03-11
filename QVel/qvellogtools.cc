/*
Copyright (C) 2002-2010 dGB Beheer B.V. All rights reserved.

This file is part of OpendTect and may be used either under the terms of:

1. The GNU General Public License version 3 or higher, as published by
the Free Software Foundation, or
2. The OpendTect Commercial License version 1 or higher, as published by
dGB Beheer B.V., or
3. The OpendTect Academic License version 1 or higher, as published by
dGB Beheer B.V.

For help to determine which license to use, please visit
http://opendtect.org/index.php/licensing.html

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*/
/*+
 * (C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
 * AUTHOR   : R.K. Singh
 * DATE     : June 2007
-*/

static const char* rcsID = "$Id: qvellogtools.cc,v 1.3 2009-07-22 16:01:27 cvsbert Exp $";

#include "qvellogtools.h"
#include "welllog.h"
#include "statruncalc.h"

QVel::LogTools::LogTools( const Well::Log& inp, Well::Log& outp )
	: inplog_(inp)
	, outplog_(outp)
{
}		

bool QVel::LogTools::runSmooth( const int inpgate )
{
    outplog_.setUnitMeasLabel( inplog_.unitMeasLabel() );

    const int gate = inpgate % 2 ? inpgate : inpgate + 1;  
    const int rad = gate / 2;
    Stats::WindowedCalc<float> wcalc(
	    		Stats::RunCalcSetup().require(Stats::Median), gate );
    const int sz = inplog_.size();
    for ( int idx=0; idx<sz+rad; idx++ )
    {
	const int cpos = idx - rad;
	if ( idx < sz )
	{
	    const float inval = inplog_.value(idx);
	    if (!mIsUdf(inval) )
		wcalc += inval;
	    if ( cpos >= rad )
		outplog_.addValue( inplog_.dah(cpos), wcalc.median() );
	}
	else
	    outplog_.addValue( inplog_.dah(cpos), inplog_.value(cpos) );

	if ( cpos<rad && cpos>=0 )
	    outplog_.addValue( inplog_.dah(cpos), inplog_.value(cpos) );
    }

    return true;
}
