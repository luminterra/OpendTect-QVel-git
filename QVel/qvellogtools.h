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
#ifndef qvellogtools_h
#define qvellogtools_h
/*+
 * (C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
 * AUTHOR   : R.K. Singh
 * DATE     : June 2007
 * ID       : $Id: qvellogtools.h,v 1.3 2009-07-22 16:01:27 cvsbert Exp $
-*/

#include "commondefs.h"

namespace Well { class Log; }

namespace QVel
{

mClass LogTools
{
public:

    			LogTools(const Well::Log& input,Well::Log& output);
			   		    

    bool		runSmooth(int gate);

protected:

    const Well::Log&	inplog_;
    Well::Log&		outplog_;

};

} // namespace

#endif
