#ifndef qveldeps_h
#define qveldeps_h

/*+
________________________________________________________________________

 (C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
 Author:	K. Tingdahl
 Date:		Oct 2013
 RCS:		$Id: moddeps.h.in 32104 2013-10-23 20:11:53Z kristofer.tingdahl@dgbes.com $
________________________________________________________________________


-*/

//
// This file is automatically generated from CMAKE. It contains includes
// to the export declaration of this module. It should normally not be included
// by the current project. Instead, the instantiate<module>export.cc is created
// for that purpose.
//

#include "qvelmod.h"

#if defined ( do_import_export )
# ifndef do_export_QVel
//Temporary allow extern declaration of extern template instantiation
#  pragma warning( push )
#  pragma warning( disable : 4231 )

#  pragma warning( pop )
# else
#  pragma message ( "Should never come here." )
// Don't include this file from it's own module.
# endif

#endif

#endif