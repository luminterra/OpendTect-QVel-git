#ifndef uiqvelmod_h
#define uiqvelmod_h

/*+
________________________________________________________________________

 (C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
 Author:	K. Tingdahl
 Date:		Mar 2006
 RCS:		$Id$
________________________________________________________________________


-*/

//
//This file is generated automatically from CMAKE. It contains export/import
//Declarations that are used on windows. It also includes deps of all modules
//that this module is dependent on.
//

#if defined( __win64__ ) || defined ( __win32__ )
# define do_import_export
#else
# ifdef do_import_export
#  undef do_import_export
# endif
#endif

#ifndef dll_export
# if defined( do_import_export )
#  define dll_export	__declspec( dllexport )
#  define dll_import	__declspec( dllimport )
#  define dll_extern	extern
# else
#  define dll_export
#  define dll_import
# endif
#endif

#if defined(uiQVel_EXPORTS) || defined(UIQVEL_EXPORTS)
# define do_export_uiQVel
#else
# if defined ( do_export_uiQVel )
#  undef do_export_uiQVel
# endif
#endif


#if defined( do_export_uiQVel )
# define Export_uiQVel	dll_export
# define Extern_uiQVel
#else
# define Export_uiQVel	dll_import
# define Extern_uiQVel	dll_extern
#endif

#if defined ( do_import_export )


#include "uiodmaindeps.h"
#include "qveldeps.h"

#endif

#endif
