#_______________________Pmake___________________________________________________
#
#	Makefile : 	Plugin 'template'
# 	Apr 2012	Bert
#	RCS :		$Id: CMakeLists.txt,v 1.1 2012/04/17 11:29:00 cvsraman Exp $
#_______________________________________________________________________________
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__win__")
SET(OD_MODULE_DEPS uiODMain QVel)
SET(OD_IS_PLUGIN yes)
SET(OD_MODULE_SOURCES
	uiqvelpi.cc
	uiqvelseistools.cc
	uiqvelwelltools.cc
	uiqvelWellmarkerdlg.cc)
SET( OD_NO_ALO_ENTRY 1 )
SET( OD_PLUGIN_ALO_EXEC ${OD_MAIN_EXEC} )
OD_INIT_MODULE()
