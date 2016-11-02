/*
 Copyright (C) 2012 LuminTerra, LLC All rights reserved.

 This file is part of OpendTect and may be used under the terms of:

 The GNU General Public License version 3 or higher, as published by
 the Free Software Foundation.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

 Ver 1.2 JB West 4/2012

*/

#include "emhorizon3d.h"
#include "emmanager.h"
#include "emobject.h"
#include "emsurfacetr.h"
namespace EM { class Horizon3D; }
#include "qvelvolumebuilder.h"
#include "qvelseistools.h"
#include "odplugin.h"
#include "tutmod.h"


mDefODPluginEarlyLoad(QVel)
mDefODPluginInfo(QVel)
{
    mDefineStaticLocalObject( PluginInfo, retpi,(
	"QVel Velocity Model Builder Base",
	"QVel",
	"LuminTerra, LLC",
	"1.2",
    	"Back-end for the LuminTerra QVel plugin."
    	"\nNon-ui QVel velocity model builder computational engine.") );
    return &retpi;
}


mDefODInitPlugin(QVel)
{
    return 0;
}
