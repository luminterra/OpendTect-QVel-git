/*
Copyright (C) 2012 LuminTerra, LLC All rights reserved.

This file is part of OpendTect and may be used under the terms of:

The GNU General Public License version 3 or higher, as published by
the Free Software Foundation.

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Ver 1.1 JB West 4/2012

*/

#include "uilistbox.h"

#include "uimenu.h"
#include "uimsg.h"
#include "uiodmenumgr.h"
#include "uivismenuitemhandler.h"
#include "uivispartserv.h"
#include "viswelldisplay.h"

#include "ioman.h"
#include "ioobj.h"
#include "ptrman.h"
#include "seistype.h"
#include "survinfo.h"

#include "odplugin.h"
#include "uiqvelwelltools.h"
#include "uiqvelseistools.h"
#ifdef uiQVel_EXPORTS
//#pragma message("uiQVel_EXPORTS")
#define TUT_EXPORTS 1
#endif
#include "uitutmod.h"
#include "tutmod.h"
static const int cQVelIdx = -1101;


mDefODPluginInfo(uiQVel)
{

    mDefineStaticLocalObject( PluginInfo, retpi,(\
	"QVel Velocity Model Builder",
	"QVel",
	"LuminTerra, LLC",
	"1.0",
	"Quick post-stack seismic velocity model builder."
	"\nFrom Horizons, T/D curves/markers and stacking velocities") );
	
   
	return &retpi;
}


class uiQVelMgr :  public CallBacker
{
public:
    uiQVelMgr(uiODMain*);

    uiODMain*	appl_;
    uiMenu*	mnuseis_;
    

    void		doSeis(CallBacker*);
    void		do2DSeis(CallBacker*);
    void		do3DSeis(CallBacker*);
    void		launchDialog(Seis::GeomType);

};


uiQVelMgr::uiQVelMgr( uiODMain* a )
    : appl_(a)
{
//    uiMenu* mnu = new uiMenu( appl_, "&QVel Velocity Model builder ..." );
    uiAction  *mnu = new uiAction("&QVel Velocity Model builder ...",
	mCB(this,uiQVelMgr,do3DSeis),"horvelanalwinicon.png");

    appl_->menuMgr().analMnu()->insertItem( mnu );

}
void uiQVelMgr::do3DSeis( CallBacker* )
{ launchDialog( Seis::Vol ); }


void uiQVelMgr::do2DSeis( CallBacker* )
{ launchDialog( Seis::Line ); }

void uiQVelMgr::doSeis( CallBacker* )
{ launchDialog( SI().has2D() ? Seis::Line : Seis::Vol ); }


void uiQVelMgr::launchDialog( Seis::GeomType tp )
{
    uiQVelSeisTools dlg( appl_, tp );
    dlg.go();
}




mDefODInitPlugin(uiQVel)
{
    static uiQVelMgr* mgr = 0; if ( mgr ) return 0;
    mgr = new uiQVelMgr( ODMainWin() );

    //uiQVelAttrib::initClass();
    return 0;
}
