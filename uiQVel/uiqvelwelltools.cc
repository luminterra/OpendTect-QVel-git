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
#include "uigeninput.h"
#include "uilistbox.h"
#include "uispinbox.h"
#include "uimsg.h"
#include "uicombobox.h"
#include "uifileinput.h"
#include "uigeninput.h"
#include "uibutton.h"
#include "uiioobjsel.h"
#include "uilistbox.h"
#include "uilabel.h"
#include "uimenu.h"
#include "uimsg.h"
#include "uiselsimple.h"
#include "uiseparator.h"
#include "uisplitter.h"
#include "uistoredattrreplacer.h"
#include "uitextedit.h"
#include "uitoolbar.h"
#include "uitaskrunner.h"
#include "bufstring.h"
#include "ioobj.h"
#include "ioman.h"
#include "strmdata.h"
#include "strmprov.h"
#include "iostrm.h"
#include "welldata.h"
#include "wellextractdata.h"
#include "wellio.h"
#include "welllog.h"
#include "welllogset.h"
#include "wellman.h"
#include "wellreader.h"
#include "wellmarker.h"
#include "uitaskrunner.h"
#include "uiqvelwelltools.h"
#include "qvelvolumebuilder.h"
#include "uiqvelWellmarkerdlg.h"

#include "tutmod.h"


uiQVelWellTools::uiQVelWellTools( uiParent* p )
    : uiDialog( p, Setup( "Well Selection",
    "Specify Wells for QVel model",
    "qvel:105.0.3") )
    , nrWellsSelected_(0)


{
    llbw_ =  new uiLabeledListBox( this, "Wells", OD::ChooseAtLeastOne );
    wellsfld_ = llbw_->box();
  
    markernms_.erase();
    mset_=NULL;;
    selectedWells_.erase() ;
    wellid_.setEmpty();
    wellobjs_.erase();
    selectedwellids_.erase();
    selectedwellnms_.erase();
    
    finaliseDone.notify( mCB(this,uiQVelWellTools,initWin) );
}

uiQVelWellTools::~uiQVelWellTools()
{
    markernms_.erase();
    
    selectedWells_.erase() ;
    wellid_.setEmpty();
    wellobjs_.erase();
    selectedwellids_.erase();
    selectedwellnms_.erase();
    // crashes for some reason in well.dll
    //delete mset_;
    //mset_=NULL;
    
}

#define mErrRet(s) { uiMSG().error(s); return false; }

bool uiQVelWellTools::acceptOK( CallBacker* )
{

    nrWellsSelected_ = wellsfld_->nrChosen();
    for ( int idx=0; idx<wellsfld_->size(); idx++ )
    {
	if ( wellsfld_->isChosen(idx) )
	{
	    selectedwellids_.add( wellobjs_[idx]->key() );
	    selectedwellnms_.add( wellobjs_[idx]->name() );
	}
    }
    return true;
}
void uiQVelWellTools::initWin( CallBacker *)
{
    wellsfld_->setEmpty();

    uiTaskRunner tr( this);
    Well::InfoCollector wic;
    if ( !tr.execute(wic) ) return;

    markernms_.erase();
    mset_ = new Well::MarkerSet();
    deepErase( wellobjs_ );
    for ( int iid=0; iid<wic.ids().size(); iid++ )
    {
	IOObj* ioobj = IOM().get( *wic.ids()[iid] );
	if ( !ioobj ) continue;
	wellobjs_ += ioobj;
	wellsfld_->addItem( ioobj->name() );

	const Well::MarkerSet& mrkrs = *wic.markers()[iid];
	for ( int imrk=0; imrk<mrkrs.size(); imrk++ ) {   
	    markernms_.addIfNew( mrkrs[imrk]->name() );
	    mset_->addIfNew((Well::Marker *)mrkrs[imrk]);
	}
    }
}

