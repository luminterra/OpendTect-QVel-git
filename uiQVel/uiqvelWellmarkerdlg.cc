/*Copyright (C) 2002-2010 dGB Beheer B.V. All rights reserved.

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
________________________________________________________________________

 (C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
 Author:	Nanne Hemstra
 Date:		October 2003
________________________________________________________________________

-*/
static const char* rcsID mUsedVar = "$Id: uiwellmarkerdlg.cc 33461 2014-02-26 06:15:35Z aneesh.tiwari@dgbes.com $";




#include "uibutton.h"
#include "uicolor.h"
#include "uifileinput.h"
#include "uigeninput.h"
#include "uilistbox.h"
#include "uimsg.h"
#include "uistratlvlsel.h"
#include "uistrattreewin.h"
#include "uitable.h"
#include "uitblimpexpdatasel.h"
#include "uitoolbutton.h"

#include "ctxtioobj.h"
#include "file.h"
#include "ioman.h"
#include "ioobj.h"
#include "iopar.h"
#include "oddirs.h"
#include "stratlevel.h"
#include "strmprov.h"
#include "survinfo.h"
#include "tabledef.h"
#include "welldata.h"
#include "wellman.h"
#include "wellimpasc.h"
#include "welltrack.h"
#include "welltransl.h"

#include "uilistbox.h"
#include "uigeninput.h"
#include "uilistbox.h"
#include "uispinbox.h"
#include "uimsg.h"
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

#include "tutmod.h"

#include "qvelvolumebuilder.h"
#include "uiqvelWellmarkerdlg.h"

uiQVelWellMarkerTools::uiQVelWellMarkerTools( uiParent* p, BufferStringSet selectedHorizons, BufferStringSet selectedMarkers, BufferStringSet associatedMarkers )
    : uiDialog( p, Setup( "Well-Marker Selection",
    "Specify Well to Marker association for QVel model",
    "qvel:105.0.3") )

{
    //tracef = fopen("/temp/traceu.out","w");
    //selectedHorizonnms_.add((const BufferStringSet&) selectedHorizons, false);
    //selectedMarkernms_.add((const BufferStringSet&)selectedMarkers, false);
   // associatedMarkers_.add((const BufferStringSet&)associatedMarkers,true);

    hzmark_ = NULL;

    llbw = (uiLabeledListBox*) new uiLabeledListBox( this, "Available Horizons", OD::ChoiceMode::ChooseOnlyOne ,  uiLabeledListBox::AboveMid );
    horizonslist_ = llbw->box();
    
    horizonslist_->selectionChanged.notify( 
	    			mCB(this,uiQVelWellMarkerTools,hznchg) );

    //uiLabel * lab  = new uiLabel(this,"<--Associated-->");
    //lab->attach(ensureRightOf,llbw);

    llbw2 = (uiLabeledListBox*) new uiLabeledListBox(this,"Associated Markers",OD::ChoiceMode::ChooseOnlyOne, uiLabeledListBox::AboveMid);
    assoclist_ = llbw2->box();
    llbw2->attach(ensureRightOf,llbw);
    assoclist_->setAllowDuplicates(true);
    
   

        // select wells
    /*
    assocbut_ = new uiPushButton( this, "Associate", true );
    assocbut_->activated.notify( mCB(this,uiQVelWellMarkerTools,assocButPush) );
    assocbut_->attach(ensureRightOf, llbw2);
    assocbut_->setToolTip("Select wells to be used for markers and/or T/D curves.");
    


    llbw3 = new uiLabeledListBox( this, "Available Markers", true ,uiLabeledListBox::AboveMid);  
    markerslist_ = llbw3->box();
    llbw3->attach(ensureRightOf,llbw2);
    */
     horizonslist_->setEmpty();
    
 //   markerslist_->setEmpty();
    
     assoclist_->setEmpty();

    for (int i = 0; i < selectedHorizons.size(); i++)
	horizonslist_->addItem(selectedHorizons[i]->str(),false,i);
    //for (int i = 0; i < selectedMarkers.size(); i++)
	//markerslist_->addItem(selectedMarkers[i]->str());
    for (int i = 0; i < associatedMarkers.size(); i++)
	assoclist_->addItem(associatedMarkers[i]->str(),false,i);

        horizonslist_->selectionChanged.notify( 
	    			mCB(this,uiQVelWellMarkerTools,hznchg) );
	assoclist_->selectionChanged.notify( 
	   			mCB(this,uiQVelWellMarkerTools,mrkchg) );

	horizonslist_->setPrefHeightInChar(associatedMarkers.size());
	assoclist_->setPrefHeightInChar(associatedMarkers.size());
 
   horizonslist_->chooseAll(false);
   assoclist_->chooseAll(false);
   assoclist_->clear();
    finaliseDone.notify( mCB(this,uiQVelWellMarkerTools,initWin) );

}

uiQVelWellMarkerTools::~uiQVelWellMarkerTools()
{
}

#define mErrRet(s) { uiMSG().error(s); return false; }

bool uiQVelWellMarkerTools::acceptOK( CallBacker* )
{

    return true;
}
void uiQVelWellMarkerTools::initWin( CallBacker *)
{
    UsrMsg("For all horizons,\nselect a horizon, then choose a marker to (dis)associate.");
}

void	uiQVelWellMarkerTools::assocButPush(CallBacker*)
{
}
void uiQVelWellMarkerTools::hznchg( CallBacker* )
{
    TypeSet<int> ahzn;
    ahzn.erase();
    //horizonslist_->getChosen(ahzn);
    horizonslist_->getChosen(ahzn);
    if (ahzn.size()<=0) return;
    hselect_ = ahzn.first();
    
    assoclist_->setChosen(ahzn);
    //horizonslist_->setItemCheckable(ahzn.first(),true);
    //horizonslist_->setItemChecked(ahzn.first(),true);
    
}
void uiQVelWellMarkerTools::mrkchg( CallBacker* )
{
       TypeSet<int> ahzn;
    ahzn.erase();
    
    horizonslist_->getChosen(ahzn);
    if (ahzn.size()<=0)
    {
     uiMSG().warning("First select a horizon, then choose a marker to (dis)associate.");
     assoclist_->clear();
     return;
    }
    int ihzn = ahzn.first();
    
    TypeSet<int> amrk;
    amrk.erase();
    assoclist_->getChosen(amrk);
    if (amrk.size()<=0) {
	uiMSG().warning("not chosen.");
	return;
    }
    int imark=amrk.first();
    
    //markerslist_->setItemCheckable(amrk.first(),true);
    //markerslist_->setItemChecked(amrk.first(),true);
    
   
    
    char * mrk ="";
    
    for(int i = 0; i < assoclist_->size(); i++) {
	
	if (i == imark) {
	    mrk= strdup(assoclist_->textOfItem(i));
	    
	}
    }

    const char * prevt = assoclist_->textOfItem(ihzn);
    

    if (!strcmp(prevt,mrk)) 
    {
	assoclist_->setItemText(ihzn,"");
	
	if (prevt && strlen(prevt) > 0)
	assoclist_->addItem(prevt);
    }
    else
    {
    assoclist_->setItemText(imark,"");
    assoclist_->setItemText(ihzn,mrk);
    
    if (prevt && strlen(prevt) > 0)
       assoclist_->addItem(prevt);
       
    }
    //horizonslist_->setAllItemsChecked(false);
    horizonslist_->chooseAll(false);
    //markerslist_->setAllItemsChecked(false);
    assoclist_->chooseAll(false);
    assoclist_->scrollToTop();
 
}
QVel::HorizonMarkersList * uiQVelWellMarkerTools::gethzmark(int & nentries)
{
    nentries = horizonslist_->size();
    if (hzmark_) delete[] hzmark_;
    hzmark_ = new QVel::HorizonMarkersList[nentries];
    for (int i = 0; i < nentries; i++)
    {
	hzmark_[i].horizons = horizonslist_->textOfItem(i);
	hzmark_[i].assocMarkers =  assoclist_->textOfItem(i);
	
    }
    
    return hzmark_;
}

