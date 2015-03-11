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

#ifndef uiqvelwellmarkerdlg_h
#define uiqvelwellmarkerdlg_h

/*+
________________________________________________________________________

(C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
Author:	Nanne Hemstra
Date:		May 2007
RCS:		$Id: uiwellmarkerdlg.h 31685 2013-09-25 14:47:06Z arnaud.huck@dgbes.com $
________________________________________________________________________

-*/

#include "uiwellmod.h"
#include "uidialog.h"
#include "wellmarker.h"

class uiGenInput;
class uiCheckBox;
class uiTable;
namespace Well { class Marker; class Track; class MarkerSet; }
class uiQVelWellMarkerTools : public uiDialog
{
public:
    uiQVelWellMarkerTools(uiParent*);
    uiQVelWellMarkerTools(uiParent*,BufferStringSet selectedHorizons, BufferStringSet selectedMarkers, BufferStringSet associatedMarkers);
    ~uiQVelWellMarkerTools();
    //const TypeSet<MultiID>& getSelectedWells() { return selectedwellids_;}
    const BufferStringSet getMarkerNames() { return associatedMarkers_;}
    int getNrWellsSelected() { return nrWellsSelected_;}
    QVel::HorizonMarkersList * gethzmark(int & nentries);
protected:
    int hselect_;
    uiListBox*      assoclist_;
    uiListBox*	    horizonslist_;
    //uiListBox*	    markerslist_;
    uiPushButton *  assocbut_;
    uiLabeledListBox* llbw;
	uiLabeledListBox* llbw2;
	uiLabeledListBox* llbw3;
    
    int             nrWellsSelected_;
    MultiID		wellid_;
    ObjectSet<IOObj>	wellobjs_;
    //TypeSet<MultiID> selectedwellids_;
    //BufferStringSet selectedHorizonnms_; 
    //BufferStringSet selectedMarkernms_;
    BufferStringSet associatedMarkers_;

    void		inpchg(CallBacker*);
    bool		acceptOK(CallBacker*);
    void        initWin( CallBacker *);
    void		assocButPush(CallBacker*);
    void    hznchg( CallBacker* );
    void    mrkchg( CallBacker* );
    QVel::HorizonMarkersList *hzmark_;
    int nhzmark_;
    FILE * tracef;

};

#endif


