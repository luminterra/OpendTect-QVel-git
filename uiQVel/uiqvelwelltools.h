/*
Copyright (C) 2012 LuminTerra, LLC All rights reserved.

This file is part of OpendTect and may be used under the terms of:

The GNU General Public License version 3 or higher, as published by
the Free Software Foundation.

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Ver 1.1 JB West 4/2012

*/

#ifndef uiqvelwelltools_h
#define uiqvelwelltools_h

#include "uidialog.h"
#include "multiid.h"
#include "uiwellmod.h"
#include "uidialog.h"
#include "wellmarker.h"


class uiGenInput;
class uiLabeledListBox;
class uiLabeledSpinBox;
namespace QVel { class LogTools; }
namespace Well { class Data; }


class uiQVelWellTools : public uiDialog
{
public:

    uiQVelWellTools(uiParent*);
    ~uiQVelWellTools();
    const TypeSet<MultiID>& getSelectedWells() { return selectedwellids_;}
    const BufferStringSet getMarkerNames() { return markernms_;}
    //Well::MarkerSet *getMarkers() { return mset_;}
    int getNrWellsSelected() { return nrWellsSelected_;}
    
protected:
    Well::MarkerSet  *mset_;
    BufferStringSet markernms_;
    uiListBox*      wellsfld_;
    uiLabeledListBox* llbw_;  
    BufferStringSet selectedWells_ ;
    int             nrWellsSelected_;
    MultiID		wellid_;
    ObjectSet<IOObj>	wellobjs_;
    TypeSet<MultiID> selectedwellids_;
    BufferStringSet selectedwellnms_;

    void		inpchg(CallBacker*);
    bool		acceptOK(CallBacker*);
    void        initWin( CallBacker *);

 
    
};


#endif
