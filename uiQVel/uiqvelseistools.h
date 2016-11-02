/*
Copyright (C) 2012 LuminTerra, LLC All rights reserved.

This file is part of OpendTect and may be used under the terms of:

The GNU General Public License version 3 or higher, as published by
the Free Software Foundation.

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Ver 1.1 JB West 4/2012

*/

#ifndef uiqvelseistools_h
#define uiqvelseistools_h

#include "uidialog.h"
#include "uichecklist.h"
#include "uiveldesc.h"
#include "uibutton.h"
#include "uilabel.h"
#include "seistype.h"
#include "uimultisurfaceread.h"
#include "emhorizon3d.h"
#include "qvelseistools.h"
#include "uiqvelWellmarkerdlg.h"

class uiSeisSel;
class uiSeisSubSel;
class uiGenInput;
class CtxtIOObj;

namespace QVel { class uiQVelSeisTools; }
//namespace QVel

class uiQVelSeisTools : public uiDialog
{
public:

    uiQVelSeisTools(uiParent*,Seis::GeomType);
    ~uiQVelSeisTools();

protected:

   // CtxtIOObj&		inctio_;
   // CtxtIOObj&		outctio_;
    Seis::GeomType	geom_;
    QVel::SeisTools&	tool_;
	uiCheckBox*     inpcheck_;
    uiVelSel*		inpfld_;
    uiSeisSubSel*	subselfld_;
    uiSeisSel*		outfld_;
    uiGroup*		scalegrp_;
    uiGenInput*		calibfld_;
    uiGenInput*		factorfld_;
    uiGenInput*		relaxfld_;
    uiGenInput*		powerfld_;
    uiGenInput*		newsdfld_;
    uiGenInput*     usemarkerfld_;
    uiGenInput*	useTDCfld_;
    uiGenInput*     gridderfld_;
    uiLabel *		welllab_;
    uiLabel *		hrzlab_;
    uiPushButton *  hrzbut_;
    uiPushButton *  wellbut_;
    uiPushButton*   mrkbut_;
    uiMultiSurfaceRead*	surfacefld_;
    TypeSet<MultiID> surfaceids_;
    BufferStringSet selectedHorizonnms_;
    BufferStringSet  associatedMarkernms_;
    BufferStringSet markerNames_;
    uiQVelWellTools* wtool_;
    uiQVelWellMarkerTools * wmt_ ;
    QVel::HorizonMarkersList * defaultHzMarkList_;


    bool		acceptOK(CallBacker*);
    void		wellButPush(CallBacker*);
    void		hrzButPush(CallBacker*);
    void		mrkButPush(CallBacker*);
    void		doProc(CallBacker*);
    void		initWin(CallBacker*);
    void		inpSel(CallBacker*);
    void		useMarkTog( CallBacker* b );
    void		useTDCTog( CallBacker* b );
    void        gridderType( CallBacker* b );
	void		inpCheck(CallBacker *);

    float		factor_;
    float       power_;
    CubeSampling *cs_;
    bool		isvel_;
    int			gridder_;
    float		relaxation_;
    bool	    useTDC_;
    
};

#endif
