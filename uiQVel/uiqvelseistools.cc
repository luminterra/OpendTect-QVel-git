/*
Copyright (C) 2012 LuminTerra, LLC All rights reserved.

This file is part of OpendTect and may be used under the terms of:

The GNU General Public License version 3 or higher, as published by
the Free Software Foundation.

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Ver 1.1 JB West 4/2012

*/

#include "cubesampling.h"
#include "uiseissel.h"
#include "uigeninput.h"
#include "uiseissubsel.h"
#include "uitaskrunner.h"
#include "uimsg.h"
#include "ioman.h"
#include "seisioobjinfo.h"
#include "seistrctr.h"
#include "seistype.h"
#include "seisselection.h"
#include "ctxtioobj.h"
#include "ioobj.h"
#include "survinfo.h"
#include "wellextractdata.h"
#include "wellmarker.h"
#include "emhorizon3d.h"
#include "emmanager.h"
#include "emobject.h"
#include "emsurfacetr.h"

#include "uiiosel.h"
#include "uiattrdesced.h"
#include "uiattrgetfile.h"
#include "uiattribfactory.h"
#include "uiattrinpdlg.h"
#include "uiattrtypesel.h"
#include "uiautoattrdescset.h"
#include "uitoolbutton.h"
#include "uicombobox.h"
#include "uifileinput.h"
#include "uigeninput.h"
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
#include "uilistbox.h"
#include "seisioobjinfo.h"
#include "veldesc.h"
#include "zdomain.h"
#include "qvelvolumebuilder.h"
#include "qvelseistools.h"
#include "tutmod.h"


namespace EM { class Horizon3D; }
namespace QVel { class SeisTools; }
#include "uiqvelwelltools.h"
#include "uiqvelseistools.h"
#include "uiqvelWellmarkerdlg.h"
#include "qvelseistools.h"

uiQVelSeisTools::uiQVelSeisTools( uiParent* p, Seis::GeomType gt )
    : uiDialog( p, Setup( "QVel Velocity Model builder",
    "Specify QVel Velocity Model parameters",
    "qvel/index.html") )
    , inctio_(*mMkCtxtIOObj(SeisTrc))
    , outctio_(*mMkCtxtIOObj(SeisTrc))
    , geom_(gt)
    , tool_(*new QVel::SeisTools)
    , wtool_(0)
    , factor_(1)
    , isvel_(false)
    , gridder_(QVel::QVelVolumeBuilder::ThinPlateSpline)
    , wmt_(0)

{
    const CallBack inpcb( mCB(this,uiQVelSeisTools,inpSel) );

    subselfld_ = uiSeisSubSel::get( this, Seis::SelSetup(geom_) );

    // The input stacking vels seismic object
    inpfld_ = new uiSeisSel( this, inctio_, uiSeisSel::Setup(geom_) );
    inpfld_->selectionDone.notify( inpcb );
    inpfld_->attach( alignedBelow, subselfld_ );
    inpfld_->setEmpty();
    inpfld_->setCaption("Optional Stacking Velocity Cube");

    // calibrate toggle
    calibfld_ = new uiGenInput( this, "Calibrate input Stacking Velocities",
	BoolInpSpec(true,"By Entire Interval","By Well Weights") );
    calibfld_->attach( alignedBelow, inpfld_ );
    calibfld_->setSensitive(false);

    // horizons
    hrzbut_ = new uiPushButton( this, "Select Horizons", true );
    hrzbut_->activated.notify( mCB(this,uiQVelSeisTools,hrzButPush) );
    hrzbut_->attach(alignedBelow, calibfld_);
    hrzbut_->setToolTip("Select the horizons to be used for Interval Velocity calculations");

    // label
    hrzlab_ = new uiLabel(this,"No horizons selected");
    hrzlab_->attach(centeredRightOf, hrzbut_);



    // select wells
    wellbut_ = new uiPushButton( this, "Select Wells", true );
    wellbut_->activated.notify( mCB(this,uiQVelSeisTools,wellButPush) );
    wellbut_->attach(alignedBelow, hrzbut_);
    wellbut_->setToolTip("Select wells to be used for markers and/or T/D curves.");

    // label
    welllab_ = new uiLabel(this,"No wells selected");
    welllab_->attach(rightOf, wellbut_);


    // use markers toggle
    usemarkerfld_ = new uiGenInput( this, "Horizon Interval Velocity Depth Control",
	BoolInpSpec(true,"From T/D Curves","From Well Markers") );
    usemarkerfld_->attach( alignedBelow, wellbut_ );
    usemarkerfld_->valuechanged.notify( mCB(this,uiQVelSeisTools,useMarkTog) );
    usemarkerfld_->setValue(true);
    usemarkerfld_->setToolTip("Choose source for depth information");

    // horizons
    mrkbut_ = new uiPushButton( this, "Associate Markers with Horizons", true );
    mrkbut_->activated.notify( mCB(this,uiQVelSeisTools,mrkButPush) );
    mrkbut_->attach(alignedBelow, usemarkerfld_);
    mrkbut_->setToolTip("Select the horizon/marker pairs for Interval Velocity calculations");
    mrkbut_->setSensitive(false);

    // use fine T/D curve toggle
    useTDCfld_ = new uiGenInput( this, "Fine Velocity Structure",
	BoolInpSpec(true,"From T/D Curves (slow)","From Seismic Velocites (if selected)") );
    useTDCfld_->attach( alignedBelow, mrkbut_);
    useTDCfld_->valuechanged.notify( mCB(this,uiQVelSeisTools,useTDCTog) );
    useTDCfld_->setValue(false);
    useTDC_ = false;
    useTDCfld_->setToolTip("Choose source for structure");

    gridderfld_ = new uiGenInput( this, "Velocity Gridder",
	BoolInpSpec(true,"Thin Plate Spline","Inverse Distance") );
    gridderfld_->attach( alignedBelow, useTDCfld_);
    gridderfld_->valuechanged.notify( mCB(this,uiQVelSeisTools,gridderType) );
    gridderfld_->setValue(true);
    gridderfld_->setToolTip("Select gridding method for interval velocity layers.");

    // Spline relaxation
    relaxation_ = 0;
    relaxfld_ = new uiGenInput( this, "Spline Flatness",
	FloatInpSpec(relaxation_) );
    relaxfld_->attach( alignedBelow, gridderfld_ );
    relaxfld_->setSensitive(true);
    relaxfld_->setToolTip("Select flatness from none [0] to high [5]");

    // well gridding radius
    factor_ = 10;
    factorfld_ = new uiGenInput( this, "Well influence radius (lines)",
	FloatInpSpec(factor_) );
    factorfld_->attach( alignedBelow, relaxfld_ );
    factorfld_->setToolTip("Select radius in grid cells (lines) where well influence is > 50%");

    // 1/r power for gridder
    power_ = 1;
    powerfld_ = new uiGenInput( this, "1/R^N power",
	FloatInpSpec(power_) );
    powerfld_->attach( alignedBelow, factorfld_ );
    powerfld_->setToolTip("Select power N for 1/(R ** N) falloff [.5--2];\ne.g. 2 means 1/r**2");

    // The output seismic object
    outctio_.ctxt.forread = false;
    outfld_ = new uiSeisSel( this, outctio_, uiSeisSel::Setup(geom_) );
    outfld_->attach( alignedBelow, powerfld_ );
    cs_ = new CubeSampling(true);

}


uiQVelSeisTools::~uiQVelSeisTools()
{
    delete inctio_.ioobj; delete &inctio_;
    delete outctio_.ioobj; delete &outctio_;
    if (wtool_)
	delete wtool_;
   delete &tool_;
}

#define mErrRet(s) { uiMSG().error(s); return false; }

#define mGetHora(varnm,fld) \
    RefMan<EM::EMObject> varnm##_emobj = \
    EM::EMM().loadIfNotFullyLoaded( (fld), &taskrunner ); \
    mDynamicCastGet(EM::Horizon3D*,varnm,varnm##_emobj.ptr()) \
    if ( !varnm ) return false;

#define mGetHor(varnm,fld) \
    varnm = (EM::Horizon3D *)\
    EM::EMM().loadIfNotFullyLoaded( (fld), &taskrunner ); \



bool uiQVelSeisTools::acceptOK( CallBacker* )
{
    // Get cubes and check

    if ( !outfld_->commitInput() )
	mErrRet("Missing Output\nPlease enter a name for the output seismic")

	// horizons
	if (surfaceids_.nrItems() < 1)
	    if (!useTDC_)
		mErrRet("Missing Input\nPlease select at least one input horizon,\n and/or use T/D curves to build model");
    // No Horizons -- we are TDC

    // wells
    if (!wtool_)
	mErrRet("Missing Input\nPlease select at least one input well with optional markers");
    if (wtool_->getNrWellsSelected() < 1)
	mErrRet("Missing Input\nPlease select at least one input well with optional markers");

    // set up all operating params
    tool_.clear();

    // processing range
    CubeSampling  cs;
    subselfld_->getSampling( cs.hrg );
    subselfld_->getZRange( cs.zrg );
    if (cs.zrg.step <19./1000. && useTDC_) 
    {
	int yn = uiMSG().question("Sampling at less than 20ms can be very slow. OK to proceed anyway?",NULL,NULL,NULL,NULL);
	if (!yn)
	    mErrRet("Model z step is too small for this type of model building.");
    }

    if (cs.zrg.start !=0 )
    {
	cs.zrg.start = 0.;
	uiMSG().warning("Model must start at t=0, start time has been set to 0.");
    }
    
    tool_.setRange( cs );

    // input file (optional)
    if ( inpfld_->commitInput() )
    {
	if (!inpfld_->isEmpty())
	    tool_.setInput( *inctio_.ioobj );
    }

    tool_.setIsVel(isvel_);

    // output
    tool_.setOutput( *outctio_.ioobj );
    bool applyall = calibfld_->getBoolValue();
    tool_.setCalibApplyAll(applyall);

    // gridder
    gridderType(NULL);
    tool_.setGridder(gridder_);
    // well factor
    float radius = factorfld_->getfValue();
    if ( mIsUdf(radius) ) radius = 10;
    tool_.setRadius( radius );

    // gridding factor
    float powfactor = powerfld_->getfValue();
    if ( mIsUdf(powfactor) ) powfactor = 1;
    tool_.setPower( powfactor);

    tool_.setHorizons(surfaceids_);

    // wells
    tool_.useMarkers(!usemarkerfld_->getBoolValue());
    tool_.setMarkerNames(wtool_->getMarkerNames());

    tool_.setWellIDs(wtool_->getSelectedWells());
    int n;
#if 1
    if (wmt_ != NULL)
	tool_.setHznMark(wmt_->gethzmark(n));
    else
	tool_.setHznMark(defaultHzMarkList_);
#endif

    tool_.setTDC(useTDC_);
    uiTaskRunner taskrunner( this );
    tool_.enableWorkControl(true);
    /*bool retval = */ taskrunner.execute( tool_ );
    uiMSG().warning(tool_.volBuilder_->getMessages());
    //return retval;
    // always indicate done with gui
	uiMSG().message("QVel job completed.\n");
    return true;
}


void uiQVelSeisTools::inpSel( CallBacker* )
{
    if ( !inpfld_->commitInput() || !inctio_.ioobj ) return;
    {
	VelocityDesc veldesc;
	isvel_ = veldesc.usePar( (*inctio_.ioobj).pars() ) &&
	    veldesc.isVelocity();
	
	const bool hasveldesc = veldesc.usePar( (*inctio_.ioobj).pars() );
	

	if (!isvel_) {
	    uiMSG().warning("Input seismic isn't tagged as a velocity type.\n","I will assume it to be Interval Velocity.");
	    isvel_ = true;
	}
	else if ( isvel_ && veldesc.type_ != VelocityDesc::Interval)
	{
	    uiMSG().warning("Input seismic is Velocity, but isn't Interval velocity.\n","Sorry, I can't convert velocity types.");
	    isvel_ = false;
	}
	bool velintime = ZDomain::isTime(( *inctio_.ioobj).pars() );
	if (isvel_ && !velintime)
	{
	    uiMSG().warning("Input seismic velocity isn't Time domain.\n","Only the trace extents will be used.");
	    isvel_ = false;
	}

    }
    calibfld_->setSensitive(isvel_);
}
void uiQVelSeisTools::hrzButPush( CallBacker* b )
{

    uiMultiSurfaceReadDlg hdlg_( this, "Horizon" );
    if ( !hdlg_.go() ) return;

    hdlg_.iogrp()->getSurfaceIds( surfaceids_ );

    uiTaskRunner taskrunner( this );

    // surfaces
    // lets' pre-load them now
    for (int ihorz=0; ihorz < surfaceids_.nrItems(); ihorz++)
    {
	EM::Horizon3D* hzn;
	mGetHor(hzn, surfaceids_[ihorz]);
	// hzn goes out of scope & is unref'd
    }
    char txt[120];
    sprintf(txt,"%d horizons selected.",(int)surfaceids_.nrItems());
    hrzlab_->setText(txt);
}


void uiQVelSeisTools::wellButPush( CallBacker* b )
{
    wtool_ = new uiQVelWellTools(this);
    wtool_->go();	
    char txt[120];
    sprintf(txt,"%d wells selected.",wtool_->getNrWellsSelected());
    welllab_->setText(txt);
}

void uiQVelSeisTools::mrkButPush(CallBacker*)
{
#if 1
    wmt_ = new uiQVelWellMarkerTools(this,selectedHorizonnms_,wtool_->getMarkerNames(),associatedMarkernms_);
    wmt_->go();
#endif
}
void uiQVelSeisTools::gridderType( CallBacker* b )
{
    // test
    if (gridderfld_->getBoolValue())
    {
	gridder_ = QVel::QVelVolumeBuilder::ThinPlateSpline;
	relaxfld_->setSensitive(true);

    }
    else
    {
	gridder_ = QVel::QVelVolumeBuilder::InverseDx;
	relaxfld_->setSensitive(false);
    }

}
void uiQVelSeisTools::useMarkTog( CallBacker* b )
{
    if (usemarkerfld_->getBoolValue())
    {
	//useTDCfld_->setValue(true);
	useTDCfld_->setSensitive(true);
	mrkbut_->setSensitive(false);
    } else {
	if (useTDCfld_->getBoolValue())
	    uiMSG().warning("I can't do that, please create a seismic volume with T/D curve structure\nand then use that as input for marker-based intervals.");
	useTDCfld_->setValue(false);
	useTDCfld_->setSensitive(false);
	mrkbut_->setSensitive(true);
    }
    if (!usemarkerfld_->getBoolValue())
    {
	if (surfaceids_.nrItems() < 1)
	{
	    uiMSG().error("Missing Input\nPlease select at least one input horizon.");
	    usemarkerfld_->setValue(true);
	    return;
	}
	if (!wtool_)
	{
	    uiMSG().error("Missing Input\nPlease select at least one input well with optional markers");
	    usemarkerfld_->setValue(true);
	    return;
	}
	if (wtool_->getNrWellsSelected() < 1)
	{
	    uiMSG().error("Missing Input\nPlease select at least one input well with optional markers");
	    usemarkerfld_->setValue(true);
	    return;
	}
	if (wtool_->getMarkerNames().size() < 1 )
	{
	    uiMSG().error("Missing Input\nPlease select wells with markers present");
	    usemarkerfld_->setValue(true);
	    return;
	}
	// match horizons to markers here

	uiTaskRunner taskrunner( this );
    int ngood = 0;
	markerNames_ = wtool_->getMarkerNames();
	BufferString emsg;
	bool haveError = false;
	selectedHorizonnms_.erase();
	associatedMarkernms_.erase();
	defaultHzMarkList_ = new QVel::HorizonMarkersList[surfaceids_.nrItems()];
	for (int ihorz = 0; ihorz < surfaceids_.nrItems(); ihorz++)
	{
	    EM::Horizon3D* hzn;
	    mGetHor(hzn, surfaceids_[ihorz]);
	    BufferString hname = hzn->name();
	    selectedHorizonnms_.add(hname);
	    //bool hit = false;
	    if (!markerNames_.isPresent(hname))			
	    {
		haveError = true;
		emsg += "Horizon <";
		emsg += hname;
		emsg += "> has no matching marker(s)\n";
		associatedMarkernms_.add("");
	    } else {
		associatedMarkernms_.add(hname);
		int idx = markerNames_.indexOf(hname);
		if (idx >=0) markerNames_.removeSingle(idx,true);
		defaultHzMarkList_[ihorz].horizons = hname.str();
		defaultHzMarkList_[ihorz].assocMarkers = hname.str();
		ngood++;
	    }
	    // hzn goes out of scope & is unref'd
	}
	for (int i=0; i < markerNames_.size(); i++)
	    associatedMarkernms_.add(markerNames_.get(i));
	if (haveError) {
	    uiMSG().warning("Missing horizon -> well marker matches:\n",emsg);
	    mrkButPush(b);
	} else {
	    emsg += ngood;
	    emsg += " horizons have been matched\n";
	    emsg += "with same-named well markers.";
	    uiMSG().message("All horizons have corresponding well markers.\n",emsg);
	}

    }
}

void uiQVelSeisTools::useTDCTog( CallBacker* b )
{
    useTDC_ = useTDCfld_->getBoolValue();
    // gridding is 1/r only
    if (useTDC_)
    {
	gridderfld_->setValue(false);
	gridderfld_->setSensitive(false);
	inpfld_->setSensitive(false);
    } else {
	gridderfld_->setSensitive(true);
	inpfld_->setSensitive(true);
    }
}
