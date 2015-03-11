/*
 Copyright (C) 2012 LuminTerra, LLC All rights reserved.

 This file is part of OpendTect and may be used under the terms of:

 The GNU General Public License version 3 or higher, as published by
 the Free Software Foundation.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

 Ver 1.1 JB West 4/2012

*/

#ifndef qvelseistools_h
#define qvelseistools_h

#include "qvelvolumebuilder.h"
#include "executor.h"
#include "cubesampling.h"
#include "samplingdata.h"
#include "multiid.h"
#include "bufstringset.h"


class IOObj;
class SeisTrc;
class SeisTrcReader;
class SeisTrcWriter;
class QVelVolumeBuilder;

namespace QVel
{

class dll_export SeisTools : public Executor
{
public:
    			SeisTools();
    virtual		~SeisTools();
    void		clear();

    const IOObj*	input() const		{ return inioobj_; }
    const IOObj*	output() const		{ return outioobj_; }
	
    void		setInput(const IOObj&);
    void		setOutput(const IOObj&);
	void        setHorizons(TypeSet<MultiID> horizons);
	void		setPower(float p) { power_ = p;}
	void		setRadius(float r) { radius_ = r;}
    void        setCalibApplyAll(bool b){ calibrateAll_ = b; }
	void		setMarkerNames(const BufferStringSet& names) {markerNames_ = names;}
	void        useMarkers(bool u) { useMarkers_ = u; }
	void		setWellIDs( const TypeSet<MultiID> names) {selectedWells_ = names;}
	void	    setHznMark(HorizonMarkersList *hzn) { hzmarklist_ = hzn;}
	void	    setTDC(bool istdc) { isTDC_ = istdc; }
	void		setIsVel(bool isvel) { isvel_ = isvel;}
	void        setRange( const CubeSampling& cs ) { cs_ = cs; }
	void		setGridder (int grid ) { gridder_ = grid; }
	void		setRelaxation (float relax) { relaxation_ = relax; }
	void		closeOutput();

    
			// Executor compliance functions
    const char*		message() const;
    od_int64		nrDone() const;		
    od_int64		totalNr() const;
    const char*		nrDoneText() const	{ return stepName_; }
			// This is where it actually happens
    int			nextStep();

	QVelVolumeBuilder *volBuilder_;

protected:

	int outputTraces(SeisTrcWriter *wrr);
	int handleTrace(int);
    IOObj*		inioobj_;
    IOObj*		outioobj_;
	TypeSet<MultiID> horizons_;
	int			numHorizons_;
    
    float		power_;
	float		radius_;
	bool        calibrateAll_;
	bool	    isTDC_;

    SeisTrcReader*	rdr_;
	bool			isvel_;
    SeisTrcWriter*	wrr_;
    SeisTrc&		trcin_;
    SeisTrc&        trcout_;
    int			nrdone_;
    mutable int		totnr_;
    uiString	errmsg_;
	BufferStringSet markerNames_;
    TypeSet<MultiID> selectedWells_ ;
    bool		createReader();
    bool		createWriter();
    void		handleTrace(bool);
	bool		first_;
	bool		useMarkers_;
	int			qvelStatus_;
	char *		stepName_;
	int			calibrationType_;
	float		calibrationFactor_;
	CubeSampling cs_;
	int			gridder_;
	float		relaxation_;
	FILE * tracef;
	HorizonMarkersList *hzmarklist_;

};

} // namespace

#endif
