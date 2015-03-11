/*
Copyright (C) 2012 LuminTerra, LLC All rights reserved.

This file is part of OpendTect and may be used under the terms of:

The GNU General Public License version 3 or higher, as published by
the Free Software Foundation.

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Ver 1.2 JB West 5/2012

*/
#include "uistring.h"
#include "cubesampling.h"
#include "seisread.h"
#include "seiswrite.h"
#include "seisioobjinfo.h"
#include "seisselectionimpl.h"
#include "seistrc.h"
#include "seisbounds.h"
#include "seistrcprop.h"
#include "ioobj.h"

#include "emhorizon3d.h"
#include "emmanager.h"
#include "emobject.h"
#include "emsurfacetr.h"
#include "ioobj.h"
#include "transl.h"
#include "array2dinterpol.h"
#include "arrayndimpl.h"
#include "cubesampling.h"
#include "horsampling.h"
#include "seistrc.h"
#include "survinfo.h"
#include "ioman.h"
#include "bufstring.h"
//#include "tutmod.h"
#ifdef __GNUC__
#include <sys/stat.h>
#else
#include <direct.h>
#endif

namespace EM { class Horizon3D; }
#include "qvelvolumebuilder.h"
#include "qvelseistools.h"

#define mErrRet(_str) \
	errmsg_ = _str; \
	return Executor::ErrorOccurred();

namespace QVel
{
	SeisTools::SeisTools()
		: Executor("QVel Interval Velocity Calculations")
		, inioobj_(0), outioobj_(0)
		, rdr_(0), wrr_(0)
		, trcin_(*new SeisTrc)
		, trcout_(*new SeisTrc)
		, first_(true)
		, volBuilder_(0)
		, stepName_((char *)"Initializing...")
		, calibrationType_(0)
		, calibrationFactor_(1)
		, isvel_(0)
		, gridder_(QVelVolumeBuilder::ThinPlateSpline)

	{
		clear();
		#ifdef __GNUC__
    mkdir("/temp",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#else
    _mkdir("/temp");
#endif
    tracef = fopen("/temp/QVelSeis.log","w");
    if (tracef == NULL) UsrMsg("Cannot create /tmp/QVel.log\n");
    if (tracef != NULL) fprintf(tracef,"BEGIN RUN\n");

	}


	QVel::SeisTools::~SeisTools()
	{
		clear();
		delete &trcin_;
		delete &trcout_;
		delete volBuilder_;
	}


	void QVel::SeisTools::clear()
	{
		delete inioobj_; inioobj_ = 0;
		delete outioobj_; outioobj_ = 0;
		delete rdr_; rdr_ = 0;
		delete wrr_; wrr_ = 0;

		totnr_ = -1; nrdone_ = 0;
		qvelStatus_ = Executor::MoreToDo();
	}

	void QVel::SeisTools::setInput( const IOObj& ioobj )
	{  delete inioobj_; inioobj_ = ioobj.clone(); }

	void QVel::SeisTools::setOutput( const IOObj& ioobj )
	{ delete outioobj_; outioobj_ = ioobj.clone(); }

	void QVel::SeisTools::setHorizons(TypeSet<MultiID> horizons)
	{ horizons_ = horizons; numHorizons_ = horizons.nrItems();}

	//Executor
	const char* QVel::SeisTools::message() const
	{
		const char * b = errmsg_.getFullString();
		return errmsg_.isEmpty() ? "Calculating..." : b;
		return "a";
	}

	// Executor
	od_int64 QVel::SeisTools::totalNr() const
	{
		if (volBuilder_ && qvelStatus_ !=Executor::Finished()) return volBuilder_->totnr_;
		return totnr_;
	}

	//Executor
	od_int64 QVel::SeisTools::nrDone() const
	{
		if (volBuilder_ && qvelStatus_ !=Executor::Finished()) return volBuilder_->nrdone_;
		return nrdone_;
	}
#if 1
	//Simple seismic trace reader
	bool QVel::SeisTools::createReader()
	{

		rdr_ = new SeisTrcReader( inioobj_ );
		Seis::RangeSelData* sd = new Seis::RangeSelData( cs_ );
		rdr_->setSelData( sd );

		rdr_->forceFloatData( true );
		if ( !rdr_->prepareWork() )
		{

			errmsg_ =  rdr_->errMsg();
			return false;
		}

		return true;
	}

	// simple seismic writer
	bool QVel::SeisTools::createWriter()
	{

		wrr_ = new SeisTrcWriter( outioobj_ );
		if (!wrr_ ) return false;

		return true;
	}

	/////////////////////////////////////////////////////

	// Executor interface

	int QVel::SeisTools::nextStep()
	{

		if (getState() == Stop) {
			closeOutput();
			return Executor::Finished();
		}
		if (first_)
		{
			first_ = false;
			nrdone_ = 0;
		}
		if ( !rdr_ )
			if (inioobj_)
				return createReader() ? Executor::MoreToDo()
				: Executor::ErrorOccurred();

		// run the Volume Builder steps until finished
		if ( !volBuilder_ )
		{
			stepName_ = (char *) "Building interval velocities...";
			errmsg_.setEmpty();

			volBuilder_ = new QVelVolumeBuilder(rdr_, wrr_, 
				selectedWells_, 
				markerNames_, 
				horizons_ , 
				numHorizons_, 
				radius_, 
				power_,
				calibrateAll_,
				useMarkers_,
				cs_,
				gridder_,
				relaxation_,
				isTDC_,
				hzmarklist_
				);
			if (!volBuilder_) return Executor::ErrorOccurred();
		}

		// more to do in volBuilder
		if (qvelStatus_ !=Executor::Finished())
		{
			qvelStatus_ = volBuilder_->nextStep();
			if (qvelStatus_ !=Executor::Finished())
				return Executor::MoreToDo();

			// drop thru on done with volBuilder steps
			stepName_ = (char *) "Writing Lines...";

			// progress meter
			nrdone_= 0;
			totnr_ = cs_.nrInl();

		}
		//if (tracef) fprintf(tracef,"Handleing traces\n");
		// Done with interval-based velocity volume building
		// handle input seismic
		if (inioobj_)
		{
			int rv = rdr_->get( trcin_.info() );
			if ( rv < 0 )
			{ 
				errmsg_ = rdr_->errMsg(); 
				closeOutput();
				return Executor::ErrorOccurred(); 
			}
			else if ( rv == 0 )
			{
				if (tracef != NULL) fprintf(tracef,"rv == 0\n");
				closeOutput();
				return Executor::Finished();
			}
			else if ( rv == 1 )
			{
				if ( !rdr_->get(trcin_) )
				{ 
					errmsg_ = rdr_->errMsg(); 
					closeOutput();
					return Executor::ErrorOccurred(); 
				}
			}

			trcout_ = trcin_;
			int trcsize = (int) (trcin_.info().sampling.start/trcin_.info().sampling.step + trcin_.size());
			//int trc2 = trcin_.size();
			trcout_.reSize(trcsize, true );
			
			// set up seismic writer 
			if ( !wrr_ )
			{
				if (!createWriter()) return Executor::ErrorOccurred();
				

					if ( !wrr_->prepareWork(trcout_) )
					{
						errmsg_ = wrr_->errMsg();
						return Executor::ErrorOccurred();
					}

			}
			if (getState() == Stop) {

				closeOutput();
			    	
				return Executor::Finished();
			}

			// Executor handle a trace at a time
			
			handleTrace(isvel_);

		} else { // no input trace, just write
			return outputTraces(wrr_);
		}

		if ( trcin_.info().nr % volBuilder_->ncols_ == 0)
			nrdone_++;

		// input trace in, and now output trace out
		if ( !wrr_->put(trcout_) )
		{ errmsg_ = wrr_->errMsg(); return Executor::ErrorOccurred(); }

		return Executor::MoreToDo();
	}

	// write out velocity model with no input stacking vels to be merged
	int QVel::SeisTools::outputTraces(SeisTrcWriter *wrr)
	{
//#ifdef DEBUG_ODT
//		static FILE *tracef=fopen("/temp/trace2.out","w");
//#else
//		static FILE * tracef = NULL;
//#endif
		// no input trace
		if (getState() == Stop) {
			closeOutput();
			return Executor::Finished();
		}
		{
			if ( !wrr_ )
			{
				wrr_ = new SeisTrcWriter( outioobj_ );
				
			}

			HorSampling inputhrg = cs_.hrg;
			const float zstep = cs_.zrg.step;
			StepInterval<int> outputzrg( mNINT32(cs_.zrg.start/zstep),
				mNINT32(cs_.zrg.stop/zstep),
				mNINT32(cs_.zrg.step/zstep) );
			if ( outputzrg.step<1 ) outputzrg.step = 1;
			int trcsz = cs_.zrg.nrSteps();
			SeisTrc * traceout;
			traceout = new SeisTrc( trcsz ) ;
			traceout->info().sampling.start = cs_.zrg.start; 
			traceout->info().sampling.step = cs_.zrg.step; 
			float dt = traceout->info().sampling.step;

			if ( !wrr_->prepareWork(*traceout) )
			{
				errmsg_ = wrr_->errMsg();
				return Executor::ErrorOccurred();
			}
			float *traceIntVels = new float[traceout->size()+1];
			float *workTrace = new float[traceout->size()+1];
			//if (tracef != NULL) fprintf(tracef," numrows %d ncols %d \n", volBuilder_->nrows_, volBuilder_->ncols_);
			float modelbot = traceout->size() * dt;
			if (isTDC_)
				volBuilder_->propagateWellIntersections(modelbot);
			// get tabularized weights
			Array2DImpl<double>* weightsTable = volBuilder_->getWeightsTable();
			for (int irow = 0; irow < volBuilder_->nrows_; irow++)
			{
				nrdone_++;
				for (int icol = 0; icol < volBuilder_->ncols_; icol++)
				{
					//if (tracef != NULL) fprintf(tracef,"    loop start irow =%d icol = %d \n",irow, icol);
					//if (tracef != NULL) fprintf(tracef,"    loop cs inl %d crl %d \n",volBuilder_->cs_.hrg.inlRange().nrSteps(), volBuilder_->cs_.hrg.crlRange().nrSteps());
					if (getState() == Stop) return Executor::Finished();
					// get interval velocities from model, in trace sampled steps
					volBuilder_->processVelocityIntervalTrace( irow, icol, 
						traceout->size(),
						traceIntVels,
						workTrace,
						dt,
						*traceout, 
						isTDC_,
						weightsTable);
					//if (tracef != NULL) fprintf(tracef,"    loop after process irow =%d icol = %d \n",irow, icol);

					int jrow = irow * (volBuilder_->cs_.hrg.step.inl()) + volBuilder_->cs_.hrg.inlRange().start;
					int jcol = icol * (volBuilder_->cs_.hrg.step.crl()) + volBuilder_->cs_.hrg.crlRange().start;

					BinID currentpos_;
					currentpos_.inl() = jrow;
					currentpos_.crl() = jcol;
					//if (tracef != NULL) fprintf(tracef,"    inl=%d crl = %d \n",jrow, jcol);
					fflush(tracef);
#if 0
					UsrMsg(" irow = ");
					UsrMsg(irow);
					UsrMsg(" jrow = ");
					UsrMsg(jrow);
					UsrMsg(" icol = ");
					UsrMsg(icol);
					UsrMsg(" jcol = ");
					UsrMsg(jcol);
#endif

					//if (tracef != NULL) fprintf(tracef,"icol %d jcol %d irow %d jrow %d\n",icol,jcol,irow,jrow);
					traceout->info().binid = currentpos_;
					traceout->info().coord = SI().transform( currentpos_ );
					trcin_.info().sampling.start = 00;
					// set output trace
					for ( int idx=0; idx<traceout->size(); idx++ )
					{
						traceout->set( idx, traceIntVels[idx], 0);
					}

					// output trace
					if ( !wrr_->put(*traceout) )
					{ errmsg_ = wrr_->errMsg(); return Executor::ErrorOccurred(); }
					traceout->info().nr++;
				}
			}

			if (tracef != NULL) fprintf(tracef,"   close \n");	
			
			closeOutput();

			delete [] traceIntVels;
			delete [] workTrace;
			delete traceout;
			if (tracef != NULL) fprintf(tracef,"   done with clean\n");
			fflush(tracef);
			return Executor::Finished();
		}
	}

	void QVel::SeisTools::closeOutput()
	{
		    if (tracef) fprintf(tracef,"Close output file\n");
			if (wrr_ != NULL)
			    wrr_->close();

			VelocityDesc veldesc(VelocityDesc::Interval);

			veldesc.fillPar( outioobj_->pars());

			if ( !IOM().commitChanges(*outioobj_) ) {
				if (tracef != NULL) fprintf(tracef,"   cannot write vel info\n");
			} else
			    if (tracef != NULL) fprintf(tracef,"   Wrote vel info\n");	

			if (tracef) fflush(tracef);
	}

	// merge incoming stacking velocity trace with computed interval model
	// adjust intervals to match either over entore interval, or well-region weighted.
	// in both cases, well picks are honored.
	void QVel::SeisTools::handleTrace(bool isvel)
	{
		float dt;
		float startTime;
		// get tabularized weights
		Array2DImpl<double>* weightsTable = volBuilder_->getWeightsTable();
		if (inioobj_)
		{
			int irow = trcin_.info().binid.inl();
			int icol = trcin_.info().binid.crl();
			dt = trcin_.info().sampling.step;
			startTime = trcin_.info().sampling.start;
			{
				irow = volBuilder_->cs_.inlIdx(irow);
				icol = volBuilder_->cs_.crlIdx(icol);

				// can't happen
				if (irow < 0 || icol < 0) return;
				// start at 0, end at end of trace
				int trcsize = (int) (trcin_.info().sampling.start/trcin_.info().sampling.step + trcin_.size());
				//int trc2 = trcin_.size();
				// if the input trace is velocity
				if (isvel)
				{
					float *traceIntVels = new float[trcsize+1];
					float *intervalAdjust = new float[volBuilder_->numHorizons()+1];
					float *modelIntVels = new float[trcsize+1];
					float *workTrace =  new float[trcsize+1];

					for (int i =0; i < volBuilder_->numHorizons()+1; i++)
						intervalAdjust[i]=1.0;
					for ( int icomp=0; icomp<trcin_.nrComponents(); icomp++ )
					{
						for ( int idx=0; idx<trcsize; idx++ )
						{
							// adjust trace to start at t=0
							float timeval = idx * dt;
							int ind0 = trcin_.info().sampling.nearestIndex<float>(timeval);
							if (ind0 < 0) ind0 = 0;
							const float v = trcin_.get( ind0, icomp );
							traceIntVels[idx] = v;
						}
						// get the model velocities as a seismic trace
						volBuilder_->processVelocityIntervalTrace( irow, icol, 
							trcsize,
							modelIntVels,
							workTrace,
							dt,
							trcin_,  isTDC_, weightsTable);

						// calculate the adjustment
						volBuilder_->computeIntervalAdjustTrace (irow, icol, 
							trcin_,
							trcsize,
							traceIntVels,
							intervalAdjust);

						// apply the adjustment
						volBuilder_->processVelocityIntervalVolumeTrace( irow, icol, 
							trcsize,
							traceIntVels,
							modelIntVels,
							dt,
							calibrateAll_,
							intervalAdjust,
							trcin_);

						// build output trace
						trcout_.info().sampling.start = 0;
						for ( int idx=0; idx < trcsize; idx++ )
						{
							trcout_.set( idx, traceIntVels[idx], icomp);
						}
					}
					delete[] traceIntVels;    
					delete [] intervalAdjust;
					delete [] modelIntVels;
				} 
				else // not stacking velocities in trcin_, just seismic -- use trace size only
				{
					// just build output trace from intervals
					float *modelIntVels = new float[trcsize+1];
					float * workTrace = new float[trcsize+1];
					volBuilder_->processVelocityIntervalTrace( irow, icol, 
						trcsize,
						modelIntVels,
						workTrace, 
						dt,
						trcin_, isTDC_, weightsTable);
					// fill output trace
					trcout_.info().sampling.start = 0;
					for ( int icomp=0; icomp<trcin_.nrComponents(); icomp++ )
						for ( int idx=0; idx<trcsize; idx++ )
						{
							trcout_.set( idx, modelIntVels[idx], icomp);
						}
						delete[] modelIntVels;
						delete[] workTrace;
				}
			}
		}

	}
#endif

} // namespace