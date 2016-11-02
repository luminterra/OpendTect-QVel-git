/*
Copyright (C) 2012 LuminTerra, LLC All rights reserved.

This file is part of OpendTect and may be used under the terms of:

The GNU General Public License version 3 or higher, as published by
the Free Software Foundation.

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Ver 1.1 JB West 4/2012

*/

#include <list>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include "executor.h"
#include "binidvalset.h"
#include "ioman.h"
#include "ioobj.h"
#include "seisbounds.h"
#include "seisread.h"
#include "seisselectionimpl.h"
#include "seistrc.h"
#include "seistrctr.h"
#include "seiswrite.h"
#include "seispacketinfo.h"
#include "sorting.h"
#include "varlenarray.h"
#include "velocitycalc.h"
#include "array2dinterpol.h"
#include "arrayndimpl.h"
#include "binidsurface.h"
#include "binidvalset.h"
#include "datapointset.h"
#include "ioobj.h"
#include "emobject.h"
#ifdef __GNUC__
#include <sys/stat.h>
#else
#include <direct.h>
#endif
#include <stdlib.h>
#include <stdio.h>

#include "mywellhorpos.h"

#undef V43
#ifdef V43
#include "mywellhorpos.h"
#else
#include "wellhorpos.h"
#endif

#include "qvelvolumebuilder.h"
#include "TpsGrid.h"
#include "survinfo.h"

static int debugCol_ = -1 ; //30;
typedef struct {
	float intV;
	float weight;
	int icol;
	int irow;
} markerStruct_;



namespace QVel
{
	QVelVolumeBuilder::QVelVolumeBuilder(SeisTrcReader* rdr, SeisTrcWriter* wrr, 
		TypeSet<MultiID> selectedWells, 
		BufferStringSet markerNames,
		TypeSet<MultiID> horizons, 
		int nrHorizons,
		float radius,
		float power,
		bool calibrateAll,
		bool useMarkers,
		CubeSampling cs,
		int griddertype,
		float relaxation, 
		bool isTDC,
		HorizonMarkersList *hzmarkersList)
		: Executor("QVel")
		, ihorz_(0)
		, topH_(0)
		, botH_(0)
		, nrTraces_(0)
		, weights_(0)
		, counts_(0)
		, firstTime_(true)
		, timeDatum_(0)
		, depthDatum_(0)
		, nrdone_(0)
		, totnr_(0)
		, isTDC_(false)

	{
		rdr_ = rdr;
		wrr_ = wrr;
		selectedWells_ = selectedWells;
		markerNames_ = markerNames;
		horizons_ = horizons;
		numHorizons_ = nrHorizons;
		hzmarkersList_ = hzmarkersList;
		radius_ = radius;
		power_ = power;
		calibrateAll_ = calibrateAll;
		useMarkers_ = useMarkers;
		intVels_ = new  Array2DImpl<double> *[numHorizons_ + 1] ;
		horizonTimes_ = new  Array2DImpl<float> *[numHorizons_ + 1] ;
		rawHorizonTimes_ = new  Array2DImpl<float> *[numHorizons_ + 1] ;
		cs_ = cs;
		gridder_ = griddertype;
		relaxation_ = relaxation;
		messages_ = "";
#ifdef __GNUC__
		mkdir("/temp",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#else
		_mkdir("/temp");
#endif
		tracef = fopen("/temp/QVel.log","w");
		if (tracef == NULL) UsrMsg("Cannot create /tmp/QVel.log\n");
		if (tracef != NULL) fprintf(tracef,"BEGIN QVEL RUN\n");
		if (tracef != NULL) fflush(tracef);
		isTDC_ = isTDC;
		zeroThicknessFound = false;
		unrealVelFound = false;
	}

	QVelVolumeBuilder::~QVelVolumeBuilder()
	{
		for (int i = 0; i < numHorizons_; i++)
		{
			delete intVels_[i];
			delete horizonTimes_[i];
			delete rawHorizonTimes_[i];
		}
		delete [] intVels_;
		delete [] horizonTimes_;
		delete [] rawHorizonTimes_;
		delete [] horizonSortV_;
		delete counts_;
		delete weights_;
		delete botH0_; // initial 0 time layer
	}

	// normalize the 1/R weights to 1.0
	void QVelVolumeBuilder::normalizeVelGrid(Array2DImpl<double>* weights,
		Array2DImpl<float>* counts,
		Array2DImpl<double>* intVels,
		int nrows,
		int ncols)
	{
		for (int irow = 0; irow < nrows; irow++)
		{
			for (int icol = 0; icol < ncols; icol++)
			{
				double vel = intVels->get(irow, icol);
				double w = weights->get(irow, icol);
				vel /= w;
				intVels->set(irow, icol, vel);
			}
		}
	}

#define mGetHor(varnm,fld) \
	varnm = (EM::Horizon3D *)\
	EM::EMM().loadIfNotFullyLoaded( (fld), 0 ) 

#define V43
#ifdef V42

	// OpenDTect V 4.2 methods
	float  QVelVolumeBuilder::intersectWell(MultiID well, int &irow, int &icol, float &zval)
	{
		irow = icol = -1;

		Well::Data* wd = Well::MGR().get( well, false );

		Well::Track timetrack = wd->track();
		// in 4.6 this wrecks it timetrack.toTime( *wd->d2TModel() );

		WellHorPos whpos( timetrack );

		EM::ObjectID emid = EM::EMM().getObjectID( horizons_[horizonSortV_[ihorz_]] );
		whpos.setHorizon( emid );
		EM::Horizon3D* tempH;
		mGetHor(tempH, horizons_[horizonSortV_[ihorz_]]);
		StepInterval<int> rowrg = tempH->geometry().rowRange();
		StepInterval<int> colrg = tempH->geometry().colRange();
		BinIDValueSet bivset(2,false);
		whpos.intersectWellHor( bivset );
		if ( bivset.nrVals() )
		{
			BinIDValueSet::Pos pos = bivset.getPos( 0 );
			BinID bid;

			bivset.get( pos, bid, zval );
			// row/col in irl/crl
			// adjust row/col to 0-start
			irow = bid.inl - cs_.hrg.start.inl;
			icol = bid.crl - cs_.hrg.start.crl;
			// and adjust to possibly decimated subset cs_
			irow/= cs_.hrg.step.inl;
			icol/= cs_.hrg.step.crl;
			float prevdah = mUdf(float);
			float z = wd->track().getDahForTVD( (zval * 1000.), prevdah );  
			return z;
		}
		return mUdf(float);
	}

#else 

	void QVelVolumeBuilder::initWellTracks()
	{
		MultiID well;
		wd_ = new Well::Data*[numWells_];
		timetrack_ = new Well::Track[numWells_];
		for (int iwell=0; iwell < numWells_; iwell++)
		{
			well = selectedWells_[iwell];
			wd_[iwell] = Well::MGR().get( well );
			timetrack_[iwell] = wd_[iwell]->track();
		}
	}
	// get Z and row/col at time tval;
	double  QVelVolumeBuilder::getWellZ(MultiID well, float & tval, int &irow, int &icol, float dt)
	{
		irow = icol = -1;

		Well::Data* wd = Well::MGR().get( well );

		Well::Track timetrack = wd->track();
		//in 4.6 this wrecks it timetrack.toTime( *wd->d2TModel() ,timetrack);

		double z = wd->d2TModel()->getDepth( tval, timetrack );
		//if (tracef != NULL) fprintf(tracef,"depth for tvle = %g\n",z);

		//int num = timetrack.nrPoints();
		//for (int j=0; j < num; j++) 
		//if (tracef != NULL) fprintf(tracef,"for J %d pos %g %g %g value %g\n",j,timetrack.pos(j).x, timetrack.pos(j).y, timetrack.pos(j).z, timetrack.dah(j));
		//float z2 = wd->d2TModel()->getDah( tval*1000., timetrack );
		//float zd = wd->d2TModel()->getDepth( tval, timetrack );

		//double v = wd->d2TModel()->getVelocityForTwt(tval, timetrack);
		//if (tracef != NULL) fprintf(tracef,"velocity is %g\n",v);

		//if (tracef != NULL) fprintf(tracef, "for tval %g z %g z2 %g zd %g \n",tval,z,z2,zd);
		//   const Coord3& crdz2 = timetrack.getPos( z2 );
		//if (tracef != NULL) fprintf(tracef,"crd z2 %g \n",crdz2.z);
		//   const Coord3 crdzd  = timetrack.getPos( zd );
		//if (tracef != NULL) fprintf(tracef,"crd zd %g \n",crdzd.z);

		const Coord3 & crd = timetrack.getPos( z );
		//if (tracef != NULL) fprintf(tracef,"crd %g \n",crd.z);
		const BinID& bid = SI().transform( crd );
		// row/col in irl/crl
		// adjust row/col to 0-start
		irow = bid.inl() - cs_.hrg.start.inl();
		icol = bid.crl() - cs_.hrg.start.crl();
		// and adjust to possibly decimated subset cs_
		irow/= cs_.hrg.step.inl();
		icol/= cs_.hrg.step.crl();

		return z;
	}
	// get Z and row/col at time tval;
	double  QVelVolumeBuilder::getWellZ(int iwell, float & tval, int &irow, int &icol, float dt)
	{
		irow = icol = -1;

		Well::Data* wd = wd_[iwell];

		Well::Track timetrack = timetrack_[iwell];
		//in 4.6 this wrecks it timetrack.toTime( *wd->d2TModel() ,timetrack);

		double z = wd->d2TModel()->getDepth( tval, timetrack );

		const Coord3 & crd = timetrack.getPos( z );
		//if (tracef != NULL) fprintf(tracef,"%g \n",crd.z);
		const BinID& bid = SI().transform( crd );
		//if (tracef != NULL) fprintf(tracef,"bid inl %d crl %d\n",bid.inl, bid.crl);
		// row/col in irl/crl
		// adjust row/col to 0-start
		irow = bid.inl() - cs_.hrg.start.inl();
		icol = bid.crl() - cs_.hrg.start.crl();
		// and adjust to possibly decimated subset cs_
		irow/= cs_.hrg.step.inl();
		icol/= cs_.hrg.step.crl();

		return z;
	}
	// Intersect well with current global active horizon
	float  QVelVolumeBuilder::intersectWell(MultiID well, int &irow, int &icol, float &tval, Array2DImpl<float>*  horz)
	{
		irow = icol = -1;
		//if (tracef != NULL) fprintf(tracef,"intersectwell\n");
		Well::Data* wd = Well::MGR().get( well );

		Well::Track timetrack = wd->track();
		//in 4.6 this wreck is timetrack.toTime( *(wd->d2TModel()) ,timetrack);

		tval = mUdf(float);    

		WellHorIntersectFinder whpos( timetrack , wd->d2TModel());

		EM::ObjectID emid = EM::EMM().getObjectID( horizons_[horizonSortV_[ihorz_]] );
		whpos.setHorizon( emid );
		EM::Horizon3D* tempH;
		mGetHor(tempH, horizons_[horizonSortV_[ihorz_]]);

		tval = whpos.findZIntersection();
		//if (tracef != NULL) fprintf(tracef,"found tz inter %g \n",tval);
		// if no intersection, this could be a truncation 
		const Coord3& crd0 = timetrack.getPos( 0. );
		const BinID& bid0 = SI().transform( crd0 );
		if (mUdf(float)==tval) 
		{
			//if (tracef != NULL) fprintf(tracef,"null intersection horiz %d\n",ihorz_);
			if (horz != NULL)
			{
				tval = horz->get(bid0.inl(), bid0.crl());
				//if (tracef != NULL) fprintf(tracef,"   null intersection fix %g\n",tval);
			}
		}

		//float z = wd->d2TModel()->getDah( tval*1000., timetrack );
		//timetrack = wd->track();
		float z = wd->d2TModel()->getDepth( tval, timetrack );
		//float zd2 = wd->d2TModel()->getDepth( tval );

		//float prevdah = mUdf(float);
		//float z2 = wd->track().getDahForTVD( tval*1000., prevdah );  
		//if (tracef != NULL) fprintf(tracef,"getDah = for tval %g = zd %f, %f tvd %f\n",tval, zd, z, z2);

		const Coord3& crd = timetrack.getPos( z );
		//if (tracef != NULL) fprintf(tracef,"getpos z=%g crd.z=%g \n",z,crd.z);

		//if (tracef != NULL) fprintf(tracef,"track nr %d \n",timetrack.nrPoints());
		//for (int kk = 0; kk < timetrack.nrPoints(); kk++) {
		//    if (tracef != NULL) fprintf(tracef,"   value %d = %g \n",kk,timetrack.value(kk));
		// }
		// if (tracef != NULL) fprintf(tracef,"track zistime %d\n",timetrack.zIsTime());

		const BinID& bid = SI().transform( crd );
		// row/col in irl/crl
		// adjust row/col to 0-start
		irow = bid.inl() - cs_.hrg.start.inl();
		icol = bid.crl() - cs_.hrg.start.crl();
		// and adjust to possibly decimated subset cs_
		irow/= cs_.hrg.step.inl();
		icol/= cs_.hrg.step.crl();
		//if (tracef != NULL) fprintf(tracef,"intersectwell %d %d\n",irow,icol);
		if (mIsUdf(tval))
		{
			return mUdf(float);
		}

		return z;
	}
#endif

	// save well intersections
	void QVelVolumeBuilder::computeWellIntersections(TypeSet<MultiID> selectedWells, int numWells, int ihorz, wellIntersectionStruct_ *wellIntersections,
		Array2DImpl<float>* topH, Array2DImpl<float>* botH,
		Array2DImpl<float>* topHextend, Array2DImpl<float>* botHextend)
	{
		int irow, icol;
		float tval;
		int jhorz = ihorz;
		//float lastT = 0.0;
		EM::ObjectID emid = EM::EMM().getObjectID( horizons_[horizonSortV_[ihorz]] );


		for (int iwell=0; iwell < numWells; iwell++)
		{
			bool fail = false;
			irow = icol = -1;

			Well::Data* wd = Well::MGR().get( selectedWells[iwell] );

			Well::Track timetrack = wd->track();
			//in 4.6 this wrecks it timetrack.toTime( *(wd->d2TModel()) ,timetrack);

			tval = mUdf(float);    

			WellHorIntersectFinder whpos( timetrack , wd->d2TModel());

			
			whpos.setHorizon( emid );
			// intersect well with horizon
			tval = whpos.findZIntersection();

			// work around bug in intersection code
			const Coord3& crd0 = timetrack.getPos( 0. );
			const BinID& bid0 = SI().transform( crd0 );
			//if (tracef != NULL) fprintf(tracef,"bid0 = %d %d \n",bid0.inl(),bid0.crl());
			if (tracef) fprintf(tracef,"\nIntersect well %d with horizon %d\n",iwell,ihorz);
#if 1
			// this can happen on a truncation above or below
			if (mUdf(float)==tval) {
				if (tracef) fprintf(tracef,"WellIntersector failed\n");
				tval = whpos.findZIntersection();
				int irow0 = bid0.inl() - cs_.hrg.start.inl();
				int icol0 = bid0.crl() - cs_.hrg.start.crl();
				// and adjust to possibly decimated subset cs_
				irow0/= cs_.hrg.step.inl();
				icol0/= cs_.hrg.step.crl();
				//if (tracef != NULL) fprintf(tracef,"cs vales %d %d %d %d \n",cs_.hrg.start.inl(),cs_.hrg.start.crl(),cs_.hrg.step.inl(),cs_.hrg.step.crl());

				irow = irow0;
				icol = icol0;

				tval = botHextend->get(irow0,icol0);

				if (tracef != NULL) fprintf(tracef,"Bottom horizon value is %g at (%d %d)\n",tval,irow0,icol0);
				if (botHextend->get(irow0, icol0) == topHextend->get(irow0,icol0)) {
					if (tracef != NULL) fprintf(tracef,"Bottom horizon truncates against top at %d %d = %g\n",irow0,icol0,topHextend->get(irow0,icol0));
				} else {
					// find deepest zero thickness
					if (tracef != NULL) fprintf(tracef,"Bottom horizon does NOT truncate against top at %d %d = %g\n",irow0,icol0,botHextend->get(irow0,icol0));
					float maxT = 0;
					for (int iih = 0; iih < numHorizons_; iih++)
					{
						if (iih != ihorz) {
							Array2DImpl<float>* testH = horizonTimes_[horizonSortV_[iih]];
							//if (tracef) fprintf(tracef,"test layer %d %g %g \n",iih,botHextend->get(irow0, icol0),testH->get(irow0, icol0));
							if (botHextend->get(irow0, icol0) == testH->get(irow0, icol0)) {
							if (tracef) fprintf(tracef,"found zero thick layer %d\n",iih);
							float dep = horizonTmax_[horizonSortV_[iih]];
							if (tracef) fprintf(tracef,"maxt %g\n",dep);
							if (dep >= maxT) {
								maxT = dep;
								jhorz = iih;
							}
						}
					}
					}
				}
				if (tracef != NULL) fprintf(tracef,"Well %d ihorz %d j horizon %d Intersect fail, fix with %g at %d %d\n",iwell,jhorz,horizonSortV_[ihorz],tval,irow0,icol0);
				if (tracef != NULL) fprintf(tracef,"TopH %g both %g \n",topH->get(irow0,icol0),botH->get(irow0,icol0));
				if (tval ==0) tval=mUdf(float);
				fail = true;
				
				
			}
			//else
#endif

#if 0
			if (tracef != NULL) 
				fprintf(tracef,"Well %d Intersect    ok, tval= %g \n",iwell, tval);
			else
				fprintf(tracef,"Well %d NO Intersect \n",iwell, tval);
			const float dah = wd->d2TModel()->getDah( 0, timetrack ) ;
			Coord3& crd2 = timetrack.getPos( dah );
			if (tracef != NULL) fprintf(tracef,"crd %g %g %g\n",crd2.x, crd2.y, crd2.z);
			if (mUdf(float)==tval) 
				if (tracef != NULL) fprintf(tracef,"  FAILED to find well/horizon \n");
				else
					if (fail) fprintf(tracef,"wellhorpos failed but I will use previous horizon!\n");

#endif

			if (mUdf(float) != tval) {
				float zz = wd->d2TModel()->getDepth( tval, timetrack );
				const Coord3&  crd = timetrack.getPos( zz );

				int irow0 = bid0.inl() - cs_.hrg.start.inl();
				int icol0 = bid0.crl() - cs_.hrg.start.crl();
				// and adjust to possibly decimated subset cs_
				irow0/= cs_.hrg.step.inl();
				icol0/= cs_.hrg.step.crl();


				const BinID& bid = SI().transform( crd );

				//if (tracef != NULL) fprintf(tracef," crd0 %d %d is %d %d \n",bid0.inl, bid0.crl, bid.inl, bid.crl);
				if (tracef != NULL) fprintf(tracef,"Horizon top %g bottom %g\n", topH->get(irow0,icol0),botH->get(irow0,icol0));


				// row/col in irl/crl
				// adjust row/col to 0-start

				irow = bid.inl() - cs_.hrg.start.inl();
				icol = bid.crl() - cs_.hrg.start.crl();
				// and adjust to possibly decimated subset cs_
				irow/= cs_.hrg.step.inl();
				icol/= cs_.hrg.step.crl();
				wellIntersections[iwell].x = irow;
				wellIntersections[iwell].y = icol;
				if (tracef != NULL) fprintf(tracef,"final intersection ij %d %d \n",irow,icol);
			}
			if (fail) {
				tval = botHextend->get(irow,icol);
				if (tracef) fprintf(tracef,"Final horizon value is %g\n",tval);
			}
			wellIntersections[iwell].z = tval;
			wellIntersections[iwell].ihorz = jhorz;
			if (tracef != NULL) fprintf(tracef,"horizon %d well %d intersec %d %d time %g horz\n",ihorz, iwell, irow, icol, tval);

			if (tracef != NULL) fprintf(tracef,"Finally: for horizon %d well %d intersection= %g %g %g \n\n",jhorz,iwell,
				wellIntersections[iwell].x, wellIntersections[iwell].y, wellIntersections[iwell].z);
			if (tracef != NULL) fflush(tracef);
		}
		
	}
	// interface for thin-plate spline code
	void QVelVolumeBuilder::addControlPoint(int wrow, int wcol,float intVel, int nrows,int ncols)
	{
		// odd requirements of the code
		// center X & Y, flip sense of Y and Z
		Vec v(wrow - nrows/2., intVel, wcol- ncols/2.);
		controlPoints_.push_back(v);
	}

	// add an intersection point to the interval velocity surface
	void QVelVolumeBuilder::addGridPoint(int wrow, int wcol, 
		Array2DImpl<double>* weights,
		Array2DImpl<float>* counts,
		Array2DImpl<double>* intVels,
		int nrows,
		int ncols,
		float intVel,
		Array2DImpl<double>* weightsTable)
	{
		// for thin-plate spline
		addControlPoint(wrow, wcol, intVel, nrows, ncols);

		// jbw 5.0 replace with table for speed

		// for 1/r type, populate entire area with 1/r weighted intersection value
		for (int irow = 0; irow < nrows; irow++)
		{
			for (int icol = 0; icol < ncols; icol++)
			{

#if 0	    //double Vel;
				double distance =sqrtf( ((float)(wrow - irow)*(float)(wrow - irow)) +
					((float)(wcol - icol)*(float)(wcol - icol)));
				double weight;

				if (distance == 0)
				{
					weight = 1.;
				}
				else
				{
					distance /= radius_;
					// 1/rsquared 
					weight = (1./pow((double)distance, (double)power_));
					//weight = (1./(distance * distance));
					//this flattens all the weighted values within 1 User Radius to 1.0
					if (weight > 1.) weight = 1.;
				}
#else
				int coldx = abs(wcol-icol);
				int rowdx = abs(wrow-irow);
				double weight = weightsTable->get(rowdx,coldx);
				//if (tracef) fprintf(tracef,"original %g table %g \n",weight,weight2);
#endif
				double ivel = intVels->get(irow, icol);
				// add in weighted value
				ivel += (intVel * weight);
				intVels->set(irow, icol, ivel);

				// accumulate weights
				double w = weights->get(irow, icol);
				weights->set(irow,icol, weight + w);

				// accumulate counts of adds
				float count = counts->get(irow,icol);
				counts->set(irow,icol, count+1);
			}
		}
	}

	// get marker or td intersecting surface row/col
	float QVelVolumeBuilder::markerDepth(Array2DImpl<float>*  horz, BufferString horizonName, MultiID well, int irow, int icol, bool useMarkers,HorizonMarkersList *hzMarkList)
	{
		int jrow, jcol;
		float tval;
		float depth, depth1;

		// if from td curves
		depth = depth1 = intersectWell(well, jrow, jcol, tval, horz);
		if (jrow < 0 || jcol < 0) depth = mUdf(float);
		if (jrow > cs_.nrInl() || jcol > cs_.nrCrl()) depth =  mUdf(float);

		if (!useMarkers)
		{
			return depth;
		}

		int hzmindex = -1;
		BufferString matchingmarker;
		for (int i=0; i < numHorizons_; i++)
			if (!strcmp(horizonName.str(), hzMarkList[i].horizons))
			{
				hzmindex = i;
				matchingmarker = hzMarkList[i].assocMarkers;
				//fprintf(tracef,"matching marker for horizon %s is %s\n",horizonName.str(), matchingmarker.str());
				break;
			}
			// directly from marker, no T/D involved
			Well::Data* wd = Well::MGR().get( well );
			for ( int idx=0; idx<wd->markers().size(); idx++ )
			{
				Well::Marker *w =  wd->markers()[idx];
				float d = w->dah();
				BufferString nm = w->name();

				// note: in ver 1, markername must == horizon name 
				// in version 2, users can override/add their own marker names

				if (nm.isEqual(matchingmarker))
				{
					const BinID& bid = SI().transform( wd->track().getPos(d) );

					// row/col in irl/crl
					// adjust row/col to 0-start
					irow = bid.inl() - cs_.hrg.start.inl();
					icol = bid.crl() - cs_.hrg.start.crl();
					// and adjust to possibly decimated subset cs_
					irow/= cs_.hrg.step.inl();
					icol/= cs_.hrg.step.crl();
					depth = d;
					return depth;
				}
			}
#if 1
			UsrMsg("Marker/Horizon intersection for horizon ");
			UsrMsg(horizonName);
			UsrMsg(" well ");
			UsrMsg(wd->name());
			UsrMsg("\n");
			UsrMsg("     not found, using T/D curve depth of ");
			MyUsrMsg(depth);
			UsrMsg("\n");
#endif
			// possibly fall thru with td curve depth loaded
			return depth;
	}

	// depth from migrated horizon
	float QVelVolumeBuilder::horizonDepth(Array2DImpl<float>* horz, int irow, int icol)
	{
		return horz->get(irow, icol);
	}

	// original time of horizon. If null, the time from the previous horizon will
	// have been migrated in by now
	float QVelVolumeBuilder::horizonTime(Array2DImpl<float>* horz, int irow, int icol)
	{
		return horz->get(irow, icol);
	}

	bool comparedx(const markerStruct_ &a,const markerStruct_ &b)
	{
		return a.weight > b.weight;
	}

	// get the interval velocity vector directly from gridding wells for a given i,j 
	double * QVelVolumeBuilder::getWellModelVelocities(int nsampls, double dt, int nrows, int ncols, int irow, int icol,
		Array2DImpl<double>* weights,Array2DImpl<float>* counts, Array2DImpl<double> *weightsTable)
	{
		double *ret = new double[nsampls];
		markerStruct_ * picks = new markerStruct_[nsampls];
		for (int i = 0; i < nsampls; i++) {
			picks[i].weight = 0.;
			picks[i].icol = 0;
			picks[i].irow = 0;
		}
		std::vector<markerStruct_> markervec;

		float tval=0.f;
		double *lastZ = new double[selectedWells_.nrItems()];
		// cache timetracks since we hit them over and over
		initWellTracks();

		// init t=0
		for (int iwell = 0; iwell < selectedWells_.nrItems();  iwell++)
		{
			int wrow = -1;
			int wcol = -1;
			float dummyt = 0.f;
			lastZ[iwell] = getWellZ(iwell, tval, wrow, wcol, dummyt);
		}
		// walk down samples
		for (int idt = 1; idt < nsampls; idt++)
		{
			tval = idt * dt;
			double TintVel = 0;
			double Tweight = 0;
			float Tcount = 0;

			// we need to intersect every time because of deviated wells
			// the distance can vary from time step to time step;
			for (int iwell = 0; iwell < selectedWells_.nrItems();  iwell++)
			{
				int wrow = -1;
				int wcol = -1;
				float dummyt = 0.f;
				double z = getWellZ(iwell, tval, wrow, wcol, dummyt);
				double intV = 2.*(z-lastZ[iwell])/dt;
				lastZ[iwell] = z;

#if 0 // use all anyway


				if (wrow < 0 || wcol < 0)  {
					picks[iwell].dx = mUdf(float);
					picks[iwell].intV = mUdf(float);
					markervec.push_back(picks[iwell]);
					continue;
				}
				if (wrow > cs_.nrInl() || wcol > cs_.nrCrl()) {
					picks[iwell].dx = mUdf(float);
					picks[iwell].intV = mUdf(float);
					markervec.push_back(picks[iwell]);
					continue;
				}
#endif

				// do a 1/r**n calculation for irow, icol for this well
				// but only if irow or icol changed in well, otherwise, well weights are the same as last ones.
				{
					int coldx, rowdx;
					if ((picks[iwell].icol != wcol ) || (picks[iwell].irow != wrow))

					{
#if 0
						// weight will have changed -- recalculate
						double distance =sqrtf( ((float)(wrow - irow)*(float)(wrow - irow)) +
							((float)(wcol - icol)*(float)(wcol - icol)));
						double weight;
						if (distance == 0)
						{
							weight = 1.;
						}
						else
						{
							distance /= radius_;
							// 1/rsquared 
							weight = (1./pow((double)distance, (double)power_));

							//this flattens all the weighted values within 1 User Radius to 1.0
							if (weight > 1.) weight = 1.;
						}
#else
						coldx = abs(wcol-icol);
						rowdx = abs(wrow-irow);
						double weight = weightsTable->get(rowdx,coldx);
#endif
						picks[iwell].weight = weight;
						picks[iwell].icol = wcol;
						picks[iwell].irow = wrow;
					}
					picks[iwell].intV = intV;
					markervec.push_back(picks[iwell]);
					//fprintf(tracef," well %d dx %g intV %g\n",iwell,distance,intV);
				}
			}
			// no need if not using nearest N
			//std::sort(markervec.begin(),markervec.end(), comparedx);
			// nearest 3 is something like a rough Nearest Neighbor
			//int iwmax = selectedWells_.nrItems() > 6 ? 6 : selectedWells_.nrItems();
			// all for now
			int iwmax = selectedWells_.nrItems();
			for (int iwell = 0; iwell < iwmax; iwell++)
			{
				float weight = markervec[iwell].weight;
				float intV = markervec[iwell].intV;

				// add in weighted value
				TintVel += (intV * weight);

				// accumulate weights

				Tweight += weight;

				// accumulate counts of adds

				Tcount++;

			}
			// done with adding points, 

			// get the int vel
			// normalize
			double vel = TintVel;
			double w = Tweight;

			vel /= w;
			ret[idt] = vel;
			markervec.clear();
		}

		return ret;
	}

	// get the interval velocity vector directly from gridding wells for a given i,j and t, within 2 horizons,
	// with conformant interpolation
	void QVelVolumeBuilder::getWellHorizonVelocities(int ihorz, double dt, int nrows, int ncols, int irow, int icol,
		float topTime, float botTime,wellIntersectionStruct_  **wellIntersections, float* velTrace, Array2DImpl<double> *weightsTable)
	{
		//if (tracef != NULL) fprintf(tracef,"getWellHorizonVelocities for horz %d row %d col %d\n",ihorz,irow,icol);
		float tval=0.f;
		double *lastZ = new double[selectedWells_.nrItems()];
		double * topZ = new double[selectedWells_.nrItems()];
		// cache timetracks since we hit them over and over
		initWellTracks();
		// real number of z steps we ned to fill
		int nz= (botTime-topTime)/dt;
		//if (tracef != NULL) fprintf(tracef,"toptime %g bot time %g\n",topTime,botTime);
		// zero thickness check
		if (nz <=0 || (botTime == mUdf(float)) || (topTime == mUdf(float)) ){
			//if (tracef != NULL) fprintf(tracef,"for BOT horz %d zero thickness for irow %d icol %d\n",ihorz,irow,icol);
			return;
		}
		//double *ret = new double[nz];
		int iTop = topTime/dt;
		int iBot = botTime/dt;
		for (int iwell = 0; iwell < selectedWells_.nrItems();  iwell++)
		{
			int wrow = -1;
			int wcol = -1;
			float wTop = topTime;
			if (ihorz > 0) wTop = wellIntersections[ihorz-1][iwell].z;
			float dummyt = 0.f;
			topZ[iwell] = lastZ[iwell] = getWellZ(iwell, wTop, wrow, wcol, dummyt);
			//if (icol == debugCol_)
			//    if (tracef != NULL) fprintf(tracef,"For Ihorz %d Well %d top z %g last z %g for well top %g\n",ihorz, iwell,topZ[iwell],lastZ[iwell],wTop);

		}
		//if (icol == debugCol_)
		//     if (tracef != NULL) fprintf(tracef, "\nhorizon toptinme %g bottime %g nz=%d\n",topTime, botTime, nz);


		for (int idt = iTop ; idt < iBot; idt++) // was <=
		{
			double TintVel = 0;
			double Tweight = 0;
			float Tcount = 0;

			//if (tracef != NULL) fprintf(tracef," for idt %d iTop %d iBot %d\n",idt,iTop,iBot);
			for (int iwell = 0; iwell < selectedWells_.nrItems();  iwell++)
			{
				int wrow = -1;
				int wcol = -1;
				float wTop = topTime;
				if (ihorz > 0) wTop = wellIntersections[ihorz-1][iwell].z;

				float wBot = wellIntersections[ihorz][iwell].z;

				float tstep = (wBot - wTop)/nz;
				if (tstep <=0) continue;
				tval = wTop + (1 + idt - iTop) * tstep;
				//if (icol == debugCol_)
				//    if (tracef != NULL) fprintf(tracef, "   for BOT horizon %d well %d well top %g bot %g tstep %g tval %g\n",ihorz, iwell, wTop, wBot, tstep, tval);
				if ((mUdf(float)==wTop) || (mUdf(float)==wBot))
				{
					//if (tracef != NULL) fprintf(tracef,"   well top or bot is missing, continue\n");
					continue;
				}
				float dummyt = 0.f;
				double z = getWellZ(iwell, tval, wrow, wcol, dummyt);
				double intV = 2.*(z-lastZ[iwell])/tstep;
				//if (tracef != NULL) fprintf(tracef,"          intv is %f\n",intV);
				lastZ[iwell] = z;

#if 0	// allow all influence   
				if (wrow < 0 || wcol < 0) {
					//if (tracef != NULL) fprintf(tracef,"well row col < 0\n");
					continue;
				}
				if (wrow > cs_.nrInl() || wcol > cs_.nrCrl()) {
					//if (tracef != NULL) fprintf(tracef,"well row col > max\n");
					continue;
				}

#endif

				// do a 1/r**n calculation for irow, icol for this well

				{

					int coldx = abs(wcol-icol);
					int rowdx = abs(wrow-irow);
					double weight = weightsTable->get(rowdx,coldx);


					//if (icol == debugCol_)
					//if (tracef != NULL) fprintf(tracef,"               vel=%g\n",intV);
					// add in weighted value
					TintVel += (intV * weight);

					// accumulate weights

					Tweight += weight;

					// accumulate counts of adds

					Tcount++;
				}
			}
			// done with adding points, 

			// get the int vel
			// normalize
			double vel = TintVel;
			double w = Tweight;
			if (Tcount <=0 || w <=0.) {
				//if (tracef != NULL) fprintf(tracef,"NO WEIGHT!\n");
				return;
			}
			vel /= w;
			velTrace[idt] = vel;
			/* if (icol == debugCol_) */
			//if (tracef != NULL) fprintf(tracef,"vel [%d] is %g at time %g\n",idt,vel, idt*dt);

		}

		velTrace[iTop] = velTrace[iTop+1];
		//if (icol == debugCol_)
		//if (tracef != NULL) fprintf(tracef,"Top vel [%d] is %g at time %g\n",iTop,velTrace[iTop], (iTop+1)*dt);

	}

	// here we actually compute the interval velocities for a single interval
	void QVelVolumeBuilder::computeIntVels(Array2DImpl<float>* topH, 
		Array2DImpl<float>* botH, 
		BufferString horizonName,
		int iBotH,
		Array2DImpl<float>* botDepths,
		Array2DImpl<float>* topDepths,
		Array2DImpl<double>* intVels,
		Array2DImpl<double>* weights,
		Array2DImpl<float>* counts,
		int nrows,
		int ncols,
		bool useMarkers,
		Array2DImpl<double>* weightsTable,
		wellIntersectionStruct_ * wellIntersections
		)
	{

		float topHrzT = 0.f;
		float topHrzZ = 0.f;
		float botHrzT = 0.f;
		float botPickZ = 0.f;
		//float lastIntv = 0.f;
		float intV     = 0.f;
		int hitrow=10;
		int hitcol=10;
		for (int iwell = 0; iwell < selectedWells_.nrItems();  iwell++)
		{
			int irow, icol;
			float tval=0.f;
			irow = icol = -1;

			if (tracef != NULL) fprintf(tracef,"\nCompute interval for horizon %s well %d",horizonName.str(), iwell);
			Well::Data* wd = Well::MGR().get( selectedWells_[iwell] );
			fprintf(tracef," %s\n",(wd->name()).str());

			//intersectWell(selectedWells_[iwell], irow, icol, tval, NULL);
			irow = wellIntersections[iwell].x;
			icol = wellIntersections[iwell].y;
			tval = wellIntersections[iwell].z;
			int jhorz = wellIntersections[iwell].ihorz;
			if (tracef != NULL) fprintf(tracef,"  intersection result %d %d time=%g hzn %d\n",irow,icol,tval,jhorz);
			hitrow = irow;
			hitcol = icol;
#if 1 // jbw new pinchout handling; care about intersection t
			if (mIsUdf(tval)) 
				continue;
#endif
			if (irow < 0 || icol < 0) 
				continue;
			if (irow > cs_.nrInl() || icol > cs_.nrCrl()) {
				if (tracef != NULL) fprintf(tracef,"  outside model %d %d\n",cs_.nrInl(),cs_.nrCrl());
				continue;
			}
			BufferString hithorizonName = horizonNames_->get(horizonSortV_[jhorz]);
			if (tracef != NULL) fprintf(tracef,"using horizon name %s for hit\n",hithorizonName.str());
			// bottom horizon time from horizon data at the well
			//botHrzT = horizonTime(botH, irow, icol);
			botHrzT = tval; // jbw new pinchout handling
			if (tracef != NULL) fprintf(tracef,"horizon time is %g\n",botHrzT);
			if (mIsUdf(botHrzT)) continue;

			// depth from pick (or T/D curve)
			botPickZ = markerDepth(botH, hithorizonName, selectedWells_[iwell], irow, icol, useMarkers,hzmarkersList_);
			if (tracef != NULL) fprintf(tracef,"ComputeIntVels: for well %d botpickz %g\n",iwell,botPickZ);
			if (mIsUdf(botPickZ)) continue;

			// top layer time
			topHrzT = horizonTime(topH,irow, icol);

			// migrated depth from prior layers & prior velocity calculations
			topHrzZ = horizonDepth(topDepths, irow, icol);
			fprintf(tracef,"At row %d col %d\n",irow,icol);
			fprintf(tracef,"Top horizon t %g bot horizon t %g\n",topHrzT, botHrzT);
			if (mIsUdf(topHrzT)) continue;
			fprintf(tracef,"Top hrz z %g bottom pick z %g\n",topHrzZ,botPickZ);

			float dz = botPickZ - topHrzZ;
			float dt = botHrzT - topHrzT;
			if (dt < 1.e-10) // zero thickness
			{
				//UsrMsg("Zero-thickness layer detected: ");
				//UsrMsg(horizonName);
				Well::Data* wd = Well::MGR().get( selectedWells_[iwell] );
				UsrMsg(" at well ");
				UsrMsg(wd->name());
				UsrMsg("\n");
				if (!zeroThicknessFound)
				{
					//logMessage("One or more zero-thickness layers have been detected.\n");
					//logMessage("See log file for more details\n");
					zeroThicknessFound = true;
				}
				continue;
			}
			if (dz < .1) // zero thickness
			{
				//UsrMsg("Zero-thickness time layer detected: ");
				//UsrMsg(horizonName);
				Well::Data* wd = Well::MGR().get( selectedWells_[iwell] );
				UsrMsg(" at well ");
				UsrMsg(wd->name());
				UsrMsg("\n");
				if (!zeroThicknessFound)
				{
					//logMessage("One or more zero-thickness layers have been detected.\n");
					//logMessage("See log file for more details\n");
					zeroThicknessFound = true;
				}
				continue;
			}

			intV = 2.*(dz/dt); // compute interval (interval average) velocity
			if (tracef != NULL) fprintf(tracef,"Interval velocity = %g\n",intV);
			// very generous sanity check (nonfatal)
			if (intV < 501. || intV > 25000)
			{
				logMessage("Discarded Unrealistic velocity calculated at ");
				logMessage(horizonName);
				Well::Data* wd = Well::MGR().get( selectedWells_[iwell] );
				logMessage(" at well ");
				logMessage(wd->name());
				logMessage(" Top horizon z ");
				logMessage(topHrzZ);
				logMessage(" Bottom  z ");
				logMessage(botPickZ);
				logMessage("\n   velocity = ");
				logMessage(intV);
				logMessage(" dz = ");
				logMessage(dz);
				logMessage(" dt = ");
				logMessage(dt);
				logMessage("\n");
				continue;
			}


			// add the calculated interval velocity to the current interval velocity grid
			addGridPoint(irow,icol, weights, counts, intVels, nrows, ncols, intV, weightsTable);
		}

		// done with adding points, now grid entore velocity grid

		// compute 1/r**n velocity grid, it "never fails"
		normalizeVelGrid(weights, counts, intVels, nrows, ncols);
		//fprintf(tracef,"after 1/r normalization vel at %d %d is %g\n",hitrow, hitcol, intVels->get(hitrow,hitcol));
		nrdone_++;

		// if asked for... and if possible, use spline fit
		if (gridder_ == QVel::QVelVolumeBuilder::ThinPlateSpline && controlPoints_.size() > 2)
		{
			TpsGrid gridder(nrows+1, ncols+1);
			gridder.setRegularization(relaxation_);
			gridder.setControlPoints(controlPoints_);

			if (gridder.calcTps() == 0)
			{
				for (int irow = 0; irow < nrows; irow++)
				{
					for (int icol = 0; icol < ncols; icol++)
					{    
						intVels->set(irow, icol, gridder.grid_[irow][icol]);

					}
				}
			}
		}
		nrdone_++;
		if (tracef != NULL) fprintf(tracef,"after surface fit vel at %d %d is %g\n",hitrow, hitcol, intVels->get(hitrow,hitcol));
		// from the velocity grid, get layer depths for next layer calc
		// if the velocity is missing, propagate top layer depth

		for (int irow = 0; irow < nrows; irow++)
		{
			for (int icol = 0; icol < ncols; icol++)
			{    
				double v = intVels->get(irow, icol);
				float dt = botH->get(irow, icol) - topH->get(irow, icol);
				float topd = topDepths->get(irow, icol);
				float dep;
				float bothrz = botH->get(irow, icol);
				if (mIsUdf(v) || mIsUdf(bothrz))
					dep = topd;
				else
					dep = (v*dt/2.)+topd;
				if ((irow==hitrow) && (icol == hitcol))
				{
					if (tracef != NULL) fprintf(tracef,"toph %g both %g\n",topH->get(irow, icol),botH->get(irow, icol));
					if (tracef != NULL) fprintf(tracef,"vel at %d %d is %g\n",irow,icol,v);
					if (tracef != NULL) fprintf(tracef,"topd dep %g new bot dep %g\n",topd,dep);
				}
				botDepths->set(irow,icol,dep);
			}
		}
		nrdone_++;
		if (tracef != NULL) fflush(tracef);
	}

	// get a vector of times for a row/col
	float * QVelVolumeBuilder::getHorizonTimes(int irow, int icol)
	{
		float * times = new float[numHorizons_];
		for (int i =0; i < numHorizons_; i++)
		{
			times[i] = (horizonTimes_[i])->get(irow, icol);
		}
		return times;
	}

	// get vector of velocities for a given row/col
	double * QVelVolumeBuilder::getHorizonVelocities(int irow, int icol)
	{
		double * vels = new double[numHorizons_];
		for (int i =0; i < numHorizons_; i++)
		{
			vels[i] = (intVels_[i])->get(irow, icol);
		}
		return vels;
	}

	// get the trace average interval velocities over an interval by
	// averaging all the finely sampled interval (instantaneous) velocites of
	// the incoming stacking vels or prior model run
	float QVelVolumeBuilder::getTraceAveIntVels(SeisTrc trcin, float traceIntVels[], float topTime, float botTime)
	{
		float velsum=0.;
		int velcount = 0;

		// our model starts at 0. Traces may not.
		int ind0 = trcin.info().sampling.nearestIndex<float>(topTime);
		if (ind0 < 0) ind0 = 0;
		int ind1 = trcin.info().sampling.nearestIndex<float>(botTime);
		if (ind1 >= trcin.size()) ind1 = trcin.size()-1;
		for (int i=ind0; i < ind1; i++)
		{
			velsum += traceIntVels[i];
			velcount ++;
		}
		return velsum / velcount;
	}

	// get a vector of adjustment multipliers such that
	// the input seismic can be calibrated with the calculated model
	bool QVelVolumeBuilder::computeIntervalAdjustTrace (int irow, int icol, 
		SeisTrc trcin, 
		int trcsize,
		float traceIntVels[],
		float intervalAdjust[])

	{
		int n = 0;
		if (numHorizons_ < 0) return true;

		try
		{
			float *horizonTimes = getHorizonTimes(irow, icol);
			double *intervalVels = getHorizonVelocities(irow, icol);
			int nBotHorizon = 0;
			float topTime = 0;
			float botTime = horizonTimes[horizonSortV_[0]];
			double intervalVelocity = intervalVels[0];
			float t = 0;
			double traceVelocity;
			float dt = trcin.info().sampling.step;
			for (n = 0; n < trcsize; n++, t += dt)
			{
				if (t <= topTime)
				{
					traceVelocity = intervalVelocity;
				}
				else if (t > botTime)
				{
					traceVelocity = getTraceAveIntVels(trcin, traceIntVels, topTime, t-dt); 
					if (mIsUdf(traceVelocity))
						intervalAdjust[nBotHorizon] = mUdf(float);
					else
						intervalAdjust[nBotHorizon] = intervalVelocity/traceVelocity;

					if (nBotHorizon < numHorizons_ -1)
					{
						topTime = botTime;
						nBotHorizon++;
						botTime = horizonTimes[horizonSortV_[nBotHorizon]];
						while (t > botTime && (nBotHorizon < numHorizons_ -1)) // missing/zero thickness
						{
							nBotHorizon++;
							botTime = horizonTimes[horizonSortV_[nBotHorizon]];
						}
						intervalVelocity = intervalVels[nBotHorizon];
					}
					else
						return true;
				}  
			}

			return true;
		}
		catch (...)
		{
			return false;
		}
	}

	// do the actual input trace calibration
	float QVelVolumeBuilder::adjustTraceVels(float inputTraceVel, 
		float inputModelVel, 
		double errorCorrection, 
		bool calibrateAll, 
		double weight,  
		float count
		)
	{
		if (mIsUdf(errorCorrection))
			return inputModelVel;

		// over entire horizon interval, just apply
		if (calibrateAll)
		{
			float ret = inputTraceVel * errorCorrection;
			return ret;
		}
		else // apply weighted calibration -- 100% at wells and dropping away from them
		{
			// normalize well 1/r sum of weight to 1
			double w = weight/count;
			float ret = (w * inputModelVel) + ((1.-w) * inputTraceVel);
			return ret;
		}
	}

	// jbw interval velocity seismic input, model input,
	//     merge the 2 using errorcorrection type & return in seismic trace
	bool QVelVolumeBuilder::processVelocityIntervalVolumeTrace(int irow, int icol, 
		int nsamples, 
		float inputTraceVels[], 
		float inputModelVels[],
		float dt, 
		bool calibrateAll, 
		float intervalAdjust[],
		SeisTrc trcin)
	{
		int n = 0;
		if (numHorizons_ <= 0) return true;

		try
		{

			float *horizonTimes = getHorizonTimes(irow, icol);
			int nBotHorizon = 0;
			float topTime = 0;
			float botTime = horizonTimes[horizonSortV_[0]];
			double errorCorrection = intervalAdjust[0];
			float t = trcin.info().sampling.start;
			//double traceVelocity;
			dt = trcin.info().sampling.step;
			for (n = 0; n < trcin.size(); n++, t += dt)
			{
				if (t <= topTime)
				{
					inputTraceVels[n] = inputModelVels[n];
				}
				else if (t > botTime)
				{
					if (nBotHorizon < numHorizons_ -1)
					{
						nBotHorizon++;
						botTime = horizonTimes[horizonSortV_[nBotHorizon]];
						while (t > botTime && (nBotHorizon < numHorizons_ -1)) // missing/zero thickness
						{
							nBotHorizon++;
							botTime = horizonTimes[horizonSortV_[nBotHorizon]];
						}
						errorCorrection = intervalAdjust[nBotHorizon];
						inputTraceVels[n] = adjustTraceVels(inputTraceVels[n], 
							inputModelVels[n], errorCorrection, calibrateAll, 
							weights_->get(irow,icol),counts_->get(irow,icol));
					}
					else // below bottom horizon -- just use trace vels
					{
						//int kk = n;
						/*
						inputTraceVels[n] = adjustTraceVels(inputTraceVels[n], 
						inputModelVels[n], errorCorrection, calibrateAll,
						weights_->get(irow,icol),counts_->get(irow,icol));
						*/
					}
				}
				else
					inputTraceVels[n] = adjustTraceVels(inputTraceVels[n], 
					inputModelVels[n], 
					errorCorrection, 
					calibrateAll,
					weights_->get(irow,icol),
					counts_->get(irow,icol));
			}
			delete [] horizonTimes;
			return true;
		}
		catch (...)
		{
			return false;
		}
	}


	// create a trace from model directly
	bool QVelVolumeBuilder::processVelocityIntervalTrace( int irow, int icol, 
		int nsampls,
		float modelTrace[],
		float workTrace[],
		float dt,
		SeisTrc trcin,
		bool isTDC,
		Array2DImpl<double> *weightsTable)
	{
		//if (tracef != NULL) fprintf(tracef," processvleocityintervaltrace row %d col %d \n",irow, icol);
		int n = 0;


		// if desired, get the gridded velocity from T/D curves at this xy;

		if (numHorizons_ <= 0) // just do time-depth-curves 1/r gridding
			try
		{
			int nrows = nrows_;
			int ncols = ncols_;
			//if (tracef != NULL) fprintf(tracef," processvleocityintervaltrace nrow %d ncol %d \n",nrows, ncols);
			double *intervalVels = getWellModelVelocities(nsampls, dt, nrows, ncols, irow, icol,
				weights_, counts_, weightsTable);

			for (n = 0; n < nsampls; n++)
			{

				modelTrace[n] = intervalVels[n];
				// test modelTrace[n] = irow;

			}
			delete [] intervalVels;
		}
		catch (...)
		{
			return false;
		}

		//if (tracef != NULL) fprintf(tracef," processvleocityintervaltrace row %d col %d \n",irow, icol);
		if (numHorizons_ <= 0) return true;

		try
		{

			float *horizonTimes = getHorizonTimes(irow, icol);
			double *intervalVels = getHorizonVelocities(irow, icol);
			int nBotHorizon = 0;
			float topTime = timeDatum_;
			float botTime = horizonTimes[horizonSortV_[0]];
			double intervalVelocity = intervalVels[0];
			float t = 0; // always start at t=0
			//double traceVelocity;
			float modelBot = nsampls * dt;


			float topH = 0;
			if (isTDC)
			{
				for (int ihorz=0; ihorz<= numHorizons_; ihorz++) 
				{
					float botH = topH;

					// propagate up
					//if (tracef != NULL) fprintf(tracef,"looking at hz %d out of %d\n",ihorz,numHorizons_);

					{
						if (ihorz < numHorizons_) {
							botH = horizonTimes[horizonSortV_[ihorz]];
							//if (tracef) fprintf(tracef,"using horizonTimes %d\n",ihorz);
						}
						else {
							botH = modelBot;
							//if (tracef) fprintf(tracef,"bot hzr is model bot\n");
						}	
					}
					// conformant interpolation
					getWellHorizonVelocities (ihorz, dt, nrows_,  ncols_,  irow, icol,
						topH, botH, wellIntersections_, modelTrace, weightsTable);
					topH = botH;
				}
			} 
			else
			{
				for (n = 0; n < nsampls; n++, t += dt)
				{
					if (t <= topTime)
					{
						modelTrace[n] = intervalVelocity;
					}
					else if (t > botTime)
					{
						if (nBotHorizon < numHorizons_ -1)
						{
							nBotHorizon++;
							botTime = horizonTimes[horizonSortV_[nBotHorizon]];
							while (t > botTime && (nBotHorizon < numHorizons_ -1)) // missing/zero thickness
							{
								nBotHorizon++;
								botTime = horizonTimes[horizonSortV_[nBotHorizon]];
							}
							intervalVelocity = intervalVels[nBotHorizon];
							modelTrace[n] = intervalVelocity;
						}
						else
							modelTrace[n] = intervalVelocity;
					}
					else
						modelTrace[n] = intervalVelocity;

				}
				delete [] horizonTimes;
				delete [] intervalVels;
				// drop thru
			}
		}
		catch (...)
		{
			return false;
		}


		return true;
	}

	void QVelVolumeBuilder::avgToIntervalVelocity(float velocities[], int nvels, float dt)
	{
		double t = -dt;
		double tp1 = 0.0;
		double z = -dt*velocities[0];
		double zp1 = 0.0;
		for (int n = 0; n < nvels-1; n++)
		{
			t = tp1;
			tp1 += dt;
			double zm1 = z;
			z = zp1;
			zp1 = tp1*velocities[n+1];
			velocities[n] = (float)((zp1 - zm1));
		}
		velocities[nvels-1] = velocities[nvels-2];
	}

	// build a vector of horizon depths from the horizon times & the model
	float* QVelVolumeBuilder::getHorizonDepths(float horizonTimes[], double intervalVels[], int nHorizons)
	{
		float topDepth = 0;
		float botDepth = 0;
		float topTime  = 0;
		float *depths = new float[nHorizons];
		for (int i = 0; i < nHorizons; i++)
		{
			if (!mIsUdf(intervalVels[i]))
				botDepth = (horizonTimes[i] - topTime) * intervalVels[i] + topDepth;
			depths[i] = botDepth;
			topDepth = botDepth;
			topTime =  horizonTimes[i];
		}
		return depths;
	}

	// if bottom horizon value missing, replace with top 
	template <class T>
	void QVelVolumeBuilder::handleMissing(Array2DImpl<T> * topH, 
		Array2DImpl<T> * botH, 
		int nrows,
		int ncols)
	{
		for ( int irow=0; irow < nrows; irow++ )
		{
			for ( int icol=0; icol < ncols; icol++ )
			{
				float z = botH->get(irow,icol);
				if (mIsUdf(z))
					botH->set(irow, icol, topH->get(irow,icol));
			}
		}
	}

	// calculate horizon min/max times for sorting horizons
	// max is used as a tiebreaker if min values are equal (zero thickness pinchout)
	// crossing horizons will cause bad results
	void QVelVolumeBuilder::scanHorizon(Array2DImpl<float> *horizonTimes, int nrows, int ncols, float & minV, float & maxV)
	{
		minV = mUdf(float);
		maxV = mUdf(float);
		for (int irow = 0; irow < nrows; irow++)
			for (int icol = 0; icol < ncols; icol++)
			{
				if (mIsUdf(minV))
					minV = horizonTimes->get(irow, icol);
				else
				{
					float t = horizonTimes->get(irow, icol);
					if (!mIsUdf(t))
						minV = minV < t ? minV : t;
				}

				if (mIsUdf(maxV))
					maxV = horizonTimes->get(irow, icol);
				else
				{
					float t = horizonTimes->get(irow, icol);
					if (!mIsUdf(t))
						maxV = maxV > t ? maxV : t;
				}
			}
	}

	// expand the horizon data to the entire survey range, filling with mUdf is needed
	// also, pack tightly by compressing out steps > 1
	void QVelVolumeBuilder::expandHorizon(Array2DImpl<float>* outH, 
		Array2D<float> *inH, 
		StepInterval<int> rowrg, 
		StepInterval<int> colrg)
	{
		//subset survey in "real" inl/crl
		int subselrowstart = cs_.hrg.start.inl();
		int subselcolstart = cs_.hrg.start.crl();
		int subselrowend = subselrowstart + cs_.hrg.nrInl()*cs_.hrg.step.inl() -1;
		int subselcolend = subselcolstart + cs_.hrg.nrCrl()*cs_.hrg.step.crl() -1;
		// horizon limits in "real" inl/crl
		int irowbeg = std::max(subselrowstart, rowrg.start);
		int icolbeg = std::max(subselcolstart, colrg.start);
		int irowend = std::min(subselrowend, rowrg.stop);
		int icolend = std::min(subselcolend, colrg.stop);
		int irb = (irowbeg - subselrowstart)/cs_.hrg.step.inl();
		int icb = (icolbeg - subselcolstart)/cs_.hrg.step.crl();
		int irend = 0;
		int icend = 0;
		for ( int row=irowbeg, ir=irb; row<=irowend; row+=cs_.hrg.step.inl(), ir++ )
		{
			int irow = rowrg.getIndex(row);
			for ( int col=icolbeg, ic=icb; col<=icolend; col+=cs_.hrg.step.crl(), ic++ )
			{
				int icol = colrg.getIndex(col);
				float t = inH->get(irow, icol);
				//if (tracef != NULL) fprintf(tracef,"hz get %d %d to %d %d -> %g\n",irow,icol,ir,ic, t);
				
				//if (ir == 130 && ic == 43)
					//if (tracef != NULL) fprintf(tracef,"hz get %d %d to %d %d -> %g\n",irow,icol,ir,ic, t);
				//if (!mIsUdf(t))
				  outH->set(ir, ic, t);
				icend = ic;
			}
			irend = ir;
		}
	}

	void QVelVolumeBuilder::filloutZTH(int fromhz, int tohz, float * horizonTmin, float * horizonTmax)
	{
		if (tracef) fprintf(tracef,"Replace missings in horizon %d with values from %d\n",fromhz,tohz);
		for (int irow=0; irow < nrows_; irow++)
		for (int icol=0; icol < ncols_; icol++)
		{
			float replval = horizonTimes_[fromhz]->get(irow, icol);
			float targetval = horizonTimes_[tohz]->get(irow, icol);
					if (!mIsUdf(replval) && mIsUdf(targetval)) {
						horizonTimes_[tohz]->set(irow, icol, replval);
					if (replval < horizonTmin[tohz])
							horizonTmin[tohz] = replval;
						if (replval > horizonTmax[tohz])
							horizonTmax[tohz] = replval;
					}
		}
	}
	// see if an adjacent cell is a zero thickness intersection
	bool QVelVolumeBuilder::fixMissing(int jhorz, int jrow, int jcol, float * horizonTmin, float * horizonTmax, int numHoriz)
	{
		bool ret = false;
		float minTimedx = .006;
		for (int ihorz=0; ihorz < numHoriz; ihorz++) {
			if (ihorz != jhorz) {
				// row + 1
				if (jrow < (nrows_ -1)) {
					float testval = horizonTimes_[ihorz]->get(jrow + 1, jcol);
					float val = horizonTimes_[jhorz]->get(jrow + 1, jcol);
					// next cell not missing, and is zero thickness
					//if (!mIsUdf(val) && !mIsUdf(testval) && tracef) fprintf(tracef,"test %g\n",	abs(val - testval));
					if (!mIsUdf(val) && !mIsUdf(testval) && (abs(val - testval) <= minTimedx)) {
						//if (tracef) fprintf(tracef,"zero thick at fix %d %d \n",jrow+1,jcol);
						//if (tracef) fprintf(tracef,"1values %g %g for h %d and %d\n",val,testval,jhorz,ihorz);
						filloutZTH(ihorz, jhorz, horizonTmin, horizonTmax);
						float replval =  horizonTimes_[ihorz]->get(jrow + 1, jcol);
						horizonTimes_[jhorz]->set(jrow,jcol, replval);
						if (replval < horizonTmin[jhorz])
							horizonTmin[jhorz] = replval;
						if (replval > horizonTmax[jhorz])
							horizonTmax[jhorz] = replval;
						
						ret = true;
					}
				}
				// col + 1
				if (jcol < (ncols_ -1)) {
					float testval = horizonTimes_[ihorz]->get(jrow , jcol + 1);
					float val = horizonTimes_[jhorz]->get(jrow , jcol + 1);
					// next cell not missing, and is zero thickness
					//if (!mIsUdf(val) && !mIsUdf(testval) && tracef) fprintf(tracef,"test %g\n",	abs(val - testval));
					if (!mIsUdf(val) && !mIsUdf(testval) && (abs(val - testval) <= minTimedx)) {
						//if (tracef) fprintf(tracef,"zero thick at fix %d %d \n",jrow,jcol+1);
						//if (tracef) fprintf(tracef,"2values %g %g for h %d and %d\n",val,testval,jhorz,ihorz);
						filloutZTH(ihorz, jhorz, horizonTmin, horizonTmax);
						float replval =  horizonTimes_[ihorz]->get(jrow, jcol + 1);
						horizonTimes_[jhorz]->set(jrow,jcol, replval);
						if (replval < horizonTmin[jhorz])
							horizonTmin[jhorz] = replval;
						if (replval > horizonTmax[jhorz])
							horizonTmax[jhorz] = replval;
						ret = true;
					}
				}
				// row -1 
				if (jrow > 0) {
					float testval = horizonTimes_[ihorz]->get(jrow - 1, jcol);
					float val = horizonTimes_[jhorz]->get(jrow - 1, jcol);
					// next cell not missing, and is zero thickness
					//if (!mIsUdf(val) && !mIsUdf(testval) && tracef) fprintf(tracef,"test %g\n",	abs(val - testval));
					if (!mIsUdf(val) && !mIsUdf(testval) && (abs(val - testval) <= minTimedx)) {
						//if (tracef) fprintf(tracef,"zero thick at fix %d %d \n",jrow-1,jcol);
						//if (tracef) fprintf(tracef,"3values %g %g for h %d and %d\n",val,testval,jhorz,ihorz);
						float replval =  horizonTimes_[ihorz]->get(jrow - 1, jcol);
						filloutZTH(ihorz, jhorz, horizonTmin, horizonTmax);
						horizonTimes_[jhorz]->set(jrow,jcol, replval);
						if (replval < horizonTmin[jhorz])
							horizonTmin[jhorz] = replval;
						if (replval > horizonTmax[jhorz])
							horizonTmax[jhorz] = replval;
						ret = true;
					}
				}
				// col - 1
				if (jcol > 0) {
					float testval = horizonTimes_[ihorz]->get(jrow , jcol - 1);
					float val = horizonTimes_[jhorz]->get(jrow , jcol -1);
					// next cell not missing, and is zero thickness
					//if (!mIsUdf(val) && !mIsUdf(testval) && tracef) fprintf(tracef,"test %g\n",	abs(val - testval));
					if (!mIsUdf(val) && !mIsUdf(testval) && (abs(val - testval) <= minTimedx)) {
						//if (tracef) fprintf(tracef,"zero thick at fix %d %d \n",jrow,jcol-1);
						//if (tracef) fprintf(tracef,"4values %g %g for h %d and %d\n",val,testval,jhorz,ihorz);
						float replval =  horizonTimes_[ihorz]->get(jrow, jcol -1);
						filloutZTH(ihorz, jhorz, horizonTmin, horizonTmax);
						horizonTimes_[jhorz]->set(jrow,jcol, replval);
						if (replval < horizonTmin[jhorz])
							horizonTmin[jhorz] = replval;
						if (replval > horizonTmax[jhorz])
							horizonTmax[jhorz] = replval;
						ret = true;
					}

				}
				
			}

		}
		return ret;
	}

	void QVelVolumeBuilder::extendZTHorizons(float *horizonTmin, 
		float *horizonTmax, 
		int   *horizonSortV, 
		int    numHoriz)
	{
		if (tracef) fprintf(tracef,"**analyze pinch-outs and zero-thickness intersections**\n");
		// inspect 
		bool changed = false;
		do {
		changed = false;
		if (tracef) fprintf(tracef, "\n*** loop\n");
		for (int ihorz=0; ihorz < numHoriz; ihorz++) {
			float tmin = horizonTmin[ihorz];
			//if (tracef) fprintf(tracef,"   *** horizon %d \n",ihorz);
			do {
			changed = false;
			int udfCount = 0;
			for (int irow = 0; irow < nrows_; irow++)
				for (int icol = 0; icol < ncols_; icol++)
				{
					float val = horizonTimes_[ihorz]->get(irow, icol);
					if (mIsUdf(val)) {
					   udfCount++;
					   if (fixMissing(ihorz, irow, icol, horizonTmin, horizonTmax, numHoriz)) {
						   changed = true;
						   //if (tracef != NULL) fprintf(tracef,"zt fix for horz %d \n",ihorz);
						   //if (tracef != NULL) fflush(tracef);
					   }

					}
				}
				if (tracef) fprintf(tracef," loop pass complete for %d number of missing %d\n",ihorz,udfCount);
			} while (changed);
			//if (tracef) fprintf(tracef,"no change 1\n");
		}
		} while (changed);
	}
	// order the horizons top to bottom
	// NB: will fail if horizons cross, or are non-truncating patches
	void QVelVolumeBuilder::orderHorizons(float *horizonTmin, 
		float *horizonTmax, 
		int   *horizonSortV, 
		int    numHoriz)
	{
		std::list<int> sortlist;
		std::list<int>::iterator it;
		sortlist.push_back(0);
		int nrItems = 1;
		float tEps = 0.003;
		for (int i = 1; i < numHoriz; i++)
		{
			bool inserted = false;
			//fprintf(tracef,"begin nrItems = %d\n",nrItems);
			it = sortlist.begin();
			for (int j = 0; j < nrItems; j++, ++it)
			{
				
				//fprintf(tracef,"begin loop j = %d\n",j);
				
				float dt = horizonTmin[*it] - horizonTmin[i];
				//fprintf(tracef,"hz %d vs %d dt is %g\n",i,*it,dt);
				if (dt > tEps)
				{
					//fprintf(tracef,"hz %d less\n",i);
					sortlist.insert(it,i);
					inserted = true;
					nrItems++;
					break;
				}
				else if ((-dt) > tEps) {
					//fprintf(tracef," try again  j = %d\n",j);
					continue;
				}
				// tie? check max time
				else if (fabs(dt) < tEps)
				{    
					//fprintf(tracef,"% tie %d vs %d\n",i,*it);
					inserted = true;
					//fprintf(tracef,"compare max %g %g \n",horizonTmax[i] , horizonTmax[*it]);
					if (horizonTmax[i] < horizonTmax[*it])
					{
						//fprintf(tracef,"insert %d\n",*it);
						sortlist.insert(it,i);
						nrItems++;
					}
					else
					{
						++it;
						//fprintf(tracef,"else insert %d\n",*it);
						sortlist.insert(it,i);
						nrItems++;
					}
					break;
				}
			
			}
			if (!inserted) { 
				//fprintf(tracef,"add on end %d\n",i);
				sortlist.push_back(i);
				nrItems++;
				if (horizonTmax[i] < horizonTmax[*it]) { 
					//fprintf(tracef,"Warn:Crossing horizons or patches, hzn %d is deeper than %d, and shallower!\n",*it,i);
					if (tracef) fprintf(tracef,"Hzn tmax %d %g less than %d %g \n",i,horizonTmax[i] ,*it,horizonTmax[*it] );
				//uiMSG().warning("Crossing horizons or disjoint patches, cannot handle properly\n");
				}
			}
		}

		int i = 0;
		for (it=sortlist.begin(); it!=sortlist.end(); ++it)
		{
			//int ival = *it;
			BufferString hn = horizonNames_->get(*it);
			fprintf(tracef,"sorted horizon list %d index %d %s min %g max %g\n",i,*it,hn.str(),horizonTmin[*it],horizonTmax[*it]);
			horizonSortV[i++] = *it;
		}
	}

	// main process loop for QVel volume building algorithm
	int    QVelVolumeBuilder::nextStep()
	{
		if (getState() == Stop) return Executor::Finished();
		if (firstTime_)
		{


			horizonSortV_ = new int[numHorizons_+1];
			horizonTmin_  = new float[numHorizons_+1];
			horizonTmax_  = new float[numHorizons_+1];
			horizonNames_ = new BufferStringSet();
			numWells_= selectedWells_.nrItems();
			wellIntersections_ = new wellIntersectionStruct_ *[numHorizons_+1];
			for (int i = 0; i < numHorizons_+1; i++) {
				wellIntersections_[i] = new wellIntersectionStruct_ [numWells_];
				for (int j =0; j < numWells_; j++)
					wellIntersections_[i][j].z=mUdf(float);
			}

			nrows_ = cs_.nrInl();
			ncols_ = cs_.nrCrl();
#if 0 // deadlocks in ODT 5.0 ?
			UsrMsg("Number of rows ");
			MyUsrMsg(nrows_);
			UsrMsg(" Number of cols ");
			MyUsrMsg(ncols_);
			UsrMsg("\n");
#endif	
			nrdone_ = 0;    
			totnr_ = numHorizons_ * (6); // for progress meter

			// calculate a table of all possible 1/r weights
			weightsTable_ = new Array2DImpl<double>( nrows_, ncols_ );
			for (int irow = 0; irow < nrows_; irow++)
			{
				for (int icol = 0; icol < ncols_; icol++)
				{
					//double Vel;
					double distance =sqrtf( ((float)(irow)*(float)(irow)) +
						((float)(icol)*(float)(icol)));
					double weight;

					if (distance == 0)
					{
						weight = 1.;
					}
					else
					{
						distance /= radius_;
						// 1/rsquared 
						weight = (1./pow((double)distance, (double)power_));
						//weight = (1./(distance * distance));
						//this flattens all the weighted values within 1 User Radius to 1.0
						if (weight > 1.) weight = 1.;
					}
					weightsTable_->set(irow,icol,weight);
				}
			}
			// we need all of the horizon times and interval velocities to build trace output
			for (int i = 0; i < numHorizons_; i++)
			{    
				intVels_[i] = new Array2DImpl<double>( nrows_, ncols_ );
				intVels_[i]->setAll( 0 );
				mGetHor(tempH_,horizons_[i]);
				Array2D<float>* tempA = tempH_->createArray2D(tempH_->sectionID(0),0);
				horizonTimes_[i] =  new Array2DImpl<float>( nrows_, ncols_ );
				rawHorizonTimes_[i] =  new Array2DImpl<float>( nrows_, ncols_ );
				BufferString nm2 = tempH_->name();
				horizonNames_->add( tempH_->name() );
				horizonTimes_[i]->setAll(mUdf(float));
				rawHorizonTimes_[i]->setAll(mUdf(float));
				StepInterval<int> rowrg = tempH_->geometry().rowRange();
				StepInterval<int> colrg = tempH_->geometry().colRange();
				expandHorizon(horizonTimes_[i], tempA, rowrg, colrg);
				expandHorizon(rawHorizonTimes_[i], tempA, rowrg, colrg);
				scanHorizon(horizonTimes_[i], nrows_, ncols_, horizonTmin_[i], horizonTmax_[i]);
				fprintf(tracef,"horizon %d %s min %g max %g\n",i,(horizonNames_->get(i)).str(),horizonTmin_[i],horizonTmax_[i]);
				delete tempA;
			}
			nrdone_++;

			if (numHorizons_ > 0) {
			// extension Zero Thickness creates zero-thickness zones
			extendZTHorizons(horizonTmin_, horizonTmax_, horizonSortV_, numHorizons_);
			// sort horizons mintime to maxtime
			orderHorizons(horizonTmin_, horizonTmax_, horizonSortV_, numHorizons_);
			}

			setName("Interval Layer");

			weights_ = new Array2DImpl<double>( nrows_, ncols_ );
			counts_ = new Array2DImpl<float>( nrows_, ncols_ );

			topH_ = NULL;

			botHextend_ = botH0_ = botH_= new Array2DImpl<float>( nrows_, ncols_ );
			botH_->setAll( 0 );

			botDepths_= new Array2DImpl<float>( nrows_, ncols_ );
			botDepths_->setAll( 0 );

			topDepths_= new Array2DImpl<float>( nrows_, ncols_ );
			topDepths_->setAll( 0 );

			firstTime_ = 0;
			for (int i = 0; i < numHorizons_; i++)
			{
				BufferString nm = horizonNames_->get(i);
				nm = horizonNames_->get(horizonSortV_[i]);
			}
		}
		// loop over horizons, working down
		if (ihorz_ < numHorizons_)
		{
			controlPoints_.erase(controlPoints_.begin(), controlPoints_.end());
			BufferString horizonName = horizonNames_->get(horizonSortV_[ihorz_]);
			weights_->setAll( 0 );
			counts_->setAll( 0 );
			// working down to last bottom is current top horizon
			topH_ = botH_;
			topHextend_ = botHextend_;
			botH_ =  rawHorizonTimes_[horizonSortV_[ihorz_]];    
			botHextend_ = horizonTimes_[horizonSortV_[ihorz_]];
			// deal with missing values in horizons
			handleMissing(topHextend_, botHextend_, nrows_, ncols_);
			nrdone_++;
			// save well intersections
			computeWellIntersections(selectedWells_, numWells_,ihorz_,wellIntersections_[ihorz_], topH_, botH_, topHextend_, botHextend_);
			// propagate using extended top horizon
			
			// here we go
			computeIntVels(topHextend_, 
				botHextend_, 
				horizonName, 
				ihorz_,  
				botDepths_, 
				topDepths_, 
				intVels_[ihorz_], 
				weights_,
				counts_,
				nrows_,
				ncols_,
				useMarkers_,
				weightsTable_,
				wellIntersections_[ihorz_]);


			delete topDepths_;
			// we don't store all the horizon depths, just the 2 we are working on
			// horizon depths can always be recalculated from the model
			topDepths_ = botDepths_;
			botDepths_ = new Array2DImpl<float>( nrows_, ncols_ );
			botDepths_->setAll( 0 );

			// propagate good velocities into any missing vels
			if (ihorz_ > 0 )
				handleMissing(intVels_[ihorz_ - 1],intVels_[ihorz_],nrows_, ncols_);
			ihorz_++;
			nrdone_++;
			return Executor::MoreToDo();
		}
		delete topDepths_;
		delete botDepths_;
		delete [] horizonTmin_;
		delete [] horizonTmax_;
		UsrMsg("Job completed.\n");
		return Executor::Finished();
	}
	void QVelVolumeBuilder::propagateWellIntersections(float modelBot)
	{
		if (numHorizons_ <= 0) return;
		// fake bottom layer
		for (int iwell = 0; iwell < numWells_; iwell++) {
			wellIntersections_[numHorizons_][iwell].z = modelBot;
			wellIntersections_[numHorizons_][iwell].x = wellIntersections_[numHorizons_ -1 ][iwell].x;
			wellIntersections_[numHorizons_][iwell].y = wellIntersections_[numHorizons_ -1 ][iwell].y;
			wellIntersections_[numHorizons_][iwell].ihorz = wellIntersections_[numHorizons_ -1 ][iwell].ihorz;
		}
		for (int iwell=0; iwell < numWells_; iwell++)
		{
			for (int ihorz=0; ihorz < numHorizons_; ihorz++)
			{
				if ( mUdf(float) == wellIntersections_[ihorz][iwell].z) {
					//if (tracef != NULL) fprintf(tracef,"for horiz %d well %d intersection is missing",ihorz, iwell);
					int thorz = ihorz+1;
					while (mUdf(float)== wellIntersections_[thorz][iwell].z)
						thorz++;
					wellIntersections_[ihorz][iwell].z = wellIntersections_[thorz][iwell].z;
					//if (tracef != NULL) fprintf(tracef," propagate up for z value= %g\n",wellIntersections_[thorz][iwell].z);
				}

			}
		}
	}
	void QVelVolumeBuilder::computeTDCVels(Array2DImpl<float> *topH,Array2DImpl<float> *botH, int ihorz, int nrows, int ncols)
	{

	}
	extern "C" void logToFile(BufferString mes) {
		std::cerr << mes <<std::endl;
	}

}; //namespace
