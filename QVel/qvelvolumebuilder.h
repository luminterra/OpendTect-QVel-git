/*
 Copyright (C) 2012 LuminTerra, LLC All rights reserved.

 This file is part of OpendTect and may be used under the terms of:

 The GNU General Public License version 3 or higher, as published by
 the Free Software Foundation.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

 Ver 1.1 JB West 4/2012

*/

#ifndef qvelvolumebuilder_h
#define qvelvolumebuilder_h

#include "executor.h"
#include "cubesampling.h"
#include "task.h"
#include "thread.h"
#include "veldesc.h"
//#include "bufstring.h"
#include "emhorizon3d.h"
#include "emmanager.h"
#include "seisread.h"
#include "seiswrite.h"
#include "ioobj.h"
#include "transl.h"
#include "arrayndimpl.h"
#include "linalg3d.h"
#include "commondefs.h"
#include "welldata.h"
#include "wellextractdata.h"
#include "wellio.h"
#include "welllog.h"
#include "welllogset.h"
#include "wellman.h"
#include "wellreader.h"
#include "wellmarker.h"
#include "welldata.h"
#include "welld2tmodel.h"
#include "welltrack.h"
#include "welllog.h"
#include "welllogset.h"
#include "wellman.h"
#include "wellmarker.h"
#include "uimsg.h"
#include "msgh.h"
extern "C" void logToFile(BufferString mes);

class IOObj;
class SeisTrc;
class SeisTrcReader;
class SeisTrcWriter;
class SeisSequentialWriter;
class Vec;



namespace QVel
{


    typedef struct {
	const char* horizons;
	const char* assocMarkers;
    }
    HorizonMarkersList;

/*!Reads in a volume Vint, and writes out a volume
    Vint. */

class dll_export QVelVolumeBuilder: public Executor
{
public:
	
	QVelVolumeBuilder(SeisTrcReader* rdr_, 
		SeisTrcWriter* wrr_, 
		TypeSet<MultiID> selectedWells_, 
		BufferStringSet markerNames_,
		TypeSet<MultiID> horizons_, 
		int nrHorizons,
		float radius,
		float power,
		bool calibrateAll,
		bool useMarkers,
		CubeSampling cs,
		int griddertype,
		float relaxation,
		bool isTDC,
		HorizonMarkersList *hzmarkersList);


	~QVelVolumeBuilder();

	int			nextStep();
	
	
    const char*		message() const	{ return "Calculating QVel velocity model..."; }  
	int   nrTraces() { return nrTraces_; }
	Array2DImpl<double>** intVels_;
	Array2DImpl<float>** horizonTimes_;
	BufferStringSet  *horizonNames_;
	Array2DImpl<double>*  weights_;
	int nrows_;
	int ncols_;
	int numHorizons() { return numHorizons_; }
	bool computeIntervalAdjustTrace (int row, int col, 
		SeisTrc trcin, 
		int trcsize,
		float traceIntVels[],
		float intervalAdjust[]);
		
	bool processVelocityIntervalVolumeTrace(int row, int col, 
		int nTrace, 
		float inputTrace[],
		float modelTrace[],
		float dt, 
		bool calib, 
		float intervalAdjust[],
		SeisTrc trcin);

	bool processVelocityIntervalTrace(  int irow, int icol, 
				int nsampls,
				float inputTrace[],
				float workTrace[],
				float dt,
				SeisTrc trcin_,
				bool isTDC,
				Array2DImpl<double> *weightsTable);

	void avgToIntervalVelocity(float velocities[], int nvels, float dt);
	CubeSampling cs_;
	int nrdone_;
	int totnr_;
	BufferString getMessages() { return messages_; }
	enum Gridder { InverseDx, ThinPlateSpline };
	double * getWellModelVelocities(int nsampls, double dt, int nrows, int ncols, int irow, int icol,
	    Array2DImpl<double>* weights,Array2DImpl<float>* counts, Array2DImpl<double>* weightsTable);
	double getWellZ(MultiID well, float & tval, int &irow, int &icol, float dt);
	double getWellZ(int iwell, float & tval, int &irow, int &icol, float dt);
	void propagateWellIntersections(float modelbot);
	Array2DImpl<double>* getWeightsTable() { return weightsTable_; }
protected:
	void computeIntVels(Array2DImpl<float>* topH, 
		Array2DImpl<float> *botH, 
		BufferString name,
		int ihorz,
		Array2DImpl<float>* botDepths,
		Array2DImpl<float>* topDepths,
		Array2DImpl<double>* intVels,
		Array2DImpl<double>* wei,
		Array2DImpl<float>* cou,
		int nrows,
		int ncols,
		bool userMarkers,
		Array2DImpl<double>* weightsTable);
	float adjustTraceVels(float inputTraceVel, 
						  float inputModelVel, 
						  double errorCorrection, 
						  bool calibrateAll, 
						  double weight, 
						  float count);

	
	void normalizeVelGrid(Array2DImpl<double>* weights,
									 Array2DImpl<float>* counts,
									 Array2DImpl<double>* intVels,
									 int nrows,
									 int ncols);
									 
	void addGridPoint(int irow, int icol, 
		Array2DImpl<double>* weights,
		Array2DImpl<float>* counts,
		Array2DImpl<double>* intVels,
		int nrows, 
		int ncols,
		float intVel,
		Array2DImpl<double>* weightsTable);
	float markerDepth(Array2DImpl<float>* horz, BufferString horizonName, MultiID well, int irow, int icol, bool useMarkers,
	    HorizonMarkersList *h);
	float horizonDepth(Array2DImpl<float>* horz, int irow, int icol);
	float horizonTime(Array2DImpl<float>* horz, int irow, int icol);
	int   getNrTraces() { return nrTraces_; }
	float intersectWell(MultiID well, int &irow, int &icol, float &zval,Array2DImpl<float>*  horz);
	void computeWellIntersections(TypeSet<MultiID> selectedWells, int numWells, int ihorz, Coord3 *wellIntersections,
	   Array2DImpl<float>* topH, Array2DImpl<float>* botH );
	void getWellHorizonVelocities(int ihorz, double dt, int nrows, int ncols, int irow, int icol,
	float topH, float botH,Coord3  **wellIntersections, float * veltrace,Array2DImpl<double>* weightsTable );

	template <class T>
	void handleMissing(Array2DImpl<T> *topH, 
					   Array2DImpl<T> *botH, 
					   int nrows,
					   int ncols);

	void expandHorizon(Array2DImpl<float>* outH, 
			Array2D<float> *inH, 
			StepInterval<int> rowrg, 
			StepInterval<int> colrg);

	float * getHorizonDepths(float horizonTimes[], double intervalVels[], int nHorizons);
	float * getHorizonTimes(int irow, int icol);
	double * getHorizonVelocities(int irow, int icol);
	float getTraceAveIntVels(SeisTrc trcin, float traceIntVels[], float topTime, float botTime);
	void scanHorizon(Array2DImpl<float> *horizonTimes, int nrows, int ncols, float &, float &);
	void orderHorizons(
					   float *horizonTmin_, 
					   float *horizontMax_, 
					   int   *horizonSortV_, 
				       int    numHorizons);
	void extendZTHorizons(
					   float *horizonTmin_, 
					   float *horizontMax_, 
					   int   *horizonSortV_, 
				       int    numHorizons);
	void insertat(int * list, int ipos, int ival, int &nrItems);
	void logMessage(BufferString mes) { messages_ += mes; }
	void logMessage(float f) { messages_ += f; }
	void MyUsrMsg(float f) { BufferString s; s+= f; UsrMsg(s);} 
	void addControlPoint(int wrow, int wcol,float intVel, int nrows,int ncols);
	void computeTDCVels(Array2DImpl<float> *topH,Array2DImpl<float> *botH, int ihorz, int nrows, int ncols);
	void initWellTracks();
	

    SeisTrcReader* rdr_;
	SeisTrcWriter* wrr_;
	TypeSet<MultiID> selectedWells_;
	BufferStringSet markerNames_;
	int			numHorizons_;
	TypeSet<MultiID> horizons_;
	int			ihorz_;
	Array2DImpl<float>* topH_;
	Array2DImpl<float>* botH_;
	Array2DImpl<float>* botH0_;
	int			nrTraces_;
	float		radius_;
	float		power_;
	bool		calibrateAll_;
	bool		useMarkers_;
	Array2DImpl<float>*  counts_;
	Array2DImpl<float>*  botDepths_;
	Array2DImpl<float>*  topDepths_;
	bool		firstTime_;
	EM::Horizon3D * tempH_;
	float		timeDatum_;
	float		depthDatum_;	
	int			*horizonSortV_;
	float		*horizonTmin_;
	float		*horizonTmax_;
	BufferString messages_;
	std::vector< Vec > controlPoints_;
	int			gridder_;
	float		relaxation_;
	int numWells_;
	Coord3  **wellIntersections_;
	FILE	*tracef;
	float *** TDCHmodel;
	bool isTDC_;
	Well::Data** wd_;
        Well::Track *timetrack_;
	HorizonMarkersList *hzmarkersList_;
	bool zeroThicknessFound;
	bool unrealVelFound;
	Array2DImpl<double> * weightsTable_;
};



}; //namespace

#endif
