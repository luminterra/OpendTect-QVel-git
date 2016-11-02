/*+
________________________________________________________________________

 (C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
 Author:        Bruno
 Date:          Jul 2010
________________________________________________________________________

-*/
static const char* rcsID mUsedVar = "$Id: wellhorpos.cc 32104 2013-10-23 20:11:53Z kristofer.tingdahl@dgbes.com $";


#include "mywellhorpos.h"

#include "emhorizon2d.h"
#include "emhorizon3d.h"
#include "emmanager.h"
#include "survinfo.h"
#include "welld2tmodel.h"
#include "welltrack.h"


myWellHorIntersectFinder::myWellHorIntersectFinder( const Well::Track& tr, 
						const Well::D2TModel* d2t )
    : track_(tr)
    , d2t_(d2t)  
{ 
}


void myWellHorIntersectFinder::setHorizon( const EM::ObjectID& emid )
{
    hor2d_ = 0; hor3d_ = 0;
    const EM::EMObject* emobj = EM::EMM().getObject( emid );
    mDynamicCastGet(const EM::Horizon2D*,hor2d,emobj)
    mDynamicCastGet(const EM::Horizon3D*,hor3d,emobj)
    hor3d_ = hor3d; hor2d_ = hor2d;
}


float myWellHorIntersectFinder::findZIntersection() const
{
    const float zstep = SI().zStep();
    const Interval<float>& dahrg = track_.dahRange();
    float zstart = d2t_ ? d2t_->getTime( dahrg.start, track_ ) : dahrg.start;
    float zstop = d2t_ ? d2t_->getTime( dahrg.stop, track_ ) : dahrg.stop;
	fprintf(tracef,"\nBegin intersection dahrg %g %g z %g %g\n",dahrg.start,dahrg.stop,zstart,zstop);
    zstart = mMAX( SI().zRange(true).start, zstart );
    zstop = mMIN( SI().zRange(true).stop, zstop );

    float zval = zstart; 
    bool isabove = true;
    bool firstvalidzfound = false;
	
    while ( zval < zstop )
    {
	
	const float dah = d2t_ ? d2t_->getDah( zval, track_ ) : zval;
	const Coord3& crd = track_.getPos( dah );
	const float horz = intersectPosHor( crd );
	fprintf(tracef,"zval %g zstop %g  horizon t %g\n",zval,zstop,horz);
	if ( mIsUdf( horz ) ) 
	{
	    zval += zstep;
	    fprintf(tracef,"  miss, new zval %g\n",zval);
	    continue;
	}
fprintf(tracef,"   hit, new zval %g\n",zval);
	if ( !firstvalidzfound )
	{
	    isabove = zval >= horz;
	    firstvalidzfound = true;
		fprintf(tracef,"first valid z\n");
	}

	if ( ( isabove && horz >= zval ) || ( !isabove && horz <= zval ) ) {
		fprintf(tracef,"return horz %g\n",horz);
	
	    return horz;
	}
	  zval += zstep;
	   fprintf(tracef,"   skip, new zval %g\n",zval);
    }
    fprintf(tracef,"miss!, new zval %g\n",zval);
    return mUdf( float );
}


float myWellHorIntersectFinder::intersectPosHor( const Coord3& pos ) const 
{
    const BinID& bid = SI().transform( pos );
    if ( !SI().isInside( bid, true ) ) {
	fprintf(tracef,"pos is NOT not inside %g %g %g\n",pos.x, pos.y, pos.z);
       return mUdf( float );	
    }

    if ( hor3d_ )
    {
	const EM::SubID subid = bid.toInt64();
	const Coord3& horpos = hor3d_->getPos( hor3d_->sectionID(0), subid ); 
	fprintf(tracef,"horpos %g %g %g\n",horpos.x, horpos.y, horpos.z);
	const BinID horbid = SI().transform( horpos );
	fprintf(tracef,"horpos %g %g %g from pos %g %g\n",horpos.x, horpos.y, horpos.z,pos.x, pos.y);
	if ( bid == horbid ) {
	    fprintf(tracef,"horizon found %g\n",horpos.z);
	    return (float)horpos.z;
	} else
	    fprintf(tracef," beid != horbid\n");
    }
    else if ( hor2d_ )
    {
	mDynamicCastGet( const Geometry::RowColSurface*, rcs, 
			 hor2d_->sectionGeometry(0));
	if ( !rcs ) return mUdf(float);

	const StepInterval<int> rowrg = rcs->rowRange();
	RowCol rc;
	for ( rc.row()=rowrg.start; rc.row()<=rowrg.stop; rc.row()+=rowrg.step )
	{
	    const StepInterval<int> colrg = rcs->colRange( rc.row() );
	    for ( rc.col()=colrg.start; rc.col()<=colrg.stop; rc.col()+=colrg.step )
	    {
		const Coord3& horpos = hor2d_->getPos( 
					hor2d_->sectionID(0), rc.toInt64() );
		const BinID horbid = SI().transform( horpos );
		if ( bid == horbid )
		    return (float) pos.z;
	    }
	}
    }
    fprintf(tracef,"mis hor\n");
    return  mUdf( float );
}
