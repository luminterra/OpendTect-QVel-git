#ifndef mywellhorpos_h
#define mywellhorpos_h

/*+
________________________________________________________________________

 (C) dGB Beheer B.V.; (LICENSE) http://opendtect.org/OpendTect_license.txt
 Author:        Bruno
 Date:          Jul 2010
 RCS:           $Id: wellhorpos.h 33397 2014-02-19 09:11:02Z mahant.mothey@dgbes.com $
________________________________________________________________________

-*/

#include "wellattribmod.h"
#include "emposid.h"
#include "namedobj.h"
#include "position.h"

namespace Well { class Track; class D2TModel; }
namespace EM { class Horizon2D; class Horizon3D; }

/*!
\brief Used to give well/horizon intersection. In theory more than one
intersection is possible ( in case of faults or deviated tracks along the
horizon ) but only one pos will be returned.
*/
static FILE *tracef = fopen("/temp/tracewh.out","w");
mExpClass(WellAttrib) myWellHorIntersectFinder
{
public:
    				myWellHorIntersectFinder(const Well::Track&,
						const Well::D2TModel* d2t=0);
				//! a d2t model is needed if z is time 

    void			setHorizon(const EM::ObjectID& emid);

    float			findZIntersection() const;
   				//return undef if not found

protected:

    const Well::Track&		track_;
    const Well::D2TModel*	d2t_;
    const EM::Horizon2D*	hor2d_;
    const EM::Horizon3D*	hor3d_;

    float			intersectPosHor(const Coord3&) const;
    
};

#endif

