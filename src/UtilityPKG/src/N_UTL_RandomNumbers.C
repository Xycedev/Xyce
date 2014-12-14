//-----------------------------------------------------------------------------
// Copyright Notice
// 
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_RandomNumbers.C,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Tom Russo, SNL
//
// Creation Date : 5/28/2014
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/05/28 18:40:03 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>
#include <N_UTL_RandomNumbers.h>

#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>
#endif


#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <time.h>

namespace Xyce {
namespace Util {
//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::gaussianRandom
// Purpose       : Provide standardized, portable source of high-quality
//                 random numbers, normally distributed with mean mu and
//                 standard deviation sigma.
// Special Notes : Method is a variant of the Box-Muller transformation.
//                 A pair of random numbers from a uniform distribution is
//                 selected such that they can be the coordinates of a
//                 unit vector.  The magnitude and phase of this vector
//                 are used in the Box-Muller transformation to create a 
//                 set of values that are normally distributed instead of
//                 uniformly distributed.
// Creator       : Tom Russo
// Creation Date : 05/27/2014
//-----------------------------------------------------------------------------
double RandomNumbers::gaussianRandom(double mu, double sigma)
{
  double x1, x2, w, y1;
  
  if (useLastGaussian_)      /* use value from previous call */
  {
    y1 = ySaveGaussian_;
    useLastGaussian_ = false;
  }
  else
  {
    do 
    {
      x1 = 2.0 * uniformRandom() - 1.0;
      x2 = 2.0 * uniformRandom() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    ySaveGaussian_ = x2 * w;
    useLastGaussian_ = true;
  }
  
  return( mu + y1 * sigma );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::uniformRandom
// Purpose       : Provide standardized, portable source of high-quality
//                 random numbers, uniformly distributed on [0,1].
// Special Notes : The code used when we don't have drand48 is based on
//                 code from Numerical Recipes, in the "Improving a bad
//                 random number generator" section.
//                 This code assumes the random number generator is already
//                 seeded.
// Creator       : Tom Russo
// Creation Date : 05/21/2014
//-----------------------------------------------------------------------------
double RandomNumbers::uniformRandom ()
{

  double res;
  
#ifdef HAVE_DRAND48
  // The rand48 package is available on most Unix-like systems, and provides
  // a high-quality random numbers.  Just use it.
  res=drand48();
#else
  
  // Otherwise, all we can count on having is the generic "rand"
  // function, which is well known to be horrid.  Let's use it to create a 
  // better choice of random numbers.
  double dum;
  int j;

  if (!randInitialized_)
  {
    randInitialized_=true;
    maxran_=RAND_MAX+1.0;
    for (j=0;j<98;++j) dum=rand();
    for (j=0;j<98;++j) randBuf_[j]=rand();
    randY_=rand();
  }

  //  j will be a number 0-97
  j=98.0*randY_/maxran_;
  randY_=randBuf_[j];
  randBuf_[j]=rand();
  res = randY_/maxran_;
#endif

  return res;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::RandomNumbers::seedRandom
// Purpose       : Initialize a random number generator
// Special Notes : 
// Creator       : Tom Russo
// Creation Date : 05/21/2014
//-----------------------------------------------------------------------------
void RandomNumbers::seedRandom(long seed)
{
  if (seed==0)
    seed=time(NULL);

#ifdef HAVE_DRAND48
  srand48(seed);
#else
  srand(seed);
#endif

}
}
}
