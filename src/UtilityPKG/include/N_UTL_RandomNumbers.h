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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_RandomNumbers.h,v $
//
// Purpose        : Provide a class with methods to generate random numbers
//                  with specific properties
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL
//
// Creation Date  :  5/28/2014
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/05/28 18:40:03 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_UTL_RandomNumbers_H
#define N_UTL_RandomNumbers_H

namespace Xyce {
namespace Util {

class RandomNumbers
{
 public:
  RandomNumbers(long seed=0)   
    : randInitialized_(false),
      useLastGaussian_(false)
    {
      seedRandom(seed);
    }

  double uniformRandom();
  double gaussianRandom(double mu, double sigma);
  void seedRandom(long seed);

 private:
  bool randInitialized_;
#ifndef HAVE_DRAND48
  // These bits of data only needed if we have no good random number generator
  double maxran_;
  int randBuf_[98];
  double randY_;
#endif

  // Internals for use by gaussian random number generator
  double ySaveGaussian_;
  bool useLastGaussian_;
};
}
}
#endif //N_UTL_RandomNumbers_H
