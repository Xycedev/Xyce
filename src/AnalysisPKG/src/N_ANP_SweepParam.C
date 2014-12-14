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
// Filename      : $RCSfile: N_ANP_SweepParam.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 9/4/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.15.2.1 $
// Revision Date  : $Date: 2014/08/26 22:31:06 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <N_ANP_SweepParam.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>


namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : SweepParam::updateCurrentVal
//
// Purpose       : Updates the values of the parameters used in a sweep.
//
// Special Notes : This is very similar to the "update" function in the
//                 class N_DEV_SweepData.  (which no longer exists).
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 10/31/2003
//-----------------------------------------------------------------------------
bool SweepParam::updateCurrentVal (int stepNumberArg)
{
  outerStepNumber  = stepNumberArg/interval;
  int inum            = outerStepNumber/maxStep;
  int localStepNumber = outerStepNumber - inum*maxStep;

  // We must keep track of whether we're at the first step of our sweep,
  // because some device manager features need to know that.
  // It is important that we only set this when localStepNumber first becoms
  // zero, not every time localStepNumber *is* zero, because an outer loop
  // might remain at 0 for quite some time.

  if (localStepNumber == 0 && localStepNumber != lastLocalStepNumber_)
  {
    sweepResetFlag_=true;
  } 
  else
  {
    sweepResetFlag_=false;
  }
  lastLocalStepNumber_=localStepNumber;

  if (type == "LIN")
  {
    currentVal = startVal + static_cast<double>(localStepNumber)*stepVal;
    ++count;
  }
  else if (type == "DEC" || type == "OCT")
  {
    currentVal = startVal*pow(stepMult, static_cast<double>(localStepNumber) );
    ++count;
  }
  else if (type == "LIST")
  {
    int size=  valList.size();
    int index = (localStepNumber < size)?localStepNumber:(size-1);
    currentVal = valList[index];
    ++count;
  }
  else
  {
    Report::DevelFatal0().in("SweepParam::updateCurrentVal") << "Unsupported type specified";
  }
#ifdef Xyce_DEBUG_TIME
  Xyce::dout() << std::endl
               << Xyce::subsection_divider << std::endl
               << "updateCurrentVal" << std::endl
               << "  name             = " << name << std::endl
               << "  stepNumberArg    = " << stepNumberArg<< std::endl
               << "  interval         = " << interval  << std::endl
               << "  outerStepNumber  = " << outerStepNumber << std::endl
               << "  localStepNumber  = " << localStepNumber << std::endl
               << "  inum             = " << inum      << std::endl
               << "  sweepResetFlag   = " << sweepResetFlag_ << std::endl
               << "  currentVal       = " << currentVal << std::endl
               << Xyce::subsection_divider << std::endl;
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SweepParam::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
std::ostream &
operator<<(std::ostream & os, const SweepParam & sp)
{
  os << "\tname            = " << sp.name
     << "\tcurrentVal      = " << sp.currentVal
     << std::endl;
  return os;
}

} // namespace Analysis
} // namespace Xyce

