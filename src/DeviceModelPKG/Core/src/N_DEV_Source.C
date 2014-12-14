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
// Filename       : $RCSfile: N_DEV_Source.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.30 $
//
// Revision Date  : $Date: 2014/04/03 19:19:24 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// N_DEV_Source
#include <N_DEV_Source.h>

// Class N_DEV_Source
#include <N_DEV_SourceData.h>

#include <N_DEV_SolverState.h>

namespace Xyce {
namespace Device {

// SourceInstance
//-----------------------------------------------------------------------------
// Function      : SourceInstance::SourceInstance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SourceInstance::SourceInstance(
  const InstanceBlock &               IB,
  ParametricData<void> &        parametric_data,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, parametric_data, factory_block),
   tranSourceData_(0),
   acSourceData_(0),
   dcSourceData_(0),
   sourceType(_SIN_DATA)
{}


//-----------------------------------------------------------------------------
// Function      : SourceInstance::~SourceInstance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SourceInstance::~SourceInstance()
{}


//-----------------------------------------------------------------------------
// Function      : SourceInstance::setFastSourceFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/27/04
//-----------------------------------------------------------------------------
void SourceInstance::setFastSourceFlag (bool value)
{
  if (tranSourceData_ != 0)
  {
    tranSourceData_->setFastTimeScaleFlag(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::getFastSourceFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/1/04
//-----------------------------------------------------------------------------
bool SourceInstance::getFastSourceFlag()
{
  bool flag=false;
  if (tranSourceData_ != 0)
  {
    flag = tranSourceData_->getFastTimeScaleFlag();
  }
  return flag;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::period
// Purpose       : Time period of periodic sources
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/31/04
//-----------------------------------------------------------------------------
double SourceInstance::period()
{
  double per = 0.0;
  if (tranSourceData_ != 0)
  {
    per = tranSourceData_->period();
  }
  return per;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::getResetFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/30/04
//-----------------------------------------------------------------------------
bool SourceInstance::getResetFlag ()
{
  if (tranSourceData_ != 0)
  {
    tranSourceData_->getResetFlag ();
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::getInstanceBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/05/06
//-----------------------------------------------------------------------------
bool SourceInstance::getInstanceBreakPoints
   (std::vector<N_UTL_BreakPoint> &breakPointTimes)
{
  bool fastSourceFlag = false;

  if (getSolverState().blockAnalysisFlag == false)
  {
    fastSourceFlag = getFastSourceFlag();
  }

  bool tmpBool = true;
  if (!fastSourceFlag && tranSourceData_ != 0)
  {
    tmpBool = tranSourceData_->getBreakPoints(breakPointTimes);
  }
  return tmpBool;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/05/06
//-----------------------------------------------------------------------------
bool SourceInstance::updateSource ()
{
  if (tranSourceData_ != 0)
  {
    tranSourceData_->updateSource ();
  }
  if (dcSourceData_ != 0)
  {
    dcSourceData_->updateSource ();
  }
  if (acSourceData_ != 0)
  {
    acSourceData_->updateSource ();
  }
  return true;
}

} // namespace Device
} // namespace Xyce
