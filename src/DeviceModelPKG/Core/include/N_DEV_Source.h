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
// Filename       : $RCSfile: N_DEV_Source.h,v $
//
// Purpose        : Source base  lasses.
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
// Revision Number: $Revision: 1.26 $
//
// Revision Date  : $Date: 2014/04/03 19:19:24 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Source_h
#define Xyce_N_DEV_Source_h

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : SourceInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
class SourceInstance : public DeviceInstance
{
public:
  SourceInstance(
     const InstanceBlock &       IB,
     ParametricData<void> &      parametric_data,
     const FactoryBlock &        factory_block);
    
  ~SourceInstance();

private:
  SourceInstance(const SourceInstance &);
  SourceInstance &operator=(const SourceInstance &);

public:
  // Additional Public Declarations
  void setFastSourceFlag (bool value);
  bool getFastSourceFlag ();

  bool getResetFlag ();

  double period();

  virtual bool getInstanceBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes);
  virtual bool updateSource ();

  virtual bool loadBVectorsforAC(double * bVecReal, double * bVecImag ) {
    return true;
  }

  protected:
  // Additional Protected Declarations
  int sourceType;  // type of source data
  SourceData *tranSourceData_;
  SourceData *acSourceData_;
  SourceData *dcSourceData_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::SourceInstance N_DEV_SourceInstance;

#endif
