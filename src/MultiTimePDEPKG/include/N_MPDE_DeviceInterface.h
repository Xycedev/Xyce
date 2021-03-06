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
// Filename       : $RCSfile: N_MPDE_DeviceInterface.h,v $
//
// Purpose        : MPDE Specific Loader
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 07/21/2008
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/08/07 23:08:54 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_MPDE_DeviceInterface_H
#define Xyce_N_MPDE_DeviceInterface_H

// ---------- Standard Includes ----------
#include <string>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_ANP_fwd.h>
#include <N_NLS_fwd.h>

// ---------- Forward declarations --------
class N_LAS_System;


//-----------------------------------------------------------------------------
// Class         : N_MPDE_DeviceInterface
// Purpose       :
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 07/21/2008
//-----------------------------------------------------------------------------
class N_MPDE_DeviceInterface
{

public:
  N_MPDE_DeviceInterface();

  void registerDeviceInterface( N_DEV_DeviceInterface *devIntPtr );
  void registerLinearSystem(N_LAS_System *tmp_system_ptr);
  void registerAnalysisManager(N_ANP_AnalysisManager *tmp_anaIntPtr);
  void registerNonlinearSolver (N_NLS_Manager *tmp_nlsMgrPtr);

  std::vector<double> getFastSourcePeriod (std::vector<std::string>& sourceNames);
  std::vector<double> registerFastSources (std::vector<std::string> & sourceNames);

  void deactivateSlowSources();
  void activateSlowSources();

  void setMPDEFlag( bool flagVal );
  void setBlockAnalysisFlag( bool flagVal );
  void setFastTime( double timeVal );

private:
  N_DEV_DeviceInterface *       devInterfacePtr_;

};

#endif

