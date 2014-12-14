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
// Filename       : $RCSfile: N_ANP_MPDE.h,v $
//
// Purpose        : MPDE analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Todd Coffey, 1414, Ting Mei 1437
//
// Creation Date  : 07/23/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.17 $
// Revision Date  : $Date: 2014/07/23 20:30:56 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_MPDE_h
#define Xyce_N_ANP_MPDE_h

// ----------   Xyce Includes   ----------
#include <N_ANP_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_StepEvent.h>
#include <N_UTL_Listener.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : MPDE
// Purpose       : MPDE analysis class
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class MPDE : public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
public:
    MPDE(AnalysisManager &analysis_manager);
    
    virtual ~MPDE() {};
    
    void notify(const StepEvent &event);

    bool getDCOPFlag();

    bool run(); 
    bool init(); 
    bool loopProcess(); 
    bool processSuccessfulDCOP(); 
    bool processFailedDCOP(); 
    bool processSuccessfulStep();
    bool processFailedStep();
    bool finish();
    bool handlePredictor();

    bool finalVerboseOutput();

private:
 
  void takeAnIntegrationStep_();

  // a flag to indicate of the simulation is paused
  bool isPaused; 

  double startDCOPtime, endTRANtime; // startTRANtime
  // Timing/loop count info
  // Device::DeviceInterface *             devInterfacePtr_;
  // N_TOP_Topology *                      topoMgrPtr_;
  // N_LOA_NonlinearEquationLoader *       nonlinearEquationLoaderPtr_;
  // N_LAS_Builder *                       appBuilderPtr_;
};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::MPDE N_ANP_MPDE;

#endif // Xyce_N_ANP_MPDE_h

