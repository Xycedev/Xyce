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
// Filename       : $RCSfile: N_DEV_CharonInterface.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/13/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.24 $
//
// Revision Date  : $Date: 2014/08/05 21:41:37 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_CharonInterface_h
#define Xyce_N_DEV_CharonInterface_h

#include <string>
#include <vector>

#include "Teuchos_RefCountPtr.hpp"

#include <N_DEV_fwd.h>
#include <N_CIR_fwd.h>
#include <N_DEV_ExternCodeInterface.h>
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>

namespace Teuchos {
  class ParameterList;
}

class N_TIA_TimeIntInfo;
class N_TIA_TwoLevelError;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : CharonInterface
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
class CharonInterface : public ExternCodeInterface
{
  public:
    CharonInterface (
      const DeviceOptions &     do1,
      const std::string &       netlist,
      const SolverState &       ss1);

    CharonInterface (const CharonInterface &right);
    virtual ~CharonInterface();

    bool initialize(N_PDS_Comm * comm = 0);

    bool simulateStep
      ( const SolverState & solState,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError
      );

    bool finalize ();
    bool run ();

    void homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals);

    void homotopyStepFailure ();

    void stepSuccess(Analysis::CurrentMode analysis);
    void stepFailure(Analysis::CurrentMode analysis);
    bool getBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes) { return true; }
    bool updateStateArrays () {return true;}
    bool startTimeStep  ( const N_TIA_TimeIntInfo & tiInfo ) {return true;}
    bool setInternalParam (std::string & name, double val) {return true;}

    bool getInitialQnorm (N_TIA_TwoLevelError & tle);

  private:

    std::string inputFileName_;
    const DeviceOptions& devOptions_;
    const SolverState & solState_;

    // The "command line" arguments
    std::vector<char*> cargs_;

    //! Input list for Charon.
    Teuchos::RefCountPtr<Teuchos::ParameterList> input_list_;

    //! Output list from Charon.
    Teuchos::RefCountPtr<Teuchos::ParameterList> output_list_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::CharonInterface N_DEV_CharonInterface;

#endif

