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
// Filename      : $RCSfile: N_IO_MeasureManager.h,v $
//
// Purpose       : This file is a class to manage measure statements in a sim.
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26.2.1 $
//
// Revision Date  : $Date: 2014/08/25 20:12:50 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef  Xyce_N_IO_MeasureManager_H
#define Xyce_N_IO_MeasureManager_H

#include <list>
#include <string>
#include <iostream>

#include <N_ANP_SweepParam.h>
#include <N_IO_Measure_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_Vector.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_Listener.h>
#include <N_ANP_StepEvent.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Class         : MeasureManager
// Purpose       : This is a manager class for handling measure statements
//                 in a simulation
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
class Manager : public Util::Listener<Analysis::StepEvent>
{
  typedef std::vector<Base *> MeasurementVector;

public:
  Manager(const std::string &netlist_filename);
  ~Manager();

  bool registerPkgOptionsMgr(PkgOptionsMgr &pkgOpt);

  void notify(const Analysis::StepEvent &step_event);

  // Return true if .measure analysis is being performed on any variables.
  bool isMeasureActive() const { return (!allMeasuresList_.empty()); }

  // add .measure line from netlist to list of things to measure
  bool addMeasure(const Util::OptionBlock & measureLine);

  void makeMeasureOps(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager);

  // Called during the simulation to update the measure objects held by this class
  // To keep things obvious, use one call for tran and another for DC
  void updateTranMeasures(Parallel::Machine comm, const double circuitTime, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
  void updateDCMeasures(Parallel::Machine comm, const std::vector<N_ANP_SweepParam> & dcParamsVec, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, const N_LAS_Vector *storeVec);
  void updateACMeasures(Parallel::Machine comm, const double frequency, const N_LAS_Vector *solnVec, const N_LAS_Vector *imaginaryVec);

  std::ostream & outputResults(std::ostream& outputStream, bool printHeader=true ) const;

  const Base *find(const std::string &name) const;

  void remeasure(N_PDS_Comm &pds_comm, const std::string &netlist_filename, const std::string &remeasure_path, bool transient_flag, OutputMgr &output_manager);

  void getMeasureValue (const std::string &name, double &value, bool &found) const;

private:
  std::string           netListFilename_;
  MeasurementVector     allMeasuresList_;
  MeasurementVector     activeMeasuresList_;
};

} // namespace Measure
} // namespace IO
} // namespace Xyce

typedef Xyce::IO::Measure::Manager N_IO_MeasureManager;

#endif  // Xyce_N_IO_MeasureManager_H
