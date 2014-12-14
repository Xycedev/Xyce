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
// Filename       : $RCSfile: N_IO_Outputter.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.68.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <list>
#include <string>
#include <vector>

#include <N_IO_Outputter.h>
#include <N_ANP_fwd.h>

#include <N_UTL_Demangle.h>
#include <N_UTL_NetlistLocation.h>
#include <N_UTL_Param.h>

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace IO {
namespace Outputter {

const static int debug = false;

void
Interface:: setAnalysisMode(
  Analysis::Analysis_Mode       analysis_mode)
{
  doSetAnalysisMode(analysis_mode);
}

void
Interface::output(
  Parallel::Machine             comm,
  const N_LAS_Vector &          solution_vector,
  const N_LAS_Vector &          state_vector,
  const N_LAS_Vector &          store_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutput" << std::endl;

  doOutputTime(comm, solution_vector, state_vector, store_vector);
}

void
Interface::outputAC(
  Parallel::Machine             comm,
  double                        frequency,
  const N_LAS_Vector &          real_solution_vector,
  const N_LAS_Vector &          imaginary_solution_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputAC" << std::endl;

  doOutputFrequency(comm, frequency, real_solution_vector, imaginary_solution_vector);
}

void Interface::outputHB(
  Parallel::Machine             comm,
  const std::vector<double> &   timePoints,
  const std::vector<double> &   freqPoints,
  const N_LAS_BlockVector &     timeDomainSolutionVec,
  const N_LAS_BlockVector &     freqDomainSolutionVecReal,
  const N_LAS_BlockVector &     freqDomainSolutionVecImaginary,
  const N_LAS_BlockVector &     timeDomainStoreVec,
  const N_LAS_BlockVector &     freqDomainStoreVecReal,
  const N_LAS_BlockVector &     freqDomainStoreVecImaginary)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputHB" << std::endl;

  doOutputHB(comm, timePoints, freqPoints,
             timeDomainSolutionVec, freqDomainSolutionVecReal, freqDomainSolutionVecImaginary,
             timeDomainStoreVec, freqDomainStoreVecReal, freqDomainStoreVecImaginary);
}

void Interface::outputMPDE(
  Parallel::Machine             comm,
  double                        time,
  const std::vector<double> &   fast_time_points,
  const N_LAS_Vector &          solution_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputMPDE" << std::endl;

  doOutputMPDE(comm, time, fast_time_points, solution_vector);
}

void Interface::outputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const N_LAS_Vector &                  solution_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputHomotopy" << std::endl;

  doOutputHomotopy(comm, parameter_names, parameter_values, solution_vector);
}

void Interface::outputSensitivity(
  Parallel::Machine             comm,
  const std::vector<double> &   objective_values,
  const std::vector<double> &   direct_values,
  const std::vector<double> &   adjoint_values,
  const std::vector<double> &   scaled_direct_values,
  const std::vector<double> &   scaled_adjoint_values,
  const N_LAS_Vector &          solution_vector,
  const N_LAS_Vector &          state_vector,
  const N_LAS_Vector &          store_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputSensitivity" << std::endl;

  doOutputSensitivity(
    comm, objective_values,
    direct_values, adjoint_values,
    scaled_direct_values, scaled_adjoint_values,
    solution_vector, state_vector, store_vector);
}

void
Interface::finishOutput()
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doFinishOutput" << std::endl;

  doFinishOutput();
}

void
Interface::startStep(
  int                           step,
  int                           max_step)
{
  doStartStep(step, max_step);
}

void
Interface::steppingComplete()
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doSteppingComplete" << std::endl;

  doSteppingComplete();
}


} // namespace Outputter
} // namespace IO
} // namespace Xyce
