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
// Filename       : $RCSfile: N_IO_OutputResponse.C,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.1 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_Message.h>
#include <N_IO_Op.h>
#include <N_IO_OutputResponse.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace IO {

OutputResponse::OutputResponse(
  OutputMgr &           output_manager)
  : outputManager_(output_manager),
    responseFileName_("response.out"),
    responseVarList_(),
    responseVarPtr_(0),
    responseNames_(),
    numResponseVars_(0)
{}

OutputResponse::~OutputResponse()
{
  for (Util::Op::OpList::iterator it = responseVarList_.begin(); it != responseVarList_.end(); ++it)
    delete *it;
}

//-----------------------------------------------------------------------------
// Function      : OutputResponse::saveResponseVarValues
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 08/11/04
//-----------------------------------------------------------------------------
void
OutputResponse::saveResponseVarValues(
  Parallel::Machine     comm,
  double                time,
  const N_LAS_Vector &  solnVecPtr)
{
  // save the independant variable
  int varNumber = 0;

  responseVarPtr_->at(varNumber) = time;

  // if (outputState_.dcSweepVector_.empty())  // DNS: is this a reliable way to determine if transient?
  // {
  //   responseVarPtr_->at(varNumber) = outputState_.circuitTime_;
  // }
  // else
  // {
  //   responseVarPtr_->at(varNumber) = outputState_.dcSweepVector_[ dcLoopNumber_ ].currentVal;
  // }
  varNumber++;

  // loop over the response variable list
  for (Util::Op::OpList::const_iterator it = responseVarList_.begin(); it != responseVarList_.end(); ++it)
  {
    double result = Util::Op::getValue(comm, *(*it), Util::Op::OpData(0, &solnVecPtr, 0, 0, 0, 0)).real();

    responseVarPtr_->at(varNumber) = result;
    varNumber++;
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputResponse::finalizeResponseVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/04/06
//-----------------------------------------------------------------------------
void
OutputResponse::finalizeResponseVars(
  double        time)
{
  // save the resuts of any end point variables if there are any
  if (numResponseVars_ != 0)
  {
    // save the independant variable
    int varNumber = 0;
    responseVarPtr_->at(varNumber) = time;

    // if (outputState_.dcSweepVector_.empty())  // DNS: is this a reliable way to determine if transient?
    // {
    //   responseVarPtr_->at(v  arNumber) = outputState_.circuitTime_;
    // }
    // else
    // {
    //   responseVarPtr_->at(varNumber) = outputState_.dcSweepVector_[ dcLoopNumber_ ].currentVal;
    // }
    varNumber++;
  }
}


// routines to tell the output manager which variables external programs will
// need as output.  By default we'll only remember the last timepoint
// or dc step unless asked to track all history.

//-----------------------------------------------------------------------------
// Function      : OutputResponse::registerResponseVars
// Purpose       : Create an objective from string submitted by external program
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/04/06
//-----------------------------------------------------------------------------

bool
OutputResponse::registerResponseVars(
  Parallel::Machine     comm,
  const std::string &   objString,
  std::vector<double> * varVectorPtr)
{
  bool result = true;
  responseVarPtr_ = varVectorPtr;

  ExtendedString sVal(objString);
  ParameterList::iterator pl_i;
  Util::Param parameter;

  sVal.toUpper();
  if (sVal.size() < 3 || sVal[1] != '(' || sVal[sVal.size()-1] != ')') {
    Report::DevelFatal0() << "OutputResponse::registerResponseVars: response var not of format V() or I(): '" << objString << "'";
  }
  numResponseVars_++;

  ParameterList pList;

  parameter.setTag(sVal.substr(0, 1));
  parameter.setVal(1.0);
  pList.push_back(parameter);

  parameter.setTag(sVal.substr(2, sVal.size()-3));
  parameter.setVal(0.0);
  pList.push_back(parameter);

  makeOps(comm, outputManager_.getOpBuilderManager(), pList.begin(), pList.end(), std::back_inserter(responseVarList_));

  return result;
}

void
OpCallback<void>::send(
  Parallel::Machine             comm,
  const N_LAS_Vector *          real_solution_vector,
  const N_LAS_Vector *          imaginary_solution_vector,
  const N_LAS_Vector *          state_vector,
  const N_LAS_Vector *          store_vector,
  const std::vector<double> *   sens_function,
  const std::vector<double> *   sens_dOdPDirect,
  const std::vector<double> *   sens_dOdPDirectScaled,
  const std::vector<double> *   sens_dOdPAdjoint,
  const std::vector<double> *   sens_dOdPAdjointScaled) const
{
  execute(Util::Op::getValue(comm, op_, Util::Op::OpData(0, real_solution_vector, imaginary_solution_vector,
                                                         state_vector, store_vector, 0, sens_function,
                                                         sens_dOdPDirect, sens_dOdPDirectScaled, sens_dOdPAdjoint, sens_dOdPAdjointScaled)));
}

void
OutputResponse::send(
  Parallel::Machine             comm,
  const N_LAS_Vector *          real_solution_vector,
  const N_LAS_Vector *          imaginary_solution_vector,
  const N_LAS_Vector *          state_vector,
  const N_LAS_Vector *          store_vector,
  const std::vector<double> *   sens_function,
  const std::vector<double> *   sens_dOdPDirect,
  const std::vector<double> *   sens_dOdPDirectScaled,
  const std::vector<double> *   sens_dOdPAdjoint,
  const std::vector<double> *   sens_dOdPAdjointScaled)
{
  for (OpCallbackVector::const_iterator it = callbacks_.begin(), end = callbacks_.end(); it != end; ++it)
  {
    const OpCallback<void> &callback = *(*it);
    callback.send(comm, real_solution_vector, imaginary_solution_vector,
                  state_vector, store_vector, sens_function,
                  sens_dOdPDirect, sens_dOdPDirectScaled, sens_dOdPAdjoint, sens_dOdPAdjointScaled);
  }
}


} // namespace IO
} // namespace Xyce
