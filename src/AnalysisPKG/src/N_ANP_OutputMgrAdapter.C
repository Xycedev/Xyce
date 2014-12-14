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
// Filename       : $RCSfile: N_ANP_OutputMgrAdapter.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:48 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_OutputMgrAdapter.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_OutputROM.h>
#include <N_IO_OutputMOR.h>
#include <N_IO_OutputResults.h>
#include <N_IO_OutputResponse.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::OutputMgrAdapter( )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
OutputMgrAdapter::OutputMgrAdapter(
  Parallel::Machine             comm,
  Util::Notifier<StepEvent> &   step_notifier)
  : Util::ListenerAutoSubscribe<StepEvent>(step_notifier),
  comm_(comm),
  outputManager_(0),
  measureManager_(0),
  fourierManager_(0),
  outputMOR_(0),
  outputResults_(0),
  outputResponse_(0),
  emptyParamVector_(),
  stepParamVector_(&emptyParamVector_),
  dcParamVector_(),
  stepAnalysisStepNumber_(0),
  stepAnalysisMaxSteps_(0),
  dcAnalysisStepNumber_(0),
  dcAnalysisMaxSteps_(0)
{}

OutputMgrAdapter::~OutputMgrAdapter()
{
  delete outputMOR_;
  delete outputResults_;
}


void
OutputMgrAdapter::notify(
  const StepEvent &     event)
{
  if (event.state_ == StepEvent::STEP_STARTED)
    stepAnalysisStepNumber_ = event.count_;
  else if (event.state_ == StepEvent::INITIALIZE)
    stepAnalysisMaxSteps_ = event.count_;
}


void OutputMgrAdapter::addOutputResults(const Util::OptionBlock & option_block)
{
  if (!outputResults_)
    outputResults_ = new IO::OutputResults(outputManager_->getNetListFilename());

  outputResults_->addResultParams(option_block);
}

void OutputMgrAdapter::addOutputResponse(const Util::OptionBlock & option_block)
{
  if (!outputResponse_)
    outputResponse_ = new IO::OutputResponse(*outputManager_);
}

void
OutputMgrAdapter::setStepParamVec(
  const std::vector<SweepParam> &       paramVec)
{
  stepParamVector_ = &paramVec;
}

void
OutputMgrAdapter::setDCParamVec(
  const std::vector<SweepParam> &       paramVec)
{
  dcParamVector_ = paramVec;
}

void
OutputMgrAdapter::tranOutput(
  double                time,
  N_LAS_Vector &        currSolutionPtr,
  N_LAS_Vector &        stateVecPtr,
  N_LAS_Vector &        storeVecPtr,
  std::vector<double> & objectiveVec_,
  std::vector<double> & dOdpVec_,
  std::vector<double> & dOdpAdjVec_,
  std::vector<double> & scaled_dOdpVec_,
  std::vector<double> & scaled_dOdpAdjVec_,
  bool                  skipPrintLineOutput)
{
  outputManager_->output(
    comm_, time, stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVector_,
    dcAnalysisStepNumber_, dcAnalysisMaxSteps_, dcParamVector_,
    currSolutionPtr, stateVecPtr, storeVecPtr, objectiveVec_,
    dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_,
    skipPrintLineOutput);

  if (outputResponse_) {
    outputResponse_->saveResponseVarValues(comm_, time, currSolutionPtr);
    outputResponse_->send(comm_,
                          &currSolutionPtr, 0, &stateVecPtr, &storeVecPtr,
                          &objectiveVec_, &dOdpVec_, &dOdpAdjVec_, &scaled_dOdpVec_, &scaled_dOdpAdjVec_);
  }
}

void
OutputMgrAdapter::dcOutput(
  int                   dcStepNumber,
  N_LAS_Vector &        currSolutionPtr,
  N_LAS_Vector &        stateVecPtr,
  N_LAS_Vector &        storeVecPtr,
  std::vector<double> & objectiveVec_,
  std::vector<double> & dOdpVec_,
  std::vector<double> & dOdpAdjVec_,
  std::vector<double> & scaled_dOdpVec_,
  std::vector<double> & scaled_dOdpAdjVec_)
{
  outputManager_->output(
    comm_, 0.0,  stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVector_,
    dcStepNumber, dcAnalysisMaxSteps_, dcParamVector_,
    currSolutionPtr, stateVecPtr, storeVecPtr, objectiveVec_,
    dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

  if (outputResponse_)
    outputResponse_->saveResponseVarValues(comm_, 0.0, currSolutionPtr);
}


void
OutputMgrAdapter::outputResult(
  const N_LAS_Vector &  currSolutionPtr,
  const N_LAS_Vector &  currStatePtr,
  const N_LAS_Vector &  currStorePtr)
{
  if (!outputResults_)
    outputResults_ = new IO::OutputResults(outputManager_->getNetListFilename());

  outputResults_->output(comm_, *outputManager_, outputManager_->getCircuitTime(), *stepParamVector_, stepAnalysisStepNumber_, currSolutionPtr, currStatePtr, currStorePtr);
}

void
OutputMgrAdapter::steppingComplete()
{
  outputManager_->steppingComplete();
  if (outputResults_)
    outputResults_->steppingComplete();
}

void
OutputMgrAdapter::finishOutput()
{
  outputManager_->finishOutput();
}

bool
OutputMgrAdapter::setupInitialConditions (
  N_LAS_Vector &        solnVec,
  N_LAS_Vector &        flagVec)
{
  return outputManager_->setupInitialConditions(comm_, solnVec, flagVec);
}

void
OutputMgrAdapter::outputDCOP(
  const N_LAS_Vector &  solution)
{
  outputManager_->outputDCOP(comm_, solution);
}


void
OutputMgrAdapter::outputMPDE (
  double                        time,
  const std::vector<double> &   fast_time_points,
  const N_LAS_Vector &          solution_vector)
{
  outputManager_->outputMPDE(comm_, time, fast_time_points, solution_vector);
}

void
OutputMgrAdapter::outputHB (
  const std::vector< double > & timePoints,
  const std::vector< double > & freqPoints,
  const N_LAS_BlockVector &     timeDomainSolutionVec,
  const N_LAS_BlockVector &     freqDomainSolutionVecReal,
  const N_LAS_BlockVector &     freqDomainSolutionVecImaginary,
  const N_LAS_BlockVector &     timeDomainStoreVec,
  const N_LAS_BlockVector &     freqDomainStoreVecReal,
  const N_LAS_BlockVector &     freqDomainStoreVecImaginary)
{
  outputManager_->outputHB(
    comm_,
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVector_,
    timePoints, freqPoints,
    timeDomainSolutionVec, freqDomainSolutionVecReal, freqDomainSolutionVecImaginary, 
    timeDomainStoreVec, freqDomainStoreVecReal, freqDomainStoreVecImaginary);
}

void
OutputMgrAdapter::outputAC (
  double                freq,
  const N_LAS_Vector &  solnVecRealPtr,
  const N_LAS_Vector &  solnVecImaginaryPtr)
{
  outputManager_->outputAC(comm_, freq, solnVecRealPtr, solnVecImaginaryPtr);
}

void
OutputMgrAdapter::outputMORTF (
  bool                                                          origSys,
  double                                                        freq,
  const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H )
{
  if (!outputMOR_)
    outputMOR_ = new IO::OutputMOR(outputManager_->getNetListFilename());

  outputMOR_->output(comm_, origSys, freq, H );
}

void
OutputMgrAdapter::resetOutputMORTF()
{
  if (!outputMOR_)
    outputMOR_ = new IO::OutputMOR(outputManager_->getNetListFilename());

  outputMOR_->reset();
}

void
OutputMgrAdapter::outputROM(
  const Teuchos::SerialDenseMatrix<int, double> &       Ghat,
  const Teuchos::SerialDenseMatrix<int, double> &       Chat,
  const Teuchos::SerialDenseMatrix<int, double> &       Bhat,
  const Teuchos::SerialDenseMatrix<int, double> &       Lhat )
{
  IO::outputROM(comm_, outputManager_->getNetListFilename(), Ghat, Chat, Bhat, Lhat );
}

void
OutputMgrAdapter::outputROM(
  const N_LAS_Matrix &                                  Ghat,
  const N_LAS_Matrix &                                  Chat,
  const Teuchos::SerialDenseMatrix<int, double> &       Bhat,
  const Teuchos::SerialDenseMatrix<int, double> &       Lhat )
{
  IO::outputROM(comm_, outputManager_->getNetListFilename(),  Ghat, Chat, Bhat, Lhat );
}

bool
OutputMgrAdapter::getOutputIntervals(
  double &                                      initialInterval,
  std::vector<std::pair< double, double > > *   intervalPairs) const
{
  return outputManager_->getOutputIntervals( initialInterval, *intervalPairs );
}


void
OutputMgrAdapter::outputHomotopy(
  const std::vector<std::string> &      paramNames,
  const std::vector<double> &           paramVals,
  N_LAS_Vector &                        solnVecPtr )
{
  outputManager_->outputHomotopy (comm_, paramNames, paramVals, solnVecPtr);
}

const Xyce::NodeNamePairMap &
OutputMgrAdapter::getAllNodes() const
{
  return outputManager_->getAllNodes();
}

} // namespace Analysis
} // namespace Xyce
