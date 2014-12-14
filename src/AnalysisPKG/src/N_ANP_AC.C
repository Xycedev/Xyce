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
//-------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_ANP_AC.C,v $
// Purpose       : AC analysis functions.
// Special Notes :
// Creator       : Ting Mei
// Creation Date :  7/11
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.67.2.1 $
// Revision Date  : $Date: 2014/08/26 22:31:06 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iomanip>

#include <N_ANP_AC.h>

#include <N_ANP_AnalysisManager.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_System.h>
#include <N_LOA_Loader.h>
#include <N_MPDE_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Timer.h>

#include <N_PDS_ParMap.h>
#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_ParComm.h>
#include <mpi.h>
#else
#include <N_PDS_SerialComm.h>
#endif

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>

#include<N_UTL_ExpressionData.h>
#include<N_NLS_ReturnCodes.h>

namespace Xyce {
namespace Analysis {
    
//-----------------------------------------------------------------------------
// Function      : AC::AC( AnalysisManager * )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
AC::AC( AnalysisManager &analysis_manager )
  : AnalysisBase(analysis_manager),
    Util::ListenerAutoSubscribe<StepEvent>(&analysis_manager),
  bVecRealPtr(0),
  bVecImagPtr(0),
  acLoopSize_(0),
  stepMult_(0.0),
  fstep_(0.0),
  currentFreq_(0.0)
{
  bVecRealPtr = linearSystem_.builder().createVector();
  bVecImagPtr = linearSystem_.builder().createVector();

  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);

}

//-----------------------------------------------------------------------------
// Function      : AC::~AC()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
AC::~AC()
{

  if (bVecRealPtr) { delete bVecRealPtr; bVecRealPtr=0; }

  if (bVecImagPtr) { delete bVecImagPtr; bVecImagPtr=0; }

}

void AC::notify(const StepEvent &event) 
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();

    bVecRealPtr->putScalar(0.0);
    bVecImagPtr->putScalar(0.0);
  }
}


//-----------------------------------------------------------------------------
// Function      : AC::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .AC statement.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 6/11
//-----------------------------------------------------------------------------
bool AC::setAnalysisParams(const N_UTL_OptionBlock & paramsBlock)
{
  Util::ParameterList::const_iterator it_tp;
  Util::ParameterList::const_iterator first = paramsBlock.getParams().begin();
  Util::ParameterList::const_iterator last = paramsBlock.getParams().end();
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag()      == "TYPE")
    {
      tiaParams_.type = it_tp->stringValue();
//      tiaParams_.tStartGiven = true;
    }
    else if (it_tp->uTag() == "NP")
    {
      tiaParams_.np = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "FSTART")
    {
      tiaParams_.fStart = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "FSTOP")
    {
      tiaParams_.fStop = it_tp->getImmutableValue<double>();
    }
  }

//  tiaParams_.debugLevel = 1;
  if (DEBUG_ANALYSIS && tiaParams_.debugLevel > 0)
  {
    dout() << section_divider << std::endl
           << "AC simulation parameters" 
           //<< Util::push << std::endl
           << std::endl
           << "number of points  = " << tiaParams_.np << std::endl
           << "starting frequency = " << tiaParams_.fStart << std::endl
           << "stop frequency = " << tiaParams_.fStop 
           //<< Util::pop << std::endl;
           << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::getDCOPFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
bool AC::getDCOPFlag()
{
  return ((getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : AC::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::run()
{
  bool bsuccess = true;

  bsuccess = bsuccess & init();
  bsuccess = bsuccess & loopProcess();

  // if processing the loop failed,
  // then skip finish step
  if( bsuccess )
  {
    bsuccess = bsuccess & finish();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AC::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::init()
{
  bool bsuccess = true;

  acLoopSize_ = setupSweepParam_();

  analysisManager_.setCurrentMode(Analysis::CURRENT_MODE_TRANOP);

  // Get set to do the operating point.
  integrationMethod_ = TIAMethod_NONE;
  workingIntgMethod_->createTimeIntegMethod(integrationMethod_);

  stepNumber            = 0;
  doubleDCOPFlag_ = loader_.getDoubleDCOPFlag ();
  doubleDCOPStep_ = tiaParams_.firstDCOPStep;

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::AC_IC));

  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loader_.setInitialGuess (analysisManager_.getTIADataStore()->nextSolutionPtr);

  // If available, set initial solution (.IC, .NODESET, etc).
  inputOPFlag_ = outputManagerAdapter_.setupInitialConditions( *(analysisManager_.getTIADataStore()->nextSolutionPtr), *(analysisManager_.getTIADataStore()->flagSolutionPtr));

  // Set a constant history for operating point calculation
  analysisManager_.getTIADataStore()->setConstantHistory();
  analysisManager_.getTIADataStore()->computeDividedDifferences();
  workingIntgMethod_->obtainCorrectorDeriv();

  // solving for DC op
  handlePredictor();
  loader_.updateSources();
  stepErrorControl_->newtonConvergenceStatus = nonlinearSolverManager_.solve();
  analysisManager_.getTIADataStore()->stepLinearCombo ();
  gatherStepStatistics_ ();
  stepErrorControl_->evaluateStepError ();

  if ( stepErrorControl_->newtonConvergenceStatus <= 0)
  {
    Report::UserError() << "Solving for DC operating point failed! Cannot continue AC analysis";
    return false;
  }

  // only output DC op if the .op was specified in the netlist
  // or if this AC analysis is called from a .step loop
  if ( analysisManager_.getStepFlag() || analysisManager_.getDotOpFlag() )
  {
    outputManagerAdapter_.dcOutput(
        stepNumber, 
        *analysisManager_.getTIADataStore()->nextSolutionPtr, 
        *analysisManager_.getTIADataStore()->nextStatePtr, 
        *analysisManager_.getTIADataStore()->nextStorePtr,
        objectiveVec_, 
          dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
    outputManagerAdapter_.finishOutput();
  }

  // This for saving the data from the DC op.  different from the above where we are
  // concerned with generating normal output.
  outputManagerAdapter_.outputDCOP(*analysisManager_.getTIADataStore()->nextSolutionPtr);

  loader_.loadBVectorsforAC (bVecRealPtr, bVecImagPtr);

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::AC_IC));

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::loopProcess()
// Purpose       : Conduct the stepping loop.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool AC::loopProcess()
{
  analysisManager_.getTIADataStore()->daeQVectorPtr->putScalar(0.0);
  analysisManager_.getTIADataStore()->daeFVectorPtr->putScalar(0.0);

  analysisManager_.getTIADataStore()->dFdxdVpVectorPtr->putScalar(0.0);
  analysisManager_.getTIADataStore()->dQdxdVpVectorPtr->putScalar(0.0);
  analysisManager_.getTIADataStore()->dQdxMatrixPtr->put(0.0);
  analysisManager_.getTIADataStore()->dFdxMatrixPtr->put(0.0);

//  integrationMethod_ = 6;

  loader_.updateState
                ((analysisManager_.getTIADataStore()->nextSolutionPtr),
               (analysisManager_.getTIADataStore()->currSolutionPtr),
               (analysisManager_.getTIADataStore()->lastSolutionPtr),
               (analysisManager_.getTIADataStore()->nextStatePtr),
               (analysisManager_.getTIADataStore()->currStatePtr),
               (analysisManager_.getTIADataStore()->lastStatePtr),
               (analysisManager_.getTIADataStore()->nextStorePtr),
               (analysisManager_.getTIADataStore()->currStorePtr),
               (analysisManager_.getTIADataStore()->lastStorePtr)
               );

  loader_.loadDAEVectors
              ((analysisManager_.getTIADataStore()->nextSolutionPtr),
               (analysisManager_.getTIADataStore()->currSolutionPtr),
               (analysisManager_.getTIADataStore()->lastSolutionPtr),
               (analysisManager_.getTIADataStore()->nextStatePtr),
               (analysisManager_.getTIADataStore()->currStatePtr),
               (analysisManager_.getTIADataStore()->lastStatePtr),
               (analysisManager_.getTIADataStore()->nextStateDerivPtr),
               (analysisManager_.getTIADataStore()->nextStorePtr),
               (analysisManager_.getTIADataStore()->currStorePtr),
               (analysisManager_.getTIADataStore()->lastStorePtr),
               (analysisManager_.getTIADataStore()->nextStoreLeadCurrQPtr),
               (analysisManager_.getTIADataStore()->daeQVectorPtr),
               (analysisManager_.getTIADataStore()->daeFVectorPtr),
               (analysisManager_.getTIADataStore()->daeBVectorPtr),
               (analysisManager_.getTIADataStore()->dFdxdVpVectorPtr),
               (analysisManager_.getTIADataStore()->dQdxdVpVectorPtr) );

  loader_.loadDAEMatrices(analysisManager_.getTIADataStore()->nextSolutionPtr,
      analysisManager_.getTIADataStore()->nextStatePtr, analysisManager_.getTIADataStore()->nextStateDerivPtr,
      analysisManager_.getTIADataStore()->nextStorePtr,
      analysisManager_.getTIADataStore()->dQdxMatrixPtr,  analysisManager_.getTIADataStore()->dFdxMatrixPtr);

  CPtr_ = rcp(analysisManager_.getTIADataStore()->dQdxMatrixPtr, false);
  GPtr_ = rcp(analysisManager_.getTIADataStore()->dFdxMatrixPtr, false);

  if (tiaParams_.debugLevel > 0)
  {
    Xyce::dout() << "dQdxMatrixPtr:" << std::endl;
    analysisManager_.getTIADataStore()->dQdxMatrixPtr->printPetraObject( Xyce::dout() );

    Xyce::dout() << "dFdxMatrixPtr:" << std::endl;
    analysisManager_.getTIADataStore()->dFdxMatrixPtr->printPetraObject( Xyce::dout() );

    Xyce::dout() << std::endl;
  }

  createLinearSystem_();

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::AC));

  for (int currentStep = 0; currentStep < acLoopSize_; ++currentStep)
  {
    updateCurrentFreq_(currentStep);

    updateLinearSystemFreq_();

    static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::AC, currentFreq_, currentStep));

    bool stepAttemptStatus = solveLinearSystem_();

    if (stepAttemptStatus)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::AC, currentFreq_, currentStep));
      processSuccessfulStep();
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::AC, currentFreq_, currentStep));
      processFailedStep();
    }
  }

  static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish(AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::AC));

  return true;
}


//-----------------------------------------------------------------------------
// Function      : AC::createLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool AC::createLinearSystem_()
{
  bool bsuccess = true;

  RCP<N_PDS_Manager> pdsMgrPtr_;
  pdsMgrPtr_ = rcp(linearSystem_.getPDSManager(), false);

  RCP<N_PDS_ParMap> baseMap;
  RCP<Epetra_CrsGraph> BaseFullGraph_;

  baseMap = rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false);
  BaseFullGraph_ = rcp( pdsMgrPtr_->getMatrixGraph("JACOBIAN"), false );

  int numBlocks = 2;

  RCP<N_PDS_ParMap> blockMap = createBlockParMap(numBlocks, *baseMap);

  // Create a block vector
  BPtr_ = rcp ( new N_LAS_BlockVector ( numBlocks, blockMap, baseMap ) );

  // -----------------------------------------------------
  // Now test block graphs.
  // -----------------------------------------------------

  std::vector<std::vector<int> > blockPattern(2);
  blockPattern[0].resize(2);
  blockPattern[0][0] = 0; blockPattern[0][1] = 1;
  blockPattern[1].resize(2);
  blockPattern[1][0] = 0; blockPattern[1][1] = 1;

  int offset = generateOffset( *baseMap );

  RCP<Epetra_CrsGraph> blockGraph = createBlockGraph( offset, blockPattern, *blockMap, *BaseFullGraph_);

  ACMatrixPtr_ = rcp ( new N_LAS_BlockMatrix( numBlocks, offset, blockPattern, *blockGraph, *BaseFullGraph_) );

  // First diagonal block
  ACMatrixPtr_->put( 0.0 ); // Zero out whole matrix
  ACMatrixPtr_->block( 0, 0 ).add(*GPtr_);
  // Second diagonal block
  ACMatrixPtr_->block( 1, 1 ).add(*GPtr_);

  BPtr_->putScalar( 0.0 );
  BPtr_->block( 0 ).addVec( 1.0, *bVecRealPtr);
  BPtr_->block( 1 ).addVec( 1.0, *bVecImagPtr);

  Amesos amesosFactory;

  XPtr_ = rcp ( new N_LAS_BlockVector (numBlocks, blockMap, baseMap) );
  XPtr_->putScalar( 0.0 );

  blockProblem = rcp(new Epetra_LinearProblem(&ACMatrixPtr_->epetraObj(), &XPtr_->epetraObj(), &BPtr_->epetraObj() ) );

  blockSolver = rcp( amesosFactory.Create( "Klu", *blockProblem ) );

  // Need to reindex the linear system because of the noncontiguous block map.
  Teuchos::ParameterList params;
  params.set( "Reindex", true );
  blockSolver->SetParameters( params );

  // Call symbolic factorization without syncronizing the values, since they are not necessary here. 
  int linearStatus = blockSolver->SymbolicFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : AC::updateLinearSystemFreq_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------

bool AC::updateLinearSystemFreq_()
{
  double omega;

  omega =  2.0 * M_PI * currentFreq_;

  ACMatrixPtr_->block( 0, 1).put( 0.0);
  ACMatrixPtr_->block( 0, 1).add(*CPtr_);
  ACMatrixPtr_->block( 0, 1).scale(-omega);

  ACMatrixPtr_->block(1, 0).put( 0.0);
  ACMatrixPtr_->block(1, 0).add(*CPtr_);
  ACMatrixPtr_->block(1, 0).scale(omega);

  // Copy the values loaded into the blocks into the global matrix for the solve.
  ACMatrixPtr_->assembleGlobalMatrix();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::solveLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :  Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/2011
//-----------------------------------------------------------------------------

bool AC::solveLinearSystem_()
{

  bool bsuccess = true;
 // Solve the block problem

  int linearStatus = blockSolver->NumericFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  linearStatus = blockSolver->Solve();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
    bsuccess = false;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AC::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::processSuccessfulStep()
{
  bool bsuccess = true;

// Output x.
  outputManagerAdapter_.outputAC (currentFreq_, XPtr_->block(0), XPtr_-> block(1));

  if ( !firstDoubleDCOPStep_() )
  {
      stepNumber += 1;
      stats_.successStepsThisParameter_ += 1;
      stats_.successfulStepsTaken_ += 1;
  }

  // This output call is for device-specific output (like from a PDE device,
  // outputting mesh-based tecplot files).  It will only work in parallel if on
  // a machine where all processors have I/O capability.
  loader_.output();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::processFailedStep()
{
  bool bsuccess = true;

  stepNumber += 1;
  acSweepFailures_.push_back(stepNumber);
  stats_.failedStepsAttempted_  += 1;
  stepErrorControl_->numberSuccessiveFailures += 1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AC::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::finish()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  Xyce::dout() << "Calling AC::finish() outputs!" << std::endl;
#endif
  outputManagerAdapter_.finishOutput();
  if (!(acSweepFailures_.empty()))
  {
    bsuccess = false;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : AC::handlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool AC::handlePredictor()
{
  analysisManager_.getTIADataStore()->setErrorWtVector();
  workingIntgMethod_->obtainPredictor();
  workingIntgMethod_->obtainPredictorDeriv();

  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  loader_.startTimeStep ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::updateCurrentFreq_(int stepNumber )
// Purpose       :
// Special Notes : Used for AC analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
bool AC::updateCurrentFreq_(int stepNumber)
{
  double fstart;
  fstart = tiaParams_.fStart;

//  tiaParams_.debugLevel = 1;

  if (tiaParams_.type=="LIN")
  {

    currentFreq_  = fstart + static_cast<double>(stepNumber)*fstep_;
  }
  else if(tiaParams_.type=="DEC" || tiaParams_.type=="OCT")
  {

    currentFreq_ = fstart*pow(stepMult_, static_cast<double>(stepNumber) );
  }
  else
  {
    Report::DevelFatal().in("AC::updateCurrentFreq_") << "AC::updateCurrentFreq_: unsupported STEP type";
  }

  if (analysisManager_.getDebugLevel() > 0)
  {
    Xyce::dout() << "currentFreq_  = " << currentFreq_ << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AC::setupSweepParam_
// Purpose       : Processes sweep parameters.
// Special Notes : Used for AC analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
int AC::setupSweepParam_()
{
  double fstart, fstop;
  double fcount = 0.0;

  fstart = tiaParams_.fStart;
  fstop = tiaParams_.fStop;

#ifdef Xyce_DEBUG_ANALYSIS
  if (analysisManager_.getDebugLevel() > 0)
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "AC::setupSweepParam_" << std::endl;
  }
#endif

    if (tiaParams_.type=="LIN")
    {
      int np = static_cast<int> (tiaParams_.np);

      if ( np == 1)
        fstep_ = 0;
      else
        fstep_  = (fstop - fstart)/(tiaParams_.np - 1);

      fcount = tiaParams_.np;
#ifdef Xyce_DEBUG_ANALYSIS
      if (analysisManager_.getDebugLevel() > 0)
      {
        Xyce::dout() << "fstep   = " << fstep_  << std::endl;
      }
#endif
    }
    else if(tiaParams_.type=="DEC")
    {
      stepMult_ = pow(10,(1/tiaParams_.np));
      fcount   = floor(fabs(log10(fstart) - log10(fstop)) * tiaParams_.np + 1);
#ifdef Xyce_DEBUG_ANALYSIS
      if (analysisManager_.getDebugLevel() > 0)
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
#endif
    }
    else if(tiaParams_.type=="OCT")
    {
      stepMult_ = pow(2,1/(tiaParams_.np));

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2=log(2.0);
      fcount   = floor(fabs(log(fstart) - log(fstop))/ln2 * tiaParams_.np + 1);
#ifdef Xyce_DEBUG_ANALYSIS
      if (analysisManager_.getDebugLevel() > 0)
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
#endif
    }
    else
    {
      Report::DevelFatal().in("AC::setupSweepParam") << "Unsupported type";
    }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return static_cast<int> (fcount);
}

} // namespace Analysis
} // namespace Xyce
