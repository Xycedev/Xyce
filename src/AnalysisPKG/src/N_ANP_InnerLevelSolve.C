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
// Filename      : $RCSfile: N_ANP_InnerLevelSolve.C,v $
// Purpose       : This file contains functions from control-algorithm that are
//                 used by the inner solve of a 2-level xyce-to-xyce solve.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 1/23/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.65.2.2 $
// Revision Date  : $Date: 2014/08/26 22:31:06 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#define Xyce_MPDE_IC
#include <N_UTL_Misc.h>

#include <iostream>
#include <ctime>

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_Transient.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_ActiveOutput.h>
#include <N_IO_CmdParse.h>
#include <N_LOA_Loader.h>
#include <N_NLS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_TwoLevelError.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::provisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/04/2009
//-----------------------------------------------------------------------------
bool AnalysisManager::provisionalStep (double maxTimeStep,  double &timeStep)
{
  bool bsuccess = true;
  bool b1 = true;
  std::string msg;
  bool dcopFlag = true;

  if (!initializeSolvers_mixedSignal_)
  {
    if (analysisMode_ == ANP_MODE_TRANSIENT)
    {
      mixedSignalAnalysisObject_ = Teuchos::rcp(new N_ANP_Transient(*this));
      mixedSignalAnalysisObject_->setAnalysisParams(tranParamsBlock);
      stepErrorControl_->resetAll();
    }
    else if (analysisMode_ == ANP_MODE_DC_SWEEP)
    {
      mixedSignalAnalysisObject_ = Teuchos::rcp(new N_ANP_DCSweep(*this));
      for (int i=0;i<dcParamsBlockVec.size();++i)
      {
        mixedSignalAnalysisObject_->setAnalysisParams(dcParamsBlockVec[i]);
      }
    }
    else
    {
      Report::DevelFatal().in("AnalysisManager::provisionalStep") << "unknown type of analysis";
    }
    mixedSignalAnalysisObject_->init();

    // Start the solvers timer.
    xyceTranTimerPtr_->resetStartTime();

    initializeSolvers_mixedSignal_=true;
  }

  primaryAnalysisObject_ = mixedSignalAnalysisObject_;

  RefCountPtr<N_ANP_Transient> mixedSignalTransientAnalysisObject =
    Teuchos::rcp_dynamic_cast<N_ANP_Transient>(mixedSignalAnalysisObject_);

  if (!Teuchos::is_null(mixedSignalTransientAnalysisObject))
  {
    dcopFlag = mixedSignalTransientAnalysisObject->getDCOPFlag();
  }

  // Now save time step info, in case this step gets rejected.

  if (!(stepErrorControl_->isFinished()))
  {
    bool stepSuccess = false;

    if (dcopFlag) // if dcop step, make one attempt.
    {
      mixedSignalAnalysisObject_->preStepDetails (maxTimeStep);
      b1 = mixedSignalAnalysisObject_->mixedSignalStep();

      // only call finalize step here if we have failed.
      if (!stepErrorControl_->stepAttemptStatus)
      {
        mixedSignalAnalysisObject_->finalizeStep();
      }
      stepSuccess = stepErrorControl_->stepAttemptStatus;
    }
    else // else, if transient step, keep re-taking the step
         // until it succeeds, or gets an unrecoverable failure,
         // such as time-step-too-small.
    {
      bool recoverableFailureFlag=true;
      while (!stepSuccess && recoverableFailureFlag)
      {
        mixedSignalAnalysisObject_->preStepDetails (maxTimeStep);
        b1 = mixedSignalAnalysisObject_->mixedSignalStep();

        // Only call finalize step here if step has failed.
        // If we succeed, we want to give Habanero the opportunity
        // to reject the step, after this function (provisionalStep)
        // exits.
        if (!stepErrorControl_->stepAttemptStatus)
        {
          recoverableFailureFlag = mixedSignalAnalysisObject_->finalizeStep ();
        }
        else
        {
          stepSuccess = true;
        }
      }
    }
    bsuccess = stepSuccess;
  }

  // get the step information.
  //

  if (dcopFlag)
  {
    timeStep = 0.0;
  }
  else
  {
    N_TIA_TimeIntInfo tiInfo;
    getTimeIntInfo(tiInfo);
    timeStep = tiInfo.nextTimeStep;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::acceptProvisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//-----------------------------------------------------------------------------
void AnalysisManager::acceptProvisionalStep ()
{
  mixedSignalAnalysisObject_->finalizeStep ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::rejectProvisionalStep
// Purpose       : Used by mixed-signal.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/2009
//-----------------------------------------------------------------------------
void AnalysisManager::rejectProvisionalStep ()
{
  stepErrorControl_->stepAttemptStatus = false;
  stepErrorControl_->updateBreakPoints();

  bool dcopFlag = false;
  RefCountPtr<N_ANP_Transient> mixedSignalTransientAnalysisObject =
    Teuchos::rcp_dynamic_cast<N_ANP_Transient>(mixedSignalAnalysisObject_);

  if (!Teuchos::is_null(mixedSignalTransientAnalysisObject))
  {
    dcopFlag = mixedSignalTransientAnalysisObject->getDCOPFlag();
  }

  if (dcopFlag)
  {
    mixedSignalAnalysisObject_->finalizeStep ();
  }
  // Transient
  else
  {
     //   bool b1 = processFailedStep();
    loader_->stepFailure(currentMode_);
    workingIntgMethod_->rejectStepForHabanero();

    mixedSignalTransientAnalysisObject->stats_.failedStepsAttempted_  += 1;
    stepErrorControl_->numberSuccessiveFailures += 1;

  } // transient

#if  0
  if (stepErrorControl_->isPauseTime())
  {
    // Failure at this point only indicates that the simulation
    // is paused and may be resumed.
    stepErrorControl_->simulationPaused();
    isPaused = true;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setExternalSolverState 
// Purpose       : 
// Special Notes : Used for multi-level Newton solves, for levels other
//                 than the top level.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/22/2014
//-----------------------------------------------------------------------------
void AnalysisManager::setExternalSolverState (const N_DEV_SolverState & ss)
{
  loader_->setExternalSolverState (ss);
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::runStep
//
// Purpose       : This function is similar to "run" except that only a single
//                 integration (or DC sweep) step will be executed.
//
// Special Notes : Used for multi-level Newton solves, for levels other
//                 than the top level.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::runStep
    (const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError)
{
  std::string msg;

  if (initializeAllFlag_ == false)
  {
    Report::DevelFatal().in("AnalysisManager::runStep") << "Call initializeAll function first";
  }

  if (analysisParamsRegistered == false)
  {
    Report::UserError() << "No analysis statement in the netlist";
    return false;
  }

  bool integration_status = false;

#ifdef Xyce_VERBOSE_TIME
  tiaParams_.printParams(Xyce::lout(), anpAnalysisModeToNLS(analysisMode_));
#endif

  solverStartTime_ = elapsedTimerPtr_->elapsedTime();

  if (stepLoopFlag_)
  {
    Report::UserError() << "Not valid to use .STEP statements in an inner solve";
  }
  else
  {
    if (analysisMode_ == ANP_MODE_TRANSIENT)
    {
#ifdef Xyce_VERBOSE_TIME
      Xyce::dout() << "AnalysisManager::runStep:" << std::endl;
      Xyce::dout() << "nextTime = " << tiInfo.nextTime << std::endl;
      Xyce::dout() << "stepSize = " << tiInfo.nextTimeStep << std::endl;
#endif
    }
    else if (analysisMode_ == ANP_MODE_DC_SWEEP)
    {
      // do nothing
    }
    else
    {
      Report::DevelFatal().in("AnalysisManager::runStep") << "Unknown type of analysis";
    }
    integration_status  = twoLevelAnalysisObject_->twoLevelStep();
    workingIntgMethod_->setupTwoLevelError(tlError);
  }

  return integration_status;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::startTimeStep_
// Purpose       : used by 2-level solves.
// Special Notes : One of the primary purposes for this function is to impose
//                 a lot of upper level information from the top level time
//                 integrator on the inner level time integrator.  This
//                 information is contained in the N_TIA_TimeIntInfo
//                 class (tiInfo here).  This includes things like the
//                 time step size, the time integration order, etc.
//
//                 In general, in a 2-level solve, an inner solver doesn't
//                 have enough information to correctly determine the
//                 step size, order, etc.  This is in part because the
//                 inner solver, while it knows about the top level solver,
//                 it cannot know about any OTHER inner solvers.
//
//                 The top level solver, however, does have enough information.
//                 It gathers break point, error analysis, and other info
//                 from all the inner solves.  So, the top level solver
//                 makes all the decisions and imposes them on the inner
//                 solves.  This function is where it does that impose.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::startTimeStep (const N_TIA_TimeIntInfo & tiInfo)
{
  // Beginning Integration flag is the only piece of data in tiInfo that is
  // currently owned by the control algorithm class.  Everything else
  // is owned by the step error control, or the integration method.
  twoLevelAnalysisObject_->setBeginningIntegrationFlag(tiInfo.beginIntegrationFlag);

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams_.debugLevel > 0)
  {
    Xyce::dout() << "AnalysisManager::startTimeStep:" << std::endl;
  }
#endif

#ifdef Xyce_VERBOSE_TIME
  if( !is_null(twoLevelAnalysisObject_) )
  {
    // an MPDE run also traverses this area so catch case where this is null
    twoLevelAnalysisObject_->printStepHeader(Xyce::lout());
  }
#endif
  if (switchIntegrator_)
    workingIntgMethod_->createTimeIntegMethod( twoLevelAnalysisObject_->getIntegrationMethod() );

  // ------------------------------------------------------------------------
  // Set the step size, current time and next time.
#ifndef Xyce_CHARON
  if ( twoLevelAnalysisObject_->getIntegrationMethod() != TIAMethod_NONE)
#endif
  {
    stepErrorControl_->updateTwoLevelTimeInfo(tiInfo);
  }


  if (twoLevelAnalysisObject_->getBeginningIntegrationFlag() && stepErrorControl_->stepAttemptStatus)
  {
    workingIntgMethod_->setTwoLevelTimeInfo(tiInfo); // new-dae only
  }

  // ------------------------------------------------------------------------
  // If we've switched the integration method, we need to obtain the
  // corrector derivative only after we've updated the TimeInfo.
  if (switchIntegrator_)
  {
    switchIntegrator_ = false;
    workingIntgMethod_->obtainCorrectorDeriv();
  }

  bool dcopFlag = true;
  RefCountPtr<N_ANP_Transient> twoLevelTransientAnalysisObject = Teuchos::rcp_dynamic_cast<N_ANP_Transient>(twoLevelAnalysisObject_);
  if (!Teuchos::is_null(twoLevelTransientAnalysisObject)) {
    dcopFlag = twoLevelTransientAnalysisObject->getDCOPFlag();
  }
#ifdef Xyce_VERBOSE_TIME
  if (!dcopFlag) stepErrorControl_->outputTimeInfo(lout());
#endif

  // ------------------------------------------------------------------------
  // Set the nonlinear solver parameters to those appropriate for the
  // transient solution, if neccessary.
  if (!dcopFlag)
  {
    nlsMgrPtr_->setAnalysisMode(anpAnalysisModeToNLS(ANP_MODE_TRANSIENT));
  }

  // Ask the method to update its coefficients
  workingIntgMethod_->updateCoeffs(); 
  twoLevelAnalysisObject_->handlePredictor();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::conductanceTest
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
void AnalysisManager::conductanceTest ()
{
  std::map<std::string,double> inputMap;
  std::vector<double> outputVector;
  std::vector< std::vector<double> > jacobian;

  // load inputMap from tiaParam.condTestDeviceNames option
  std::list< std::string >::iterator currentDeviceName = tiaParams_.condTestDeviceNames.begin();
  std::list< std::string >::iterator endDeviceName = tiaParams_.condTestDeviceNames.end();
  while( currentDeviceName != endDeviceName )
  {
    inputMap[ *currentDeviceName ] = 0.0;
    ++currentDeviceName;
  }

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams_.debugLevel > 0)
  {
    Xyce::dout() << "AnalysisManager::conductanceTest()" << std::endl;
    currentDeviceName = tiaParams_.condTestDeviceNames.begin();
    while( currentDeviceName != endDeviceName )
    {
      Xyce::dout() << "currentDeviceName = \"" << *currentDeviceName << "\" added to inputMap[ "
            << *currentDeviceName << " ] = " << inputMap[ *currentDeviceName ] << std::endl;
      ++currentDeviceName;
    }
  }
#endif

  int isize=inputMap.size();
  outputVector.resize(isize,0.0);
  jacobian.resize(isize);
  for (int i=0;i<isize;++i)
  {
    jacobian[i].resize(isize,0.0);
  }

  bool b1 = nlsMgrPtr_->obtainConductances(
        inputMap,
        outputVector,
        jacobian
    );

  int iE1, iE2;
  int numElectrodes = isize;

  FILE *fp1;
  fp1 = fopen("conductance.txt","w");

  fprintf(fp1, "%s", "Conductance array: \n");
  fprintf(fp1,"%s", "              ");
  if (b1)
  {
    std::map<std::string,double>::iterator iterM = inputMap.begin();
    std::map<std::string,double>::iterator  endM = inputMap.end  ();
    for (iE2 = 0; iE2 < numElectrodes; ++iE2,++iterM)
    {
      std::string srcname = iterM->first;
      fprintf(fp1,"\t%14s",srcname.c_str());
    }
    fprintf(fp1,"%s", "\n");

    iterM = inputMap.begin();
    for (iE1 = 0; iE1 < numElectrodes; ++iE1, ++iterM)
    {
      std::string srcname = iterM->first;
      fprintf(fp1,"%14s",srcname.c_str());
      for (iE2 = 0; iE2 < numElectrodes; ++iE2)
      {
        fprintf(fp1,"\t%14.4e",jacobian[iE1][iE2]);
      }
      fprintf(fp1,"%s", "\n");
    }
    fprintf(fp1,"%s", "\n");
  }
  else
  {
    fprintf(fp1,"%s", "\nConductance calculation failed!\n");
  }

  fclose(fp1);

}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::startupSolvers
// Purpose       :
// Special Notes : Used only for 2-level solves.
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 3/10/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::startupSolvers()
{
  bool bsuccess = true;
  std::string msg;

  if (analysisMode_ == ANP_MODE_TRANSIENT)
  {
    nlsMgrPtr_->resetAll(Nonlinear::DC_OP);
    twoLevelAnalysisObject_ = Teuchos::rcp(new N_ANP_Transient(*this));
    twoLevelAnalysisObject_->setAnalysisParams(tranParamsBlock);
    stepErrorControl_->resetAll();
  }
  else if (analysisMode_ == ANP_MODE_DC_SWEEP)
  {
    stepErrorControl_->setTIAParams();
    twoLevelAnalysisObject_ = Teuchos::rcp(new N_ANP_DCSweep(*this));
    for (int i=0;i<dcParamsBlockVec.size();++i)
    {
      twoLevelAnalysisObject_->setAnalysisParams(dcParamsBlockVec[i]);
    }
  }
  else
  {
    Report::UserError() << "Multi-Level Newton solves only supports DC and Transient analysis";
    return false;
  }
  
  primaryAnalysisObject_ = twoLevelAnalysisObject_;
  twoLevelAnalysisObject_->init();

  // outputManagerAdapter_->prepareOutput(analysisMode_);
  activeOutput_ = new IO::ActiveOutput(*outMgrPtr_);
  activeOutput_->add(pdsMgrPtr_->getPDSComm()->comm(), analysisMode_);

  // Reset the solvers timer.
  xyceTranTimerPtr_->resetStartTime();

  // Hardwire the erroption parameter to 1, which will force the inner
  // solve (initiated by this function) to only use the Newton success/failure
  // as step criteria.  Predictor-corrector information is handled in the
  // upper level solver.
  tiaParams_.errorAnalysisOption = 1;
  // tscoffe 03/07/08  With the new errorAnalysisOption = 1 features of nlmin
  // and nlmax, we should enforce the original mode where nlmin=nlmax=maxNLSteps.
  //int maxNLSteps = ???;
  //tiaParams_.NLmax = maxNLSteps;
  //tiaParams_.NLmin = maxNLSteps;

 return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::finishSolvers
// Purpose       :
// Special Notes : Used only for 2-level solves.
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/10/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::finishSolvers ()
{
  bool bsuccess = true;
  std::string msg;

  twoLevelAnalysisObject_->finish();

  delete activeOutput_;
  activeOutput_ = 0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::homotopyStepSuccess
// Purpose       : Lower-level processing of a successful homotopy step,
//                 which was controlled from the upper level of a 2-level solve.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void AnalysisManager::homotopyStepSuccess
    ( const std::vector<std::string> & paramNames,
      const std::vector<double> & paramVals)
{
#ifdef Xyce_DEBUG_ANALYSIS
  std::string netListFile = commandLine_.getArgumentValue("netlist");
  Xyce::dout() << "\n " << netListFile;
  Xyce::dout() << " AnalysisManager::homotopyStepSuccess " << std::endl;
#endif
  // output:
  outputManagerAdapter_->outputHomotopy( paramNames, paramVals, *getTIADataStore()->nextSolutionPtr );

  // update the data arrays:
  getTIADataStore()->updateSolDataArrays();

  // pass info to the next level down, if it exists.
  loader_->homotopyStepSuccess (paramNames,paramVals);

  return;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::homotopyStepFailure
//
// Purpose       : Lower-level processing of a failed homotopy step,
//                 which was controlled from the upper level of a
//                 2-level solve.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void AnalysisManager::homotopyStepFailure ()
{
#ifdef Xyce_DEBUG_ANALYSIS
  std::string netListFile = commandLine_.getArgumentValue("netlist");
  Xyce::dout() << "\n " << netListFile;
  Xyce::dout() << " AnalysisManager::homotopyStepFailure " << std::endl;
#endif

  // The solutions currently in place represent failure.  Get rid of them.
  getTIADataStore()->usePreviousSolAsPredictor ();

  // pass info to the next level down, if it exists.
  loader_->homotopyStepFailure ();

  return;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void AnalysisManager::stepSuccess(CurrentMode analysisUpper)
{

#ifdef Xyce_DEBUG_ANALYSIS
  std::string netListFile = commandLine_.getArgumentValue("netlist");
  Xyce::dout() << "\n " << netListFile;
  Xyce::dout() << " AnalysisManager::stepSuccess " << std::endl;
#endif

  setCurrentMode(analysisUpper);
  stepErrorControl_->stepAttemptStatus = true;
  switch (analysisUpper)
  {
    case Analysis::CURRENT_MODE_TRANOP:
      {
        N_ANP_Transient * twoLevelTransientAnalysisObject = dynamic_cast<N_ANP_Transient*>(&*twoLevelAnalysisObject_);
        if (twoLevelTransientAnalysisObject)
        {
          twoLevelTransientAnalysisObject->processSuccessfulDCOP();
        }
        else
        {
          Report::DevelFatal().in("AnalysisManager::stepSuccess") << "Failed dynamic_cast of twoLevelAnalysisObject to N_ANP_Transient.";
        }
      }
      break;
    case Analysis::CURRENT_MODE_TRANSIENT:
      twoLevelAnalysisObject_->processSuccessfulStep();
      break;
    case Analysis::CURRENT_MODE_DC_SWEEP:
      twoLevelAnalysisObject_->processSuccessfulStep();
      break;
    default:
      Report::DevelFatal().in("AnalysisManager::stepSuccess") << "Unknown type of analysis";
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void AnalysisManager::stepFailure(Analysis::CurrentMode analysisUpper)
{
  setCurrentMode(analysisUpper);
  stepErrorControl_->stepAttemptStatus = false;
  switch (analysisUpper)
  {
    case Analysis::CURRENT_MODE_TRANOP:
      {
        N_ANP_Transient * twoLevelTransientAnalysisObject = dynamic_cast<N_ANP_Transient*>(&*twoLevelAnalysisObject_);
        if (twoLevelTransientAnalysisObject)
        {
          twoLevelTransientAnalysisObject->processFailedDCOP();
        }
        else
        {
          Report::DevelFatal().in("AnalysisManager::stepFailure") << "Failed dynamic_cast of twoLevelAnalysisObject to N_ANP_Transient.";
        }
      }
      break;
    case Analysis::CURRENT_MODE_TRANSIENT:
      twoLevelAnalysisObject_->processFailedStep();
      break;
    case Analysis::CURRENT_MODE_DC_SWEEP:
      twoLevelAnalysisObject_->processFailedStep();
      break;
    default:
      Report::DevelFatal().in("AnalysisManager::stepFailure") << "Unknown type of analysis";
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitialQnorm
// Purpose       : Used for 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool AnalysisManager::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  bool bsuccess = true;
  workingIntgMethod_->getInitialQnorm (tle);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBreakPoints
// Purpose       : Used for 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool AnalysisManager::getBreakPoints(std::vector<N_UTL_BreakPoint> &breakPointTimes)
{
  bool bsuccess = true;
  loader_->getBreakPoints(breakPointTimes);
  return bsuccess;
}

} // namespace Analysis
} // namespace Xyce
