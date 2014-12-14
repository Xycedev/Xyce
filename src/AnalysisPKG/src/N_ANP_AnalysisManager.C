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
// Filename      : $RCSfile: N_ANP_AnalysisManager.C,v $
// Purpose       : This file contains the functions which define the time
//                 domain & integration algorithm classes.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.144.2.5 $
// Revision Date  : $Date: 2014/08/28 21:00:43 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <ctime>
#include <iostream>
#include <sstream>

#include <N_ANP_AnalysisManager.h>

#include <N_ANP_AC.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_Dakota.h>
#include <N_ANP_HB.h>
#include <N_ANP_MOR.h>
#include <N_ANP_MPDE.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_Report.h>
#include <N_ANP_Step.h>
#include <N_ANP_Transient.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_ActiveOutput.h>
#include <N_IO_RestartMgr.h>
#include <N_LAS_LAFactory.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LOA_CktLoader.h>
#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_MPDE_Manager.h>
#include <N_NLS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_MPDEInterface.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_TwoLevelError.h>
#include <N_TOP_Topology.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Stats.h>
#include <N_UTL_Timer.h>

namespace Xyce {
namespace Analysis {

namespace {

// .RESULTS
struct AnalysisManager_ResultOptionsReg : public IO::PkgOptionsReg
{
  AnalysisManager_ResultOptionsReg(AnalysisManager &analysis_manager)
    : analysisManager_(analysis_manager)
  {}

  bool operator()(const Util::OptionBlock & option_block)
  {
    analysisManager_.getOutputManagerAdapter().addOutputResults(option_block);

    return true;
  }

private:
  AnalysisManager &     analysisManager_;
};

  // .TRAN
  struct AnalysisManager_TransAnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_TransAnalysisReg( AnalysisManager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTranAnalysisParams( options ); }

    AnalysisManager * Mgr;
  };

  // .options TIMEINT
  struct AnalysisManager_TranOptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_TranOptionsReg( AnalysisManager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setTranOptions( options ); }

    AnalysisManager * Mgr;
  };

  // .DC
  struct AnalysisManager_DCAnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_DCAnalysisReg( AnalysisManager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setDCAnalysisParams( options ); }

    AnalysisManager * Mgr;
  };

  // .options OP_IO
  struct AnalysisManager_DCOPOptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_DCOPOptionsReg( AnalysisManager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setDCOPRestartParams( options ); }

    AnalysisManager * Mgr;
  };


  // .OP
  struct AnalysisManager_OPAnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_OPAnalysisReg( AnalysisManager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setOPAnalysisParams( options ); }

    AnalysisManager * Mgr;
  };

  // .STEP
  struct AnalysisManager_STEPAnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_STEPAnalysisReg( AnalysisManager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setSTEPAnalysisParams( options ); }

    AnalysisManager * Mgr;
  };

  // .options SAVE
  struct AnalysisManager_SaveOptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_SaveOptionsReg( AnalysisManager * mgr )
    : Mgr(mgr)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return Mgr->setSaveOptions( options ); }

    AnalysisManager * Mgr;
  };

  // .MPDE
  struct AnalysisManager_MPDE_AnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_MPDE_AnalysisReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMPDEAnalysisParams( options ); }

    AnalysisManager * anpInt_;
  };

  // .HB
  struct AnalysisManager_HB_AnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_HB_AnalysisReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setHBAnalysisParams( options ); }

    AnalysisManager * anpInt_;
  };

  // .options HBINT
  struct AnalysisManager_HB_OptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_HB_OptionsReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setHBOptions( options ); }

    AnalysisManager * anpInt_;
  };

  // .options LINSOL
  struct AnalysisManager_LinSolReg : public IO::PkgOptionsReg
  {
    AnalysisManager_LinSolReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setLinSol( options ); }

    AnalysisManager * anpInt_;
  };

  // .options LINSOL-HB
  struct AnalysisManager_HB_LinSolReg : public IO::PkgOptionsReg
  {
    AnalysisManager_HB_LinSolReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setHBLinSol( options ); }

    AnalysisManager * anpInt_;
  };

  // .options MPDEINT
  struct AnalysisManager_MPDE_OptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_MPDE_OptionsReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMPDEOptions( options ); }

    AnalysisManager * anpInt_;
  };

  // .options TIMEINT-MPDE
  struct AnalysisManager_MPDE_TranMPDEOptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_MPDE_TranMPDEOptionsReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setTRANMPDEOptions( options ); }

    AnalysisManager * anpInt_;
  };

  // .AC
  struct AnalysisManager_AC_AnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_AC_AnalysisReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setACAnalysisParams( options ); }

    AnalysisManager * anpInt_;
  };

  // .MOR
  struct AnalysisManager_MOR_AnalysisReg : public IO::PkgOptionsReg
  {
    AnalysisManager_MOR_AnalysisReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMORAnalysisParams( options ); }

    AnalysisManager * anpInt_;
  };

  // .options MOR_OPTS
  struct AnalysisManager_MOR_OptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_MOR_OptionsReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setMOROptions( options ); }

    AnalysisManager * anpInt_;
  };

  
  // .options SENS 
  struct AnalysisManager_SensOptionsReg : public IO::PkgOptionsReg
  {
    AnalysisManager_SensOptionsReg( AnalysisManager * anpInt )
    : anpInt_(anpInt)
    {}

    bool operator()( const N_UTL_OptionBlock & options )
    { return anpInt_->setSensOptions( options ); }

    AnalysisManager * anpInt_;
  };


} // namespace <unnamed>


const char *
analysisModeName(Analysis_Mode mode) 
{
  static const char * const mode_names[] = {"Invalid", "DC OP", "DC Sweep", "Transient", "MPDE", "HB", "AC", "MOR"};

  if (mode < sizeof(mode_names)/sizeof(mode_names[0]))
    return mode_names[mode];
  else
    return mode_names[0];
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::AnalysisManager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------

AnalysisManager::AnalysisManager(IO::CmdParse & cp, Stats::Stat root_stat)
  : Util::Notifier<StepEvent>(),
  Util::Notifier<AnalysisEvent>(),
  Util::ListenerAutoSubscribe<StepEvent>(this),
  Util::ListenerAutoSubscribe<AnalysisEvent>(this),
  commandLine_(cp),
  tiaParams_(cp),
  workingIntgMethod_(0),
  stepErrorControl_(0),
  linearSystem_(0),
  nlsMgrPtr_(0),
  loader_(0),
  cktLoaderPtr_(0),
  restartPtr_(0),
  pkgOptMgrPtr_(0),
  nonlinearEquationLoaderPtr_(0),
  devInterfacePtr_(0),
  topoMgrPtr_(0),
  outMgrPtr_(0),
  appBuilderPtr_(0),
  pdsMgrPtr_(0),
  outputManagerAdapter_(0),
  tiaDataStore_(0),
  analysisMode_(ANP_MODE_TRANSIENT),
  analysisParamsRegistered(false),
  firstTime(true),
  oldPercentComplete(0.0),
  startSimTime(-1.0),
  calledBeforeTwoLevelTran_(false),
  switchIntegrator_(false),
  initializeAllFlag_(false),
  startTRANtime_(0.0),
  stepLoopFlag_(false),
  stepLoopInitialized_(false),
  dcLoopInitialized_(false),
  gui_(false),
  daeStateDerivFlag_(true),
  initializeSolvers_mixedSignal_(false),
  dcLoopSize_(0),
  sweepSourceResetFlag_(true),
  progressFlag_(true),
  rootStat_("Analysis Manager", root_stat),
  xyceTranTimerPtr_(0),
  elapsedTimerPtr_(0),
  solverStartTime_(0.0),
  dakotaRunFlag_(false),
  dakotaIterationNumber_(0),
  saveTime_(0.0),
  saveTimeGiven_(false),
  saveFlag_(false),
  savedAlready_(false),
  dcopRestartFlag_(false),
  dotOpSpecified_(false),
  initialOutputInterval_(0.0),
  outputIntervals_(),
  nextOutputTime_(0.0),
  initialRestartInterval_(0.0),
  restartIntervals_(),
  nextRestartSaveTime_(0.0),
  blockAnalysisFlag_(false),
  hbFlag_(false),
  mpdeFlag_(false),
  sensFlag_(false),
  analysisObject_(),
  stepAnalysisTarget_(),
  dakotaAnalysisTarget_(),
  primaryAnalysisObject_(),
  mpdeMgrPtr_(0),
  tiaMPDEIfacePtr_(),
  twoLevelAnalysisObject_(),
  mixedSignalAnalysisObject_(),
  breakPointRestartStep(0),
  currentMode_(CURRENT_MODE_TRANOP)
{
  gui_ = commandLine_.argExists("-gui");

  // check for maximum order on the command line
  if ( commandLine_.argExists ("-maxord") )
  {
    tiaParams_.maxOrder =
      atoi( commandLine_.getArgumentValue( "-maxord" ).c_str() );

    if ( tiaParams_.maxOrder < 1) tiaParams_.maxOrder = 1;
    if ( tiaParams_.maxOrder > 5) tiaParams_.maxOrder = 5;
  }

  // Create the MPDEIface object and register it with tiaControl
  tiaMPDEIfacePtr_ = rcp(new N_TIA_MPDEInterface(tiaParams_));
  tiaMPDEIfacePtr_->registerTIAControl(this);

  return;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::~AnalysisManager
//
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
AnalysisManager::~AnalysisManager()
{
  delete cktLoaderPtr_;
  delete nonlinearEquationLoaderPtr_;
  delete outputManagerAdapter_;
  delete mpdeMgrPtr_;
  delete workingIntgMethod_;
  delete tiaDataStore_;
  delete stepErrorControl_;
  delete xyceTranTimerPtr_;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Jul 14 12:15:30 2014
//-----------------------------------------------------------------------------
///
/// Notification that there is a StepEvent.
///
/// @param step_event   information about the event
///
void
AnalysisManager::notify(
  const StepEvent &     step_event)
{
  if (step_event.state_ == StepEvent::STEP_STARTED) {
    tiaDataStore_->setZeroHistory();
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Jul 14 12:15:30 2014
//-----------------------------------------------------------------------------
///
/// Notification that there is a AnalysisEvent.
///
/// @param time_integrator_event   information about the event
///
void
AnalysisManager::notify(
  const AnalysisEvent &     analysis_event)
{
//  Xyce::dout() << Dump<AnalysisEvent>(analysis_event) << std::endl;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::resetAll
//
// Purpose       : just like a destructor without the death
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
void AnalysisManager::resetAll()
{
  // tiaDataStore_ is created in initializeAll
  delete tiaDataStore_;
  tiaDataStore_ = 0;

  // stepErrorControl_ is created in initializeAll
  delete stepErrorControl_;
  stepErrorControl_ = 0;

  // wimPtr is created in initializeAll
  delete workingIntgMethod_;
  workingIntgMethod_ = 0;

  // Reset step statistics to zero.
  primaryAnalysisObject_->resetAll();

  initializeAllFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBlockAnalysisFlag
// Purpose       : "get" function for MPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getBlockAnalysisFlag () const
{
  if (!mpdeMgrPtr_)
  {
    return blockAnalysisFlag_;
  }
  return mpdeMgrPtr_->blockAnalysisFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getMPDEFlag
// Purpose       : "get" function for MPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getMPDEFlag ()
{
  if (!mpdeMgrPtr_)
  {
    return mpdeFlag_;
  }
  return mpdeMgrPtr_->getMPDEFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getMPDEStartupFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getMPDEStartupFlag ()
{
  if (!mpdeMgrPtr_)
  {
    return false;
  }
  return mpdeMgrPtr_->getMPDEStartupFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getMPDEIcFlag
// Purpose       : get function for MPDE initial condition flag (true if MPDE & IC
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getMPDEIcFlag()
{
  if (!mpdeMgrPtr_)
  {
    return false;
  }
  return mpdeMgrPtr_->getMPDEIcFlag();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getWaMPDEFlag ()
// Purpose       : "get" function for WaMPDE flag. (true if not IC)
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/31/08
//-----------------------------------------------------------------------------
bool AnalysisManager::getWaMPDEFlag ()
{
  if (!mpdeMgrPtr_)
  {
    return false;
  }
  return mpdeMgrPtr_->getWaMPDEFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::isPaused
//
// Purpose       : return the true if the simulation is currently paused
//
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical and Microsystems Simulation
// Creation Date : 02/18/2008
//-----------------------------------------------------------------------------
bool AnalysisManager::isPaused()
{
  return stepErrorControl_->isPauseTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::initializeAll
// Purpose       : This function performs the final initializations.  Mainly,
//                 it initializes the two data container classes, DataStore and
//                 StepErrorControl.  It also registers the neccessary vectors
//                 with the LAS system class.
// Special Notes : This function should *only* be called after all the
//                 registrations have all been performed.  In particular, the
//                 N_LAS_System class *must* be registered before this function
//                 is called.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
bool AnalysisManager::initializeAll(N_LOA_Loader * tmpLoaderPtr)
{
  static Stats::Stat initialize_stat("Initialize", rootStat_);

  Stats::TimeBlock x(initialize_stat);

  if (!linearSystem_)
  {
    Report::DevelFatal0().in("AnalysisManager::initializeAll")
      << "Register LAS system first.";
  }

  tiaParams_.solutionSize = linearSystem_->getSolutionSize();
  tiaParams_.stateSize    = linearSystem_->getStateSize();

  // allocate data store class, which will allocate all the vectors.
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams_.debugLevel > 0)
  {
    Xyce::dout() << "AnalysisManager::initializeAll.  " ;
    Xyce::dout() << "  maximum order = " << tiaParams_.maxOrder << std::endl;
  }
#endif

  // Allocate circuit loader: (must be before step error control allocation)
  if (tmpLoaderPtr==0)
  {
    loader_ = cktLoaderPtr_ = new N_LOA_CktLoader();
    cktLoaderPtr_->registerDeviceInterface(devInterfacePtr_);
  }
  else
  {
    loader_ = tmpLoaderPtr;
  }

  delete tiaDataStore_;
  tiaDataStore_ = new N_TIA_DataStore(&tiaParams_, linearSystem_);
  nlsMgrPtr_->registerTIADataStore(tiaDataStore_);

  delete stepErrorControl_;
  stepErrorControl_ = new N_TIA_StepErrorControl(commandLine_, *this, tiaParams_);

  if (mpdeMgrPtr_)
  {
    mpdeMgrPtr_->registerApplicationLoader(loader_);
  }

  // Now that data store has been created, we can also create the
  // working integration method object.
  delete workingIntgMethod_;
  workingIntgMethod_ = new N_TIA_WorkingIntegrationMethod(tiaParams_, *stepErrorControl_, *tiaDataStore_);
  stepErrorControl_->registerWIMPtr(workingIntgMethod_);

  // allocate nls equation loader (must be after sec and wim allocation)
  delete nonlinearEquationLoaderPtr_;
  nonlinearEquationLoaderPtr_ = new N_LOA_NonlinearEquationLoader(*tiaDataStore_, *loader_, *workingIntgMethod_, *pdsMgrPtr_, daeStateDerivFlag_);
  nonlinearEquationLoaderPtr_->registerDeviceInterface(devInterfacePtr_); 

  nlsMgrPtr_->registerLoader(nonlinearEquationLoaderPtr_);

  tiaMPDEIfacePtr_->registerTIADataStore(&*tiaDataStore_);
  tiaMPDEIfacePtr_->registerTIAStepErrorControl(&*stepErrorControl_);
  // The final time has to be forced to be a "paused" breakpoint,
  // after any registerTIAParams call.
  // RLS: May need to conditionally call this only when this is an MPDE run
  // not sure if we know that when this routine is called.
  //setPauseTime(tiaParams_.finalTime);

  registerOutputIntervals();
  registerRestartIntervals();

  // register current:
  linearSystem_->registerCurrStaVector(&(tiaDataStore_->currStatePtr));
  linearSystem_->registerCurrSolVector(&(tiaDataStore_->currSolutionPtr));

  // register next:
  linearSystem_->registerNextStaVector(&(tiaDataStore_->nextStatePtr));
  linearSystem_->registerNextSolVector(&(tiaDataStore_->nextSolutionPtr));

  // register last:
  linearSystem_->registerLastStaVector(&(tiaDataStore_->lastStatePtr));
  linearSystem_->registerLastSolVector(&(tiaDataStore_->lastSolutionPtr));

  linearSystem_->registerFlagSolVector(&(tiaDataStore_->flagSolutionPtr));

  // register temporaries:
  linearSystem_->registerTmpSolVector(&(tiaDataStore_->tmpSolVectorPtr));
  linearSystem_->registerTmpStaVector(&(tiaDataStore_->tmpStaVectorPtr));
  linearSystem_->registerTmpStaDerivVector(&(tiaDataStore_->tmpStaDerivPtr));
  linearSystem_->registerTmpStaDivDiffVector(&(tiaDataStore_->tmpStaDivDiffPtr));

  // register next derivatives:
  linearSystem_->registerNextSolDerivVector(&(tiaDataStore_->nextSolutionDerivPtr));
  linearSystem_->registerNextStaDerivVector(&(tiaDataStore_->nextStateDerivPtr));

  // register the device mask
  linearSystem_->registerDeviceMaskVector(&(tiaDataStore_->deviceMaskPtr));

  // Get the RHS and the Jacobian
  tiaDataStore_->JMatrixPtr    = linearSystem_->getJacobianMatrix();
  tiaDataStore_->RHSVectorPtr  = linearSystem_->getRHSVector();

  // DAE formulation vectors
  linearSystem_->registerDAEQVector     ( tiaDataStore_->daeQVectorPtr );
  linearSystem_->registerDAEFVector     ( tiaDataStore_->daeFVectorPtr );
  linearSystem_->registerDAEBVector     ( tiaDataStore_->daeBVectorPtr );

  // DAE formulation matrices
  linearSystem_->registerDAEdQdxMatrix  ( tiaDataStore_->dQdxMatrixPtr );
  linearSystem_->registerDAEdFdxMatrix  ( tiaDataStore_->dFdxMatrixPtr );

  // Get the limiter vectors
  tiaDataStore_->dFdxdVpVectorPtr = linearSystem_->getdFdxdVpVector ();
  tiaDataStore_->dQdxdVpVectorPtr = linearSystem_->getdQdxdVpVector ();

  tiaDataStore_->limiterFlag = loader_->getLimiterFlag ();

  // This should probably be moved elsewhere later.  If the user has
  // specified that steps should only be accepted when the nonlinear solver
  // truly converges, then the "nearConvergence" return code needs to be
  // negative.  ERK. 7/02/03
  if ( !(tiaParams_.nlNearConvFlag) )
  {
    N_NLS_ReturnCodes retCodes;
    retCodes.nearConvergence = -3;
    nlsMgrPtr_->setReturnCodes(retCodes);
  }

  // same for the "small update" case.
  if ( !(tiaParams_.nlSmallUpdateFlag) )
  {
    N_NLS_ReturnCodes retCodes;
    retCodes.smallUpdate = -4;
    nlsMgrPtr_->setReturnCodes(retCodes);
  }

  // check if analysis was specified.  If not, but .OP was specified,
  // then set up a DC calculation.
  if ( !analysisParamsRegistered )
  {
    if ( dotOpSpecified_ )
    {
      analysisMode_ = ANP_MODE_DC_SWEEP;
      analysisParamsRegistered = true;
    }
    else // flag an error.
    {
      Report::UserError0() << "No analysis statement in the netlist";
      return false;
    }
  }

  // Allocate analysis objects, and also set up params.
  allocateAnalysisObject_ ();

  initializeAllFlag_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::run
// Purpose       : Execute the control loop for the set analysis type.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
bool AnalysisManager::run()
{
  std::string msg;

  if (initializeAllFlag_ == false)
  {
    Report::DevelFatal0().in("AnalysisManager::run")
      << "Call the initializeAll function first";
  }

  if (!analysisParamsRegistered)
  {
    Report::UserError0() << "No analysis statement in the netlist";
    return false;
  }

  // now that the problem is set up, we can tell the datastore what the var type are
  // this info can be used by the time integrator to control error in a way that is
  // appropriate for different variable types (like V's and I's).
  std::vector<char> varTypes;
  topoMgrPtr_->returnVarTypeVec( varTypes );
  tiaDataStore_->varTypeVec = varTypes;

  bool runStatus = false;
  {
    // This prepares the outputters for this analysis mode
    IO::ActiveOutput active(outputManagerAdapter_->getOutputManager());

    active.setStepSweep(outputManagerAdapter_->getStepParamVec());
    active.setDCSweep(outputManagerAdapter_->getDCParamVec());
    active.add(pdsMgrPtr_->getPDSComm()->comm(), analysisMode_);

#ifdef Xyce_VERBOSE_TIME
    tiaParams_.printParams(Xyce::lout(), anpAnalysisModeToNLS(analysisMode_));
#endif

#ifndef Xyce_NO_MASKED_WRMS_NORMS
    if (loader_->loadDeviceMask())
    {
      //Xyce::dout() << " Nontrivial mask! " << std::endl;
    }
    else
    {
      // Xyce::dout() << " Trivial mask! " << std::endl;
    }
#endif

    // 6/24/2010: ERK Note: this needs a refactor!
    // For HB and MPDE, it is necessary to reallocate different analysis types as the simulation
    // progresses fr//                     om the initial condition phase to the block analysis phase.  For example,
    // in MPDE, it is common to solve a DC, then a series of transients, and then the full MPDE
    // simulation.
    //
    // I recently moved some stuff from the analysis manager down into specific analysis type
    // classes, which is where they really belong.  However, doing this required that the
    // order of setup change somewhat.  DC and sweep parameters now cannot be processed until
    // the DCSweep or Step classes are allocated.  They are now primarily allocated
    // in the intializeAll function.  They need to be allocated prior to this ::run function,
    // because they need to happen before restart files are read in.
    //
    // However, to preserve MPDE and HB, there has to be an option of reallocating analysis
    // classes here.  So, that is still the case, but this should be refactored to make
    // it cleaner.
    if ( getBlockAnalysisFlag () )
    {
      allocateAnalysisObject_ ();
    }

    Report::safeBarrier(pdsMgrPtr_->getPDSComm()->comm());

    solverStartTime_ = elapsedTimerPtr_->elapsedTime();

    // Start the solvers timers.
    xyceTranTimerPtr_->resetStartTime();

    runStatus = analysisObject_->run();

    if (tiaParams_.condTestFlag)
    {
      conductanceTest ();
    }
  }

  return runStatus;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::allocateAnalysisObject_
// Purpose       : Allocate analysis objects, and also setup params.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/10
//-----------------------------------------------------------------------------
void AnalysisManager::allocateAnalysisObject_ ()
{
  std::string msg("");

  if( !tiaParams_.resume )
  {
    if (analysisMode_ == ANP_MODE_TRANSIENT )
    {
      analysisObject_ = Teuchos::rcp(new Transient(*this));
      analysisObject_->setAnalysisParams(tranParamsBlock);
      stepErrorControl_->resetAll();
    }
    else if (analysisMode_ == ANP_MODE_DC_SWEEP)
    {
      analysisObject_ = Teuchos::rcp(new DCSweep(*this));
      for (int i=0;i<(int)dcParamsBlockVec.size();++i)
      {
        analysisObject_->setAnalysisParams(dcParamsBlockVec[i]);
      }
    }
    else if (analysisMode_ == ANP_MODE_MPDE)
    {
      analysisObject_ = Teuchos::rcp(new MPDE(*this));
      analysisObject_->setAnalysisParams(mpdeParamsBlock);
    }
    else if (analysisMode_ == ANP_MODE_HB)
    {
      analysisObject_ = Teuchos::rcp(new HB(*this));
      analysisObject_->setAnalysisParams(hbParamsBlock);
      Teuchos::rcp_dynamic_cast<HB>(analysisObject_)->setHBOptions(hbOptionsBlock);
      Teuchos::rcp_dynamic_cast<HB>(analysisObject_)->setHBLinSol(hbLinSolBlock);
      Teuchos::rcp_dynamic_cast<HB>(analysisObject_)->setLinSol(linSolBlock);
      setHBFlag( true );
    }
    else if ( analysisMode_ == ANP_MODE_AC)
    {
      analysisObject_ = Teuchos::rcp(new AC(*this));
      analysisObject_->setAnalysisParams(acParamsBlock);
    }
    else if ( analysisMode_ == ANP_MODE_MOR)
    {
      analysisObject_ = Teuchos::rcp(new MOR(*this));
      analysisObject_->setAnalysisParams(morParamsBlock);
    }
    else
    {
      Report::UserError0() <<  "Unknown type of analysis";
      return;
    }

    if( !is_null(analysisObject_) )
    {
      analysisObject_->setParamsWithOutputMgrAdapter ( *outputManagerAdapter_ );
    }

    // need to throw an error here as this really shouldn't happen.
    if( is_null( analysisObject_ ) )
    {
      Report::DevelFatal0().in("AnalysisManager::initializeAll")
        << "Unable to allocate analysis type";
    }

    // ultimately sensitivity calculation should get its own analysis object.
    if (sensFlag_)
    {
      analysisObject_->setSensFlag();
    }

    // The primary analysis object is stored in case there is an outer
    // loop analysis such as .STEP or .DAKOTA.   Many accessors such
    // as "getStepNumber" need to access this object rather than the
    // step object.
    primaryAnalysisObject_ = analysisObject_;

    if( stepLoopFlag_ )
    {
      stepAnalysisTarget_ = analysisObject_;
      analysisObject_ = Teuchos::rcp(new Step(*this, *stepAnalysisTarget_.get()));
      for (int i=0;i<(int)stepParamsBlockVec.size();++i)
      {
        analysisObject_->setAnalysisParams(stepParamsBlockVec[i]);
      }
      analysisObject_->setParamsWithOutputMgrAdapter ( *outputManagerAdapter_ );
    }

    if (dakotaRunFlag_ && is_null(dakotaAnalysisTarget_) )
    {
      dakotaAnalysisTarget_ = analysisObject_;
      analysisObject_ = Teuchos::rcp(new Dakota(*this, *dakotaAnalysisTarget_.get()));
      analysisObject_->setAnalysisParams(dakotaParamsBlock);
      analysisObject_->setParamsWithOutputMgrAdapter ( *outputManagerAdapter_ );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::printLoopInfo
// Purpose       : Prints out time loop information.
// Special Notes : Prints stats from save point start to save point finish.
//                 Special case 0,0 is entire run to this point
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
bool AnalysisManager::printLoopInfo(int start, int finish)
{
  return primaryAnalysisObject_->printLoopInfo (start,finish);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::testDCOPOutputTime_
// Purpose       : Similar to testOutputTime, except that this is for
//                 DCOP restart files.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::testDCOPOutputTime_()
{
  bool flag(true);

  if( !dcopRestartFlag_ )
  {
    flag = false;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::testSaveOutputTime_
// Purpose       : Similar to testOutputTime, except that this is for
//                 .SAVE files.
// Special Notes : Only outputs 1x.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::testSaveOutputTime_()
{
  bool flag(true);

  if( !saveFlag_ )
  {
    flag = false;
  }
  else if (stepErrorControl_->currentTime < saveTime_)
  {
    flag = false;
  }
  else if (savedAlready_)
  {
    flag = false;
  }

  if (flag==true)
  {
    savedAlready_ = true;
    Xyce::dout() <<"Calling SAVE outputs!" <<std::endl;
  }

  return flag;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::testRestartSaveTime_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel ComputationalSciences.
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::testRestartSaveTime
()
{
  bool flag;

#ifdef Xyce_DEBUG_RESTART
  Xyce::dout() << "TESTING FOR RESTART SAVE" << std::endl
               << Xyce::subsection_divider << std::endl
               << "stepErrorControl_->currentTime: " << stepErrorControl_->currentTime << std::endl
               << "nextSaveTime: " << nextRestartSaveTime_ << std::endl
               << "initialRestartInterval_: " << initialRestartInterval_ << std::endl;
  if (!(restartIntervals_.empty()))
  {
    Xyce::dout() << "First restart interval: " << restartIntervals_[0].first << std::endl;
  }
  else
  {
    Xyce::dout() << "restartIntervals_ is empty" << std::endl;
  }
#endif

  if (initialRestartInterval_ == 0.0)
  {
    flag = false;
  }
  else if (stepErrorControl_->currentTime < nextRestartSaveTime_)
  {
    flag = false;
  }
  else if (restartIntervals_.empty())
  {
    while (nextRestartSaveTime_ <= stepErrorControl_->currentTime)
    {
      nextRestartSaveTime_ += initialRestartInterval_;
    }
    flag = true;
  }
  else if (stepErrorControl_->currentTime < restartIntervals_[0].first)
  {
    while (nextRestartSaveTime_ <= stepErrorControl_->currentTime)
    {
      nextRestartSaveTime_ += initialRestartInterval_;
    }
    if (nextRestartSaveTime_ > restartIntervals_[0].first)
    {
      nextRestartSaveTime_ = restartIntervals_[0].first;
    }
    flag = true;
  }
  else
  {
    std::pair<double, double> currInterval, nextInterval;
    int size = restartIntervals_.size();
    for (int i = 0; i < size; ++i)
    {
      if (restartIntervals_[i].first <= stepErrorControl_->currentTime)
      {
        currInterval = restartIntervals_[i];
        if ((i+1) < (int)restartIntervals_.size())
        {
          nextInterval = restartIntervals_[i+1];
        }
      }
    }
    int step = static_cast <int> ((stepErrorControl_->currentTime-currInterval.first) /
                                  currInterval.second);
    nextRestartSaveTime_ = currInterval.first + (step+1)*currInterval.second;

    if (nextInterval.first && (nextInterval.first!=currInterval.first)
        && (nextRestartSaveTime_>=nextInterval.first))
    {
      nextRestartSaveTime_ = nextInterval.first;
    }
    flag = true;
  }

#ifdef Xyce_DEBUG_RESTART
  Xyce::dout() << "new nextSaveTime: " << nextRestartSaveTime_ << std::endl
               << "restart flag: " << flag << std::endl
               << Xyce::subsection_divider << std::endl;
#endif

  return flag;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::partialTimeDerivative
// Purpose       : Returns the current partial time derivative for either the
//                 solution or state vector.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
double AnalysisManager::partialTimeDerivative()
{
  double partDT = workingIntgMethod_->partialTimeDeriv();

  // Add a check here to try to prevent capacitive spiral of death.
  // This is still a "research" option, so is not on by default.
  if (tiaParams_.jacLimitFlag)
  {
#ifdef Xyce_DEBUG_ANALYSIS
    if (tiaParams_.debugLevel > 0)
    {
      Xyce::dout() << "AnalysisManager::partialTimeDerivative.";
      Xyce::dout() << "   Using jac limit = " << tiaParams_.jacLimit << std::endl;
    }
#endif
    if (partDT > tiaParams_.jacLimit)
    {
      partDT = tiaParams_.jacLimit;
    }
  }

  return partDT;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBreakpointTol
// Purpose       : Returns the breakpoint tolerance.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/7/02
//-----------------------------------------------------------------------------
double AnalysisManager::getBreakpointTol()
{
  return N_UTL_BreakPoint::getBPTol();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setBreakpointTol
// Purpose       : Sets the breakpoint tolerance.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/7/02
//-----------------------------------------------------------------------------
void AnalysisManager::setBreakpointTol(double bptol)
{
  N_UTL_BreakPoint::setBPTol(bptol);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerRestartIntervals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::registerRestartIntervals()
{
  double initInt;
  std::vector< std::pair<double,double> > intPairs;
  restartPtr_->getRestartIntervals(initInt, intPairs);
  initialRestartInterval_ = initInt;
  nextRestartSaveTime_    = stepErrorControl_->initialTime;
  restartIntervals_       = intPairs;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerOutputIntevals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::registerOutputIntervals()
{
  double initInt;
  std::vector< std::pair<double,double> > intPairs;
  outputManagerAdapter_->getOutputIntervals( initInt, & intPairs );
  initialOutputInterval_ = initInt;
  nextOutputTime_ = stepErrorControl_->initialTime;
  outputIntervals_ = intPairs;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTranAnalysisParams
// Purpose       : Sets transient analysis parameters (from .TRAN)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool AnalysisManager::setTranAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysisMode_ = ANP_MODE_TRANSIENT;
  tranParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDCAnalysisParams
// Purpose       : Sets the DC sweep calculation parameters (from .DC)
//
// Special Notes : This function will be called multiple times if there is more
//                 than one sweep variable.  The parser separates each variable
//                 into separate option blocks prior to calling this function.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::setDCAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysisMode_ = ANP_MODE_DC_SWEEP;
  // before saving a copy of paramsBlock, check and see if it already
  // is in the dcParamsBlockVec.  If so then replace the original copy.
  // This replacement is necessary if the first copy had an expression
  // element that was resolved in later parsing.

  bool foundMatch = false;
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecItr = dcParamsBlockVec.begin();
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecEnd = dcParamsBlockVec.end();
  while( paramsBlockVecItr != paramsBlockVecEnd )
  {
    if( paramsBlockVecItr->compareParamLists( paramsBlock ) )
    {
      // these are the same
      foundMatch = true;
      break;
    }
    paramsBlockVecItr++;
  }

  if( foundMatch )
  {
    // replace the existing one with the new one
    *paramsBlockVecItr = paramsBlock;
  }
  else
  {
    // save the new one.
    dcParamsBlockVec.push_back (paramsBlock); // save a copy for later.
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setOPAnalysisParams
// Purpose       : Handle OP statement. (.OP)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
bool AnalysisManager::setOPAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  dotOpSpecified_ = true;
  opParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSTEPAnalysisParams
// Purpose       : Sets the STEP calculation parameters. (from .STEP)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/03
//-----------------------------------------------------------------------------
bool AnalysisManager::setSTEPAnalysisParams(
  const N_UTL_OptionBlock & paramsBlock)
{
  stepLoopFlag_ = true;
  // before saving a copy of paramsBlock, check and see if it already
  // is in the stepParamsBlockVec.  If so then replace the original copy.
  // This replacement is necessary if the first copy had an expression
  // element that was resolved in later parsing.
  bool foundMatch = false;
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecItr = stepParamsBlockVec.begin();
  std::vector<N_UTL_OptionBlock>::iterator paramsBlockVecEnd = stepParamsBlockVec.end();
  while( paramsBlockVecItr != paramsBlockVecEnd )
  {
    if( paramsBlockVecItr->compareParamLists( paramsBlock ) )
    {
      // these are the same
      foundMatch = true;
      break;
    }
    paramsBlockVecItr++;
  }

  if( foundMatch )
  {
    // replace the existing one with the new one
    *paramsBlockVecItr = paramsBlock;
  }
  else
  {
    // save the new one.
    stepParamsBlockVec.push_back (paramsBlock); // save a copy for later.
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSaveOptions
// Purpose       : Sets the Save parameters.
// Special Notes : Most of these parameters are handled in the output manager,
//                 rather than here.  So, most params are a noop here, except
//                 for "TIME".
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::setSaveOptions(
  const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams_.debugLevel > 0)
  {
    Xyce::dout() << "In AnalysisManager::setSaveOptions" << std::endl;
  }
#endif

  saveFlag_ = true;

  Util::ParameterList::const_iterator iterPL = OB.getParams().begin();
  Util::ParameterList::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
    Xyce::dout() << "iterPL->tag = " << iterPL->tag() << std::endl;
#endif
    if (iterPL->tag() == "TYPE")
    {
      // noop
    }
    else if (iterPL->tag() == "FILE")
    {
      // noop
    }
    else if (iterPL->tag() == "TIME")
    {
      saveTime_ = iterPL->getImmutableValue<double>();
      saveTimeGiven_ = true;
    }
    else if (iterPL->tag() == "LEVEL")
    {
      // noop
    }
    else
    {
      // noop.  If it gets here there is an error on the .SAVE line.
      // However, the IO manager has a trap for this, so do nothing
      // here.
    }

    ++iterPL;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setACAnalysisParams
// Purpose       : Sets the AC sweep calculation parameters (from .AC)
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::setACAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysisMode_ = ANP_MODE_AC;
  acParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setMORAnalysisParams
// Purpose       : Sets the MOR calculation parameters (from .MOR)
//
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 05/29/12
//-----------------------------------------------------------------------------
bool AnalysisManager::setMORAnalysisParams
  (const N_UTL_OptionBlock & paramsBlock)
{
  analysisParamsRegistered = true;
  analysisMode_ = ANP_MODE_MOR;
  morParamsBlock = paramsBlock; // save a copy for later.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setMOROptions
// Purpose       :
// Special Notes : These are from '.options mor'
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 05/31/12
//-----------------------------------------------------------------------------
bool AnalysisManager::setMOROptions(const N_UTL_OptionBlock & OB)
{
  Util::ParameterList::const_iterator it_tpL;
  Util::ParameterList::const_iterator first = OB.getParams().begin();
  Util::ParameterList::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="METHOD")
    {
      ExtendedString stringVal ( it_tpL->stringValue() );
      stringVal.toUpper();
      tiaParams_.morMethod = stringVal;
    }
    else if (it_tpL->uTag()=="SAVEREDSYS")
    {
      tiaParams_.morSaveRedSys = static_cast<bool>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="COMPORIGTF")
    {
      tiaParams_.morCompOrigTF = static_cast<bool>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="COMPREDTF")
    {
      tiaParams_.morCompRedTF = static_cast<bool>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="COMPTYPE")
    {
      ExtendedString stringVal ( it_tpL->stringValue() );
      stringVal.toUpper();
      tiaParams_.morCompType = stringVal;
    }
    else if (it_tpL->uTag()=="COMPNP")
    {
      tiaParams_.morCompNP = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="COMPFSTART")
    {
      tiaParams_.morCompFStart = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="COMPFSTOP")
    {
      tiaParams_.morCompFStop = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="EXPPOINT")
    {
      tiaParams_.morExpPoint = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="SCALETYPE")
    {
      tiaParams_.morScaleType = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="SCALEFACTOR")
    {
      tiaParams_.morScaleFactor = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="SCALEFACTOR1")
    {
      tiaParams_.morScaleFactor1 = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="SPARSIFICATIONTYPE")
    {
      tiaParams_.morSparsificationType = it_tpL->getImmutableValue<int>();
    }
    else
    {
      Report::UserError0() << it_tpL->uTag() << " is not a recognized model-order reduction option.";
    }
  }

  // If we are computing the transfer function, make sure the frequency range is valid.
  if (tiaParams_.morCompOrigTF || tiaParams_.morCompRedTF)
  {
    if (tiaParams_.morCompFStop < tiaParams_.morCompFStart)
    {
      Report::UserError() << ".options mor COMPFSTART = " << tiaParams_.morCompFStart << " > " << tiaParams_.morCompFStop << " = COMPFSTOP!";
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDCOPRestartParams
// Purpose       : Sets the DCOP restart parameters.
// Special Notes : Most of the dcop restart parameters are used and handled by
//                 the IO::OutputMgr class, so this function here doens't need
//                 to do very much.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisManager::setDCOPRestartParams(const N_UTL_OptionBlock & OB)
{
  dcopRestartFlag_ = true;

  Util::ParameterList::const_iterator iterPL = OB.getParams().begin();
  Util::ParameterList::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
    Xyce::dout() << "iterPL->tag = " << iterPL->tag() << std::endl;
#endif

    if (iterPL->tag() == "INPUT")
    {
      // do nothing, this is handled in the output manager.
    }
    else if (iterPL->tag() == "OUTPUT")
    {
      // do nothing, this is handled in the output manager.
    }
    else if (iterPL->tag() == "TIME")
    {
      saveTime_ = iterPL->getImmutableValue<double>();
      saveTimeGiven_ = true;
    }
    else
    {
      // noop.  If it gets here there is an error on the .SAVE line.
      // However, the IO manager has a trap for this, so do nothing
      // here.
    }

    ++iterPL;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTranOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::setTranOptions(const N_UTL_OptionBlock & OB)
{
  Util::ParameterList::const_iterator it_tpL;
  Util::ParameterList::const_iterator first = OB.getParams().begin();
  Util::ParameterList::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag()=="METHOD")
    {
//      tiaParams_.integrationMethod=it_tpL->iVal();

      if (it_tpL->isInteger())
        tiaParams_.integrationMethod=it_tpL->getImmutableValue<int>();
      else
      {

        ExtendedString stringVal ( it_tpL->stringValue() );
        stringVal.toUpper();

        if (stringVal == "TRAP" || stringVal == "TRAPEZOIDAL")
          tiaParams_.integrationMethod = 7;
        else if (stringVal == "BDF")
          tiaParams_.integrationMethod = 6;
        else if (stringVal == "GEAR")
          tiaParams_.integrationMethod = 8;
        else
        {
          Report::UserError0() << "Unsupported transient method type";
        }
      }

    }
#ifdef Xyce_DEBUG_ANALYSIS
    else if (it_tpL->uTag()=="CONSTSTEP")
    {
      tiaParams_.constantStepSize
        =  static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
#endif
    else if (it_tpL->uTag()=="USEDEVICEMAX")
    {
      tiaParams_.useDeviceTimeStepMax =  static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="RELTOL")
    {
      tiaParams_.relErrorTol=it_tpL->getImmutableValue<double>();
      tiaParams_.relErrorTolGiven = true;
    }
    else if (it_tpL->uTag()=="ABSTOL")
    {
      tiaParams_.absErrorTol=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="DOUBLEDCOPSTEP")
    {
      tiaParams_.doubleDCOPStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="FIRSTDCOPSTEP")
    {
      tiaParams_.firstDCOPStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="LASTDCOPSTEP")
    {
      tiaParams_.lastDCOPStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="BPENABLE" )
    {
      tiaParams_.bpEnable =  static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="RESTARTSTEPSCALE" )
    {
      tiaParams_.restartTimeStepScale=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="EXITTIME" )
    {
      tiaParams_.exitTime=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="EXITSTEP" )
    {
      tiaParams_.exitStep=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag() == "MINTIMESTEPSBP")
    {
      tiaParams_.minTimeStepsBP = it_tpL->getImmutableValue<int>();
      tiaParams_.minTimeStepsBPGiven = true;
    }
    else if (it_tpL->uTag()=="ERROPTION" )
    {
      tiaParams_.errorAnalysisOption=it_tpL->getImmutableValue<int>();
      // tscoffe/tmei 12/7/07:  If error option = 1 (NO LTE) then make sure minTimeStepBP is enabled.
      if (tiaParams_.errorAnalysisOption == 1)
      {
        if (!tiaParams_.minTimeStepsBPGiven)
        {
          tiaParams_.minTimeStepsBPGiven = true;
        }
      }
    }
    else if (it_tpL->uTag()=="NLNEARCONV" )
    {
      tiaParams_.nlNearConvFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="NLSMALLUPDATE" )
    {
      tiaParams_.nlSmallUpdateFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="JACLIMITFLAG" )
    {
      tiaParams_.jacLimitFlag= static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="JACLIMIT" )
    {
      tiaParams_.jacLimit = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="DAESTATEDERIV" )
    {
      daeStateDerivFlag_ = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="MAXORD" )
    {
      tiaParams_.maxOrder = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="MINORD" )
    {
      tiaParams_.minOrder = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="TIMESTEPSREVERSAL" )
    {
      tiaParams_.timestepsReversal =it_tpL->getImmutableValue<bool>();
    }
    else if (it_tpL->uTag()=="TESTFIRSTSTEP" )
    {
      tiaParams_.testFirstStep = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="DELMAX" )
    {
      tiaParams_.delmax =it_tpL->getImmutableValue<double>();
      tiaParams_.delmaxGiven = true;
//      tiaParams_.delmax =it_tpL->iVal();
    }
    else if (it_tpL->uTag()=="NLMIN" )
    {
      tiaParams_.NLmin=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="NLMAX" )
    {
      tiaParams_.NLmax=it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="OUTPUTINTERPMPDE")
    {
      tiaParams_.outputInterpMPDE = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="NEWLTE")
    {
      tiaParams_.newLte =it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="NEWBPSTEPPING")
    {
      tiaParams_.newBPStepping = it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="INTERPOUTPUT")
    {
      tiaParams_.interpOutputFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="CONDTEST")
    {
      tiaParams_.condTestFlag = static_cast<bool> (it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag()=="CONDTESTDEVICENAME")
    {
      tiaParams_.condTestDeviceNames.push_back(it_tpL->stringValue() );
    }
    else if (it_tpL->uTag() == "DTMIN")
    {
      tiaParams_.userSpecMinTimeStep = it_tpL->getImmutableValue<double>();
      tiaParams_.userSpecMinTimeStepGiven = true;
    }
    else if (it_tpL->uTag() == "PASSNLSTALL")
    {
      tiaParams_.passNLStall = it_tpL->getImmutableValue<bool>();
    }
    else if (it_tpL->uTag() == "FASTTESTS")
    {
      tiaParams_.fastTests = it_tpL->getImmutableValue<bool>();
    }
    else if (it_tpL->uTag() == "MINTIMESTEPRECOVERY")
    {
      tiaParams_.minTimeStepRecoveryCounter = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag()=="VOLTZEROTOL" )
    {
      tiaParams_.voltZeroTol=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="CURRZEROTOL" )
    {
      tiaParams_.currZeroTol=it_tpL->getImmutableValue<double>();
    }
    else if (it_tpL->uTag()=="DEBUGLEVEL" )
    {
#ifdef Xyce_DEBUG_TIME
      // set (or override) debug levels based on command line options
      if ( commandLine_.argExists( "-tdl" ) )
      {
        tiaParams_.debugLevel = atoi( commandLine_.getArgumentValue( "-tdl" ).c_str() );
      }

      else
      {
        tiaParams_.debugLevel = it_tpL->getImmutableValue<int>();
      }
#endif
    }
    else if (it_tpL->uTag()=="HISTORYTRACKINGDEPTH" )
    {
      tiaParams_.historyTrackingDepth = it_tpL->getImmutableValue<int>();
    }
    else
    {
      Report::UserError() << it_tpL->uTag() << " is not a recognized time integration option";
    }
  }

  if (tiaParams_.NLmin > tiaParams_.NLmax)
  {
    Report::UserError() << ".options timeint NLMIN = " << tiaParams_.NLmin << " > " << tiaParams_.NLmax << " = NLMAX!";
  }

  if (tiaParams_.firstDCOPStep < 0) tiaParams_.firstDCOPStep = 0;
  if (tiaParams_.firstDCOPStep > 1) tiaParams_.firstDCOPStep = 1;
  if (tiaParams_.lastDCOPStep < 0)  tiaParams_.lastDCOPStep  = 0;
  if (tiaParams_.lastDCOPStep > 1)  tiaParams_.lastDCOPStep  = 1;

  // check for maximum order on the command line
  if ( commandLine_.argExists ("-maxord") )
  {

    tiaParams_.maxOrder =
      atoi( commandLine_.getArgumentValue( "-maxord" ).c_str() );
  }

  if ( tiaParams_.maxOrder < 1) tiaParams_.maxOrder = 1;
  if ( tiaParams_.maxOrder > 5) tiaParams_.maxOrder = 5;

  if (tiaParams_.newLte == true)
  {
    if (tiaParams_.relErrorTolGiven != true)
      tiaParams_.relErrorTol = 1.0e-3;
  }

  return true;

}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setMPDEAnalysisParams
// Purpose       : Sets the MPDE sweep calculation parameters (from .MPDE)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setMPDEAnalysisParams(const N_UTL_OptionBlock & OB)
{
  if (!mpdeMgrPtr_)
  {
    setupMPDEMgr_();
  }
  bool bsuccess = mpdeMgrPtr_->setMPDEAnalysisParams(OB);
  analysisMode_ = ANP_MODE_MPDE;
  analysisParamsRegistered = true;
  mpdeParamsBlock = OB;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setMPDEOptions
// Purpose       :
// Special Notes : from '.options mpdeint'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setMPDEOptions(const N_UTL_OptionBlock & OB)
{
  if (!mpdeMgrPtr_)
  {
    setupMPDEMgr_();
  }
  mpdeMgrPtr_->setMPDEOptions(OB);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setHBAnalysisParams
// Purpose       : Sets the HB sweep calculation parameters (from .HB)
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setHBAnalysisParams(const N_UTL_OptionBlock & OB)
{
  Util::ParameterList::const_iterator it_tpL;
  Util::ParameterList::const_iterator first = OB.getParams().begin();
  Util::ParameterList::const_iterator last  = OB.getParams().end();

  for (it_tpL = first; it_tpL != last; ++it_tpL)
  {
    if (it_tpL->uTag() == "FREQ")
    {

      tiaParams_.freqs = it_tpL->getValue<std::vector<double> >();
      tiaParams_.freqGiven = true;
    }
  }

  if (tiaParams_.freqs[0] <= 0.0 )
  {
    Report::UserError() << "Frequency of oscillation " << tiaParams_.freqs[0] << " is less than or equal to zero, invalid .HB specification";
  }

  if ((DEBUG_ANALYSIS & tiaParams_.debugLevel) > 0)
  {
    dout() << section_divider << std::endl
           << "HB transient simulation parameters" 
           //<< Util::push << std::endl
           << std::endl
           << "HB frequency = " << tiaParams_.freqs[0] << std::endl
           //<< Util::pop << std::endl;
           << std::endl;
  }

  setBlockAnalysisFlag( true );
  devInterfacePtr_->setBlockAnalysisFlag( true );
  analysisMode_ = ANP_MODE_HB;
  analysisParamsRegistered = true;
  hbParamsBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setHBOptions
// Purpose       :
// Special Notes : from '.options hb'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setHBOptions(const N_UTL_OptionBlock & OB)
{
  hbOptionsBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
bool AnalysisManager::setLinSol(const N_UTL_OptionBlock & OB)
{
  linSolBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setHBLinSol(const N_UTL_OptionBlock & OB)
{
  hbLinSolBlock = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTRANMPDEOptions
// Purpose       :
// Special Notes : from '.options timeint-mpde'
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::setTRANMPDEOptions(const N_UTL_OptionBlock & OB)
{
  if (!mpdeMgrPtr_)
  {
    setupMPDEMgr_();
  }
  mpdeMgrPtr_->registerTranMPDEOptions(OB);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/5/13
//-----------------------------------------------------------------------------
bool AnalysisManager::setSensOptions(const N_UTL_OptionBlock & OB)
{
  sensFlag_=true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::completeOPStartStep
// Purpose       : Call to rotate next state to current state following a
//               : constrained DCOP solve when using a previous operating point
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 06/27/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::completeOPStartStep  ( )
{
  bool bsuccess = true;

  bsuccess = tiaDataStore_->updateStateDataArrays ();
  tiaDataStore_->setConstantHistory();
  tiaDataStore_->equateTmpVectors();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::completeHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 successful homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::completeHomotopyStep
    ( const std::vector<std::string> & paramNames,
      const std::vector<double> & paramVals,
      N_LAS_Vector * solnVecPtr )
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  std::string netListFile = commandLine_.getArgumentValue("netlist");
  Xyce::dout() << "\n " << netListFile;
  Xyce::dout() << " AnalysisManager::completeHomotopyStep " << std::endl;
#endif

  // Rotate the data vectors:
  bool bs1 = tiaDataStore_->updateStateDataArrays ();    bsuccess = bsuccess && bs1;
  tiaDataStore_->setConstantHistory();
  tiaDataStore_->equateTmpVectors();

  // Pass info in to the lower level solver:
  loader_->homotopyStepSuccess (paramNames,paramVals);

  // Call output
  outputManagerAdapter_->outputHomotopy( paramNames, paramVals, *solnVecPtr );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::failHomotopyStep
// Purpose       : Call to do final cleanup after a single,
//                 failed homotopy step.
// Special Notes : Called from the Xyce-LOCA interface.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//-----------------------------------------------------------------------------
bool AnalysisManager::failHomotopyStep ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_ANALYSIS
  std::string netListFile = commandLine_.getArgumentValue("netlist");
  Xyce::dout() << "\n " << netListFile;
  Xyce::dout() << " AnalysisManager::failHomotopyStep " << std::endl;
#endif

  // Pass info in to the lower level solver:
  loader_->homotopyStepFailure ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setupMPDEMgr_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
void AnalysisManager::setupMPDEMgr_()
{
  // Allocate new MPDE Manager
  mpdeMgrPtr_ = new N_MPDE_Manager(commandLine_);

  // Register MPDE manager with everything
  bool bs1 = true;
  bool bsuccess = true;
  bs1 = doMPDEregistrations_();
  bsuccess = bsuccess && bs1;

  if (!bsuccess)
  {
    Report::DevelFatal0().in("AnalysisManager::setupMPDEMgr_") << "Registration function failed";
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;

  // this work was in the constructor, but now we don't know the pkgOptMgrPtr_ until
  // it is registered.  So, do this work now.
  std::string netListFile = "";
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT", netListFile, new AnalysisManager_TranOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TRAN", netListFile, new AnalysisManager_TransAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "RESULT", netListFile, new AnalysisManager_ResultOptionsReg( *this ) );

  pkgOptMgrPtr_->submitRegistration(
      "DC", netListFile, new AnalysisManager_DCAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "OP", netListFile, new AnalysisManager_OPAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "STEP", netListFile, new AnalysisManager_STEPAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "OP_IO", netListFile, new AnalysisManager_DCOPOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SAVE", netListFile, new AnalysisManager_SaveOptionsReg( this ) );

  // MPDE specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MPDE", netListFile, new AnalysisManager_MPDE_AnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "MPDEINT", netListFile, new AnalysisManager_MPDE_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT-MPDE", netListFile, new AnalysisManager_MPDE_TranMPDEOptionsReg( this ) );

  // HB Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "HB", netListFile, new AnalysisManager_HB_AnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "HBINT", netListFile, new AnalysisManager_HB_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "LINSOL-HB", netListFile, new AnalysisManager_HB_LinSolReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "LINSOL", netListFile, new AnalysisManager_LinSolReg( this ) );

  // AC Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "AC", netListFile, new AnalysisManager_AC_AnalysisReg( this ) );

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR", netListFile, new AnalysisManager_MOR_AnalysisReg( this ) );

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR_OPTS", netListFile, new AnalysisManager_MOR_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SENS", netListFile, new AnalysisManager_SensOptionsReg( this ) );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::doMPDEregistrations_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter
// Creation Date : 5/23/2014
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisManager::doMPDEregistrations_()
{
  bool bs=true;

  if (mpdeMgrPtr_)
  {
    mpdeMgrPtr_->registerAnalysisManager(this);
    mpdeMgrPtr_->registerDeviceInterface(devInterfacePtr_);
    mpdeMgrPtr_->registerParallelManager(pdsMgrPtr_);
    mpdeMgrPtr_->registerTopology(topoMgrPtr_);
    mpdeMgrPtr_->registerRestartManager(restartPtr_);
    mpdeMgrPtr_->registerOutputManager(outMgrPtr_);
    mpdeMgrPtr_->registerApplicationLoader(loader_);
    mpdeMgrPtr_->registerNonlinearEquationLoader(nonlinearEquationLoaderPtr_);
    mpdeMgrPtr_->registerApplicationBuilder(appBuilderPtr_);
    mpdeMgrPtr_->registerLinearSystem(linearSystem_);

    mpdeMgrPtr_->registerTIAMPDEInterface(tiaMPDEIfacePtr_);
  }
  else
  {
    bs=false;
  }
  return bs;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::initializeTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
void AnalysisManager::initializeTransientModel()
{
  workingIntgMethod_->createTimeIntegMethod(TIAMethod_BACKWARD_DIFFERENTIATION_15);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::evalTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//-----------------------------------------------------------------------------
bool AnalysisManager::evalTransientModel(
    double t,
    N_LAS_Vector * SolVectorPtr,
    N_LAS_Vector * CurrSolVectorPtr,
    N_LAS_Vector * LastSolVectorPtr,
    N_LAS_Vector * StaVectorPtr,
    N_LAS_Vector * CurrStaVectorPtr,
    N_LAS_Vector * LastStaVectorPtr,
    N_LAS_Vector * StaDerivVectorPtr,
    N_LAS_Vector * StoVectorPtr,
    N_LAS_Vector * CurrStoVectorPtr,
    N_LAS_Vector * LastStoVectorPtr,
    N_LAS_Vector * stoLeadCurrQVectorPtr,
    N_LAS_Vector * QVectorPtr,
    N_LAS_Vector * FVectorPtr,
    N_LAS_Vector * BVectorPtr,
    N_LAS_Vector * dFdxdVpVectorPtr,
    N_LAS_Vector * dQdxdVpVectorPtr,
    N_LAS_Matrix * dQdxMatrixPtr,
    N_LAS_Matrix * dFdxMatrixPtr
    )
{
  // This is F,Q load:
  bool bsuccess = loader_->loadDAEVectors(
        SolVectorPtr,
        CurrSolVectorPtr,
        LastSolVectorPtr,
        StaVectorPtr,
        CurrStaVectorPtr,
        LastStaVectorPtr,
        StaDerivVectorPtr,
        StoVectorPtr,
        CurrStoVectorPtr,
        LastStoVectorPtr,
        stoLeadCurrQVectorPtr,
        QVectorPtr,
        FVectorPtr,
        BVectorPtr,
        dFdxdVpVectorPtr,
        dQdxdVpVectorPtr
        );
  // This is dQdx, dFdx load:
  bsuccess = bsuccess && loader_->loadDAEMatrices(
      SolVectorPtr,
      StaVectorPtr,
      StaDerivVectorPtr,
      StoVectorPtr,
      dQdxMatrixPtr,
      dFdxMatrixPtr);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::evalTransientModelState
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : private
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
bool AnalysisManager::evalTransientModelState(
    double t,
    N_LAS_Vector * SolVectorPtr,
    N_LAS_Vector * StaVectorPtr,
    N_LAS_Vector * StoVectorPtr
    )
{
  // This is part of state vector load:
  stepErrorControl_->currentTime = t;
  stepErrorControl_->nextTime = t;
  loader_->updateSources(); // updates source values wrt time t.
  bool bsuccess = loader_->updateState(
      SolVectorPtr,
      SolVectorPtr,
      SolVectorPtr,
      StaVectorPtr,
      StaVectorPtr,
      StaVectorPtr,
      StoVectorPtr,
      StoVectorPtr,
      StoVectorPtr
      );
  return bsuccess;
}

// ***** Accessor methods *****
//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setBeginningIntegrationFlag(bool bif)
{
  primaryAnalysisObject_->setBeginningIntegrationFlag(bif);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getBeginningIntegrationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getBeginningIntegrationFlag()
{
  return primaryAnalysisObject_->getBeginningIntegrationFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setIntegrationMethod(int im)
{
  primaryAnalysisObject_->setIntegrationMethod (im);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
unsigned int AnalysisManager::getIntegrationMethod ()
{
  return primaryAnalysisObject_->getIntegrationMethod ();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalLinearSolutionTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalLinearSolutionTime() const
{
  return primaryAnalysisObject_->getTotalLinearSolutionTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalResidualLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalResidualLoadTime() const
{
  return primaryAnalysisObject_->getTotalResidualLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTotalJacobianLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
double AnalysisManager::getTotalJacobianLoadTime() const
{
  return primaryAnalysisObject_->getTotalJacobianLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDoubleDCOPEnabled ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getDoubleDCOPEnabled ()
{
  return primaryAnalysisObject_->getDoubleDCOPEnabled ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::isSimulationComplete
// Purpose       : return boolean signifying whether simulation complete or
//                 not.
//
// Special Notes :THIS VERSION IS ONLY VALID FOR TRANSIENT RUNS, where
//                 completion of the simulation means integration to final
//                 time.
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//-----------------------------------------------------------------------------
bool AnalysisManager::isSimulationComplete()
{
  if (analysisMode_ == ANP_MODE_TRANSIENT)
  {
    return (stepErrorControl_->isFinished());
  }
  else
  {
    Report::DevelFatal0().in("AnalysisManager::simulationComplete") << "Called for non-transient run, not currently valid";
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getPauseTime
//
// Purpose       : return the time at which the simulation will pause
//
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
double AnalysisManager::getPauseTime()
{
  return(tiaParams_.pauseTime);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
int AnalysisManager::getStepNumber ()
{
  int number=0;
  if( !is_null( primaryAnalysisObject_ ) )
  {
    number = primaryAnalysisObject_->getStepNumber();
  }
  return number;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
int AnalysisManager::getTranStepNumber ()
{
  int number=0;
  if (analysisMode_ == ANP_MODE_TRANSIENT && !is_null(primaryAnalysisObject_) )
  {
    number = primaryAnalysisObject_->getTranStepNumber();
  }
  return number;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void AnalysisManager::setStepNumber (int step)
{
  if( !is_null( primaryAnalysisObject_ ) )
  {
    primaryAnalysisObject_->setStepNumber(step);
  }
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setTranStepNumber
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
void AnalysisManager::setTranStepNumber (int step)
{
  if( !is_null( primaryAnalysisObject_ ) )
  {
    primaryAnalysisObject_->setTranStepNumber(step);
  }
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitTranFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 12/16/2010
//-----------------------------------------------------------------------------
bool AnalysisManager::getInitTranFlag()
{
  int stepNum =  getStepNumber();
  return (stepNum <= 0);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTime
// Purpose       : Gets the next time value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getTime() const
{
  return stepErrorControl_->nextTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getCurrentTime
// Purpose       : Gets the current time value.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 8/21/2009
//-----------------------------------------------------------------------------
double AnalysisManager::getCurrentTime() const
{
  return stepErrorControl_->currentTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getFinalTime
// Purpose       : Gets the final time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getFinalTime() const
{
  return stepErrorControl_->finalTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInitialTime
// Purpose       : Gets the initial time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getInitialTime() const
{
  return stepErrorControl_->initialTime;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStartingTimeStep
// Purpose       : Gets the starting time-step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisManager::getStartingTimeStep ()
{
  return stepErrorControl_->startingTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTimeIntMode
// Purpose       : Gets the time-integration method.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
int AnalysisManager::getTimeIntMode()
{
  return primaryAnalysisObject_->getIntegrationMethod();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDCOPFlag
// Purpose       : Gets a flag indicating we are in steady state.
//                  (steady=true)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getDCOPFlag ()
{
  if ( !is_null(primaryAnalysisObject_) )
  {
    return (primaryAnalysisObject_->getDCOPFlag());
  }
  return false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getInputOPFlag
// Purpose       : Gets a flag indicating we are starting from a previous OP
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/13/06
//-----------------------------------------------------------------------------
bool AnalysisManager::getInputOPFlag ()
{
  return primaryAnalysisObject_->getInputOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTranOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a transient simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getTranOPFlag ()
{
  return ((analysisMode_ == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
     && (primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getACOPFlag
// Purpose       : Gets a flag indicating we are in a DCOP calculation that is
//                 the initialization for a transient simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/16/12
//-----------------------------------------------------------------------------
bool AnalysisManager::getACOPFlag ()
{
  return ((analysisMode_ == ANP_MODE_AC || primaryAnalysisObject_->isAnalysis(ANP_MODE_AC))
     && (primaryAnalysisObject_->getIntegrationMethod())==TIAMethod_NONE);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDCSweepFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getDCSweepFlag ()
{
  return ((analysisMode_ == ANP_MODE_DC_SWEEP || primaryAnalysisObject_->isAnalysis(ANP_MODE_DC_SWEEP)));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTransientFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/05
//-----------------------------------------------------------------------------
bool AnalysisManager::getTransientFlag () const
{
  return (((analysisMode_ == ANP_MODE_TRANSIENT || primaryAnalysisObject_->isAnalysis(ANP_MODE_TRANSIENT))
     && (primaryAnalysisObject_->getIntegrationMethod()) !=TIAMethod_NONE));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDoubleDCOPStep
// Purpose       : Gets the double DC Operating Point step value.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/01
//-----------------------------------------------------------------------------
int AnalysisManager::getDoubleDCOPStep()
{
  return primaryAnalysisObject_->getDoubleDCOPStep();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getCurrentStepSize
// Purpose       : Returns the "current" time step size.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/15/01
//-----------------------------------------------------------------------------
double AnalysisManager::getCurrentStepSize()
{
  return stepErrorControl_->currentTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getLastStepSize
// Purpose       : Returns the "last" time step size.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/15/01
//-----------------------------------------------------------------------------
double AnalysisManager::getLastStepSize()
{
  return stepErrorControl_->lastTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::updateDerivs
// Purpose       : Calls the time  int. method to update the corrector
//                 derivatives.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool AnalysisManager::updateDerivs()
{
  workingIntgMethod_->obtainCorrectorDeriv();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::updateDivDiffs
// Purpose       : Updates the divided difference values.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
bool AnalysisManager::updateDivDiffs()
{
  tiaDataStore_->computeDividedDifferences();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::updateDerivsBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
bool AnalysisManager::updateDerivsBlock(
  const std::list<index_pair> & solGIDList,
  const std::list<index_pair> & staGIDList)
{
  tiaDataStore_->computeDivDiffsBlock(solGIDList,staGIDList);
  workingIntgMethod_->updateDerivsBlock (solGIDList, staGIDList);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::equateTmpVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool AnalysisManager::equateTmpVectors()
{
  return tiaDataStore_->equateTmpVectors();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerElapsedTimer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/24/06
//-----------------------------------------------------------------------------
bool AnalysisManager::registerElapsedTimer(Util::Timer * et)
{
  elapsedTimerPtr_ = et;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::restartDataSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
int AnalysisManager::restartDataSize( bool pack )
{
  return stepErrorControl_->restartDataSize( pack );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::dumpRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::dumpRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  return stepErrorControl_->dumpRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::restoreRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::restoreRestartData
  (char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack)
{
  return stepErrorControl_->restoreRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getSolnVarData( const int & gid,
                                             std::vector<double> & varData )
{
  return tiaDataStore_->getSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::getStateVarData( const int & gid,
                                              std::vector<double> & varData )
{
  return tiaDataStore_->getStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getStoreVarData( const int & gid,
                                              std::vector<double> & varData )
{
  return tiaDataStore_->getStoreVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getOrder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 02/15/07
//-----------------------------------------------------------------------------
int AnalysisManager::getOrder ()
{
  return workingIntgMethod_->getOrder();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNumberOfSteps
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int AnalysisManager::getNumberOfSteps ()
{
  return workingIntgMethod_->getNumberOfSteps();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getUsedOrder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int AnalysisManager::getUsedOrder ()
{
  return workingIntgMethod_->getUsedOrder();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getNscsco
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Computational Sciences
// Creation Date : 03/01/07
//-----------------------------------------------------------------------------
int AnalysisManager::getNscsco ()
{
  return workingIntgMethod_->getNscsco();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::setSolnVarData( const int & gid,
                                             const std::vector<double> & varData )
{
  return tiaDataStore_->setSolnVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisManager::setStateVarData( const int & gid,
                                              const std::vector<double> & varData )
{
  return tiaDataStore_->setStateVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::setStoreVarData( const int & gid,
                                              const std::vector<double> & varData )
{
  return tiaDataStore_->setStoreVarData( gid, varData );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 09/06/01
//-----------------------------------------------------------------------------
bool AnalysisManager::registerTIAParams(const N_TIA_TIAParams & tia_params)
{
  tiaParams_ = tia_params;

  if (stepErrorControl_)
    stepErrorControl_->setTIAParams();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerLinearSystem
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerLinearSystem(N_LAS_System * linear_system)
{
  linearSystem_ = linear_system;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerNLSManager
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/17/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerNLSManager(N_NLS_Manager * nlsMgrPtr_tmp)
{
  nlsMgrPtr_ = nlsMgrPtr_tmp;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerLoader
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerLoader(N_LOA_Loader * loader)
{
  loader_ = loader;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerOutputMgr
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerOutputMgr(IO::OutputMgr * output_manager)
{
  outMgrPtr_ = output_manager;
  outputManagerAdapter_ = new OutputMgrAdapter(pdsMgrPtr_->getPDSComm()->comm(), *this);

  outputManagerAdapter_->registerOutputMgr( outMgrPtr_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerRestartMgr
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerRestartMgr(
  IO::RestartMgr * restartPtr_tmp)
{
  restartPtr_ = restartPtr_tmp;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerDeviceInterface ( N_DEV_DeviceInterface * devInterfacePtr_tmp )
{
  devInterfacePtr_ = devInterfacePtr_tmp;
  return true;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerTopology( N_TOP_Topology * topoMgrPtr_tmp )
{
  topoMgrPtr_ = topoMgrPtr_tmp;
  
  return true;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerApplicationBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisManager::registerApplicationBuilder( N_LAS_Builder * appBuilderPtr_tmp )
{
  appBuilderPtr_ = appBuilderPtr_tmp;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::registerParallelServices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/19/02
//-----------------------------------------------------------------------------
bool AnalysisManager::registerParallelServices(N_PDS_Manager * pds)
{
  pdsMgrPtr_ = pds;

  xyceTranTimerPtr_ = new Util::Timer(*pds->getPDSComm());

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setNextSolVectorPtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 04/22/2003
//-----------------------------------------------------------------------------
bool AnalysisManager::setNextSolVectorPtr (N_LAS_Vector * solVecPtr)
{
  return tiaDataStore_->setNextSolVectorPtr (solVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setPauseTime
//
// Purpose       : Set the time at which to pause the simulation
//
// Special Notes : This is a temporary implementation that I'm using to
//                 begin the process of hiding time integrator internals from
//                 the "N_CIR_Xyce::simulateUntil" method, so that I can
//                 ultimately change those internals without the simulateUntil
//                 method
//                 In the zero order version, just sets tiaParams_.pauseTime
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
void AnalysisManager::setPauseTime(double pauseTime)
{
  stepErrorControl_->setBreakPoint(N_UTL_BreakPoint(pauseTime, Util::PAUSE_BREAKPOINT));
}



//-----------------------------------------------------------------------------
// Function      : AnalysisManager::resumeSimulation
//
// Purpose       : set flag to signify that simulation is continuation of
//                 previously paused simulation.
//
// Special Notes : This is a temporary implementation that I'm using to
//                 begin the process of hiding time integrator internals from
//                 the "N_CIR_Xyce::simulateUntil" method, so that I can
//                 ultimately change those internals without the simulateUntil
//                 method
//                 In the zero order version, just sets tiaParams_.resume=true
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/2004
//-----------------------------------------------------------------------------
void AnalysisManager::resumeSimulation()
{
  tiaParams_.resume=true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::unset_resumeSimulation
// Purpose       : UN-set flag to signify that simulation is continuation of
//                 previously paused simulation.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 06/12/2013
//-----------------------------------------------------------------------------
void AnalysisManager::unset_resumeSimulation()
{
  tiaParams_.resume=false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::outputIntervalSpecified_
// Purpose       : Return true if user has specified output control options
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems modeling
// Creation Date : 2/12/07
//-----------------------------------------------------------------------------
bool AnalysisManager::outputIntervalSpecified()
{
  // This is necessary and sufficient ---
  return (initialOutputInterval_ > 0.0);
}

// routines to get/set Dakota run flags and actually run a Dakota iteration
//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisManager::getDakotaRunFlag()
{
  return dakotaRunFlag_;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setDakotaRunFlag( bool flag )
{
  dakotaRunFlag_ = flag;
}


//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
int AnalysisManager::getDakotaIteration()
{
  return dakotaIterationNumber_;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::setDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::setDakotaIteration( int iterNumber )
{
    dakotaIterationNumber_ = iterNumber;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::getTimeIntInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisManager::getTimeIntInfo (N_TIA_TimeIntInfo & tiInfo)
{
  tiInfo.currentOrder         = getOrder ();
  tiInfo.numberOfSteps        = getNumberOfSteps ();
  tiInfo.usedOrder            = getUsedOrder ();
  tiInfo.nscsco               = getNscsco ();

  tiInfo.pdt                  = partialTimeDerivative(); // alpha/DT
  tiInfo.nextTimeStep         = getCurrentStepSize();
  tiInfo.currTimeStep         = getLastStepSize();
  tiInfo.currentTime          = stepErrorControl_->currentTime;
  tiInfo.nextTime             = getTime();
  tiInfo.finalTime            = getFinalTime();
  tiInfo.startingTimeStep     = getStartingTimeStep();
  tiInfo.bpTol                = getBreakpointTol();

  tiInfo.dcopFlag             = getDCOPFlag ();
  tiInfo.inputOPFlag          = getInputOPFlag ();
  tiInfo.tranopFlag           = getTranOPFlag ();
  tiInfo.acopFlag             = getACOPFlag ();
  tiInfo.transientFlag        = getTransientFlag ();
  tiInfo.dcsweepFlag          = getDCSweepFlag ();
  tiInfo.sweepSourceResetFlag = getSweepSourceResetFlag ();

  tiInfo.timeStepNumber       = getStepNumber();
  tiInfo.initTranFlag         = getInitTranFlag ();
  tiInfo.beginIntegrationFlag = getBeginningIntegrationFlag();

  tiInfo.doubleDCOPStep       = getDoubleDCOPStep();
  tiInfo.doubleDCOPEnabled    = getDoubleDCOPEnabled ();

  tiInfo.stepLoopIter         = 0;
  if( stepLoopFlag_ )
  {
    analysisObject_->getStepIter();
  }

  tiInfo.timeIntMode          = getTimeIntMode();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::silenceProgress
// Purpose       : Shut up "Percent Complete" noises
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void AnalysisManager::silenceProgress()
{
  progressFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisManager::enableProgress
// Purpose       : stop shutting up "Percent Complete" noises
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void AnalysisManager::enableProgress()
{
  progressFlag_ = true;
}

N_PDS_Comm &
AnalysisManager::getPDSComm() {
  return *pdsMgrPtr_->getPDSComm();
}

//-----------------------------------------------------------------------------
// Function      : anpAnalysisModeToNLS
// Purpose       : Converte between N_NLS_Manager.h AnalysisMode enum and
//               : expanded ANP_AnalysisManager.h ANP_Analysis_Mode enum.
// Special Notes :
// Scope         :
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
Nonlinear::AnalysisMode anpAnalysisModeToNLS(Analysis_Mode mode)
{
  Nonlinear::AnalysisMode outMode;
  if (mode == ANP_MODE_TRANSIENT)
  {
    outMode = Nonlinear::TRANSIENT;
  }
  else if (mode == ANP_MODE_DC_OP)
  {
    outMode = Nonlinear::DC_OP;
  }
  else if (mode == ANP_MODE_DC_SWEEP)
  {
    outMode = Nonlinear::DC_SWEEP;
  }
  else if (mode == ANP_MODE_HB)
  {
    outMode = Nonlinear::HB_MODE;
  }
  else
  {
    outMode = Nonlinear::NUM_MODES; // Should be this be TRANSIENT?
  }
  return outMode;
}

} // namespace Analysis
} // namespace Xyce
