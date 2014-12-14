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
// Filename      : $RCSfile: N_ANP_HB.C,v $
// Purpose       : HB analysis functions.
// Special Notes :
// Creator       : Todd Coffey, 1414, Ting Mei, 1437
// Creation Date : 07/23/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.94.2.2 $
// Revision Date  : $Date: 2014/08/28 21:00:43 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_UTL_APFT.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_HB.h>
#include <N_ANP_Transient.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_Report.h>
#include <N_MPDE_Manager.h>
#include <N_MPDE_Discretization.h>

#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_DataStore.h>

#include <N_LOA_HBLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_PrecondFactory.h>
#include <N_LAS_System.h>
#include <N_LAS_BlockSystemHelpers.h>

#include <N_NLS_Manager.h>

// #include <N_IO_OutputMgr.h>
#include <N_IO_ActiveOutput.h>

#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Timer.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseHelpers.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : HB::HB( AnalysisManager * )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
HB::HB( AnalysisManager &analysis_manager )
  : AnalysisBase(analysis_manager),
    Util::ListenerAutoSubscribe<StepEvent>(&analysis_manager),
  debugLevel(0),
  isPaused(false),
  startDCOPtime(0.0),
  endTRANtime(0.0),
  isTransient_(false),
  isDCSweep_(false),
  test_(false),
  size_(21),
  period_(1.0),
  startUpPeriods_(0),
  startUpPeriodsGiven_(false),
  startUpPeriodsFinished_(false),
  saveIcData_(false), 
  transTiaParams_( AnalysisBase::tiaParams_),
  voltLimFlag_(1),
  intmodMax_(0),
  method_("APFT"),
  taHB_(1),
  intmodMaxGiven_(false),
  fastTimeDisc_(0),
  fastTimeDiscOrder_(1),
  resetForStepCalledBefore_(false),
  hbLoaderPtr_(0)
{
  devInterfacePtr_ = analysisManager_.getDeviceInterface();
  nonlinearEquationLoaderPtr_ = analysisManager_.getNonlinearEquationLoader();
  appBuilderPtr_ = analysisManager_.getAppBuilder();
  pdsMgrPtr_ = analysisManager_.getPDSManager();
}

HB::~HB()
{
  delete hbLoaderPtr_;
}

void HB::notify(const StepEvent &event) 
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();

    if (resetForStepCalledBefore_)
    {
      goodSolutionVec_.clear();
      goodStateVec_.clear();
      goodQVec_.clear();
      goodStoreVec_.clear();

      stepErrorControl_->resetAll();

      analysisManager_.setNextOutputTime(0.0);
      analysisManager_.registerLinearSystem(&linearSystem_);
      analysisManager_.resumeSimulation();

      nonlinearSolverManager_.resetAll(Nonlinear::DC_OP);
      nonlinearSolverManager_.setMatrixFreeFlag( false );
      nonlinearSolverManager_.registerLinearSystem(&linearSystem_);
      nonlinearSolverManager_.registerLoader( nonlinearEquationLoaderPtr_);
      nonlinearSolverManager_.setLinSolOptions( saved_lsOB_ );

      devInterfacePtr_->registerLinearSystem(&linearSystem_);

      // un-set the fast source flag in the devices
      std::vector<std::string> srcVec;
      devInterfacePtr_->deRegisterFastSources( srcVec );

      analysisManager_.initializeAll(&loader_);

      devInterfacePtr_->initializeAll();
      devInterfacePtr_->setMPDEFlag( false );

      nonlinearSolverManager_.initializeAll();

      analysisManager_.unset_resumeSimulation();

      analysisManager_.getXyceTranTimer().resetStartTime();
    }

    resetForStepCalledBefore_=true;
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::getStepNumber()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
int HB::getStepNumber ()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    return analysisObject_->getStepNumber();
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : HB::setStepNumber()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
void HB::setStepNumber (int step)
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    analysisObject_->setStepNumber( step );
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::setBeginningIntegrationFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
void HB::setBeginningIntegrationFlag(bool bif)
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    analysisObject_->setBeginningIntegrationFlag( bif );
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::getBeginningIntegrationFlag()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
bool HB::getBeginningIntegrationFlag()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    return analysisObject_->getBeginningIntegrationFlag();
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
void HB::setIntegrationMethod(int im)
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    analysisObject_->setIntegrationMethod( im );
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::getIntegrationMethod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/20/13
//-----------------------------------------------------------------------------
unsigned int HB::getIntegrationMethod ()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    return analysisObject_->getIntegrationMethod ();
  }
  return TIAMethod_NONE;
}


//-----------------------------------------------------------------------------
// Function      : HB::getDoubleDCOPStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
int HB::getDoubleDCOPStep()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    return analysisObject_->getDoubleDCOPStep();
  }
  else
  {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::getDCOPFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/24/2014
//-----------------------------------------------------------------------------
bool HB::getDCOPFlag ()
{
  if ( !Teuchos::is_null( analysisObject_ ) )
  {
    if ( isTransient_ )
    {
      return analysisObject_->getDCOPFlag ();
    }
    // DCSweep is a special case in the HB analysis type.  The HB calculation is
    // performed through a DC analysis object.  However, this function (getDCOPFlag)
    // is called (ultimately) from the device package to determine a bunch of
    // state-dependent load decisions.  For an HB calculation, it should NOT
    // do a DCOP load.  It needs to do a full transient load.
    else 
    {
      return false;
    }
  }
  else
  {
    // the assumption here is that if there is no analysis object, then
    // the simulation hasn't started yet.  When it *does* start, the first
    // analysis will be DCOP.
    //return true;
    return false;
  }
}

//-----------------------------------------------------------------------------
// Function      : HB::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::run()
{

  // initializeAll_();
  //   get TIAParams from AnalysisManager
  //   create HBBuilder
  //     generateMaps
  //     generateStateMaps
  //   create vectors for IC
  //
  // if (test_) {
  //   runTests_();
  // } else {
  //
  //   computeInitialCondition_();
  //     run startup periods
  //     integrate one period
  //     interpolate to evenly spaced points
  //     set this as IC for HB nonlinear problem
  //
  //   setupHBProblem_();
  //
  //   runHBProblem_();
  //
  // }

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
// Function      : HB::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::init()
{
  bool returnValue=true;

  Xyce::lout() << " ***** Running HB initial conditions....\n" << std::endl;
#ifdef Xyce_DEBUG_HB
  Xyce::dout() << std::endl
         << section_divider << std::endl
         << "  HB::init()" << std::endl;
#endif // Xyce_DEBUG_HB

  //Store copy of transient TIAParams for HB run
  transTiaParams_ = analysisManager_.getTIAParams();
  doubleDCOPFlag_ = loader_.getDoubleDCOPFlag();
  doubleDCOPStep_ = transTiaParams_.firstDCOPStep;

  setFreqPoints_();
  //
  // If it was requested, advance the solution a fixed number of startup periods

  period_ =   1.0/freqPoints_[ (size_-1)/2 + 1 ]; 

#ifdef Xyce_DEBUG_HB 
  Xyce::dout() << "HB period =" <<  period_ << std::endl;
#endif

  if (transTiaParams_.freqs.size() == 1)
    setInitialGuess_();
  else
  {
    setTimePoints_();

    fastTimes_.resize(size_+1);
    goodTimePoints_.resize(size_+1);
    
    goodTimePoints_ = fastTimes_;

  }

  // now that we have size_, continue with the initialization of objects for HB
  mpdeDiscPtr_ = rcp(new N_MPDE_Discretization(
     static_cast<N_MPDE_Discretization::Type>(fastTimeDisc_), fastTimeDiscOrder_ ));

  hbBuilderPtr_ = rcp(new N_LAS_HBBuilder( size_, mpdeDiscPtr_ ));

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::init():  Generate Maps\n";
#endif // Xyce_DEBUG_HB
  hbBuilderPtr_->generateMaps( rcp(pdsMgrPtr_->getParallelMap( "SOLUTION" ), false),
                               rcp(pdsMgrPtr_->getParallelMap( "SOLUTION_OVERLAP_GND" ), false) );
  hbBuilderPtr_->generateStateMaps( rcp(pdsMgrPtr_->getParallelMap( "STATE" ),false) );
  hbBuilderPtr_->generateStoreMaps( rcp(pdsMgrPtr_->getParallelMap( "STORE" ),false) );
  hbBuilderPtr_->generateGraphs( *pdsMgrPtr_->getMatrixGraph( "JACOBIAN" ));

  HBICVectorPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  HBICStateVectorPtr_ = hbBuilderPtr_->createTimeDomainStateBlockVector();
  HBICStoreVectorPtr_ = hbBuilderPtr_->createTimeDomainStoreBlockVector();
  HBICQVectorPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();

  HBICVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();

//  HBICStateVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeStateBlockVector();
//  HBICQVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();
//  HBICStoreVectorFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeStoreBlockVector();

  transTiaParams_.finalTime =  period_;
  analysisManager_.registerTIAParams (transTiaParams_);

  // set the fast source flag in the devices
  std::vector<std::string> srcVec;
  devInterfacePtr_->registerFastSources( srcVec );

  //analysisManager_.anaIntPtr->getnlHBOptions(saved_nlHBOB_);
  nonlinearSolverManager_.getHBOptions(saved_nlHBOB_);

  // Create HB Loader.
  delete hbLoaderPtr_;
  hbLoaderPtr_ = new N_LOA_HBLoader( mpdeState_, mpdeDiscPtr_ );
  hbLoaderPtr_->registerHBBuilder(hbBuilderPtr_);
  hbLoaderPtr_->registerAppBuilder(appBuilderPtr_);

  // Create DFT for HB Loader
  // NOTE:  For single-tone HB the DFT will probably be a FFT, for multi-tone a specialized 
  //        implementation of the Util::DFTInterfaceDecl will need to be made and registered with
  //        the HB loader.

  ftInData_.resize( size_ );
  ftOutData_.resize( size_ +1 );
  iftInData_.resize( size_  +1 );
  iftOutData_.resize( size_ );
  if (transTiaParams_.freqs.size() == 1)
  {
    if (ftInterface_ == Teuchos::null)
    {
      ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( size_ ) );
      ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
    } 
    else if (ftInterface_->getFFTInterface()->getSignalLength() != size_)
    {
      ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( size_ ) );
      ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
    }
    hbLoaderPtr_->registerDFTInterface( ftInterface_->getFFTInterface() );
  }
  else
  {
    createFT_();

    dftInterface_ = Teuchos::rcp( new N_UTL_APFT<std::vector<double> >( idftMatrix_, dftMatrix_ ) );
    dftInterface_->registerVectors( Teuchos::rcp( &ftInData_, false ), Teuchos::rcp( &ftOutData_, false ), Teuchos::rcp( &iftInData_, false ), Teuchos::rcp( &iftOutData_, false ) );
//    ftInterface_-> registerMatrices(idftMatrix_, dftMatrix_);

    hbLoaderPtr_->registerDFTInterface( dftInterface_ );
  } 

  if (taHB_==1)
  {
  // Pick up IC data from the initial transient.
  for (int i=0 ; i<size_ ; ++i)
  {
#ifdef Xyce_DEBUG_HB
    if( debugLevel > 0 )
    {
      Xyce::dout() << "HB::init():  Loading initial condition data from time: fastTimes_["
                << i << "] = " << fastTimes_[i] << std::endl;
    }
#endif // Xyce_DEBUG_HB

    HBICVectorPtr_->block(i) = *(goodSolutionVec_[i]);
    HBICStateVectorPtr_->block(i) = *(goodStateVec_[i]);
    HBICQVectorPtr_->block(i) = *(goodQVec_[i]);
    HBICStoreVectorPtr_->block(i) = *(goodStoreVec_[i]);
  }

  hbLoaderPtr_->permutedFFT(*HBICVectorPtr_, &*HBICVectorFreqPtr_);

  }
#ifdef Xyce_DEBUG_HB
  if ( debugLevel > 1 )
  {
    Xyce::dout() << "HB Initial Condition Solution!\n";
    HBICVectorPtr_->printPetraObject(std::cout);
    Xyce::dout() << "HB Initial Condition State Vector!\n";
    HBICStateVectorPtr_->printPetraObject(std::cout);
    Xyce::dout() << "HB Initial Condition Store Vector!\n";
    HBICStoreVectorPtr_->printPetraObject(std::cout);
  }
#endif // Xyce_DEBUG_HB

  //Destroy Solvers, etc. from IC phase and prepare for HB
  //-----------------------------------------

//  if ( voltLimFlag_ == 0 )
//    devInterfacePtr_->setVoltageLimiterFlag (false);
  devInterfacePtr_->setVoltageLimiterFlag (voltLimFlag_);

  devInterfacePtr_->setMPDEFlag( true );
  analysisManager_.resetAll();

  //-----------------------------------------

  //Finish setup of HB Loader
  //-----------------------------------------
  hbLoaderPtr_->registerAppLoader( rcp(&loader_, false) );
  hbLoaderPtr_->registerDeviceInterface( devInterfacePtr_ );
  goodTimePoints_.resize(size_+1);
  goodTimePoints_[size_] = period_;

  for( int i = 0; i < size_; ++i )
  {
    goodTimePoints_[i] = goodTimePoints_[i] - transTiaParams_.initialTime;
  }

  hbLoaderPtr_->setFastTimes(goodTimePoints_);

  //-----------------------------------------
  //Construct Solvers, etc. for HB Phase
  //-----------------------------------------
  lasHBSysPtr_ = rcp(new N_LAS_System());
  //-----------------------------------------

  analysisManager_.registerTIAParams( transTiaParams_ );
  analysisManager_.registerLinearSystem( &*lasHBSysPtr_ );

  //hack needed by TIA initialization currently
  hbBuilderPtr_->registerPDSManager( pdsMgrPtr_ );

  lasHBSysPtr_->registerAnalysisManager(&analysisManager_);
  lasHBSysPtr_->registerPDSManager( pdsMgrPtr_ );
  lasHBSysPtr_->registerBuilder( &*hbBuilderPtr_ );

  //need to cut out unnecessary stuff from this call for new dae
  lasHBSysPtr_->initializeSystem();

  // Give NLS Manager the same old nonlinearEquationLoader as it just calls the TIA loader in newDAE
  nonlinearSolverManager_.registerLinearSystem( &*lasHBSysPtr_ );
  nonlinearSolverManager_.setLinSolOptions( saved_lsHBOB_ );
  nonlinearSolverManager_.setMatrixFreeFlag( true );

  // Let the HB loader know that the application of the operator is matrix free
  hbLoaderPtr_->setMatrixFreeFlag( true );

  if (Teuchos::is_null( precFactory_ ))
  {
    // Generate the HB preconditioner factory.
    precFactory_ = rcp( new N_LAS_HBPrecondFactory( saved_lsHBOB_ ) );
  }

  // Register application loader with preconditioner factory
  RCP<N_LAS_HBPrecondFactory> tmpPrecFactory
    = rcp_dynamic_cast<N_LAS_HBPrecondFactory>( precFactory_ );

  tmpPrecFactory->registerAppBuilder( appBuilderPtr_ );
  tmpPrecFactory->registerHBLoader( rcp(hbLoaderPtr_, false) );
  tmpPrecFactory->registerHBBuilder( hbBuilderPtr_ );
  tmpPrecFactory->setFastTimes( goodTimePoints_ );
  tmpPrecFactory->setTimeSteps( timeSteps_ );

  nonlinearSolverManager_.registerPrecondFactory( precFactory_ );
  //-----------------------------------------

  //Initialization of Solvers, etc. for HB Phase
  //-----------------------------------------
  //Dummy call to setup time integrator for transient
  analysisManager_.resumeSimulation();
  analysisManager_.initializeAll( hbLoaderPtr_ );
  analysisManager_.setMPDEFlag( true );

  nonlinearSolverManager_.initializeAll();
  nonlinearSolverManager_.setAnalysisMode(anpAnalysisModeToNLS(ANP_MODE_HB));

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << section_divider << std::endl;
#endif // Xyce_DEBUG_HB

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::loopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::loopProcess()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Beginning full HB simulation....\n" << std::endl;

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << std::endl
         << section_divider << std::endl
         << "  HB::loopProcess" << std::endl;
#endif // Xyce_DEBUG_HB

  N_TIA_DataStore * dsPtr = analysisManager_.getTIADataStore();
  *(dsPtr->nextSolutionPtr) = *(HBICVectorFreqPtr_.get());
  *(dsPtr->nextStatePtr) = *(HBICStateVectorPtr_.get());
  *(dsPtr->nextStorePtr) = *(HBICStoreVectorPtr_.get());

  // try to run the problem
  analysisObject_ = Teuchos::rcp(new DCSweep(analysisManager_));
  returnValue = analysisObject_->run();

  // Add in simulation times
  accumulateStatistics_();

  // print out analysis info
  Xyce::lout() << " ***** Harmonic Balance Computation Summary *****" << std::endl;
  analysisObject_->printLoopInfo( 0, 0 );

#ifdef Xyce_DEBUG_HB
  dout() << section_divider << std::endl;
#endif // Xyce_DEBUG_HB

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::processSuccessfulDCOP()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processSuccessfulDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processSuccessfulStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processFailedStep()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::processFailedDCOP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::processFailedDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::finish()
{
  stats_ = hbStatCounts_;

  return true;
}

bool HB::handlePredictor()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::finalVerboseOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 09/18/2008
//-----------------------------------------------------------------------------
bool HB::finalVerboseOutput()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HB::setHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 05/13/13
//-----------------------------------------------------------------------------
bool HB::setHBOptions(const Util::OptionBlock & OB)
{
  Util::ParameterList::const_iterator iterPL = OB.getParams().begin();
  Util::ParameterList::const_iterator  endPL = OB.getParams().end();

  for( ; iterPL != endPL; ++iterPL )
  {
    ExtendedString tag = iterPL->tag();
    tag.toUpper();

    if (std::string(tag,0,7) == "NUMFREQ" ) 
    {
      size_ = iterPL->getImmutableValue<int>();
      numFreqs_.push_back(size_);
    }
    else if ( tag == "STARTUPPERIODS" )
    {
      startUpPeriods_ = iterPL->getImmutableValue<int>();

      if (startUpPeriods_ > 0)
        startUpPeriodsGiven_ = true;
    }
    else if( tag == "SAVEICDATA" )
    {
      saveIcData_ = true;
    }
    else if( tag == "TEST" )
    {
      test_     = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if (tag == "DEBUGLEVEL" )
    {
      debugLevel = iterPL->getImmutableValue<int>();
    }
    else if ( tag == "TAHB" )
    {
      taHB_ = iterPL->getImmutableValue<int>();
    }
    else if ( tag == "VOLTLIM" )
    {
      voltLimFlag_ = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if ( tag == "INTMODMAX" )
    {
      intmodMax_ = iterPL->getImmutableValue<int>();

      if ( intmodMax_  > 0)
        intmodMaxGiven_ = true;
    }
    else if ( tag == "METHOD" )
    {
      ExtendedString stringVal ( iterPL->stringValue() );
      stringVal.toUpper();
      method_ = stringVal;
    }
    else
    {
      UserWarning(*this) << "Unrecognized HBINT option " << tag;
    }
  }

  if (numFreqs_.size() != 0 && numFreqs_.size() != transTiaParams_.freqs.size() ) 
  {
    Report::UserError() << "The size of numFreq does not match the number of tones in .hb!";
  }      


  if (numFreqs_.size() == 0 )
  {

    numFreqs_.resize( transTiaParams_.freqs.size() ) ;
    for (int i=0; i < transTiaParams_.freqs.size(); i++ )
    {
      numFreqs_[i] = size_;
    }

  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : HB::setLinSol
// Purpose       : this is needed for .STEP to work with HB
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/12/2013
//-----------------------------------------------------------------------------
bool HB::setLinSol(const Util::OptionBlock & OB)
{
  // Save the non-HB linear solver option block
  saved_lsOB_ = OB;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/13/13
//-----------------------------------------------------------------------------
bool HB::setHBLinSol(const Util::OptionBlock & OB)
{
  // Save the HB linear solver option block
  saved_lsHBOB_ = OB;

  // Generate the HB preconditioner factory.
  precFactory_ = rcp( new N_LAS_HBPrecondFactory( OB ) );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::isAnalysis
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 5/13/13
// Notes         : Alternatively, we could try to cast the analysis object
//               : However, this method is called a lot.
//-----------------------------------------------------------------------------
bool HB::isAnalysis( int analysis_type )
{
  bool returnValue = false;

  if ( analysis_type == ANP_MODE_TRANSIENT )
  {
    returnValue = isTransient_;
  }
  if ( analysis_type == ANP_MODE_DC_SWEEP )
  {
    returnValue = isDCSweep_;
  }
  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::prepareHBOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 08/20/07
//-----------------------------------------------------------------------------
void HB::prepareHBOutput(
    N_LAS_Vector & solnVecPtr,
    std::vector<double> & timePoints,
    std::vector<double> & freqPoints,
    RCP<N_LAS_BlockVector> & timeDomainSolnVec,
    RCP<N_LAS_BlockVector> & freqDomainSolnVecReal,
    RCP<N_LAS_BlockVector> & freqDomainSolnVecImag,
    RCP<N_LAS_BlockVector> & timeDomainStoreVec,
    RCP<N_LAS_BlockVector> & freqDomainStoreVecReal,
    RCP<N_LAS_BlockVector> & freqDomainStoreVecImag
    ) const
{
  N_LAS_BlockVector & blockSolVecPtr = dynamic_cast<N_LAS_BlockVector &>(solnVecPtr);

  Teuchos::RCP<N_LAS_BlockVector> bStoreVecFreqPtr_ = hbLoaderPtr_->getStoreVecFreqPtr();

  timeDomainStoreVec = hbBuilderPtr_->createTimeDomainStoreBlockVector();

  if (bStoreVecFreqPtr_->blockCount() > 0 )
  {
    hbLoaderPtr_->permutedIFT(*bStoreVecFreqPtr_, &*timeDomainStoreVec);
  } 

  //TD solution
  timeDomainSolnVec = hbBuilderPtr_->createTimeDomainBlockVector(); 
  int blockCount = timeDomainSolnVec->blockCount();
  int N = timeDomainSolnVec->block(0).globalLength(); 

  timePoints.resize(size_);

  for( int i = 0; i < size_; ++i )
  {
    timePoints[i] = fastTimes_[i] - transTiaParams_.initialTime;
  }

  freqPoints = freqPoints_;

  // Create individual block vectors to store the real and imaginary parts separately.
  Teuchos::RCP<N_PDS_ParMap> baseMap = Teuchos::rcp_const_cast<N_PDS_ParMap>( hbBuilderPtr_->getBaseSolutionMap() );
  Teuchos::RCP<N_PDS_ParMap> globalMap = createBlockParMap( blockCount, *baseMap );
  freqDomainSolnVecReal = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalMap, baseMap ) );
  freqDomainSolnVecImag = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalMap, baseMap ) );

  hbLoaderPtr_->permutedIFT(blockSolVecPtr, &*timeDomainSolnVec); 

  // Now copy over the frequency domain solution, real and imaginary parts separately, into the output vectors.
  for (int j=0; j<N; j++)
  {
    // See if this time-domain solution variable is owned by the local processor.
    // If so, this processor owns the entire j-th block of the blockSolVecPtr vector,
    // and the j-th entry of every block in the freqDomainSolnVec[Real/Imag] vector.
    int lid = baseMap->globalToLocalIndex( j );
    N_LAS_Vector& solBlock = blockSolVecPtr.block( j );

    N_LAS_Vector& realVecRef =  freqDomainSolnVecReal->block((blockCount-1)/2);
    N_LAS_Vector& imagVecRef =  freqDomainSolnVecImag->block((blockCount-1)/2);

    if (lid >= 0)
    { 
      realVecRef[lid] = solBlock[0];  
      imagVecRef[lid] = solBlock[1];
    }

    for (int i=1; i <= (blockCount-1)/2; ++i)
    {
      N_LAS_Vector& realVecRef_neg =  freqDomainSolnVecReal->block((blockCount-1)/2 - i);
      N_LAS_Vector& imagVecRef_neg =  freqDomainSolnVecImag->block((blockCount-1)/2 - i);
      N_LAS_Vector& realVecRef_pos =  freqDomainSolnVecReal->block((blockCount-1)/2 + i);
      N_LAS_Vector& imagVecRef_pos =  freqDomainSolnVecImag->block((blockCount-1)/2 + i);

      if (lid >= 0)
      {
        realVecRef_neg[lid] = solBlock[ 2*(blockCount-i) ];
        imagVecRef_neg[lid] = solBlock[ 2*(blockCount-i) + 1 ];
        realVecRef_pos[lid] = solBlock[ 2*i ];
        imagVecRef_pos[lid] = solBlock[ 2*i+1 ];
      }
    }
  } 

  // proceed to store variables
  Teuchos::RCP<N_PDS_ParMap> baseStoreMap = Teuchos::rcp_const_cast<N_PDS_ParMap>( hbBuilderPtr_->getBaseStoreMap() );
  Teuchos::RCP<N_PDS_ParMap> globalStoreMap = createBlockParMap( blockCount, *baseStoreMap );
  freqDomainStoreVecReal = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalStoreMap, baseStoreMap ) );
  freqDomainStoreVecImag = Teuchos::rcp( new N_LAS_BlockVector( blockCount, globalStoreMap, baseStoreMap ) );

  N = timeDomainStoreVec->block(0).globalLength(); 

  for (int j=0; j<N; j++)
  {
    // See if this time-domain solution variable is owned by the local processor.
    // If so, this processor owns the entire j-th block of the blockSolVecPtr vector,
    // and the j-th entry of every block in the freqDomainSolnVec[Real/Imag] vector.
    int lid = baseStoreMap->globalToLocalIndex( j );
    N_LAS_Vector& storeBlock =  bStoreVecFreqPtr_->block( j );

    N_LAS_Vector& realVecRef =  freqDomainStoreVecReal->block((blockCount-1)/2);
    N_LAS_Vector& imagVecRef =  freqDomainStoreVecImag->block((blockCount-1)/2);

    if (lid >= 0)
    { 
      realVecRef[lid] = storeBlock[0];  
      imagVecRef[lid] = storeBlock[1];
    }

    for (int i=1; i <= (blockCount-1)/2; ++i)
    {
      N_LAS_Vector& realVecRef_neg =  freqDomainStoreVecReal->block((blockCount-1)/2 - i);
      N_LAS_Vector& imagVecRef_neg =  freqDomainStoreVecImag->block((blockCount-1)/2 - i);
      N_LAS_Vector& realVecRef_pos =  freqDomainStoreVecReal->block((blockCount-1)/2 + i);
      N_LAS_Vector& imagVecRef_pos =  freqDomainStoreVecImag->block((blockCount-1)/2 + i);

      if (lid >= 0)
      {
        realVecRef_neg[lid] = storeBlock[ 2*(blockCount-i) ];
        imagVecRef_neg[lid] = storeBlock[ 2*(blockCount-i) + 1 ];
        realVecRef_pos[lid] = storeBlock[ 2*i ];
        imagVecRef_pos[lid] = storeBlock[ 2*i+1 ];
      }
    }
  } 

#ifdef Xyce_DEBUG_HB
//  Xyce::dout() << "HB X Vector FD" << std::endl;
  freqDomainSolnVecReal->printPetraObject(std::cout);
  freqDomainSolnVecImag->printPetraObject(std::cout);

  Xyce::dout() << "HB Store Vector FD" << std::endl;

  freqDomainStoreVecReal->printPetraObject(std::cout); 
  freqDomainStoreVecImag->printPetraObject(std::cout);
#endif // Xyce_DEBUG_HB


}


//-----------------------------------------------------------------------------
// Function      : HB::accumulateStatistics()
// Purpose       : Add in the statistics from the current analysis object
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, 1355, Electrical Models & Simulation
// Creation Date : 05/29/13
//-----------------------------------------------------------------------------
void HB::accumulateStatistics_()
{
  hbStatCounts_ += analysisObject_->stats_;
}


//-----------------------------------------------------------------------------
// Function      : HB::runTol_
// Purpose       : Conducts transient run to determine right tolerance
//                 parameters for IC calculation
// Special Notes :
// Scope         : private
// Creator       : T. Mei, 1437, Electrical and Micro Modeling
// Creation Date : 02/23/09
//-----------------------------------------------------------------------------
bool HB::runTol_()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Computing tolerance parameters for HB IC calculation....\n" << std::endl;

  N_TIA_TIAParams & tiaParams = analysisManager_.getTIAParams();
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  // now try the real thing
  tiaParams.initialTime = 0;
  tiaParams.finalTime = period_;
  tiaParams.pauseTime = tiaParams.finalTime;
  tiaParams.resume = false;
  tiaParams.maxOrder = 1;

  // If start up periods are not being run then we can use this transient to compute HB ICs.
  if (!startUpPeriodsGiven_)
  {
    tiaParams.saveTimeStepsFlag = true;
  }

  // register new tiaParams with time integrator
  analysisManager_.registerTIAParams (tiaParams);

  // Create a transient analysis object for this section.
  isTransient_ = true;
  analysisObject_ = Teuchos::rcp(new Transient(analysisManager_));
  analysisObject_->setAnalysisParams( Util::OptionBlock() );
  Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
  returnValue = analysisObject_->run();

  if (!returnValue)
  {
    Report::UserError() << "Calculation of tolerance parameters failed for relErrorTol = " << tiaParams.relErrorTol;
    return false;
  }

  int numPoints = analysisManager_.getStepNumber();

  while((numPoints < (1.2*size_)) && (tiaParams.relErrorTol>= 1e-6))
  {
    {
      Report::UserWarning() << "Tolerance parameters refined, re-running with relErrorTol = " << tiaParams.relErrorTol/10;
    }

    if (!startUpPeriodsGiven_)
    {
      // Clear the fast time data storage before performing the next transient
      N_TIA_DataStore * dsPtr = analysisManager_.getTIADataStore();
      dsPtr->resetFastTimeData();
    }

    tiaParams.relErrorTol =  tiaParams.relErrorTol/10;
    analysisManager_.registerTIAParams (tiaParams);

    // Create a transient analysis object for this section.
    analysisObject_ = Teuchos::rcp(new Transient(analysisManager_));
    analysisObject_->setAnalysisParams( Util::OptionBlock() );
    Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
    bool retV = analysisObject_->run();

    if (!retV)
    {
      Report::UserError() << "Calculation of tolerance parameters failed for relErrorTol = " << tiaParams.relErrorTol;
      return false;
    }
    returnValue = retV && returnValue;

    numPoints = analysisManager_.getStepNumber();
  }

  // Add in simulation times
  accumulateStatistics_();

  // Reset parameters
  tiaParamsSave.relErrorTol  = tiaParams.relErrorTol;
  analysisManager_.registerTIAParams (tiaParamsSave);

  // Reset transient flag before exiting.
  isTransient_ = false;

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::runStartupPeriods_()
// Purpose       : Runs normal transient problem through the requested
//                 number of startup periods
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool HB::runStartupPeriods_()
{
  bool returnValue = true;

  Xyce::lout() << "  ***** Computing " << startUpPeriods_ << " start up periods for HB IC calculation...." << std::endl;

  // need to advance time by startUpPeriods_ * period_
  N_TIA_TIAParams & tiaParams = analysisManager_.getTIAParams();
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  // set DAE initial time = 0.0
  tiaParams.initialTime = 0.0;
  tiaParams.finalTime = startUpPeriods_ * period_;
  tiaParams.pauseTime = tiaParams.finalTime;

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::runStartupPeriods_():  Advancing time through "
               << startUpPeriods_ << " startup periods"
               << " initialTime = " << tiaParams.initialTime
               << " finalTime = " << tiaParams.finalTime << std::endl;

  Xyce::dout() << "HB::runStartupPeriods_():  Double DCOP tiaParams:"
               << " doubleDCOPStep = " << tiaParams.doubleDCOPStep
               << " firstDCOPStep = " << tiaParams.firstDCOPStep
               << " lastDCOPStep = " << tiaParams.lastDCOPStep << std::endl;

  Xyce::dout() << "HB::runStartupPeriods_():  Double DCOP tiaParamsSave:"
               << " doubleDCOPStep = " << tiaParamsSave.doubleDCOPStep
               << " firstDCOPStep = " << tiaParamsSave.firstDCOPStep
               << " lastDCOPStep = " << tiaParamsSave.lastDCOPStep << std::endl;

#endif // Xyce_DEBUG_HB

  {
    Xyce::IO::ActiveOutput x(*analysisManager_.getOutputManager());
    x.add(Xyce::IO::PrintType::HB_STARTUP, ANP_MODE_HB);

    // register new tiaParams with time integrator
    analysisManager_.registerTIAParams (tiaParams);

    // Create a transient analysis object for this section.
    isTransient_ = true;
    analysisObject_ = Teuchos::rcp(new Transient(analysisManager_));
    analysisObject_->setAnalysisParams( Util::OptionBlock() );
    Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
    returnValue = analysisObject_->run();
    isTransient_ = false;

    // Add in simulation times
    accumulateStatistics_();

    analysisManager_.getOutputManagerAdapter().finishOutput();
  }

  // reset the output filename suffix
  // analysisManager_.outMgrPtr->setOutputFilenameSuffix("");

  // put the dsPtr->currentSolutionPtr into dcOpSol and State Vec so that it
  // is used as our initial condition for the pending fast time scale runs
  N_TIA_DataStore * dsPtr = analysisManager_.getTIADataStore();
  dcOpSolVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->currSolutionPtr) ));
  dcOpStateVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->currStatePtr) ));
  dcOpQVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->daeQVectorPtr) ));
  dcOpStoreVecPtr_ = rcp( new N_LAS_Vector( *(dsPtr->currStorePtr) ));

  // tell HB to start after this startup period
  tiaParamsSave.initialTime = startUpPeriods_ * period_;
  transTiaParams_.initialTime = startUpPeriods_ * period_;
  analysisManager_.registerTIAParams (tiaParamsSave);
  startUpPeriodsFinished_ = true;

  return returnValue;
}


//  bool setfreqPoints_();

//-----------------------------------------------------------------------------
// Function      : HB::setFreqPoints_
// Purpose       : Set frequency spectrum for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 03/03/2014
//-----------------------------------------------------------------------------
bool HB::setFreqPoints_()
{
  if ( !intmodMaxGiven_)
  {

    int maxValue = 0;
    if (numFreqs_.size() != 0 )
    // find the max of numFreqs
    {
      maxValue = (numFreqs_[0] - 1)/2;
      for (int i=1; i<numFreqs_.size(); ++i)
      {
        if ((numFreqs_[i] - 1)/2 > maxValue)
          maxValue = (numFreqs_[i] - 1)/2 ;

      }
    } 
    else
    {
      maxValue = (size_ - 1)/2;
    }
    
    intmodMax_ = maxValue;
 
  }


//  std::vector<int> numPosFreqs;

  std::vector<int> k;
 
  int numAnalysisFreqs = transTiaParams_.freqs.size();

  numPosFreqs.resize(numAnalysisFreqs);

  k.resize(numAnalysisFreqs);

  k[0] = 1;

  if (numFreqs_.size() != 0 )
    numPosFreqs[0] = (numFreqs_[0] - 1)/2;     
  else
    numPosFreqs[0]= (size_ - 1)/2;

  int numTotalFrequencies;

  int numExtraFreqs = 0;

  if (numPosFreqs[0] > intmodMax_) 
  {
    numExtraFreqs += ( numPosFreqs[0] - intmodMax_ );
    numFreqs_[0] = (intmodMax_*2 + 1);
  }

  numTotalFrequencies = numFreqs_[0];

  for (int i=1; i < numAnalysisFreqs; i++)
  {
    numPosFreqs[i] = (numFreqs_[i] - 1)/2;
   
    if (numPosFreqs[i] > intmodMax_)
    {
      numExtraFreqs += ( numPosFreqs[i] - intmodMax_ );
      numFreqs_[i] = (intmodMax_*2 + 1);
    }

    k[i] = k[i-1] * numFreqs_[i -1 ];

    numTotalFrequencies *= numFreqs_[i];  
  }

#ifdef Xyce_DEBUG_HB
  for (int i=0; i< numAnalysisFreqs; i++)
  {
    Xyce::dout() << "HB index i" << i << std::endl;
    Xyce::dout() << "HB numPosFreqs =" << numPosFreqs[i] << std::endl;
    Xyce::dout() << "HB k =" << k[i] << std::endl;
  }
  Xyce::dout() << "HB numTotalFrequencies =" << numTotalFrequencies<< std::endl;

  Xyce::dout() << "HB numextrafreqs =" << numExtraFreqs << std::endl;

#endif

  int numIndex = numTotalFrequencies;

  Teuchos::SerialDenseMatrix<int,double> indexMatrix(numAnalysisFreqs, numTotalFrequencies);
//  Teuchos::SerialDenseMatrix<int,int>  indexMatrix(numAnalysisFreqs, numTotalFrequencies);


#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB intmodMax =" <<  intmodMax_ << std::endl;  
#endif

  int nextIndex; 

  int idxMod,  idxValues;
  int sumIndex;

  std::vector<int> goodIndex;
  for (int i=0; i < numIndex; i++)      //  column
  {
    nextIndex = i;
    sumIndex = 0;

    for (int j= (numAnalysisFreqs - 1); j >= 0;  j-- )       // row 
    {
      idxMod = nextIndex%k[j];
      idxValues =  (nextIndex - idxMod)/k[j];

      indexMatrix (j, i) = static_cast<double>(idxValues - (numFreqs_[j] - 1)/2 );
      nextIndex = idxMod;
      sumIndex += abs(idxValues - (numFreqs_[j] - 1)/2 );

    }

    if( sumIndex <= intmodMax_) 
      goodIndex.push_back(i);

  }

  int diaindexSize = goodIndex.size();
 
  Teuchos::SerialDenseMatrix<int,double> diaindexMatrix( numAnalysisFreqs, (diaindexSize + numExtraFreqs) );
  diaindexMatrix.putScalar(0.0); 

  for (int i=0; i < diaindexSize; i++)       
  {
    for (int j= (numAnalysisFreqs - 1); j >= 0;  j-- )
      diaindexMatrix (j, i) = indexMatrix (j, goodIndex[i]);
  }

#ifdef Xyce_DEBUG_HB 
  for (int i=0; i< diaindexSize; i++)
  {
    std::cout << "good index i = " << i << goodIndex[i]  << std::endl;
  }


  std::cout <<  " checking diamond indexMatrix" << std::endl;
  diaindexMatrix.print(std::cout);

  std::cout <<  " checking indexMatrix" << std::endl;
  indexMatrix.print(std::cout);

#endif 
  
  int extraIndexPos = diaindexSize;
  for (int i=0; i < numAnalysisFreqs ; i++)      //  column
  {

    if (numPosFreqs[i] > intmodMax_)
    {

      for (int j=0; j < (numPosFreqs[i] - intmodMax_ ) ; j++)
        diaindexMatrix (i, extraIndexPos + j ) = static_cast<double>(intmodMax_ + j + 1 );

      extraIndexPos += (numPosFreqs[i] - intmodMax_ );

    }  

  }


#ifdef Xyce_DEBUG_HB 
  std::cout <<  " checking diamond indexMatrix after axis" << std::endl;
  diaindexMatrix.print(std::cout);
#endif 
 
// get the positive frequencies
  
  std::vector<double> posfreqPoints_;     
   
  int posindexSize = (diaindexSize - 1)/2;
  posfreqPoints_.resize(posindexSize + numExtraFreqs); 

  Teuchos::SerialDenseMatrix<int,double> currindexMatrix( Teuchos::View, diaindexMatrix, numAnalysisFreqs, (posindexSize + numExtraFreqs), 0, posindexSize+1 );
  Teuchos::SerialDenseVector<int,double> currfreqPoints( Teuchos::View, &posfreqPoints_[0], (posindexSize + numExtraFreqs ) );

  Teuchos::SerialDenseVector<int,double> hbFreqs( Teuchos::View, &transTiaParams_.freqs[0], numAnalysisFreqs);
//    Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i], oversampleRate*size_-(i+1) );
  currfreqPoints.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, currindexMatrix, hbFreqs, 0.0 );


#ifdef Xyce_DEBUG_HB 
  std::cout << "checking positive frequencies" << std::endl;
  currindexMatrix.print(std::cout);
  hbFreqs.print(std::cout); 
  currfreqPoints.print(std::cout);
#endif 

  for (int i=0; i < posindexSize; i++)       
  {
 
    if (posfreqPoints_[i] < 0.0)
      posfreqPoints_[i] = fabs( posfreqPoints_[i]);
      
  }
  
//  size_ = (posindexSize + numExtraFreqs ) *2 + 1;

  std::sort(posfreqPoints_.begin(), posfreqPoints_.end() );


#ifdef Xyce_DEBUG_HB
  for (int i=0; i< posfreqPoints_.size(); i++)
  {
    std::cout << "pos frequency point " <<  posfreqPoints_[i] << std::endl;
  }
#endif
  posfreqPoints_.erase(std::unique(posfreqPoints_.begin(), posfreqPoints_.end() ), posfreqPoints_.end() );


  if (abs( posfreqPoints_[0]) < 2.0*Util::MachineDependentParams::MachinePrecision() )
    posfreqPoints_.erase( posfreqPoints_.begin()); 

  size_ = ( posfreqPoints_.size() ) *2 + 1;
#ifdef Xyce_DEBUG_HB
  for (int i=0; i<posfreqPoints_.size(); i++)
  {
    std::cout << "pos frequency point after " <<  posfreqPoints_[i] << std::endl;
  }

  Xyce::dout() << "HB size =" << size_ << std::endl;
#endif
  
  int i=0;
  freqPoints_.resize(size_);

  for( i = 0; i < size_; ++i )
  {
    if (i < (size_-1)/2)
      freqPoints_[i] = - posfreqPoints_[ (size_-1)/2 - i - 1 ];
    else if (i > (size_-1)/2)
      freqPoints_[i] =  posfreqPoints_[ i - (size_-1)/2 - 1 ]; 
    else
      freqPoints_[i] = 0.0;
  } 

#ifdef Xyce_DEBUG_HB
  for (int i=0; i<freqPoints_.size(); i++)
  {
    std::cout << " frequency point " << freqPoints_[i] << std::endl;
  }
#endif


  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::setTimePoints_
// Purpose       : Set time points for multi-tone HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 03/05/2014
//-----------------------------------------------------------------------------
bool HB::setTimePoints_()
{
  // NOTE:  Need to make this parallel safe.

  int posFreq = (size_-1)/2;
  int oversampleRate = 2; 
  int periodSampleMultiplier = 1;
//  int periodSampleMultiplier = 1;

  Teuchos::BLAS<int,double> blas;
  std::vector<double> testPoints(oversampleRate*size_);

//  std::cout << "Period is" << period_ << std::endl;
  for (int i=0; i<oversampleRate*size_; ++i)
  {
    testPoints[i] = periodSampleMultiplier*period_*((Teuchos::ScalarTraits<double>::random()+1)/2);
 //    testPoints[i] = periodSampleMultiplier*period_/(oversampleRate*size_)* static_cast<double>(i);
//    std::cout << "testPoints" << i << "=" << testPoints[i]  << std::endl;
  }

  Teuchos::SerialDenseMatrix<int,double> testMatrix(size_,oversampleRate*size_);
  // NOTE: i represents frequency, j represents sample time
  // Set DC values first.
  for (int j=0; j<oversampleRate*size_; ++j)
  {
    testMatrix(0,j) = 1.0;
  }
  // Set rest of frequency values
  for (int i=1; i<=posFreq; i++)
  {
//    for (int j=0; j<oversampleRate*size_; j+2)
    for (int j=0; j<oversampleRate*size_; j++)
    {
      testMatrix(2*i-1,j) = cos(2*M_PI*freqPoints_[posFreq+i]*testPoints[j]);
      testMatrix(2*i,j) = sin(2*M_PI*freqPoints_[posFreq+i]*testPoints[j]);
    }
  }
  // Check matrix is correct
//  std::cout << "Checking testMatrix" << std::endl;
//  testMatrix.print(std::cout);

  // Now for orthogonalization (yipeee)
  std::vector<double> weightVector(oversampleRate*size_);
  for (int i=0; i<size_; ++i)
  {
    // Find column with largest norm, choose as next vector for orthogonalization
    int maxIndex = 0;
    double maxValue = 0.0;
    for (int j=i; j<oversampleRate*size_; ++j)
    {
      Teuchos::SerialDenseMatrix<int,double> tempVector( Teuchos::View, testMatrix, size_, 1, 0, j );
      weightVector[j] = tempVector.normFrobenius();
    
      if (weightVector[j] > maxValue)
      {
        maxIndex = j;
        maxValue = weightVector[j];
      }
    }  
   

    // Swap time and vector.
    std::swap( testPoints[i], testPoints[maxIndex] );
    Teuchos::SerialDenseVector<int,double> newSwapVector2 = Teuchos::getCol<int,double>( Teuchos::Copy, testMatrix, maxIndex );
    Teuchos::SerialDenseVector<int,double> newSwapVector = Teuchos::getCol<int,double>( Teuchos::View, testMatrix, i );
    Teuchos::setCol<int,double>( newSwapVector, maxIndex, testMatrix );
    Teuchos::setCol<int,double>( newSwapVector2, i, testMatrix );

//    std::cout << "Checking testMatrix for i= " << i << std::endl;
//    testMatrix.print(std::cout);   
 
    // Compute inner product with vector from last time point with rest of time points
    Teuchos::SerialDenseMatrix<int,double> currTestMatrix( Teuchos::View, testMatrix, size_, oversampleRate*size_-(i+1), 0, i+1 );
    Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i+1], oversampleRate*size_-(i+1) );
//    Teuchos::SerialDenseVector<int,double> currWeightVector( Teuchos::View, &weightVector[i], oversampleRate*size_-(i+1) );
    currWeightVector.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, currTestMatrix, newSwapVector, 0.0 );

//    std::cout << "The current norm is" << std::endl;
//    currWeightVector.print(std::cout);

    // Subtract off scaled vector from rest of time points.
    // NOTE:  maxValue is the norm of the last time point.
    for (int j=i+1; j<oversampleRate*size_; ++j)
    {
       Teuchos::SerialDenseMatrix<int,double> currVector( Teuchos::View, testMatrix, size_, 1, 0, j );
       blas.AXPY( size_, -(currWeightVector[j-(i+1)]/(maxValue*maxValue) ), newSwapVector.values(), 1, currVector.values(), 1 );
    }

//    std::cout << "Checking testMatrix after orthogonalization  for i = " << i << std::endl;
//    testMatrix.print(std::cout);

  }
//  std::cout << "Checking testMatrix after orthogonalization." << std::endl;
//  testMatrix.print(std::cout);

  // Sort the chosen test points and then copy them into the fastTimes_ vector.
  std::sort( testPoints.begin(), testPoints.begin()+size_ );

  fastTimes_.resize (size_); 

  for (int i=0; i<size_; ++i)
  {
    fastTimes_[i] = testPoints[i];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HB::createFT_
// Purpose       : Create the DFT and IFT matrices using the time points for multi-tone HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 03/05/2014
//-----------------------------------------------------------------------------
bool HB::createFT_()
{
  int posFreq = (size_-1)/2;
  idftMatrix_.reshape(size_,size_);

  // NOTE: i represents frequency, j represents sample time
  // Set DC values first.
  for (int i=0; i<size_; ++i)
  {
    idftMatrix_(i,0) = 1.0;
  }
  // Set rest of frequency values
  for (int i=0; i<size_; i++)
  {
    for (int j=1; j<=posFreq; j++)
    {
      idftMatrix_(i,2*j-1) = cos(2*M_PI*freqPoints_[posFreq+j]*fastTimes_[i]);
      idftMatrix_(i,2*j) = sin(2*M_PI*freqPoints_[posFreq+j]*fastTimes_[i]);
    }
  }  

  std::cout << "Checking IDFTmatrix" << std::endl;
  idftMatrix_.print(std::cout);

  // Compute DFT matrix.
  dftMatrix_ = idftMatrix_;
  Teuchos::SerialDenseSolver<int,double> ftSolver;
  ftSolver.setMatrix( Teuchos::rcp( &dftMatrix_, false ) );
  ftSolver.invert();

  std::cout << "Checking DFTmatrix" << std::endl;
  dftMatrix_.print(std::cout);

  std::cout << "Checking IDFTmatrix after inverse:" << std::endl;
  idftMatrix_.print(std::cout);

  return true;

}

//-----------------------------------------------------------------------------
// Function      : HB::setInitialGuess_
// Purpose       : Set initial guess for HB analysis.
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 03/03/2014
//-----------------------------------------------------------------------------
bool HB::setInitialGuess_()
{

  bool returnValue = true;  

  if (taHB_ == 1)
  {
  bool retTol1 = runTol_(); returnValue = returnValue && retTol1;

  // Start up periods need to be run before the initial condition is computed, otherwise
  // just used the solution from the tolerance calculation.
  if( startUpPeriodsGiven_ )
  {
    bool startupPeriodsSuccess = runStartupPeriods_();
    if (!startupPeriodsSuccess)
    {
      Report::UserError() << "Failed to calculate the startup periods";
      return false;
    }
    returnValue = returnValue && startupPeriodsSuccess;

    bool icSuccess = runTransientIC_();
    if (!icSuccess)
    {
      Report::UserError() << "Initial HB Transient failed";
      return false;
    }
    returnValue = returnValue && icSuccess;
  }

  interpolateIC_();

  }
  else
  {

    double TimeStep = period_/static_cast<double>(size_);
    timeSteps_.push_back( TimeStep );
    fastTimes_.resize(size_+1);

    goodTimePoints_.resize(size_+1);
    for( int i = 0; i <= size_; ++i )
    {
      fastTimes_[i] = transTiaParams_.initialTime + static_cast<double>(i) * TimeStep;
    }
    goodTimePoints_ = fastTimes_;

  }

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::runTransientIC_
// Purpose       : Conducts a regular transient run for HB initial conditions
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 10/03/2008
//-----------------------------------------------------------------------------
bool HB::runTransientIC_()
{
  bool returnValue = true;

  Xyce::lout() << " ***** Running transient to compute HB initial condition....\n" << std::endl;

  // this prevents extra DC op data from being printed.
  devInterfacePtr_->setMPDEFlag( true );

  if(saveIcData_)
  {
    // Keep the initial condition data
    // analysisManager_.outMgrPtr->setOutputFilenameSuffix( ".hb_ic" );
  }

  // use an initial transient run to create a set of time points for the fast time scale
  N_TIA_TIAParams & tiaParams = analysisManager_.getTIAParams();
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  //tiaParams.initialTime = 0; // should be tiaParams_.initialTime;
  tiaParams.initialTime = transTiaParams_.initialTime;
  //tiaParams.finalTime = period_; // should be tiaParams_.initialTime + period_;
  tiaParams.finalTime = transTiaParams_.initialTime + period_;
  tiaParams.pauseTime = tiaParams.finalTime;
  tiaParams.saveTimeStepsFlag = true;
  tiaParams.maxOrder = 1;

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::runTransientIC_():  Advancing time from"
       << " initialTime = " << tiaParams.initialTime
       << " finalTime = " << tiaParams.finalTime << std::endl;
#endif // Xyce_DEBUG_HB

  // Initial conditions will be set if startup periods were run.
  if ( startUpPeriodsGiven_ )
  {
    tiaParams.NOOP = true;

    N_TIA_DataStore * dsPtr = analysisManager_.getTIADataStore();
    *(dsPtr->nextSolutionPtr) = *(dcOpSolVecPtr_.get());
    *(dsPtr->nextStatePtr) = *(dcOpStateVecPtr_.get());
    *(dsPtr->daeQVectorPtr) = *(dcOpQVecPtr_.get());
    *(dsPtr->nextStorePtr) = *(dcOpStoreVecPtr_.get());
  }

  // register new tiaParams with time integrator
  analysisManager_.registerTIAParams (tiaParams);

  // Create a transient analysis object for this section.
  isTransient_ = true;
  analysisObject_ = Teuchos::rcp(new Transient(analysisManager_));
  analysisObject_->setAnalysisParams( Util::OptionBlock() );
  Teuchos::rcp_dynamic_cast<Transient>(analysisObject_)->resetForHB();
  returnValue = analysisObject_->run();
  isTransient_ = false;

  // Add in simulation times
  accumulateStatistics_();

  if(saveIcData_)
  {
    // reset suffix
    // analysisManager_.outMgrPtr->setOutputFilenameSuffix( "" );
  }

  tiaParamsSave.initialTime += period_;  // start HB problem after this transient init.
  analysisManager_.registerTIAParams (tiaParamsSave);
  devInterfacePtr_->setMPDEFlag( false );

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : HB::interpolateIC_()
// Purpose       : Tries to filter the fast time points from a transient run
//                 so that points are not too close together
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool HB::interpolateIC_()
{
  Xyce::lout() << " ***** Interpolating transient solution for IC calculation....\n" << std::endl;

  N_TIA_DataStore *dsPtr = analysisManager_.getTIADataStore();
  int numPoints = dsPtr->timeSteps.size();

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "HB::interpolateIC_(): Initial transient run produced " << numPoints << " points." << std::endl;
#endif

  std::vector<int> goodIndicies;
  goodTimePoints_.resize(size_);

  double TimeStep = period_/static_cast<double>(size_);
  timeSteps_.push_back( TimeStep );
  for( int i = 0; i < size_; ++i )
  {
    goodTimePoints_[i] = transTiaParams_.initialTime + static_cast<double>(i) * TimeStep;
  }
  fastTimes_ = goodTimePoints_;
  fastTimes_.resize(size_+1);
  fastTimes_[size_] = transTiaParams_.initialTime + period_;

  int breakpoints = 0;          // need to keep track of how many breakpoints there are
  int startIndex = 0;

  // always keep first point
  goodIndicies.push_back(startIndex);
  int GoodTimePointIndex = startIndex + 1;

  for( int i=startIndex; i < numPoints - 1 ; i++ )
  {
    // count up breakpoints
    if( dsPtr->timeStepsBreakpointFlag[i] == true )
    {
      breakpoints++;
    }

#ifdef Xyce_DEBUG_HB
    if( debugLevel > 0 )
    {
      Xyce::dout() << "\t\t timeStep[ " << i << " ] = " << dsPtr->timeSteps[i];
      if( dsPtr->timeStepsBreakpointFlag[i] == true )
      {
        Xyce::dout() << "  Breakpoint";
      }
      Xyce::dout() << std::endl;
    }
#endif
    while( ( GoodTimePointIndex < size_ )  && (dsPtr->timeSteps[i] <= goodTimePoints_[GoodTimePointIndex]) && (goodTimePoints_[GoodTimePointIndex] < dsPtr->timeSteps[i+1]))
    {
      // found a good point so save the index
      goodIndicies.push_back( i );
      GoodTimePointIndex =  GoodTimePointIndex+1;
    }
  }

  for(int i=0; i<size_; i++ )
  {
    int currentIndex = goodIndicies[i];
    N_LAS_Vector * firstSolVecPtr = dsPtr->fastTimeSolutionVec[currentIndex];
    N_LAS_Vector * secondSolVecPtr = dsPtr->fastTimeSolutionVec[currentIndex+1];

    N_LAS_Vector * firstStateVecPtr = dsPtr->fastTimeStateVec[currentIndex];
    N_LAS_Vector * secondStateVecPtr = dsPtr->fastTimeStateVec[currentIndex+1];

    N_LAS_Vector * firstQVecPtr = dsPtr->fastTimeQVec[currentIndex];
    N_LAS_Vector * secondQVecPtr = dsPtr->fastTimeQVec[currentIndex+1];

    N_LAS_Vector * firstStoreVecPtr = dsPtr->fastTimeStoreVec[currentIndex];
    N_LAS_Vector * secondStoreVecPtr = dsPtr->fastTimeStoreVec[currentIndex+1];

    double fraction = (goodTimePoints_[i] -  dsPtr->timeSteps[currentIndex])/(dsPtr->timeSteps[currentIndex+1] -  dsPtr->timeSteps[currentIndex]);

    RCP<N_LAS_Vector> InterpICSolVecPtr = rcp( new N_LAS_Vector( *secondSolVecPtr ) );
    RCP<N_LAS_Vector> InterpICStateVecPtr = rcp( new N_LAS_Vector( *secondStateVecPtr ) );
    RCP<N_LAS_Vector> InterpICQVecPtr = rcp( new N_LAS_Vector( *secondQVecPtr ) );
    RCP<N_LAS_Vector> InterpICStoreVecPtr = rcp( new N_LAS_Vector( *secondStoreVecPtr ) );

    InterpICSolVecPtr->putScalar(0.0);
    InterpICStateVecPtr->putScalar(0.0);
    InterpICQVecPtr->putScalar(0.0);
    InterpICStoreVecPtr->putScalar(0.0);

    InterpICSolVecPtr->linearCombo(-1.0, *firstSolVecPtr, 1.0, *secondSolVecPtr );
    InterpICSolVecPtr->linearCombo(1.0, *firstSolVecPtr, fraction , *InterpICSolVecPtr);

    InterpICStateVecPtr->linearCombo(-1.0, *firstStateVecPtr, 1.0, *secondStateVecPtr );
    InterpICStateVecPtr->linearCombo(1.0, *firstStateVecPtr, fraction , *InterpICStateVecPtr);

    InterpICQVecPtr->linearCombo(-1.0, *firstQVecPtr, 1.0, *secondQVecPtr );
    InterpICQVecPtr->linearCombo(1.0, *firstQVecPtr, fraction , *InterpICQVecPtr);

    InterpICStoreVecPtr->linearCombo(-1.0, *firstStoreVecPtr, 1.0, *secondStoreVecPtr );
    InterpICStoreVecPtr->linearCombo(1.0, *firstStoreVecPtr, fraction , *InterpICStoreVecPtr);

    goodSolutionVec_.push_back(InterpICSolVecPtr);
    goodStateVec_.push_back(InterpICStateVecPtr);
    goodQVec_.push_back(InterpICQVecPtr);
    goodStoreVec_.push_back(InterpICStoreVecPtr);
  }

  // Clean up the fast time data since we are finished computing the initial condition.
  // The fast time data can take a considerable amount of memory for large problems.
  dsPtr->resetFastTimeData();

  return true;
}

} // namespace Analysis
} // namespace Xyce

