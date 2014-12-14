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
// Filename      : $RCSfile: N_MPDE_Manager.h,v $
//
// Purpose       : This file defines the manager class for the MPDE package
//
// Special Notes :
//
// Creator       : Robert Hoekstra, 9233, Computational Sciences
//
// Creation Date : 3/11/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.90 $
//
// Revision Date  : $Date: 2014/08/08 20:08:07 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_MANAGER_H
#define Xyce_MPDE_MANAGER_H

#include <string>
#include <map>

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_IO_fwd.h>
#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_TOP_fwd.h>

#include <N_LAS_BlockVector.h>

#include <N_IO_PkgOptionsMgr.h>

#include <N_MPDE_State.h>
#include <N_TIA_TIAParams.h>
#include <N_MPDE_WarpedPhaseCondition.h>
#include <N_MPDE_DeviceInterface.h>
#include <N_LAS_System.h>
#include <N_NLS_Manager.h>

#include <Teuchos_RCP.hpp>

class N_MPDE_Loader;
class N_MPDE_Builder;
class N_MPDE_Discretization;

class N_TIA_MPDEInterface;
class N_PDS_Manager;

class N_LOA_Loader;
class N_LOA_NonlinearEquationLoader;
class N_LAS_Builder;
// class N_LAS_System;

class N_LAS_PrecondFactory;

//-----------------------------------------------------------------------------
// Class         : N_MPDE_Manager
// Purpose       : MPDE Manager Class
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
class N_MPDE_Manager
{
 public:

  enum MPDE_IC
  {
    MPDE_IC_DCOP,         // 0  -- use DC sweep to approx IC
    MPDE_IC_SAWTOOTH,     // 1  -- use many small transients to approx IC.
    WAMPDE_IC_TRANSIENT,  // 2  -- use a transient sim as an approx IC.
                          //       New algorithms are added to get better IC for 
                          //       complex WaMPDE circuits. genie 042914.
    MPDE_IC_TRAN_MAP,     // 3  -- use many transient sims map out an IC.
    MPDE_IC_TWO_PERIOD,    // 4  -- use two periods to interpolate IC.
    MPDE_IC_TRANSIENT    // 5  -- use a transient sim as an approx IC. 
                          //       This is the old ic=2. genie 042914
  };


  // Default constructor
  N_MPDE_Manager(N_IO_CmdParse & cp);

  // Destructor
  ~N_MPDE_Manager()
  {}

  //Runs MPDE analysis
  bool run();

  // MPDE analysis flag
  bool blockAnalysisFlag() const { return blockAnalysisFlag_; }

  //Registrations
  void registerAnalysisManager( N_ANP_AnalysisManager *anaIntPtr );

  void registerTIAMPDEInterface( Teuchos::RCP<N_TIA_MPDEInterface> tiaMPDEIfacePtr );

  void registerNonlinearSolver( Teuchos::RCP<N_NLS_Manager> nlsMgrPtr );

  void registerDeviceInterface ( N_DEV_DeviceInterface *devInterfacePtr );

  void registerParallelManager( N_PDS_Manager *pdsMgrPtr );

  void registerTopology( N_TOP_Topology *topoMgrPtr );

  void registerRestartManager( N_IO_RestartMgr *resMgrPtr );

  void registerOutputManager( N_IO_OutputMgr *outMgrPtr );

  //registrations of application system builder and loader
  void registerApplicationLoader( N_LOA_Loader *appLoaderPtr );

  void registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader *appLoaderPtr );

  void registerApplicationBuilder( N_LAS_Builder *appBuilderPtr );

  void registerLinearSystem(N_LAS_System *lasSysPtr);

  // Method to register the utility options.
  bool setMPDEAnalysisParams(const N_UTL_OptionBlock & OB);

  // Method to register the utility options.
  bool setMPDEOptions(const N_UTL_OptionBlock & OB);

  // Method to register the utility options.
  bool registerTranMPDEOptions(const N_UTL_OptionBlock & OB);

  // "get" function for the fast time scale.
  //double getFastTime ();

  // "get" function for MPDE flag. (true if not IC)
  bool getMPDEFlag ();

  bool getMPDEStartupFlag();
  // get function for MPDE initial condition flag (true if MPDE & IC)
  bool getMPDEIcFlag();

  // "get" function for WaMPDE flag. (true if not IC)
  bool getWaMPDEFlag ();

  const std::vector<double> & getFastTimePoints () const;

    const std::vector<double> & getFreqPoints () const;

  // "get" function for GID of phi variable in Warped MPDE case
  int getPhiGID ();

 private :
  void printParams_ ();

 public:
  int debugLevel;

 private:

  //Run Stages
  bool initializeAll_();
  bool initializeMPDEAll_(int size); //genie 111413
  bool initializeOscOut_();
  bool runInitialCondition_();
  bool runDCOP_();
  bool runStartupPeriods_();
  bool runTransientIC_();
  bool filterFastTimePoints_();
  bool setupMPDEProblem_();
  bool runMPDEProblem_();
  bool runTests_();

  void setMPDEOnFlag( bool flagVal );
  void setRunFlag( bool flagVal );
  void setBlockAnalysisFlag_( bool flagVal );

  // diagnostic function
  double checkPeriodicity_();

  N_IO_CmdParse & commandLine_;

  N_ANP_AnalysisManager *               anaIntPtr_;

  Teuchos::RCP<N_TIA_MPDEInterface> tiaMPDEIfacePtr_;
  N_TIA_TIAParams tiaMPDEParams_;
  N_UTL_OptionBlock saved_tranMPDEOB_;

  N_NLS_Manager                         nlsMgrPtr_;
  N_LAS_System                          lasMPDESysPtr_;
  N_DEV_DeviceInterface *               devInterfacePtr_;
  Teuchos::RCP< N_MPDE_DeviceInterface > mpdeDeviceInterfacePtr_;
  N_PDS_Manager *                       pdsMgrPtr_;
  N_TOP_Topology *                      topoMgrPtr_;
  N_IO_RestartMgr *                     resMgrPtr_;
  N_IO_OutputMgr *                      outMgrPtr_;
  N_LOA_Loader *                        appLoaderPtr_;
  N_LOA_NonlinearEquationLoader *       nonlinearEquationLoaderPtr_;
  N_LAS_Builder *                       appBuilderPtr_;
  N_LAS_System *                        lasSysPtr_;
//  Teuchos::RCP<N_LAS_System> lasMPDESysPtr_;

  N_MPDE_State mpdeState_;
  Teuchos::RCP<N_MPDE_Loader> mpdeLoaderPtr_;
  Teuchos::RCP<N_MPDE_Builder> mpdeBuilderPtr_;
  Teuchos::RCP<N_MPDE_Discretization> mpdeDiscPtr_;

  Teuchos::RCP<N_LAS_Vector> dcOpSolVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpStateVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpQVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpStoreVecPtr_;

  Teuchos::RCP<N_LAS_Vector> endIcSolVecPtr_;
  Teuchos::RCP<N_LAS_Vector> endIcStateVecPtr_;
  Teuchos::RCP<N_LAS_Vector> endIcQVecPtr_;
  Teuchos::RCP<N_LAS_Vector> endIcStoreVecPtr_;

  // Do we have a .MPDE analysis line
  bool blockAnalysisFlag_;

  //Testing Flag
  bool test_;

  //MPDE Problem Size Factor
  int size_;
  bool tranRunForSize_;      // use an initial transient run to calculate size_
  int maxCalcSize_;          // max size to use from transient run.
  bool maxCalcSizeGiven_;

  //MPDE Fast Src driving the fast time scale oscillations
  std::string fastSrc_;
  bool fastSrcGiven_;
  std::vector<std::string> srcVec_;

  // Independent variable for warped MPDE.
  std::string oscOut_;
  bool oscOutGiven_;

  // Number of steps to take during the start of the full MPDE
  // calculation before turning on truncation errror control
  int nonLteSteps_;
  bool nonLteStepsGiven_;

  //MPDE Fast Time Scale Period
  double period_;
  bool periodGiven_;

  //MPDE number of fast time periods to integrate over and IGNORE before
  // getting initial conditions for MPDE.  Default is zero.
  int startUpPeriods_;
  bool startUpPeriodsGiven_;
  bool startUpPeriodsFinished_;
  bool saveIcData_;

  // MPDE initial condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICVectorPtr_;

  // MPDE initial state condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICStateVectorPtr_;

  // MPDE initial Q condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICQVectorPtr_;

  // MPDE initial store condition
  Teuchos::RCP<N_LAS_BlockVector> mpdeICStoreVectorPtr_;

  //MPDE fast time points
  // 12/5/06 tscoffe:  Note, the period T2 is the last element in fastTimes.
  // This means that the number of fast time points is fastTimes_.size()-1
  std::vector<double> fastTimes_;

  std::vector<double> freqPoints_;

  // Fast time discretization
  int fastTimeDisc_;
  int fastTimeDiscOrder_;

  //if we pull data directly from an initial transient run, keep
  // a list of the indices we used so that we can pull out the solution
  // and state data too.
  std::vector<int> indicesUsed_;
  std::vector<bool> nonPeriodicFlags;
  // warped MPDE setting in netlist
  bool warpMPDE_;

  // WaMPDE phase condition class
  Teuchos::RCP<N_MPDE_WarpedPhaseCondition> warpMPDEPhasePtr_;

  // WaMPDE OSCOUT GID:
  int warpMPDEOSCOUT_;

  // WaMPDE Phase equation:
  int warpPhase_;
  bool warpPhaseGiven_;

  // WaMPDE Phase equation constant:
  double warpPhaseCoeff_;
  bool warpPhaseCoeffGiven_;

  // frequency domain flag
  bool fftFlag_;

  // Number of fast periods used for initial condition.
  int icPer_;

  // Initial condition strategy.
  int initialCondition_;

  // MPDE mode flag.  if false, run initial condition, if true run MPDE.
  bool mpdeOnFlag_;

  // MPDE IC flag; if true then we're calculating ic's for MPDE
  bool mpdeICFlag_;

  // myhsieh 081213
  // WaMPDE initial transient IC flag; if true then the initial condition for the 
  // initial fast transient run is given (e.g. .IC is pecified in the netlist to
  // kickoff the oscillation). This parameter is specifically for WaMPDE since
  // MPDE circuits have oscillations coming from OSCSRC.
  bool warpMPDEICFlag_;


  // debug flags:
  bool dcopExitFlag_;
  bool icExitFlag_;
  int  exitSawtoothStep_;

  // An analysis-dependent preconditioner factory.
  Teuchos::RCP<N_LAS_PrecondFactory> precFactory_;

};

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getMPDEFlag()
{
  return mpdeOnFlag_;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getMPDEStartupFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getMPDEStartupFlag()
{
    return startUpPeriodsFinished_;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getMPDEIcFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437, Elec. Modeling
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getMPDEIcFlag()
{
  return mpdeICFlag_;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::getWaMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline bool N_MPDE_Manager::getWaMPDEFlag()
{
  return warpMPDE_;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerAnalysisManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerAnalysisManager( N_ANP_AnalysisManager *anaIntPtr )
{
  anaIntPtr_ = anaIntPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerTIAMPDEInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerTIAMPDEInterface( Teuchos::RCP<N_TIA_MPDEInterface> tiaMPDEIfacePtr )
{
  tiaMPDEIfacePtr_ = tiaMPDEIfacePtr;
}

// //-----------------------------------------------------------------------------
// // Function      : N_MPDE_Manager::registerNonlinearSolver
// // Purpose       :
// // Special Notes :
// // Scope         : public
// // Creator       : Robert Hoekstra, 9233, Computational Sciences
// // Creation Date : 03/12/04
// //-----------------------------------------------------------------------------
// inline void N_MPDE_Manager::registerNonlinearSolver( N_NLS_Manager *nlsMgrPtr )
// {
//   nlsMgrPtr_ = nlsMgrPtr;
// }

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerDeviceInterface( N_DEV_DeviceInterface *devInterfacePtr )
{
  devInterfacePtr_ = devInterfacePtr;
  mpdeDeviceInterfacePtr_ = Teuchos::rcp( new N_MPDE_DeviceInterface() );
  mpdeDeviceInterfacePtr_->registerDeviceInterface( devInterfacePtr );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerParallelManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerParallelManager(N_PDS_Manager *pdsMgrPtr )
{
  pdsMgrPtr_ = pdsMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerTopology(N_TOP_Topology *topoMgrPtr )
{
  topoMgrPtr_ = topoMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerRestartManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerRestartManager( N_IO_RestartMgr *resMgrPtr )
{
  resMgrPtr_ = resMgrPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerOutputManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerOutputManager( N_IO_OutputMgr *outMgrPtr )
{
  outMgrPtr_ = outMgrPtr;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerApplicationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerApplicationLoader( N_LOA_Loader *appLoaderPtr )
{
  appLoaderPtr_ = appLoaderPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/06/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader *nonlinearEquationLoaderPtr )
{
  nonlinearEquationLoaderPtr_ = nonlinearEquationLoaderPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerApplicationBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerApplicationBuilder( N_LAS_Builder *appBuilderPtr )
{
  appBuilderPtr_ = appBuilderPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::registerLinearSystem( N_LAS_System *lasSysPtr )
{
  lasSysPtr_ = lasSysPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getFreqTimePoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/02/08
//-----------------------------------------------------------------------------
inline const std::vector<double> & N_MPDE_Manager::getFreqPoints() const
{
  return freqPoints_;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getFastTimePoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 04/06/04
//-----------------------------------------------------------------------------
inline const std::vector<double> & N_MPDE_Manager::getFastTimePoints() const
{
  return fastTimes_;
}



//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::getPhiGID
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey 1414
// Creation Date : 01/10/07
//-----------------------------------------------------------------------------
inline int N_MPDE_Manager::getPhiGID()
{
  return warpMPDEPhasePtr_->getPhiGID();
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setMPDEOnFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey 1414
// Creation Date : 01/10/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::setMPDEOnFlag( bool flagVal )
{
  mpdeOnFlag_ = flagVal;
  mpdeDeviceInterfacePtr_->setMPDEFlag( flagVal );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey 1414
// Creation Date : 01/10/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Manager::setBlockAnalysisFlag_( bool flagVal )
{
  blockAnalysisFlag_ = flagVal;
  mpdeDeviceInterfacePtr_->setBlockAnalysisFlag( flagVal );
}


#endif //Xyce_MPDE_MANAGER_H
