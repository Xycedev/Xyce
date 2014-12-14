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
// Filename      : $RCSfile: N_ANP_AnalysisManager.h,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.78.2.4 $
//
// Revision Date  : $Date: 2014/08/28 21:00:43 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AnalysisManager_h
#define Xyce_N_ANP_AnalysisManager_h

#include <list>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_UTL_Xyce.h>
#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_LOA_fwd.h>
#include <N_LAS_fwd.h>
#include <N_NLS_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_StepEvent.h>
#include <N_ANP_AnalysisEvent.h>
#include <N_NLS_Manager.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_Listener.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Stats.h>
#include <N_UTL_Timer.h>

class N_MPDE_Manager;

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Class         :
// Purpose       : This function converts between Nonlinear::AnalysisMode
//               : and AnalysisManager.h ANP_Analysis_Mode
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
Nonlinear::AnalysisMode anpAnalysisModeToNLS(Analysis_Mode mode);

const char *analysisModeName(Analysis_Mode mode);

//-----------------------------------------------------------------------------
// Class         : AnalysisManager
//
// Purpose       : This class manages, allocates, and sets up the
//                 various analysis types, such as DC, Trn,a HB, etc.
//
// Special Notes : Some of this class was once in the N_TIA_ControlAlgorithm
//                 class, which was set up back when Xyce only did transient
//                 simulations.  As we added analysis types it became necessary
//                 to refactor the code so that each analysis type had its own
//                 set of classes.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 6/01/00 (N_TIA_ControlAlgorithm, now deprecated)
// Creation Date : 1/24/08 (for this version of the class. Date is approximate)
//-----------------------------------------------------------------------------
class AnalysisManager : public Util::Notifier<StepEvent>,
                        public Util::Notifier<AnalysisEvent>,
                        public Util::ListenerAutoSubscribe<StepEvent>,
                        public Util::ListenerAutoSubscribe<AnalysisEvent>
{
  public:
    // Default constructor.
    AnalysisManager(IO::CmdParse & cp, Stats::Stat root_stat);

    // Destructor
    ~AnalysisManager();

    void notify(const StepEvent &step_event);
    void notify(const AnalysisEvent &analysis_event);

    bool registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr );

    // Execution functions:
    void resetAll();

    // Execute the control loop for the set analysis type.
    bool run();

    // This function performs the final initializations.  Mainly, it initializes
    // the two data container classes, DataStore and StepErrorControl.  It also
    // registers the neccessary vectors with the LAS system class.
    bool initializeAll(N_LOA_Loader * tmpLoaderPtr = 0);

    // Returns the current partial time derivative
    double partialTimeDerivative();

    // Gets the next time-step value.
    double getTime() const;

    // Gets the current time-step value.
    double getCurrentTime() const;

    // Gets the final time-step value.
    double getFinalTime() const;

    // Gets the initial time-step value.
    double getInitialTime() const;

    // Gets the starting time-step value.
    double getStartingTimeStep();

    // Updates the divided difference values.
    bool updateDivDiffs();

    // Calls the time int. method to update the corrector derivatives.
    bool updateDerivs();

    // updates state vectors following initial solve with previous operating point constraints
    bool completeOPStartStep();

    // updates the State vectors. This function is called from the LOCA interface.
    bool completeHomotopyStep
      ( const std::vector<std::string> & paramNames,
        const std::vector<double> & paramVals,
        N_LAS_Vector * solnVecPtr );

    bool failHomotopyStep ();

    // This function equates the 6 temporary vectors with their "next" vector
    // equivalents.  This function is neccessary for the nonlinear solver damping
    // loop.
    bool equateTmpVectors();

    // Compute an estimate of the error in the integration step.
    bool updateDerivsBlock(const std::list< index_pair > & solGIDList,
                           const std::list< index_pair > & staGIDList);

    // Prints out time loop information.
    bool printLoopInfo(int start, int finish);

    // Gets the time-integration method.
    int getTimeIntMode();

    // Get the steady-state flag (true if time int mode is none)
    bool getDCOPFlag ();

    // Get the dcop flag for transient.(true if doing DCOP for transient initialization)
    bool getTranOPFlag ();

    // Get the dcop flag for AC.(true if doing DCOP for AC initialization)
    bool getACOPFlag ();

    bool getDCSweepFlag ();

    bool getDotOpFlag() {return dotOpSpecified_;};
    bool getStepFlag() {return stepLoopFlag_;};

    bool getSweepSourceResetFlag () {return sweepSourceResetFlag_;};
    void setSweepSourceResetFlag (bool ssrf) { sweepSourceResetFlag_=ssrf;} ;

    bool getTransientFlag () const;

    // Is the doubleDCOP algorithm enabled?
    bool getDoubleDCOPEnabled ();

    void setCurrentMode(CurrentMode current_mode) {
      currentMode_ = current_mode;
    }

    CurrentMode getCurrentMode() const {
      return currentMode_;
    }

    // Get block analysis information for HB
    void setBlockAnalysisFlag( bool flagVal ) { blockAnalysisFlag_ = flagVal; }
    bool getBlockAnalysisFlag() const;

    void setHBFlag( bool flagVal ) { hbFlag_ = flagVal; }
    bool getHBFlag () {return hbFlag_; }

    // gets the index of the DCOP step.
    // 0 = nonlinear poisson, 1=full TCAD
    int getDoubleDCOPStep();

    // Gets/sets the step number.
    int getStepNumber();
    int getTranStepNumber();

    void setStepNumber(int step);
    void setTranStepNumber(int step);

    // This is true only at the beginning of integration, not at a breakpoint.
    bool getInitTranFlag();

    // Gets the step number.
    //int getStepLoopIter () { return stepLoopIter_; }

    // Returns the time integration order
    int getOrder ();

    int getNumberOfSteps ();
    int getUsedOrder ();
    int getNscsco ();

    const IO::CmdParse &getCommandLine() const {
      return commandLine_;
    }

    // Returns the "current" time step size.
    double getCurrentStepSize();

    // Returns the "last" time step size.
    double getLastStepSize();

    // Returns the breakpoint tolerance.
    double getBreakpointTol();

    // Sets the breakpoint tolerance.
    void setBreakpointTol(double bptol);

    // returns whether transient analysis is completed
    bool isSimulationComplete();

    // Gets the size of the restart data (bytes?).
    int restartDataSize( bool pack );

    // Sets the transient calculations parameters
    bool setTranAnalysisParams(const Util::OptionBlock & paramsBlock);

    // Sets the DC sweep calculation parameters.
    bool setDCAnalysisParams(const Util::OptionBlock & paramsBlock);

    // Method to handle OP statements.
    bool setOPAnalysisParams(const Util::OptionBlock & paramsBlock);

    // Sets the STEP calculation parameters.
    bool setSTEPAnalysisParams(const Util::OptionBlock & paramsBlock);

    // Sets the SAVE parameters.
    bool setSaveOptions(const Util::OptionBlock & OB);

    // Sets the DCOP restart parameters.
    bool setDCOPRestartParams(const Util::OptionBlock & OB);

    // Method to register the AC Analysis options.
    bool setACAnalysisParams(const Util::OptionBlock & OB);

    // Method to register the MOR Analysis options.
    bool setMORAnalysisParams(const Util::OptionBlock & OB);

    // Method to register the MOR utility options.
    bool setMOROptions(const Util::OptionBlock & OB);

    // sets a time at which to pause the simulation
    void setPauseTime(double pauseTime);
    // returns time at which to pause the simulation
    double getPauseTime();

    // returns true if the simulation is currently paused
    bool isPaused();

    // signal that simulation is being resumed from previously paused state
    void resumeSimulation();
    // reset the resume flag to false.
    void unset_resumeSimulation();

    // Registers the options block reference.
    bool setTranOptions(const Util::OptionBlock & OB);

    // Method to register the MPDE Analysis options.
    bool setMPDEAnalysisParams(const Util::OptionBlock & OB);

    // Method to register the MPDE utility options.
    bool setMPDEOptions(const Util::OptionBlock & OB);

    // Method to register the HB Analysis options.
    bool setHBAnalysisParams(const Util::OptionBlock & OB);

    // Method to register the HB utility options.
    bool setHBOptions(const Util::OptionBlock & OB);

    // Method to register the linear solver / preconditioning options.
    bool setLinSol(const Util::OptionBlock & OB);

    // Method to register the HB linear solver / preconditioning options.
    bool setHBLinSol(const Util::OptionBlock & OB);

    // Method to register the MPDE utility options.
    bool setTRANMPDEOptions(const Util::OptionBlock & OB);

    // Method to register the sensitivity options.
    bool setSensOptions(const Util::OptionBlock & OB);

    // Registers the TIA parameters block reference.
    bool registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp);

    // Registers the linear system pointer.
    bool registerLinearSystem(N_LAS_System * linearSystem_tmp);

    // Registers the nonlinear system manager pointer.
    bool registerNLSManager(Nonlinear::Manager * nlsMgrPtr_tmp);

    // Registers the nonlinear loader pointer.
    bool registerLoader(N_LOA_Loader * loader_tmp);

    // Registers the output manager pointer.
    bool registerOutputMgr(IO::OutputMgr * outputPtr_tmp);

    // Registers the restart manager pointer.
    bool registerRestartMgr(IO::RestartMgr * restartPtr_tmp);

    // Method to register the Device inerfaace pointer
    bool registerDeviceInterface ( N_DEV_DeviceInterface * devInterfacePtr );

    // Method to register the Topology pointer
    bool registerTopology( N_TOP_Topology * topoMgrPtr );

    // Method to register the Restart Manager pointer
    bool registerRestartManager( IO::RestartMgr * resMgrPtr );

    // Method to register the Output Manager pointer
    bool registerOutputManager( IO::OutputMgr * outMgrPtr );

    // Method to register the Application Builder pointer
    bool registerApplicationBuilder( N_LAS_Builder * appBuilderPtr );

    // Registers the parallel services manager pointer.
    bool registerParallelServices(N_PDS_Manager * pds_tmp);

    // Registers the restart intervals.
    bool registerRestartIntervals();

    // Registers the restart output intervals.
    bool registerOutputIntervals();

    // Registers the elapsed time timer
    bool registerElapsedTimer(Util::Timer *);

    // Writes-out the restart data.
    bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack);

    // Restores the restart data.
    bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

    IO::RestartMgr &getRestartManager() {
      return *restartPtr_;
    }

    // Gets the solution variable data.
    bool getSolnVarData(const int & gid, std::vector< double > & varData);

    // Gets the state variable data.
    bool getStateVarData(const int & gid, std::vector< double > & varData);

    // Gets the store variable data.
    bool getStoreVarData(const int & gid, std::vector< double > & varData);

    // Sets the solution variable data.
    bool setSolnVarData(const int & gid, const std::vector< double > & varData);

    // Sets the state variable data.
    bool setStateVarData(const int & gid, const std::vector< double > & varData);

    // Sets the store variable data.
    bool setStoreVarData(const int & gid, const std::vector< double > & varData);

    // set/get beginning integration flag
    void setBeginningIntegrationFlag(bool bif);
    bool getBeginningIntegrationFlag();

    // set/get integration method
    void setIntegrationMethod(int im);
    unsigned int getIntegrationMethod();

    // Gets the total time spent in the linear solvers.
    double getTotalLinearSolutionTime() const;

    // Gets the total time spent in the residual load calculations.
    double getTotalResidualLoadTime() const;

    // Gets the total time spent in the Jacobian load calculations.
    double getTotalJacobianLoadTime() const;

    // set the next solution vector pointer.  Needed for NOX...
    bool setNextSolVectorPtr (N_LAS_Vector * solVecPtr);

    // Habanero API mixed signal functions:
    bool provisionalStep (double maxTimeStep, double &currTimeStep);
    void acceptProvisionalStep ();
    void rejectProvisionalStep ();

    // Two-level Newton API functions:
    // Execute the control loop for the set analysis type,
    // for a set number of steps.
    void setExternalSolverState (const N_DEV_SolverState & ss);
    bool runStep
      (const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError);

    void conductanceTest ();
    bool startupSolvers ();
    bool finishSolvers ();

    void homotopyStepSuccess
      ( const std::vector<std::string> & paramNames,
        const std::vector<double> & paramVals);

    void homotopyStepFailure ();

    void stepSuccess(CurrentMode analysisUpper);
    void stepFailure(CurrentMode analysisUpper);
    bool getInitialQnorm (N_TIA_TwoLevelError & tle);
    bool getBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes);
    bool startTimeStep (const N_TIA_TimeIntInfo & tiInfo);


    // routines to get/set Dakota run flags and actually run a Dakota iteration
    bool getDakotaRunFlag();
    void setDakotaRunFlag( bool flag );
    int getDakotaIteration();
    void setDakotaIteration( int iterNumber );
    void getTimeIntInfo (N_TIA_TimeIntInfo & tiInfo);

    void initializeTransientModel();
    bool evalTransientModel(
        double t,
        N_LAS_Vector * SolVectorPtr,
        N_LAS_Vector * CurrSolVectorPtr,
        N_LAS_Vector * LasSolVectorPtr,
        N_LAS_Vector * StaVectorPtr,
        N_LAS_Vector * CurrStaVectorPtr,
        N_LAS_Vector * LasStaVectorPtr,
        N_LAS_Vector * StaDerivVectorPtr,
        N_LAS_Vector * StoVectorPtr,
        N_LAS_Vector * CurrStoVectorPtr,
        N_LAS_Vector * LasStoVectorPtr,
        N_LAS_Vector * stoLeadCurrQVectorPtr,
        N_LAS_Vector * QVectorPtr,
        N_LAS_Vector * FVectorPtr,
        N_LAS_Vector * BVectorPtr,
        N_LAS_Vector * dFdxdVpVectorPtr,
        N_LAS_Vector * dQdxdVpVectorPtr,
        N_LAS_Matrix * dQdxMatrixPtr,
        N_LAS_Matrix * dFdxMatrixPtr
        );
    bool evalTransientModelState(
        double t,
        N_LAS_Vector * SolVectorPtr,
        N_LAS_Vector * StaVectorPtr,
        N_LAS_Vector * StoVectorPtr
        );

    //-----------------------------------------------------------------------------
    // Function      : AnalysisManager::getTIAParams
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : Eric R. Keiter, SNL, Computational Sciences
    // Creation Date : 11/05/04
    //-----------------------------------------------------------------------------
    N_TIA_TIAParams &getTIAParams()
    {
      return tiaParams_;
    }

  private :

    bool getInputOPFlag ();

    // Allocate and register the MPDE Manager
    void setupMPDEMgr_();
    bool doMPDEregistrations_();

    // allocate analysis objects:
    void allocateAnalysisObject_();

    void initializeIntegrationProcess_();

    void computeDividedDifferences_();

    // Sets the nonlinear solver solution parameters.
    void setNLSParams_();

  public:
    // Queries about output and restart times
    bool outputIntervalSpecified();

    bool testRestartSaveTime();
    // Queries to the MPDE Manager.
    // keep there like this so we don't have to
    // expose the MPDE Manager pointer
    void setMPDEFlag( bool flagVal ) { mpdeFlag_ = flagVal; }
    bool getMPDEFlag ();    // "get" function for MPDE flag. (true if not IC)

    bool getMPDEIcFlag();   // get function for MPDE initial condition flag (true if MPDE & IC)
    bool getMPDEStartupFlag();   // True if you have done an initial transient simulation before starting MPDE IC calculation
    bool getWaMPDEFlag ();  // "get" function for WaMPDE flag. (true if not IC)

    // For DCOP restart and/or .SAVE files.
    bool testDCOPOutputTime_();
    bool testSaveOutputTime_();

  public:
    N_MPDE_Manager *getMPDEManager() const { return mpdeMgrPtr_; }

    N_TIA_DataStore *getTIADataStore() {
      return tiaDataStore_;
    }

    const AnalysisBase &getAnalysisObject() const {
      return *primaryAnalysisObject_;
    }

    void silenceProgress();

    void enableProgress();

  Device::DeviceInterface *getDeviceInterface() {
    return devInterfacePtr_;
  }

  N_LOA_NonlinearEquationLoader *getNonlinearEquationLoader() {
    return nonlinearEquationLoaderPtr_;
  }

  N_LAS_Builder *getAppBuilder() {
    return appBuilderPtr_;
  }

  Nonlinear::Manager *getNonlinearSolverManager() const {
    return nlsMgrPtr_;
  }

  N_TIA_StepErrorControl *&getStepErrorControl() {
    return stepErrorControl_;
  }

  N_TIA_WorkingIntegrationMethod *&getWorkingIntgMethod() {
    return workingIntgMethod_;
  }

  N_PDS_Manager *getPDSManager() const {
    return pdsMgrPtr_;
  }

  IO::OutputMgr *getOutputManager() const {
    return outMgrPtr_;
  }

  N_LAS_System *getLinearSystem() const {
    return linearSystem_;
  }

  bool getSwitchIntegrator() const {
    return switchIntegrator_;
  }

  void setSwitchIntegrator(bool switch_itegrator) {
    switchIntegrator_ = switch_itegrator;
  }

  void setNextOutputTime(double next_output_time) {
    nextOutputTime_ = next_output_time;
  }

  double getNextOutputTime() const {
    return nextOutputTime_;
  }

  double getInitialOutputInterval() const {
    return initialOutputInterval_;
  }

  void setStepLoopInitialized(bool step_loop_initialized) {
    stepLoopInitialized_ = step_loop_initialized;
  }

  const std::vector< std::pair < double, double > > &getOutputIntervals() const {
    return outputIntervals_;
  }

  bool isStepLoopInitialized() const {
    return stepLoopInitialized_;
  }
  
  Util::Timer &getXyceTranTimer() {
    return *xyceTranTimerPtr_;
  }

  OutputMgrAdapter &getOutputManagerAdapter() const {
    return *outputManagerAdapter_;
  }

  N_TIA_StepErrorControl *getSECPtr() const {
    return stepErrorControl_;
  }

  int getDebugLevel() const {
    return tiaParams_.debugLevel;
  }

  Topo::Topology &getTopology() {
    return *topoMgrPtr_;
  }

  N_TIA_WorkingIntegrationMethod &getWorkingIntegrationMethod() {
    return *workingIntgMethod_;
  }

  N_LOA_Loader &getLoader() {
    return *loader_;
  }

  void setAnalysisMode(Analysis_Mode analysis_mode) {
    analysisMode_ = analysis_mode;
  }

  Analysis_Mode getAnalysisMode() const {
    return analysisMode_;
  }

  double getSolverStartTime() const {
    return solverStartTime_;
  }

  double getStartTranTime() const {
    return startTRANtime_;
  }
  
  void setStartTranTime(double start_tran_time) {
    startTRANtime_ = start_tran_time;
  }

  bool getProgressFlag() const {
    return progressFlag_;
  }

  double getSaveTime() const {
    return saveTime_;
  }

  N_PDS_Comm &getPDSComm();

  private:
    IO::ActiveOutput *    activeOutput_;

    IO::CmdParse &                      commandLine_;                   ///< Command line object

    N_TIA_TIAParams                     tiaParams_;                     ///< Current time-integration method parameters
    N_TIA_WorkingIntegrationMethod *    workingIntgMethod_;             ///< Working intergration method
    N_TIA_StepErrorControl *            stepErrorControl_;              ///< Pointer to the TIA step-error control object.
    N_LAS_System *                      linearSystem_;                  ///< Pointer to the linear system information and containers
    Nonlinear::Manager *                nlsMgrPtr_;                     ///< Pointer to the nonlinear solver manager.
    N_LOA_Loader *                      loader_;                        ///< Pointer to the nonlinear loader object.
    N_LOA_CktLoader *                   cktLoaderPtr_;                  ///< 'real' pointer to the ckt-loader.
    IO::RestartMgr *                    restartPtr_;                    ///< Pointer to the restart manager.
    IO::PkgOptionsMgr *                 pkgOptMgrPtr_;                  ///< package options manager
    N_LOA_NonlinearEquationLoader *     nonlinearEquationLoaderPtr_;    ///< Pointer to the application loader
    Device::DeviceInterface *           devInterfacePtr_;               ///< Pointer to the device interface
    Topo::Topology *                    topoMgrPtr_;                    ///< Pointer to the topology manager
    IO::OutputMgr *                     outMgrPtr_;                     ///< Pointer to the output manager
    N_LAS_Builder *                     appBuilderPtr_;                 ///< Pointer to the applicaiton builder
    N_PDS_Manager *                     pdsMgrPtr_;                     ///< Pointer to the parallel services manager.
    OutputMgrAdapter *                  outputManagerAdapter_;          ///< Output manager adapter
    N_TIA_DataStore *                   tiaDataStore_;                  ///< TIA data store object.

    Analysis::Analysis_Mode             analysisMode_;

    bool analysisParamsRegistered;

    bool firstTime;
    double oldPercentComplete;
    double startSimTime;

    bool calledBeforeTwoLevelTran_;

    // Switch the integration flag.
    bool switchIntegrator_;

    // DC Operating Point flag.
    bool initializeAllFlag_;      // true if the initializeAll function has been
                                  // called once.

    double startTRANtime_;

    bool stepLoopFlag_;           // true if there is an external
                                  // parameter stepper loop around the dcop or
                                  // transient simulation loops.

    bool stepLoopInitialized_;    // true if the step loop has been set up.
    bool dcLoopInitialized_;      // true if the dc sweep loop has been set up.
    bool gui_;                    // command line arg -gui is present


    bool daeStateDerivFlag_;   // true if running new-DAE and need
                              // state derivative.  true by default.
                              // If set to true, it breaks MPDE.

    bool initializeSolvers_mixedSignal_;

    int dcLoopSize_;

    bool sweepSourceResetFlag_;

    Stats::Stat         rootStat_;

    // Flag to decide whether to print progress
    bool progressFlag_;

    // Xyce timing utility for timing the transient simulation CPU time.
    Util::Timer *       xyceTranTimerPtr_;

    // Xyce timing utility for timing elapsed run time
    Util::Timer *       elapsedTimerPtr_;

    double solverStartTime_;

    bool dakotaRunFlag_;
    int dakotaIterationNumber_;

    // for .SAVE and/or DCOP restart.
    double saveTime_;
    bool saveTimeGiven_;
    bool saveFlag_;
    bool savedAlready_;
    bool dcopRestartFlag_;

    // .OP flag(s)
    bool dotOpSpecified_;

    // output and restart interval info
    double initialOutputInterval_;
    std::vector< std::pair < double, double > > outputIntervals_;
    double nextOutputTime_;

    double initialRestartInterval_;
    std::vector< std::pair < double, double > > restartIntervals_;
    double nextRestartSaveTime_;

    // for HB, MPDE, or any other block analysis type
    bool blockAnalysisFlag_;
    bool hbFlag_;
    bool mpdeFlag_;

    // sensitivity flag(s)
    bool sensFlag_;

    // ref counted pointers for various analyses. Not all are used in every simulation
    Teuchos::RefCountPtr<AnalysisBase> analysisObject_;
    Teuchos::RefCountPtr<AnalysisBase> stepAnalysisTarget_;
    Teuchos::RefCountPtr<AnalysisBase> dakotaAnalysisTarget_;
    Teuchos::RefCountPtr<AnalysisBase> primaryAnalysisObject_;

    // MPDE stuff:
    N_MPDE_Manager *                    mpdeMgrPtr_;
    RefCountPtr<N_TIA_MPDEInterface> tiaMPDEIfacePtr_;

    // Two level Newton API:
    RefCountPtr<AnalysisBase>           twoLevelAnalysisObject_;

    // Habanero mixed-signal API:
    RefCountPtr<AnalysisBase>           mixedSignalAnalysisObject_;

    // parameter objects for each analysis type.  These are saved until the
    // analysis object of choice is allocated.
    Util::OptionBlock tranParamsBlock;
    Util::OptionBlock opParamsBlock;
    Util::OptionBlock acParamsBlock;
    Util::OptionBlock morParamsBlock;

    // Different block passed in for each sweep variable, so need to have
    // these containers be vectors.
    std::vector<Util::OptionBlock> dcParamsBlockVec;
    std::vector<Util::OptionBlock> stepParamsBlockVec;

    Util::OptionBlock mpdeParamsBlock;
    Util::OptionBlock hbParamsBlock;
    Util::OptionBlock hbOptionsBlock;
    Util::OptionBlock hbLinSolBlock;
    Util::OptionBlock linSolBlock;
    Util::OptionBlock dakotaParamsBlock;

    CurrentMode currentMode_;

  public:
    unsigned int breakPointRestartStep;
};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::AnalysisManager N_ANP_AnalysisManager;

#endif // Xyce_N_ANP_AnalysisManager_h
