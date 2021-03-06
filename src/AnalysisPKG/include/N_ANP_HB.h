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
// Filename       : $RCSfile: N_ANP_HB.h,v $
//
// Purpose        : HB analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Todd Coffey, 1414, Ting Mei 1437
//
// Creation Date  : 07/23/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.41.2.1 $
// Revision Date  : $Date: 2014/08/28 21:00:43 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_HB_h
#define Xyce_N_ANP_HB_h

#include <vector>

#include <N_ANP_fwd.h>
#include <N_LOA_fwd.h>
#include <N_LAS_fwd.h>
#include <N_ANP_AnalysisBase.h>
#include <N_ANP_StepEvent.h>
#include <N_MPDE_State.h>

#include <N_UTL_DFTInterfaceDecl.hpp>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Listener.h>

class N_MPDE_Discretization;

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : HB
// Purpose       : HB analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class HB : public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
public:
    HB( AnalysisManager &anaManagerPtr );
    virtual ~HB();

    void notify(const StepEvent &event);

  // Method to set HB options
    bool setHBOptions(const N_UTL_OptionBlock & OB);

    // Method to set HB linear solver / preconditioning options
    bool setHBLinSol(const N_UTL_OptionBlock & OB);

    // Method to set non-HB linear solver / preconditioning options (needed for .STEP)
    bool setLinSol(const N_UTL_OptionBlock & OB);

    // Override these methods using the current analysisObject_.
    int getStepNumber ();
    void setStepNumber (int step);
    void setBeginningIntegrationFlag(bool bif);
    bool getBeginningIntegrationFlag();
    void setIntegrationMethod (int im);
    unsigned int getIntegrationMethod ();

    int getDoubleDCOPStep();

    bool getDCOPFlag();

    bool run(); 
    bool init(); 
    bool loopProcess(); 
    bool processSuccessfulDCOP(); 
    bool processFailedDCOP();
    bool processSuccessfulStep();
    bool processFailedStep();
    bool finish();
    bool handlePredictor();

    bool finalVerboseOutput();

  // Utility function for AnalysisManager to determine if a complex analysis type (like HB)
  // is performing another type of analysis under the hood.  This is necessary for populating
  // the N_TIA_TimeIntInfo struct.  For straightforward analysis types, this method is not
  // needed because the AnalysisManager already knows which analysis is being performed.
  bool isAnalysis( int analysis_type );

  // Transform the current solution vector for time domain and frequency domain output
  void prepareHBOutput(N_LAS_Vector & solnVecPtr,
                       std::vector<double> & timePoints,
                       std::vector<double> & freqPoints,
                       Teuchos::RCP<N_LAS_BlockVector> & timeDomainSolnVec,
                       Teuchos::RCP<N_LAS_BlockVector> & freqDomainSolnVecReal,
                       Teuchos::RCP<N_LAS_BlockVector> & freqDomainSolnVecImaginary,
                       Teuchos::RCP<N_LAS_BlockVector> & timeDomainStoreVec,
                       Teuchos::RCP<N_LAS_BlockVector> & freqDomainStoreVecReal,
                       Teuchos::RCP<N_LAS_BlockVector> & freqDomainStoreVecImaginary) const;

  int debugLevel;

private:

  // Add in solver info and timing info from current analysisObject_
  void accumulateStatistics_();

  bool runTol_();
  bool runStartupPeriods_();
  bool runTransientIC_();
  bool interpolateIC_();

  bool setFreqPoints_();

  bool setInitialGuess_();

  bool setTimePoints_();

  bool createFT_(); 
 
  // Flag to indicate of the simulation is paused
  bool isPaused;

  // Timing/loop count info
  double startDCOPtime, endTRANtime; // startTRANtime

  Device::DeviceInterface *             devInterfacePtr_;
  N_LOA_NonlinearEquationLoader *       nonlinearEquationLoaderPtr_;
  N_LAS_Builder *                       appBuilderPtr_;
  N_PDS_Manager *                       pdsMgrPtr_;
  Teuchos::RCP<AnalysisBase>            analysisObject_;

  // Current analysis state flags.
  bool isTransient_, isDCSweep_;

  //Testing Flag
  bool test_;

  // Problem Size
  int size_;

  std::vector<int> numPosFreqs;
  std::vector<int> numFreqs_;

  // Periodicity Information
  double period_;

  // Number of fast time periods to integrate over and IGNORE before
  // getting initial conditions for HB.  Default is zero.
  int startUpPeriods_;
  bool startUpPeriodsGiven_;

  bool startUpPeriodsFinished_;
  bool saveIcData_;

  // Stored copy of transient TIAParams
  N_TIA_TIAParams transTiaParams_;

  // Transient assisted HB.
  int taHB_;

  bool voltLimFlag_;
  int intmodMax_;      
  std::string method_;

  bool intmodMaxGiven_;

  // HB loader, builder, system, and DFT
  N_LOA_HBLoader *              hbLoaderPtr_;
  Teuchos::RCP<N_LAS_HBBuilder> hbBuilderPtr_;
  Teuchos::RCP<N_LAS_System> lasHBSysPtr_;
  Teuchos::RCP<N_UTL_FFTInterface<std::vector<double> > > ftInterface_;
  std::vector<double> ftInData_, ftOutData_, iftInData_, iftOutData_;

  // Time discretization
  int fastTimeDisc_;
  int fastTimeDiscOrder_;
  std::vector<double> fastTimes_;
  std::vector<double> timeSteps_;
  std::vector<double> freqPoints_;
  Teuchos::RCP<N_MPDE_Discretization> mpdeDiscPtr_;
  N_MPDE_State mpdeState_;

  // Fourier matrices

  Teuchos::RCP<N_UTL_DFTInterfaceDecl<std::vector<double> > > dftInterface_;
  Teuchos::SerialDenseMatrix<int,double> idftMatrix_, dftMatrix_;

  // Linear solver and nonlinear solver options
  N_UTL_OptionBlock saved_lsHBOB_;
  N_UTL_OptionBlock saved_lsOB_;
  N_UTL_OptionBlock saved_nlHBOB_;

  // An analysis-dependent preconditioner factory.
  Teuchos::RCP<N_LAS_PrecondFactory> precFactory_;

  // Local storage vectors
  Teuchos::RCP<N_LAS_Vector> dcOpSolVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpStateVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpQVecPtr_;
  Teuchos::RCP<N_LAS_Vector> dcOpStoreVecPtr_;

  std::vector<double> goodTimePoints_;
  std::vector<Teuchos::RCP<N_LAS_Vector> > goodSolutionVec_;
  std::vector<Teuchos::RCP<N_LAS_Vector> > goodStateVec_;
  std::vector<Teuchos::RCP<N_LAS_Vector> > goodQVec_;
  std::vector<Teuchos::RCP<N_LAS_Vector> > goodStoreVec_;

  // HB initial condition
  Teuchos::RCP<N_LAS_BlockVector> HBICVectorPtr_;
  Teuchos::RCP<N_LAS_BlockVector> HBICVectorFreqPtr_;

  // HB initial state condition
  Teuchos::RCP<N_LAS_BlockVector> HBICStateVectorPtr_;
//  Teuchos::RCP<N_LAS_BlockVector> HBICStateVectorFreqPtr_;

  // HB initial Q condition
  Teuchos::RCP<N_LAS_BlockVector> HBICQVectorPtr_;
 //  Teuchos::RCP<N_LAS_BlockVector> HBICQVectorFreqPtr_;

  // HB initial store condition
  Teuchos::RCP<N_LAS_BlockVector> HBICStoreVectorPtr_;

  // HB statistics
  StatCounts         hbStatCounts_;

  // int hbTotalNumberSuccessfulStepsTaken_;
  // int hbTotalNumberFailedStepsAttempted_;
  // int hbTotalNumberJacobiansEvaluated_;
  // int hbTotalNumberIterationMatrixFactorizations_;
  // int hbTotalNumberLinearSolves_;
  // int hbTotalNumberFailedLinearSolves_;
  // int hbTotalNumberLinearIters_;
  // int hbTotalNumberResidualEvaluations_;
  // int hbTotalNonlinearConvergenceFailures_;
  // double hbTotalResidualLoadTime_;
  // double hbTotalJacobianLoadTime_;
  // double hbTotalLinearSolutionTime_;

  bool resetForStepCalledBefore_;
};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::Transient N_ANP_Transient;

#endif // Xyce_N_ANP_HB_h

