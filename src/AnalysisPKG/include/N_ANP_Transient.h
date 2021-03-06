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
// Filename       : $RCSfile: N_ANP_Transient.h,v $
//
// Purpose        : Transient analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.37 $
// Revision Date  : $Date: 2014/08/05 21:41:35 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Transient_h
#define Xyce_N_ANP_Transient_h

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_StepEvent.h>
#include <N_UTL_FixedQueue.h>
#include <N_UTL_Listener.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : Transient
// Purpose       : Transient analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Transient : public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
  public:
    Transient(AnalysisManager &analysis_manager);

    virtual ~Transient() {}

  void notify(const StepEvent &event);
  
    bool setAnalysisParams(const Util::OptionBlock & paramsBlock);

    bool run();
    bool init();
    bool loopProcess();
    bool processSuccessfulDCOP();
    bool processFailedDCOP();
    bool processSuccessfulStep();
    bool processFailedStep();
    bool finish();
    bool handlePredictor();
    bool resetForHB();
    void outputQueuedData();
    bool finalVerboseOutput();

    void printStepHeader(std::ostream &os);
    void printProgress(std::ostream &os);


    // mixed-signal specific
    void preStepDetails (double maxTimeStepFromHabanero);
    bool mixedSignalStep();
    bool finalizeStep ();

    // Two Level specific
    bool twoLevelStep();

    void setDCOPFlag(bool flag) { dcopFlag_ = flag; }
    bool getDCOPFlag() { return dcopFlag_; }

    int getDCStats() { return dcStats; }
    int getTranStats() { return tranStats; }

  private:
    Parallel::Machine     comm_;

    // 05/26/09 Coffey,Schiek,Mei:  We moved these functions from AnalysisManager  ----v
    unsigned int initialIntegrationMethod_;
    bool firstTranOutput_;
    std::vector<double> outputInterpolationTimes_;
    void computeOutputInterpolationTimes_(double currTime);
    void updateOutputTime_(double currTime);
    bool testOutputTime_();

    void noopOutputs ();
    void tranopOutputs ();
    void tranStepOutputs ();
    // 05/26/09 Coffey,Schiek,Mei:  We moved these functions from AnalysisManager  ----^

    void takeAnIntegrationStep_();

    bool retakeAndAcceptTimeStep( double aTimeStep );

    // a flag to indicate of the simulation is paused
    bool isPaused;
    bool dcopFlag_;               // true if this is a DCOP calculation.

    double startDCOPtime, endTRANtime; // startTRANtime
    bool gui_;                    // command line arg -gui is present

    bool historyTrackingOn_;      // bool to indicate if history tracking is on.

    // These are used to track the minimum estimated error over tol. of failed
    // time steps so that if we're going to exit with a time step too small error,
    // we can have the option of accepting whatever step had the minimum error.
    double minEstErrorOverTol;
    int    stepNumberAtMinEstErrorOverTol;
    double timeStepAtMinEstErrorOverTol;
    // Timing/loop count info
    // int dcStats;
    //
    // for handling expressions given for a time dependent max time step
    bool maxTimeStepExpressionGiven_;
    std::string maxTimeStepExpressionAsString_;
    RefCountPtr<Util::ExpressionData> maxTimeStepExpressionRCPtr_;

    // here we store stats on the last few time steps
    // to report if Xyce fails. A user can use this info
    // to figure how how to make the simulation work.
    // The number of items saved is set in the constructor
    int queueSize_;
    Util::FixedQueue<double> timeQueue_;
    Util::FixedQueue<double> timeStepQueue_;
    Util::FixedQueue<int> stepStatusQueue_;
    Util::FixedQueue<double> estErrorOverTolQueue_;
    Util::FixedQueue<int> nonlinearSolverStatusQueue_;
    Util::FixedQueue<int> nonlinearSolverNumIterationsQueue_;
    Util::FixedQueue<double> nonlinearSolverMaxNormQueue_;
    Util::FixedQueue<double> nonlinearSolverMaxNormIndexQueue_;
    // Util::FixedQueue<double> nonlinearSolverNormQueue_;

    // solution variable names vector, used in outputting Queued data.
    std::vector<std::string> nameVec_;

    bool firstTime;
    double oldPercentComplete;
    double startSimTime;
    int dcStats;
    int tranStats;

    // Sensitivity data.  This will need a refactor to work in transient, but for
    // now this needs to be here to compile.
    std::vector<double> objectiveVec_; 
    std::vector<double> dOdpVec_; 
    std::vector<double> dOdpAdjVec_;
    std::vector<double> scaled_dOdpVec_;
    std::vector<double> scaled_dOdpAdjVec_;
};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::Transient N_ANP_Transient;

#endif // Xyce_N_ANP_Transient_h
