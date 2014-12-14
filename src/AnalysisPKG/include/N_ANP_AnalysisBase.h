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
// Filename       : $RCSfile: N_ANP_AnalysisBase.h,v $
//
// Purpose        : Base class for Analysis types
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
// Revision Number: $Revision: 1.31.2.1 $
//
// Revision Date  : $Date: 2014/08/28 21:00:43 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AnalysisBase_h
#define Xyce_N_ANP_AnalysisBase_h

#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_IO_fwd.h>
#include <N_TIA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_UTL_JSON.h>

class N_TIA_Assembler;
class N_LAS_System;
class N_LOA_Loader;

namespace Xyce {
namespace Analysis {

struct StatCounts 
{
  StatCounts();
  StatCounts &operator+=(const StatCounts &stats);

  unsigned int        successfulStepsTaken_;        ///< Number of consecutive successful time-integration steps.
  unsigned int        successStepsThisParameter_;
  unsigned int        failedStepsAttempted_;        ///< Total number of failed time-integration steps.
  unsigned int        jacobiansEvaluated_;
  unsigned int        iterationMatrixFactorizations_;
  unsigned int        linearSolves_;
  unsigned int        failedLinearSolves_;
  unsigned int        linearIters_;
  unsigned int        residualEvaluations_;
  unsigned int        nonlinearConvergenceFailures_;
  double              linearSolutionTime_;
  double              residualLoadTime_;
  double              jacobianLoadTime_;
};

Util::JSON &operator<<(Util::JSON &json, const StatCounts &s);

//-------------------------------------------------------------------------
// Class         : AnalysisBase
// Purpose       : Base class for common analysis functions
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class AnalysisBase
{
  public: 
    AnalysisBase(AnalysisManager &analysis_manager);
    virtual ~AnalysisBase();

  public:
    virtual bool setAnalysisParams(const Util::OptionBlock & paramsBlock) {return true;}
    virtual bool outputFailureStats() {return true;}

    virtual void setParamsWithOutputMgrAdapter(OutputMgrAdapter & outputManagerAdapter) {}

    virtual int getStepIter () { return 0; }
    virtual int getStepNumber () { return stepNumber; }
    virtual void setStepNumber (int step) { stepNumber=step; }

    virtual void setTranStepNumber (int step) { tranStepNumber=step; }
    virtual int getTranStepNumber () { return tranStepNumber; }

    virtual void setSensFlag() {sensFlag_=true;}
    virtual bool getDCOPFlag() = 0;

    virtual bool run() = 0;
    virtual bool init() = 0;
    virtual bool loopProcess() = 0;
    virtual bool processSuccessfulStep() = 0;
    virtual bool processFailedStep() = 0;
    virtual bool finish() = 0;
    virtual bool handlePredictor() = 0;


    //-----------------------------------------------------------------------------
    // Function      : AnalysisBase::printStepHeader()
    // Purpose       : Prints out step information.
    // Special Notes : 
    // Scope         : public
    // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
    // Creation Date : 6/26/00
    //-----------------------------------------------------------------------------
    virtual void printStepHeader(std::ostream &os) 
    {}

    //-----------------------------------------------------------------------------
    // Function      : AnalysisBase::printProgress()
    // Purpose       : Outputs run completion percentage and estimated
    //                 time-to-completion.
    // Special Notes : 
    // Scope         : public
    // Creator       : Scott A. Hutchinson, SNL, Computational Sciences
    // Creation Date : 06/07/2002
    //-----------------------------------------------------------------------------
    virtual void printProgress(std::ostream &os) 
    {}

    // mixed-signal
    virtual void preStepDetails(double maxTimeStepFromHabanero) {}
    virtual bool mixedSignalStep() { return true; }
    virtual bool finalizeStep() { return true; }

    // Two Level specific
    virtual bool twoLevelStep() { return true; }

    // Utility function for AnalysisManager to determine if a complex analysis type (like HB)
    // is performing another type of analysis under the hood.  This is necessary for populating
    // the N_TIA_TimeIntInfo struct.  For straightforward analysis types, this method is not
    // needed because the AnalysisManager already knows which analysis is being performed.
    virtual bool isAnalysis( int analysis_type ) { return false; }

    virtual bool printLoopInfo(int start, int finish);

    virtual void setBeginningIntegrationFlag(bool bif) {beginningIntegration = bif;}
    virtual bool getBeginningIntegrationFlag() {return beginningIntegration;}

    virtual void setIntegrationMethod (int im) {integrationMethod_= im;}
    virtual unsigned int getIntegrationMethod () {return integrationMethod_;}

    virtual bool getInputOPFlag(){return inputOPFlag_;}

    bool resetForStepAnalysis();
    void resetAll();
    int saveLoopInfo ();

    // step statistic functions
    void gatherStepStatistics_();
    double getTotalLinearSolutionTime() const;
    double getTotalResidualLoadTime() const;
    double getTotalJacobianLoadTime() const;

    // functions related to PDE problems requiring two DCOP solves
    bool getDoubleDCOPEnabled();
    virtual int getDoubleDCOPStep();
    bool firstDoubleDCOPStep_();

  const StatCounts &getStatCounts(int index = -1) const
  {
    if (index == -1)
      return saveStatCountsVector_.back();
    else
      return saveStatCountsVector_[index];
  }

protected:
    AnalysisManager &                           analysisManager_;
    N_LAS_System &                              linearSystem_;
    N_LOA_Loader &                              loader_;
    N_LOA_Loader &                              nonlinearEquationLoader_;
    Nonlinear::Manager &                        nonlinearSolverManager_;
    OutputMgrAdapter &                          outputManagerAdapter_;
    N_TIA_StepErrorControl *&                   stepErrorControl_;              ///< Ref to pointer since AnalysisManager will change pointer
    N_TIA_WorkingIntegrationMethod *&           workingIntgMethod_;             ///< Ref to pointer since AnalysisManager will change pointer
    N_TIA_TIAParams &                           tiaParams_;

    bool                                        beginningIntegration;

    unsigned int                                integrationMethod_;             ///< Current time-integration method flag.

    // NOTE: For now tranStepNumber is the same as stepNumber, but later I will
    //       change it.  stepNumber will later include both dcop and tran steps.
    //       I haven't changed it yet b/c I need to check what devices call
    //       getStepNumber, and what they expect to get.

    unsigned int                                stepNumber;                     ///< Time-integration step number counter.
    unsigned int                                tranStepNumber;

    bool                                        doubleDCOPFlag_;                ///< true if doing a double-DCOP is possible.
    int                                         doubleDCOPStep_;                ///< current step in the DCOP loop.

    bool                                        sensFlag_;
    bool                                        inputOPFlag_;                   ///< true if starting from an initial condition.

    std::vector<StatCounts>                     saveStatCountsVector_;

public:
    StatCounts                                  stats_;
};

StatCounts operator-(const StatCounts &s0, const StatCounts &s1);

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalLinearSolutionTime() const
{
  return stats_.linearSolutionTime_;
}

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalResidualLoadTime() const
{
  return stats_.residualLoadTime_;
}

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalJacobianLoadTime() const
{
  return stats_.jacobianLoadTime_;
}

//-----------------------------------------------------------------------------
inline bool AnalysisBase::getDoubleDCOPEnabled () 
{
  return doubleDCOPFlag_;
}

//-----------------------------------------------------------------------------
inline int AnalysisBase::getDoubleDCOPStep () 
{
  return doubleDCOPStep_ ;
}

bool updateSweepParams(N_LOA_Loader &loader, AnalysisManager &analysis_manager, int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);
int setupSweepLoop(N_LOA_Loader &loader, int debug_level, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);

} // namespace Analysis
} // namespace Xyce


typedef Xyce::Analysis::AnalysisBase N_ANP_AnalysisBase;

#endif

