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
// Filename      : $RCSfile: N_ANP_AnalysisBase.C,v $
// Purpose       : Base class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.48.2.1 $
// Revision Date  : $Date: 2014/08/28 21:00:43 $
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisBase.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_SweepParam.h>
#include <N_ERH_Message.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_OptionBlock.h>
#include <N_LOA_NonlinearEquationLoader.h>

namespace Xyce {
namespace Analysis {

StatCounts::StatCounts()
  : successfulStepsTaken_(0),
    successStepsThisParameter_(0),
    failedStepsAttempted_(0),
    jacobiansEvaluated_(0),
    iterationMatrixFactorizations_(0),
    linearSolves_(0),
    failedLinearSolves_(0),
    linearIters_(0),
    residualEvaluations_(0),
    nonlinearConvergenceFailures_(0),
    linearSolutionTime_(0.0),
    residualLoadTime_(0.0),
    jacobianLoadTime_(0.0)
{}

StatCounts &
StatCounts::operator+=(
  const StatCounts &   stats)
{
  successfulStepsTaken_ += stats.successfulStepsTaken_;
  successStepsThisParameter_ += stats.successStepsThisParameter_;
  failedStepsAttempted_ += stats.failedStepsAttempted_;
  jacobiansEvaluated_ += stats.jacobiansEvaluated_;
  iterationMatrixFactorizations_ += stats.iterationMatrixFactorizations_;
  linearSolves_ += stats.linearSolves_;
  failedLinearSolves_ += stats.failedLinearSolves_;
  linearIters_ += stats.linearIters_;
  residualEvaluations_ += stats.residualEvaluations_;
  nonlinearConvergenceFailures_ += stats.nonlinearConvergenceFailures_;
  linearSolutionTime_ += stats.linearSolutionTime_;
  residualLoadTime_ += stats.residualLoadTime_;
  jacobianLoadTime_ += stats.jacobianLoadTime_;

  return *this;
}


Util::JSON &operator<<(Util::JSON &json, const StatCounts &s) 
{
  json << Util::JSON::open
       << Util::nameValuePair("successfulStepsTaken", s.successfulStepsTaken_) << Util::JSON::sep
       << Util::nameValuePair("failedStepsAttempted", s.failedStepsAttempted_) << Util::JSON::sep
       << Util::nameValuePair("jacobiansEvaluated", s.jacobiansEvaluated_) << Util::JSON::sep
       << Util::nameValuePair("iterationMatrixFactorizations", s.iterationMatrixFactorizations_) << Util::JSON::sep
       << Util::nameValuePair("linearSolves", s.linearSolves_) << Util::JSON::sep
       << Util::nameValuePair("failedLinearSolves", s.failedLinearSolves_) << Util::JSON::sep
       << Util::nameValuePair("linearIters", s.linearIters_) << Util::JSON::sep
       << Util::nameValuePair("residualEvaluations", s.residualEvaluations_) << Util::JSON::sep
       << Util::nameValuePair("nonlinearConvergenceFailures", s.nonlinearConvergenceFailures_) << Util::JSON::sep
       << Util::nameValuePair("residualLoadTime", s.residualLoadTime_) << Util::JSON::sep
       << Util::nameValuePair("jacobianLoadTime", s.jacobianLoadTime_) << Util::JSON::sep
       << Util::nameValuePair("linearSolutionTime", s.linearSolutionTime_)
       << Util::JSON::close;

  return json;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::AnalysisBase
// Purpose       : Constructor
// Special Notes : 
// Scope         : public
// Creator       : Todd S. Coffey, SNL.
// Creation Date : 01/29/08
//-----------------------------------------------------------------------------
AnalysisBase::AnalysisBase( AnalysisManager &analysis_manager )
  : analysisManager_(analysis_manager),
    linearSystem_(*analysis_manager.getLinearSystem()),
    loader_(analysis_manager.getLoader()),
    nonlinearEquationLoader_(*analysis_manager.getNonlinearEquationLoader()),
    nonlinearSolverManager_(*analysis_manager.getNonlinearSolverManager()),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    stepErrorControl_(analysis_manager.getStepErrorControl()),
    workingIntgMethod_(analysis_manager.getWorkingIntgMethod()),
    tiaParams_(analysis_manager.getTIAParams()),
    beginningIntegration(true),
    integrationMethod_(TIAMethod_NONE),
    stepNumber(0),
    tranStepNumber(0),
    doubleDCOPFlag_(false),
    doubleDCOPStep_(0),
    sensFlag_(false),
    inputOPFlag_(false),
    saveStatCountsVector_(),
    stats_()
{}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::~AnalysisBase()
// Purpose       : destructor
// Special Notes : 
// Scope         : 
// Creator       : Todd S. Coffey, SNL.
// Creation Date : 01/29/08
//-----------------------------------------------------------------------------
AnalysisBase::~AnalysisBase()
{}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::resetForStepAnalysis() 
// Purpose       : When doing a .STEP sweep, some data must be reset to its 
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool AnalysisBase::resetForStepAnalysis()
{
  stats_.successStepsThisParameter_ = 0;
  stepNumber = 0;
  beginningIntegration = true;

  return true;
}



//-----------------------------------------------------------------------------
// Function      : AnalysisBase::resetAll
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/16/10
//-----------------------------------------------------------------------------
void AnalysisBase::resetAll()
{
  stepNumber = 0;
  tranStepNumber = 0;

  stats_ = StatCounts();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::saveLoopInfo
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/16/10
//-----------------------------------------------------------------------------
int AnalysisBase::saveLoopInfo ()
{
  if (saveStatCountsVector_.empty()) // push back empty stats as sentinel
    saveStatCountsVector_.push_back(StatCounts());

  saveStatCountsVector_.push_back(stats_);

  return saveStatCountsVector_.size() - 1;
}

StatCounts
operator-(
  const StatCounts &   s0,
  const StatCounts &   s1)
{
  StatCounts s;

  s.successfulStepsTaken_ = s0.successfulStepsTaken_ - s1.successfulStepsTaken_;
  s.failedStepsAttempted_ = s0.failedStepsAttempted_ - s1.failedStepsAttempted_;
  s.jacobiansEvaluated_ = s0.jacobiansEvaluated_ - s1.jacobiansEvaluated_;
  s.iterationMatrixFactorizations_ = s0.iterationMatrixFactorizations_ - s1.iterationMatrixFactorizations_;
  s.linearSolves_ = s0.linearSolves_ - s1.linearSolves_;
  s.failedLinearSolves_ = s0.failedLinearSolves_ - s1.failedLinearSolves_;

  if (s0.linearIters_ > s1.linearIters_)
    s.linearIters_ = s0.linearIters_ - s1.linearIters_;

  s.residualEvaluations_ = s0.residualEvaluations_ - s1.residualEvaluations_;
  s.nonlinearConvergenceFailures_ = s0.nonlinearConvergenceFailures_ - s1.nonlinearConvergenceFailures_;
  s.residualLoadTime_ = s0.residualLoadTime_ - s1.residualLoadTime_;
  s.jacobianLoadTime_ = s0.jacobianLoadTime_ - s1.jacobianLoadTime_;
  s.linearSolutionTime_ = s0.linearSolutionTime_ - s1.linearSolutionTime_;

  return s;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::printLoopInfo
// Purpose       : Prints out time loop information.
// Special Notes : Prints stats from save point start to save point finish.
//                 Special case 0,0 is entire run to this point
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
bool AnalysisBase::printLoopInfo(int t1, int t2)
{
  bool bsuccess = true;

  if (t1 == 0 && t2 == 0)
  {
    t2 = saveLoopInfo();
  }

  lout() << "\tNumber Successful Steps Taken:\t\t" << saveStatCountsVector_[t2].successfulStepsTaken_ - saveStatCountsVector_[t1].successfulStepsTaken_ << std::endl
         << "\tNumber Failed Steps Attempted:\t\t" << saveStatCountsVector_[t2].failedStepsAttempted_ - saveStatCountsVector_[t1].failedStepsAttempted_ << std::endl
         << "\tNumber Jacobians Evaluated:\t\t" << saveStatCountsVector_[t2].jacobiansEvaluated_ - saveStatCountsVector_[t1].jacobiansEvaluated_ << std::endl
         << "\tNumber Iteration Matrix Factorizations:\t" << saveStatCountsVector_[t2].iterationMatrixFactorizations_ - saveStatCountsVector_[t1].iterationMatrixFactorizations_ << std::endl
         << "\tNumber Linear Solves:\t\t\t" << saveStatCountsVector_[t2].linearSolves_ - saveStatCountsVector_[t1].linearSolves_ << std::endl
         << "\tNumber Failed Linear Solves:\t\t" << saveStatCountsVector_[t2].failedLinearSolves_ - saveStatCountsVector_[t1].failedLinearSolves_ << std::endl;

  if (saveStatCountsVector_[t2].linearIters_ > saveStatCountsVector_[t1].linearIters_)
  {
    lout() << "\tNumber Linear Solver Iterations:\t" << saveStatCountsVector_[t2].linearIters_ - saveStatCountsVector_[t1].linearIters_ << std::endl;
  }
  lout() << "\tNumber Residual Evaluations:\t\t" << saveStatCountsVector_[t2].residualEvaluations_ - saveStatCountsVector_[t1].residualEvaluations_ << std::endl
         << "\tNumber Nonlinear Convergence Failures:\t" << saveStatCountsVector_[t2].nonlinearConvergenceFailures_ - saveStatCountsVector_[t1].nonlinearConvergenceFailures_ << std::endl
         << "\tTotal Residual Load Time:\t\t" << saveStatCountsVector_[t2].residualLoadTime_ - saveStatCountsVector_[t1].residualLoadTime_ << " seconds" << std::endl
         << "\tTotal Jacobian Load Time:\t\t" << saveStatCountsVector_[t2].jacobianLoadTime_ - saveStatCountsVector_[t1].jacobianLoadTime_ << " seconds" << std::endl
         << "\tTotal Linear Solution Time:\t\t" << saveStatCountsVector_[t2].linearSolutionTime_ - saveStatCountsVector_[t1].linearSolutionTime_ << " seconds" << std::endl << std::endl;

  return bsuccess;
}


void AnalysisBase::gatherStepStatistics_ ()
{
  if (stepErrorControl_->newtonConvergenceStatus <= 0)
  {
    ++stats_.nonlinearConvergenceFailures_;
  }

  stats_.jacobiansEvaluated_      += nonlinearSolverManager_.getNumJacobianLoads();
  stats_.linearSolves_            += nonlinearSolverManager_.getNumLinearSolves();
  stats_.failedLinearSolves_      += nonlinearSolverManager_.getNumFailedLinearSolves();
  stats_.linearIters_             += nonlinearSolverManager_.getTotalNumLinearIters();
  stats_.residualEvaluations_     += nonlinearSolverManager_.getNumResidualLoads();
  stats_.iterationMatrixFactorizations_ += nonlinearSolverManager_.getNumJacobianFactorizations();
  stats_.linearSolutionTime_            += nonlinearSolverManager_.getTotalLinearSolveTime();
  stats_.residualLoadTime_              += nonlinearSolverManager_.getTotalResidualLoadTime();
  stats_.jacobianLoadTime_              += nonlinearSolverManager_.getTotalJacobianLoadTime();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisBase::firstDoubleDCOPStep_
// Purpose       : If the current step is the first step of
//                  a "doubleDCOP", then return "true".
//
//  Explanation:
//
//  If there are PDE semiconductor devices as part of this problem,
//  there may need to be a "double-pass"
//
//  first pass  = nonlinear poisson solution
//  second pass = drift diffusion solution
//
// Special Notes : Only PDE problems can ever return true.
// Scope         : 
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/21/04
//-----------------------------------------------------------------------------
bool AnalysisBase::firstDoubleDCOPStep_ ()
{
  return (getDoubleDCOPEnabled() && getDoubleDCOPStep() != tiaParams_.lastDCOPStep);
}


//-----------------------------------------------------------------------------
// Function      : setupSweepLoop
// Purpose       : Processes sweep parameters.
// Special Notes : Used for DC and STEP analysis classes.
// Scope         : public
// Creator       : Eric R. Keiter, SNL.
// Creation Date : 08/21/04
//-----------------------------------------------------------------------------
int setupSweepLoop(N_LOA_Loader &loader, int debug_level, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end)
{
  // loop over the param containers, and check that all the params exist.
  // (the device package will complain if it can't find the param)
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    SweepParam &sweep_param = (*it);

    loader.getParamAndReduce(sweep_param.name);
  }

#ifdef Xyce_DEBUG_ANALYSIS
  if (debug_level > 0)
  {
    Xyce::dout() << std::endl << std::endl
                 << Xyce::subsection_divider << std::endl
                 << "AnalysisManager::setupSweepLoop" << std::endl;
  }
#endif

  double pinterval = 1.0;
  double pcount = 0.0, pstart, pstop, pstep;

  // loop over the param containers:
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    SweepParam &sweep_param = (*it);

#ifdef Xyce_DEBUG_ANALYSIS
    if (debug_level > 0)
    {
      Xyce::dout() << "name = " << sweep_param.name << std::endl;
    }
#endif
    // set interval:
    sweep_param.interval = static_cast<int> (pinterval);

    // This stuff should probably be moved up into the SweepParam class.

    // obtain next pinterval:
    if (sweep_param.type=="LIN")
    {
      pstart = sweep_param.startVal;
      pstop  = sweep_param.stopVal;
      pstep  = sweep_param.stepVal;
      // ----------
      // pcount = floor(((pstop - pstart)/pstep) + 1.0);
      // The computation of "pcount" above is notoriously prone to roundoff
      // error, especially if the "floor" function is an in-line function
      // and is subject to high levels of optimization on x86 processors.
      // The symptom is that the last step of a DC sweep (or other sweep)
      // gets lost for very specific combinations of start/stop/step.
      // The next few lines are an attempt to mitigate this roundoff issue,
      // which was present in Xyce for years, and was inherited from SPICE3F5,
      // from which the above expression was taken verbatim.

      // Compute the number of steps of size pstep between pstart and pstop
      pcount = floor(((pstop - pstart)/pstep));
      // Here we're checking that adding one more step doesn't pass pstop
      // by more than machine precision.  If we're within machine precision
      // of pstop by taking one more step, that must mean we undercounted 
      // due to roundoff in the division --- if we hadn't undercounted, we'd 
      // exceed pstop by (nearly) a full pstep.
      if ( fabs(pstop-(pstart+(pcount+1.0)*pstep)) < 2.0*Util::MachineDependentParams::MachinePrecision())
      {
        pcount += 1.0;
      }
      
      // Pcount is now the exact number of steps of size pstep between pstart
      // and pstop, with roundoff handled mostly cleanly.

      // finally, because our actual loop does a loop from zero to maxStep-1,
      // we have to pad maxStep (as is done in the original pcount expression
      // above) to get the full range.
      pcount += 1.0;

      // done this way, we should no longer miss final steps of DC sweeps.
      // Fixed 31 Jul 2012.  This was bug 695 in Bugzilla, and had plagued
      // us since Xyce was first ported to Linux with GCC.
      // ----------

      sweep_param.maxStep = static_cast<int>(pcount);
#ifdef Xyce_DEBUG_ANALYSIS
      if (debug_level > 0)
      {
        Xyce::dout() << "pstart  = " << pstart << std::endl;
        Xyce::dout() << "pstop   = " << pstop  << std::endl;
        Xyce::dout() << "pstep   = " << pstep  << std::endl;
        Xyce::dout() << "pstop-pstart/pstep = " << ((pstop - pstart)/pstep) << std::endl;
        Xyce::dout() << "floor ()= " << floor(((pstop - pstart)/pstep)+1.0) << std::endl;
        Xyce::dout() << "pcount  = " << pcount << std::endl;
        Xyce::dout() << "maxStep = " << sweep_param.maxStep << std::endl;
      }
#endif
    }
    else if(sweep_param.type=="DEC")
    {
      double numSteps = static_cast<double>(sweep_param.numSteps);
      // stepMult could also be calculated as pow(10,(1/numSteps))
      double stepMult = exp(log(10.0)/numSteps);
      sweep_param.stepMult = stepMult;

      pstart   = sweep_param.startVal;
      pstop    = sweep_param.stopVal;
      pcount   = floor(fabs(log10(pstart) - log10(pstop)) * numSteps + 1);
      sweep_param.maxStep = static_cast<int>(pcount);
    }
    else if(sweep_param.type=="OCT")
    {
      double numSteps = static_cast<double>(sweep_param.numSteps);
      // stepMult could also be calculated as pow(2,1/(numSteps))
      double stepMult = exp(log(2.0)/numSteps);

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw 
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2=log(2.0);

      sweep_param.stepMult = stepMult;
      pstart   = sweep_param.startVal;
      pstop    = sweep_param.stopVal;
      pcount   = floor(fabs(log(pstart) - log(pstop))/ln2 * numSteps + 1);
      sweep_param.maxStep = static_cast<int>(pcount);
    }
    else if(sweep_param.type=="LIST")
    {
      pcount = sweep_param.valList.size();
      sweep_param.maxStep = sweep_param.valList.size();
    }
    else
    {
      Report::UserError0() << " Unsupported STEP type";
    }
    pinterval *= pcount;

#ifdef Xyce_DEBUG_ANALYSIS
    if (debug_level > 0)
    {
      Xyce::dout() << "parameter = " << sweep_param.name << std::endl;
      Xyce::dout() << "pcount    = " << pcount << std::endl;
      Xyce::dout() << "pinterval = " << pinterval << std::endl;
    }
#endif
  }

  // At this point, pinterval equals the total number of steps 
  // for the step loop.
  return static_cast<int>(pinterval);
}

//-----------------------------------------------------------------------------
// Function      : updateSweepParams
// Purpose       : Update parameters either for DC or STEP sweeps
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool updateSweepParams(N_LOA_Loader &loader, AnalysisManager &analysis_manager, int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end)
{
  bool resetFlag = false;

  // set parameter(s)
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    (*it).updateCurrentVal(step_count);
    resetFlag = resetFlag || (*it).getSweepResetFlag();
    loader.setParam ((*it).name, (*it).currentVal);
  }

  // Tell the manager if any of our sweeps are being reset in this loop iteration.
  analysis_manager.setSweepSourceResetFlag(resetFlag);

  return true;
}

} // namespace Analysis
} // namespace Xyce

