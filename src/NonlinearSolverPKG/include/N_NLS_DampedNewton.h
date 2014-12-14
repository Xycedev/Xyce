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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_DampedNewton.h,v $
//
// Purpose        : Specification file for the implemenation of the Newton
//                  trust-region related methods.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences.
//
// Creation Date  : 04/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.77 $
//
// Revision Date  : $Date: 2014/08/07 23:08:54 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_DampedNewton_h
#define Xyce_N_NLS_DampedNewton_h

#include <N_IO_fwd.h>
#include <N_NLS_fwd.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_NLParams.h>
#include <N_NLS_ParamMgr.h>
#include <N_NLS_ReturnCodes.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : DampedNewton
// Purpose       : This class implements a damped Newton nonlinear solver.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
class DampedNewton : public NonLinearSolver
{
public:

   DampedNewton(N_IO_CmdParse & cp);
   ~DampedNewton();

  bool setOptions(const N_UTL_OptionBlock& OB);
  bool setTranOptions(const N_UTL_OptionBlock& OB);
  bool setHBOptions(const N_UTL_OptionBlock& OB);
  
  bool initializeAll();

  int solve (NonLinearSolver * nlsTmpPtr = NULL);

  int takeFirstSolveStep (NonLinearSolver * nlsTmpPtr = NULL);
  int takeOneSolveStep   ();

  int getNumIterations() const;
#ifdef Xyce_DEBUG_NONLINEAR
  int getDebugLevel() const;
  bool getScreenOutputFlag () const;
  double getDebugMinTime() const;
  double getDebugMaxTime() const;
  int getDebugMinTimeStep() const;
  int getDebugMaxTimeStep() const;
  bool getMMFormat () const;
#endif

  int getContinuationStep() const;
  int getParameterNumber() const;

  bool isFirstContinuationParam() const;
  bool isFirstSolveComplete() const;

  void setAnalysisMode(AnalysisMode mode);

  double getMaxNormF() const
  { return maxNormRHS_; };

  int getMaxNormFindex() const
  { return maxNormRHSindex_; };

protected:

private:

  void updateWeights_();

  void printHeader_(std::ostream &os);
  void printFooter_(std::ostream &os);
  void printStepInfo_(std::ostream &os, int step);

  bool rhs_();
  bool newton_();

  void direction_();
  void updateX_();
  bool computeStepLength_();

  bool   divide_();
  bool   backtrack_();
  bool   fullNewton_();
  bool   spiceNewton_();
  bool   bankRose_();
  bool   descent_();
  bool   simpleBacktrack_();
  bool   simpleBt_(double gsinit, double finit);

  double constrain_();
  void   setForcing_(const double);
  void   evalModNewton_();
  int    converged_();

  void resetCountersAndTimers_();

public:

protected:

private:
  //! Convergence Rate
  double resConvRate_;

  //! Weighted convergence Rate
  double wtUpdateConvRate_;

  //! Constraint object pointer
  ConstraintBT* nlConstraintPtr_;

  // Current RHS norms:
  double normRHS_;
  double maxNormRHS_;
  int maxNormRHSindex_;

  double normRHS_init_;  // used in calculating the relative norm.

  // Current dx norm:
  double normDX_;

  // Current weighted dx norm:
  double wtNormDX_;

  // Current relative RHS norm:
  double normRHS_rel_;

  // Current solution norm:
  double normSoln_;

  // step length:
  double stepLength_;

  // backtracking upper and lower bounds and percentage change (0..1):
  double BTUpper_;
  double BTLower_;

  // constraint factor:
  double constraintFactor_;

  // current nonlinear solver step:
  unsigned nlStep_;

  // current Newton step:
  unsigned newtonStep_;

  // current modified Newton step:
  unsigned modNewtonStep_;

  // current steepest descent step:
  unsigned descentStep_;

  // current line-search step:
  unsigned searchStep_;

  //! Pointer to direction vector
  /*! \todo Is this vector internal or external?? Is it allocated or just a pointer?? */
  N_LAS_Vector* searchDirectionPtr_;

  // Needed for some reason
  N_LAS_Vector* tmpVectorPtr_;

  int iNumCalls_;

  double delta_;

  // Flag for controlling the loading of the Jacobian matrix.
  bool loadJacobianFlag_;

  // Flag for determining if this solver has been called
  // before, and deltaXTol.  The first time the solver is
  // called we may want to tighten up the convergence tolerance
  // a la Petzold, et al.  These were static variables
  // down in the "solve" function of this class.  ERK.
  bool firstTime;
  double initialDeltaXTol;

  NLParams nlParams;

  double etaOld;
  double nlResNormOld;
  double tmpConvRate;

  int count;
};

//---------------------------------------------------------------------------
// Function      : DampedNewton::getNumIterations
//
// Return Type   : Integer (current number of nonlinear iterations)
//---------------------------------------------------------------------------
inline int DampedNewton::getNumIterations() const
{
  return nlStep_;
}

#ifdef Xyce_DEBUG_NONLINEAR
//---------------------------------------------------------------------------
// Function      : DampedNewton::getDebugLevel
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int DampedNewton::getDebugLevel() const
{
  return nlParams.getDebugLevel();
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::getScreenOutputFlag 
//
// Return Type   : int
//---------------------------------------------------------------------------
inline bool DampedNewton::getScreenOutputFlag () const
{
  return nlParams.getScreenOutputFlag ();
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::getDebugMinTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double DampedNewton::getDebugMinTime() const
{
  return nlParams.getDebugMinTime();
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::getDebugMaxTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double DampedNewton::getDebugMaxTime() const
{
  return nlParams.getDebugMaxTime();
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::getDebugMinTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int DampedNewton::getDebugMinTimeStep() const
{
  return nlParams.getDebugMinTimeStep();
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::getDebugMaxTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int DampedNewton::getDebugMaxTimeStep() const
{
  return nlParams.getDebugMaxTimeStep();
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::getMMFormat
//
// Return Type   : bool
//---------------------------------------------------------------------------
inline bool DampedNewton::getMMFormat () const
{
  return nlParams.getMMFormat ();
}

#endif // debug nonlin

//---------------------------------------------------------------------------
// Function      : DampedNewton::setAnalysisMode
//
// Purpose       : Specify the analysis mode to be used by the nonlinear
//                 solver in the next call to solve(). This *may* affect
//                 the parameters used by the solver. 
//
// See Also      : setOptions, setTranOptions
//
// - Input Arguments -
//
//    mode       : Mode to be used in the next nonlinear solve.
//---------------------------------------------------------------------------
inline void DampedNewton::setAnalysisMode(AnalysisMode mode)
{
  nlpMgrPtr_->setAnalysisMode(mode);
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::resetCountersAndTimers_
// Purpose       :  Reset the counters and timers
//---------------------------------------------------------------------------
inline void DampedNewton::resetCountersAndTimers_()
{
  NonLinearSolver::resetCountersAndTimers_();
  resConvRate_ = 0.0;
  wtUpdateConvRate_ = 1.0;
}

//---------------------------------------------------------------------------
// Function      : DampedNewton::resetCountersAndTimers_
// Purpose       :  Reset the counters and timers
//---------------------------------------------------------------------------
inline bool DampedNewton::isFirstSolveComplete() const
{
  return true;
}

} // namespace Nonlinear
} // namespace Xyce

typedef Xyce::Nonlinear::DampedNewton N_NLS_DampedNewton;

#endif

