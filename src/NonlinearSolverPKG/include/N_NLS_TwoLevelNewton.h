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
// Filename       : $RCSfile: N_NLS_TwoLevelNewton.h,v $
//
// Purpose        : Defines N_NLS_TwoLevelNewton class.
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/20/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.81 $
//
// Revision Date  : $Date: 2014/08/07 23:08:54 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_TwoLevelNewton_h
#define Xyce_N_NLS_TwoLevelNewton_h

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_NLS_NLParams.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_TwoLevelEnum.h>

// Loader for RHS and Jacobian
class N_LOA_Loader;

// Linear Algebra Support
class N_LAS_Matrix;
class N_LAS_Vector;
class N_LAS_System;

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : TwoLevelNewton
// Purpose       : This class is the manager for the two level newton
//                 algorithm.  Mostly, it will contain a control loop
//                 which repeatedly calls the "true" nonlinear solver
//                 (either DampedNewton's solve, or NOX, depending
//                 on the options specified).
//
//                 The idea of two level newton is divide the problem
//                 up into sub-problems and solve them separately
//                 for part or all of the solve.
//
//                 The control loop contained in this class repeatedly
//                 loops over all the sub-problems, and also determines
//                 if the total problem has converged.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------

class TwoLevelNewton : public NonLinearSolver
{
public:
  TwoLevelNewton(bool noxFlag, bool noxFlagInner, N_IO_CmdParse & cp);
  ~TwoLevelNewton();

  int getNumIterations() const;
#ifdef Xyce_DEBUG_NONLINEAR
  int getDebugLevel () const;
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

  int solve (NonLinearSolver * nlsTmpPtr=NULL);
  bool setOptions     (const N_UTL_OptionBlock & OB);
  bool setTranOptions (const N_UTL_OptionBlock & OB);
  bool setHBOptions (const N_UTL_OptionBlock & OB);
  bool setTwoLevelOptions     (const N_UTL_OptionBlock & OB);
  bool setTwoLevelTranOptions (const N_UTL_OptionBlock & OB);
  bool setLocaOptions (const N_UTL_OptionBlock& OB);
  bool setTwoLevelLocaOptions (const N_UTL_OptionBlock& OB);
  bool registerLinearSystem(N_LAS_System* ptr);
  bool registerAnalysisManager (N_ANP_AnalysisManager * tiaPtr_tmp);
  bool registerLoader(N_LOA_Loader* ptr);
  bool registerOutputMgr (N_IO_OutputMgr * ptr);
  bool initializeAll ();
  TwoLevelNewtonMode  getCouplingMode ();
  void setAnalysisMode (AnalysisMode mode);
  bool setPetraOptions (const N_UTL_OptionBlock & OB);
  int  getContStepNumber ();
  bool enableSensitivity ();

  // solver statistics:
  virtual int getNumResidualLoads();
  virtual int getNumJacobianLoads();
  virtual int getNumLinearSolves();
  virtual int getNumFailedLinearSolves();
  virtual int getNumJacobianFactorizations();
  virtual unsigned int getTotalNumLinearIters();
  virtual double getTotalLinearSolveTime();
  virtual double getTotalResidualLoadTime();
  virtual double getTotalJacobianLoadTime();

  double getMaxNormF() const;
  int getMaxNormFindex() const;
private:
  TwoLevelNewton();

  void printStepInfo_(int step, int success, TwoLevelNewtonMode solveType);

  void zeroInnerLoopStatistics_ ();
  void calcInnerLoopStatistics_ ();

  void calcOuterLoopStatistics_ ();

  bool calcCouplingTerms_ ();

  int continuationLoop_ ();
  int locaLoop_ ();

  int algorithm0_();
  int algorithm1_();
  int algorithm2_();
  int algorithm3_();
  int algorithm4_();
  int algorithm5_();

public:

private:
  // Pointer to the outer loop nonlinear solver that is being used.
  NonLinearSolver *nlsOuterPtr_;

  // Pointer to the inner loop nonlinear solver that is being used.
  NonLinearSolver *nlsInnerPtr_;

  // Pointer to the time integrator that is being used.
  N_ANP_AnalysisManager *tiaPtr_;

  // number of iterations for the "outer" loop, if it exists. (algorithm>0)
  int maxOuterSteps_;

  // number of iterations for the inner "voltage stepper" loop.
  // This parameter is kind of irrelevant, if running with variable step
  // size.
  int maxContSteps_;
  int maxContStepsTran_;
  int contStep_;

  double increaseContScalar_;
  double decreaseContScalar_;

  // this variable determines which variant of two level newton
  // is being used.
  //
  // Here's the key:
  //
  // 0 = full newton , one level (as though two-level weren't being used)
  //      This is probably the best algorithm for transient.
  //
  // 1 = full newton outter loop, device pde newton inner loop.
  //
  // 2 = full newton outter loop, with pde device inner loop and
  //     continuation.  Sort of a 3-level.
  //
  // 3 = ckt only outter loop, with pde device inner loop, which uses
  //     continuation.  This is sort of a 3-level Newton.  This is
  //     the best algorithm for DCOP.
  //
  int twoLevelAlgorithm_;
  int twoLevelAlgorithmTran_;  // same, but for transient.

  // this flag indicates if the current solve iteration is of the
  // outer loop or the inner loop.
  bool outerLoopActiveFlag_;

  AnalysisMode externalAnalysisMode;

  bool setupOuterLoopParamsFlag_;
  bool setupTranParamsFlag_;
  bool noxFlag_;
  bool noxFlagInner_;

  // inner loop statistics vars:
  int numResidualLoads_;
  int numJacobianLoads_;
  int numLinearSolves_;
  int numFailedLinearSolves_;
  int numJacobianFactorizations_;
  unsigned int totalNumLinearIters_;
  double totalLinearSolveTime_;
  double totalResidualLoadTime_;
  double totalJacobianLoadTime_;

  bool numInterfaceNodesSetup_;

  std::vector<int> numInterfaceNodes_;

  TwoLevelNewtonMode twoLevelCouplingMode_;

  N_LAS_Vector *savedRHSPtr_;
  N_LAS_Vector *savedNextSolPtr_;
  N_LAS_Vector *jdxpVectorPtr_;

  int numSubProblems_;
  int continuationType_;
  bool innerLoopFailFatal_;
  bool totalSolveFailFatal_;
  bool doFullNewtonFinalEnforcement_;

  NonLinearSolver * nlsPassingPtr_;

  bool firstDCOPFlag_;

  bool continuationCalledBefore_;

  // continuation parameters used for algorithm 4.  (not the
  // continuationLoop_ function, that's different)
  std::vector<std::string> paramNameList;
  std::vector<double> paramFinalVal;
  std::vector<double> paramCurrentVal;

  // These options are saved in case they have to be modified
  // over the course of the solve.
  N_UTL_OptionBlock innerSolverOptions_;
  N_UTL_OptionBlock innerLocaOptions_;
  N_UTL_OptionBlock outerLocaOptions_;

  // voltage limiter tolerance:
  double voltLimTol_;
};

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getContStepNumber
// Purpose       : returns the current continuation step number.
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Compuational Sciences
// Creation Date : 10/24/02
//-----------------------------------------------------------------------------
inline int TwoLevelNewton::getContStepNumber ()
{
  return contStep_;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::isFirstSolveComplete
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/22/06
//-----------------------------------------------------------------------------
inline bool TwoLevelNewton::isFirstSolveComplete () const
{
  return true;
}

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getDebugLevel
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
inline int TwoLevelNewton::getDebugLevel () const
{
  return -100;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
inline bool TwoLevelNewton::getScreenOutputFlag () const
{
  return false;
}

//---------------------------------------------------------------------------
// Function      : TwoLevelNewton::getDebugMinTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double TwoLevelNewton::getDebugMinTime() const
{
  return 0.0;
}

//---------------------------------------------------------------------------
// Function      : TwoLevelNewton::getDebugMaxTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double TwoLevelNewton::getDebugMaxTime() const
{
  return N_UTL_MachineDependentParams::DoubleMax();
}

//---------------------------------------------------------------------------
// Function      : TwoLevelNewton::getDebugMinTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int TwoLevelNewton::getDebugMinTimeStep() const
{
  return 0;
}

//---------------------------------------------------------------------------
// Function      : TwoLevelNewton::getDebugMaxTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int TwoLevelNewton::getDebugMaxTimeStep() const
{
  return N_UTL_MachineDependentParams::IntMax();
}

//---------------------------------------------------------------------------
// Function      : TwoLevelNewton::getMMFormat
//
// Return Type   : bool
//---------------------------------------------------------------------------
inline bool TwoLevelNewton::getMMFormat () const
{
  return false;
}

#endif // debug nonlin

} // namespace Nonlinear
} // namespace Xyce

typedef Xyce::Nonlinear::TwoLevelNewton N_NLS_TwoLevelNewton;

#endif
