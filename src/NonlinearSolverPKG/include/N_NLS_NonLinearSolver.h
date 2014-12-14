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
// Filename       : $RCSfile: N_NLS_NonLinearSolver.h,v $
//
// Purpose        : Specification file which declares an interface common to
//                  all supported nonlinear solver algorithms.  The Manager
//                  class uses this interface to call a concrete algorithm.
//
// Special Notes  : This is the "Strategy" class in the Strategy design
//                  pattern.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.118 $
//
// Revision Date  : $Date: 2014/08/07 23:08:54 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NonLinearSolver_h
#define Xyce_N_NLS_NonLinearSolver_h

// ---------- Standard Includes ----------

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_NLS_Manager.h> // for definition of AnalysisMode
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_NLS_fwd.h>
#include <N_NLS_ReturnCodes.h>
#include <N_NLS_NonLinInfo.h>
#include <N_TOP_Topology.h>

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_Matrix;
class N_LAS_System;

class N_LAS_Solver;
class N_LAS_Problem;

class N_LOA_Loader;

class N_LAS_PrecondFactory;

class N_PDS_Manager;

// ---------- Using Declarations ----------
using Teuchos::RefCountPtr;

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : NonLinearSolver
// Purpose       : Nonlinear Solver Abstract Class
// Special Notes : Many of the virtual functions should not, in general,
//                 be redefined in derived classes. Check the Virtual Notes
//                 on each function.
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
class NonLinearSolver
{

public:

  NonLinearSolver(N_IO_CmdParse & cp);
  virtual ~NonLinearSolver();

  virtual bool setOptions(const N_UTL_OptionBlock& OB) = 0;
  virtual bool setTranOptions(const N_UTL_OptionBlock& OB) = 0;
  virtual bool setHBOptions(const N_UTL_OptionBlock& OB) = 0;
  virtual bool setLocaOptions(const N_UTL_OptionBlock& OB);
  virtual bool setTwoLevelLocaOptions(const N_UTL_OptionBlock& OB);
  virtual bool setTwoLevelOptions    (const N_UTL_OptionBlock& OB);
  virtual bool setTwoLevelTranOptions(const N_UTL_OptionBlock& OB);
  virtual bool setPetraOptions(const N_UTL_OptionBlock& OB);
  virtual bool setDCOPRestartOptions(const N_UTL_OptionBlock& OB);
  virtual bool setICOptions(const N_UTL_OptionBlock& OB);
  virtual bool setNodeSetOptions(const N_UTL_OptionBlock& OB);

  virtual bool registerRHSVector(N_LAS_Vector* ptr);
  virtual bool registerLoader(N_LOA_Loader* ptr);
  virtual bool registerLinearSystem(N_LAS_System* ptr);
  virtual bool registerTwoLevelSolver (TwoLevelNewton * ptr);
  virtual bool registerParamMgr (ParamMgr * ptr);
  virtual bool registerTopology(N_TOP_Topology * ptr);
  virtual bool registerPrecondFactory(const RefCountPtr<N_LAS_PrecondFactory>& ptr);
  virtual bool registerParallelMgr(N_PDS_Manager * pdsMgrPtr);
  virtual bool registerAnalysisManager(N_ANP_AnalysisManager* tmp_anaIntPtr);
  virtual bool registerOutputMgr (N_IO_OutputMgr * outPtr);
  virtual bool registerTIADataStore(N_TIA_DataStore * tiaDSPtr);

  virtual bool initializeAll();

  virtual int solve (NonLinearSolver * nlsTmpPtr = NULL) = 0;
  virtual inline int takeFirstSolveStep (NonLinearSolver * nlsTmpPtr = NULL);
  virtual inline int takeOneSolveStep ();

  virtual int getNumIterations() const = 0;
#ifdef Xyce_DEBUG_NONLINEAR
  virtual int getDebugLevel() const = 0;
  virtual bool getScreenOutputFlag () const=0;
  virtual double getDebugMinTime() const = 0;
  virtual double getDebugMaxTime() const = 0;
  virtual int getDebugMinTimeStep() const = 0;
  virtual int getDebugMaxTimeStep() const = 0;
  virtual bool getMMFormat () const = 0;
#endif

  virtual bool isFirstContinuationParam() const = 0;
  virtual bool isFirstSolveComplete() const = 0;
  virtual int getContinuationStep() const = 0;
  virtual int getParameterNumber() const = 0;

  virtual bool getLocaFlag ();

  virtual inline int getNumResidualLoads();
  virtual inline int getNumJacobianLoads();
  virtual inline int getNumLinearSolves();
  virtual inline int getNumFailedLinearSolves();
  virtual inline int getNumJacobianFactorizations();
  virtual inline unsigned int getTotalNumLinearIters();
  virtual inline double getTotalLinearSolveTime();
  virtual inline double getTotalResidualLoadTime();
  virtual inline double getTotalJacobianLoadTime();

  virtual TwoLevelNewtonMode getCouplingMode ();

  virtual void setAnalysisMode(AnalysisMode mode) = 0;
  virtual void resetAll (AnalysisMode mode);
  virtual void setReturnCodes (const ReturnCodes & retCodesTmp);
  virtual bool enableSensitivity () {return true;}
  virtual bool getMatrixFreeFlag();
  virtual void setMatrixFreeFlag(bool matrixFreeFlag);

  virtual double getMaxNormF() const = 0;
  virtual int getMaxNormFindex () const = 0;

#ifdef Xyce_DEBUG_NONLINEAR
  // use for debugging:
  void debugOutput1 (N_LAS_Matrix & jacobian, N_LAS_Vector & rhs);
  void debugOutput3 (N_LAS_Vector & dxVector, N_LAS_Vector & xVector);

  void debugOutputDAE();

#ifdef Xyce_DEBUG_VOLTLIM
  void debugOutputJDX_VOLTLIM ();
#endif

  void  setDebugFlags ();
#endif

  virtual bool applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result);

protected:

  virtual void resetCountersAndTimers_();
  virtual bool setX0_();
  virtual bool rhs_();
  virtual bool jacobian_();
  virtual bool newton_();
  virtual bool gradient_();

protected:

  std::string netlistFileName_;
  N_LAS_Vector** nextSolVectorPtrPtr_;
  N_LAS_Vector** currSolVectorPtrPtr_;
  N_LAS_Vector** tmpSolVectorPtrPtr_;
  N_LAS_Vector* rhsVectorPtr_;

#ifdef Xyce_DEBUG_VOLTLIM
  N_LAS_Matrix* jacTestMatrixPtr_;
  N_LAS_Matrix* dFdxTestMatrixPtr_;
  N_LAS_Matrix* dQdxTestMatrixPtr_;
  N_LAS_Vector* dxVoltlimVectorPtr_;
  N_LAS_Vector* jdxVLVectorPtr_; // old-DAE
  N_LAS_Vector* fdxVLVectorPtr_; // new-DAE
  N_LAS_Vector* qdxVLVectorPtr_; // new-DAE
#endif

  N_LAS_Matrix* jacobianMatrixPtr_;
  N_LAS_Vector* gradVectorPtr_;
  N_LAS_Vector* NewtonVectorPtr_;
  N_LAS_Vector* solWtVectorPtr_;
  N_LAS_System* lasSysPtr_;
  N_LAS_Solver * lasSolverPtr_;
  RefCountPtr<N_LAS_Problem> lasProblemRCPtr_;
  RefCountPtr<N_LAS_PrecondFactory> lasPrecPtr_;
  N_UTL_OptionBlock* petraOptionBlockPtr_;
  N_LOA_Loader* loaderPtr_;
  N_ANP_AnalysisManager* anaIntPtr_;
  TwoLevelNewton * tlnPtr_;
  ParamMgr * nlpMgrPtr_;
  N_IO_OutputMgr * outMgrPtr_;
  Teuchos::RefCountPtr<N_TOP_Topology> topologyRcp_;
  N_PDS_Manager * pdsMgrPtr_;
  N_TIA_DataStore * dsPtr_;

  int numJacobianLoads_;
  int numJacobianFactorizations_;
  int numLinearSolves_;
  int numFailedLinearSolves_;
  int numResidualLoads_;
  unsigned int totalNumLinearIters_;
  double totalLinearSolveTime_;
  double totalResidualLoadTime_;
  double totalJacobianLoadTime_;
  ReturnCodes retCodes_;
  bool matrixFreeFlag_;

  N_IO_CmdParse & commandLine_;

  friend class ConductanceExtractor;
  friend class Sensitivity;
  friend class TwoLevelNewton;
  friend class Manager;

  int outputStepNumber_;  // this is either the time step number or the dc
                          // sweep step number, depending on the mode.   It
			  // is only used in setting up output file names.
  bool debugTimeFlag_;

  int contStep_;
};

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumResidualLoads
// Return Type   : Integer (Get the total number of residual loads)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumResidualLoads()
{
  return numResidualLoads_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumJacobianLoads
// Return Type   : Integer (Get the total number of Jacobian loads)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumJacobianLoads()
{
  return numJacobianLoads_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumLinearSolves
// Return Type   : Integer (total number of successful linear solves)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumLinearSolves()
{
  return numLinearSolves_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumFailedLinearSolves
// Return Type   : Integer (total number of failed linear solves)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumFailedLinearSolves()
{
  return numFailedLinearSolves_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumJacobianFactorizations
// Return Type   : Integer (total number of Jacobian factorizations)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumJacobianFactorizations()
{
  return numJacobianFactorizations_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalNumLinearIters
// Return Type   : unsigned int (total number of iterative linear solver
//                 iterations)
//---------------------------------------------------------------------------
inline unsigned int NonLinearSolver::getTotalNumLinearIters()
{
  return totalNumLinearIters_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalLinearSolveTime
// Return Type   : double (total linear solve time in seconds)
//---------------------------------------------------------------------------
inline double NonLinearSolver::getTotalLinearSolveTime()
{
  return totalLinearSolveTime_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalResidualLoadTime
// Return Type   : double (total residual load (claculation) time in seconds)
//---------------------------------------------------------------------------
inline double NonLinearSolver::getTotalResidualLoadTime()
{
  return totalResidualLoadTime_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalJacobianLoadTime
// Return Type   : double (total Jacobian load (claculation) time in seconds)
//---------------------------------------------------------------------------
inline double NonLinearSolver::getTotalJacobianLoadTime()
{
  return totalJacobianLoadTime_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::takeFirstSolveStep
// Return Type   : int
//---------------------------------------------------------------------------
inline int NonLinearSolver::takeFirstSolveStep
  (NonLinearSolver * nlsTmpPtr)
{
  return -1;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::takeOneSolveStep
// Return Type   : int
//---------------------------------------------------------------------------
inline int NonLinearSolver::takeOneSolveStep ()
{
  return -1;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::resetAll
// Return Type   : void
//---------------------------------------------------------------------------
inline void NonLinearSolver::resetAll (AnalysisMode mode)
{
  setAnalysisMode(mode);
}
//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setReturnCodes
// Purpose       : Allows the user to set return codes.
//
// Special Notes : This was put in place mainly to help with continuation
//                 solves.  If running continuation, I don't want the
//                 solver to consider "nearConvergence" to be a success.
//                 There are other circumstances that I may wish to turn
//                 this off as well.  For example, if running with the
//                 predictor/corrector step error control turned off (this
//                 is an option), I also would not want "nearConverged" to
//                 be good enough.
//
//                 The rational for having a "nearConverged" option  is the
//                 idea that we can just let the time integrator handle
//                 things and let it use the predictor/corrector stuff to
//                 determine whether or not to reject the step.  If that
//                 stuff is turned off, you otherwise unavailable, you
//                 don't want to do this.
//
//                 This way, the calling code can, if it wants, decide for
//                 itself what senarios will be considered a solver success.
//
// Scope         : public
// Creator       : Eric R. Keiter, 9233
// Creation Date : 10/30/04
//-----------------------------------------------------------------------------
inline void NonLinearSolver::setReturnCodes
  (const ReturnCodes & retCodesTmp)
{
  retCodes_ = retCodesTmp;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::getLocaFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/09/05
//-----------------------------------------------------------------------------
inline bool NonLinearSolver::getLocaFlag ()
{
  return false;
}

} // namespace Nonlinear
} // namespace Xyce

typedef Xyce::Nonlinear::NonLinearSolver N_NLS_NonLinearSolver;

#endif

