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
// Filename       : $RCSfile: N_LOA_Loader.h,v $
//
// Purpose        : This file contains class definitions for the base
//                  loader class.
//
// Special Notes  : This is a potentially confusing set of classes.  As of
//                  this writing (2/11/2007), loaders are used in 2 different
//                  interface layers:
//
//                  1) between the nonlinear solver and time integrator:
//                      Nonlinear Equation loader.
//
//                  2) between the time integrator and the device package:
//                      Ckt loader, MPDE, and SawTooth loaders
//
//                  This gets confusing, because 2 different purposes are being
//                  served.
//
//                  The Nonlinear Equation loader exists to insulate the nonlinear
//                  solver from having to know what kind of problem (in a math
//                  sense) this is.  For DCOP, the rhs consists of F(x) and B(t),
//                  but for transient F(x,t) , dQdt(x,t) and B(t). The Nonlinear Equation
//                  loader hides this decision from the nonlinear solver.
//
//                  The other layer (Ckt, MPDE and sawtooth) does 2 things.
//                  (1) it insulates the solvers from the specific physics.
//                  ie, the solvers (theoretically) don't know this is a
//                  circuit code. (2) We can hide MPDE details behind this
//                  interface, so the rest of the code can mostly be unaware
//                  if we are running an MPDE simulation or not.
//
//                  This structure is not perfect, and because it serves several
//                  different purposes, confusing.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.95.2.1 $
//
// Revision Date  : $Date: 2014/09/02 22:49:48 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_Loader_H
#define Xyce_LOA_Loader_H

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>
#include <N_DEV_fwd.h>
#include <N_ANP_fwd.h>

// ----------   Forward Declarations  ----------
class N_LAS_Vector;
class N_LAS_Matrix;
class N_TIA_TwoLevelError;

//-----------------------------------------------------------------------------
// Class         : N_LOA_Loader
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
class N_LOA_Loader
{
public:

  // Defat constructor.
  N_LOA_Loader();

  // Destructor.
  virtual ~N_LOA_Loader();

  virtual bool loadJacobian () { return true; };

  // Virtual nonlinear Jacobian apply method.
  virtual bool applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result)
  { return false; }

  virtual bool loadRHS () { return true; }

  virtual bool loadSensitivityResiduals () { return true; }

  virtual bool loadFinalSensitivityDerivatives () { return true; }

  // return value indicates whether mask is nontrivial
  // base class method does nothing.
  virtual bool loadDeviceMask() {return false;}

  // Virtual nonlinear new-DAE Matrices  load method.
  virtual bool loadDAEMatrices  (N_LAS_Vector * tmpSolVectorPtr,
                                 N_LAS_Vector * tmpStaVectorPtr,
                                 N_LAS_Vector * tmpStaDerivVectorPtr,
                                 N_LAS_Vector * tmpStoVectorPtr,
                                 N_LAS_Matrix * tmpdQdxMatrixPtr,
                                 N_LAS_Matrix * tmpdFdxMatrixPtr)
  { return false; }

  // Virtual nonlinear new-DAE vectors load method.
  virtual bool loadDAEVectors   (N_LAS_Vector * nextSolVectorPtr,
                                 N_LAS_Vector * currSolVectorPtr,
                                 N_LAS_Vector * lastSolVectorPtr,
                                 N_LAS_Vector * nextStaVectorPtr,
                                 N_LAS_Vector * currStaVectorPtr,
                                 N_LAS_Vector * lastStaVectorPtr,
                                 N_LAS_Vector * StaDerivVectorPtr,
                                 N_LAS_Vector * nextStoVectorPtr,
                                 N_LAS_Vector * currStoVectorPtr,
                                 N_LAS_Vector * lastStoVectorPtr,
                                 N_LAS_Vector * stoLeadCurrQVectorPtr,
                                 N_LAS_Vector * QVectorPtr,
                                 N_LAS_Vector * FVectorPtr,
                                 N_LAS_Vector * BVectorPtr,
                                 N_LAS_Vector * dFdxdVpVectorPtr,
                                 N_LAS_Vector * dQdxdVpVectorPtr)
  { return false; }

  // Virtual nonlinear new-DAE Matrices  apply method.
  // tmpdQdxVecVector = tmpdQdxMatrix * tmpVecVector
  // tmpdFdxVecVector = tmpdFdxMatrix * tmpVecVector
  virtual bool applyDAEMatrices  (N_LAS_Vector * tmpSolVectorPtr,
                                  N_LAS_Vector * tmpStaVectorPtr,
                                  N_LAS_Vector * tmpStaDerivVectorPtr,
                                  N_LAS_Vector * tmpStoVectorPtr,
                                  const N_LAS_Vector & tmpVecVectorPtr,
                                  N_LAS_Vector * tmpdQdxVecVectorPtr,
                                  N_LAS_Vector * tmpdFdxVecVectorPtr)
  { return false; }

  virtual bool updateState      (N_LAS_Vector * nextSolVectorPtr,
                                 N_LAS_Vector * currSolVectorPtr,
                                 N_LAS_Vector * lastSolVectorPtr,
                                 N_LAS_Vector * nextStaVectorPtr,
                                 N_LAS_Vector * currStaVectorPtr,
                                 N_LAS_Vector * lastStaVectorPtr,
                                 N_LAS_Vector * nextStoVectorPtr,
                                 N_LAS_Vector * currStoVectorPtr,
                                 N_LAS_Vector * lastStoVectorPtr
                                 )
  { return false; }

  virtual bool loadBVectorsforAC (N_LAS_Vector * bVecRealPtr,
                                  N_LAS_Vector * bVecImagPtr)
  { return false; }

  virtual bool getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec, std::vector<int>& bMatPosEntriesVec)
  { return false; }

//  virtual bool getInductorsEntriesforMOR(std::vector<int>& inductorEntriesVec)
//  { return false; }

  // Virtual function for setting the initial guess.
  virtual bool setInitialGuess(N_LAS_Vector * solVectorPtr)
  { return false; }

  // Virtual function for determining if an analytic sensitivity (df/dp or dqdp) is available.
  virtual bool analyticSensitivitiesAvailable (std::string & name) = 0;

  virtual void getAnalyticSensitivities
    (std::string & name, 
     std::vector<double> & dfdpVec, 
     std::vector<double> & dqdpVec,
     std::vector<double> & dbdpVec,
     std::vector<int> & FindicesVec,
     std::vector<int> & QindicesVec,
     std::vector<int> & BindicesVec
     ) = 0;

  // Virtual function for setting a single parameter value.
  virtual bool setParam (std::string & name, double val) = 0; 

  // Virtual function for getting a single parameter value.
  virtual double getParamAndReduce (const std::string & name) = 0; 
  virtual bool getParamAndReduce (const std::string & name, double & val) = 0; 

  int getVoltageDropRow (std::string & srcName) { return -1; }
  virtual int getVposRow (std::string & srcName) { return -1; }

  // Virtual method which is called to update the sources.
  virtual bool updateSources() { return false; }

  // Virtual method which initializes the nonlinear problem.
  virtual bool initializeProblem( N_LAS_Vector * nextSolVectorPtr,
                          N_LAS_Vector * currSolVectorPtr,
                          N_LAS_Vector * lastSolVectorPtr,
                          N_LAS_Vector * nextStaVectorPtr,
                          N_LAS_Vector * currStaVectorPtr,
                          N_LAS_Vector * lastStaVectorPtr,
                          N_LAS_Vector * StateDerivVectorPtr,
                          N_LAS_Vector * nextStoVectorPtr,
                          N_LAS_Vector * currStoVectorPtr,
                          N_LAS_Vector * lastStoVectorPtr,
                          N_LAS_Vector * QVectorPtr,
                          N_LAS_Vector * FVectorPtr,
                          N_LAS_Vector * BVectorPtr,
                          N_LAS_Vector * dFdxdVpVectorPtr,
                          N_LAS_Vector * dQdxdVpVectorPtr)  = 0;
  //{ return false; }

  virtual bool getLinearSystemFlag() { return false; }

  // This is somewhat circuit-specific, unfortunately.
  virtual bool getLimiterFlag() { return false; }

  // Virtual method which gets the double DC Operating Point flag - used for
  // PDE devices.
  virtual bool getDoubleDCOPFlag() { return false; }

  virtual int  enablePDEContinuation () {return -1;}
  virtual bool disablePDEContinuation () {return false;}

  // functions related to the two-level newton:
  virtual void getNumInterfaceNodes (std::vector<int> & numINodes) { return; }
  virtual bool loadCouplingRHS (int iSubProblem, int iCouple, N_LAS_Vector * dfdvPtr) { return false; }
  virtual bool calcCouplingTerms (int iSubProblem, int iCouple, const N_LAS_Vector * dxdvPtr) { return false; }
  virtual bool raiseDebugLevel (int increment){return false;}

  virtual bool output();
  virtual bool finishOutput();

  // Accessors for timing info.

  // Virtual method which gets the nonlinear residual load time.
  virtual double getResidualTime() { return -1.0; }

  // Virtual method which gets the nonlinear Jacobian load time.
  virtual double getJacobianTime() { return -1.0; }

  // Virtual method which gets the time integration required breakpoint times
  // (in a vector).
  virtual bool getBreakPoints(std::vector< N_UTL_BreakPoint > & breakPointTimes)
  { return false; }

  // Virtual accessor which returns the maximum time step size (in seconds).
  virtual double getMaxTimeStepSize() { return 0.0; }

  // Get block size from device options for block gainscale homotopy
  virtual int getHomotopyBlockSize() const { return 1; }

  // Get convergence info from devices
  virtual bool allDevsConverged() { return true; };

  // Get convergence info from inner-solves
  virtual bool innerDevsConverged() { return true; };

  // Functions needed by the NEW (power node) 2-level algorithm:
  virtual void homotopyStepSuccess
    (const std::vector<std::string> & paramNames,
     const std::vector<double> & paramVals) {}

  virtual void homotopyStepFailure () {}

  virtual void stepSuccess(Xyce::Analysis::CurrentMode analysis) {}
  virtual void stepFailure(Xyce::Analysis::CurrentMode analysis) {}


  //TVR: Method to be called when time integrator accepts a step, before any
  // tinkering with times or vectors
  virtual void acceptStep() {}

  virtual bool getInitialQnorm (std::vector<N_TIA_TwoLevelError> & tleVec )
  {return false;}
  virtual bool getInnerLoopErrorSums (std::vector<N_TIA_TwoLevelError> & tleVec )
  {return false;};

  virtual bool updateStateArrays () {return true;}
  virtual bool startTimeStep () {return true;}

  virtual void setExternalSolverState (const N_DEV_SolverState & ss)
  {}

protected:
private :

};

#endif
