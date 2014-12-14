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
// Filename       : $RCSfile: N_LOA_NonlinearEquationLoader.C,v $
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11.2.3 $
//
// Revision Date  : $Date: 2014/09/02 22:49:48 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Xyce Includes   ----------
#include <N_TIA_TimeIntegrationMethods.h>
#include <N_TIA_DataStore.h>

#include <N_LOA_NonlinearEquationLoader.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Timer.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::N_LOA_NonlinearEquationLoader
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
N_LOA_NonlinearEquationLoader::N_LOA_NonlinearEquationLoader( N_TIA_DataStore & ds,
                                   N_LOA_Loader & loader,
                                   N_TIA_WorkingIntegrationMethod & wim,
                                   N_PDS_Manager & pds,
                                   bool daeStateDerivFlag
                                   )
: ds_(ds),
  loader_(loader),
  wim_(wim),
  pdsMgr(pds),
  residualTimerPtr_(0),
  jacobianTimerPtr_(0),
  daeStateDerivFlag_(daeStateDerivFlag)
{
  residualTimerPtr_ = new N_UTL_Timer(*(pdsMgr.getPDSComm()));
  jacobianTimerPtr_ = new N_UTL_Timer(*(pdsMgr.getPDSComm()));
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::~N_LOA_NonlinearEquationLoader
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
N_LOA_NonlinearEquationLoader::~N_LOA_NonlinearEquationLoader()
{
  if (residualTimerPtr_ != 0)  delete residualTimerPtr_;
  if (jacobianTimerPtr_ != 0)  delete jacobianTimerPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::loadRHS
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the residual (RHS).
//
// Special Notes : All the contributions to the RHS come (in the new DAE
//                 form) from the device package, as Q, F, and B.  The 
//                 RHS needs dQdt + F - B.  As dQdt is determined by the
//                 time integration package, the final summation should be
//                 managed from here.
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 03/04/04
//-----------------------------------------------------------------------------
bool N_LOA_NonlinearEquationLoader::loadRHS ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  residualTimerPtr_->resetStartTime();

  ds_.daeQVectorPtr->putScalar(0.0);
  ds_.daeFVectorPtr->putScalar(0.0);
  ds_.daeBVectorPtr->putScalar(0.0);

  ds_.dFdxdVpVectorPtr->putScalar(0.0);
  ds_.dQdxdVpVectorPtr->putScalar(0.0);

  // Update the state. Note - underneath this call, most of the calculations
  // pertaining to the currents, conductances, etc. will happen.
  tmpBool = loader_.updateState 
              ((ds_.nextSolutionPtr), 
               (ds_.currSolutionPtr), 
               (ds_.lastSolutionPtr), 
               (ds_.nextStatePtr),
               (ds_.currStatePtr),
               (ds_.lastStatePtr),
               (ds_.nextStorePtr),
               (ds_.currStorePtr),
               (ds_.lastStorePtr)
               );

  bsuccess = bsuccess && tmpBool;

  if (daeStateDerivFlag_)
  {
    wim_.updateStateDeriv ();
  }

  // first load the 2 components: Q, F and B
  tmpBool = loader_.loadDAEVectors   
              ((ds_.nextSolutionPtr), 
               (ds_.currSolutionPtr), 
               (ds_.lastSolutionPtr), 
               (ds_.nextStatePtr),
               (ds_.currStatePtr),
               (ds_.lastStatePtr),
               (ds_.nextStateDerivPtr),
               (ds_.nextStorePtr),
               (ds_.currStorePtr),
               (ds_.lastStorePtr),
               (ds_.nextStoreLeadCurrQPtr),
               (ds_.daeQVectorPtr),
               (ds_.daeFVectorPtr),
               (ds_.daeBVectorPtr),
               (ds_.dFdxdVpVectorPtr),
               (ds_.dQdxdVpVectorPtr) );
  bsuccess = bsuccess && tmpBool;
  wim_.updateLeadCurrent();
  // Now determine dQdt:
  // now sum them all together, to create the total.
  // f(x) is given by:
  //
  //    f(x) = dQ/dt + F(x) = 0
  //
  // if running with SEPARATE_F_AND_B then instead:
  //    f(x) = dQ/dt + F(x) - B(t)= 0
  //
  // Note, the nonlinear solver is expecting the RHS vector to
  // contain -f(x).  Or, possibly -f(x) + J*dx, if voltage 
  // limiting is on.
  wim_.obtainResidual();

  // Update the total load time
  residualTime_ = residualTimerPtr_->elapsedTime();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::loadJacobian
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the Jacobian.
//
// Special Notes : All the contributions to the Jacobian 
//                 come (in the new DAE form)
//                 from the device package, as dQdx and dFdx. 
//
//                 The Jacobian is: 
//
//                 J = df/dx = d(dQdt)/dx + dF/dx 
//
//                 As dQdt is determined by the time integration package, 
//                 the final summation should be managed from here.
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 03/04/04
//-----------------------------------------------------------------------------
bool N_LOA_NonlinearEquationLoader::loadJacobian ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  jacobianTimerPtr_->resetStartTime();

  ds_.dQdxMatrixPtr->put(0.0);
  ds_.dFdxMatrixPtr->put(0.0);

  // first load the 2 components: dQdx and dFdx
  tmpBool = loader_.loadDAEMatrices 
              ((ds_.nextSolutionPtr),
               (ds_.nextStatePtr),
               (ds_.nextStateDerivPtr),
               (ds_.nextStorePtr),
               (ds_.dQdxMatrixPtr),
               (ds_.dFdxMatrixPtr));
  bsuccess = bsuccess && tmpBool;

  // Now determine the d(dQdt)/dx stuff:

  // now sum them all together, to create the total:
  // J = alpha/dt * dQ/dx + dF/dx
  wim_.obtainJacobian();

  // Update the total load time
  jacobianTime_ = jacobianTimerPtr_->elapsedTime();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::applyJacobian
//
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool N_LOA_NonlinearEquationLoader::applyJacobian 
  (const N_LAS_Vector& input, N_LAS_Vector& result)
{
  bool bsuccess = true;
  bool tmpBool = true;

  // Start the timer...
  jacobianTimerPtr_->resetStartTime();

  // first load the 2 components: dQdx and dFdx
  tmpBool = loader_.applyDAEMatrices 
              ((ds_.nextSolutionPtr),
               (ds_.nextStatePtr),
               (ds_.nextStateDerivPtr),
               (ds_.nextStorePtr),
               input,
               (ds_.dQdxVecVectorPtr),
               (ds_.dFdxVecVectorPtr));
  bsuccess = bsuccess && tmpBool;

  // Now determine the d(dQdt)/dx stuff:

  // now sum them all together, to create the total:
  // J = alpha/dt * dQ/dx + dF/dx
  wim_.applyJacobian(input, result);

  // Update the total load time
  jacobianTime_ = jacobianTimerPtr_->elapsedTime();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::loadSensitivityResiduals
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the sensitivity residual (RHS).
//
// Special Notes : All the contributions to the vector come (in the new DAE
//                 form) from the device package, as dQp, dFp, and dBdp.  
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_LOA_NonlinearEquationLoader::loadSensitivityResiduals ()
{
  bool bsuccess = true;
  wim_.obtainSensitivityResiduals ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_NonlinearEquationLoader::loadSensitivityResiduals
//
// Purpose       : This function manages the various function calls necessary
//                 to assemble the DAE form of the sensitivity residual (RHS).
//
// Special Notes : All the contributions to the vector come (in the new DAE
//                 form) from the device package, as dQp, dFp, and dBdp.  
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 09/02/2014
//-----------------------------------------------------------------------------
bool N_LOA_NonlinearEquationLoader::loadFinalSensitivityDerivatives () 
{
  bool bsuccess = true;
  wim_.loadFinalSensitivityDerivatives ();
  return bsuccess;
}

