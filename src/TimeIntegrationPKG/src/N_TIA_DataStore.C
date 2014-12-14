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
// Filename      : $RCSfile: N_TIA_DataStore.C,v $
//
// Purpose       : This file creates & initializes the data arrays needed for
//                 the time integration algorithms.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.152.2.9 $
//
// Revision Date  : $Date: 2014/09/03 11:14:50 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>

#include <N_TIA_DataStore.h>

#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_fwd.h>

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::N_TIA_DataStore
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_DataStore::N_TIA_DataStore(N_TIA_TIAParams * tiaPtr, N_LAS_System * lsPtr)
  : limiterFlag(false),
    lasSysPtr(lsPtr),
    solutionSize(tiaPtr->solutionSize),
    stateSize(tiaPtr->stateSize),
    tiaParamsPtr_(tiaPtr),
    nextSolPtrSwitched(false),
    tmpSolVectorPtr(0),
    tmpStaVectorPtr(0),
    tmpStaDerivPtr(0),
    tmpStaDivDiffPtr(0),
    tmpStoVectorPtr(0),
    xn0Ptr (0),
    currSolutionPtr(0),
    lastSolutionPtr(0),
    oldeSolutionPtr(0),
    nextSolutionPtr(0),
    flagSolutionPtr(0),
    savedNextSolutionPtr(0),
    currStatePtr(0),
    lastStatePtr(0),
    oldeStatePtr(0),
    nextStatePtr(0),
    currStorePtr(0),
    lastStorePtr(0),
    oldeStorePtr(0),
    nextStorePtr(0),
    currStoreLeadCurrQPtr(0),
    lastStoreLeadCurrQPtr(0),
    oldeStoreLeadCurrQPtr(0),
    nextStoreLeadCurrQPtr(0),
    currStoreLeadCurrQDerivPtr(0),
    lastStoreLeadCurrQDerivPtr(0),
    oldeStoreLeadCurrQDerivPtr(0),
    nextStoreLeadCurrQDerivPtr(0),
    currSolutionDerivPtr(0),
    lastSolutionDerivPtr(0),
    oldeSolutionDerivPtr(0),
    nextSolutionDerivPtr(0),
    currStateDerivPtr(0),
    lastStateDerivPtr(0),
    oldeStateDerivPtr(0),
    nextStateDerivPtr(0),
    resMatVecPtr(0),
    currSolutionDivDiffPtr(0),
    lastSolutionDivDiffPtr(0),
    oldeSolutionDivDiffPtr(0),
    nextSolutionDivDiffPtr(0),
    currStateDivDiffPtr(0),
    lastStateDivDiffPtr(0),
    oldeStateDivDiffPtr(0),
    nextStateDivDiffPtr(0),
    errWtVecPtr(0),
    absErrTolPtr(0),
    relErrTolPtr(0),
    JMatrixPtr(0),
    RHSVectorPtr(0),
#ifdef Xyce_DEBUG_DEVICE
    JdxpVectorPtr(0),
#endif
    newtonCorrectionPtr(0),
    deviceMaskPtr(0),
    indexVecsInitialized(false),

    qErrWtVecPtr(0),
    daeQVectorPtr(0),
    daeFVectorPtr(0),
    daeBVectorPtr(0),
    dQdxMatrixPtr(0),
    dFdxMatrixPtr(0),
    //saved_dQdxMatrixPtr(0),
    dQdxVecVectorPtr(0),
    dFdxVecVectorPtr(0),
    dFdxdVpVectorPtr(0),
    dQdxdVpVectorPtr(0),
    qn0Ptr(0),
    qpn0Ptr(0),
    sn0Ptr(0),
    spn0Ptr(0),
    ston0Ptr(0),
    stopn0Ptr(0),
    qNewtonCorrectionPtr(0),
    sNewtonCorrectionPtr(0),
    stoNewtonCorrectionPtr(0),
    stoLeadCurrQNewtonCorrectionPtr(0),
    delta_x(0),
    delta_q(0),
    tmpXn0APtr(0),
    tmpXn0BPtr(0)
{
  // temporary vectors:
  tmpSolVectorPtr = lasSysPtr->builder().createVector();
  tmpStaVectorPtr = lasSysPtr->builder().createStateVector();
  tmpStaDerivPtr = lasSysPtr->builder().createStateVector();
  tmpStaDivDiffPtr = lasSysPtr->builder().createStateVector();
  tmpStoVectorPtr = lasSysPtr->builder().createStoreVector();

  xn0Ptr = lasSysPtr->builder().createVector();

  // solution vectors:
  currSolutionPtr = lasSysPtr->builder().createVector();
  lastSolutionPtr = lasSysPtr->builder().createVector();
  nextSolutionPtr = lasSysPtr->builder().createVector();
  flagSolutionPtr = lasSysPtr->builder().createVector();

  // state vectors:
  currStatePtr    = lasSysPtr->builder().createStateVector();
  lastStatePtr    = lasSysPtr->builder().createStateVector();
  nextStatePtr    = lasSysPtr->builder().createStateVector();

  // store vectors:
  currStorePtr    = lasSysPtr->builder().createStoreVector();
  lastStorePtr    = lasSysPtr->builder().createStoreVector();
  nextStorePtr    = lasSysPtr->builder().createStoreVector();
  currStoreLeadCurrQPtr = lasSysPtr->builder().createStoreVector();
  lastStoreLeadCurrQPtr = lasSysPtr->builder().createStoreVector();
  nextStoreLeadCurrQPtr = lasSysPtr->builder().createStoreVector();

  // solution derivative vectors:
  currSolutionDerivPtr = lasSysPtr->builder().createVector();
  lastSolutionDerivPtr = lasSysPtr->builder().createVector();
  nextSolutionDerivPtr = lasSysPtr->builder().createVector();

  // state derivative vectors:
  currStateDerivPtr = lasSysPtr->builder().createStateVector();
  lastStateDerivPtr = lasSysPtr->builder().createStateVector();
  nextStateDerivPtr = lasSysPtr->builder().createStateVector();

  // store derivative vec for lead currents.
  currStoreLeadCurrQDerivPtr = lasSysPtr->builder().createStoreVector();
  lastStoreLeadCurrQDerivPtr = lasSysPtr->builder().createStoreVector();
  nextStoreLeadCurrQDerivPtr = lasSysPtr->builder().createStoreVector();

  // solution scaled divided differences:
  currSolutionDivDiffPtr = lasSysPtr->builder().createVector();
  lastSolutionDivDiffPtr = lasSysPtr->builder().createVector();
  nextSolutionDivDiffPtr = lasSysPtr->builder().createVector();

  // state scaled divided differences:
  currStateDivDiffPtr = lasSysPtr->builder().createStateVector();
  lastStateDivDiffPtr = lasSysPtr->builder().createStateVector();
  nextStateDivDiffPtr = lasSysPtr->builder().createStateVector();

  // error vectors:
  errWtVecPtr      = lasSysPtr->builder().createVector();
  absErrTolPtr     = lasSysPtr->builder().createVector();
  relErrTolPtr     = lasSysPtr->builder().createVector();

  errWtVecPtr->putScalar(1.0);

  deviceMaskPtr     = lasSysPtr->builder().createVector();
  deviceMaskPtr->putScalar(1.0);

  absErrTolPtr->putScalar(tiaParamsPtr_->absErrorTol);
  relErrTolPtr->putScalar(tiaParamsPtr_->relErrorTol);

  // nonlinear solution vectors:
  newtonCorrectionPtr = lasSysPtr->builder().createVector();

  // new-DAE stuff:
  // Error Vectors
  qErrWtVecPtr      = lasSysPtr->builder().createVector();

  // DAE formulation vectors
  daeQVectorPtr      = lasSysPtr->builder().createVector();
  daeFVectorPtr      = lasSysPtr->builder().createVector();
  daeBVectorPtr      = lasSysPtr->builder().createVector();

  // DAE formulation matrices
  dQdxMatrixPtr = lasSysPtr->builder().createMatrix();
  dFdxMatrixPtr = lasSysPtr->builder().createMatrix();
  //saved_dQdxMatrixPtr = lasSysPtr->builder().createMatrix();

  dQdxVecVectorPtr = lasSysPtr->builder().createVector();
  dFdxVecVectorPtr = lasSysPtr->builder().createVector();

  // History arrays
  int sizeOfHistory = tiaParamsPtr_->maxOrder+1;
  for (int i=0;i<sizeOfHistory;++i)
  {
    xHistory.push_back(lasSysPtr->builder().createVector());
    qHistory.push_back(lasSysPtr->builder().createVector());
    sHistory.push_back(lasSysPtr->builder().createStateVector());
    stoHistory.push_back(lasSysPtr->builder().createStoreVector());
    stoLeadCurrQHistory.push_back(lasSysPtr->builder().createStoreVector());
  }

  // Predictors
  qn0Ptr = lasSysPtr->builder().createVector();
  qpn0Ptr = lasSysPtr->builder().createVector();
  sn0Ptr = lasSysPtr->builder().createStateVector();
  spn0Ptr = lasSysPtr->builder().createStateVector();
  ston0Ptr = lasSysPtr->builder().createStoreVector();
  stopn0Ptr = lasSysPtr->builder().createStoreVector();
  stoQn0Ptr = lasSysPtr->builder().createStoreVector();
  stoQpn0Ptr = lasSysPtr->builder().createStoreVector();

  // Nonlinear solution vector:
  qNewtonCorrectionPtr = lasSysPtr->builder().createVector();
  sNewtonCorrectionPtr = lasSysPtr->builder().createStateVector();
  stoNewtonCorrectionPtr = lasSysPtr->builder().createStoreVector();
  stoLeadCurrQNewtonCorrectionPtr = lasSysPtr->builder().createStoreVector();

  // Step-size selection temporary vectors
  delta_x = lasSysPtr->builder().createVector();
  delta_q = lasSysPtr->builder().createVector();

  // Temporary vector for MPDE & WaMPDE interpolation
  tmpXn0APtr = lasSysPtr->builder().createVector();
  tmpXn0BPtr = lasSysPtr->builder().createVector();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::N_TIA_DataStore
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
N_TIA_DataStore::~N_TIA_DataStore()
{
  delete tmpSolVectorPtr;
  delete tmpStaVectorPtr;
  delete tmpStaDerivPtr;
  delete tmpStaDivDiffPtr;
  delete tmpStoVectorPtr;
  delete xn0Ptr;

  delete currSolutionPtr;
  delete lastSolutionPtr;

  delete nextSolutionPtr;
  delete flagSolutionPtr;

  delete currStatePtr;
  delete lastStatePtr;

  delete nextStatePtr;

  delete currStorePtr;
  delete lastStorePtr;
  delete nextStorePtr;
  delete currStoreLeadCurrQPtr;
  delete lastStoreLeadCurrQPtr;
  delete nextStoreLeadCurrQPtr;
  delete currStoreLeadCurrQDerivPtr;
  delete lastStoreLeadCurrQDerivPtr;
  delete nextStoreLeadCurrQDerivPtr;

  delete currSolutionDerivPtr;
  delete lastSolutionDerivPtr;

  delete nextSolutionDerivPtr;

  delete currStateDerivPtr;
  delete lastStateDerivPtr;

  delete nextStateDerivPtr;

  delete currSolutionDivDiffPtr;
  delete lastSolutionDivDiffPtr;

  delete nextSolutionDivDiffPtr;

  delete currStateDivDiffPtr;
  delete lastStateDivDiffPtr;

  delete nextStateDivDiffPtr;

  delete errWtVecPtr;
  delete absErrTolPtr;
  delete relErrTolPtr;

  delete deviceMaskPtr;

  delete newtonCorrectionPtr;

  //new-DAE:
  // Error Vectors
  delete qErrWtVecPtr;

  // DAE formulation vectors
  delete daeQVectorPtr;
  delete daeFVectorPtr;
  delete daeBVectorPtr;

  // DAE formulation matrices
  delete dQdxMatrixPtr;
  delete dFdxMatrixPtr;
  //delete saved_dQdxMatrixPtr;

  // HB temporary vectors
  delete dQdxVecVectorPtr;
  delete dFdxVecVectorPtr;

  Xyce::deleteList(xHistory.begin(), xHistory.end());
  Xyce::deleteList(qHistory.begin(), qHistory.end());
  Xyce::deleteList(sHistory.begin(), sHistory.end());
  Xyce::deleteList(stoHistory.begin(), stoHistory.end());
  Xyce::deleteList(stoLeadCurrQHistory.begin(), stoLeadCurrQHistory.end());

  delete qn0Ptr;
  delete qpn0Ptr;
  delete sn0Ptr;
  delete spn0Ptr;
  delete ston0Ptr;
  delete stopn0Ptr;
  delete stoQn0Ptr;
  delete stoQpn0Ptr;

  // Nonlinear solution vector:
  delete qNewtonCorrectionPtr;
  delete sNewtonCorrectionPtr;
  delete stoNewtonCorrectionPtr;
  delete stoLeadCurrQNewtonCorrectionPtr;

  // Step-size selection temporary vectors
  delete delta_x;
  delete delta_q;

  // Temporary vector for WaMPDE interpolation
  delete tmpXn0APtr;
  delete tmpXn0BPtr;

  // Delete data in the fast time storage for HB and MPDE
  Xyce::deleteList(fastTimeSolutionVec.begin(), fastTimeSolutionVec.end());
  Xyce::deleteList(fastTimeStateVec.begin(), fastTimeStateVec.end());
  Xyce::deleteList(fastTimeQVec.begin(), fastTimeQVec.end());
  Xyce::deleteList(fastTimeStoreVec.begin(), fastTimeStoreVec.end());

  // delete the sensitivity related stuff, if necessary
  deleteSensitivityArrays();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::deleteSensitivityArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/6/2014
//-----------------------------------------------------------------------------
void N_TIA_DataStore::deleteSensitivityArrays()
{
  int numParams = currDfdpPtrVector.size();

  if (resMatVecPtr != 0) delete resMatVecPtr;

  for (int ip=0;ip<numParams;++ip)
  {
    delete sensRHSPtrVector[ip];

    delete currDfdpPtrVector[ip];
    delete lastDfdpPtrVector[ip];
    delete nextDfdpPtrVector[ip];

    delete currDqdpPtrVector[ip];
    delete lastDqdpPtrVector[ip];
    delete nextDqdpPtrVector[ip];

    delete currDbdpPtrVector[ip];
    delete lastDbdpPtrVector[ip];
    delete nextDbdpPtrVector[ip];

    delete currDXdpPtrVector[ip];
    delete lastDXdpPtrVector[ip];
    delete nextDXdpPtrVector[ip];

    delete currDQdxDXdpPtrVector[ip];
    delete lastDQdxDXdpPtrVector[ip];
    delete nextDQdxDXdpPtrVector[ip];

    delete currDQdxDXdpDerivPtrVector[ip];
    delete lastDQdxDXdpDerivPtrVector[ip];
    delete nextDQdxDXdpDerivPtrVector[ip];

    delete currDqdpDerivPtrVector[ip];
    delete lastDqdpDerivPtrVector[ip];
    delete nextDqdpDerivPtrVector[ip];

    // predictors:
    delete dfdpn0PtrVector[ip];
    delete dqdpn0PtrVector[ip];
    delete dbdpn0PtrVector[ip];
    delete dXdpn0PtrVector[ip];

    delete dfdppn0PtrVector[ip];
    delete dqdppn0PtrVector[ip];
    delete dbdppn0PtrVector[ip];
    delete dXdppn0PtrVector[ip];

    // newton correction:
    delete dfdpNewtonCorrectionPtrVector[ip];
    delete dqdpNewtonCorrectionPtrVector[ip];
    delete dbdpNewtonCorrectionPtrVector[ip];
    delete dXdpNewtonCorrectionPtrVector[ip];

    // history
    int sizeOfHistory = tiaParamsPtr_->maxOrder+1;
    for (int i=0;i<sizeOfHistory;++i)
    {
      delete dfdpHistory[ip][i];
      delete dqdpHistory[ip][i];
      delete dbdpHistory[ip][i];
      delete dXdpHistory[ip][i];
      delete dQdxdXdpHistory[ip][i];
    }
    dfdpHistory[ip].clear();
    dqdpHistory[ip].clear();
    dbdpHistory[ip].clear();
    dXdpHistory[ip].clear();
    dQdxdXdpHistory[ip].clear();
  }

  sensRHSPtrVector.clear();

  currDfdpPtrVector.clear();
  lastDfdpPtrVector.clear();
  nextDfdpPtrVector.clear();

  currDqdpPtrVector.clear();
  lastDqdpPtrVector.clear();
  nextDqdpPtrVector.clear();

  currDbdpPtrVector.clear();
  lastDbdpPtrVector.clear();
  nextDbdpPtrVector.clear();

  currDXdpPtrVector.clear();
  lastDXdpPtrVector.clear();
  nextDXdpPtrVector.clear();

  currDQdxDXdpPtrVector.clear();
  lastDQdxDXdpPtrVector.clear();
  nextDQdxDXdpPtrVector.clear();

  currDQdxDXdpDerivPtrVector.clear();
  lastDQdxDXdpDerivPtrVector.clear();
  nextDQdxDXdpDerivPtrVector.clear();

  currDqdpDerivPtrVector.clear();
  lastDqdpDerivPtrVector.clear();
  nextDqdpDerivPtrVector.clear();

  
  // predictors:
  dfdpn0PtrVector.clear();
  dqdpn0PtrVector.clear();
  dbdpn0PtrVector.clear();
  dXdpn0PtrVector.clear();

  dfdppn0PtrVector.clear();
  dqdppn0PtrVector.clear();
  dbdppn0PtrVector.clear();
  dXdppn0PtrVector.clear();

  // newton correction:
  dfdpNewtonCorrectionPtrVector.clear();
  dqdpNewtonCorrectionPtrVector.clear();
  dbdpNewtonCorrectionPtrVector.clear();
  dXdpNewtonCorrectionPtrVector.clear();

  // history
  dfdpHistory.clear();
  dqdpHistory.clear();
  dbdpHistory.clear();
  dXdpHistory.clear();
  dQdxdXdpHistory.clear();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::allocateSensitivityArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/6/2014
//-----------------------------------------------------------------------------
void N_TIA_DataStore::allocateSensitivityArrays(int numParams)
{
  N_LAS_Builder & bld_ = lasSysPtr->builder();

  resMatVecPtr = bld_.createVector();

  sensRHSPtrVector.resize(numParams);

  currDfdpPtrVector.resize(numParams);
  lastDfdpPtrVector.resize(numParams);
  oldeDfdpPtrVector.resize(numParams);
  nextDfdpPtrVector.resize(numParams);

  currDqdpPtrVector.resize(numParams);
  lastDqdpPtrVector.resize(numParams);
  oldeDqdpPtrVector.resize(numParams);
  nextDqdpPtrVector.resize(numParams);

  currDbdpPtrVector.resize(numParams);
  lastDbdpPtrVector.resize(numParams);
  oldeDbdpPtrVector.resize(numParams);
  nextDbdpPtrVector.resize(numParams);

  currDXdpPtrVector.resize(numParams);
  lastDXdpPtrVector.resize(numParams);
  oldeDXdpPtrVector.resize(numParams);
  nextDXdpPtrVector.resize(numParams);

  currDQdxDXdpPtrVector.resize(numParams);
  lastDQdxDXdpPtrVector.resize(numParams);
  oldeDQdxDXdpPtrVector.resize(numParams);
  nextDQdxDXdpPtrVector.resize(numParams);

  currDQdxDXdpDerivPtrVector.resize(numParams);
  lastDQdxDXdpDerivPtrVector.resize(numParams);
  oldeDQdxDXdpDerivPtrVector.resize(numParams);
  nextDQdxDXdpDerivPtrVector.resize(numParams);

  currDqdpDerivPtrVector.resize(numParams);
  lastDqdpDerivPtrVector.resize(numParams);
  oldeDqdpDerivPtrVector.resize(numParams);
  nextDqdpDerivPtrVector.resize(numParams);

  
  // predictors:
  dfdpn0PtrVector.resize(numParams);
  dqdpn0PtrVector.resize(numParams);
  dbdpn0PtrVector.resize(numParams);
  dXdpn0PtrVector.resize(numParams);

  dfdppn0PtrVector.resize(numParams);
  dqdppn0PtrVector.resize(numParams);
  dbdppn0PtrVector.resize(numParams);
  dXdppn0PtrVector.resize(numParams);

  // newton correction:
  dfdpNewtonCorrectionPtrVector.resize(numParams);
  dqdpNewtonCorrectionPtrVector.resize(numParams);
  dbdpNewtonCorrectionPtrVector.resize(numParams);
  dXdpNewtonCorrectionPtrVector.resize(numParams);

  // history
  dfdpHistory.resize(numParams);
  dqdpHistory.resize(numParams);
  dbdpHistory.resize(numParams);
  dXdpHistory.resize(numParams);
  dQdxdXdpHistory.resize(numParams);

  for (int ip=0;ip<numParams;++ip)
  {
    sensRHSPtrVector[ip] = bld_.createVector();

    currDfdpPtrVector[ip] = bld_.createVector();
    lastDfdpPtrVector[ip] = bld_.createVector();
    nextDfdpPtrVector[ip] = bld_.createVector();

    currDqdpPtrVector[ip] = bld_.createVector();
    lastDqdpPtrVector[ip] = bld_.createVector();
    nextDqdpPtrVector[ip] = bld_.createVector();

    currDbdpPtrVector[ip] = bld_.createVector();
    lastDbdpPtrVector[ip] = bld_.createVector();
    nextDbdpPtrVector[ip] = bld_.createVector();

    currDXdpPtrVector[ip] = bld_.createVector();
    lastDXdpPtrVector[ip] = bld_.createVector();
    nextDXdpPtrVector[ip] = bld_.createVector();

    currDQdxDXdpPtrVector[ip] = bld_.createVector();
    lastDQdxDXdpPtrVector[ip] = bld_.createVector();
    nextDQdxDXdpPtrVector[ip] = bld_.createVector();

    currDQdxDXdpDerivPtrVector[ip] = bld_.createVector();
    lastDQdxDXdpDerivPtrVector[ip] = bld_.createVector();
    nextDQdxDXdpDerivPtrVector[ip] = bld_.createVector();

    currDqdpDerivPtrVector[ip] = bld_.createVector();
    lastDqdpDerivPtrVector[ip] = bld_.createVector();
    nextDqdpDerivPtrVector[ip] = bld_.createVector();

    // predictors:
    dfdpn0PtrVector[ip] = bld_.createVector();
    dqdpn0PtrVector[ip] = bld_.createVector();
    dbdpn0PtrVector[ip] = bld_.createVector();
    dXdpn0PtrVector[ip] = bld_.createVector();

    dfdppn0PtrVector[ip] = bld_.createVector();
    dqdppn0PtrVector[ip] = bld_.createVector();
    dbdppn0PtrVector[ip] = bld_.createVector();
    dXdppn0PtrVector[ip] = bld_.createVector();

    // newton correction:
    dfdpNewtonCorrectionPtrVector[ip] = bld_.createVector();
    dqdpNewtonCorrectionPtrVector[ip] = bld_.createVector();
    dbdpNewtonCorrectionPtrVector[ip] = bld_.createVector();
    dXdpNewtonCorrectionPtrVector[ip] = bld_.createVector();

    // history
    int sizeOfHistory = tiaParamsPtr_->maxOrder+1;
    for (int i=0;i<sizeOfHistory;++i)
    {
      dfdpHistory[ip].push_back(bld_.createVector());
      dqdpHistory[ip].push_back(bld_.createVector());
      dbdpHistory[ip].push_back(bld_.createVector());
      dXdpHistory[ip].push_back(bld_.createVector());
      dQdxdXdpHistory[ip].push_back(bld_.createVector());
    }
  }
}

// ----------------------------------------------------------------------------
// -----------------------  DataStore Class Functions -------------------------
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::printOutPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/23/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::printOutPointers()
{
  Xyce::dout() << "olde ptr = " << oldeSolutionPtr << std::endl;
  Xyce::dout() << "last ptr = " << lastSolutionPtr << std::endl;
  Xyce::dout() << "curr ptr = " << currSolutionPtr << std::endl;
  Xyce::dout() << "next ptr = " << nextSolutionPtr << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setConstantHistory
// Purpose       : This function is called  after the operating point
//                 calculation has been called.  Once the operating point
//                 solution has been obtained, the code should regard that
//                 solution as having been the existing constant solution since
//                 the dawn of time.
// Special Notes : The most recent solution, etc., are in the "next" vectors.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/23/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::setConstantHistory()
{
#ifdef Xyce_DEBUG_TIME
  Xyce::dout() << "\nN_TIA_DataStore::setConstantHistory" << std::endl;
#endif

  // Solutions:
  *(lastSolutionPtr) = *(nextSolutionPtr);
  *(currSolutionPtr) = *(nextSolutionPtr);

  // Derivative of Solutions:
  *(lastSolutionDerivPtr) = *(nextSolutionDerivPtr);
  *(currSolutionDerivPtr) = *(nextSolutionDerivPtr);

  // Scaled Divided Differences of solutions:
  *(lastSolutionDivDiffPtr) = *(nextSolutionDivDiffPtr);
  *(currSolutionDivDiffPtr) = *(nextSolutionDivDiffPtr);

  // States:
  *(lastStatePtr) = *(nextStatePtr);
  *(currStatePtr) = *(nextStatePtr);

  // Stores:
  *(lastStorePtr) = *(nextStorePtr);
  *(currStorePtr) = *(nextStorePtr);

  *(lastStoreLeadCurrQPtr) = *(nextStoreLeadCurrQPtr);
  *(currStoreLeadCurrQPtr) = *(nextStoreLeadCurrQPtr);

  // Derivative of States:
  *(lastStateDerivPtr) = *(nextStateDerivPtr);
  *(currStateDerivPtr) = *(nextStateDerivPtr);

  // lead current derivative info in store
  *(lastStoreLeadCurrQDerivPtr) = *(nextStoreLeadCurrQDerivPtr);
  *(currStoreLeadCurrQDerivPtr) = *(nextStoreLeadCurrQDerivPtr);

  // Scaled Divided Differences of states:
  *(lastStateDivDiffPtr) = *(nextStateDivDiffPtr);
  *(currStateDivDiffPtr) = *(nextStateDivDiffPtr);

  setConstantSensitivityHistory();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setConstantSensitivityHistory
// Purpose       : 
// Special Notes : The most recent solution, etc., are in the "next" vectors.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void N_TIA_DataStore::setConstantSensitivityHistory()
{
  // Sensitivities:
  for (int ip=0;ip<nextDfdpPtrVector.size();++ip)
  {
    *(lastDfdpPtrVector[ip]) = *(nextDfdpPtrVector[ip]);
    *(currDfdpPtrVector[ip]) = *(nextDfdpPtrVector[ip]);

    *(lastDqdpPtrVector[ip]) = *(nextDqdpPtrVector[ip]);
    *(currDqdpPtrVector[ip]) = *(nextDqdpPtrVector[ip]);

    *(lastDbdpPtrVector[ip]) = *(nextDbdpPtrVector[ip]);
    *(currDbdpPtrVector[ip]) = *(nextDbdpPtrVector[ip]);

    *(lastDqdpDerivPtrVector[ip]) = *(nextDqdpDerivPtrVector[ip]);
    *(currDqdpDerivPtrVector[ip]) = *(nextDqdpDerivPtrVector[ip]);

    *(lastDXdpPtrVector[ip]) = *(nextDXdpPtrVector[ip]);
    *(currDXdpPtrVector[ip]) = *(nextDXdpPtrVector[ip]);

    *(lastDQdxDXdpPtrVector[ip]) = *(nextDQdxDXdpPtrVector[ip]);
    *(currDQdxDXdpPtrVector[ip]) = *(nextDQdxDXdpPtrVector[ip]);

    *(lastDQdxDXdpDerivPtrVector[ip]) = *(nextDQdxDXdpDerivPtrVector[ip]);
    *(currDQdxDXdpDerivPtrVector[ip]) = *(nextDQdxDXdpDerivPtrVector[ip]);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::resetAll
//
// Purpose       : This function resets everything so that a transient loop
//                 can be started from the beginning.
//
// Special Notes : This function was needed for HB.
//
// Scope         : public
// Creator       : T. Mei, SNL
// Creation Date : 02/26/09
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::resetAll ()
{
  absErrTolPtr->putScalar(tiaParamsPtr_->absErrorTol);
  relErrTolPtr->putScalar(tiaParamsPtr_->relErrorTol);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::resetFastTimeData
//
// Purpose       : This function deletes the information from all the vectors
//                 that store fast time data for HB and MPDE
//
// Special Notes : This function was needed for HB.
//
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 06/05/13
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::resetFastTimeData()
{
  // Clear the time step vectors
  timeSteps.clear();
  timeStepsBreakpointFlag.clear();

  // Delete any stored up any solution or state info
  Xyce::deleteList(fastTimeSolutionVec.begin(), fastTimeSolutionVec.end());
  Xyce::deleteList(fastTimeStateVec.begin(), fastTimeStateVec.end());
  Xyce::deleteList(fastTimeQVec.begin(), fastTimeQVec.end());
  Xyce::deleteList(fastTimeStoreVec.begin(), fastTimeStoreVec.end());

  fastTimeSolutionVec.clear();
  fastTimeStateVec.clear();
  fastTimeQVec.clear();
  fastTimeStoreVec.clear();

  return true;
}



//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::updateSolDataArrays
// Purpose       : Update the necessary integration data arrays for
//                 preparation of the next integration step. This is done after
//                 a successful step has been taken.
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::updateSolDataArrays()
{
#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    Xyce::dout() << "\nN_TIA_DataStore::updateSolDataArrays " << std::endl;
  }
#endif

  // if the next solution has been switched out (probably because of NOX),
  // then reset it
  if (nextSolPtrSwitched) unsetNextSolVectorPtr();

  // Solutions:
  oldeSolutionPtr = lastSolutionPtr;
  lastSolutionPtr = currSolutionPtr;
  currSolutionPtr = nextSolutionPtr;
  nextSolutionPtr = oldeSolutionPtr;

  // Derivative of Solutions:
  oldeSolutionDerivPtr = lastSolutionDerivPtr;
  lastSolutionDerivPtr = currSolutionDerivPtr;
  currSolutionDerivPtr = nextSolutionDerivPtr;
  nextSolutionDerivPtr = oldeSolutionDerivPtr;

  // Scaled Divided Differences of solutions:
  oldeSolutionDivDiffPtr = lastSolutionDivDiffPtr;
  lastSolutionDivDiffPtr = currSolutionDivDiffPtr;
  currSolutionDivDiffPtr = nextSolutionDivDiffPtr;
  nextSolutionDivDiffPtr = oldeSolutionDivDiffPtr;

  // States:
  oldeStatePtr = lastStatePtr;
  lastStatePtr = currStatePtr;
  currStatePtr = nextStatePtr;
  nextStatePtr = oldeStatePtr;

  // Stores:
  oldeStorePtr = lastStorePtr;
  lastStorePtr = currStorePtr;
  currStorePtr = nextStorePtr;
  nextStorePtr = oldeStorePtr;
  oldeStoreLeadCurrQPtr = lastStoreLeadCurrQPtr;
  lastStoreLeadCurrQPtr = currStoreLeadCurrQPtr;
  currStoreLeadCurrQPtr = nextStoreLeadCurrQPtr;
  nextStoreLeadCurrQPtr = oldeStoreLeadCurrQPtr;

  // Derivative of States:
  oldeStateDerivPtr = lastStateDerivPtr;
  lastStateDerivPtr = currStateDerivPtr;
  currStateDerivPtr = nextStateDerivPtr;
  nextStateDerivPtr = oldeStateDerivPtr;

  // lead curent component of store
  oldeStoreLeadCurrQDerivPtr = lastStoreLeadCurrQDerivPtr;
  lastStoreLeadCurrQDerivPtr = currStoreLeadCurrQDerivPtr;
  currStoreLeadCurrQDerivPtr = nextStoreLeadCurrQDerivPtr;
  nextStoreLeadCurrQDerivPtr = oldeStoreLeadCurrQDerivPtr;

  // Scaled Divided Differences of states:
  oldeStateDivDiffPtr = lastStateDivDiffPtr;
  lastStateDivDiffPtr = currStateDivDiffPtr;
  currStateDivDiffPtr = nextStateDivDiffPtr;
  nextStateDivDiffPtr = oldeStateDivDiffPtr;

  // Sensitivities, if they were requested:
  int iparamSizeF = currDfdpPtrVector.size();
  for (int ip=0;ip<iparamSizeF;++ip)
  {
    oldeDfdpPtrVector[ip] = lastDfdpPtrVector[ip];
    lastDfdpPtrVector[ip] = currDfdpPtrVector[ip];
    currDfdpPtrVector[ip] = nextDfdpPtrVector[ip];
    nextDfdpPtrVector[ip] = oldeDfdpPtrVector[ip];
  }

  int iparamSizeQ = currDqdpPtrVector.size();
  for (int ip=0;ip<iparamSizeQ;++ip)
  {
    oldeDqdpPtrVector[ip] = lastDqdpPtrVector[ip];
    lastDqdpPtrVector[ip] = currDqdpPtrVector[ip];
    currDqdpPtrVector[ip] = nextDqdpPtrVector[ip];
    nextDqdpPtrVector[ip] = oldeDqdpPtrVector[ip];
  }

  int iparamSizeQderiv = currDqdpDerivPtrVector.size();
  for (int ip=0;ip<iparamSizeQderiv;++ip)
  {
    oldeDqdpDerivPtrVector[ip] = lastDqdpDerivPtrVector[ip];
    lastDqdpDerivPtrVector[ip] = currDqdpDerivPtrVector[ip];
    currDqdpDerivPtrVector[ip] = nextDqdpDerivPtrVector[ip];
    nextDqdpDerivPtrVector[ip] = oldeDqdpDerivPtrVector[ip];
  }

  int iparamSizeB = currDbdpPtrVector.size();
  for (int ip=0;ip<iparamSizeB;++ip)
  {
    oldeDbdpPtrVector[ip] = lastDbdpPtrVector[ip];
    lastDbdpPtrVector[ip] = currDbdpPtrVector[ip];
    currDbdpPtrVector[ip] = nextDbdpPtrVector[ip];
    nextDbdpPtrVector[ip] = oldeDbdpPtrVector[ip];
  }

  int iparamSizeX = currDXdpPtrVector.size();
  for (int ip=0;ip<iparamSizeX;++ip)
  {
    oldeDXdpPtrVector[ip] = lastDXdpPtrVector[ip];
    lastDXdpPtrVector[ip] = currDXdpPtrVector[ip];
    currDXdpPtrVector[ip] = nextDXdpPtrVector[ip];
    nextDXdpPtrVector[ip] = oldeDXdpPtrVector[ip];
  }

  int iparamSizeQX = currDQdxDXdpPtrVector.size();
  for (int ip=0;ip<iparamSizeQX;++ip)
  {
    oldeDQdxDXdpPtrVector[ip] = lastDQdxDXdpPtrVector[ip];
    lastDQdxDXdpPtrVector[ip] = currDQdxDXdpPtrVector[ip];
    currDQdxDXdpPtrVector[ip] = nextDQdxDXdpPtrVector[ip];
    nextDQdxDXdpPtrVector[ip] = oldeDQdxDXdpPtrVector[ip];
  }


  int iparamSizeQXd = currDQdxDXdpDerivPtrVector.size();
  for (int ip=0;ip<iparamSizeQXd;++ip)
  {
    oldeDQdxDXdpDerivPtrVector[ip] = lastDQdxDXdpDerivPtrVector[ip];
    lastDQdxDXdpDerivPtrVector[ip] = currDQdxDXdpDerivPtrVector[ip];
    currDQdxDXdpDerivPtrVector[ip] = nextDQdxDXdpDerivPtrVector[ip];
    nextDQdxDXdpDerivPtrVector[ip] = oldeDQdxDXdpDerivPtrVector[ip];
  }

  // copy contents of "curr" into "next".  This is to insure
  // that at a minimum, the initial guess for the Newton solve
  // will at least be the results of the previous Newton solve.
  *(nextSolutionPtr) = *(currSolutionPtr);
  *(nextStatePtr)    = *(currStatePtr);
  *(nextStorePtr)    = *(currStorePtr);

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::updateStateDataArrays
//
//
// Purpose       : Same as updateSolDataArrays, but this function only
//                 advances the state vector, and leaves the
//                 solution alone.
//
// Special Notes : The main usage of this function is LOCA.
//                 After each continuation step, LOCA needs to call
//                 this function.  LOCA keeps track of solution vectors
//                 on its own, which is why updateSolDataArrays
//                 is inappropriate for LOCA.
//
//                 This is necessary for voltage limiting to be
//                 consistent with LOCA.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, 9233.
// Creation Date : 3/06/05
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::updateStateDataArrays()
{
#ifdef Xyce_DEBUG_TIME
  Xyce::dout() << "\nN_TIA_DataStore::updateStateDataArrays " << std::endl;
#endif
  // States:
  oldeStatePtr = lastStatePtr;
  lastStatePtr = currStatePtr;
  currStatePtr = nextStatePtr;
  nextStatePtr = oldeStatePtr;

  // Stores:
  oldeStorePtr = lastStorePtr;
  lastStorePtr = currStorePtr;
  currStorePtr = nextStorePtr;
  nextStorePtr = oldeStorePtr;

  // q component of lead current
  oldeStoreLeadCurrQPtr = lastStoreLeadCurrQPtr;
  lastStoreLeadCurrQPtr = currStoreLeadCurrQPtr;
  currStoreLeadCurrQPtr = nextStoreLeadCurrQPtr;
  nextStoreLeadCurrQPtr = oldeStoreLeadCurrQPtr;

  // Derivative of States:
  oldeStateDerivPtr = lastStateDerivPtr;
  lastStateDerivPtr = currStateDerivPtr;
  currStateDerivPtr = nextStateDerivPtr;
  nextStateDerivPtr = oldeStateDerivPtr;

  // Scaled Divided Differences of states:
  oldeStateDivDiffPtr = lastStateDivDiffPtr;
  lastStateDivDiffPtr = currStateDivDiffPtr;
  currStateDivDiffPtr = nextStateDivDiffPtr;
  nextStateDivDiffPtr = oldeStateDivDiffPtr;


  // Now, make the "next" stuff the same as the "curr" stuff.
  // This is done because at the end of the tranop, but before
  // the transient phase starts, the function setConstantHistory
  // will be called.  When it is called, the most recent values
  // for state variables need to be in the "next" vectors.
  //
  // As long as LOCA solves are never used for transient, and only
  // for DC and tranop solves, this is OK.  This is, of course,
  // a bit of a kludge.
  *(nextStatePtr) = *(currStatePtr);
  *(nextStorePtr) = *(currStorePtr);

  return true;
}

#ifndef OLD_PLOT
//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputSolDataArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::outputSolDataArrays(std::ostream &os)
{
  char tmp[256];
  os << std::endl
     << Xyce::section_divider << std::endl
     << std::endl
     << "  Solution Vectors:\n          Current               Last               Olde               Error" << std::endl;

  for (unsigned int k = 0; k < solutionSize; ++k)
  {
    os << (*(currSolutionPtr))[k] << (*(lastSolutionPtr))[k] << (*(oldeSolutionPtr))[k] << (*(newtonCorrectionPtr))[k] << std::endl;
  }

  os << Xyce::section_divider << std::endl;
}

#else

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputSolDataArrays
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::outputSolDataArrays(std::ostream &os)
{
  for (unsigned int k = 0; k < solutionSize; ++k)
  {
    os << "\t" <<  (*(currSolutionPtr))[0][k] << std::endl;
  }
  os << std::endl;
}

#endif

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::enableOrderOneStart
// Purpose       : Initialize arrays for "re-starting" integration at order 1.
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::enableOrderOneStart()
{
  // solution vectors:
  *lastSolutionPtr = *currSolutionPtr;
  *oldeSolutionPtr = *currSolutionPtr;

  nextSolutionDerivPtr->putScalar(0.0);
  currSolutionDerivPtr->putScalar(0.0);
  currSolutionDivDiffPtr->putScalar(0.0);

  // state vectors:
  *lastStatePtr = *currStatePtr;
  *oldeStatePtr = *currStatePtr;

  nextStateDerivPtr->putScalar(0.0);
  currStateDerivPtr->putScalar(0.0);
  currStateDivDiffPtr->putScalar(0.0);

  // store vectors:
  *lastStorePtr = *currStorePtr;
  *oldeStorePtr = *currStorePtr;
  *lastStoreLeadCurrQPtr = *currStoreLeadCurrQPtr;
  *oldeStoreLeadCurrQPtr = *currStoreLeadCurrQPtr;
  nextStoreLeadCurrQDerivPtr->putScalar(0.0);
  currStoreLeadCurrQDerivPtr->putScalar(0.0);

  // sensitivities
  int iparamSizeF = currDfdpPtrVector.size();
  for (int ip=0;ip<iparamSizeF;++ip)
  {
    *(oldeDfdpPtrVector[ip]) = *(currDfdpPtrVector[ip]);
    *(lastDfdpPtrVector[ip]) = *(currDfdpPtrVector[ip]);
    currDfdpPtrVector[ip]->putScalar(0.0);
    nextDfdpPtrVector[ip]->putScalar(0.0);
  }

  int iparamSizeQ = currDqdpPtrVector.size();
  for (int ip=0;ip<iparamSizeQ;++ip)
  {
    *(oldeDqdpPtrVector[ip]) = *(currDqdpPtrVector[ip]);
    *(lastDqdpPtrVector[ip]) = *(currDqdpPtrVector[ip]);
    currDqdpPtrVector[ip]->putScalar(0.0);
    nextDqdpPtrVector[ip]->putScalar(0.0);
  }

  int iparamSizeQd = currDqdpDerivPtrVector.size();
  for (int ip=0;ip<iparamSizeQ;++ip)
  {
    *(oldeDqdpDerivPtrVector[ip]) = *(currDqdpDerivPtrVector[ip]);
    *(lastDqdpDerivPtrVector[ip]) = *(currDqdpDerivPtrVector[ip]);
    currDqdpDerivPtrVector[ip]->putScalar(0.0);
    nextDqdpDerivPtrVector[ip]->putScalar(0.0);
  }

  int iparamSizeB = currDbdpPtrVector.size();
  for (int ip=0;ip<iparamSizeB;++ip)
  {
    *(oldeDbdpPtrVector[ip]) = *(currDbdpPtrVector[ip]);
    *(lastDbdpPtrVector[ip]) = *(currDbdpPtrVector[ip]);
    currDbdpPtrVector[ip]->putScalar(0.0);
    nextDbdpPtrVector[ip]->putScalar(0.0);
  }

  int iparamSizeX = currDXdpPtrVector.size();
  for (int ip=0;ip<iparamSizeX;++ip)
  {
    *(oldeDXdpPtrVector[ip]) = *(currDXdpPtrVector[ip]);
    *(lastDXdpPtrVector[ip]) = *(currDXdpPtrVector[ip]);
    currDXdpPtrVector[ip]->putScalar(0.0);
    nextDXdpPtrVector[ip]->putScalar(0.0);
  }

  int iparamSizeQX = currDQdxDXdpPtrVector.size();
  for (int ip=0;ip<iparamSizeQX;++ip)
  {
    *(oldeDQdxDXdpPtrVector[ip]) = *(currDQdxDXdpPtrVector[ip]);
    *(lastDQdxDXdpPtrVector[ip]) = *(currDQdxDXdpPtrVector[ip]);
    currDQdxDXdpPtrVector[ip]->putScalar(0.0);
    nextDQdxDXdpPtrVector[ip]->putScalar(0.0);
  }


  int iparamSizeQXd = currDQdxDXdpDerivPtrVector.size();
  for (int ip=0;ip<iparamSizeQX;++ip)
  {
    *(oldeDQdxDXdpDerivPtrVector[ip]) = *(currDQdxDXdpDerivPtrVector[ip]);
    *(lastDQdxDXdpDerivPtrVector[ip]) = *(currDQdxDXdpDerivPtrVector[ip]);
    currDQdxDXdpDerivPtrVector[ip]->putScalar(0.0);
    nextDQdxDXdpDerivPtrVector[ip]->putScalar(0.0);
  }

}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputPredictedSolution
// Purpose       : Set ErrorEstimate array values to zero.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::outputPredictedSolution(std::ostream &os)
{
  os << Xyce::subsection_divider << std::endl;

  os << "  Predicted Solution:" << std::endl;

  for (unsigned int k = 0; k< solutionSize; ++k)
  {
    os << (*(nextSolutionPtr))[k] << std::endl;
  }

  os << Xyce::subsection_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::outputPredictedDerivative
// Purpose       : Set ErrorEstimate array values to zero.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::outputPredictedDerivative(std::ostream &os)
{
  os << Xyce::subsection_divider << std::endl;

  os << "  Predicted Derivative:" << std::endl;

  for (unsigned int k = 0; k < solutionSize; ++k)
  {
    os << (*(nextSolutionDerivPtr))[k] << std::endl;
  }

  os << Xyce::subsection_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialErrorNormSum
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialErrorNormSum ()
{
  double errorNorm = 0.0;
  newtonCorrectionPtr->wRMSNorm(*errWtVecPtr, &errorNorm);
  double sum = errorNorm*errorNorm;
  double length = newtonCorrectionPtr->globalLength();

  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::globalLength
// Purpose       : Needed by 2-level solves.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double N_TIA_DataStore::globalLength ()
{
  return newtonCorrectionPtr->globalLength();
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::computeDividedDifferences
// Purpose       : Compute the scaled (by current stepsize) divided difference
//                 approximation to the derivative.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
void N_TIA_DataStore::computeDividedDifferences()
{
  // First do the solution divided difference:
  nextSolutionDivDiffPtr->linearCombo(1.0, *(nextSolutionPtr), -1.0, *(currSolutionPtr));

  // Now do the state divided difference:
  nextStateDivDiffPtr->linearCombo(1.0, *(nextStatePtr), -1.0, *(currStatePtr));

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::computeDivDiffsBlock
// Purpose       : Compute the scaled (by current stepsize) divided difference
//                 approximation to the derivative.  This is the same as
//                 the  function "computeDividedDifferences", except that
//                 the  operation is only performed on a sub-block of the
//                 vectors.
//
// Special Notes : Not done yet...
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
void N_TIA_DataStore::computeDivDiffsBlock (
         const std::list<index_pair> & solGIDList,
         const std::list<index_pair> & staGIDList)
{}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::equateTmpVectors
// Purpose       : This function equates the 6 temporary vectors  with
//                 their "next" vector equivalents.  This function
//                 is neccessary for the nonlinear solver damping loop.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::equateTmpVectors()
{
  // next solution vector:
  *(tmpSolVectorPtr) = *(nextSolutionPtr);

  // next state vector:
  *(tmpStaVectorPtr) = *(nextStatePtr);

  // next store vector:
  *(tmpStoVectorPtr) = *(nextStorePtr);

  // next state divided difference vector:
  *(tmpStaDivDiffPtr) = *(nextStateDivDiffPtr);

  // next state derivative vector:
  *(tmpStaDerivPtr)  = *(nextStateDerivPtr);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::usePreviousSolAsPredictor
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::usePreviousSolAsPredictor ()
{
  bool bsuccess = true;

  *(nextSolutionPtr) = *(currSolutionPtr);
  *(nextStatePtr)    = *(currStatePtr);
  *(nextStorePtr)    = *(currStorePtr);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setNextSolVectorPtr
// Purpose       :
// Special Notes : Only needed for NOX, and other solvers that prefer to
//                 own the solution vector.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/03
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setNextSolVectorPtr (N_LAS_Vector * solVecPtr)
{
  // only save the old pointer if it hasn't been switched yet.
  if (!nextSolPtrSwitched)
  {
    savedNextSolutionPtr = nextSolutionPtr;
    nextSolPtrSwitched = true;
  }
  nextSolutionPtr = solVecPtr;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::unsetNextSolVectorPtr
// Purpose       :
// Special Notes : This also copies over the solution, in addition to the
//                 pointer.
//
//                 This is only called when it is time to rotate the
//                 pointers for the next time step.  If we have been
//                 running with NOX, or with some other solver that prefers
//                 to own the solution vector, this is neccessary, and the
//                 switch flag should be set.  Otherwise, not.
//
//                 Basically, if we've been running with NOX, then the next
//                 solution vector ptr has probably been switched out at least
//                 once.  We need to maintain the history, so we make a
//                 copy of this switched solution vector, and the
//                 restore the old pointer.
//
//                 This is kludgy, but will have to do for now.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/22/03
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::unsetNextSolVectorPtr ()
{
  if (nextSolPtrSwitched)
  {
    *(savedNextSolutionPtr) = *(nextSolutionPtr);
    nextSolutionPtr = savedNextSolutionPtr;
    nextSolPtrSwitched = false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setZeroHistory
// Purpose       : Sets everything to zero.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/24/07
//-----------------------------------------------------------------------------
void N_TIA_DataStore::setZeroHistory()
{
#ifdef Xyce_DEBUG_TIME
  Xyce::dout() << "\nN_TIA_DataStore::setZeroHistory" << std::endl;
#endif
  // Solutions:
  nextSolutionPtr->putScalar(0.0);
  nextSolutionDerivPtr->putScalar(0.0);
  nextSolutionDivDiffPtr->putScalar(0.0);

  // States:
  nextStatePtr->putScalar(0.0);
  nextStateDerivPtr->putScalar(0.0);
  nextStateDivDiffPtr->putScalar(0.0);
  nextStorePtr->putScalar(0.0);
  nextStoreLeadCurrQPtr->putScalar(0.0);

  for (int ip=0;ip<nextDfdpPtrVector.size();++ip)
  {
    nextDfdpPtrVector[ip]->putScalar(0.0);
    nextDqdpPtrVector[ip]->putScalar(0.0);
    nextDqdpDerivPtrVector[ip]->putScalar(0.0);
    nextDbdpPtrVector[ip]->putScalar(0.0);
    nextDXdpPtrVector[ip]->putScalar(0.0);
    nextDQdxDXdpPtrVector[ip]->putScalar(0.0);
    nextDQdxDXdpDerivPtrVector[ip]->putScalar(0.0);
  }

  qErrWtVecPtr->putScalar(0.0);

  // DAE formulation vectors
  daeQVectorPtr->putScalar(0.0);
  daeFVectorPtr->putScalar(0.0);
  daeBVectorPtr->putScalar(0.0);

  // Predictors
  xn0Ptr->putScalar(0.0);
  qn0Ptr->putScalar(0.0);
  qpn0Ptr->putScalar(0.0);
  sn0Ptr->putScalar(0.0);
  spn0Ptr->putScalar(0.0);
  stoQn0Ptr->putScalar(0.0);
  stoQpn0Ptr->putScalar(0.0);

  for (int ip=0;ip<dfdpn0PtrVector.size();++ip)
  {
    dfdpn0PtrVector[ip]->putScalar(0.0);
    dqdpn0PtrVector[ip]->putScalar(0.0);
    dbdpn0PtrVector[ip]->putScalar(0.0);
    dXdpn0PtrVector[ip]->putScalar(0.0);

    dfdppn0PtrVector[ip]->putScalar(0.0);
    dqdppn0PtrVector[ip]->putScalar(0.0);
    dbdppn0PtrVector[ip]->putScalar(0.0);
    dXdppn0PtrVector[ip]->putScalar(0.0);
  }

  // Nonlinear solution vector:
  qNewtonCorrectionPtr->putScalar(0.0);
  sNewtonCorrectionPtr->putScalar(0.0);
  stoNewtonCorrectionPtr->putScalar(0.0);
  stoLeadCurrQNewtonCorrectionPtr->putScalar(0.0);

  for (int ip=0;ip<dfdpNewtonCorrectionPtrVector.size();++ip)
  {
    dfdpNewtonCorrectionPtrVector[ip]->putScalar(0.0);
    dqdpNewtonCorrectionPtrVector[ip]->putScalar(0.0);
    dbdpNewtonCorrectionPtrVector[ip]->putScalar(0.0);
    dXdpNewtonCorrectionPtrVector[ip]->putScalar(0.0);
  }

  // Step-size selection temporary vectors
  delta_x->putScalar(0.0);
  delta_q->putScalar(0.0);

  // Temporary vector for WaMPDE interpolation
  tmpXn0APtr->putScalar(0.0);
  tmpXn0BPtr->putScalar(0.0);

  // This just sets the "oldDAE" history vectors to zero.
  setConstantHistory ();

  // new-DAE history:
  int sizeOfHistory = xHistory.size();
  for (int i = 0; i < sizeOfHistory; ++i)
  {
    xHistory[i]->putScalar(0.0);
    qHistory[i]->putScalar(0.0);
    sHistory[i]->putScalar(0.0);
    stoHistory[i]->putScalar(0.0);
    stoLeadCurrQHistory[i]->putScalar(0.0);

    for (int ip=0;ip<dfdpHistory.size();++ip)
    {
      std::vector<N_LAS_Vector *> dfdp = dfdpHistory[ip];
      std::vector<N_LAS_Vector *> dqdp = dqdpHistory[ip];
      std::vector<N_LAS_Vector *> dbdp = dbdpHistory[ip];
      std::vector<N_LAS_Vector *> dXdp = dXdpHistory[ip];
      std::vector<N_LAS_Vector *> dQdxdXdp = dQdxdXdpHistory[ip];

      dfdp[i]->putScalar(0.0);
      dqdp[i]->putScalar(0.0);
      dbdp[i]->putScalar(0.0);
      dXdp[i]->putScalar(0.0);
      dQdxdXdp[i]->putScalar(0.0);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setErrorWtVector_
// Purpose       : Set the Error Weight Vector (defined in terms of the
//                 solution approximation and error tolerances).
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL,Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
void N_TIA_DataStore::setErrorWtVector()
{
  // to avoid several conditionals within a loop we traverse many times, we'll
  // presort the variables into types so we can loop without also making conditional
  // checks
  if( ! indexVecsInitialized )
  {
    // figure out which unknowns are V's, I's or masked.
    numVVars = 0;
    numIVars = 0;
    numMaskedVars = 0;

    bool nTDMF=lasSysPtr->getNonTrivialDeviceMaskFlag();

    // first count how many of each type we have
    for (int k = 0; k < solutionSize; ++k)
    {
      if (nTDMF && (*(deviceMaskPtr))[k] == 0.0)
      {
        numMaskedVars++;
      }
      else if( (tiaParamsPtr_->fastTests == true) && (varTypeVec[k] == 'I') )
      {
        // we only count I vars if we're doing fast tests, otherwise we treat them as v vars
        numIVars++;
      }
      else
      {
        numVVars++;
      }
    }

    // set up our storage
    indexVVars.resize(numVVars);
    indexIVars.resize(numIVars);
    indexMaskedVars.resize(numMaskedVars);

    // now fill the index arrays
    for (int k = 0, currI = 0, currV = 0, currM = 0 ; k < solutionSize; ++k)
    {
      if (nTDMF && (*(deviceMaskPtr))[k] == 0.0)
      {
        indexMaskedVars[currM] = k;
        currM++;
      }
      else if( (tiaParamsPtr_->fastTests == true) && (varTypeVec[k] == 'I') )
      {
        indexIVars[currI] = k;
        currI++;
      }
      else
      {
        indexVVars[currV] = k;
        currV++;
      }
    }

    indexVecsInitialized = true;
  }

#ifdef Xyce_DEBUG_TIME
  char tmp[256];

  if (tiaParamsPtr_->debugLevel > 1)
  {
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "N_TIA_DataStore::setErrorWtVector" << std::endl << std::endl
                 << "   errorWtVector    currSolution     relErrorTol    absErrorTol" << std::endl
                 << "   --------------  --------------  --------------  --------------" << std::endl;
  }
#endif

  if (tiaParamsPtr_->newLte == true)
//  if (tiaParamsPtr_->integrationMethod == 7 && tiaParamsPtr_->newLte == true)
  {
    double currMaxValue = 0.0;
    currSolutionPtr->infNorm(&currMaxValue);

    errWtVecPtr->putScalar(currMaxValue);

#ifdef Xyce_DEBUG_TIME
    if (tiaParamsPtr_->debugLevel > 0)
    {
      std::vector<int> index(1, -1);
      currSolutionPtr->infNormIndex( &index[0] );
      Xyce::dout() << "currMaxValue = " << currMaxValue << ", currMaxValueIndex = " << index[0] << std::endl;
    }
#endif // Xyce_DEBUG_TIME
  }
  else
  {
    errWtVecPtr->absValue(*currSolutionPtr);

#ifdef Xyce_DEBUG_TIME
    if (tiaParamsPtr_->debugLevel > 0)
      {
        double currMaxValue = 0.0;
        currSolutionPtr->infNorm(&currMaxValue);
        std::vector<int> index(1, -1);
        currSolutionPtr->infNormIndex( &index[0] );
        Xyce::dout() << "currMaxValueoldLte = " << currMaxValue << ", currMaxValueIndex = " << index[0]  << std::endl;
      }
#endif
  }

  qErrWtVecPtr->absValue(*daeQVectorPtr);

  if( tiaParamsPtr_->fastTests == true )
  {
    // Voltage variables
    for (int k = 0; k < numVVars; ++k)
    {
//      errWtVecPtr->putScalar(currMaxValue);
      if( (*(errWtVecPtr))[indexVVars[k]] < tiaParamsPtr_->voltZeroTol )
      {
        (*(errWtVecPtr))[indexVVars[k]] = N_UTL_MachineDependentParams::MachineBig();
        (*(qErrWtVecPtr))[indexVVars[k]] = N_UTL_MachineDependentParams::MachineBig();
      }
      else
      {
        (*(errWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(errWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];
        (*(qErrWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(qErrWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];
      }
    }


    // Current variables
    for (int k = 0; k < numIVars; ++k)
    {
//      errWtVecPtr->absValue(*currSolutionPtr);
      if( (*(errWtVecPtr))[indexIVars[k]] < tiaParamsPtr_->currZeroTol )
      {
        (*(errWtVecPtr))[indexIVars[k]] = N_UTL_MachineDependentParams::MachineBig();
        (*(qErrWtVecPtr))[indexIVars[k]] = N_UTL_MachineDependentParams::MachineBig();
      }
      else
      {
        (*(errWtVecPtr))[indexIVars[k]] =  (*(errWtVecPtr))[indexIVars[k]] + (*(absErrTolPtr))[indexIVars[k]];
        (*(qErrWtVecPtr))[indexIVars[k]] = (*(qErrWtVecPtr))[indexIVars[k]] + (*(absErrTolPtr))[indexIVars[k]];
      }
    }

    // Masked variables
    for (int k = 0; k < numMaskedVars; ++k)
    {
      (*(errWtVecPtr))[indexMaskedVars[k]] = (*(qErrWtVecPtr))[indexMaskedVars[k]] = N_UTL_MachineDependentParams::MachineBig();
    }
  }
  else
  {
     // Voltage variables
    for (int k = 0; k < numVVars; ++k)
    {
//      errWtVecPtr->putScalar(currMaxValue);
      (*(errWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(errWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];
      (*(qErrWtVecPtr))[indexVVars[k]] = (*(relErrTolPtr))[indexVVars[k]] * (*(qErrWtVecPtr))[indexVVars[k]] + (*(absErrTolPtr))[indexVVars[k]];
    }

    // Current variables
    // if fastTests == false, then I vars are treated with V vars above.

    // Masked variables
    for (int k = 0; k < numMaskedVars; ++k)
    {
      (*(errWtVecPtr))[indexMaskedVars[k]] = (*(qErrWtVecPtr))[indexMaskedVars[k]] = N_UTL_MachineDependentParams::MachineBig();
    }

  }

#ifdef Xyce_DEBUG_TIME
  for (int k = 0; k < solutionSize; ++k)
  {
    if (tiaParamsPtr_->debugLevel > 1)
    {
      sprintf(tmp,"%16.6e%16.6e%16.6e%16.6e",
              (*(errWtVecPtr))     [k],
              (*(currSolutionPtr)) [k],
              (*(relErrTolPtr))    [k],
              (*(absErrTolPtr))    [k]);

      Xyce::dout() << tmp << std::endl;
    }
  }
#endif

#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    Xyce::dout() << ""  << std::endl
                 << Xyce::section_divider << std::endl;
  }
#endif
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::WRMS_errorNorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/25/05
//-----------------------------------------------------------------------------
double N_TIA_DataStore::WRMS_errorNorm()
{
  double errorNorm = 0.0, qErrorNorm = 0.0;
  newtonCorrectionPtr->wRMSNorm(*errWtVecPtr, &errorNorm);
  qNewtonCorrectionPtr->wRMSNorm(*qErrWtVecPtr, &qErrorNorm);

#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    Xyce::dout() << "N_TIA_DataStore::errorNorm = " << errorNorm << std::endl;
    Xyce::dout() << "N_TIA_DataStore::qErrorNorm = " << qErrorNorm << std::endl;
  }
#endif

  // This is for handling the total errorNorm for 2-level solves.
  //
  // Note:  This version of the function assumes that the error norm
  // and q-error norm are the same size.  (they have to be...)
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;
    double totalQSum = qErrorNorm*qErrorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum;
      double innerQSum = innerErrorInfoVec[i].qErrorSum;
      double innerSize = innerErrorInfoVec[i].innerSize;

#ifdef Xyce_DEBUG_TIME
      Xyce::dout() << "DSdae:innerSum["<<i<<"] = " << innerSum <<std::endl;
      Xyce::dout() << "DSdae:innerQSum["<<i<<"] = " << innerQSum <<std::endl;
      Xyce::dout() << "DSdae:innerSize["<<i<<"] = " << innerSize <<std::endl;
#endif

      totalSize += innerSize;
      totalSum += innerSum;
      totalQSum += innerQSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
    qErrorNorm = sqrt(recip*totalQSum);

#ifdef Xyce_DEBUG_TIME
    Xyce::dout() << "DSdae:upperSize = " << upperSize << std::endl;
    Xyce::dout() << "DSdae:totalSum = " << totalSum << std::endl;
    Xyce::dout() << "DSdae:totalQSum = " << totalQSum << std::endl;
    Xyce::dout() << "DSdae:totalSize = " << totalSize << std::endl;
    Xyce::dout() << "DSdae:2-level errorNorm = " << errorNorm << std::endl;
    Xyce::dout() << "DSdae:2-level qErrorNorm = " << qErrorNorm << std::endl;
#endif
  }

#ifndef Xyce_USE_Q_NORM
  //errorNorm = errorNorm;
#else
  errorNorm = sqrt(0.5*errorNorm*errorNorm+0.5*qErrorNorm*qErrorNorm);
#endif
  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialQErrorNormSum
// Purpose       : Needed by 2-level solves.  This is the Q-vector version
//                 of this function.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialQErrorNormSum ()
{
  double qErrorNorm = 0.0;
  qNewtonCorrectionPtr->wRMSNorm(*qErrWtVecPtr, &qErrorNorm);
  double sum = qErrorNorm*qErrorNorm;
  double length = qNewtonCorrectionPtr->globalLength();
  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_m1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_m1 (int currentOrder)
{
  double sum = 0.0;

  if (currentOrder>1)
  {
    delta_x->linearCombo(1.0,*(xHistory[currentOrder]),1.0,*newtonCorrectionPtr);
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_m2
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_m2 (int currentOrder)
{
  double sum = 0.0;

  if (currentOrder>2)
  {
    delta_x->linearCombo(1.0,*(xHistory[currentOrder]),1.0,*newtonCorrectionPtr);
    delta_x->linearCombo(1.0,*(xHistory[currentOrder-1]),1.0,*delta_x);
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_p1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_p1 (int currentOrder, int maxOrder)
{
  double sum = 0.0;

  if (currentOrder<maxOrder)
  {
    delta_x->linearCombo(1.0,*newtonCorrectionPtr,-1.0,*(xHistory[currentOrder+1]));
    double norm = 0.0;
    delta_x->wRMSNorm(*errWtVecPtr, &norm);
    sum = norm*norm;
    double length = newtonCorrectionPtr->globalLength();
    sum *= length;
  }

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::partialSum_q1
// Purpose       : Needed by 2-level solves.
//
// Special Notes : A weighted RMS norm is this:
//
//                  norm = sqrt ( 1/n * sum (x/w)^2 )
//
//                  What this function returns is: sum (x/w)^2
//
//                  It will later be summed with other patial sums
//                  to get the complete WRMS value.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::partialSum_q1 ()
{
  double sum = 0.0;

  double norm = 0.0;
  (qHistory[1])->wRMSNorm(*qErrWtVecPtr, &norm);

  sum = norm*norm;
  double length = qNewtonCorrectionPtr->globalLength();
  sum *= length;

  return sum;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_m1
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_m1()
{
  double errorNorm = 0.0;
  delta_x->wRMSNorm(*errWtVecPtr, &errorNorm);

  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum_m1;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_m2
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_m2()
{
  double errorNorm = 0.0;
  delta_x->wRMSNorm(*errWtVecPtr, &errorNorm);

  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum_m2;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_p1
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_p1()
{
  double errorNorm = 0.0;
  delta_x->wRMSNorm(*errWtVecPtr, &errorNorm);

  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].xErrorSum_p1;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::delta_x_errorNorm_q1
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 03/15/07
//-----------------------------------------------------------------------------
double N_TIA_DataStore::delta_x_errorNorm_q1()
{
  double errorNorm = 0.0;
  (qHistory[1])->wRMSNorm(*qErrWtVecPtr, &errorNorm);

  // This is for handling the total errorNorm for 2-level solves.
  // It really should get a separate overloaded function.
  if ( !(innerErrorInfoVec.empty()) )
  {
    double upperSize = newtonCorrectionPtr->globalLength();
    int sumSize = innerErrorInfoVec.size();

    double totalSize = upperSize;
    double totalSum = errorNorm*errorNorm*upperSize;

    for (int i=0;i<sumSize;++i)
    {
      double innerSum = innerErrorInfoVec[i].q1HistorySum;
      double innerSize = innerErrorInfoVec[i].innerSize;

      totalSize += innerSize;
      totalSum += innerSum;
    }

    double recip = 1.0/totalSize;
    errorNorm = sqrt(recip*totalSum);
  }

  return errorNorm;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::stepLinearCombo
// Purpose       : setup the newtonCorrection vectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 02/20/07
//-----------------------------------------------------------------------------
void N_TIA_DataStore::stepLinearCombo()
{
  // 03/16/04 tscoffe:  update the newton correction.  Note:  this should be
  // available from NOX, but for now I'm going to do the difference anyway.
  newtonCorrectionPtr->linearCombo (1.0,*nextSolutionPtr,-1.0,*xn0Ptr);

  // We need to compute the correction in Q here
  // I'm assuming dsDaePtr_->daeQVectorPtr will be fresh from the end of the
  // nonlinear solve.
  qNewtonCorrectionPtr->linearCombo (1.0,*daeQVectorPtr,-1.0,*qn0Ptr);

  // We also need a State correction between the time steps
  sNewtonCorrectionPtr->linearCombo (1.0,*nextStatePtr,-1.0,*sn0Ptr);

  // We also need a Store correction between the time steps
  stoNewtonCorrectionPtr->linearCombo (1.0,*nextStorePtr,-1.0,*ston0Ptr);
  stoLeadCurrQNewtonCorrectionPtr->linearCombo (1.0, *nextStoreLeadCurrQPtr, -1.0, *stoQn0Ptr);

  // correction on sensitivities:
  for (int ip=0;ip<dfdpNewtonCorrectionPtrVector.size();++ip)
  {
    (dfdpNewtonCorrectionPtrVector[ip])->linearCombo (1.0,*(nextDfdpPtrVector[ip]),-1.0,*(dfdpn0PtrVector[ip]));
    (dqdpNewtonCorrectionPtrVector[ip])->linearCombo (1.0,*(nextDqdpPtrVector[ip]),-1.0,*(dqdpn0PtrVector[ip]));
    (dbdpNewtonCorrectionPtrVector[ip])->linearCombo (1.0,*(nextDbdpPtrVector[ip]),-1.0,*(dbdpn0PtrVector[ip]));
    (dXdpNewtonCorrectionPtrVector[ip])->linearCombo (1.0,*(nextDXdpPtrVector[ip]),-1.0,*(dXdpn0PtrVector[ip]));
  }

#ifdef Xyce_DEBUG_TIME
  if (tiaParamsPtr_->debugLevel > 1)
  {
    Xyce::dout() << "\n newtonCorrection: \n" << std::endl;
    newtonCorrectionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
    Xyce::dout() << "\n qNewtonCorrection: \n" << std::endl;
    qNewtonCorrectionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << "\n sNewtonCorrection: \n" << std::endl;
    sNewtonCorrectionPtr->printPetraObject(Xyce::dout());
    Xyce::dout() << std::endl;
  }
#endif // Xyce_DEBUG_TIME

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::getSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::getSolnVarData( const int & gid,
				      std::vector<double> & varData )
{
  varData.resize(23);
  int i=0;
  varData[i++] = tmpSolVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = currSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextSolutionPtr->getElementByGlobalIndex( gid );
  varData[i++] = currSolutionDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastSolutionDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextSolutionDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = currSolutionDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastSolutionDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextSolutionDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = errWtVecPtr->getElementByGlobalIndex( gid );
  varData[i++] = absErrTolPtr->getElementByGlobalIndex( gid );
  varData[i++] = relErrTolPtr->getElementByGlobalIndex( gid );
  varData[i++] = newtonCorrectionPtr->getElementByGlobalIndex( gid );
  varData[i++] = qErrWtVecPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeQVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeFVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = daeBVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = xn0Ptr->getElementByGlobalIndex ( gid );
  varData[i++] = qn0Ptr->getElementByGlobalIndex ( gid );
  varData[i++] = qpn0Ptr->getElementByGlobalIndex ( gid );
  varData[i++] = dFdxdVpVectorPtr->getElementByGlobalIndex ( gid );
  varData[i++] = dQdxdVpVectorPtr->getElementByGlobalIndex ( gid );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::getStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::getStateVarData( const int & gid,
			 	       std::vector<double> & varData )
{
  int i=0;
  varData.resize( 14 );
  varData[i++] = tmpStaVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = tmpStaDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = tmpStaDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStatePtr->getElementByGlobalIndex( gid );
  varData[i++] = currStateDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStateDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStateDerivPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStateDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStateDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStateDivDiffPtr->getElementByGlobalIndex( gid );
  varData[i++] = sn0Ptr->getElementByGlobalIndex( gid );
  varData[i++] = spn0Ptr->getElementByGlobalIndex( gid );
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::getStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date :
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::getStoreVarData( const int & gid,
			 	       std::vector<double> & varData )
{
  int i=0;
  varData.resize( 6 );
  varData[i++] = tmpStoVectorPtr->getElementByGlobalIndex( gid );
  varData[i++] = currStorePtr->getElementByGlobalIndex( gid );
  varData[i++] = lastStorePtr->getElementByGlobalIndex( gid );
  varData[i++] = nextStorePtr->getElementByGlobalIndex( gid );
  varData[i++] = ston0Ptr->getElementByGlobalIndex( gid );
  varData[i++] = stopn0Ptr->getElementByGlobalIndex( gid );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setSolnVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setSolnVarData( const int & gid,
				      const std::vector<double> & varData )
{
  int i=0;
  tmpSolVectorPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextSolutionPtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currSolutionDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  lastSolutionDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  nextSolutionDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  currSolutionDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  lastSolutionDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  nextSolutionDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  errWtVecPtr->setElementByGlobalIndex           ( gid, varData[i++] );
  absErrTolPtr->setElementByGlobalIndex          ( gid, varData[i++] );
  relErrTolPtr->setElementByGlobalIndex          ( gid, varData[i++] );
  newtonCorrectionPtr->setElementByGlobalIndex   ( gid, varData[i++] );
  qErrWtVecPtr->setElementByGlobalIndex          ( gid, varData[i++] );
  daeQVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  daeFVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  daeBVectorPtr->setElementByGlobalIndex         ( gid, varData[i++] );
  xn0Ptr->setElementByGlobalIndex                ( gid, varData[i++] );
  qn0Ptr->setElementByGlobalIndex                ( gid, varData[i++] );
  qpn0Ptr->setElementByGlobalIndex               ( gid, varData[i++] );
  dFdxdVpVectorPtr->setElementByGlobalIndex      ( gid, varData[i++] );
  dQdxdVpVectorPtr->setElementByGlobalIndex      ( gid, varData[i++] );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setStateVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setStateVarData( const int & gid,
			 	       const std::vector<double> & varData )
{
  int i=0;
  tmpStaVectorPtr->setElementByGlobalIndex    ( gid, varData[i++] );
  tmpStaDerivPtr->setElementByGlobalIndex     ( gid, varData[i++] );
  tmpStaDivDiffPtr->setElementByGlobalIndex   ( gid, varData[i++] );
  currStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextStatePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  currStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  lastStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  nextStateDerivPtr->setElementByGlobalIndex  ( gid, varData[i++] );
  currStateDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  lastStateDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );
  nextStateDivDiffPtr->setElementByGlobalIndex( gid, varData[i++] );

  sn0Ptr->setElementByGlobalIndex             ( gid, varData[i++] );
  spn0Ptr->setElementByGlobalIndex            ( gid, varData[i++] );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_TIA_DataStore::setStoreVarData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/27/06
//-----------------------------------------------------------------------------
bool N_TIA_DataStore::setStoreVarData( const int & gid,
			 	       const std::vector<double> & varData )
{
  int i=0;
  tmpStoVectorPtr->setElementByGlobalIndex    ( gid, varData[i++] );
  currStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  lastStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  nextStorePtr->setElementByGlobalIndex       ( gid, varData[i++] );
  ston0Ptr->setElementByGlobalIndex           ( gid, varData[i++] );
  stopn0Ptr->setElementByGlobalIndex          ( gid, varData[i++] );
  return true;
}

