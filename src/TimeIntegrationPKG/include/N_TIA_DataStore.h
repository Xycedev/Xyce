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
// Filename      : $RCSfile: N_TIA_DataStore.h,v $
//
// Purpose       : This file handles the class that defines the data arrays
//                 needed for the time integration algorithms.
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
// Revision Number: $Revision: 1.90.2.4 $
//
// Revision Date  : $Date: 2014/09/02 22:49:49 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_DATA_STORE_H
#define Xyce_N_TIA_DATA_STORE_H

// ---------- Standard Includes ----------

#include <list>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_LAS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_TIA_TwoLevelError.h>

//-----------------------------------------------------------------------------
// Class         : N_TIA_DataStore
// Purpose       : This is the class for defining data arrays needed in the
//                 time integration algorithms.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_DataStore
{
  public:
    N_TIA_DataStore(N_TIA_TIAParams * tiaPtr, N_LAS_System * lsPtr);
    ~N_TIA_DataStore();

    void allocateSensitivityArrays(int numParams);
    void deleteSensitivityArrays();

  private:
    N_TIA_DataStore(const N_TIA_DataStore& rhs);
    N_TIA_DataStore &operator=(const N_TIA_DataStore& rhs);
    N_TIA_DataStore();

  public:
    // TIA Arrays (pointers) for  Integration Solution Process:
    unsigned int solutionSize;
    unsigned int stateSize;

    // temporary vectors:
    N_LAS_Vector * tmpSolVectorPtr;
    N_LAS_Vector * tmpStaVectorPtr;
    N_LAS_Vector * tmpStaDerivPtr;
    N_LAS_Vector * tmpStaDivDiffPtr;

    N_LAS_Vector * tmpStoVectorPtr;

    // Predictors
    N_LAS_Vector * xn0Ptr;

    // Solutions:
    N_LAS_Vector * currSolutionPtr;
    N_LAS_Vector * lastSolutionPtr;
    N_LAS_Vector * oldeSolutionPtr;
    N_LAS_Vector * nextSolutionPtr;
    N_LAS_Vector * flagSolutionPtr;

    N_LAS_Vector * savedNextSolutionPtr;

    // States:
    N_LAS_Vector * currStatePtr;
    N_LAS_Vector * lastStatePtr;
    N_LAS_Vector * oldeStatePtr;
    N_LAS_Vector * nextStatePtr;

    // Storage:
    N_LAS_Vector * currStorePtr;
    N_LAS_Vector * lastStorePtr;
    N_LAS_Vector * oldeStorePtr;
    N_LAS_Vector * nextStorePtr;
    // for lead current calculations.  F component is
    // held in the store vector, Q component is here
    N_LAS_Vector * currStoreLeadCurrQPtr;
    N_LAS_Vector * lastStoreLeadCurrQPtr;
    N_LAS_Vector * oldeStoreLeadCurrQPtr;
    N_LAS_Vector * nextStoreLeadCurrQPtr;

    // sensitivity vectors:
    std::vector<N_LAS_Vector*> sensRHSPtrVector;

    std::vector<N_LAS_Vector*> currDfdpPtrVector;
    std::vector<N_LAS_Vector*> lastDfdpPtrVector;
    std::vector<N_LAS_Vector*> oldeDfdpPtrVector;
    std::vector<N_LAS_Vector*> nextDfdpPtrVector;

    std::vector<N_LAS_Vector*> currDqdpPtrVector;
    std::vector<N_LAS_Vector*> lastDqdpPtrVector;
    std::vector<N_LAS_Vector*> oldeDqdpPtrVector;
    std::vector<N_LAS_Vector*> nextDqdpPtrVector;

    std::vector<N_LAS_Vector*> currDbdpPtrVector;
    std::vector<N_LAS_Vector*> lastDbdpPtrVector;
    std::vector<N_LAS_Vector*> oldeDbdpPtrVector;
    std::vector<N_LAS_Vector*> nextDbdpPtrVector;

    std::vector<N_LAS_Vector*> currDXdpPtrVector;
    std::vector<N_LAS_Vector*> lastDXdpPtrVector;
    std::vector<N_LAS_Vector*> oldeDXdpPtrVector;
    std::vector<N_LAS_Vector*> nextDXdpPtrVector;

    // matvecs:
    std::vector<N_LAS_Vector*> currDQdxDXdpPtrVector;
    std::vector<N_LAS_Vector*> lastDQdxDXdpPtrVector;
    std::vector<N_LAS_Vector*> oldeDQdxDXdpPtrVector;
    std::vector<N_LAS_Vector*> nextDQdxDXdpPtrVector;

    std::vector<N_LAS_Vector*> currDQdxDXdpDerivPtrVector;
    std::vector<N_LAS_Vector*> lastDQdxDXdpDerivPtrVector;
    std::vector<N_LAS_Vector*> oldeDQdxDXdpDerivPtrVector;
    std::vector<N_LAS_Vector*> nextDQdxDXdpDerivPtrVector;

    // Derivatives of Solutions:
    N_LAS_Vector * currSolutionDerivPtr;
    N_LAS_Vector * lastSolutionDerivPtr;
    N_LAS_Vector * oldeSolutionDerivPtr;
    N_LAS_Vector * nextSolutionDerivPtr;

    // Derivatives of States:
    N_LAS_Vector * currStateDerivPtr;
    N_LAS_Vector * lastStateDerivPtr;
    N_LAS_Vector * oldeStateDerivPtr;
    N_LAS_Vector * nextStateDerivPtr;

    // Derivatives of Store for lead curent calculations
    N_LAS_Vector * currStoreLeadCurrQDerivPtr;
    N_LAS_Vector * lastStoreLeadCurrQDerivPtr;
    N_LAS_Vector * oldeStoreLeadCurrQDerivPtr;
    N_LAS_Vector * nextStoreLeadCurrQDerivPtr;

    // Derivatives of dq/dp for sensitivity calculations
    std::vector<N_LAS_Vector*> currDqdpDerivPtrVector;
    std::vector<N_LAS_Vector*> lastDqdpDerivPtrVector;
    std::vector<N_LAS_Vector*> oldeDqdpDerivPtrVector;
    std::vector<N_LAS_Vector*> nextDqdpDerivPtrVector;

    // Sensitivity residual matvec term:
    N_LAS_Vector * resMatVecPtr;

    // Scaled Divided Differences
    N_LAS_Vector * currSolutionDivDiffPtr;
    N_LAS_Vector * lastSolutionDivDiffPtr;
    N_LAS_Vector * oldeSolutionDivDiffPtr;
    N_LAS_Vector * nextSolutionDivDiffPtr;

    // Scaled Divided Differences
    N_LAS_Vector * currStateDivDiffPtr;
    N_LAS_Vector * lastStateDivDiffPtr;
    N_LAS_Vector * oldeStateDivDiffPtr;
    N_LAS_Vector * nextStateDivDiffPtr;

    // Error Vectors
    N_LAS_Vector * errWtVecPtr;
    N_LAS_Vector * absErrTolPtr;
    N_LAS_Vector * relErrTolPtr;

    // Jacobian and RHS
    N_LAS_Matrix * JMatrixPtr;
    N_LAS_Vector * RHSVectorPtr;
#ifdef Xyce_DEBUG_DEVICE
    N_LAS_Vector * JdxpVectorPtr; // only used in the device manager for debug purposes these days.
#endif

    // NonLinear Solution Vectors
    N_LAS_Vector * newtonCorrectionPtr;

    // Mask for error norms (to allow some equations not to take part in
    // weighted norms)
    N_LAS_Vector * deviceMaskPtr;
    // TVR: I toyed with the idea of having this flag here, but went with
    // keeping it in the LAS_System instead --- we call an accessor method
    // to set and get the flag when we create or use the mask.
    // flag showing whether mask is trivial or not
    //   bool nonTrivialDeviceMask;

    // this is a simple vector indicating var types.  Now
    // we just handle V and I differently, but we could do more with this
    std::vector<char> varTypeVec;

    // To remove conditionals from setErrorWtVector() we'll create
    // lists of indexes of unknows that are handled in different ways
    std::vector<int> indexVVars;
    std::vector<int> indexIVars;
    std::vector<int> indexMaskedVars;
    int numVVars, numIVars, numMaskedVars;
    bool indexVecsInitialized;

    // limiter flag:
    bool limiterFlag;

    // 2-level information:
    std::vector<N_TIA_TwoLevelError> innerErrorInfoVec;

    // new-DAE data (originally from the new-DAE derrived class)
    // Error Vectors
    N_LAS_Vector * qErrWtVecPtr;

    // DAE formulation vectors
    N_LAS_Vector * daeQVectorPtr;
    N_LAS_Vector * daeFVectorPtr;
    N_LAS_Vector * daeBVectorPtr;

    // voltage limiting vectors
    N_LAS_Vector * dFdxdVpVectorPtr;
    N_LAS_Vector * dQdxdVpVectorPtr;

    // DAE formulation matrices
    N_LAS_Matrix * dQdxMatrixPtr;
    N_LAS_Matrix * dFdxMatrixPtr;

    // temporary stop-gap, needed for sensitivity analysis.
    //N_LAS_Matrix * saved_dQdxMatrixPtr;

    // HB temporary Matvec storage vectors
    N_LAS_Vector * dQdxVecVectorPtr;
    N_LAS_Vector * dFdxVecVectorPtr;

    // History arrays
    std::vector<N_LAS_Vector*> xHistory;
    std::vector<N_LAS_Vector*> qHistory;
    std::vector<N_LAS_Vector*> sHistory;    // state history
    std::vector<N_LAS_Vector*> stoHistory;  // store history
    std::vector<N_LAS_Vector*> stoLeadCurrQHistory;  // store history for lead current Q component.

    // sensitivity histories
    std::vector< std::vector<N_LAS_Vector*> > dfdpHistory;
    std::vector< std::vector<N_LAS_Vector*> > dqdpHistory;
    std::vector< std::vector<N_LAS_Vector*> > dbdpHistory;
    std::vector< std::vector<N_LAS_Vector*> > dXdpHistory;

    std::vector< std::vector<N_LAS_Vector*> > dQdxdXdpHistory;

    // Predictors
    N_LAS_Vector * qn0Ptr;
    N_LAS_Vector * qpn0Ptr;

    N_LAS_Vector * sn0Ptr;
    N_LAS_Vector * spn0Ptr;

    N_LAS_Vector * ston0Ptr;
    N_LAS_Vector * stopn0Ptr;

    N_LAS_Vector * stoQn0Ptr;
    N_LAS_Vector * stoQpn0Ptr;

    // For sensitivities
    std::vector<N_LAS_Vector*> dfdpn0PtrVector;
    std::vector<N_LAS_Vector*> dqdpn0PtrVector;
    std::vector<N_LAS_Vector*> dbdpn0PtrVector;
    std::vector<N_LAS_Vector*> dXdpn0PtrVector;

    std::vector<N_LAS_Vector*> dfdppn0PtrVector;
    std::vector<N_LAS_Vector*> dqdppn0PtrVector;
    std::vector<N_LAS_Vector*> dbdppn0PtrVector;
    std::vector<N_LAS_Vector*> dXdppn0PtrVector;

    // Nonlinear solution vector:
    N_LAS_Vector * qNewtonCorrectionPtr;
    N_LAS_Vector * sNewtonCorrectionPtr;
    N_LAS_Vector * stoNewtonCorrectionPtr;
    N_LAS_Vector * stoLeadCurrQNewtonCorrectionPtr;

    // For sensitivities
    std::vector<N_LAS_Vector*> dfdpNewtonCorrectionPtrVector;
    std::vector<N_LAS_Vector*> dqdpNewtonCorrectionPtrVector;
    std::vector<N_LAS_Vector*> dbdpNewtonCorrectionPtrVector;
    std::vector<N_LAS_Vector*> dXdpNewtonCorrectionPtrVector;

    // Step-size selection temporary vectors
    N_LAS_Vector * delta_x;
    N_LAS_Vector * delta_q;

    // Temporary vectors for WaMPDE interpolation
    N_LAS_Vector * tmpXn0APtr;
    N_LAS_Vector * tmpXn0BPtr;

    // These are for MPDE fast time scale points
    std::vector<double> timeSteps;
    std::vector<bool> timeStepsBreakpointFlag;
    std::vector<N_LAS_Vector*> fastTimeSolutionVec;
    std::vector<N_LAS_Vector*> fastTimeStateVec;
    std::vector<N_LAS_Vector*> fastTimeQVec;
    std::vector<N_LAS_Vector*> fastTimeStoreVec;

    std::vector<double> objectiveVec_;
    std::vector<double> dOdpVec_;
    std::vector<double> dOdpAdjVec_;
    std::vector<double> scaled_dOdpVec_;
    std::vector<double> scaled_dOdpAdjVec_;
    std::vector<double> paramOrigVals_; 

  // DataStore Functions
  public:
    void initializeDataArrays();
    void enableOrderOneStart();

    void updateSolDataArrays();
    bool updateStateDataArrays();

    void outputSolDataArrays(std::ostream &os);
    void setConstantHistory();
    void setConstantSensitivityHistory();
    void setZeroHistory();
    void setErrorWtVector();
    double WRMS_errorNorm();

    double partialErrorNormSum();
    double partialQErrorNormSum();

    double partialSum_m1(int currentOrder);
    double partialSum_m2(int currentOrder);

    double globalLength ();
    void computeDividedDifferences();

    void computeDivDiffsBlock(const std::list<index_pair> & solGIDList,
                              const std::list<index_pair> & staGIDList);

    void printOutPointers ();
    bool equateTmpVectors ();
    bool usePreviousSolAsPredictor ();

    void outputPredictedSolution(std::ostream &os);
    void outputPredictedDerivative(std::ostream &os);

    void stepLinearCombo ();

    double partialSum_p1(int currentOrder, int maxOrder);
    double partialSum_q1();

    double delta_x_errorNorm_m1();
    double delta_x_errorNorm_m2();
    double delta_x_errorNorm_p1();
    double delta_x_errorNorm_q1();

    bool getSolnVarData( const int & gid, std::vector<double> & varData );
    bool getStateVarData( const int & gid, std::vector<double> & varData );
    bool setSolnVarData( const int & gid, const std::vector<double> & varData );
    bool setStateVarData( const int & gid, const std::vector<double> & varData );
    bool getStoreVarData( const int & gid, std::vector<double> & varData );
    bool setStoreVarData( const int & gid, const std::vector<double> & varData );

    N_LAS_System * lasSysPtr;

    bool setNextSolVectorPtr (N_LAS_Vector * solVecPtr);
    bool unsetNextSolVectorPtr ();

    bool resetAll ();
    bool resetFastTimeData ();

  protected:
    N_TIA_TIAParams * tiaParamsPtr_;
    double          * dataBlockPtr;
    bool nextSolPtrSwitched;
};

#endif // Xyce_N_TIA_DATA_STORE_H

