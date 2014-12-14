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
// Filename       : $RCSfile: N_NLS_Sensitivity.C,v $
//
// Purpose        : Body for the sensitivity class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.95.2.8 $
//
// Revision Date  : $Date $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

#include <N_UTL_Misc.h>
#include <algorithm>
#include <sstream>
#include <stdexcept>

// ----------   Xyce Includes   ----------

#include <N_NLS_Sensitivity.h>
#include <N_NLS_Manager.h>

#include <N_LOA_Loader.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>

#include <N_ERH_ErrorMgr.h>
#include <N_ANP_AnalysisManager.h>

#include <N_TOP_Topology.h>

#include <N_UTL_Expression.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_OptionBlock.h>

#include <N_TIA_DataStore.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>

// ----------   Static Declarations ----------

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : Sensitivity::Sensitivity
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
Sensitivity::Sensitivity (NonLinearSolver & nls,
                                      N_TOP_Topology & topTmp,
                                      N_IO_CmdParse & cp)
    : NonLinearSolver(cp),
      debugLevel_(1),
      solutionSize_(0),
      solveDirectFlag_(false),
      solveAdjointFlag_(true),
      outputScaledFlag_(false),
      outputUnscaledFlag_(true),
      maxParamStringSize_(0),
      stdOutputFlag_(true),
      fileOutputFlag_(false),
      dakotaFileOutputFlag_(false),
      forceFD_(false),
      numSolves_(0),
      difference(SENS_FWD),
      objFuncGiven_(false),
      objFuncGIDsetup_(false),
      expNumVars_(0),
      expVal_(0.0),
      objFuncString_(""),
      curValue_(0.0),
      objFuncEval_(0.0),
      dOdp_(0.0),
      sqrtEta_(1.0e-8),
      sqrtEtaGiven_(false),
      dOdXVectorPtr_(0),
      lambdaVectorPtr_(0),
      savedRHSVectorPtr_(0),
      savedNewtonVectorPtr_(0),

      origFVectorPtr_(0),
      pertFVectorPtr_(0),
      origQVectorPtr_(0),
      pertQVectorPtr_(0),
      origBVectorPtr_(0),
      pertBVectorPtr_(0),

      nls_(nls),
      top_(topTmp),
      expPtr_(0),
      numSensParams_(0)
{
  // if the base nonlinear solver class had a copy constructor, I could use
  // that here... maybe I'll set that up later. ERK 11/15/02.
  lasSysPtr_    = nls_.lasSysPtr_;
  anaIntPtr_    = nls_.anaIntPtr_;
  loaderPtr_    = nls_.loaderPtr_;
  rhsVectorPtr_ = nls_.rhsVectorPtr_;

  NewtonVectorPtr_     = nls_.NewtonVectorPtr_;
  lasSolverPtr_        = nls_.lasSolverPtr_;
  jacobianMatrixPtr_   = nls_.jacobianMatrixPtr_;
  nextSolVectorPtrPtr_ = nls_.nextSolVectorPtrPtr_;
  currSolVectorPtrPtr_ = nls_.currSolVectorPtrPtr_;

  dOdXVectorPtr_        = lasSysPtr_->builder().createVector();
  savedRHSVectorPtr_    = lasSysPtr_->builder().createVector();
  savedNewtonVectorPtr_ = lasSysPtr_->builder().createVector();

  origFVectorPtr_ = lasSysPtr_->builder().createVector();
  pertFVectorPtr_ = lasSysPtr_->builder().createVector();

  origQVectorPtr_ = lasSysPtr_->builder().createVector();
  pertQVectorPtr_ = lasSysPtr_->builder().createVector();

  origBVectorPtr_ = lasSysPtr_->builder().createVector();
  pertBVectorPtr_ = lasSysPtr_->builder().createVector();

  solutionSize_ = lasSysPtr_->getSolutionSize();
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::~Sensitivity
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
Sensitivity::~Sensitivity()
{
  delete dOdXVectorPtr_;
  dOdXVectorPtr_ = 0;

  delete savedRHSVectorPtr_;
  savedRHSVectorPtr_ = 0;

  delete savedNewtonVectorPtr_;
  savedNewtonVectorPtr_ = 0;

  delete origFVectorPtr_;
  origFVectorPtr_ = 0;

  delete pertFVectorPtr_;
  pertFVectorPtr_ = 0;

  delete origQVectorPtr_;
  origQVectorPtr_ = 0;

  delete pertQVectorPtr_;
  pertQVectorPtr_ = 0;


  delete origBVectorPtr_;
  origBVectorPtr_ = 0;

  delete pertBVectorPtr_;
  pertBVectorPtr_ = 0;

  if (lambdaVectorPtr_) // only allocated for adjoint solves
  {
    delete lambdaVectorPtr_;
    lambdaVectorPtr_ = 0;
  }

  if (expPtr_)
  {
    delete expPtr_;
    expPtr_ = 0;
  }

  // For all the stuff that is to be deleted in the nonlinear solver
  // base class destructor, just set those pointers to zero because
  // they will have been deleted already.
  //
  // This is the consequence of having the sensitivity class derived
  // off of the nonlinear solver base class, but having it use all
  // the same linear solver objects, etc., as the nonlinear solver used
  // to solve the problem.
  NewtonVectorPtr_     = 0;
  gradVectorPtr_       = 0;
  solWtVectorPtr_      = 0;
  petraOptionBlockPtr_ = 0;
  lasSolverPtr_        = 0;
}



//-----------------------------------------------------------------------------
// Function      : Sensitivity::stdOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/18/2014
//-----------------------------------------------------------------------------
void Sensitivity::stdOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities)
{
  // Send the sensitivity information to the screen:
#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  int myPID = pdsCommPtr->procID();
  if (myPID==0)
#endif
  {
    Xyce::dout() << "\n"<<idString << " Sensitivities of objective function:" 
        << objFuncString_ << " = " << objFuncEval_ << std::endl;

    Xyce::dout() << std::setw(maxParamStringSize_)<<"Name";
    Xyce::dout() << "\t"<<std::setw(13)<<"Value";
    Xyce::dout() << "\t"<<std::setw(13)<<"Sensitivity";
    Xyce::dout() << "\t"<<std::setw(13)<<"Normalized"<<std::endl;

    for (int iparam=0; iparam< numSensParams_; ++iparam)
    {
      Xyce::dout() << std::setw(maxParamStringSize_)<<paramNameVec_[iparam];

      Xyce::dout() << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << paramVals[iparam];

      Xyce::dout() << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << sensitivities[iparam];

      Xyce::dout() << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << scaled_sensitivities[iparam] << std::endl;
    }
  }
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr->barrier();
#endif
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::fileOutput
// Purpose       : Dump sensitivity information out to a file.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/18/2014
//-----------------------------------------------------------------------------
void Sensitivity::fileOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities)
{

#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  int myPID = pdsCommPtr->procID();
  if (myPID==0)
#endif
  {
    std::ostringstream numSolvesOStr;
    numSolvesOStr << numSolves_;
    std::string dodpFileName = netlistFileName_ + numSolvesOStr.str() + "_dodp" + idString +".txt";
    FILE *fp = fopen(dodpFileName.c_str(),"w");
    for (int iparam=0;iparam< numSensParams_; ++iparam)
    {
      fprintf(fp,"\t%16.8e\n", sensitivities[iparam]);
    }
    fclose(fp);
  }
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr->barrier();
#endif
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::dakOutput
// Purpose       : Dump sensitivity information out to a dakota-style file.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/19/2014
//-----------------------------------------------------------------------------
void Sensitivity::dakOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities)
{
#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  int myPID = pdsCommPtr->procID();
  if (myPID==0)
#endif
  {
    // write a file format that can be used by dakota :
    // Note that currently, this will simply overwrite the same 
    // file every time this function is called.
    std::string dakotaFileName = netlistFileName_ + "_dodp" + idString + "_all.txt";
    FILE *fp2 = fopen(dakotaFileName.c_str(),"w");
    fprintf(fp2,"%16.8e", objFuncEval_ );
    fprintf(fp2,"%s","\n[\n");
    for (int iparam=0;iparam< numSensParams_; ++iparam)
    {
      fprintf(fp2,"\t%16.8e\n", sensitivities[iparam]);
    }
    fprintf(fp2,"%s","]\n");
    fclose(fp2);
  }
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr->barrier();
#endif
}


//-----------------------------------------------------------------------------
// Function      : Sensitivity::icSensitivity
// Purpose       : This function is called when there is a NOOP or UIC.  The
//                 dqdp vector needs to be set up for the history to be correct.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Sensitivity::icSensitivity ( 
     std::vector<double> & objectiveVec,
     std::vector<double> & dOdpVec, 
     std::vector<double> & dOdpAdjVec,
     std::vector<double> & scaled_dOdpVec, 
     std::vector<double> & scaled_dOdpAdjVec)
{
  if (!solveDirectFlag_ && !solveAdjointFlag_) return 1;

  N_TIA_DataStore & ds = *dsPtr_;

  ds.dOdpVec_.clear();
  ds.dOdpAdjVec_.clear();

  ds.scaled_dOdpVec_.clear();
  ds.scaled_dOdpAdjVec_.clear();

  // first get the derivatives of the RHS vector w.r.t. the
  // user-specified optimization parameters.
  loadSensitivityResiduals ();

  calcObjFuncDerivs ();

  objectiveVec.clear();
  objectiveVec.push_back(objFuncEval_);
  ds.objectiveVec_.clear();
  ds.objectiveVec_.push_back(objFuncEval_);

  if (solveDirectFlag_) 
  {
    ds.dOdpVec_.resize(numSensParams_,0.0);
    ds.scaled_dOdpVec_.resize(numSensParams_,0.0);

    if (outputUnscaledFlag_)
    {
      dOdpVec = ds.dOdpVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpVec = ds.scaled_dOdpVec_;
    }

    if (stdOutputFlag_)
    {
      stdOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_);
    }
  }

  if (solveAdjointFlag_) 
  {
    ds.dOdpAdjVec_.resize(numSensParams_,0.0);
    ds.scaled_dOdpAdjVec_.resize(numSensParams_,0.0);

    if (outputUnscaledFlag_)
    {
      dOdpAdjVec = ds.dOdpAdjVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpAdjVec = ds.scaled_dOdpAdjVec_;
    }
    if (stdOutputFlag_)
    {
      stdOutput(std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::solve
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/21/02
//-----------------------------------------------------------------------------
int Sensitivity::solve ( 
     std::vector<double> & objectiveVec,
     std::vector<double> & dOdpVec, 
     std::vector<double> & dOdpAdjVec,
     std::vector<double> & scaled_dOdpVec, 
     std::vector<double> & scaled_dOdpAdjVec)
{
  if (!solveDirectFlag_ && !solveAdjointFlag_) return 1;

  N_TIA_DataStore & ds = *dsPtr_;

  ds.dOdpVec_.clear();
  ds.dOdpAdjVec_.clear();

  ds.scaled_dOdpVec_.clear();
  ds.scaled_dOdpAdjVec_.clear();

  // It may now be neccessary to re-load the jacobian and rhs vectors.
  // It is necccessary, for example, if the two-level Newton solver is the
  // solver being used.
  nls_.enableSensitivity ();

  // first get the derivatives of the RHS vector w.r.t. the
  // user-specified optimization parameters.
  loadSensitivityResiduals ();

  calcObjFuncDerivs ();

  objectiveVec.clear();
  objectiveVec.push_back(objFuncEval_);
  ds.objectiveVec_.clear();
  ds.objectiveVec_.push_back(objFuncEval_);

  if (solveDirectFlag_) 
  {
    solveDirect ();

    if (outputUnscaledFlag_)
    {
      dOdpVec = ds.dOdpVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpVec = ds.scaled_dOdpVec_;
    }
  }

  if (solveAdjointFlag_) 
  {
    solveAdjoint ();
    if (outputUnscaledFlag_)
    {
      dOdpAdjVec = ds.dOdpAdjVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpAdjVec = ds.scaled_dOdpAdjVec_;
    }
  }

  numSolves_++;

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::solveDirect
//
// Purpose       : This function calculates the direct sensitivities for
//                 the user specified parameters.
//
//                 The ultimate goal of this function is to obtain dO/dp,
//                 where O is the objective function, and p is a
//                 user-defined optimization parameter.
//
//                 This is a bit confusing because the DAKOTA folks use
//                 a different naming convention than I have tended to
//                 use.  In Dakota's documentation, dO/dp would be referred
//                 to as df/du.  In Xyce, f is already considered to be the
//                 right-hand-side residual vector.  For clarity, this is
//                 the key between my notation and DAKOTA's:
//                        DAK     Xyce
//                         f       O
//                         y       x
//                         u       p
//                         c       f
//
//                 To obtain dOdp, the device manager is first called, and
//                 told to calculate dfdp, which is the derivative of the
//                 residual vector w.r.t. user-defined params.  It does
//                 this by numerical differentiation.  Then, after that,
//                 this function solves the equation:
//
//                 J dx/dp = df/dp    ->   dx/dp = J^-1 df/dp
//
//                 for each p.  
//
//                 J is the jacobian matrix, so J=df/dx.  (dc/dy)
//
//                 After performing these linear solves, dO/dp is to be
//                 obtained using the chain rule by:
//
//                 dO/dp = - dO/dx * dx/dp  + dO/dp
//
//                 The O in the dO/dp on the left side of this equation
//                 should have a hat over it "^", to indicate that it is
//                 different than the O on the right hand side.
//
//                 Note, this method is best if you have lots of objective
//                 functions, and a small number of parameters.  For adjoint
//                 calculations it is the other way around.
//
//                 11/19/02.
//                 It is assumed (for now) that dO/dp on the right hand side
//                 is zero, i.e., there is no direct
//                 dependence on p by O.  So, for example, if a user
//                 defined parameter to be used is the length of a MOSFET,
//                 the MOSFET length will NOT appear in the analytical
//                 expression for O.
//
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
int Sensitivity::solveDirect ()
{
#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "In Sensitivity::solveDirect" << std::endl;
  }
#endif

  N_TIA_DataStore & ds = *dsPtr_;

  int iparam;
  std::vector<N_LAS_Vector *> & dfdpPtrVector_ = ds.nextDfdpPtrVector;
  std::vector<N_LAS_Vector *> & dXdpPtrVector_ = ds.nextDXdpPtrVector;
  std::vector<N_LAS_Vector *> & currdXdpPtrVector_ = ds.currDXdpPtrVector;
  std::vector<N_LAS_Vector *> & sensRHSPtrVector = ds.sensRHSPtrVector;

  std::vector<N_LAS_Vector *> & DQdxDXdpPtrVector_ = ds.nextDQdxDXdpPtrVector;

  // first save a copy of the rhs and newton vectors, in case we want them later.
  savedRHSVectorPtr_->update(1.0, *(rhsVectorPtr_), 0.0);
  savedNewtonVectorPtr_->update(1.0, *(NewtonVectorPtr_), 0.0);

  // via the loader, setup all the sensitivity residuals.
  loaderPtr_->loadSensitivityResiduals();

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::string matrixFile = netlistFileName_ + "_directMatrix.txt";
    jacobianMatrixPtr_->writeToFile(const_cast<char *>(matrixFile.c_str()));
  }
#endif

  // Now solve the series of linear systems to get dXdp.
  for (iparam=0; iparam< numSensParams_; ++iparam)
  {
    // copy the sensitivity residual into the rhs vector location,
    // as that is the RHS vector that the linear solver expects.
    rhsVectorPtr_->update(1.0, *(sensRHSPtrVector[iparam]), 0.0);

    lasSolverPtr_->solve();

    // copy the result of the linear solve into the dxdp data structure.
    (dXdpPtrVector_[iparam])->update(1.0, *(NewtonVectorPtr_), 0.0);

#ifdef Xyce_DEBUG_NONLINEAR
    // do debug output.
    if (debugLevel_ > 0)
    {
      Xyce::dout() << "iparam="<<iparam << "\t" << paramNameVec_[iparam] <<std::endl;
      for (int k = 0; k < solutionSize_; ++k)
      {
        Xyce::dout() << "dXdp[" << std::setw(3) << k << "] = "<< std::setw(15)<< std::scientific
          << std::setprecision(8)<< (*(dXdpPtrVector_[iparam]))[k]
          <<std::endl;
      }

      std::ostringstream filename; 
      filename << netlistFileName_ << "_dxdp";
      filename << std::setw(3) << std::setfill('0') << iparam;
      filename << ".txt";
      dXdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));
    }
#endif

    // Now store the DQdx*dXdp matvec
    N_LAS_Matrix & dQdx = *(ds.dQdxMatrixPtr);
    bool Transpose = false;
    (DQdxDXdpPtrVector_[iparam])->putScalar(0.0);
    dQdx.matvec( Transpose , *(dXdpPtrVector_[iparam]), *(DQdxDXdpPtrVector_[iparam]) );

  }// end of param for loop

  // Now update the time derivative of the matvec (this is needed for trap).
  loaderPtr_->loadFinalSensitivityDerivatives ();

  // Restore the RHS and Newton vectors.
  rhsVectorPtr_->update(1.0, *(savedRHSVectorPtr_),0.0);
  NewtonVectorPtr_->update(1.0, *(savedNewtonVectorPtr_),0.0);

  // Now get the final dOdp's (one for each param).
  for (iparam=0; iparam< numSensParams_; ++iparam)
  {
    double tmp = dOdXVectorPtr_->dotProduct( (*(dXdpPtrVector_[iparam]))  );
    tmp += dOdp_;

    ds.dOdpVec_.push_back(tmp);
 
    // get scaled value.  dO/dp*(p/100)
    double normalize = ds.paramOrigVals_[iparam]/100.0;
    tmp *= normalize;
    ds.scaled_dOdpVec_.push_back(tmp);
  }

  if (stdOutputFlag_)
  {
    stdOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_);
  }

  if (fileOutputFlag_)
  {
    fileOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_);
  }

  if (dakotaFileOutputFlag_)
  {
    dakOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_);
  }

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::calcObjFuncDerivs
//
// Purpose       : This function assumes an objective function, O, and
//                 uses it to calculate other stuff, including dO'/du,
//                 dO/dy, dO/du.
//
//                 The objective function simply assumes that you are
//                 trying to match a solution variable with a
//                 user-specified number.   (objective=1)
//
//                 Or that you are trying to match your final solution
//                 with a specified solution. (objective=2)
//
// Special Notes : This is mostly for testing purposes.  Usually, the
//                 objective function would be set and managed from
//                 the DAKOTA side.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/02
//-----------------------------------------------------------------------------
bool Sensitivity::calcObjFuncDerivs ()
{
  bool bsuccess = true;
  int i;

  dOdXVectorPtr_->putScalar(0.0);

  // first check if the user specified a function via the expressions
  // class:
  int found(0);
  int found2(0);
  bool foundLocal(false);
  bool foundLocal2(false);

  N_PDS_Comm & comm = *(pdsMgrPtr_->getPDSComm());
  int myPID = comm.procID();

  if (objFuncGiven_)
  {
    if (!objFuncGIDsetup_)
    {
      // set up the gid's:
      expVarGIDs_.resize( expNumVars_, -1 );
      for (i = 0; i < expNumVars_; ++i)
      {
        std::list<int> svGIDList1, dummyList;
        char type1;
        foundLocal = top_.getNodeSVarGIDs(NodeID(expVarNames_[i], Xyce::_VNODE), svGIDList1, dummyList, type1);
        found = static_cast<int>(foundLocal);
        Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found, 1);

        foundLocal2 = false;
        if (!found)// if looking for this as a voltage node failed, try a "device" (i.e. current) node.
        {
          foundLocal2 = top_.getNodeSVarGIDs(NodeID(expVarNames_[i], Xyce::_DNODE), svGIDList1, dummyList, type1);
        }
        found2 = static_cast<int>(foundLocal2);
        Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found2, 1);

        if (!found && !found2)
        {
          static std::string tmp = "objective function variable not found!\n";
          N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, tmp);
        }

        if (found || found2)
        {
          int tmpGID=-1;
          if(svGIDList1.size()==1)
          {
            tmpGID = svGIDList1.front();
          }
          expVarGIDs_[i] = tmpGID;
        }
      }
      objFuncGIDsetup_ = true;
    }

    // obtain the expression variable values.  It will only grab this value if it is owned on this processor.
    expVarVals_.resize (expNumVars_, 0.0);
    expVarDerivs_.resize (expNumVars_, 0.0);

    comm.barrier();
    for (i = 0; i < expNumVars_; ++i)
    {
      int tmpGID=expVarGIDs_[i];
      double tmpVal=0.0;
      int root=-1;
      if (tmpGID >= 0)
      {
        tmpVal = (*nextSolVectorPtrPtr_)->getElementByGlobalIndex(tmpGID, 0);
        root = myPID;
      }
      Xyce::Parallel::AllReduce(comm.comm(), MPI_MAX, &root, 1);

      comm.bcast( &tmpVal, 1, root );

      expVarVals_[i] = tmpVal;
    }
    comm.barrier();
    
    //get expression value and partial derivatives
    expPtr_->evaluate( expVal_, expVarDerivs_, expVarVals_ );
    objFuncEval_ = expVal_;
    dOdXVectorPtr_->putScalar(0.0);
    for (i=0;i<expNumVars_;++i)
    {
      int tmpGID = expVarGIDs_[i];
      double tmpDODX = expVarDerivs_[i];

#ifdef Xyce_DEBUG_NONLINEAR
      if (debugLevel_ > 0)
      {
        Xyce::dout() << "i="<<i<<"  gid = " << tmpGID << "  dodx = "<< tmpDODX << std::endl;
      }
#endif

      if (tmpGID >= 0)
      {
        dOdXVectorPtr_->setElementByGlobalIndex(tmpGID, tmpDODX, 0);
      }
    }

    // Assuming this is zero:
    dOdp_ = 0.0;
  }
  else // give a warning.
  {
    static std::string tmp = " ***** WARNING: objective function was not specified.\n\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, tmp);
  }

  dOdXVectorPtr_->fillComplete();

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::string filename = netlistFileName_ + "_dodx.txt";
    dOdXVectorPtr_->writeToFile(const_cast<char *>(filename.c_str()));
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::solveAdjoint
// Purpose       : Solves first for the vector, lambda, of adjoint variables.
//                 Afterwards, it solves for dO/dp.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
int Sensitivity::solveAdjoint ()
{
#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "In Sensitivity::solveAdjoint" << std::endl;
  }
#endif

  N_TIA_DataStore & ds = *dsPtr_;

  std::vector<N_LAS_Vector *> & dfdpPtrVector_ = ds.nextDfdpPtrVector;
  std::vector<N_LAS_Vector *> & dXdpPtrVector_ = ds.nextDXdpPtrVector;

  // first save a copy of the rhs vector, in case we want it later.
  savedRHSVectorPtr_->update(1.0, *(rhsVectorPtr_),0.0);
  savedNewtonVectorPtr_->update(1.0, *(NewtonVectorPtr_),0.0);

  lambdaVectorPtr_ = lasSysPtr_->builder().createVector();
  bool useTran = jacobianMatrixPtr_->useTranspose ();
  jacobianMatrixPtr_->setUseTranspose (true);

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::string matrixFile = netlistFileName_ + "_adjointMatrix.txt";
    jacobianMatrixPtr_->writeToFile(const_cast<char *>(matrixFile.c_str()));
  }
#endif

  // copy the current dOdx vector into the rhs vector data structure.
  rhsVectorPtr_->update(1.0, *(dOdXVectorPtr_),0.0);

  int status = lasSolverPtr_->solve();
  if (status!=0)
  {
    std::string msg("Sensitivity::solveAdjoint.  Solver failed\n");
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // allocate the dxdp vector for this param, and
  // copy the resulting deltax vector into the dxdp data structure.
  lambdaVectorPtr_->update(1.0, *(NewtonVectorPtr_),0.0);

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::string filename = netlistFileName_ + "_lambda.txt";
    lambdaVectorPtr_->writeToFile(const_cast<char *>(filename.c_str()));
  }
#endif

  // Now that we have lambda, get the dOdp's by doing dot products of
  // lambda * df/dp.

  // do the final dot products, one for each param.
  for (int iparam=0; iparam< numSensParams_; ++iparam)
  {
    double tmp = -1.0 * lambdaVectorPtr_->dotProduct(*(dfdpPtrVector_[iparam]));
    ds.dOdpAdjVec_.push_back(tmp);

    // get scaled value.  dO/dp*(p/100)
    double normalize = ds.paramOrigVals_[iparam]/100.0;
    tmp *= normalize;
    ds.scaled_dOdpAdjVec_.push_back(tmp);
  }


  // restore the useTranspose flag to the original setting. (probably
  // false)
  jacobianMatrixPtr_->setUseTranspose (useTran);

  // Restore the RHS and Newton vectors.
  rhsVectorPtr_->update(1.0, *(savedRHSVectorPtr_),0.0);
  NewtonVectorPtr_->update(1.0, *(savedNewtonVectorPtr_),0.0);

  // set the sensitivity information to the screen:
  if (stdOutputFlag_)
  {
    stdOutput(std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_);
  }
  if (fileOutputFlag_)
  {
    fileOutput(std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_);
  }

  if (dakotaFileOutputFlag_)
  {
    dakOutput(std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_);
  }

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::loadSensitivityResiduals
//
// Purpose       : This function is to be called after a calculation has
//                 converged.  It computes the sensitivies of different components
//                 of the residual.  Recall the DAE form of the residual:
//                
//                 F = dq/dt + f - b
//
//                 This function sets up the following derivatives:  
//
//                 df/dp, dq/dp, db/dp where p is the parameter.
//
//                 Which can ultimately be assembled to give dF/dp:
//
//                         d(dq/dp)    
//                 dF/dp = ------- +  df/dp - db/dp
//                           dt
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
bool Sensitivity::loadSensitivityResiduals ()
{
  int iparam;
  std::string msg;
  N_LAS_System & lasSys_ = (*lasSysPtr_);

  N_TIA_DataStore & ds = *dsPtr_;

  std::vector<N_LAS_Vector *> & dfdpPtrVector_ = ds.nextDfdpPtrVector;
  std::vector<N_LAS_Vector *> & dqdpPtrVector_ = ds.nextDqdpPtrVector;
  std::vector<N_LAS_Vector *> & dbdpPtrVector_ = ds.nextDbdpPtrVector;
  std::vector<N_LAS_Vector *> & dXdpPtrVector_ = ds.nextDXdpPtrVector;

  ds.paramOrigVals_.clear();

  // it is necessary to load the Jacobian here to make sure we have the most
  // up-to-date matrix.  The Jacobian is not loaded for the final 
  // evaluation of the residual in the Newton solve.
  loaderPtr_->loadJacobian ();

  // Loop over the vector of parameters.  For each parameter, find the
  // device entity (a model or an instance) which corresponds to it, and
  // perform the finite difference calculation.
  std::vector<std::string>::iterator firstParam = paramNameVec_.begin ();
  std::vector<std::string>::iterator lastParam  = paramNameVec_.end ();
  std::vector<std::string>::iterator iterParam;
  for ( iterParam=firstParam, iparam=0;
        iterParam!=lastParam; ++iterParam, ++iparam )
  {
    std::string paramName(*iterParam);

    // get the original value of this parameter, to be used for scaling and/or 
    // numerical derivatives later.
    double paramOrig = 0.0;
    bool found = loaderPtr_->getParamAndReduce(paramName, paramOrig);
    if (!found)
    {
      std::string msg("Sensitivity::loadSensitivityResiduals: cannot find parameter ");
      msg += paramName;
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::DEV_FATAL, msg);
    }
    ds.paramOrigVals_.push_back(paramOrig);

    // check if derivative is available analytically.  If not, take FD.
    bool analyticAvailable = false;
    if (!forceFD_)
    {
      analyticAvailable = loaderPtr_->analyticSensitivitiesAvailable (paramName);
    }
    if (analyticAvailable)
    {
      std::vector<double> dfdpVec;
      std::vector<double> dqdpVec;
      std::vector<double> dbdpVec;

      std::vector<int> FindicesVec;
      std::vector<int> QindicesVec;
      std::vector<int> BindicesVec;

      loaderPtr_->getAnalyticSensitivities(paramName,
          dfdpVec,dqdpVec,dbdpVec,
          FindicesVec, QindicesVec, BindicesVec);

      dfdpPtrVector_[iparam]->putScalar(0.0);
      dqdpPtrVector_[iparam]->putScalar(0.0);
      dbdpPtrVector_[iparam]->putScalar(0.0);

      int Fsize=FindicesVec.size();
      for (int i=0;i<Fsize;++i)
      {
        N_LAS_Vector & dfdpRef = *(dfdpPtrVector_[iparam]);
        dfdpRef[FindicesVec[i]]  = dfdpVec[i];
      }

      int Qsize=QindicesVec.size();
      for (int i=0;i<Qsize;++i)
      {
        N_LAS_Vector & dqdpRef = *(dqdpPtrVector_[iparam]);
        dqdpRef[QindicesVec[i]]  = dqdpVec[i];
      }

      int Bsize=BindicesVec.size();
      for (int i=0;i<Bsize;++i)
      {
        N_LAS_Vector & dbdpRef = *(dbdpPtrVector_[iparam]);
        dbdpRef[BindicesVec[i]]  = dbdpVec[i];
      }

      dfdpPtrVector_[iparam]->fillComplete();
      dqdpPtrVector_[iparam]->fillComplete();
      dbdpPtrVector_[iparam]->fillComplete();

#ifdef Xyce_DEBUG_NONLINEAR
      if (debugLevel_ > 0)
      {
        Xyce::dout() << *iterParam << ": ";
        Xyce::dout().setf(std::ios::scientific);
        Xyce::dout() << std::endl;

        for (int k1 = 0; k1 < solutionSize_; ++k1)
        {
          Xyce::dout() 
            <<"dfdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(dfdpPtrVector_[iparam]))[k1]
            <<" dqdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(dqdpPtrVector_[iparam]))[k1]
            <<" dbdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(dbdpPtrVector_[iparam]))[k1]
            <<std::endl;
        }
      }
#endif
    }
    else
    {
#ifdef Xyce_DEBUG_NONLINEAR
      if (debugLevel_ > 0)
      {
        Xyce::dout() << std::endl << "  Calculating numerical df/dp, dq/dp and db/dp for: ";
        Xyce::dout() << *iterParam << std::endl;
      }
#endif

      // save a copy of the DAE vectors
      origFVectorPtr_->update(1.0, *(ds.daeFVectorPtr), 0.0);
      origQVectorPtr_->update(1.0, *(ds.daeQVectorPtr), 0.0);
      origBVectorPtr_->update(1.0, *(ds.daeBVectorPtr), 0.0);

      // now perturb the value of this parameter.
      double dp = sqrtEta_ * (1.0 + fabs(paramOrig));
      double paramPerturbed = paramOrig;

      if (difference==SENS_FWD)
      {
        paramPerturbed += dp;
      }
      else if (difference==SENS_REV)
      {
        paramPerturbed -= dp;
      }
      else if (difference==SENS_CNT)
      {
        static std::string tmp = "difference=central not supported.\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
      }
      else
      {
        static std::string tmp = "difference not recognized!\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
      }

#ifdef Xyce_DEBUG_NONLINEAR
      if (debugLevel_ > 0)
      {
        Xyce::dout() << std::setw(maxParamStringSize_)<< *iterParam
          << " dp = " << std::setw(11)<< std::scientific<< std::setprecision(4) << dp 
          << " original value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramOrig 
          << " modified value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramPerturbed 
          <<std::endl;
      }
#endif
      loaderPtr_->setParam (paramName, paramPerturbed);

      // Now that the parameter has been perturbed,
      // calculate the numerical derivative.

      // Load F,Q and B.
      loaderPtr_->loadRHS();

      // save the perturbed DAE vectors
      pertFVectorPtr_->update(1.0, *(ds.daeFVectorPtr), 0.0);
      pertQVectorPtr_->update(1.0, *(ds.daeQVectorPtr), 0.0);
      pertBVectorPtr_->update(1.0, *(ds.daeBVectorPtr), 0.0);

      // calculate the df/dp vector.  
      double rdp=1/dp;
      dfdpPtrVector_[iparam]->putScalar(0.0);
      dfdpPtrVector_[iparam]->addVec (+1.0, *(pertFVectorPtr_)); //+Fperturb
      dfdpPtrVector_[iparam]->addVec (-1.0, *(origFVectorPtr_)); //-Forig
      dfdpPtrVector_[iparam]->scale(rdp);

      // calculate the dq/dp vector.  
      dqdpPtrVector_[iparam]->putScalar(0.0);
      dqdpPtrVector_[iparam]->addVec (+1.0, *(pertQVectorPtr_)); //+Fperturb
      dqdpPtrVector_[iparam]->addVec (-1.0, *(origQVectorPtr_)); //-Forig
      dqdpPtrVector_[iparam]->scale(rdp);

      // calculate the db/dp vector.  
      dbdpPtrVector_[iparam]->putScalar(0.0);
      dbdpPtrVector_[iparam]->addVec (+1.0, *(pertBVectorPtr_)); //+Fperturb
      dbdpPtrVector_[iparam]->addVec (-1.0, *(origBVectorPtr_)); //-Forig
      dbdpPtrVector_[iparam]->scale(rdp);

#ifdef Xyce_DEBUG_NONLINEAR
      if (debugLevel_ > 0)
      {
        Xyce::dout() << *iterParam << ": ";
        Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
        Xyce::dout() << "deviceSens_dp = " << dp << std::endl;

        for (int k1 = 0; k1 < solutionSize_; ++k1)
        {

          Xyce::dout() 
            <<"fpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertFVectorPtr_))[k1]
            <<" forig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origFVectorPtr_))[k1]
            <<" dfdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(dfdpPtrVector_[iparam]))[k1]
            <<std::endl;
        }

        Xyce::dout() << std::endl;
        for (int k1 = 0; k1 < solutionSize_; ++k1)
        {
          Xyce::dout() 
            <<"qpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertQVectorPtr_))[k1]
            <<" qorig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origQVectorPtr_))[k1]
            <<" dqdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(dqdpPtrVector_[iparam]))[k1]
            <<std::endl;
        }

        Xyce::dout() << std::endl ;
        for (int k1 = 0; k1 < solutionSize_; ++k1)
        {
          Xyce::dout() 
            <<"bpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertBVectorPtr_))[k1]
            <<" borig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origBVectorPtr_))[k1]
            <<" dbdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(dbdpPtrVector_[iparam]))[k1]
            <<std::endl;

        }

        std::ostringstream filename; 
        filename << netlistFileName_ << "_dfdp";
        filename << std::setw(3) << std::setfill('0') << iparam;
        filename << ".txt";
        dfdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

        filename.str("");
        filename << netlistFileName_ << "_fpert";
        filename << std::setw(3) << std::setfill('0') << iparam;
        filename << ".txt";
        pertFVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));

        filename.str("");
        filename << netlistFileName_ << "_dqdp";
        filename << std::setw(3) << std::setfill('0') << iparam;
        filename << ".txt";
        dqdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

        filename.str("");
        filename << netlistFileName_ << "_qpert";
        filename << std::setw(3) << std::setfill('0') << iparam;
        filename << ".txt";
        pertQVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));

        filename.str("");
        filename << netlistFileName_ << "_dbdp";
        filename << std::setw(3) << std::setfill('0') << iparam;
        filename << ".txt";
        dbdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

        filename.str("");
        filename << netlistFileName_ << "_bpert";
        filename << std::setw(3) << std::setfill('0') << iparam;
        filename << ".txt";
        pertBVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));
      }
#endif

      // now reset the parameter and rhs to previous values.
      loaderPtr_->setParam (paramName, paramOrig);
      
      //rhsVectorPtr_->update(-1.0, *(origFVectorPtr_), 0.0);
      ds.daeFVectorPtr->update(1.0, *(origFVectorPtr_), 0.0);
      ds.daeQVectorPtr->update(1.0, *(origQVectorPtr_), 0.0);
      ds.daeBVectorPtr->update(1.0, *(origBVectorPtr_), 0.0);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setOptions
//
// Purpose       : This function processes the .SENS line
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/02
//-----------------------------------------------------------------------------
bool Sensitivity::setOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;
  std::list<N_UTL_Param>::const_iterator iter = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator end   = OB.getParams().end();

  numSensParams_ = 0;
  for ( ; iter != end; ++ iter)
  {
    if (iter->uTag() == "OBJFUNC")
    {
      objFuncString_ = iter->stringValue();
      expPtr_ = new N_UTL_Expression(iter->stringValue());
      objFuncGiven_ = true;
    }
    else if ( std::string( iter->uTag() ,0,5) == "PARAM") // this is a vector
    {
      ExtendedString tag = iter->stringValue();
      tag.toUpper();
      // set up the initial skeleton of the maps:
      ++numSensParams_;
      paramNameVec_.push_back(tag);
      int sz = tag.size();
      if (sz > maxParamStringSize_)
      {
        maxParamStringSize_ = sz;
      }
    }
    else
    {
      Xyce::Report::UserWarning() << iter->uTag() 
        << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  // parse the expression now, so if there are any errors, they will come
  // up early in the simulation.
  if (objFuncGiven_)
  {
    // setup the names:
    expVarNames_.clear();

    std::vector<std::string> nodes;
    expPtr_->get_names(XEXP_NODE, nodes);
    std::vector<std::string> instances;
    expPtr_->get_names(XEXP_INSTANCE, instances);

    // Make the current (instance) strings all upper case.
    // The topology directory apparently requires this.
    std::vector<std::string>::iterator iter;
    for (iter=instances.begin();iter!=instances.end();++iter)
    {
      ExtendedString tmpString = *iter;
      tmpString.toUpper ();
      *iter  = tmpString;
    }

    expVarNames_.insert(expVarNames_.end(), nodes.begin(), nodes.end());
    expVarNames_.insert(expVarNames_.end(), instances.begin(), instances.end());

    // Order the names in the expression so that it agrees with the order
    // in expVarNames.
    if ( !expVarNames_.empty() )
    {
      expPtr_->order_names( expVarNames_ );
    }

    expNumVars_ = expVarNames_.size();
  }

#ifdef Xyce_DEBUG_NONLINEAR
  if (debugLevel_ > 0)
  {
    std::vector<std::string>::iterator iter;
    std::vector<std::string>::iterator begin = paramNameVec_.begin();
    std::vector<std::string>::iterator end   = paramNameVec_.end  ();

    for (iter=begin;iter!=end;++iter)
    {
      Xyce::dout() << *iter<<std::endl;
    }
  }
#endif

  dsPtr_->allocateSensitivityArrays(numSensParams_);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setSensitivityOptions
// Purpose       : This function processes the .options SENSITIVITY line
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/20/2013
//-----------------------------------------------------------------------------
bool Sensitivity::setSensitivityOptions(const N_UTL_OptionBlock &OB)
{
  bool bsuccess = true;
  std::list<N_UTL_Param>::const_iterator it  = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator end = OB.getParams().end();
  for ( ; it != end; ++ it)
  {
    if ((*it).uTag() == "ADJOINT")
    {
      solveAdjointFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DIRECT")
    {
      solveDirectFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUTSCALED")
    {
      outputScaledFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUTUNSCALED")
    {
      outputUnscaledFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "STDOUTPUT")
    {
      stdOutputFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DAKOTAFILE")
    {
      dakotaFileOutputFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DIAGNOSTICFILE")
    {
      fileOutputFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "FORCEFD")
    {
      forceFD_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DIFFERENCE")
    {
      ExtendedString sval=(*it).stringValue();
      sval.toUpper();
      if(sval=="FORWARD")
      {
        difference=SENS_FWD;
      }
      else if(sval=="REVERSE")
      {
        difference=SENS_REV;
      }
      else if(sval=="CENTRAL")
      {
        difference=SENS_CNT;
        static std::string tmp = "difference=central not supported.\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
      }
      else
      {
        static std::string tmp = "difference not recognized!\n";
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, tmp);
      }
    }
    else if ((*it).uTag() == "SQRTETA")
    {
      sqrtEta_ = (*it).getImmutableValue<double>();
      sqrtEtaGiven_ = true;
    }
#ifdef Xyce_DEBUG_NONLINEAR
    else if ((*it).uTag() == "DEBUGLEVEL")
    {
      debugLevel_ = (*it).getImmutableValue<int>();
    }
#endif
    else
    {
      Xyce::Report::UserWarning() << (*it).uTag() 
        << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setTranOptions
//
// Purpose       : Not used yet.  Same as setOptions, but for transient
//                 mode, if and when the sensitivity stuff is ever set
//                 up to work in transient.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
bool Sensitivity::setTranOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool Sensitivity::setHBOptions(const N_UTL_OptionBlock& OB)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::getMaxNormF
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and MEMS Modeling
// Creation Date : 9/28/2009
//-----------------------------------------------------------------------------
double Sensitivity::getMaxNormF() const
{
  return nls_.getMaxNormF();
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::getMaxNormFindex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and MEMS Modeling
// Creation Date : 9/28/2009
//-----------------------------------------------------------------------------
int Sensitivity::getMaxNormFindex () const
{
  return nls_.getMaxNormFindex ();
}

} // namespace Nonlinear
} // namespace Xyce
