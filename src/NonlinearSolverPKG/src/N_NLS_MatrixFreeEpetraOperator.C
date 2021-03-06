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
// Filename       : $RCSfile: N_NLS_MatrixFreeEpetraOperator.C,v $
//
// Purpose        :
//
// Creator        : Todd Coffey, 1414
//
// Creation Date  : 09/04/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9 $
//
// Revision Date  : $Date: 2014/08/07 23:08:54 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_MatrixFreeEpetraOperator.h>
#include <N_LAS_Vector.h>
#include <N_PDS_ParMap.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : matrixFreeEpetraOperator
// Purpose       : non-member constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
RefCountPtr<MatrixFreeEpetraOperator> matrixFreeEpetraOperator(
    RefCountPtr<NonLinearSolver> nonlinearSolver,
    RefCountPtr<N_LAS_Vector> solVector,
    RefCountPtr<N_LAS_Vector> rhsVector,
    RefCountPtr<const Epetra_Map> solutionMap
    )
{
  RefCountPtr<MatrixFreeEpetraOperator> epetraOperator =
    rcp(new MatrixFreeEpetraOperator);
  epetraOperator->initialize(nonlinearSolver,
      solVector,
      rhsVector,
      solutionMap
      );
  return epetraOperator;
}


//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::MatrixFreeEpetraOperator
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
MatrixFreeEpetraOperator::MatrixFreeEpetraOperator()
{
  isInitialized_ = false;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::MatrixFreeEpetraOperator
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
MatrixFreeEpetraOperator::~MatrixFreeEpetraOperator()
{
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::initialize
// Purpose       : Initialization
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
void MatrixFreeEpetraOperator::initialize(
      RefCountPtr<NonLinearSolver> nonlinearSolver,
      RefCountPtr<N_LAS_Vector> solVector,
      RefCountPtr<N_LAS_Vector> rhsVector,
      RefCountPtr<const Epetra_Map> solutionMap
    )
{
  nonlinearSolverRCPtr_ = nonlinearSolver;
  solVectorRCPtr_ = solVector;
  rhsVectorRCPtr_ = rhsVector;
  solutionMap_ = solutionMap;
  isInitialized_ = true;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::SetUseTranspose
// Purpose       : Define if transpose Apply and ApplyInverse is to be used.
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::SetUseTranspose(bool UseTranspose)
{
  // This is not supported for the HB load layers.
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Apply
// Purpose       : Apply matrix free operator with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::Apply(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  // Convert these to N_LAS_MultiVectors and call the other Apply

  // Cast away the const until the Apply which will enforce it.
  // This is necessary because there is no const view function in N_LAS_MultiVector
//  Epetra_MultiVector* Xptr = const_cast<Epetra_MultiVector*>(&X);
//  N_LAS_MultiVector las_X(Xptr);  // This is the wrong thing to do, when it goes out of scope, it deletes the Epetra_MultiVector Ptr.
//  N_LAS_MultiVector las_Y(&Y);

  // COPY the multi-vector data into new objects on the stack.
  Epetra_MultiVector* Xcopy = new Epetra_MultiVector(X); // This gets deleted by the N_LAS_MultiVector below
  Epetra_MultiVector* Ycopy = new Epetra_MultiVector(Y); // This gets deleted by the N_LAS_MultiVector below
  N_LAS_MultiVector las_X(Xcopy); // this co-ops the Epetra_MultiVector and uses (and owns) its memory
  N_LAS_MultiVector las_Y(Ycopy); // this co-ops the Epetra_MultiVector and uses (and owns) its memory
  int status = Apply(las_X,las_Y);
  // COPY the Ycopy data back into Y
  Y = las_Y.epetraObj();
  return(status);
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Apply
// Purpose       : Apply matrix free operator with N_LAS_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::Apply(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  if (!isInitialized_)
  {
    std::string msg = "MatrixFreeEpetraOperator::Apply:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  bool status = true;
  for (int i=0 ; i<X.numVectors() ; ++i)
  {
    const N_LAS_Vector x(X.epetraVector(i));
    N_LAS_Vector y(Y.epetraVector(i));
    bool localStatus = nonlinearSolverRCPtr_->applyJacobian(x,y);
    status = status && localStatus;
  }
  if (status)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}
//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::ApplyInverse
// Purpose       : Apply inverse of matrix free operator with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::ApplyInverse(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  std::string msg = "MatrixFreeEpetraOperator::ApplyInverse is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::ApplyInverse
// Purpose       : Apply inverse of matrix free operator with N_LAS_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int MatrixFreeEpetraOperator::ApplyInverse(
  const N_LAS_MultiVector& X,
  N_LAS_MultiVector& Y
  ) const
{
  std::string msg = "MatrixFreeEpetraOperator::ApplyInverse is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::NormInf
// Purpose       : Norm Inf of matrix
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
double MatrixFreeEpetraOperator::NormInf() const
{
  std::string msg = "MatrixFreeEpetraOperator::NormInf is not supported!";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  return -1.0;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Label
// Purpose       : Label for operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const char * MatrixFreeEpetraOperator::Label() const
{
  return "Matrix Free Harmonic Balance Epetra Operator";
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::UseTranspose
// Purpose       : Query for useTranspose setting
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
bool MatrixFreeEpetraOperator::UseTranspose() const
{
  // Use Transpose is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::HasNormInf
// Purpose       : Query for normInf support
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
bool MatrixFreeEpetraOperator::HasNormInf() const
{
  // Norm Inf is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::Comm
// Purpose       : Return Epetra_Comm object
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Comm & MatrixFreeEpetraOperator::Comm() const
{
  if (!isInitialized_)
  {
    std::string msg = "MatrixFreeEpetraOperator::Comm:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(rhsVectorRCPtr_->epetraObj().Comm());
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::OperatorDomainMap
// Purpose       : Return Epetra_Map corresponding to domain of operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Map & MatrixFreeEpetraOperator::OperatorDomainMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "MatrixFreeEpetraOperator::OperatorDomainMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  return(*solutionMap_);
}

//-----------------------------------------------------------------------------
// Function      : MatrixFreeEpetraOperator::OperatorRangeMap
// Purpose       : Return Epetra_Map corresponding to range of operator
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
const Epetra_Map & MatrixFreeEpetraOperator::OperatorRangeMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "MatrixFreeEpetraOperator::OperatorRangeMap:  I'm not initialized!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  const Epetra_Map* emap = dynamic_cast<const Epetra_Map*>(&rhsVectorRCPtr_->epetraObj().Map());
  return(*solutionMap_);
}

} // namespace Nonlinear
} // namespace Xyce
