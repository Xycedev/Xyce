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
// Filename       : $$
//
// Purpose        : Constraint Backtracking Class.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 01/26/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------

#include <N_NLS_ConstraintBT.h>

#include <N_LAS_Vector.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_NLS_NLParams.h>

#ifdef Xyce_DEBUG_NONLINEAR
#include <N_ERH_ErrorMgr.h>
#endif

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::ConstraintBT
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
ConstraintBT::ConstraintBT()
 : constraintMinVector_(0),
   constraintMaxVector_(0),
   constraintChangeVector_(0),
   constraintTempVector_(0)
{

  // Initialize protected data
  resetThetaBoundNeg();
  resetThetaBoundPos();
  resetThetaChange();

}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::ConstraintBT
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
ConstraintBT::ConstraintBT(const ConstraintBT & right)
{

}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::~ConstraintBT
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
ConstraintBT::~ConstraintBT()
{
  if (constraintMinVector_ != 0)    delete constraintMinVector_;
  if (constraintMaxVector_ != 0)    delete constraintMaxVector_;
  if (constraintChangeVector_ != 0) delete constraintChangeVector_;
  if (constraintTempVector_ != 0)   delete constraintTempVector_;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::operator=
// Purpose       : Assignment operator
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
ConstraintBT & ConstraintBT::operator=(
  const ConstraintBT & right)
{
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::operator==
// Purpose       : Equal operator
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
int ConstraintBT::operator==(const ConstraintBT & right) const

{
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::operator!=
// Purpose       : Not-Equal operator
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
int ConstraintBT::operator!=(const ConstraintBT & right) const

{
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::initializeAll
// Purpose       : Not-Equal operator
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/29/02
//-----------------------------------------------------------------------------
bool ConstraintBT::initializeAll
    (N_LAS_System * lasSysPtr, const NLParams & nlParams)
{
  // create and initialize constraint backtracking vectors:
  constraintMinVector_ = lasSysPtr->builder().createVector();
  constraintMaxVector_ = lasSysPtr->builder().createVector();

  constraintMinVector_->putScalar(nlParams.getGlobalBTMin());
  constraintMaxVector_->putScalar(nlParams.getGlobalBTMax());

  constraintChangeVector_ = lasSysPtr->builder().createVector();

  constraintChangeVector_->putScalar (nlParams.getGlobalBTChange());

  constraintTempVector_ = lasSysPtr->builder().createVector();

  constraintTempVector_->putScalar(0.0);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::updateThetaBoundNeg
// Purpose       : Updates the minimum bound value for the backtracking
//                 algorithm.
// Special Notes : This is implemented according to a internal
//                 communication with John Shadid (SNL) on their MPSalsa
//                 constraint backtracking implementation.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
void ConstraintBT::updateThetaBoundNeg(const N_LAS_Vector * oldSoln,
                                             const N_LAS_Vector * solnUpdate)

{
  N_LAS_Vector * solnPtr;
  N_LAS_Vector * solnUpdatePtr;

  // Initialize
  solnPtr       = const_cast<N_LAS_Vector *> (oldSoln);
  solnUpdatePtr = const_cast<N_LAS_Vector *> (solnUpdate);

  // First, form a vector of constraints...
  for (int i = 0; i < solnPtr->localLength(); ++i)
  {
    if ((*(solnUpdatePtr))[i] < 0.0)
      (*(constraintTempVector_))[i] = ((*(constraintMinVector_))[i] -
      (*(solnPtr))[i]) / (*(solnUpdatePtr))[i];
    else
      (*(constraintTempVector_))[i] = N_UTL_MachineDependentParams::DoubleMax();

    if (Xyce::DEBUG_NONLINEAR)
      Xyce::dout() << "ConstraintBT::updateThetaBoundNeg: min: " << (*(constraintMinVector_))[i] << "\n"
                   << "ConstraintBT::updateThetaBoundNeg: soln: " << (*(solnPtr))[i] << "\n"
                   << "ConstraintBT::updateThetaBoundNeg: solnUpdate: " << (*(solnUpdatePtr))[i] << "\n"
                   << "ConstraintBT::updateThetaBoundNeg: constraint: " << (*(constraintTempVector_))[i];
  }

  // Find minimum
  constraintTempVector_->minValue(&thetaBoundNeg_);

  if (Xyce::DEBUG_NONLINEAR)
    Xyce::dout() << "ConstraintBT::updateThetaBoundNeg: thetaBoundNeg_: " << thetaBoundNeg_;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::updateThetaBoundPos
// Purpose       : Updates the maximum bound value for the backtracking
//                 algorithm.
// Special Notes : This is implemented according to a internal
//                 communication with John Shadid (SNL) on their MPSalsa
//                 constraint backtracking implementation.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
void ConstraintBT::updateThetaBoundPos(const N_LAS_Vector * oldSoln,
                                             const N_LAS_Vector * solnUpdate)

{
  N_LAS_Vector * solnPtr;
  N_LAS_Vector * solnUpdatePtr;

  // Initialize
  solnPtr       = const_cast<N_LAS_Vector *> (oldSoln);
  solnUpdatePtr = const_cast<N_LAS_Vector *> (solnUpdate);

  // First, form a vector of constraints...
  for (int i = 0; i < solnPtr->localLength(); ++i)
  {
    if ((*(solnUpdatePtr))[i] > 0.0)
      (*(constraintTempVector_))[i] = ((*(constraintMaxVector_))[i] -
      (*(solnPtr))[i]) / (*(solnUpdatePtr))[i];
    else
      (*(constraintTempVector_))[i] = N_UTL_MachineDependentParams::DoubleMax();

    if (Xyce::DEBUG_NONLINEAR)
      Xyce::dout() << "ConstraintBT::updateThetaBoundPos: max: " << (*(constraintMaxVector_))[i] << "\n"
                   << "ConstraintBT::updateThetaBoundPos: soln: " << (*(solnPtr))[i] << "\n"
                   << "ConstraintBT::updateThetaBoundPos: solnUpdate: " << (*(solnUpdatePtr))[i] << "\n"
                   << "ConstraintBT::updateThetaBoundPos: constraint: " << (*(constraintTempVector_))[i];
  }

  // Find minimum
  constraintTempVector_->minValue(&thetaBoundPos_);

  if (Xyce::DEBUG_NONLINEAR)
    Xyce::dout() << "ConstraintBT::updateThetaBoundPos: thetaBoundPos_: " << thetaBoundPos_;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::updateThetaChange
// Purpose       : Updates the percentage-change bound value for the
//                 backtracking algorithm.
// Special Notes : This is implemented according to a internal
//                 communication with John Shadid (SNL) on their MPSalsa
//                 constraint backtracking implementation.
//                 This function returns:
//
//                 theta_u = min {gamma_i | oldSoln_i | / | solnUpdate_i | }
//                            i
//
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
void ConstraintBT::updateThetaChange(const N_LAS_Vector * oldSoln,
                                           const N_LAS_Vector * solnUpdate)

{
  N_LAS_Vector * solnPtr;
  N_LAS_Vector * solnUpdatePtr;

  // Initialize
  solnPtr       = const_cast<N_LAS_Vector *> (oldSoln);
  solnUpdatePtr = const_cast<N_LAS_Vector *> (solnUpdate);

  // First, form a vector of constraints...
  for (int i = 0; i < solnPtr->localLength(); ++i)
  {
    if (fabs((*(solnUpdatePtr))[i]) > N_UTL_MachineDependentParams::DoubleMin() &&
        fabs((*(solnPtr))[i]) > 0.0)
      (*(constraintTempVector_))[i] = (*(constraintChangeVector_))[i] *
        fabs((*(solnPtr))[i]) / fabs((*(solnUpdatePtr))[i]);

    else
      (*(constraintTempVector_))[i] = N_UTL_MachineDependentParams::DoubleMax();

    if (Xyce::DEBUG_NONLINEAR)
      Xyce::dout() << "ConstraintBT::updateThetaChange: change: " << (*(constraintChangeVector_))[i] << "\n"
                   << "ConstraintBT::updateThetaChange: constraint: " << (*(constraintTempVector_))[i];
  }

  // Find minimum
  constraintTempVector_->minValue(&thetaChange_);

  if (Xyce::DEBUG_NONLINEAR)
    Xyce::dout() << "ConstraintBT::updateThetaChange: thetaChange_: " << thetaChange_;
}

} // namespace Nonlinear
} // namespace Xyce
