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

#ifndef Xyce_N_NLS_ConstraintBT_h
#define Xyce_N_NLS_ConstraintBT_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_NLS_fwd.h>

// ----------   Fwd Declarations  ----------

class N_LAS_Vector;
class N_LAS_System;

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : ConstraintBT
// Purpose       : Supports constraint backtracking (damping) for the nonlinear
//                 solver.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------

class ConstraintBT
{
public:

  // ***** Constructors and destructor and general operators *****

  ConstraintBT();

  ~ConstraintBT();

  ConstraintBT & operator=(const ConstraintBT & right);

  bool initializeAll (N_LAS_System * lasSysPtr, 
		      const NLParams & nlParams);

  // ***** Methods *****

  void updateThetaBoundNeg(const N_LAS_Vector * oldSoln,
                           const N_LAS_Vector * solnUpdate);
  void updateThetaBoundPos(const N_LAS_Vector * oldSoln,
                           const N_LAS_Vector * solnUpdate);
  void updateThetaChange(const N_LAS_Vector * oldSoln,
                         const N_LAS_Vector * solnUpdate);

  // ***** Accessor Methods *****

  inline void   setThetaBoundNeg(double value);
  inline void   resetThetaBoundNeg();
  inline double getThetaBoundNeg() const;

  inline void   setThetaBoundPos(double value);
  inline void   resetThetaBoundPos();
  inline double getThetaBoundPos() const;

  inline void   setThetaChange(double value);
  inline void   resetThetaChange();
  inline double getThetaChange() const;

protected:

  // Global bounds for the constraint backtracking - the names correspond
  // roughly to those used by John Shadid (SNL) in his writeup on the MPSalsa
  // implementation.
  double thetaBoundNeg_;
  double thetaBoundPos_;
  double thetaChange_;

private:

  // Copy constructor - private until implemented...
  ConstraintBT(const ConstraintBT & right);

  int operator==(const ConstraintBT & right) const;
  int operator!=(const ConstraintBT & right) const;

  // Constraint backtracking vectors.
  N_LAS_Vector  * constraintMinVector_;
  N_LAS_Vector  * constraintMaxVector_;
  N_LAS_Vector  * constraintChangeVector_;
  N_LAS_Vector  * constraintTempVector_;

};

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::setThetaBoundNeg
// Purpose       : Accessor method to set the negative bound.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline void ConstraintBT::setThetaBoundNeg(double value)
{
  thetaBoundNeg_ = value;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::resetThetaBoundNeg
// Purpose       : Accessor method to reset the negative bound to the default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline void ConstraintBT::resetThetaBoundNeg()
{
  thetaBoundNeg_ = 1.0;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::getThetaBoundNeg
// Purpose       : Accessor method which returns the negative bound constraint.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline double ConstraintBT::getThetaBoundNeg() const
{
  return thetaBoundNeg_;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::setThetaBoundPos
// Purpose       : Accessor method to set the positive bound.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline void ConstraintBT::setThetaBoundPos(double value)
{
  thetaBoundPos_ = value;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::resetThetaBoundPos
// Purpose       : Accessor method to reset the positive bound to the default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline void ConstraintBT::resetThetaBoundPos()
{
  thetaBoundPos_ = 1.0;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::getThetaBoundPos
// Purpose       : Accessor method which returns the positive bound constraint.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline double ConstraintBT::getThetaBoundPos() const
{
  return thetaBoundPos_;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::setThetaChange
// Purpose       : Accessor method to set the percentage change bound.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline void ConstraintBT::setThetaChange(double value)
{
  thetaChange_ = value;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::resetThetaChange
// Purpose       : Accessor method to reset the percentage change to the
//                 default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline void ConstraintBT::resetThetaChange()
{
  thetaChange_ = 1.0;
}

//-----------------------------------------------------------------------------
// Function      : ConstraintBT::getThetaChange
// Purpose       : Accessor method which returns the percentage change
//                 constraint.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/26/01
//-----------------------------------------------------------------------------
inline double ConstraintBT::getThetaChange() const
{
  return thetaChange_;
}

} // namespace Nonlinear
} // namespace Xyce

typedef Xyce::Nonlinear::ConstraintBT N_NLS_ConstraintBT;

#endif
