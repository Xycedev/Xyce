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
// Filename       : $RCSfile: N_NLS_Sensitivity.h,v $
//
// Purpose        : This file contains the sensitivity class.   It mostly
//                  manages the calculations of direct (and possibly later,
//                  adjoint) sensitivities.
//
// Special Notes  : The main reason that this class is derived from
//                  N_NLS_NonLinearSolver is that this class needs to
//                  do a series of linear solves, using the jacobian
//                  matrix.  This seemed similar enough to the requirements
//                  of a nonlinear solver to have one derived off the other.
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/30/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.56.2.1 $
//
// Revision Date  : $Date: 2014/08/19 16:28:07 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_Sensitivity_h
#define Xyce_N_NLS_Sensitivity_h

#include<vector>

#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>
#include <N_NLS_NLParams.h>
#include <N_NLS_NonLinearSolver.h>

namespace Xyce {
namespace Nonlinear {

enum sensDiffMode
{
  SENS_FWD,
  SENS_REV,
  SENS_CNT,
  NUM_DIFF_MODES
};

//-----------------------------------------------------------------------------
// Class         : Sensitivity
// Purpose       :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------

class Sensitivity : public NonLinearSolver
{
public:
  Sensitivity ( NonLinearSolver & nls_, 
                      N_TOP_Topology & top_,
                      N_IO_CmdParse & cp);

  ~Sensitivity ();

   bool icSensitivity (
       std::vector<double> & objectiveVec,
       std::vector<double> & dOdpVec, 
       std::vector<double> & dOdpAdjVec,
       std::vector<double> & scaled_dOdpVec, 
       std::vector<double> & scaled_dOdpAdjVec);

   int solve (NonLinearSolver * nlsTmpPtr=NULL) {return -1;};
   int solve (
       std::vector<double> & objectiveVec,
       std::vector<double> & dOdpVec, 
       std::vector<double> & dOdpAdjVec,
       std::vector<double> & scaled_dOdpVec, 
       std::vector<double> & scaled_dOdpAdjVec);

   int solveDirect  ();
   int solveAdjoint ();

   void stdOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities
       );

   void fileOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities
       );

   void dakOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities
       );

   bool loadSensitivityResiduals ();

   bool calcObjFuncDerivs ();

   bool setOptions(const N_UTL_OptionBlock& OB);
   bool setSensitivityOptions(const N_UTL_OptionBlock& OB);
   bool setTranOptions(const N_UTL_OptionBlock& OB);
   bool setHBOptions(const N_UTL_OptionBlock& OB);
   
   // Note, many of the following are here b/c they are purely
   // virtual functions of the nonlinear solver class.
   // They don't have much meaning here.  It may turn out that
   // having this class derive off of the NonLinearSolver
   // class doesn't make much sense.  If so, I'll change it later. ERK

   int getNumIterations() const;
#ifdef Xyce_DEBUG_NONLINEAR
   int getDebugLevel() const;
   bool getScreenOutputFlag() const;
   double getDebugMinTime() const;
   double getDebugMaxTime() const;
   int getDebugMinTimeStep() const;
   int getDebugMaxTimeStep() const;
   bool getMMFormat () const;
#endif
   double getMaxNormF() const;
   int getMaxNormFindex() const;

   int getContinuationStep() const;
   int getParameterNumber() const;
   bool isFirstContinuationParam() const;
   bool isFirstSolveComplete() const;
   void setAnalysisMode(AnalysisMode mode);

protected:
private:

public:
protected:
private:

  //bool allocateddXVec_;
  int debugLevel_;
  int solutionSize_;
  bool solveDirectFlag_;
  bool solveAdjointFlag_;
  bool outputScaledFlag_; // include scaled sensitivities in IO 
  bool outputUnscaledFlag_; // include unscaled sensitivities in IO
  int maxParamStringSize_;

  bool stdOutputFlag_;
  bool fileOutputFlag_;
  bool dakotaFileOutputFlag_;
  bool forceFD_;
  int numSolves_;

  // expression related stuff:
  int difference;
  bool objFuncGiven_;
  bool objFuncGIDsetup_;
  int            expNumVars_;
  std::vector<std::string> expVarNames_;
  std::vector<int>    expVarGIDs_;
  std::vector<int>    expVarLocal_;
  std::vector<double> expVarVals_;
  std::vector<double> expVarDerivs_;
  double         expVal_;
  std::string objFuncString_;

  double curValue_;   // current value of the variable.
  double objFuncEval_;// value of the evaluated objective function.
  double dOdp_;
  double sqrtEta_;
  bool sqrtEtaGiven_;

  N_LAS_Vector* dOdXVectorPtr_; // size of solution vector.

  N_LAS_Vector * lambdaVectorPtr_;
  N_LAS_Vector * savedRHSVectorPtr_;
  N_LAS_Vector * savedNewtonVectorPtr_;

  N_LAS_Vector * origFVectorPtr_;
  N_LAS_Vector * pertFVectorPtr_;

  N_LAS_Vector * origQVectorPtr_;
  N_LAS_Vector * pertQVectorPtr_;

  N_LAS_Vector * origBVectorPtr_;
  N_LAS_Vector * pertBVectorPtr_;

  NonLinearSolver & nls_;

  N_TOP_Topology & top_;

  N_UTL_Expression * expPtr_;

  int numSensParams_;
  std::vector<std::string> paramNameVec_;
};

//-----------------------------------------------------------------------------
// Function      : Sensitivity::getNumIterations
// Purpose       : doesn't do anything, is just a placeholder.
// Special Notes : This one may be needed later, I'm not sure.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
inline int Sensitivity::getNumIterations() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline int Sensitivity::getContinuationStep() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline int Sensitivity::getParameterNumber() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline bool Sensitivity::isFirstContinuationParam() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline bool Sensitivity::isFirstSolveComplete() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setAnalysisMode
// Purpose       : doesn't do anything, is just a placeholder.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
inline void Sensitivity::setAnalysisMode(AnalysisMode mode)
{

}

#ifdef Xyce_DEBUG_NONLINEAR
//-----------------------------------------------------------------------------
// Function      : Sensitivity::getDebugLevel
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
inline int Sensitivity::getDebugLevel() const
{
  return -100;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::getScreenOutputFlag
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
inline bool Sensitivity::getScreenOutputFlag () const
{
  return false;
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMinTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double Sensitivity::getDebugMinTime() const
{
  return 0.0;
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMaxTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double Sensitivity::getDebugMaxTime() const
{
  return N_UTL_MachineDependentParams::DoubleMax();
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMinTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int Sensitivity::getDebugMinTimeStep() const
{
  return 0;
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMaxTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int Sensitivity::getDebugMaxTimeStep() const
{
  return N_UTL_MachineDependentParams::IntMax();
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getMMFormat
//
// Return Type   : bool
//---------------------------------------------------------------------------
inline bool Sensitivity::getMMFormat () const
{
  return false;
}

#endif // debug nonlin

} // namespace Nonlinear
} // namespace Xyce

typedef Xyce::Nonlinear::Sensitivity N_NLS_Sensitivity;

#endif // Xyce_N_NLS_Sensitivity_h

