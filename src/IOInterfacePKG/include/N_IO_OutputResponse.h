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
// Filename       : $RCSfile: N_IO_OutputResponse.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2014/08/04 17:23:44 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputResponse_h
#define Xyce_N_IO_OutputResponse_h

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace IO {

template <class T>
class OpCallback;

template <>
class OpCallback<void>
{
public:
  OpCallback(const Util::Op::Operator &op)
    : op_(op)
  {}

  virtual ~OpCallback()
  {}

  void send(
    Parallel::Machine             comm,
    const N_LAS_Vector *          real_solution_vector,
    const N_LAS_Vector *          imaginary_solution_vector,
    const N_LAS_Vector *          state_vector,
    const N_LAS_Vector *          store_vector,
    const std::vector<double> *   sens_function,
    const std::vector<double> *   sens_dOdPDirect,
    const std::vector<double> *   sens_dOdPDirectScaled,
    const std::vector<double> *   sens_dOdPAdjoint,
    const std::vector<double> *   sens_dOdPAdjointScaled) const;

  virtual void execute(complex value) const = 0;

private:  
  const Util::Op::Operator &  op_;
};

template <class T>
class OpCallback : OpCallback<void>
{
public:
  typedef void(T::*Action)(complex value) const;

  OpCallback(const T &t, Action action, Util::Op::Operator &op)
    : OpCallback<void>(op),
      t_(t),
      action_(action)
  {}

  virtual ~OpCallback()
  {}

  virtual void execute(complex value) const
  {
    (t_.*action_)(value);
  }

private:
  const T &             t_;
  Action                action_;
};

typedef std::vector<OpCallback<void> *> OpCallbackVector;

class OutputResponse
{
public:
  OutputResponse(OutputMgr &output_manager);

  ~OutputResponse();

  // routines to get/set variables that external programs such as Dakota want to query the value of.
  void setResponseFilename( std::string aFileName )
  {
    responseFileName_ = aFileName;
  }

  std::string getResponseFilename() const
  {
    return responseFileName_;
  }

  void finalizeResponseVars(double time);
  bool registerResponseVars(Parallel::Machine comm, const std::string &objString, std::vector<double> *variable_vector);
  void saveResponseVarValues(Parallel::Machine comm, double time, const N_LAS_Vector &solnVecPtr);

  template<class T>
  void addOpCallback(const T &t, typename OpCallback<T>::Action action)
  {
    callbacks_.push_back(new OpCallback<T>(t, action));
  }

  void send(
    Parallel::Machine             comm,
    const N_LAS_Vector *          real_solution_vector,
    const N_LAS_Vector *          imaginary_solution_vector,
    const N_LAS_Vector *          state_vector,
    const N_LAS_Vector *          store_vector,
    const std::vector<double> *   sens_function,
    const std::vector<double> *   sens_dOdPDirect,
    const std::vector<double> *   sens_dOdPDirectScaled,
    const std::vector<double> *   sens_dOdPAdjoint,
    const std::vector<double> *   sens_dOdPAdjointScaled);

private:
  OutputMgr &                   outputManager_;
  std::string                   responseFileName_;
  Util::Op::OpList              responseVarList_;
  std::vector<double> *         responseVarPtr_;
  std::vector<std::string>      responseNames_;
  int                           numResponseVars_;
  OpCallbackVector              callbacks_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputResponse_h
