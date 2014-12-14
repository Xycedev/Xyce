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
// Filename       : $RCSfile: N_IO_OutputterLocal.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.1 $
//
// Revision Date  : $Date: 2014/08/14 23:17:20 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterLocal_h
#define Xyce_N_IO_OutputterLocal_h

#include <list>
#include <string>
#include <vector>

#include <N_IO_Outputter.h>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Demangle.h>
#include <N_UTL_Param.h>

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace IO {

namespace Outputter {

class TimeInterface : public Interface
{
public:
  virtual void doOutputFrequency(
    Parallel::Machine           comm,
    double                      frequency,
    const N_LAS_Vector &        real_solution_vector,
    const N_LAS_Vector &        imaginary_solution_vector)
  {}

  virtual void doOutputHB(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const std::vector<double> & freqPoints,
    const N_LAS_BlockVector &   timeDomainSolutionVec,
    const N_LAS_BlockVector &   freqDomainSolutionVecReal,
    const N_LAS_BlockVector &   freqDomainSolutionVecImaginary,
    const N_LAS_BlockVector &   timeDomainStoreVec,
    const N_LAS_BlockVector &   freqDomainStoreVecReal,
    const N_LAS_BlockVector &   freqDomainStoreVecImaginary)
  {}

  virtual void doOutputMPDE(
    Parallel::Machine           comm,
    double                      time,
    const std::vector<double> & fast_time_points,
    const N_LAS_Vector &        solution_vector)
  {}

  virtual void doOutputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names,
    const std::vector<double> &         param_values,
    const N_LAS_Vector &                solution_vector)
  {}

  virtual void doOutputSensitivity(
    Parallel::Machine           comm,
    const std::vector<double> & objective_values,
    const std::vector<double> & direct_values,
    const std::vector<double> & adjoint_values,
    const std::vector<double> & scaled_direct_values,
    const std::vector<double> & scaled_adjoint_values,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector)
  {}
};

class FrequencyInterface : public Interface
{
public:
  virtual void doOutputTime(
    Parallel::Machine           comm,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector)
  {}

  virtual void doOutputHB(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const std::vector<double> & freqPoints,
    const N_LAS_BlockVector &   timeDomainSolutionVec,
    const N_LAS_BlockVector &   freqDomainSolutionVecReal,
    const N_LAS_BlockVector &   freqDomainSolutionVecImaginary,
    const N_LAS_BlockVector &   timeDomainStoreVec,
    const N_LAS_BlockVector &   freqDomainStoreVecReal,
    const N_LAS_BlockVector &   freqDomainStoreVecImaginary)
  {}

  virtual void doOutputMPDE(
    Parallel::Machine           comm,
    double                      time,
    const std::vector<double> & fast_time_points,
    const N_LAS_Vector &        solution_vector)
  {}

  virtual void doOutputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names,
    const std::vector<double> &         param_values,
    const N_LAS_Vector &                solution_vector)
  {}

  virtual void doOutputSensitivity(
    Parallel::Machine           comm,
    const std::vector<double> & objective_values,
    const std::vector<double> & direct_values,
    const std::vector<double> & adjoint_values,
    const std::vector<double> & scaled_direct_values,
    const std::vector<double> & scaled_adjoint_values,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector)
  {}
};

class HBInterface : public Interface
{
public:
  virtual void doOutputFrequency(
    Parallel::Machine           comm,
    double                      frequency,
    const N_LAS_Vector &        real_solution_vector,
    const N_LAS_Vector &        imaginary_solution_vector)
  {}

  virtual void doOutputTime(
    Parallel::Machine           comm,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector)
  {}

  virtual void doOutputMPDE(
    Parallel::Machine           comm,
    double                      time,
    const std::vector<double> & fast_time_points,
    const N_LAS_Vector &        solution_vector)
  {}

  virtual void doOutputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names,
    const std::vector<double> &         param_values,
    const N_LAS_Vector &                solution_vector)
  {}

  virtual void doOutputSensitivity(
    Parallel::Machine           comm,
    const std::vector<double> & objective_values,
    const std::vector<double> & direct_values,
    const std::vector<double> & adjoint_values,
    const std::vector<double> & scaled_direct_values,
    const std::vector<double> & scaled_adjoint_values,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector)
  {}
};


void fixupColumns(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager, PrintParameters &print_parameters, Util::Op::OpList &op_list);

std::string outputFilename(const std::string &filename, const std::string &default_extension, const std::string &suffix, const std::string &net_list_filename);

std::ostream &printHeader(std::ostream &os, const PrintParameters &print_parameters);
std::ostream &printHeader(std::ostream &os, const Table::ColumnList &column_list, const std::string &delimiter);
std::ostream &printHeader(std::ostream &os, const Table::Column &column);
std::ostream &printValue(std::ostream &os, const Table::Column &column,
                         const std::string &delimiter, const int column_index, double value);

//-----------------------------------------------------------------------------
// Function      : filter
// Purpose       : Applies a filter to a double value, returning 0 if the
//                 absolute value of the data is less than the filter.
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 11/25/2013
//-----------------------------------------------------------------------------
inline double filter(double value, double filter)
{
  return std::abs(value) < filter ? 0.0 : value;
}

} // namespace Outputter

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterLocal_h
