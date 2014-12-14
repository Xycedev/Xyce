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
// Filename       : $RCSfile: N_IO_OutputterHomotopyTecplot.h,v $
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
// Revision Number: $Revision: 1.7.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:50 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterHomotopyTecplot_h
#define Xyce_N_IO_OutputterHomotopyTecplot_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

class HomotopyTecPlot : public Interface
{
public:
  HomotopyTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~HomotopyTecPlot();

private:
  HomotopyTecPlot(const HomotopyTecPlot &);
  HomotopyTecPlot &operator=(const HomotopyTecPlot &);

public:
  virtual void doSetAnalysisMode(Analysis::Analysis_Mode analysis_mode) {
    printParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doOutputTime(
    Parallel::Machine           comm, 
    const N_LAS_Vector &        solution_vector, 
    const N_LAS_Vector &        state_vector, 
    const N_LAS_Vector &        store_vector)
  {}

  virtual void doFinishOutput();
  virtual void doStartStep(int current_step, int number_of_step);
  virtual void doSteppingComplete();

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
    Parallel::Machine    comm,
    double               time,
    const std::vector<double> & fast_time_points,
    const N_LAS_Vector & solution_vector)
  {}

  virtual void doOutputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names, 
    const std::vector<double> &         param_values, 
    const N_LAS_Vector &                solution_vector);

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

private:
  void homotopyHeader(
    const std::vector<std::string> &    parameter_names, 
    const std::vector<double> &         param_values, 
    const N_LAS_Vector &                solution_vector);

private:
  OutputMgr &           outputManager_;
  PrintParameters       printParameters_;
  std::string           outFilename_;
  std::ostream *        os_;
  int                   index_;
  int                   currentStep_;
  int                   numberOfSteps_;

  Table::ColumnList     columnList_;
  Util::Op::OpList      opList_;
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterHomotopyTecplot_h
