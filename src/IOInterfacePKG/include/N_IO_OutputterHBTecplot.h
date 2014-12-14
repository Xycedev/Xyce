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
// Filename       : $RCSfile: N_IO_OutputterHBTecplot.h,v $
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
// Revision Number: $Revision: 1.4.2.3 $
//
// Revision Date  : $Date: 2014/08/25 20:12:50 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterHBTecplot_h
#define Xyce_N_IO_OutputterHBTecplot_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

class HBTecPlot : public HBInterface
{
public:
  HBTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters);

  virtual ~HBTecPlot();

private:
  HBTecPlot(const HBTecPlot &);
  HBTecPlot &operator=(const HBTecPlot &);

public:
  virtual void doSetAnalysisMode(Analysis::Analysis_Mode analysis_mode) {
    freqPrintParameters_.analysisMode_ = analysis_mode;
    timePrintParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doOutputHB(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const std::vector<double> & freqPoints,
    const N_LAS_BlockVector &   timeDomainSolutionVec,
    const N_LAS_BlockVector &   freqDomainSolutionVecReal,
    const N_LAS_BlockVector &   freqDomainSolutionVecImaginary,
    const N_LAS_BlockVector &   timeDomainStoreVec,
    const N_LAS_BlockVector &   freqDomainStoreVecReal,
    const N_LAS_BlockVector &   freqDomainStoreVecImaginary);

  virtual void doFinishOutput();

  virtual void doStartStep(int current_step, int number_of_step);

  virtual void doSteppingComplete();

private:
  OutputMgr &           outputManager_;
  PrintParameters       freqPrintParameters_;
  PrintParameters       timePrintParameters_;
  std::string           timeFilename_;
  std::string           freqFilename_;
  std::ostream *        tos_;
  std::ostream *        fos_;
  int                   index_;
  int                   currentStep_;
  int                   numberOfSteps_;

  Util::Op::OpList      timeOpList_;
  Util::Op::OpList      freqOpList_;
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterHBTecplot_h
