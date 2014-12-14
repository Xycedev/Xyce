//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_IO_OutputterSensitivity.C,v $
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
// Revision Number: $Revision: 1.2.2.4 $
//
// Revision Date  : $Date: 2014/08/26 22:46:45 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputterSensitivity.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputterFrequencyPrn.h>
#include <N_IO_OutputterFrequencyCSV.h>
#include <N_IO_OutputterSensitivityPrn.h>
#include <N_IO_OutputterSensitivityTecplot.h>
#include <N_IO_OutputterSensitivityDakota.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void enableSensitivityOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Analysis_Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::SENS);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) 
    {
      PrintParameters sensitivity_print_parameters = (*it);

      if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
        sensitivity_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      if (sensitivity_print_parameters.printIndexColumn_)
        sensitivity_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

      output_manager.fixupPrintParameters(comm, sensitivity_print_parameters);

      Outputter::Interface *outputter;
      if (sensitivity_print_parameters.format_ == Format::STD) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
        outputter = new Outputter::SensitivityPrn(comm, output_manager, sensitivity_print_parameters);
      }
      else if (sensitivity_print_parameters.format_ == Format::TECPLOT) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".SENS.dat";
        outputter = new Outputter::SensitivityTecPlot(comm, output_manager, sensitivity_print_parameters);
      }
      else if (sensitivity_print_parameters.format_ == Format::DAKOTA) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".SENS.txt";
        outputter = new Outputter::SensitivityDakota(comm, output_manager, sensitivity_print_parameters);
      }
      else
      {
        Report::UserWarning0()
          << "Sensitivity output cannot be written in requested format, using standard format";

        sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
        outputter = new Outputter::SensitivityPrn(comm, output_manager, sensitivity_print_parameters);
      }

      output_manager.addOutputter(PrintType::SENS, outputter);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
