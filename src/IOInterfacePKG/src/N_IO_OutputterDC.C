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
// Filename       : $RCSfile: N_IO_OutputterDC.C,v $
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
// Revision Number: $Revision: 1.5.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputterDC.h>
#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterTimeProbe.h>
#include <N_IO_OutputterTimeRaw.h>
#include <N_IO_OutputterTimeRawASCII.h>
#include <N_IO_OutputterTimeTecplot.h>

namespace Xyce {
namespace IO {
namespace Outputter {

void enableDCOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Analysis_Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::DC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it) {
      PrintParameters dc_print_parameters = (*it);


      if (dc_print_parameters.printIndexColumn_)
        dc_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

      output_manager.fixupPrintParameters(comm, dc_print_parameters);

      Outputter::Interface *outputter_prn;
      if (dc_print_parameters.format_ == Format::STD) {
        outputter_prn = new Outputter::TimePrn(comm, output_manager, dc_print_parameters);
      }
      else if (dc_print_parameters.format_ == Format::CSV) {
        outputter_prn = new Outputter::TimeCSV(comm, output_manager, dc_print_parameters);
      }
      else if (dc_print_parameters.format_ == Format::RAW) {
        if (analysis_mode == Analysis::ANP_MODE_DC_SWEEP)
          dc_print_parameters.variableList_.push_front(Util::Param("sweep", 0.0));

        outputter_prn = new Outputter::TimeRaw(comm, output_manager, dc_print_parameters);
      }
      else if (dc_print_parameters.format_ == Format::RAW_ASCII) {
        if (analysis_mode == Analysis::ANP_MODE_DC_SWEEP)
          dc_print_parameters.variableList_.push_front(Util::Param("sweep", 0.0));

        outputter_prn = new Outputter::TimeRawAscii(comm, output_manager, dc_print_parameters);
      }
      else if (dc_print_parameters.format_ == Format::TECPLOT) {
        outputter_prn = new Outputter::TimeTecplot(comm, output_manager, dc_print_parameters);
      }
      else if (dc_print_parameters.format_ == Format::PROBE) {
        outputter_prn = new Outputter::TimeProbe(comm, output_manager, dc_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "DC output cannot be written in " << dc_print_parameters.format_ << " format, using standard format";

        outputter_prn = new Outputter::TimePrn(comm, output_manager, dc_print_parameters);
      }

      output_manager.addOutputter(PrintType::TRAN, outputter_prn);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
