//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-Transient04-94AL85000 with Sandia Corporation, the U.S.
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
// Filename       : $RCSfile: N_IO_OutputterHB.C,v $
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
// Revision Number: $Revision: 1.2.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputterTransient.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputterFrequencyPrn.h>
#include <N_IO_OutputterFrequencyCSV.h>
#include <N_IO_OutputterHBPrn.h>
#include <N_IO_OutputterHBCSV.h>
#include <N_IO_OutputterHBTecplot.h>
#include <N_IO_OutputterTimeTecplot.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

void enableHBOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Analysis_Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result1 = output_manager.findOutputParameter(OutputType::HB_FD);
  std::pair<OutputParameterMap::const_iterator, bool> result2 = output_manager.findOutputParameter(OutputType::HB_TD);
  if (result1.second && result2.second)
  {
    for (std::vector<PrintParameters>::const_iterator it1 = (*result1.first).second.begin(), end1 = (*result1.first).second.end(), it2 = (*result2.first).second.begin();
         it1 != end1; ++it1, ++it2) {
      PrintParameters freq_print_parameters = (*it1);
      PrintParameters time_print_parameters = (*it2);

      output_manager.fixupPrintParameters(comm, freq_print_parameters);
      output_manager.fixupPrintParameters(comm, time_print_parameters);

      Outputter::Interface *outputter_hb;
      if (freq_print_parameters.format_ == Format::STD) {
        outputter_hb = new Outputter::HBPrn(comm, output_manager, freq_print_parameters, time_print_parameters);
      }
      else if (freq_print_parameters.format_ == Format::CSV) {
        outputter_hb = new Outputter::HBCSV(comm, output_manager, freq_print_parameters, time_print_parameters);
      }
      else if (freq_print_parameters.format_ == Format::TECPLOT) {
        outputter_hb = new Outputter::HBTecPlot(comm, output_manager, freq_print_parameters, time_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "HB output cannot be written in " << freq_print_parameters.format_ << " format, using standard format";

        outputter_hb = new Outputter::HBPrn(comm, output_manager, freq_print_parameters, time_print_parameters);
      }

      output_manager.addOutputter(PrintType::HB, outputter_hb);
    }
  }

  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::HB_IC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters hb_ic_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, hb_ic_print_parameters);

      Outputter::Interface *outputter_init;
      if (hb_ic_print_parameters.format_ == Format::STD) {
        outputter_init = new Outputter::TimePrn(comm, output_manager, hb_ic_print_parameters);
      }
      else if (hb_ic_print_parameters.format_ == Format::CSV) {
        outputter_init = new Outputter::TimeCSV(comm, output_manager, hb_ic_print_parameters);
      }
      else if (hb_ic_print_parameters.format_ == Format::TECPLOT) {
        outputter_init = new Outputter::TimeTecplot(comm, output_manager, hb_ic_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "HB output cannot be written in " << hb_ic_print_parameters.format_ << " format, using standard format";

        outputter_init = new Outputter::TimePrn(comm, output_manager, hb_ic_print_parameters);
      }

      output_manager.addOutputter(PrintType::HB_IC, outputter_init);
    }
  }

  result = output_manager.findOutputParameter(OutputType::HB_STARTUP);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters hb_startup_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, hb_startup_print_parameters);

      Outputter::Interface *outputter_startup;
      if (hb_startup_print_parameters.format_ == Format::STD) {
        outputter_startup = new Outputter::TimePrn(comm, output_manager, hb_startup_print_parameters);
      }
      else if (hb_startup_print_parameters.format_ == Format::CSV) {
        outputter_startup = new Outputter::TimeCSV(comm, output_manager, hb_startup_print_parameters);
      }
      else if (hb_startup_print_parameters.format_ == Format::TECPLOT) {
        outputter_startup = new Outputter::TimeTecplot(comm, output_manager, hb_startup_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "HB output cannot be written in " << hb_startup_print_parameters.format_ << " format, using standard format";

        outputter_startup = new Outputter::TimePrn(comm, output_manager, hb_startup_print_parameters);
      }

      output_manager.addOutputter(PrintType::HB_STARTUP, outputter_startup);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
