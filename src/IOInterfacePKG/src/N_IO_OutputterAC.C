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
// Filename       : $RCSfile: N_IO_OutputterAC.C,v $
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
// Revision Number: $Revision: 1.4.2.2 $
//
// Revision Date  : $Date: 2014/08/25 20:12:51 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputterAC.h>
#include <N_IO_OutputterFrequencyCSV.h>
#include <N_IO_OutputterFrequencyPrn.h>
#include <N_IO_OutputterFrequencyProbe.h>
#include <N_IO_OutputterFrequencyRaw.h>
#include <N_IO_OutputterFrequencyRawASCII.h>
#include <N_IO_OutputterFrequencyTecplot.h>
#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterTimeProbe.h>
#include <N_IO_OutputterTimeRaw.h>
#include <N_IO_OutputterTimeRawASCII.h>
#include <N_IO_OutputterTimeTecplot.h>

namespace Xyce {
namespace IO {
namespace Outputter {

void enableACOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Analysis_Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::AC_IC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it) {
      PrintParameters ac_ic_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, ac_ic_print_parameters);

      Outputter::Interface *outputter_prn;
      if (ac_ic_print_parameters.format_ == Format::STD) {
        ac_ic_print_parameters.defaultExtension_ = ".TD.prn";
        outputter_prn = new Outputter::TimePrn(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::CSV) {
        ac_ic_print_parameters.defaultExtension_ = ".TD.csv";
        outputter_prn = new Outputter::TimeCSV(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::PROBE) {
        ac_ic_print_parameters.defaultExtension_ = ".TD.csd";
        outputter_prn = new Outputter::TimeProbe(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::TECPLOT) {
        ac_ic_print_parameters.defaultExtension_ = ".TD.dat";
        outputter_prn = new Outputter::TimeTecplot(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::RAW)
      {
        ac_ic_print_parameters.defaultExtension_ = ".raw";
        outputter_prn = new Outputter::TimeRaw(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::RAW_ASCII)
      {
        ac_ic_print_parameters.defaultExtension_ = ".raw";
        outputter_prn = new Outputter::TimeRawAscii(comm, output_manager, ac_ic_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "AC output cannot be written in " << ac_ic_print_parameters.format_ << " format, using standard format";

        outputter_prn = new Outputter::TimePrn(comm, output_manager, ac_ic_print_parameters);
      }

      output_manager.addOutputter(PrintType::AC_IC, outputter_prn);
    }
  }

  result = output_manager.findOutputParameter(OutputType::AC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it) {
      PrintParameters ac_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, ac_print_parameters);

      Outputter::Interface *outputter_fd;
      if (ac_print_parameters.format_ == Format::STD) {
        outputter_fd = new Outputter::FrequencyPrn(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::CSV) {
        outputter_fd = new Outputter::FrequencyCSV(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::PROBE) {
        outputter_fd = new Outputter::FrequencyProbe(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::TECPLOT) {
        outputter_fd = new Outputter::FrequencyTecplot(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::RAW) {
          outputter_fd = new Outputter::FrequencyRaw(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::RAW_ASCII) {
        outputter_fd = new Outputter::FrequencyRawAscii(comm, output_manager, ac_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "AC output cannot be written in " << ac_print_parameters.format_ << " format, using standard format";

        outputter_fd = new Outputter::FrequencyPrn(comm, output_manager, ac_print_parameters);
      }

      output_manager.addOutputter(PrintType::AC, outputter_fd);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
