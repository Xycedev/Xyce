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
// Filename       : $RCSfile: N_IO_OutputterMPDE.C,v $
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
#include <N_IO_OutputterMPDEPrn.h>
#include <N_IO_OutputterMPDETecplot.h>
#include <N_IO_OutputterTimeTecplot.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

void enableMPDEOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Analysis_Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::MPDE);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters mpde_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, mpde_print_parameters);

      Outputter::Interface *outputter_mpde;
      if (mpde_print_parameters.format_ == Format::STD) {
        outputter_mpde = new Outputter::MPDEPrn(comm, output_manager, mpde_print_parameters);
      }
      else if (mpde_print_parameters.format_ == Format::TECPLOT) {
        outputter_mpde = new Outputter::MPDETecplot(comm, output_manager, mpde_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "MPDE output cannot be written in " << mpde_print_parameters.format_ << " format, using standard format";

        outputter_mpde = new Outputter::MPDEPrn(comm, output_manager, mpde_print_parameters);
      }

      output_manager.addOutputter(PrintType::MPDE, outputter_mpde);
    }
  }
  

  result = output_manager.findOutputParameter(OutputType::MPDE_IC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters mpde_ic_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, mpde_ic_print_parameters);

      Outputter::Interface *outputter_mpde_ic;
      if (mpde_ic_print_parameters.format_ == Format::STD) {
        outputter_mpde_ic = new Outputter::TimePrn(comm, output_manager, mpde_ic_print_parameters);
      }
      else if (mpde_ic_print_parameters.format_ == Format::TECPLOT) {
        outputter_mpde_ic = new Outputter::TimeTecplot(comm, output_manager, mpde_ic_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "MPDE output cannot be written in " << mpde_ic_print_parameters.format_ << " format, using standard format";

        outputter_mpde_ic = new Outputter::TimePrn(comm, output_manager, mpde_ic_print_parameters);
      }

      output_manager.addOutputter(PrintType::MPDE_IC, outputter_mpde_ic);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
