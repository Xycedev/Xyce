//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-RawOverride04-94AL85000 with Sandia Corporation, the U.S.
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
// Filename       : $RCSfile: N_IO_OutputterRawOverride.C,v $
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

#include <N_IO_OutputterRawOverride.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputterFrequencyPrn.h>
#include <N_IO_OutputterFrequencyCSV.h>
#include <N_IO_OutputterOverrideRaw.h>
#include <N_IO_OutputterOverrideRawASCII.h>
// #include <N_IO_OutputterTimeProbe.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

void enableRawOverrideOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Analysis_Mode analysis_mode)
{
  PrintParameters raw_print_parameters = output_manager.getDefaultPrintParameters();

  output_manager.fixupPrintParameters(comm, raw_print_parameters);

  Outputter::Interface *outputter;
  if (raw_print_parameters.format_ == Format::RAW)
    outputter = new Outputter::OverrideRaw(comm, output_manager, raw_print_parameters);
  else
    outputter = new Outputter::OverrideRawAscii(comm, output_manager, raw_print_parameters);

  output_manager.addOutputter(PrintType::RAW_OVERRIDE, outputter);
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
