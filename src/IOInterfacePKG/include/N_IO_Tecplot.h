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
// Filename       : $RCSfile: N_IO_Tecplot.h,v $
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
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/07/10 12:49:42 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_Tecplot_h
#define Xyce_N_IO_Tecplot_h

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : getTecplotTimeDateStamp
// Purpose       : Get current date and time and format for .PRINT output
// Special Notes : tecplot version of getTimeDateStamp.
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 6/14/2013
//-----------------------------------------------------------------------------
std::string getTecplotTimeDateStamp();

void tecplotTimeHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager);
void tecplotFreqHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager);

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Tecplot_h
