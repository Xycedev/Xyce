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
// Filename       : $RCSfile: N_IO_Tecplot.C,v $
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
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2014/07/10 12:49:43 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_Outputter.h>
#include <N_IO_OutputMgr.h>

#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : getTecplotTimeDateStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jun 25 13:20:41 2014
//-----------------------------------------------------------------------------
///
/// 
///
/// @invariant
///
///
/// @return 
///
///
std::string getTecplotTimeDateStamp()
{
  const time_t now = time( NULL);
  char timeDate[ 40 ];

  // format for output
  strftime( timeDate, 40, "TIME= \" %I:%M:%S %p %b %d, %Y \" ", localtime( &now));

  return std::string( timeDate);
}

//-----------------------------------------------------------------------------
// Function      : tecplotTimeHeader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jun 25 13:20:35 2014
//-----------------------------------------------------------------------------
///
/// 
///
/// @invariant
///
/// @param os 
/// @param print_title 
/// @param title 
/// @param op_list 
/// @param output_manager 
///
///
void tecplotTimeHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);
  os.precision(2);

  if (print_title)
  {
    os << "TITLE = \"" << title << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the user-specified solution vars:
    for (Util::Op::OpList::const_iterator it = op_list.begin() ; it != op_list.end(); ++it)
      os << "\" " << (*it)->getName() << "\" " << std::endl;

    // output some AUXDATA
    os << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl;

    if (!output_manager.getTempSweepFlag())
    {
      os << "DATASETAUXDATA TEMP = \"" << output_manager.getCircuitTemp() << " \"" << std::endl;
    }

  } // print header calls=0

  os << "ZONE F=POINT ";


  if (output_manager.getStepParamVec().empty())
  {
    os << "T=\"Xyce data\" ";
  }
  else
  {
    os << "T= \" ";
    for (std::vector<N_ANP_SweepParam>::const_iterator it = output_manager.getStepParamVec().begin(); it != output_manager.getStepParamVec().end(); ++it)
    {
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }
  os << std::endl;

  // put in the various sweep parameters as auxdata:
  if (!output_manager.getStepParamVec().empty())
  {
    for (std::vector<N_ANP_SweepParam>::const_iterator it = output_manager.getStepParamVec().begin(); it != output_manager.getStepParamVec().end(); ++it)
    {
      // convert any ":" or "%" in the name to a "_", so as not to confuse tecplot.
      std::string name(it->name);
      std::replace(name.begin(), name.end(), '%', '_');
      std::replace(name.begin(), name.end(), ':', '_');
      os << "AUXDATA " << name << " = " << "\" " << it->currentVal << "\" ";
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : tecplotFreqHeader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Jun 25 13:21:00 2014
//-----------------------------------------------------------------------------
///
/// 
///
/// @invariant
///
/// @param os 
/// @param print_title 
/// @param title 
/// @param op_list 
/// @param output_manager 
///
///
void tecplotFreqHeader(std::ostream &os, bool print_title, const std::string title, const Util::Op::OpList &op_list, const OutputMgr &output_manager)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);
  os.precision(2);

  if (print_title)
  {
    os << " TITLE = \" Xyce Frequency Domain data, " << title << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the user-specified solution vars:
    for (Util::Op::OpList::const_iterator it = op_list.begin() ; it != op_list.end(); ++it)
    {
      os << "\" ";
      if ( (*it)->getName() == "FREQUENCY" )
      {
        os << "FREQ";
      }
      else
      {
        os << (*it)->getName() ;
      }
      os << "\" " << std::endl;
    }
    os << "DATASETAUXDATA ";
    os << getTecplotTimeDateStamp();
    os << std::endl;

    if (!output_manager.getTempSweepFlag())
    {
      os << "DATASETAUXDATA TEMP = \"" << output_manager.getCircuitTemp() << " \"" << std::endl;
    }
  }

  // output some AUXDATA
  os << "ZONE F=POINT  ";

  if (output_manager.getStepParamVec().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    for (std::vector<N_ANP_SweepParam>::const_iterator it = output_manager.getStepParamVec().begin(); it != output_manager.getStepParamVec().end(); ++it)
    {
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;

  // put in the various sweep parameters as auxdata:
  if (!output_manager.getStepParamVec().empty())
  {
    for (std::vector<N_ANP_SweepParam>::const_iterator iterParam = output_manager.getStepParamVec().begin();
    iterParam != output_manager.getStepParamVec().end();
    ++iterParam)
    {
      // convert any ":" or "%" in the name to a "_", so as not to confuse tecplot.
      std::string tmpName(iterParam->name);
      replace(tmpName.begin(), tmpName.end(), '%', '_');
      replace(tmpName.begin(), tmpName.end(), ':', '_');
      os << "AUXDATA " << tmpName << " = " << "\" " << iterParam->currentVal << "\" ";
    }
    os << std::endl;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
