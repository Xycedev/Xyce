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
// Filename       : $RCSfile: N_IO_InitialConditions.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4 $
//
// Revision Date  : $Date: 2014/08/05 21:41:38 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_InitialConditions_h
#define Xyce_N_IO_InitialConditions_h

namespace Xyce {
namespace IO {

struct InitialConditionsData
{
    enum {IC_TYPE_UNDEFINED, IC_TYPE_DCOP_RESTART, IC_TYPE_IC, IC_TYPE_NODESET};

    InitialConditionsData()
      : icType_(IC_TYPE_UNDEFINED),
        output_op_(false),
        output_op_file_(""),
        input_op_(false),
        input_op_file_(""),
        ICflag_(false),
        nodesetflag_(false),
        saveFlag_(false),
        saveFileLevel_("ALL"),
        saveFileType_(".NODESET"),
        saveOutputFile_(""),
        total_soln_(0),
        op_found_(0)
    {}

    int icType_;
    bool output_op_;
    std::string output_op_file_;
    bool input_op_;
    std::string input_op_file_;
    bool ICflag_;
    bool nodesetflag_;
    bool saveFlag_;
    std::string saveFileLevel_;
    std::string saveFileType_;
    std::string saveOutputFile_;
    int total_soln_;
    int op_found_;

    NodeNamePairMap opData_;
};

void writeDCOPRestart(Parallel::Machine comm, std::ofstream &os, const NodeNamePairMap &all_nodes);
bool readDCOPRestart(Parallel::Machine comm, std::ifstream &is, const NodeNamePairMap &all_nodes, N_LAS_Vector & solnVec, N_LAS_Vector & flagVec, NodeNamePairMap &opData_, int &op_found_, int &total_soln_);

bool setupIC_or_NODESET(Parallel::Machine comm, const NodeNamePairMap &all_nodes, N_LAS_Vector & solnVec, N_LAS_Vector & flagVec, int icType, std::vector<Util::OptionBlock> & initBlockVec, NodeNamePairMap &opData_, int &op_found_, int &total_soln_);

void outputIC_or_NODESET (Parallel::Machine comm, std::ofstream &os, const std::string &saveFileType_, NodeNamePairMap &all_nodes, const N_LAS_Vector & solution);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_InitialConditions_h
