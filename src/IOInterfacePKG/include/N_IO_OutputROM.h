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
// Filename       : $RCSfile: N_IO_OutputROM.h,v $
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
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/06/26 14:13:02 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputROM_h
#define Xyce_N_IO_OutputROM_h

#include <string>

#include <Teuchos_SerialDenseMatrix.hpp>

#include <N_PDS_fwd.h>

class N_LAS_Matrix;

namespace Xyce {
namespace IO {

  void outputROM(Parallel::Machine comm, const std::string &netlist_filename,
    const Teuchos::SerialDenseMatrix<int, double>& Ghat,
    const Teuchos::SerialDenseMatrix<int, double>& Chat,
    const Teuchos::SerialDenseMatrix<int, double>& Bhat,
    const Teuchos::SerialDenseMatrix<int, double>& Lhat);

  void outputROM(Parallel::Machine comm, const std::string &netlist_filename,
    const N_LAS_Matrix& Ghat,
    const N_LAS_Matrix& Chat,
    const Teuchos::SerialDenseMatrix<int, double>& Bhat,
    const Teuchos::SerialDenseMatrix<int, double>& Lhat);



} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputROM_h
