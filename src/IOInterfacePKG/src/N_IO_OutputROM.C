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
// Filename       : $RCSfile: N_IO_OutputROM.C,v $
//
// Purpose        : Output Manager
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
// Revision Date  : $Date: 2014/06/26 14:13:03 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#undef HAVE_LIBPARMETIS
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_CrsMatrix.h>

#include <N_IO_OutputROM.h>
#include <N_IO_mmio.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputROM
// Purpose       : Output reduced order model to file in Matrix Market format.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 6/7/12
//-----------------------------------------------------------------------------
void outputROM(Parallel::Machine comm, const std::string &netlist_filename,
      const Teuchos::SerialDenseMatrix<int, double>& Ghat,
      const Teuchos::SerialDenseMatrix<int, double>& Chat,
      const Teuchos::SerialDenseMatrix<int, double>& Bhat,
      const Teuchos::SerialDenseMatrix<int, double>& Lhat
     )
{
  if (Parallel::rank(comm) == 0)
  {
    // Eliminate the ".cir" suffix from the filename.
    // Leaving ".cir" in the netlist_filename causes issues when reading the files back in.
    int pos = netlist_filename.find(".cir");
    std::string baseoutputname;
    if (pos > 0)
    {
      baseoutputname = netlist_filename.substr(0, pos);
    }
    else
    {
      baseoutputname = netlist_filename;
    }

    // Open files for Ghat, Chat, Bhat, and Lhat
    FILE *c_file, *g_file, *b_file, *l_file;
    MMIO::MM_typecode matcode;
    std::string cfile = baseoutputname + ".Chat";
    std::string gfile = baseoutputname + ".Ghat";
    std::string bfile = baseoutputname + ".Bhat";
    std::string lfile = baseoutputname + ".Lhat";
    c_file = fopen(cfile.c_str(), "w");
    g_file = fopen(gfile.c_str(), "w");
    b_file = fopen(bfile.c_str(), "w");
    l_file = fopen(lfile.c_str(), "w");
    if (c_file == NULL || g_file == NULL || b_file == NULL || l_file == NULL)
    {
      Report::DevelFatal0() << "Cannot open one of the ROM files for output: "
        << cfile << ", " << gfile << ", " << bfile << ", " << lfile;
    }
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_general(&matcode);
    mm_set_real(&matcode);

    int ret = 0;

    // Write the headers
    ret = MMIO::mm_write_banner(g_file, matcode);
    ret = MMIO::mm_write_banner(c_file, matcode);
    ret = MMIO::mm_write_banner(b_file, matcode);
    ret = MMIO::mm_write_banner(l_file, matcode);

    // Write the matrix array sizes
    ret = MMIO::mm_write_mtx_array_size(g_file, Ghat.numRows(), Ghat.numCols());
    ret = MMIO::mm_write_mtx_array_size(c_file, Chat.numRows(), Chat.numCols());
    ret = MMIO::mm_write_mtx_array_size(b_file, Bhat.numRows(), Bhat.numCols());
    ret = MMIO::mm_write_mtx_array_size(l_file, Lhat.numRows(), Lhat.numCols());

    // Write Ghat
    for (int j=0; j<Ghat.numCols(); j++) {
      for (int i=0; i<Ghat.numRows(); i++) {
        fprintf(g_file, "%22.16e\n", Ghat(i, j));
      }
    }

    // Write Chat
    for (int j=0; j<Chat.numCols(); j++) {
      for (int i=0; i<Chat.numRows(); i++) {
        fprintf(c_file, "%22.16e\n", Chat(i, j));
      }
    }

    // Write Bhat
    for (int j=0; j<Bhat.numCols(); j++) {
      for (int i=0; i<Bhat.numRows(); i++) {
        fprintf(b_file, "%22.16e\n", Bhat(i, j));
      }
    }

    // Write Lhat
    for (int j=0; j<Lhat.numCols(); j++) {
      for (int i=0; i<Lhat.numRows(); i++) {
        fprintf(l_file, "%22.16e\n", Lhat(i, j));
      }
    }

    // Close the files
    fclose(g_file);
    fclose(c_file);
    fclose(b_file);
    fclose(l_file);
  } // end proc0 check
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputROM
// Purpose       : Output reduced order model to file in Matrix Market format.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 6/7/12
//-----------------------------------------------------------------------------
void outputROM(Parallel::Machine comm, const std::string &netlist_filename,
  const N_LAS_Matrix& Ghat,
  const N_LAS_Matrix& Chat,
  const Teuchos::SerialDenseMatrix<int, double>& Bhat,
  const Teuchos::SerialDenseMatrix<int, double>& Lhat)
{
  // Eliminate the ".cir" suffix from the filename.
  // Leaving ".cir" in the netlist_filename causes issues when reading the files back in.
  int pos = netlist_filename.find(".cir");
  std::string baseoutputname;
  if (pos > 0)
  {
    baseoutputname = netlist_filename.substr(0, pos);
  }
  else
  {
    baseoutputname = netlist_filename;
  }

  // Create file string for Chat and Ghat
  std::string gfile = baseoutputname + ".Ghat";
  std::string cfile = baseoutputname + ".Chat";

  // Get Epetra_CrsMatrix objects from Ghat and Chat
  Epetra_CrsMatrix& epetraGhat =(const_cast<N_LAS_Matrix*>(&Ghat))->epetraObj();
  Epetra_CrsMatrix& epetraChat =(const_cast<N_LAS_Matrix*>(&Chat))->epetraObj();

  // Write out objects using EpetraExt
  EpetraExt::RowMatrixToMatrixMarketFile(gfile.c_str(), epetraGhat);
  EpetraExt::RowMatrixToMatrixMarketFile(cfile.c_str(), epetraChat);

  // Write out Bhat and Lhat.
  // NOTE:  Only do this on one processor if running in parallel.

  if (Parallel::rank(comm) == 0)
  {

    // Open files for Bhat, and Lhat
    FILE *b_file, *l_file;
    MMIO::MM_typecode matcode;

    std::string bfile = baseoutputname + ".Bhat";
    std::string lfile = baseoutputname + ".Lhat";

    b_file = fopen(bfile.c_str(), "w");
    l_file = fopen(lfile.c_str(), "w");
    if (b_file == NULL || l_file == NULL)
    {
      Report::DevelFatal0() << "Cannot open one of the ROM files for output: "
        << bfile << ", " << lfile;
    }
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_general(&matcode);
    mm_set_real(&matcode);

    int ret = 0;

    // Write the headers
    ret = MMIO::mm_write_banner(b_file, matcode);
    ret = MMIO::mm_write_banner(l_file, matcode);

    // Write the matrix array sizes
    ret = MMIO::mm_write_mtx_array_size(b_file, Bhat.numRows(), Bhat.numCols());
    ret = MMIO::mm_write_mtx_array_size(l_file, Lhat.numRows(), Lhat.numCols());

    // Write Bhat
    for (int j=0; j<Bhat.numCols(); j++) {
      for (int i=0; i<Bhat.numRows(); i++) {
        fprintf(b_file, "%22.16e\n", Bhat(i, j));
      }
    }

    // Write Lhat
    for (int j=0; j<Lhat.numCols(); j++) {
      for (int i=0; i<Lhat.numRows(); i++) {
        fprintf(l_file, "%22.16e\n", Lhat(i, j));
      }
    }

    // Close the files
    fclose(b_file);
    fclose(l_file);
  } // end proc0 check
}

} // namespace IO
} // namespace Xyce
