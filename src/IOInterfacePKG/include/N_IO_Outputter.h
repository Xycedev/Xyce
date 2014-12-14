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
// Filename       : $RCSfile: N_IO_Outputter.h,v $
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
// Revision Number: $Revision: 1.50.2.3 $
//
// Revision Date  : $Date: 2014/08/26 18:54:29 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Outputter_h
#define Xyce_N_IO_Outputter_h

#include <list>
#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_ANP_fwd.h>

#include <N_UTL_NetlistLocation.h>
#include <N_UTL_Param.h>

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace IO {

typedef std::list<Util::Param> ParameterList;

namespace Format {
enum Format {STD, TECPLOT, PROBE, CSV, RAW, RAW_ASCII, DAKOTA};
}

//-----------------------------------------------------------------------------
// Class         : Table
//
// Purpose       : This struct manages the header data
//
// Special Notes :
//
// Creator       : Todd Coffey, 1414
// Creation Date : 9/19/08
//-----------------------------------------------------------------------------
struct Table
{
  // Enum for column justification (left/center/right) in std header output
  enum Justification
    {
      JUSTIFICATION_LEFT,
      JUSTIFICATION_CENTER,
      JUSTIFICATION_RIGHT,
      JUSTIFICATION_NONE
    };

  struct Column
  {
    Column()
      : name_(),
        format_(std::ios_base::scientific),
        width_(17),
        precision_(9),
        justification_(JUSTIFICATION_LEFT)
    {}

    Column(const Column &column)
      : name_(column.name_),
        format_(column.format_),
        width_(column.width_),
        precision_(column.precision_),
        justification_(column.justification_)
    {}

    Column(const std::string &name, std::ios_base::fmtflags format, int width, int precision, Justification justification)
      : name_(name),
        format_(format),
        width_(width),
        precision_(precision),
        justification_(justification)
    {}

    std::string             name_;
    std::ios_base::fmtflags format_;
    int                     width_;
    int                     precision_;
    Justification           justification_;
  };

  typedef std::vector<Column> ColumnList;

  Table()
  {}

  Table(const Table &table)
    : columnList_(table.columnList_.begin(), table.columnList_.end())
  {}

  Table &operator=(const Table &table)
  {
    columnList_.assign(table.columnList_.begin(), table.columnList_.end());

    return *this;
  }

  virtual ~Table()
  {}

  void addColumn(std::string name, std::ios_base::fmtflags format, int width, int precision, Justification justification)
  {
    columnList_.push_back(Column(name, format, width, precision, justification));
  }

  void addColumn(std::string name, int width, int precision, Justification justification)
  {
    columnList_.push_back(Column(name, std::ios_base::scientific, width, precision, justification));
  }

  ColumnList          columnList_;
};


//-----------------------------------------------------------------------------
// Class         : Table
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jul 17 15:16:48 2014
//-----------------------------------------------------------------------------
///
/// A fallback print parameters are generated from a none specific print
/// type.  I.E. .PRINT TRAN creates a fallback .PRINT HOMOTOPY print
/// parameters, while .PRINT HOMOTOPY does not.  Fallback print
/// parameters are removed if a specific print parameters is provided.
///
struct PrintParameters
{
  PrintParameters()
    : fallback_(false),
      filename_(),
      suffix_(),
      defaultExtension_(),
      analysisMode_(Analysis::ANP_MODE_INVALID),
      overrideRaw_(false),
      asciiRaw_(false),
      format_(Format::STD),
      printIndexColumn_(true),
      variableList_(),
      table_(),
      streamWidth_(17),
      streamPrecision_(9),
      timeWidth_(8),
      delimiter_(),
      outputTimeScaleFactor_(1.0),
      filter_(0.0),
      expandComplexTypes_(false)
  {}

  PrintParameters(const PrintParameters &print_parameters)
    : fallback_(print_parameters.fallback_),
      filename_(print_parameters.filename_),
      suffix_(print_parameters.suffix_),
      defaultExtension_(print_parameters.defaultExtension_),
      netlistLocation_(print_parameters.netlistLocation_),
      analysisMode_(print_parameters.analysisMode_),
      overrideRaw_(print_parameters.overrideRaw_),
      asciiRaw_(print_parameters.asciiRaw_),
      format_(print_parameters.format_),
      printIndexColumn_(print_parameters.printIndexColumn_),
      variableList_(print_parameters.variableList_.begin(), print_parameters.variableList_.end()),
      table_(print_parameters.table_),
      streamWidth_(print_parameters.streamWidth_),
      streamPrecision_(print_parameters.streamPrecision_),
      timeWidth_(print_parameters.timeWidth_),
      delimiter_(print_parameters.delimiter_),
      outputTimeScaleFactor_(print_parameters.outputTimeScaleFactor_),
      filter_(print_parameters.filter_),
      expandComplexTypes_(print_parameters.expandComplexTypes_)
  {}

  PrintParameters &operator=(const PrintParameters &print_parameters)
  {
    fallback_ = print_parameters.fallback_;
    filename_ = print_parameters.filename_;
    suffix_ = print_parameters.suffix_;
    defaultExtension_ = print_parameters.defaultExtension_;
    netlistLocation_ = print_parameters.netlistLocation_;
    analysisMode_ = print_parameters.analysisMode_;
    overrideRaw_ = print_parameters.overrideRaw_;
    asciiRaw_ = print_parameters.asciiRaw_;
    format_ = print_parameters.format_;
    printIndexColumn_ = print_parameters.printIndexColumn_;
    variableList_.assign(print_parameters.variableList_.begin(), print_parameters.variableList_.end());
    table_ = print_parameters.table_;
    streamWidth_ = print_parameters.streamWidth_;
    streamPrecision_ = print_parameters.streamPrecision_;
    timeWidth_ = print_parameters.timeWidth_;
    delimiter_ = print_parameters.delimiter_;
    outputTimeScaleFactor_ = print_parameters.outputTimeScaleFactor_;
    filter_ = print_parameters.filter_;
    expandComplexTypes_ = print_parameters.expandComplexTypes_;

    return *this;
  }

public:
  virtual ~PrintParameters()
  {}

  bool                          fallback_;                      ///< True if the creating PrintType is generic
  std::string                   filename_;                      ///< Output filename
  std::string                   suffix_;                        ///< Suffix (placed between filename and extension)
  std::string                   defaultExtension_;              ///< Extension if none provided

  NetlistLocation               netlistLocation_;               ///< Location in netlist creating this PrintParameters
  Analysis::Analysis_Mode       analysisMode_;                  ///< Analysis mode that triggered the current print
  bool                          overrideRaw_;                   ///< True if -r specified on command line
  bool                          asciiRaw_;                      ///< True if -a specified on command line
  Format::Format                format_;                        ///< Format specified
  bool                          printIndexColumn_;              ///< True if INDEX column is to be printed
  ParameterList                 variableList_;                  ///< Description of variables to be printed
  Table                         table_;                         ///< Formatting of table for print
  int                           streamWidth_;                   ///< Column width
  int                           streamPrecision_;               ///< Column precision
  int                           timeWidth_;                     ///< Really?
  std::string                   delimiter_;                     ///< Delimiter
  double                        outputTimeScaleFactor_;         ///< Output time in something other than seconds (such as milli-seconds)
  double                        filter_;
  bool                          expandComplexTypes_;
};

namespace Outputter {

class Interface
{
public:
  Interface()
  {}

  virtual ~Interface()
  {}

  void setAnalysisMode(Analysis::Analysis_Mode analysis_mode);

  void startStep(int step, int max_step);

  void output(
    Parallel::Machine           comm,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector);

  void outputAC(
    Parallel::Machine           comm,
    double                      frequency,
    const N_LAS_Vector &        real_solution_vector,
    const N_LAS_Vector &        imaginary_solution_vector);

  void outputHB(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const std::vector<double> & freqPoints,
    const N_LAS_BlockVector &   timeDomainSolutionVec,
    const N_LAS_BlockVector &   freqDomainSolutionVecReal,
    const N_LAS_BlockVector &   freqDomainSolutionVecImaginary,
    const N_LAS_BlockVector &   timeDomainStoreVec,
    const N_LAS_BlockVector &   freqDomainStoreVecReal,
    const N_LAS_BlockVector &   freqDomainStoreVecImaginary);

  void outputMPDE(
    Parallel::Machine           comm,
    double                      time,
    const std::vector<double> & fast_time_points,
    const N_LAS_Vector &        solution_vector);

  void outputHomotopy(
    Parallel::Machine                   comm,
    const std::vector<std::string> &    parameter_names,
    const std::vector<double> &         parameter_values,
    const N_LAS_Vector &                solution_vector);

  void outputSensitivity(
    Parallel::Machine                   comm,
    const std::vector<double> &         objective_values,
    const std::vector<double> &         direct_values,
    const std::vector<double> &         adjoint_values,
    const std::vector<double> &         scaled_direct_values,
    const std::vector<double> &         scaled_adjoint_values,
    const N_LAS_Vector &                solution_vector,
    const N_LAS_Vector &                state_vector,
    const N_LAS_Vector &                store_vector);

  void finishOutput();

  void steppingComplete();

private:
  virtual void doSetAnalysisMode(Analysis::Analysis_Mode analysis_mode) = 0;

  virtual void doOutputTime(
    Parallel::Machine           comm,
    const N_LAS_Vector &        solution_vector,
    const N_LAS_Vector &        state_vector,
    const N_LAS_Vector &        store_vector) = 0;

  virtual void doOutputFrequency(
    Parallel::Machine           comm,
    double                      frequency,
    const N_LAS_Vector &        real_solution_vector,
    const N_LAS_Vector &        imaginary_solution_vector) = 0;

  virtual void doOutputHB(
    Parallel::Machine           comm,
    const std::vector<double> & timePoints,
    const std::vector<double> & freqPoints,
    const N_LAS_BlockVector &   timeDomainSolutionVec,
    const N_LAS_BlockVector &   freqDomainSolutionVecReal,
    const N_LAS_BlockVector &   freqDomainSolutionVecImaginary,
    const N_LAS_BlockVector &   timeDomainStoreVec,
    const N_LAS_BlockVector &   freqDomainStoreVecReal,
    const N_LAS_BlockVector &   freqDomainStoreVecImaginary) = 0;

  virtual void doOutputMPDE(
    Parallel::Machine           comm,
    double                      time,
    const std::vector<double> & fast_time_points,
    const N_LAS_Vector &        solution_vector) = 0;

  virtual void doOutputHomotopy(
    Parallel::Machine                   comm,
     const std::vector<std::string> &   parameter_names,
     const std::vector<double> &        param_values,
     const N_LAS_Vector &               solution_vector) = 0;

  virtual void doOutputSensitivity(
    Parallel::Machine                    comm,
     const std::vector<double> &         objective_values,
     const std::vector<double> &         direct_values,
     const std::vector<double> &         adjoint_values,
     const std::vector<double> &         scaled_direct_values,
     const std::vector<double> &         scaled_adjoint_values,
     const N_LAS_Vector &                solution_vector,
     const N_LAS_Vector &                state_vector,
     const N_LAS_Vector &                store_vector) = 0;

  virtual void doFinishOutput() = 0;

  virtual void doStartStep(int step, int max_step) = 0;

  virtual void doSteppingComplete() = 0;
};

} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Outputter_h
