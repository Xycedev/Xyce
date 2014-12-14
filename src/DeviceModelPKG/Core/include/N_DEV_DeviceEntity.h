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
//  (at your option) any later version.
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
// Filename       : $RCSfile: N_DEV_DeviceEntity.h,v $
//
// Purpose        : This file contains the device entity base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/03/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.127.2.1 $
//
// Revision Date  : $Date: 2014/08/13 20:36:35 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceEntity_h
#define Xyce_N_DEV_DeviceEntity_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_NetlistLocation.h>
#include <N_DEV_Pars.h>
#include <N_DEV_InstanceName.h>

class N_LAS_Vector;

namespace Xyce {
namespace Device {

typedef std::map<std::string, std::vector<Param>, LessNoCase> CompositeParamMap;

void populateParams(const ParameterMap &parameter_map, std::vector<Param> &param_list, CompositeParamMap &composite_param_map);
void setParameters(CompositeParam &composite_param, const std::string &pName, const Param &ndParam );
void setParameters(DeviceEntity &entity, std::vector<Param>::const_iterator begin, std::vector<Param>::const_iterator end, const DeviceOptions &device_options);

struct Depend
{
  std::string                 name;
  Util::Expression *          expr;
  union resUnion
  {
    double *                result;
    std::vector<double> *   resVec;
  } resultU;
  int                         vectorIndex;
  std::vector<double>         vals;
  std::vector<std::string>    global_params;
  int                         n_vars, lo_var;
};

//-----------------------------------------------------------------------------
// Class         : DeviceEntity
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/11/02
//-----------------------------------------------------------------------------
class DeviceEntity : public ParameterBase
{
public:
  DeviceEntity(
     ParametricData<void> &    parametric_data,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options,
     const std::string &       netlist_path,
     int                       netlist_line);

private:
  DeviceEntity(const DeviceEntity &);                   ///< No copying
  DeviceEntity &operator=(const DeviceEntity &);        ///< No assignment

public:
  virtual ~DeviceEntity();

  virtual bool processParams() = 0;

  virtual bool processInstanceParams() = 0;

  virtual CompositeParam *constructComposite(const std::string &composite_name, const std::string &param_name) 
  {
    return NULL;
  }

  bool setDefaultParam(double val);
  double getDefaultParam() const;

  bool scaleParam(const std::string & paramName, double val, double val0);
  bool scaleParam(const std::string & paramName, double val);
  bool scaleDefaultParam(double val);

  bool analyticSensitivityAvailable (const std::string & paramName);
  bool getAnalyticSensitivity ( const std::string & paramName,
                                std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec );

  bool setParam(const std::string & paramName, double val);
  bool getParam(const std::string & paramName, double & result) const;
  bool getParamBreakpoints( std::vector<Util::BreakPoint> & );

  bool updateDependentParameters(N_LAS_Vector & vars);
  bool updateDependentParameters(double temp_tmp);
  bool updateGlobalParameters(GlobalParameterMap &);
  bool updateDependentParameters();

  double setDependentParameter(Util::Param &, double *, ParameterType::ExprAccess);
  double setDependentParameter(Util::Param &, std::vector<double> *, int , ParameterType::ExprAccess);
  void setDependentParameter(Util::Param & par, Depend & dependentParam, ParameterType::ExprAccess depend);

  void setDefaultParams() 
  {
    setDefaultParameters(*this, getParameterMap().begin(), getParameterMap().end(), devOptions_);
  }

  void setParams(const std::vector<Param> & params) 
  {
    setParameters(*this, params.begin(), params.end(), devOptions_);
  }

public:
  bool given(const std::string & parameter_name) const;

  virtual std::ostream &printName(std::ostream &os) const = 0;

  void setDefaultParamName(const std::string &default_param_name) 
  {
    defaultParamName_ = default_param_name;
  }

  const std::vector<Depend> &getDependentParams() 
  {
    return dependentParams_;
  }

  void addDependentParameter(const Depend &param) {
    dependentParams_.push_back(param);
  }

  const DeviceOptions &getDeviceOptions() const 
  {
    return devOptions_;
  }

  const SolverState &getSolverState() const 
  {
    return solState_;
  }

  const NetlistLocation &netlistLocation() const 
  {
    return netlistLocation_;
  }

  const ParameterMap &getParameterMap() const 
  {
    return parametricData_.getMap();
  }

private:
  void escape(std::string &) const;
  void checkDepend(ParameterType::ExprAccess &);

private:
  std::string                 defaultParamName_;
  ParametricData<void> &      parametricData_;
  NetlistLocation             netlistLocation_;

  const SolverState &         solState_;
  Globals &                   globals_;

  const DeviceOptions &       devOptions_;
  std::vector<Depend>         dependentParams_;

protected:
  std::vector<int>            expVarGIDs;
  std::vector<int>            expVarLIDs;
  std::vector<std::string>    expVarNames;
  std::vector<double>         expVarVals;
  std::vector<double>         eVarVals;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Depend Depend;

#endif // Xyce_N_DEV_DeviceEntity_h
