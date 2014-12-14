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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_ArtificialParameters.h,v $
//
// Purpose        : Implement the MOSFET Level 1 static model
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/04/28 21:48:23 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ArtificialParameters_h
#define Xyce_N_DEV_ArtificialParameters_h

#include <Xyce_config.h>

#include <vector>

#include <N_DEV_fwd.h>

namespace Xyce {
namespace Device {
namespace ArtificialParameters {

typedef std::vector<DeviceInstance *> InstanceVector;
typedef std::map<ModelTypeId, InstanceVector> ModelTypeInstanceVectorMap;

struct ArtificialParameter
{
  SolverState &getSolverState(DeviceMgr &device_manager);
  const SolverState &getSolverState(const DeviceMgr &device_manager) const;
  DeviceOptions &getDeviceOptions(DeviceMgr &device_manager);
  const DeviceOptions &getDeviceOptions(const DeviceMgr &device_manager) const;
  ModelTypeInstanceVectorMap &getModelTypeInstanceVectorMap(DeviceMgr &device_manager);
  const ModelTypeInstanceVectorMap &getModelTypeInstanceVectorMap(const DeviceMgr &device_manager) const;
  InstanceVector &getInstanceVector(DeviceMgr &device_manager);
  const InstanceVector &getInstanceVector(const DeviceMgr &device_manager) const;

  virtual bool setValue(DeviceMgr &device_manager, double value) = 0;
  virtual double getValue(const DeviceMgr &device_manager) const = 0;
};

struct MOSFETGainScaleParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct MOSFETNLTermScaleParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct MOSFETLParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct MOSFETWParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct MOSFETSizeScaleParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct MOSFETTOXParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct BJTBFParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct BJTNFParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct BJTNRParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct BJTExpOrdParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct DiodeNParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct VsrcScaleParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct PDEAlphaParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct PDEBetaParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct PDEChargeAlphaParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct GSteppingParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct GMinParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

struct TempParam : public ArtificialParameter
{
  virtual bool setValue(DeviceMgr &device_manager, double value);
  virtual double getValue(const DeviceMgr &device_manager) const;
};

} // namespace ArtificialParameters
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_ArtificialParameters_h
