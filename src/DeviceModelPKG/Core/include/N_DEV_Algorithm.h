//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Algorithm.h,v $
//
// Purpose        :
//
//
//
// Special Notes  :
//
//
// Creator        :
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10 $
//
// Revision Date  : $Date: 2014/04/23 13:32:01 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Algorithm_h
#define Xyce_N_DEV_Algorithm_h

#include <iterator>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_DEV_Device.h>

namespace Xyce {
namespace Device {


// Ultimately making std::for_each work would be most awesome.  But that will take some doing.

//-----------------------------------------------------------------------------
// Function      : forEachInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:04:01 2014
//-----------------------------------------------------------------------------
///
/// Call forEachInstance on object d, passing operator op
///
/// @param d Object to call forEachInstance on
/// @param op Operator to call passing in DeviceInstance pointer
///
template <class D, class Op>
void forEachInstance(const D &d, Op op) {
  d.forEachInstance(op);
}


//-----------------------------------------------------------------------------
// Function      : forEachModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:04:01 2014
//-----------------------------------------------------------------------------
///
/// Call forEachModel on object d, passing operator op
///
/// @param d Object to call forEachModel on
/// @param op Operator to call passing in DeviceModel pointer
///
template <class D, class Op>
void forEachModel(const D &d, Op op) {
  d.forEachModel(op);
}


//-----------------------------------------------------------------------------
// Function      : getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 11:12:37 2014
//-----------------------------------------------------------------------------
///
/// Returns the name of the specified object
///
/// @param c object to return associated name
///
/// @return name associated with the object
///
template <class C>
const std::string &getName(const C *c);


//-----------------------------------------------------------------------------
// Class         : DeviceInstanceOutIteratorOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 10:58:03 2014
//-----------------------------------------------------------------------------
///
/// Operator to populate a container via an output iterator with device instances.
///
template <class Out>
struct DeviceInstanceOutIteratorOp: public DeviceInstanceOp
{
    //-----------------------------------------------------------------------------
    // Function      : DeviceInstanceOutIteratorOp
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 10:59:27 2014
    //-----------------------------------------------------------------------------
    ///
    /// Copies output iterator into operator
    ///
    DeviceInstanceOutIteratorOp(Out it)
      : it_(it)
    {}

    //-----------------------------------------------------------------------------
    // Function      : operator()
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 11:00:21 2014
    //-----------------------------------------------------------------------------
    ///
    /// Copies instance pointer through output iterator
    ///
    virtual bool operator()(DeviceInstance *instance) {
      (*it_)++ = instance;

      return true;
    }

    Out it_;                    ///< Output iterator
};


//-----------------------------------------------------------------------------
// Class         : DeviceModelOutIteratorOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 10:58:03 2014
//-----------------------------------------------------------------------------
///
/// Operator to populate a container via an output iterator with device models.
///
template <class Out>
struct DeviceModelOutIteratorOp: public DeviceModelOp
{
    //-----------------------------------------------------------------------------
    // Function      : DeviceModelOutIteratorOp
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 10:59:27 2014
    //-----------------------------------------------------------------------------
    ///
    /// Copies output iterator into operator
    ///
    DeviceModelOutIteratorOp(Out it)
      : it_(it)
    {}

    //-----------------------------------------------------------------------------
    // Function      : operator()
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 11:00:21 2014
    //-----------------------------------------------------------------------------
    ///
    /// Copies model pointer through output iterator
    ///
    virtual bool operator()(DeviceModel *model) {
      (*it_)++ = model;

      return true;
    }

    Out it_;                    ///< Output iterator
};


//-----------------------------------------------------------------------------
// Class         : DeviceInstanceNameOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 10:58:03 2014
//-----------------------------------------------------------------------------
///
/// Operator to populate a container of strings with the names of the device instances.
///
template<class Out>
struct DeviceInstanceNameOp: public DeviceInstanceOp
{
    //-----------------------------------------------------------------------------
    // Function      : DeviceInstanceNameOp
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 10:59:27 2014
    //-----------------------------------------------------------------------------
    ///
    /// Copies output iterator to it_
    ///
    DeviceInstanceNameOp(Out it)
      : it_(it)
    {}

    //-----------------------------------------------------------------------------
    // Function      : operator
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 12:08:10 2014
    //-----------------------------------------------------------------------------
    ///
    /// Write the result of getName(instance) to the output iterator
    ///
    virtual bool operator()(DeviceInstance *device_instance) {
      (*it_)++ = getName(device_instance);

      return true;
    }

    Out         it_;            ///< Output iterator
};


//-----------------------------------------------------------------------------
// Function      : MapOp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:19:13 2014
//-----------------------------------------------------------------------------
///
/// operator to populate a map from device instance name to device instance pointer
///
template<class C>
struct MapOp : public DeviceInstanceOp
{
    //-----------------------------------------------------------------------------
    // Function      : MapOp
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 12:19:57 2014
    //-----------------------------------------------------------------------------
    ///
    /// Destination map reference into the operator
    ///
    /// @param map map to insert name to instance mapping
    ///
    MapOp(std::map<std::string, C *> &map)
      : map_(map)
    {}

    //-----------------------------------------------------------------------------
    // Function      : operator()
    // Purpose       :
    // Special Notes :
    // Scope         : public
    // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
    // Creation Date : Tue Apr 22 12:21:47 2014
    //-----------------------------------------------------------------------------
    ///
    /// 
    ///
    virtual bool operator()(DeviceInstance *device_instance)
    {
      map_[getName(device_instance)] = static_cast<C *>(device_instance);
      return true;
    }

    std::map<std::string, C *> &        map_;
};


//-----------------------------------------------------------------------------
// Function      : getDeviceInstances
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 11:02:16 2014
//-----------------------------------------------------------------------------
///
/// Calls forEachInstance() on d which iterates through all the instances copies them to the output iterator.
///
/// @param d D to call forEachInstance
/// @param it output iterator
///
template <class D, class Out>
void getDeviceInstances(const D &d, Out it) {
  forEachInstance(d, DeviceInstanceOutIteratorOp<Out>(it));
}


//-----------------------------------------------------------------------------
// Function      : getDeviceModels
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 11:02:16 2014
//-----------------------------------------------------------------------------
///
/// Calls forEachModel() on d which iterates through all the models copies them to the output iterator.
///
/// @param d D to call forEachModel
/// @param it output iterator
///
template <class D, class Out>
void getDeviceModels(const D &d, Out it)
{
  forEachModel(d, DeviceModelOutIteratorOp<Out>(it));
}


//-----------------------------------------------------------------------------
// Function      : getDeviceInstanceNames
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:22:16 2014
//-----------------------------------------------------------------------------
///
/// Calls forEachInstance() on d which iterates through all the instance and copies the name to the output iterator.
///
/// @param d D to call forEachInstance
/// @param it output iterator
///
template<class D, class Out>
void getDeviceInstanceNames(const D &d, Out it)
{
  forEachInstance(d, DeviceInstanceNameOp<Out>(it));
}


//-----------------------------------------------------------------------------
// Function      : mapDeviceInstances
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Tue Apr 22 12:22:22 2014
//-----------------------------------------------------------------------------
///
/// Calls forEachInstance() on d which iterates through all the instances and inserts the device instance name to pointer mapping.
///
/// @param d D to call forEachInstance
/// @param it output iterator
///
template<class D, class X>
void mapDeviceInstances(const D &d, std::map<std::string, X *> &map) {
  forEachInstance(d, MapOp<X>(map));
}

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Algorithm_h
