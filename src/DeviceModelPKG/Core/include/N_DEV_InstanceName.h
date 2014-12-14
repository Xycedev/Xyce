//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_InstanceName.h,v $
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/05/19 20:00:59 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_InstanceName_h
#define Xyce_N_DEV_InstanceName_h

#include <string>

#include <N_DEV_fwd.h>
#include <N_UTL_NoCase.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : InstanceName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon May 12 15:53:34 2014
//-----------------------------------------------------------------------------
///
/// Devices and models are each named.  Models are not encoded, so
/// simple string representation is sufficient.  However Devices are a
/// different lot.  They are encoded as
///
/// [s:]*xname
/// [s:]*Ytype!name
/// [s:]*Utype!name!count
///
/// where s is a subciruit name, x is a single letter device type, type
/// is a multiletter device type (no Y or U prefix) and count is a
/// special input count for the U device.
///
/// Currently encoded names are accepted and then decoded using the
/// getter's.  In the future, these will be stored in component form and
/// then encoded onyl as needed.
///
class InstanceName
{
private:
  InstanceName(char device_letter, const std::string &device_type, const std::string &subcircuit_name, const std::string &device_name, int num_inputs)
    : deviceLetter_(device_letter),
      deviceType_(device_type),
      subcircuitName_(subcircuit_name),
      deviceName_(device_name),
      numInputs_(num_inputs),
      name_()
  {
    encode();
  }

public:
  //-----------------------------------------------------------------------------
  // Function      : InstanceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 15:57:34 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates an empty entity name.
  ///
  ///
  InstanceName()
    : deviceLetter_(0),
      deviceType_(),
      subcircuitName_(),
      deviceName_(),
      numInputs_(0),
      name_()
  {}

  //-----------------------------------------------------------------------------
  // Function      : InstanceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 15:58:19 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates an entity
  ///
  /// @invariant
  ///
  /// @param name model or encoded device entity name
  ///
  /// @return 
  ///
  ///
  explicit InstanceName(const std::string &name)
    : deviceLetter_(0),
      deviceType_(),
      deviceName_(),
      subcircuitName_(),
      numInputs_(0),
      name_(name)
  {
    decode();
  }

  //-----------------------------------------------------------------------------
  // Function      : InstanceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 15:58:19 2014
  //-----------------------------------------------------------------------------
  ///
  /// Creates an entity
  ///
  /// @invariant
  ///
  /// @param name model or encoded device entity name
  ///
  /// @return 
  ///
  ///
  InstanceName &operator=(const InstanceName &entity_name) {
    deviceLetter_ = entity_name.deviceLetter_;
    deviceType_ = entity_name.deviceType_;
    deviceName_ = entity_name.deviceName_;
    subcircuitName_ = entity_name.subcircuitName_;
    numInputs_ = entity_name.numInputs_;
    name_ = entity_name.name_;
    return *this;
  }

  //-----------------------------------------------------------------------------
  // Function      : getDeviceLetter
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 16:19:04 2014
  //-----------------------------------------------------------------------------
  ///
  /// Return the first letter of the device specification after the
  /// subcircuit.
  ///
  /// The device letter is Y or U for multiletter device names or the
  /// letter that identifies the device for sinlge letter names.
  ///
  /// @return single letter associated with the device
  ///
  char getDeviceLetter() const {
    return deviceLetter_;
  }

  //-----------------------------------------------------------------------------
  // Function      : getDeviceType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 16:02:24 2014
  //-----------------------------------------------------------------------------
  ///
  /// Decodes the device type.
  ///
  /// The device type is a string containing the first letter for a
  /// single letter device name, or the string starting after the Y or
  /// U up to but noe including the exclamation point (!).
  ///
  /// Subcircuit prefixes are not included.
  ///
  /// @return string representing the device type
  ///
  ///
  std::string getDeviceType() const {
    return deviceType_;
  }


  //-----------------------------------------------------------------------------
  // Function      : getDeviceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 16:04:37 2014
  //-----------------------------------------------------------------------------
  ///
  /// Decodes the device name. 
  ///
  /// The device name is a string containing the first letter for a
  /// single letter device name, or the string starting after the Y or
  /// U up to but noe including the exclamation point (!).
  ///
  /// Subcircuit prefixes are not included.
  ///
  /// @return string representing the device name
  ///
  ///
  std::string getDeviceName() const {
    return deviceName_;
  }


  //-----------------------------------------------------------------------------
  // Function      : getDeviceName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 16:04:37 2014
  //-----------------------------------------------------------------------------
  ///
  /// Decodes the device name. 
  ///
  /// The device name is a string containing the first letter for a
  /// single letter device name, or the string starting after the Y or
  /// U up to but noe including the exclamation point (!).
  ///
  /// @return string representing the device name
  ///
  ///
  std::string getSubcircuitName() const {
    return subcircuitName_;
  }


  //-----------------------------------------------------------------------------
  // Function      : getNumInputs
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 12 16:23:42 2014
  //-----------------------------------------------------------------------------
  ///
  /// For the U device, return the number of inputs which have been
  /// encoded into the device name.  (These probably should be encoded
  /// into the device parameters, not the name.)
  ///
  /// @invariant
  ///
  ///
  /// @return 
  ///
  ///
  int getNumInputs() const {
    return numInputs_;
  }

  //-----------------------------------------------------------------------------
  // Function      : getEncodedName
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Mon May 19 07:40:55 2014
  //-----------------------------------------------------------------------------
  ///
  /// Return the instance name encoded as:
  ///   [s:]*xname
  ///   [s:]*Ytype!name
  ///   [s:]*Utype!name!count
  ///
  /// @return encoded instance name
  ///
  ///
  const std::string &getEncodedName() const
  {
    return name_;
  }

private:
  void decode();
  void encode();

  char decodeDeviceLetter() const;
  std::string decodeDeviceType() const;
  std::string decodeDeviceName() const;
  std::string decodeSubcircuitName() const;
  int decodeNumInputs() const;

private:
  char                deviceLetter_;          ///< Device letter (Y or U included)
  std::string         deviceType_;            ///< Device type from netlist (does NOT include Y or U prefix)
  std::string         deviceName_;            ///< Device name from netlist (subcircuit NOT included)
  std::string         subcircuitName_;        ///< Device subcircuit from netlist
  int                 numInputs_;             ///< Hack for U type device (needs to become a parameter)

  std::string         name_;                  ///< Complete encoded name
};

inline bool operator==(const InstanceName &entity_name, const std::string &name) {
  return equal_nocase(entity_name.getEncodedName(), name);
}

inline bool operator<(const InstanceName &entity_name, const std::string &name) {
  return less_nocase(entity_name.getEncodedName(), name);
}

inline bool operator!=(const InstanceName &entity_name, const std::string &name) {
  return !equal_nocase(entity_name.getEncodedName(),  name);
}

std::string setupOutputName(const InstanceName &name);

std::ostream &operator<<(std::ostream &os, const InstanceName &entity_name);

void spiceInternalName(std::string &entity_name);
std::string spiceInternalName(const InstanceName &entity_name, const std::string &lead);
std::string spiceStoreName(const InstanceName &entity_name, const std::string &lead);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_InstanceName_h

