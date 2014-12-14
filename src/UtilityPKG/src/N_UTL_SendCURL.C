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
// Filename       : $RCSfile: N_UTL_SendCURL.C,v $
//
// Purpose        :
//
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.6 $
//
// Revision Date  : $Date: 2014/08/30 04:35:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_SendCURL.h>

#ifdef Xyce_USE_CURL
#include <curl/curl.h>
#endif

namespace Xyce {
namespace Util {

namespace {
//-----------------------------------------------------------------------------
// Function      : write_callback
// Purpose       : Discard output from remote http responses
// Special Notes : 
// Scope         : file-local (unnamed namespace)
// Creator       : Dave Baur
// Creation Date : 8/27/14
//-----------------------------------------------------------------------------
size_t write_callback(char *ptr, size_t size, size_t nmemb, void *userdata)
{
  return size*nmemb;
}
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::sendTrackingData
// Purpose       : Send a packet of tracking info to a tracking URL
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 3/26/14
//-----------------------------------------------------------------------------
bool
sendTrackingData(
  const char *          url,
  const char *          proxy,
  const std::string &   message)
{
#ifdef Xyce_USE_CURL
  CURLcode result = CURLE_OK;

  if (url && url[0] != '\0') {
    curl_global_init(CURL_GLOBAL_ALL);

    CURL *curl = curl_easy_init();
    if (!curl)
      return false;

    struct curl_slist *slist=NULL;
    slist = curl_slist_append(slist,"Expect:");
    slist = curl_slist_append(slist,"Content-Type: application/json; charset=utf-8");

    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);                 // 1 to enabled verbose logging
    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);          // 1 to verify SSL peer

    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_URL, url);
    if (proxy && proxy[0] != '\0')
      result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_PROXY, proxy);

    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_TIMEOUT, 2L);
    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_POST, 1L);
    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_POSTFIELDS, message.c_str());
    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, message.size());
    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_HTTPHEADER, slist);
    // Set the write function to discard server response
    result = result != CURLE_OK ? result : curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    result = result != CURLE_OK ? result : curl_easy_perform(curl);

    curl_slist_free_all(slist);
    curl_easy_cleanup(curl);
  }

  return result == CURLE_OK;

#else
  return true;
#endif  // Xyce_USE_CURL
}

} // namespace Util
} // namespace Xyce
