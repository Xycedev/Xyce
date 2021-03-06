//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2003, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : Xyce.C
//
// Purpose        : front end for standalone Xyce executable
//
// Special Notes  :
//
// Creator        : Eric Rankin
//
// Creation Date  : 01/28/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2014/05/13 14:50:45 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------


#include <N_CIR_Xygra.h>
#include <N_ERH_ErrorMgr.h>
#include <iostream>
#include <fstream>

// Function to be called if memory runs out:
void _new_handler (void)
{
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, "OUT OF MEMORY (error in 'new')");
}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       : front end for standalone Xyce executable
// Special Notes :
// Scope         : 
// Creator       : Eric Rankin
// Creation Date : 01/28/2004
//-----------------------------------------------------------------------------
int main( int iargs, char *cargs[] )
{
  // Set divide by zero, and invalid operation handling on linux
  // Set out of memory detection on all systems
  set_new_handler (&_new_handler);
 
  N_CIR_Xygra * XycePtr = new N_CIR_Xygra();

  bool bsuccess = XycePtr->initialize(iargs, cargs);
  vector <string> deviceNames;
  vector <double> vN;

  if (bsuccess) 
    bsuccess = XycePtr->getDeviceNames("YXYGRA",deviceNames);

  ifstream infile("s-and-k.out",ifstream::in);

  if (bsuccess)
    bsuccess=infile.good();

  if (bsuccess)
  {
    vector<vector<int> >coilWindings;
    vector<vector<string> >coilNames;
    vector<int> numNodes;
    vector<int> numWindings;
    coilWindings.resize(deviceNames.size());
    coilNames.resize(deviceNames.size());
    numNodes.resize(deviceNames.size());
    numWindings.resize(deviceNames.size());
    for (int i=0; i < deviceNames.size(); ++i)
    {
      XycePtr->xygraGetCoilWindings(deviceNames[i],coilWindings[i]);
      XycePtr->xygraGetCoilNames(deviceNames[i],coilNames[i]);
      numNodes[i]=XycePtr->xygraGetNumNodes(deviceNames[i]);
      numWindings[i]=XycePtr->xygraGetNumWindings(deviceNames[i]);
      
      cout << " Xygra device " << deviceNames[i] << " has " 
           << coilWindings[i].size() << " coils " << endl;
      for (int j=0; j<coilWindings[i].size(); j++)
      {
        cout << "    coil " << j << " is named " << coilNames[i][j] << " and has " << coilWindings[i][j] 
             << "windings" << endl;
      }
      cout << "     for a total of " << numNodes[i] << " nodes and " << numWindings[i] << " windings." << endl;
    }

    if (deviceNames.size() != 1)
    {
      cerr << " Sorry, this test is designed to work only with one Xygra device. " << endl;
    }
    else
    {
      
      double completedTime, timeStep;
      completedTime = 0.0;
      bool opComplete = false;

      while (!(XycePtr->simulationComplete()) && bsuccess)
      {
        cout << "Simulation incomplete, completedTime is " << completedTime 
             <<"." << endl;


        // We will read each line of the input file sequentially to get
        // the target final time, the s vector and K matrix for that target
        // time.
        double targetTime;
        vector<double> sV;
        vector<vector<double> > kM;
        sV.resize(numWindings[0]); // we only have 1 device for sure
        kM.resize(numWindings[0]);

        infile >> targetTime;
        if (infile.eof()) break; // won't know this till after we try reading
/********
        for (int i=0; i<numWindings[0]; i++)
          infile >> sV[i];
*/
        for (int i=0; i<numWindings[0]; i++)
        {
          kM[i].resize(numWindings[0]);
          for (int j=0; j<numWindings[0]; j++)
            infile >> kM[i][j];
        }
        for (int i=0; i<numWindings[0]; i++)
          infile >> sV[i];

        // Dump the vector and matrices:

        cout << "S vector for time " << targetTime << ":" << endl;
        for (int i=0; i<numWindings[0]; i++)
          cout << sV[i] << " " ;
        cout << endl;

        cout << "K matrix for time " << targetTime << ":" << endl;
        cout << "--------" << endl;
        for (int i=0; i<numWindings[0]; i++)
        {
          for (int j=0; j<numWindings[0]; j++)
            cout << kM[i][j] << " " ;
          cout << endl; 
        }
        cout << "--------" << endl;

        // By leaving off time argument, we don't interpolate.
        XycePtr->xygraSetSources(deviceNames[0],sV);
        XycePtr->xygraSetK(deviceNames[0],kM);

        if (targetTime > 0.0)
        {
          cout << "Calling simulateUntil with requested time " << targetTime << endl;

          bsuccess = XycePtr->simulateUntil(targetTime,completedTime);
          cout << "Simulated to " << completedTime << endl;
          for (int i=0; i<deviceNames.size(); i++)
          {
            int offset=0;
            XycePtr->xygraGetVoltages(deviceNames[i], vN);
            cout << " Nodal voltages for device " << deviceNames[i] << endl;
            for (int coil=0; coil<coilWindings[i].size(); coil++)
            {
              cout << "   Coil " << coil << ":" << endl;
              for (int node=0; node<coilWindings[i][coil]+1;node++)
              {
                cout << "    node " << node << " voltage " << vN[offset++]
                     << endl;
              }
            }
          }
        }
      }
    }
  }
  delete XycePtr;

  (bsuccess) ? exit(0) : exit(-1);
}

