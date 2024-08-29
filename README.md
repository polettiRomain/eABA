## Basics
This library extends the Articulated-Body Algorithm implemented in OpenFOAM v2206. 
It allows to enforce the kinematics of some joints of a multibody system while the rest is driven by the Newton-Euler equations. 
Complete information are found in the paper: 'insert doi'

## Folder description
This GitHub contains the library and a test case. The test case is a body and wing drone for which an extensive description is provided in the paper. 

## Installation
1. Copy the given rigidBodyDynamics and rigidBodyMeshMotion in your local src folder
2. Compile with wmake 
3. Include these libraries in the controlDict of your own test case:
   libs           ("liboverset.so" "libMyrigidBodyDynamics2.so" "libMyrigidBodyMeshMotion2.so");

## Running
1. Load locally OpenFOAM v2206
2. Run the with the Allrun script in the TEST_CASE folder 
