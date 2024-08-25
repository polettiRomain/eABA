This code extends the Articulated-Body Algorithm implemented in OpenFOAM v2206. .
The complete information are given in ...

**How to use it? **

1. Copy the given rigidBodyDynamics and rigidBodyMeshMotion in your local src folder
2. Compile with wmake 
3. Include these libraries in the controlDict of your own test case:
   libs           ("liboverset.so" "libMyrigidBodyDynamics2.so" "libMyrigidBodyMeshMotion2.so");
