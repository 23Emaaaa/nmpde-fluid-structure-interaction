// Gmsh project created on Mon Nov 24 23:29:58 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.03};
//+
Point(2) = {1, 0, 0, 0.03};
//+
Point(3) = {0, 1, 0, 0.03};
//+
Point(4) = {1, 1, 0, 0.03};
//+
Point(5) = {0.3, 0.3, 0, 0.03};
//+
Point(6) = {0.7, 0.3, 0, 0.03};
//+
Point(7) = {0.3, 0.7, 0, 0.03};
//+
Point(8) = {0.7, 0.7, 0, 0.03};
//+
Point(9) = {0, 0.3, 0, 0.03};
//+
Point(10) = {1, 0.3, 0, 0.03};
//+
Point(11) = {0.5, 1, 0, 0.03};
//+
Line(1) = {1, 2};
//+
Line(2) = {1, 9};
//+
Line(3) = {9, 5};
//+
Line(4) = {5, 7};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 6};
//+
Line(7) = {6, 10};
//+
Line(8) = {10, 2};
//+
Line(9) = {10, 4};
//+
Line(10) = {3, 11};
//+
Line(11) = {11, 4};
//+
Line(12) = {3, 9};
//+
Curve Loop(1) = {10, 11, -9, -7, -6, -5, -4, -3, -12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 5, 6, 7, 8, -1, 2, 3};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Fluid", 12) = {1};
//+
Physical Surface("Solid", 13) = {2};
//+
Physical Curve("Inflow", 14) = {10};
//+
Physical Curve("Outflow", 15) = {11};
//+
Physical Curve("BottomWall", 16) = {1};
//+
Physical Curve("LeftWallSolid", 17) = {2};
//+
Physical Curve("LeftWallFluid", 18) = {12};
//+
Physical Curve("RightWallSolid", 19) = {8};
//+
Physical Curve("RightWallFluid", 20) = {9};
//+




