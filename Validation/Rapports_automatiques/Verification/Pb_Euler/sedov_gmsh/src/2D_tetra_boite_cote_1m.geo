
eps = 0.007;

Point(1) = {-0.5, -0.5, 0 ,eps};
//+
Point(2) = {-0.5, 0.5, 0, eps};
//+
Point(3) = {0.5, 0.5, 0, eps};
//+
Point(4) = {0.5, -0.5, 0, eps};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Surface(1) = {1};
//+
Physical Curve("left", 5) = {1};
//+
Physical Curve("up", 6) = {2};
//+
Physical Curve("right", 7) = {3};
//+
Physical Curve("down", 8) = {4};
//+
Physical Surface("domain", 9) = {1};
//+
Transfinite Surface {1} = {1, 2, 3, 4}; //structure
//+
Recombine Surface {1} = 1;


