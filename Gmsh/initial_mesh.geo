//GMSH file basic. Deal.ii understand mesh generated from this file. 
//You need to build on this
// Outer point
Point(3) = {0, 0, 0, 1.0};
Point(4) = {1, 0, 0, 1.0};
Point(5) = {1, 1, 0, 1.0};
Point(6) = {0, 1, 0, 1.0};
// Cylinder point
Point(7) = {0.7,0.3,0, 1.0}; 
Point(8) = {0.7,0.7,0, 1.0}; 
Point(9) = {0.3,0.7,0, 1.0};
Point(10) = {0.3,0.3,0, 1.0};
// Center of a Circle 
Point(11) = {0.5,0.5,0, 1.0};

//Construct block lines  
Line(1) = {3, 4}; Transfinite Curve {1} = 4 Using Progression 1;
Line(2) = {4, 5}; Transfinite Curve {2} = 4 Using Progression 1;
Line(3) = {5, 6}; Transfinite Curve {3} = 4 Using Progression 1;
Line(4) = {6, 3}; Transfinite Curve {4} = 4 Using Progression 1;
// Construct Circle 
Circle(5) = {7,11,10}; Transfinite Curve {5} = 4 Using Progression 1;
Circle(6) = {8,11,7}; Transfinite Curve {6} = 4 Using Progression 1;
Circle(7) = {9,11,8}; Transfinite Curve {7} = 4 Using Progression 1;
Circle(8) = {10,11,9}; Transfinite Curve {8} = 4 Using Progression 1;
// Link Between Circle and block 
Line(9) =  {3, 10}; Transfinite Curve {9} = 4 Using Progression 1;
Line(10) = {9, 6}; Transfinite Curve {10} = 4 Using Progression 1; 
Line(11) = {5, 8}; Transfinite Curve {11} = 4 Using Progression 1;
Line(12) = {7, 4}; Transfinite Curve {12} = 4 Using Progression 1;



// Boundary Condition : You don't need it  
//Physical Curve("cylinder") = {7, 8, 5, 6};
//Physical Curve("inlet") = {1};
//Physical Curve("outlet") = {3};
//Physical Curve("top") = {2};
//Physical Curve("bottom") = {4};

//Define boundary indicator in deal.II

Physical Line(0) = {1, 2, 3, 4};
Physical Line(1) = {5, 8, 7, 6};



// Surface Selection 
Curve Loop(1) = {5, -9, 1, -12};
Plane Surface(1) = {1};
Curve Loop(2) = {6, 12, 2, 11};
Plane Surface(2) = {2};
Curve Loop(3) = {7, -11, 3, -10};
Plane Surface(3) = {3};
Curve Loop(4) = {8, 10, 4, 9};
Plane Surface(4) = {4};


// you need the physical surface, because that is what deal.II reads in
// designate plance surfaces 
Physical Surface(17) = {1,2,3,4};

Transfinite Surface {1}; 
Transfinite Surface {2}; 
Transfinite Surface {3}; 
Transfinite Surface {4};



// some parameters for the meshing: deal.ii specific 
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.09;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";

