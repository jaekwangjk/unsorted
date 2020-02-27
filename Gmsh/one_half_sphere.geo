// Outer point

a1=0.025; 
a2=0.025;
d=0.2; 
L=0.8;
H=0.1; 

//Zone1
Point(1) = {-L/2, 0, 0, 1.0};
Point(2) = {-d/2-d/4, 0,0, 1.0};
Point(3) = {-d/2-d/4, H,0, 1.0};
Point(4) = {-L/2, H, 0,1.0};

Line(1) = {1, 2}; 
Line(2) = {2, 3};
Line(3) = {3, 4}; 
Line(4) = {4, 1}; 

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};


//Zone2
Point(5) = {-d/2 - a1, 0,0,1.0}; 
Point(6) = {-d/2 - a1/Sqrt(2), a1/Sqrt(2),0,1.0};
Point(7) = {-d/2 + a1/Sqrt(2), a1/Sqrt(2),0,1.0}; 
Point(8) = {-d/2 + a1, 0,0,1.0}; 
Point(9) = {-d/2 +d/4, H,0,1.0};
Point(10) = {-d/2 +d/4, 0,0,1.0}; 

Point(99) = {-d/2 , 0,0,1.0}; 

Line(5) = {2, 5}; 
Circle(6) = {5,99,6}; 
Circle(7) = {6,99,7}; 
Line(8) = {6,3}; 
Circle(9) = {7,99,8}; 
Line(10) = {8,10}; 
Line(11) = {10,9}; 
Line(12) = {9,3}; 
Line(13) = {7,9}; 


Curve Loop(2) = {5,6,8,-2};
Plane Surface(2) = {2};

Curve Loop(3) = {7,13,12,-8};
Plane Surface(3) = {3};

Curve Loop(4) = {9,10,11,-13};
Plane Surface(4) = {4};


//Zone3_OneSphereBall_ends_Here

Point(11) = {0.2 , 0,0,1.0}; 
Point(12) = {0.2 , H,0,1.0}; 
Line(14) = {10,11}; 
Line(15) = {11,12};
Line(16) = {12,9}; 

Curve Loop(5) = {14,15,16,-11};
Plane Surface(5) = {5};

//Tranfinte Curves
Transfinite Curve {1} = 4 Using Progression 1;
Transfinite Curve {2} = 4 Using Progression 1;
Transfinite Curve {3} = 4 Using Progression 1;
Transfinite Curve {4} = 4 Using Progression 1;
Transfinite Curve {5} = 4 Using Progression 1;
Transfinite Curve {6} = 4 Using Progression 1;
Transfinite Curve {7} = 4 Using Progression 1;
Transfinite Curve {8} = 4 Using Progression 1;
Transfinite Curve {9} = 4 Using Progression 1;
Transfinite Curve {10} = 4 Using Progression 1;
Transfinite Curve {11} = 4 Using Progression 1;
Transfinite Curve {12} = 4 Using Progression 1;
Transfinite Curve {13} = 4 Using Progression 1;
Transfinite Curve {14} = 4 Using Progression 1;
Transfinite Curve {15} = 4 Using Progression 1;
Transfinite Curve {16} = 4 Using Progression 1;



//Define boundary indicator in deal.II
Physical Line(10) = {6,7,9}; //sphere
Physical Line(12) = {4}; //inlet
Physical Line(13) = {15}; //outlet
Physical Line(14) = {3,12,16}; //wall
Physical Line(15) = {1,5,10,14}; //axis


// you need the physical surface, because that is what deal.II reads in
// designate plance surfaces 
Physical Surface(1) = {1,2,3,4,5};

Transfinite Surface {1};
Transfinite Surface {2}; 
Transfinite Surface {3}; 
Transfinite Surface {4};
Transfinite Surface {5};


// some parameters for the meshing:deal.ii specific 
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.09;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";

