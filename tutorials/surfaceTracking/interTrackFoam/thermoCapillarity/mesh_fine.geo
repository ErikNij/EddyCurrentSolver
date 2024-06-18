ps = 1.5;

d = 0.1;
dByTwo = d/2;

Point(1)  = {0   ,0   ,-dByTwo, ps};
Point(2)  = {0   ,0.15,-dByTwo, ps};
Point(3)  = {0   ,4.00,-dByTwo, ps};
Point(4)  = {4.00,4.00,-dByTwo, ps};
Point(5)  = {4.00,0   ,-dByTwo, ps};
Point(6)  = {0.15,0   ,-dByTwo, ps};

Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Circle(6) = {6, 1, 2};

Line Loop(7) = {2, 3, 4, 5, 6};
Plane Surface(7) = {7};

Extrude {0, 0, d} {
    Surface{7};
    Layers{1};
    Recombine;
}

Physical Volume("fluid") = {1};

Physical Surface("inner") = {33};
Physical Surface("side_x") = {25};
Physical Surface("side_y") = {21};
Physical Surface("symm_x") = {17};
Physical Surface("symm_y") = {29};
Physical Surface("empty_z") = {7, 34};

size = 0.01;

arc_size = 0.25;
arc_n = arc_size/size*2;

inner_size = 3.85;
inner_n = inner_size/size/5;
inner_e = 0.95;

outer_size = 4;
outer_n = outer_size/size/14;
outer_e = 1;

Transfinite Line {6, 13} = arc_n Using Progression 1;
Transfinite Line {5, 12} = inner_n Using Progression inner_e;
Transfinite Line {2, 9} = inner_n Using Progression 1/inner_e;

Transfinite Line {4, 11} = outer_n Using Progression outer_e;
Transfinite Line {3, 10} = outer_n Using Progression 1/outer_e;