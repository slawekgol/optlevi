Include "coil.geo";
Include "stale.geo";

// calculating centers of each turn in the coil

// startind = 3;
// number_of_turns = 8;
// number_of_turns_inactive = 6;


// indx() = {};
// indy() = {};

// // the round part of the coil gutter
// For i In {0:number_of_turns-1}
//   x = coil_x + r * Cos(turn_start + turn_dist * (i-(number_of_turns-1)/2));
//   y = coil_y + r * Sin(turn_start + turn_dist * (i-(number_of_turns-1)/2));
//   If (i == 0)
//     lx = x;
//     ly = y;
//   EndIf
//   If (i == number_of_turns - 1)
//     rx = x;
//     ry = y;
//   EndIf
//   indx() += x;
//   indy() += y;
// EndFor

// // the left (closer to the axis) wall of the coil cutter
// number_of_turns_l = 0; // the number of turns
// number_of_turns_lstep = 0.007; // distance between the centers of the turns

// For i In {1:number_of_turns_l}
//   x = lx;
//   y = ly + i * number_of_turns_lstep;
//   indx() = {x, indx()};
//   indy() = {y, indy()};
// EndFor

// // the right (further form the axis) wall of the coil cutter

// number_of_turns_r=0;
// number_of_turns_rstep=0.014;

// For i In {1:number_of_turns_r}
//   x = rx;
//   y = ry + i * number_of_turns_rstep;
//   indx() += x;
//   indy() += y;
// EndFor

// number_of_turns = #indx();

rind = 0.003; // radius of conductor (turn intersection)



chy = 0.000;  //  y-coordinate of charge
chx = charge_R; //  left x-coordinate of charge

rint =  20 * chx; // internal radius of infinite region
rext = 1.2 * rint; // external radius of infinite region

//Printf("number_of_turns = %g;", number_of_turns) > "param.txt";
//Printf("rint = %g;", rint) >> "param.txt";
//Printf("rext = %g;", rext) >> "param.txt";

lcind = rind / 4; // mesh size in the coil
lcair = 0.05; // mesh size in the air

lcch = 0.0005;


p = 100;
indl() = {};

ind_s=1000;
ind_in()={};
For i In {1:#indx[]}
  Point(p) = {indx(i-1), indy(i-1), 0, lcind};
  Point(p + 1) = {indx(i - 1) - rind, indy(i - 1), 0, lcind};
  Point(p + 2) = {indx(i - 1), indy(i - 1) - rind, 0, lcind};
  Point(p + 3) = {indx(i - 1) + rind, indy(i - 1), 0, lcind};
  Point(p + 4) = {indx(i - 1), indy(i - 1) + rind, 0, lcind};
  Circle(p) = {p + 1, p, p + 2};
  Circle(p + 1) = {p + 2, p, p + 3};
  Circle(p + 2) = {p + 3, p, p + 4};
  Circle(p + 3) = {p + 4, p, p + 1};
  Curve Loop(p) = {p : p + 3};
  indl() += {p};

  Plane Surface(ind_s+i) = {p};
  Physical Surface(Str(Sprintf("Ind_in %g", i)), IND_IN + i - 1) = {ind_s+i};
  p += 10;
EndFor



Point(1) = {0, 0, 0, lcair};
Point(2) = {0, -rint, 0, lcair};
Point(3) = {rint, 0, 0, lcair};
Point(4) = {0, rint, 0, lcair};
Point(5) = {0, -rext, 0, lcair};
//Point(6) = {rext, 0, 0, lcair};
//Point(7) = {0, rext, 0, lcair};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
//Circle(3) = {5, 1, 6};
//Circle(4) = {6, 1, 7};
//Line(5) = {5, 2};
Line(6) = {2, 1};
Line(7) = {1, 4};
//Line(8) = {4, 7};

If (distribution==0)
  Point(9) ={chx,chy,0,lcch};
  Point(10) = {chx - charge_r  , chy  , 0, lcch};
  Point(11) = {chx, chy + charge_r , 0, lcch};
  Point(12) = {chx + charge_r  , chy , 0, lcch};
  Point(13) = {chx, chy - charge_r , 0, lcch};
  Circle(10) = {10, 9, 11};
  Circle(11) = {11, 9, 12};
  Circle(12) = {12, 9, 13};
  Circle(13) = {13, 9, 10};
  Curve Loop(2) = {10, 11, 12, 13};
  Plane Surface(1) = {2};
  chargePlane=1;
  
EndIf

If(distribution>0)
  Point(9) ={chx-area_r,chy-area_y,0,lcch};
  Point(10) = {chx - area_r  , chy+area_y  , 0, lcch};
  Line(11) = {9,10};
  
  
  tmp[]=Extrude{2*area_r,0,0} {Curve{11};Layers{(2*area_r/lcch)};Recombine;};
  Printf("tmp = %g %g %g %g", tmp[0], -tmp[2],-11,tmp[3]);
  Curve Loop(2) = {tmp[0], -tmp[2],-11, -tmp[3]};
  chargePlane = tmp[1];
EndIf

Curve Loop(1) = {1, 2, -7, -6};

//Curve Loop(3) = {2, 8, -4, -3, 5, 1};


Plane Surface(2) = {1, 2, indl()};
//Plane Surface(3) = {3};

Physical Surface("Charge_in", CHARGE_IN) = chargePlane;
Physical Surface("Air_in",  AIR_IN) = 2;
//Physical Surface("Air inf_in", 3) = 3;

Physical Curve("Axis", AXIS) = {6, 7};
Physical Curve("Inf", INF) =  {1,2};

