cl1=1;
Point(2) = {877.6, 221.7, cl1};
Point(4) = {754.9, 276.3, cl1};
Point(5) = {584.9, 362.3, cl1};
Point(6) = {420.6, 400.8, cl1};
Point(7) = {318.6, 438.3, cl1};
Point(8) = {302.9, 407.6, cl1};
Point(9) = {333.1, 364.7, cl1};
Point(10) = {361.6, 327.6, cl1};
Point(11) = {383.9, 293.5, cl1};
Point(12) = {365.6, 263.6, cl1};
Point(13) = {432.9, 252, cl1};
Point(14) = {442.7, 217.7, cl1};
Point(15) = {511.4, 216.7, cl1};
Point(16) = {548.6, 240, cl1};
Point(17) = {569.9, 248.4, cl1};
Point(18) = {570.2, 231.1, cl1};
Point(19) = {566.3, 209, cl1};
Point(20) = {566.8, 141.8, cl1};
Point(21) = {591.4, 131, cl1};
Point(22) = {635.2, 139.6, cl1};
Point(23) = {622.6, 141.9, cl1};
Point(24) = {694.3, 136, cl1};
Point(25) = {640.8, 138.4, cl1};
Point(26) = {717.4, 124.2, cl1};
Point(27) = {763.1, 150.7, cl1};
Point(28) = {814.8, 178.1, cl1};
Point(29) = {764.8, 154.6, cl1};
Point(30) = {873.2, 226, cl1};
Point(31) = {748.9, 279.3, cl1};
Point(32) = {651.8, 325.7, cl1};
Point(33) = {581.8, 362.1, cl1};Line(1) = {2, 30};
Line(2) = {30, 4};
Line(3) = {4, 31};
Line(4) = {31, 32};
Line(5) = {32, 5};
Line(6) = {5, 33};
Line(7) = {33, 6};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 9};
Line(11) = {9, 10};
Line(12) = {10, 11};
Line(13) = {11, 12};
Line(14) = {12, 13};
Line(15) = {13, 14};
Line(16) = {14, 15};
Line(17) = {15, 16};
Line(18) = {16, 17};
Line(19) = {17, 18};
Line(20) = {18, 19};
Line(21) = {19, 20};
Line(22) = {20, 21};
Line(23) = {21, 23};
Line(24) = {23, 22};
Line(25) = {22, 25};
Line(26) = {25, 24};
Line(27) = {24, 26};
Line(28) = {26, 27};
Line(29) = {27, 29};
Line(30) = {29, 28};
Line(31) = {28, 2};
Line Loop(32) = {31, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
Plane Surface(33) = {32};
Physical Point(34) = {2, 30, 4, 31, 32, 5, 33, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 22, 25, 24, 26, 27, 29, 28};
Physical Line(35) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
Physical Surface(36) = {33};
Recombine Surface {33};
