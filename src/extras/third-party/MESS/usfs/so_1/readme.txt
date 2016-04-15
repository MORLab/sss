Function Handles for Second Order LQR Problem!
Linear Beam 1D Example from 
http://portal.uni-freiburg.de/imteksimulation/downloads/benchmark.
Second Order System
Mx" + D x' + K x 	= B u
y 	= C x
Transformed to 
First Order System
|-K 0||x' |= |0 -K||x | + |0|
|0  M||x''|  |-K -D||x'|   |B|u
   E   x'  =  Ax        +  B
Attention the Matrix M D K are symmetric and quadratic.
K is a fullrank Matrix.
The fieldnames have to end with _  to indicate that the Data 
are inputdata for the Algorithm.
eqn.M_
eqn.K_
eqn.D_
eqn.B