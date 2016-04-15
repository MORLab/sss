The second order system

   M x"(t) + D x'(t) + K x(t) = B u(t)
                         y(t) = C x(t)
       
is transformed to the first order system

   E z'(t) = A z(t) + G u(t)
      y(t) = L z(t)
  
where

      | D  M|
   E= | M  0|
   
      |-K  0|
   A= | 0  M|
   
      | B |
   G= | 0 |
   
   L= [C  0]
   
         | x(t)  |
   z(t)= | x'(t) | .
   
Matrices M, D, K are assumed to be quadratic, symmetric and positive definit.

The fieldnames have to end with _  to indicate that the data 
are inputdata for the algorithm.
eqn.M_
eqn.K_
eqn.D_
eqn.B

[Benner, Kuerschner, Saak: An improved numerical method for balanced truncation for symmetric second-order systems, 2013]
