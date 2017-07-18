%% M.E.S.S.
% Solves large sparse matrix equations
%
%% Description
% M.E.S.S. is the successor to the LYAPACK Toolbox for MATLAB. It is 
% intended for solving large sparse matrix equations. The new version has 
% been rewriten in large parts to fit the drastic upgrades in the Matlab 
% releases since 2000. Additionally new solvers for differential Riccati 
% equations extend the functionality and many enhancements upgrade the 
% efficiency and runtime behaviour enlarging the number of unknowns that 
% can now be computed. 
%
% Amongst other things, M.E.S.S. can solve Lyapunov and Riccati equations, 
% and perform model reduction of systems in state space and structured 
% differential algebraic form. M.E.S.S. has been implemented in MATLAB as 
% well as C, with bindings also via MEX and to Python. M.E.S.S. is 
% therefore not restricted to the solution of "academic toy problems". 
% Several measures have been taken to enhance the computational performance
% of MESS routines in both implementations. To put this into the right 
% perspective, Lyapunov equations of order 20 000 were solved by MESS 
% within less than a minute on a regular laptop computer. On a 64bit 
% computeserver algebraic Riccati equations of order 250 000 can be solved 
% in well below an hour and solutions to Lyapunov equations for 3d 
% multiphysics applications with roughly 500 000 DOFs have been computed in
% only a few hours. When using standard (dense) methods, supercomputers are
% needed to solve problems of this size in reasonable time.
%
% For more toolbox info, have a look at the 
% <http://www.mpi-magdeburg.mpg.de/projects/mess>
%
%% References
% * *[1] Saak, J. and Koehler, M. and Benner, P. (2016)*, M-M.E.S.S.-1.0 -- 
% The Matrix Equations Sparse Solvers library
% * *[2] Penzl (2000)*, LYAPACK - A MATLAB Toolbox for Large Lyapunov
% and Riccati Equations, Model Reduction Problems, and
% Linear-Quadratic Optimal Control Problems.
%

%%
% <html>
%   <hr>
%   <p class="copy">&copy; 2015-2017 RT Technische Universit&auml;t M&uuml;nchen
%        <tt class="minicdot">&#149;</tt>
%        <a href="https://www.rt.mw.tum.de/?sss">Website</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:txts/LICENSE.txt">License</a>
%        <tt class="minicdot">&#149;</tt>
%        <a href="file:txts/RELEASE.txt">Release Notes</a>
%   </p>
% <div>
% <table>
%  <tr>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%             <img src="img/logo_sss_long.png" alt="sss_Logo" height="40px">
%      </td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/Logo_Textzusatz_rechts_engl_Chair.png" alt="RT_Logo" height="40px"></td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/logoTFD_SR_RGB.png" alt="TFD_Logo" height="40px"></td>
%   <td style="background-color:#ffffff; border:0; width:25%; vertical-align:middle; text-align:center">
%      <img src="img/TUM-logo.png" alt="TUM_Logo" height="40px"></td>
%  </tr>
% </table>
% </div>
% </html>
