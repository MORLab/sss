%% Functions Reference
% UI-functions of the sss toolbox sorted by categories

%% Sparse State-Space
%
% * <.\appendhelp.html |append|>  --  Appends a set of sparse LTI systems (sss)
% * <.\bodehelp.html |bode|>  --  Plots the bode diagram of an LTI system
% * <.\bodemaghelp.html |bodemag|>  --  Bode magnitude plot of an LTI system
% * <.\bodeplothelp.html |bodeplot|>  --  Bode plot of an LTI system
% * <.\c2dhelp.html |c2d|>  --  Converts a sss object from continues to discrete
% * <.\clearhelp.html |clear|>  --  Deletes all state-space dimension related properties
% * <.\connecthelp.html |connect|>  --  Block diagram interconnections of dynamic systems sparse LTI system (sss)
% * <.\connectSsshelp.html |connectSss|>  --  Connects an appended sparse state space  LTI system (sss) with feedback matrix K
% * <.\dcgainhelp.html |dcgain|>  --  Computes the dcgain of a LTI-system
% * <.\decayTimehelp.html |decayTime|>  --  Computes the time period in which a sparse LTI system levels off
% * <.\diaghelp.html |diag|>  --  Transforms an LTI system to (block)-diagonal representation
% * <.\disphelp.html |disp|>  --  Displays information about a sparse state-space model
% * <.\eighelp.html |eig|>  --  Compute eigenvalues and eigenvectors of a sparse state-space model
% * <.\eigshelp.html |eigs|>  --  Compute eigenvalues of the sparse state space system using sparse matrices.
% * <.\frdhelp.html |frd|>  --  Convert to frequency-response data model
% * <.\freqresphelp.html |freqresp|>  --  Frequency response of sparse state-space systems.
% * <.\getfrdhelp.html |getfrd|>  --  Get frd-object(s) for frequency-response functions
% * <.\impulsehelp.html |impulse|>  --  Computes and/or plots the impulse response of a sparse LTI system
% * <.\issdhelp.html |issd|>  --  Check strict dissipativity of sparse LTI system
% * <.\isstablehelp.html |isstable|>  --  Check stability of sparse LTI system
% * <.\lsimhelp.html |lsim|>  --  Simulate time response of sparse dynamic system to arbitrary inputs
% * <.\lyapcholhelp.html |lyapchol|>  --  Solve Lyapunov equations
% * <.\minushelp.html |minus|>  --  Computes difference of two sparse LTI systems.
% * <.\mtimeshelp.html |mtimes|>  --  Computes the product of two LTI systems.
% * <.\normhelp.html |norm|>  --  Computes the p-norm of an sss LTI system
% * <.\polehelp.html |pole|>  --  Compute largest poles of an LTI system
% * <.\plushelp.html |plus|>  --  Computes sum of two sparse LTI systems.
% * <.\pzmaphelp.html |pzmap|>  --  Pole-zero plot of sparse state-space system
% * <.\residuehelp.html |residue|>  --  Computes residues, poles and feedthrough of an LTI system
% * <.\sigmahelp.html |sigma|>  --  Plots the singular values of the frequency response of a sparse LTI system
% * <.\simhelp.html |sim|>  --  Simulates a sss system using iddata input time series
% * <.\sizehelp.html |size|>  --  Computes the size of a sparse LTI system (sss)
% * <.\spyhelp.html |spy|>  --  Plot sparsity pattern of sss system
% * <.\sshelp.html |ss|>  --  Converts sparse LTI system (sss) to MATLAB built-in ss object
% * <.\ssshelp.html |sss|>  --  Create sparse state-space (sss) model, convert to sss model
% * <.\stephelp.html |step|>  --  Computes and/or plots the step response of a sparse LTI system
% * <.\truncatehelp.html |truncate|>  --  Truncates a sparse LTI system (sss)
% * <.\zerohelp.html |zero|>  --  Compute largest invariant zeros of an LTI system
% * <.\zpkhelp.html |zpk|>  --  Compute largest poles and zeros or zpk object of an LTI system

%% Extras
%
% * <.\bodeToFighelp.html |bodeToFig|>  --  Transforms MATLAB bode/step/impulse plot into normal figure
% * <.\ispdhelp.html |ispd|>  --  Determines if a matrix is positive definite
% * <.\iterativeRefinementhelp.html |iterativeRefinement|>  --  Use iterative LSE solver to refine LSE solution
% * <.\loadSsshelp.html |loadSss|>  --  Creates an sss-object from .mat file data
% * <.\parseOptshelp.html |parseOpts|>  --  Creates an optional struct after parsing
% * <.\partitionhelp.html |partition|>  --  Partitions a matrix into submatrices
% * <.\second2firsthelp.html |second2first|>  --  Convert a 2nd order system to state space representation
% * <.\solveLsehelp.html |solveLse|>  --  Solve linear system of equations

%% Sim
%
% * <.\simBackwardEulerhelp.html |simBackwardEuler|>  --  Integrates sss model using backward (implicit) Euler
% * <.\simDiscretehelp.html |simDiscrete|>  --  Integrates discrete time model
% * <.\simForwardEulerhelp.html |simForwardEuler|>  --  Integrates sss model using forward (explicit) Euler
% * <.\simInithelp.html |simInit|>  --  Initialization for simulation functions
% * <.\simRK4help.html |simRK4|>  --  Integrates sss model using Runge-Kutta 4
% * <.\simUpdatehelp.html |simUpdate|>  --  Update for simulation functions

%% Third Party
%
% * <.\MESShelp.html |M.E.S.S.|>  --  Solves large sparse matrix equations

%% Demos
%
% * <.\sss_gettingStartedhelp.html |sss_gettingStarted|>  --  Introductory demo to sss toolbox
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
