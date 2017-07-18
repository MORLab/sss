%% What is sss?
%
% *sss* is a sparse state-space and system analysis toolbox for
% MATLAB designed to work with large-scale dynamical systems in
% state-space.
%
% *sss* extends the capabilities of the Control System Toolbox by defining
% sparse state-space (sss) objects and implementing large-scale and
% sparsity exploiting versions of common system and control functions (such
% as bode, step, pzmap, eig,...).
%
% By using *sss*, it is possible to define and analyze dynamic
% system objects with state-space dimensions higher than
% $\mathcal{O}(10^4)$, which is generally the limit for standard built-in
% ss and dss objects.
%
% *sss* is a MATLAB toolbox developed at the Model Order Reduction Lab (MORLAB)
% of the Chair of Automatic Control in collaboration with the Chair of 
% Thermofluid Dynamics, TU M&uuml;nchen.
%
% The *sss* toolbox is also comprised in the <https://www.rt.mw.tum.de/?sssMOR
% sssMOR toolbox>, which in addition to the sss functions further contains classic 
% and state-of-the-art model order reduction (MOR) techniques to capture the 
% dynamics of large-scale systems in reduced order models of significantly smaller size.
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
