classdef testSim < sssTest
% testSim - Testing of sim.m
%
% Description:
%   The function sim.m is tested on:
%    + Simulation of an artifitialy constructed, stable SSS, SISO benchmark
%    % + Simulation of a SSS, SISO benchmark system (building).
%    % + Simulation of a SSS, SISO benchmark system (heat-cond).
%    % + Simulation of a DSSS, SISO benchmark system (SpiralInductorPeec).
%    % + Simulation of a SSS, MIMO benchmark system (iss).
%    % + Simulation of a SSS, MIMO benchmark system (CDplayer).
%    % + Simulation of a DSSS, MIMO benchmark system (rail_1357).
%    % + Verification of the obtained results.
%
%------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> morlab@rt.mw.tum.de <-
%------------------------------------------------------------------
% Authors:      Maria Cruz Varona
% Last Change:  05 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    methods(Test)
        
        function testSISOartificial(testCase)
           
            % create a artificial systems which has all eigenvalues on the 
            % left side of the imaginary axis
            lambda  = ones(100,1) + rand(100,1)*5;
            A       = diag(-lambda);
            B       = ones(100,1);
            C       = ones(1,100);
            sys     = sss(A,B,C);
            sys_ss  = ss(sys);
            
            % create an input signal u as a iddata-object
            Ts      = 1e-4;
            t       = 0:Ts:10;
            u       = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
            datau   = iddata([],u',Ts);
            
            % simulate with the sss-toolbox functions
            dataSparseRK4       = sim(sys,          datau,'RK4');
            dataSparseFoEuler   = sim(sys,          datau,'forwardEuler');
            dataSparseBaEuler   = sim(sys,          datau,'backwardEuler');
            dataSparseDiscrete  = sim(c2d(sys,Ts),  datau,'discrete');
            
            % simulate with the ss-toolbox function (for verification)
            actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
            expSim = lsim(sys_ss,u,t);
            
%             % plot the result
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
%             legend('RK4','ForwEuler','BackEuler','Discrete','lsim');
            
            % check if all solutions are equal
            for iCase=1:4
                verification(testCase, actSim(:,iCase), expSim);
            end
            
            
        end
        
        function testLsim(testCase)
            sys = sss('building');
            
            u       = [0:.1:100].';
            ts      = 0.01; %s
            tVec    = 0:ts:ts*(length(u)-1);

            % sss
            [y1,t1,x1]      = lsim(sys,u,ts);   %call with ts
            [y2,t2,x2]      = lsim(sys,u,tVec); %call with tVec
            % ss
            [yM, tM, xM]    = lsim(ss(sys),u,tVec);

            % check sizes
            verifySize(testCase,y1,size(yM));
            verifySize(testCase,y2,size(yM));
            verifySize(testCase,t1,size(tM));
            verifySize(testCase,t2,size(tM));
            verifySize(testCase,x1,size(xM));
            verifySize(testCase,x2,size(xM));
            
            % check time vector limits
            tol = 1e-10;
            verifyEqual(testCase,t1(1),     tM(1),  'AbsTol',tol);
            verifyEqual(testCase,t2(1),     tM(1),  'AbsTol',tol);
            verifyEqual(testCase,t1(end),   tM(end),'AbsTol',tol);
            verifyEqual(testCase,t2(end),   tM(end),'AbsTol',tol);
            
            tol = 1e-3;
            % verify simulation results
            verifyLessThan(testCase,...
                    [norm(tM-t1),norm(tM-t2),...
                     norm(yM-y1),norm(yM-y2),...
                     norm(xM-x1),norm(xM-x2)],...
                     tol*ones(1,6))
        end
        
%         function testSISObench1(testCase)
%             load('building.mat');
%             sysSparse=sss(A,B,C);
%             sys=ss(full(A),full(B),full(C),zeros(1,1));
%             nOutputs = sysSparse.p;
%             
% %             isstable = isstable(sysSparse)
%             
%             Ts = 1e-4;
%             t = 0:Ts:10;
%             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
%             datau = iddata([],u',Ts);
%             
%             dataSparseRK4 = sim(sysSparse,datau,'RK4');
%             dataSparseFoEuler = sim(sysSparse,datau,'forwardEuler');
%             dataSparseBaEuler = sim(sysSparse,datau,'backwardEuler');
%             dataSparseDiscrete = sim(c2d(sysSparse,Ts),datau,'discrete');
%             
%             actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
%             expSim = lsim(sys,u,t);
%             
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
%             legend('RK4','ForwEuler','BackEuler','Discrete','lsim');
%             
%             for iCase=1:4
%                 verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
%             end
%         end
%         
%         function testSISObench2(testCase)
%             load('heat-cont.mat');
%             sysSparse=sss(A,B,C);
%             sys=ss(full(A),full(B),full(C),zeros(1,1));
%             nOutputs = sysSparse.p;
%             
% %             isstable = isstable(sysSparse)
%             
%             Ts = 1e-4;
%             t = 0:Ts:10;
% %             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
%             u = double(t>=Ts);
%             datau = iddata([],u',Ts);
%             
%             dataSparseRK4 = sim(sysSparse,datau,'RK4');
%             dataSparseFoEuler = sim(sysSparse,datau,'forwardEuler');
%             dataSparseBaEuler = sim(sysSparse,datau,'backwardEuler');
%             dataSparseDiscrete = sim(c2d(sysSparse,Ts),datau,'discrete');
%             
%             actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
%             expSim = lsim(sys,u,t);
%             
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
%             legend('RK4','ForwEuler','BackEuler','Discrete','lsim');
%             
%             for iCase=1:4
%                 verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
%             end
%         end
% 
%         function testDSSSISObench(testCase)
%             load('SpiralInductorPeec.mat');
%             sysSparse=sss(A,B,C,[],E);
%             sys=ss(sysSparse);
%             nOutputs = sysSparse.p;
%             
% %             isstable = isstable(sysSparse)
%             
%             Ts = 1e-10; % time constant of the system lies approx. by Tmax = 5e-8s
%             % Due to the Shannon-theorem: Ts = 1/2 * Tmax
%             t = 0:Ts:3e-7;
% %             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
%             u = double(t>=Ts);
%             datau = iddata([],u',Ts);
%             
%             dataSparseRK4 = sim(sysSparse,datau,'RK4'); %bad results; needs probably a smaller sample time Ts
%             dataSparseFoEuler = sim(sysSparse,datau,'forwardEuler'); %bad results; needs probably a smaller sample time Ts
%             dataSparseBaEuler = sim(sysSparse,datau,'backwardEuler'); %good results
%             dataSparseDiscrete = sim(c2d(sysSparse,Ts),datau,'discrete'); %bad results; needs probably a smaller sample time Ts
%             % Sample time of Ts = 1e-18 still yield bad results
%             
%             actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
%             expSim = lsim(sys,u,t);
%             
% %             figure; plot(t,actSim); hold on; plot(t,expSim); 
% %             xlabel('Time [s]'); ylabel('Amplitude');
% %             legend('RK4','ForwEuler','BackEuler','Discrete','lsim');
% 
%             figure; plot(t,dataSparseBaEuler.y); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
%             legend('BackEuler','lsim');
%             
%             for iCase=1:4
%                 verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
%             end
%         end
% 
%         function testMIMObench1(testCase)
%             load('iss.mat');
%             sysSparse=sss(A,B,C);
%             sys=ss(full(A),full(B),full(C),zeros(3,3));
%             nOutputs = sysSparse.p;
%             
% %             isstable = isstable(sysSparse);
% 
%             Ts = 1e-4;
%             t = 0:Ts:1;
% %             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
%             u = double(t>=Ts);
%                         
%             U = repmat(u,size(B,2),1);
%             dataU = iddata([],U',Ts);
%             
%             dataSparseRK4 = sim(sysSparse,dataU,'RK4'); %good results
%             dataSparseFoEuler = sim(sysSparse,dataU,'forwardEuler'); %good results
%             dataSparseBaEuler = sim(sysSparse,dataU,'backwardEuler'); %good results
%             dataSparseDiscrete = sim(c2d(sysSparse,Ts),dataU,'discrete'); %good results
%             
%             actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
%             expSim = lsim(sys,U,t);
%             
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
% 
%             for iCase=1:4
%                 verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
%             end
%         end
% 
%         function testMIMObench2(testCase)
%             load('CDplayer.mat');
%             sysSparse=sss(A,B,C);
%             sys=ss(full(A),full(B),full(C),zeros(2,2));
%             nOutputs = sysSparse.p;
%             
% %             isstable = isstable(sysSparse);
% 
%             Ts = 1e-7;
%             t = 0:Ts:1;
% %             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
%             u = double(t>=Ts);
%                         
%             U = repmat(u,size(B,2),1);
%             dataU = iddata([],U',Ts);
%             
%             dataSparseRK4 = sim(sysSparse,dataU,'RK4'); %good results
%             dataSparseFoEuler = sim(sysSparse,dataU,'forwardEuler'); %good results
%             dataSparseBaEuler = sim(sysSparse,dataU,'backwardEuler'); %good results
%             dataSparseDiscrete = sim(c2d(sysSparse,Ts),dataU,'discrete'); %good results
%             
%             actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
%             expSim = lsim(sys,U,t);
%             
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
% 
%             for iCase=1:4
%                 verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
%             end
%         end
% 
%         function testMIMObench3(testCase)
%             load('rail_1357.mat');
%             sysSparse=sss(A,B,C,[],E);
%             sys=dss(full(A),full(B),full(C),zeros(size(C,1),size(B,2)),full(E));
%             nOutputs = sysSparse.p;
%             
% %             isstable = isstable(sysSparse)
% 
%             Ts = 1e-4;
%             t = 0:Ts:2;
% %             u = idinput(length(t),'rgs',[0 0.5/(1/2/Ts)])';
%             u = double(t>=Ts);
%                         
%             U = repmat(u,size(B,2),1);
%             dataU = iddata([],U',Ts);
%             
%             dataSparseRK4 = sim(sysSparse,dataU,'RK4'); %good results
%             dataSparseFoEuler = sim(sysSparse,dataU,'forwardEuler'); %good results
%             dataSparseBaEuler = sim(sysSparse,dataU,'backwardEuler'); %good results
%             dataSparseDiscrete = sim(c2d(sysSparse,Ts),dataU,'discrete'); %good results
%             
%             actSim = [dataSparseRK4.y,dataSparseFoEuler.y,dataSparseBaEuler.y,dataSparseDiscrete.y];
%             expSim = lsim(sys,U,t);
%             
%             figure; plot(t,actSim); hold on; plot(t,expSim); 
%             xlabel('Time [s]'); ylabel('Amplitude');
% 
%             for iCase=1:4
%                 verification(testCase, actSim(:,(iCase*nOutputs-nOutputs+1):iCase*nOutputs), expSim);
%             end
%         end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution, expSolution, 'RelTol', 1e-1, 'AbsTol', 1e-2, ...
    'Difference between actual and expected exceeds relative tolerance');
end
