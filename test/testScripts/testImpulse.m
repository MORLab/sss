classdef testImpulse < sssTest
    % testImpulse - testing of impulse.m   
    %
    % ------------------------------------------------------------------
    %   This file is part of sssMOR, a Sparse State Space, Model Order
    %   Reduction and System Analysis Toolbox developed at the Institute
    %   of Automatic Control, Technische Universitaet Muenchen.
    %   For updates and further information please visit www.rt.mw.tum.de
    %   For any suggestions, submission and/or bug reports, mail us at
    %                     -> sssMOR@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto, Maria Cruz Varona
    %               Jorge Luiz Moreira Silva, Stefan Jaensch
    % Last Change:  18 Apr 2016
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function testImpulse1(testCase)
            filePath = fileparts(mfilename('fullpath'));
            toolPath = fileparts(filePath);
            load([toolPath '/../benchmarks/iss.mat'])
            
            sys = ss(full(A),full(B),full(C),[]);
            [sys,g] = balreal(sys);
            sys = modred(ss(sys),g<1e-3)+1;
            sysSss = sss(sys);
            sysSssE = sss(sys.A*2,sys.B*2,sys.C,sysSss.D,eye(size(sys.A))*2);                       
           
            tFinal = 15;
            [actH,actT]=impulse(ss(sysSss),tFinal);
            expH = {}; expT = {};
            ODEset = odeset;
            ODEset.AbsTol = 1e-7;
            ODEset.RelTol = 1e-7;
            [expH{end+1},expT{end+1}]=impulse(sysSss,actT);
            [expH{end+1},expT{end+1}]=impulse(sysSss,actT,struct('nMin',0,'ode','ode113','odeset',ODEset));
            [expH{end+1},expT{end+1}]=impulse(sysSss,actT,struct('nMin',0,'ode','ode15s','odeset',ODEset));
            [expH{end+1},expT{end+1}]=impulse(sysSss,actT,struct('nMin',0,'ode','ode23','odeset',ODEset));
            [expH{end+1},expT{end+1}]=impulse(sysSss,actT,struct('nMin',0,'ode','ode45','odeset',ODEset));
            [expH{end+1},expT{end+1}]=impulse(sysSssE,actT,struct('nMin',0));            
            
            actSolution={actH+1,actT};
            for i = 1:length(expH)
                expSolution={expH{i}+1,expT{i}};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actH, 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actT , 'double', 'Instances not matching');
                verifySize(testCase, actH, size(expH{i}), 'Size not matching');
                verifySize(testCase, actT, size(expT{i}), 'Size not matching');
            end
        end
        function testImpulsePlot(testCase)
            filePath = fileparts(mfilename('fullpath'));
            toolPath = fileparts(filePath);
            load([toolPath '/../benchmarks/iss.mat'])
            
            sys = ss(full(A),full(B),full(C),[]);
            [sys,g] = balreal(sys);
            sys = modred(ss(sys),g<1e-3)+1;
            sys.InputName = {'u1';'u2';'u3';};
            sys.OutputName = {'y1';'y2';'y3';};
            sys = sss(sys);
            sysSss = sys;
            sysSssE = sss(sys.A*2,sys.B*2,sys.C,sysSss.D,eye(size(sys.A))*2);                                   
            sys1_1 = sys(1,1);
            sys12_3 = sys(1:2,3);
            sysSS = ss(sys);
            tFinal = 15;
            figure            
            impulse(sysSss,sysSssE,sys1_1,sys12_3,sysSS,tFinal)
            title('testImpulsePlot');
        end
        function testImpulseBasic(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    ODEset = odeset;
                    ODEset.AbsTol = 1e-12;
                    ODEset.RelTol = 1e-12;
                    [actSolution,t]=impulse(sys,struct('nMin',0,'odeset',ODEset));

                    expSolution=impulse(ss(sys),t);
                    
                    % impulse response oscillates strongly, so a few
                    % failing indices are expected (RelTol 0.05, AbsTol
                    % 1e-4)
                    failIndices=zeros(size(actSolution,1),1);
                    for j=1:size(actSolution,1)
                        for k=1:size(actSolution,2)
                            for l=1:size(actSolution,3)
                                if abs(actSolution(j,k,l)-expSolution(j,k,l))/abs(expSolution(j,k,l)) > 0.05 && abs(actSolution(j,k,l)-expSolution(j,k,l))>1e-4
                                    failIndices(j)=1;
                                end
                            end
                        end
                    end
                    verifyLessThanOrEqual(testCase,nnz(failIndices)/size(actSolution,1),0.16);
                end
            end
        end
        function testImpulseTime(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    
                    ODEset = odeset;
                    ODEset.AbsTol = 1e-10;
                    ODEset.RelTol = 1e-10;
                    
                    % time vector
                    t=0.1:0.1:0.5;
                    [actSolution]=impulse(sys,t,struct('nMin',0,'odeset',ODEset));
                    expSolution=impulse(ss(sys),t);
                    verification(testCase,actSolution(2:end-1,:),expSolution(2:end-1,:));

                    % final time
                    Tfinal=1;
                    [~,t]=impulse(sys,Tfinal,struct('nMin',0,'odeset',ODEset));
                    verifyEqual(testCase,t(end),Tfinal,'AbsTol',0.06);
                end
            end
        end
        function testImpulseMultiSys(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    sys2=loadSss('building');
                    sys3=loadSss('CDplayer');
                    Tfinal=0.5;
                    
                    % test call
                    impulse(sys,'b-',ss(sys2),'r--',sys3,'g:',Tfinal,struct('nMin',0));
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',0.16,'AbsTol',5e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end