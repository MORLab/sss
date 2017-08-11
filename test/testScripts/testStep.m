classdef testStep < sssTest
    % testStep - testing of step.m   
    %
    % ------------------------------------------------------------------
    %   This file is part of sssMOR, a Sparse State Space, Model Order
    %   Reduction and System Analysis Toolbox developed at the Institute
    %   of Automatic Control, Technische Universitaet Muenchen.
    %   For updates and further information please visit www.rt.mw.tum.de
    %   For any suggestions, submission and/or bug reports, mail us at
    %                     -> morlab@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto, Maria Cruz Varona
    %               Jorge Luiz Moreira Silva, Stefan Jaensch
    % Last Change:  18 Apr 2016
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function testStep1(testCase)
            filePath = fileparts(mfilename('fullpath'));
            toolPath = fileparts(filePath);
            load([toolPath '/../benchmarks/iss.mat'])
            
            sys = ss(full(A),full(B),full(C),[]);
            [sys,g] = balreal(sys);
            sys = modred(ss(sys),g<1e-3)+1;
            sysSss = sss(sys);
            sysSssE = sss(sys.A*2,sys.B*2,sys.C,sysSss.D,eye(size(sys.A))*2);                                   
            
            tFinal = 15;
            [actH,actT]=step(ss(sysSss),tFinal);
            expH = {}; expT = {};
            ODEset = odeset;
            ODEset.AbsTol = 1e-7;
            ODEset.RelTol = 1e-7;
            [expH{end+1},expT{end+1}]=step(sysSss,actT);
            [expH{end+1},expT{end+1}]=step(sysSss,actT,struct('nMin',0,'ode','ode113','odeset',ODEset));
            [expH{end+1},expT{end+1}]=step(sysSss,actT,struct('nMin',0,'ode','ode15s','odeset',ODEset));
            [expH{end+1},expT{end+1}]=step(sysSss,actT,struct('nMin',0,'ode','ode23','odeset',ODEset));
            [expH{end+1},expT{end+1}]=step(sysSss,actT,struct('nMin',0,'ode','ode45','odeset',ODEset));
            [expH{end+1},expT{end+1}]=step(sysSssE,actT,struct('nMin',0));            
            
            actSolution={actH,actT};
            for i = 1:length(expH)
                expSolution={expH{i},expT{i}};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actH, 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actT , 'double', 'Instances not matching');
                verifySize(testCase, actH, size(expH{i}), 'Size not matching');
                verifySize(testCase, actT, size(expT{i}), 'Size not matching');
            end
        end
        function testStepPlot(testCase)
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
            step(sysSss,sysSssE,sys1_1,sys12_3,sysSS,tFinal)
            title('testStepPlot');
        end
        function testStepBasic(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && ~strcmp(sys.Name,'CDplayer')
                    ODEset = odeset;
                    ODEset.AbsTol = 1e-10;
                    ODEset.RelTol = 1e-10;
                    [expSolution,t]=step(ss(sys));
                    actSolution=step(sys,t,struct('nMin',0,'odeset',ODEset));
                    verification(testCase,actSolution,expSolution);
                end
            end
        end
        function testStepTime(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    ODEset = odeset;
                    ODEset.AbsTol = 1e-8;
                    ODEset.RelTol = 1e-8;
                    
                    % time vector
                    t=0.1:0.1:0.3;
                    [actSolution]=step(sys,t,struct('nMin',0,'odeset',ODEset));
                    expSolution=step(ss(sys),t);
                    verification(testCase,actSolution,expSolution);

                    % final time
                    Tfinal=0.3;
                    [~,t]=step(sys,Tfinal,struct('nMin',0,'odeset',ODEset));
                    verifyEqual(testCase,t(end),Tfinal,'AbsTol',0.05);
                end
            end
        end
        function testStepMultiSys(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    close all
                    sys2=sss('building');
                    t=0.01:0.01:0.03;
                    Tfinal=0.01;
                    
                    % test different calls
                    step(sys,sys2,t,struct('nMin',0,'tsMin',1e-3));
                    step(sys,ss(sys2),Tfinal,struct('nMin',0,'tsMin',1e-3));
                    step(sys,'b-',sys2,'r--',Tfinal,struct('nMin',0,'tsMin',1e-3));
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1.1e-3,'AbsTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end