classdef testNorm < sssTest
    % testNorm - testing of norm.m
    %
    % Description:
    %   The function norm.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs (norm(sys),
    %       norm(sys,inf), [n,fpeak]=norm(sys,inf), norm(sys,2)
    %
    % ------------------------------------------------------------------
    %   This file is part of sssMOR, a Sparse State Space, Model Order
    %   Reduction and System Analysis Toolbox developed at the Institute
    %   of Automatic Control, Technische Universitaet Muenchen.
    %   For updates and further information please visit www.rt.mw.tum.de
    %   For any suggestions, submission and/or bug reports, mail us at
    %                     -> morlab@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto
    %               Jorge Luiz Moreira Silva
    % Last Change:  26 Out 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function testNorm1(testCase)
            warning('off','sss:isstable:EigsFailed');
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                if ~sys_sss.isDae
                    sys_ss=ss(sys_sss);
                    
                    actNorm1=norm(sys_sss);
                    [actNorm2,actFreq]=norm(sys_sss,inf);
                    actNorm = [actNorm1,actNorm2,actFreq];
                    
                    expNorm1=norm(sys_ss);
                    [expNorm2,expFreq]=norm(sys_ss,inf);
                    expNorm = [expNorm1,expNorm2,expFreq];
                    verification(testCase, actNorm, expNorm);
                end
            end
        end
        function testLyapchol(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.n>300
                        % options for mess
                        % eqn struct: system data
                        eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'type','N','haveE',sys.isDescriptor);

                        % opts struct: mess options
                        messOpts.adi=struct('shifts',struct('l0',20,'kp',50,'km',25,'b0',ones(sys.n,1),...
                            'info',0),'maxiter',300,'restol',0,'rctol',1e-12,...
                            'info',0,'norm','fro');
                        oper = operatormanager('default');
                        
                        % get adi shifts
                        [messOpts.adi.shifts.p,~,~,~,~,~,~,eqn]=mess_para(eqn,messOpts,oper);

                        % low rank adi
                        [R,~,eqn]=mess_lradi(eqn,messOpts,oper);

                        actNrmAdi1=norm(R'*eqn.C','fro');
                        Opts.lyapchol='adi';
                        actNrmAdi2=norm(sys,Opts);
                        Opts.lyapchol='builtIn';
                        actNrmBuiltIn=norm(sys,Opts);
                        expNrm=norm(ss(sys));
                        
                        actSolution={actNrmAdi1, actNrmAdi2, actNrmBuiltIn};
                        expSolution={expNrm, expNrm, expNrm};
                        verifyEqual(testCase, actSolution, expSolution,'RelTol',1e-3,...
                        'Difference between actual and expected solution.');
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end