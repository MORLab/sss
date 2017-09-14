classdef testSigma < sssTest
    % testSigma - testing of sigma.m
    %
    % Description:
    %   The function sigma.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case just if it runs for the syntax
    %    [mag,omega]=sigma(sys) and sigma(sys)
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
        function testBench(testCase)
            for i=1:length(testCase.sysCell)
                sysSparse=testCase.sysCell{i};
                if ~sysSparse.isDae
                    [expMag,omega]=sigma(ss(sysSparse));
                    actMag=sigma(sysSparse, omega');
                    verification(testCase, actMag, expMag);
                end
            end
        end
        function inputFunctionality(testCase)
            % input of frequnecy range
            
            for i=1:length(testCase.sysCell)
                sys  = testCase.sysCell{i};
                if ~sys.isDae
                    w = {10,100};

                    [~,omega] = sigma(sys,w);
                    verifyEqual(testCase,omega(1),w{1},  'Wrong frequency returned');
                    verifyEqual(testCase,omega(end),w{2},'Wrong frequency returned');
                end
            end
        end
        
        function plotFunctionalitySISO(testCase)
            % verify the correct plot compared to built-in when omega is
            % not passed
            sys1 = sss('building');
            sys2 = sss('beam');
            sys3 = sss('eady');
            sys4 = sss('fom');
            sys5 = sss('iss'); sys5 = sys5(1,1);
            
            figure; sigma(sys1,sys2,'r-',sys3,'k--',sys4,sys5);
            figure; sigma(ss(sys1),ss(sys2),'r-',ss(sys3),'k--',ss(sys4),ss(sys5));
            figure; sigma(sys1,sys2,'r-',sys3,'k--',sys4,sys5,{1e-3,1e6}); % frequency range is passed
        end
        
        function plotFunctionalityMIMO(testCase)
            % verify the correct plot compared to built-in when omega is
            % not passed
            sys1 = sss('CDPlayer');
            sys2 = sss('iss'); sys2 = sys2(1:2,1:2);
            
            figure; sigma(sys1,sys2,'r-');
            figure; sigma(ss(sys1),ss(sys2),'r-');
            figure; sigma(sys1,sys2,'r-',{1e-3,1e6}); % frequency range is passed
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-2,'AbsTol',0.005,...
    'Difference between actual and expected exceeds relative tolerance');
end