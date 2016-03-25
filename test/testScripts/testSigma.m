classdef testSigma < sssTest
    % testSigma - testing of signa.m
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
    %                     -> sssMOR@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto
    %               Jorge Luiz Moreira Silva
    % Last Change:  26 Out 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function testBench(testCase)
            for i=1:length(testCase.sysCell)
                if testCase.sysCell{i}.isSiso
                    sysSparse=testCase.sysCell{i};
                    sys=ss(sysSparse);
                    sigma(sysSparse);
                    [expMag,omega]=sigma(sys);
                    mag=sigma(sysSparse, omega');
                    actMag=zeros(1,size(mag,3));
                    actMag(:)=mag(1,1,:);
                    verification(testCase, sort(actMag), sort(expMag(1,:)));
                    close all;
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-3,'AbsTol',1e-4,...
    'Difference between actual and expected exceeds relative tolerance');
end