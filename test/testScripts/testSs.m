classdef testSs < sssTest
    % testSs - testing of ss.m
    %
    % Description:
    %   The function norm.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs: sys=ss(sysSparse)
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
                sys=ss(sysSparse);
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution(1:4),  expSolution(1:4),'RelTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end