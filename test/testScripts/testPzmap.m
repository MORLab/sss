classdef testPzmap < sssTest
    % testPzmap - testing of pzmap.m
    %
    % Description:
    %   The function pzmap.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs:
    %    pzmap(sys),[p,z]=pzmap(sys)
    %
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
    %               Jorge Luiz Moreira Silva
    % Last Change:  05 Nov 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
   
    methods(Test)
        function test1(testCase)
            for i=1:length(testCase.sysCell)
                sysSparse=testCase.sysCell{i};
                sys=ss(sysSparse);
                [actP,actZ]=pzmap(sysSparse);
                actPZ={sort(actP),sort(actZ)};
                pzmap(sysSparse);
                [expP,expZ]=pzmap(sys);
                expPZ={sort(expP),sort(expZ)};
                verification(testCase, actPZ, expPZ);
                close all
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',0.1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end