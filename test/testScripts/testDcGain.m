classdef testDcGain < sssTest
    % testDcGain - testing of dcgain.m   
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
    %               Lisa Jeschek
    % Last Change:  20 May 2016
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function testDcGain1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    actSolution=dcgain(sys);
                    expSolution=dcgain(ss(sys));
                    verification (testCase, actSolution, expSolution);
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1.1e-3,'AbsTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end