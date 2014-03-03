classdef BccSettings < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        debug = false;
        
        %INFERENCE/MODELLING ----------------------------------------------
        
        useLikelihood = false; %defer to likelihood only (ignore class proportions)
        priorAdjLevel = 100;
        
        %SCORE MAPPING ----------------------------------------------------
        zeroScoreMap = 0;
        IbccMapFunction = [];%@mapScoresDoNothing; 
        %IbccMapFunction = @mapGZSNScoresToVotes; %for making responses into binary votes        
                
        scoreMap = [1 2 0; 0 1 -1];
        minScore = 0;
        maxScore = 1;        
        %scoreMap = [1 2 3 0; -1 3 1 0];
        
        %HYPERPARAMETERS --------------------------------------------------
        nu = {[3 3]};
        iNu = 1;
        
        Alpha = [0.5 0.3 0.05; 0.18 0.36 0.41]; %*used in paper*
        % Alpha = [2.5 0.5 1.5; 1.5 3.5 3];            
        % Alpha = [1 1 1; 1 1 1] ./ 10;
        % Alpha = [56 28 7; 42 77 71] ./ 100; %based on true proportions for all GZSN data.
        % Alpha = [10 10 2; 10 10 10];
        % Alpha = [0.7 0.26 0.04; 0.1 0.35 0.55];
        % Alpha = [0.5 0.3 0.05; 0.18 0.36 0.41] .* 20;
        % Alpha = [0.7 0.3; 0.2 0.8] .* 10;
        %     Alpha = [0.2 0.2 0.2; 0.2 0.2 0.2];
        % Alpha = [0.2 0.2 1; 0.5 0.7 0.1];
                
        AlphaDiff = [0.5 0.3 0.05; 0.3 0.3 0.4; 0.18 0.36 0.41; 0.4 0.3 0.3];

        %EXPERT LABEL HANDLING --------------------------------------------
        
        targetsAsAgents = [];
        
        %Treat the ptfType labels as data from a reliable agent. Otherwise they
        %will be ignored. For test data points, these values will be ignored, so
        %this only makes a difference to the results if includeUnlabelled is true.         
        screenLabelsAsAgent = true;
        trustedAlpha = [];
        trustedAgentProb = 0.99; %make this a column vector if we have more trusted agents
        
        %ITERATIONS -------------------------------------------------------
        convThreshold = 10^-2;
        maxIt = 1000;
        fixedNIt = 0; %250; %use a fixed number of iterations
        convIt = 2;
        
        %DYNAMIC IBCC -----------------------------------------------------
        changeRateMod = 1;
        
        %LEGACY -----------------------------------------------------------
        %for IBCC models using a distribution over the Alphas
        propKnown = [0 0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.3 0.4];
        iPropKnown = 1;
        
        lambdaSym = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
        iLambdaSym = 1;
        
        lambdaMag = [1 10 50 100];
        iLambdaMag = 1;
        
        trustedLambdaSym = [];
        trustedLambdaMag = [];
    
        %detailed runs use these settings
        %propKnown = [0 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.2];
        %lambdaSym = [0.5, 0.6, 0.7, 0.8, 0.9];
        %lambdaMag = [0.1 0.5 1 5];
    end
    
    methods
        function obj = BccSettings()
            %each cell contains posteriors: each set of results is a row,
            %each row is a repeat with same settings.
            obj.targetsAsAgents = [];
        end   
    end
    
end

