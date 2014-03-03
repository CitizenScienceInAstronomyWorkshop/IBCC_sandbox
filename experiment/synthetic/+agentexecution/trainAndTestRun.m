function [agents basePost baseLogodds agentConfidence] = trainAndTestRun(...
    synthSet, labelledTrainingData, labelledTestData, reuseSubs)
%generate synthetic data, create logistic regression agents with subsets of
%data inputs, train and classify new data points. Cluster the agents based
%on similarity of results.

%experimental_settings
N = size(labelledTestData, 1);%synthSet.nSamples;
A = synthSet.nAgents;
S = synthSet.nSensors();

%allocate sensors to agents

subscriptions = synthSet.subs;

if isempty(subscriptions) || reuseSubs == false
    subscriptions = agentexecution.subscribeSensors(A, S, synthSet);
end

%train agents on the dataset
agents = {A};

Ntr = synthSet.nTrainingSamples;
for i=1:A
    if isempty(labelledTrainingData)
        D = [];
        T = [];
    else
        trainingBatch = i;
        if ~isempty(synthSet.trainingBatch)
            trainingBatch = synthSet.trainingBatch(i);
        end
        D = labelledTrainingData( ((trainingBatch-1)*Ntr+1):(trainingBatch*Ntr), 1:S);
        T = labelledTrainingData( ((trainingBatch-1)*Ntr+1):(trainingBatch*Ntr), S+1);
    end
    
    s = subscriptions{i};
    
    if strcmp(synthSet.agentType, 'logistic')
        import baselearners.logisticRegressionAgent;
        agent = logisticRegressionAgent(s, i, synthSet.lengthHistory, synthSet.cal(i), synthSet.bias(i));
    elseif strcmp(synthSet.agentType, 'dlr')  
        import baselearners.dlrAgent;        
        agent = dlrAgent(s, i, false);
    end
    
    agent.train(D, T);
    agents{i} = agent;
end

basePost = zeros(A, N);
baseLogodds = zeros(A, N);

% classifications = zeros(N, 1);
test_data = zeros(N, S);

% sensor_credit = zeros(S,N);
% sensor_weight = zeros(S,N);

%load a matrix for storing the test data and its labels and changepoint
%labels
%labelledTestData = dlmread('test_data.dat');


test_data(1:N, :) = labelledTestData(1:N, 1:S);

credit_duration = 10.0;
credit_range = 1:credit_duration;
historic_weights = credit_range ./ (credit_duration * ones(1, credit_duration));
historic_weights = ones(S,1) * historic_weights;

for i=1:A
    agent = agents{i};
    if exist('labels_for_agents', 'var')
        agent_labels = labels_for_agents(:, i);
    else
        agent_labels = labelledTestData(:, S+1);
    end
    agent.init_test(test_data, agent_labels);
end
display(['testing ' num2str(N) ' samples']);
for n=1:N
    
    %display(sprintf('testing:  %i', n));
    
    point = labelledTestData(n, 1:S);
    classification = labelledTestData(n, S+1);
                
    sum_agent_sensors = ones(S, 1);
        
    for i=1:A    
        agent = agents{i};
                
        [post, a] = agent.classify(n);
                     
        basePost(i, n) = post;
        baseLogodds(i, n) = a;
        
%         %train the agent on this data point as well
%         agent.retrain(n);
%         
%         if classification ~= -1
%             agent_credit = (0.5 - abs(classification-post)) ./ sum(agent.s);
%         
%             sum_agent_sensors = sum_agent_sensors + ...
%                 ((classification-post).^2) .* transpose(agent.s);
%         
%             sensor_credit(:, n) = sensor_credit(:, n) + ...
%                 (agent_credit .* transpose(agent.s) );
%                 
%             weights = agent.W(2:agent.num_active_sensors+1);
%             sensor_weight(agent.s==1, n) = sensor_weight(agent.s==1, n) + ...
%                 (agent_credit .* weights ./ sum(abs(weights)) );
%         end
        
    end
    
%     sensor_credit(:, n) = sensor_credit(:, n) ./ sum_agent_sensors;
%     sensor_weight(:, n) = sensor_weight(:, n) ./ sum_agent_sensors;
%     
%     c = n-credit_duration+1;
%     this_duration = credit_duration;
%     if c <= 0 
%         c = 1;
%         this_duration = n;
%     end
%     
%     credit_range = credit_duration-this_duration+1:credit_duration;
%     
%     sensor_credit(:, n) = sensor_credit(:, n) + sum(historic_weights(:, credit_range) .* sensor_credit(:, c:n), 2);
%     sensor_weight(:, n) = sensor_weight(:, n) + sum(historic_weights(:, credit_range) .* sensor_weight(:, c:n), 2);
%     
end

% sensor_acc_credit = scoring.calculate_sensor_credit(agents, A, S, N, classifications, basePost);

display('Trained and tested agents.');
% 
% if strcmp(synthSet.clusterData, 'logodds')
%     data = baseLogodds;
% elseif strcmp(synthSet.clusterData, 'post')
%     data = basePost;
% else
%     display('Defaulting to agent confidence based on logodds: unrecognised setting');
%     data = baseLogodds;
% end

% labels = labelledTestData(1:N, S+1);
% test_data = labelledTestData(1:N, 1:S);
% analyser = scoring.classification_analyser(data, labels, test_data);   
agentConfidence = zeros(A, N);
% for n=1:N
%     agentConfidence(1:A, n) = analyser.agentConfidence(1:A, n);
% end
% 
% display('Done. Calculated agent confidences');
end
