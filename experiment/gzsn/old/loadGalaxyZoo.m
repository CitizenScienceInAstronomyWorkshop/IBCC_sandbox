fid = fopen('/homes/49/edwin/matlab/combination/data/galaxyZoo3/resultset2.csv');
snData = textscan(fid, '%d%d%d%d%d%d%s%s%d%d%d%d%d%s', 'Delimiter',',');
fclose(fid);

%importdata('/homes/49/edwin/matlab/combination/data/galaxyZoo3/resultset_sample.csv', ',');
%importdata('/homes/49/edwin/matlab/combination/data/galaxyZoo2/gz.dat', ' ');
%importdata('/homes/49/edwin/matlab/combination/data/galaxyZoo/resultset_deduped.csv', ',');

%data cleansing - hopefully not necessary now
if exist('cleanData', 'var') && cleanData

    %set values that were not scored for answers 2, 4, and 9 to -1 as they should be (data cleaning),
    %leave empty entries as 0.
    snData.data(snData{4}==2, 3) = -1;
    snData.data(snData{4}==4, 3) = -1;
    snData.data(snData{4}==9, 3) = -1;

    display(sprintf('No. duplicate -1 scores: %d', sum(snData.data(:, 3)< -1)));
    snData.data(snData{3}< -1, 3) = -1;
    snData.data( isnan(snData{3}), 3) = 0;

    %make entries with null zooniverse user id to a new user id, i.e. one
    %larger than the maximum existing one.
    snData.data((snData{2} == 0), 2) = max(snData.data(:, 2)) + 1;


    %make the data into integers
    %snData.data(:, 2) = int32(snData.data(:, 2));

    %remove invalid task/answer annotations (invalid as you can't repeat the
    %same tree path; we are assuming that if the classification id is the same,
    %then the duplicate annotations represent data that should have been
    %deleted when the user pressed back. If the classification ID is different,
    %the user has classified the asset more than once, which is valid.)
    for dupeTask=1:max(snData{12})
        dupeTaskIdxs = find(snData{12}==dupeTask);
        firstMatrix = sparse(ones(length(dupeTaskIdxs),1), double(snData{2}(dupeTaskIdxs)), double(snData{12}(dupeTaskIdxs)));

        invClassIds = find(firstMatrix>1);
        invAnnIds = [];
        for c=invClassIds
            annIds = snData{1}(snData{2}==c & snData{12}==dupeTask);
            if isempty(annIds)
                continue;
            end
            validAnnId = annIds(end);
            %just set to zero and ignore for now if another task was completed that
            %should supersede the first
            snData{13}(snData{2}==c & snData{1}<validAnnId & snData{12}>=dupeTask) = 0; 
            snData{9}(snData{2}==c) = sum(snData{13}(snData{2}==c));
        end
    end

end

save('/home/edwin/work_overflow/galaxyZoo3/snData.mat', 'snData');

% snMatrix = sparse(double(snData{3}), double(snData{4}), double(snData{9}));
% %snMatrix = sparse(snData.data(:, 2), snData.data(:, 1), snData.data(:, 3));
%  
% if exist('cleanData', 'var') && cleanData
%     %Some entries appear to have followed an invalid tree path when answering
%     %questions. Let's just hope it's only those that ended up with invalid
%     %scores - set them to sensible scores.
%     snMatrix(snMatrix == -2) = -1;
%     %these appear to be missing the data about how they reached the final task.
%     %The answer to this should result in a score of 3 as 1 more would be
%     %awarded by the previous question.
%     snMatrix(snMatrix == 2) = 3;
% 
%     %remove rows and columns of all zeros
%     snMatrix = snMatrix(:, sum(abs(snMatrix), 1)~=0);
%     snMatrix = snMatrix(sum(abs(snMatrix), 2)~=0, :);
% end
% 
% % snMatrix(snMatrix < -1) = -1;

%use reloadGalaxyZoo instead of the stuff below

% import settings.*
% vbIbccPaper.exp2seq;
% 
% %pick the combination methods to use.
% combMethods = {...
%     combiners.MeanDecision.shortLabel,...
%     combiners.SimpleMajorityVoting.shortLabel, ...
%     combiners.WeightedMajority.subShortLabel, ...
%     combiners.WeightedSum.shortLabel,...
%     combiners.IbccVb.shortLabel,...    
%     %combiners.Ibcc.shortLabel, ...
%     };
% 
% snBaseOutputs = cell(3,1);
% snBaseOutputs(1) = snData(3);
% snBaseOutputs(2) = snData(4);
% snBaseOutputs(3) = snData(9);
% 
% runner = gzRunner(expSettings, snBaseOutputs, []);
% runner.runCombineAll(false, 'dontTouch', true, combMethods, true); %new_data, runBaseAgents, drawGraphs, combMethods)
