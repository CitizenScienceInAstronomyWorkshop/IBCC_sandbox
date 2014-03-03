% 
% aucs_greedy = cell(1,10);
% lam_greedy = cell(1,10);
% 
% aucs_clus = cell(1,10);
% lam_clus = cell(1,10);
% 
% aucs_rand = cell(1,10);
% lam_rand = cell(1,10);
% 
% aucs_mostUncert = cell(1,10);
% lam_mostUncert = cell(1,10);

for repeat=1:1
    display(['repeats: ' num2str(repeat)]);
%     clear initLabelIds
%         
% %     selectMethod = 'uncertGreedy';    
% %     
% %     display('Two-stage, greedy')
% %     
% %     classifierMethod = 'two-stage';
% % 
% %     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
% %     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% % 
% %     runSelP_selMethods
% %     
% %     aucs_greedy{repeat} = aucs;
% %     lam_greedy{repeat} = lam_uncert;    
    
    selectMethod = 'uncert';    
    
    display('Two-stage, uncert')
    
    classifierMethod = 'two-stage';

    featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
    fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';

    runSelP_selMethods
    
    aucs_mostUncert{repeat} = aucs;
    lam_mostUncert{repeat} = lam_uncert;    

    selectMethod = 'random';    
    
    display('Two-stage, random')
    
    classifierMethod = 'two-stage';

    featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
    fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';

    runSelP_selMethods
    
    aucs_rand{repeat} = aucs;
    lam_rand{repeat} = lam_uncert;     
        
    selectMethod = 'uncertClust';    
    
    display('Two-stage, clus')
    
    classifierMethod = 'two-stage';

    featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
    fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';

    runSelP_selMethods
    
    aucs_clus{repeat} = aucs;
    lam_clus{repeat} = lam_uncert;    
%     
%     display('Two-stage')
%     
%     classifierMethod = 'two-stage';
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP_selMethods

%     display('VB: uncertain worker labels')
%     
%     classifierMethod = 'VB-workerUncert';
%     crowdTrust = [20 80];
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP_selMethods
%     
%     display('VB: uncertain worker labels')
%     
%     classifierMethod = 'VB-workerUncert';
%     crowdTrust = [10 90];
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP_selMethods
% %     
%     display('VB: uncertain worker labels')
%     
%     classifierMethod = 'VB-workerUncert';
%     crowdTrust = [4 16];
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP_selMethods
% 
%     display('VB: uncertain labels, crowd as unit')
%     
%     classifierMethod = 'VB-uncertLabels';
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP_selMethods
%     
%     display('VB: assume correct labels')
%     
% 
%     classifierMethod = 'VB';
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP_selMethods
%     
%     display('1500 features');
% 
%     featureFile = '/homes/49/edwin/data/trec/trec8_sample/matrices_1500_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/trec8_sample/filemap.csv';
% 
%     runSelP_selMethods
end

fclose('all');