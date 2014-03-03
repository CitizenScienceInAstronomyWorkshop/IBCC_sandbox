for repeats=1:1
    display(['repeats: ' num2str(repeats)]);
    clear initLabelIds
    selectMethod = 'uncertClust';    
%     
    
    display('Two-stage-G')
    
    classifierMethod = 'two-stage-G';

    featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
    fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';

    runSelP
    title('Two-stage-G Classifier TREC 7');
    
    display('Two-stage')
    
    classifierMethod = 'two-stage';

    featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
    fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';

    runSelP
    title('Two-stage Classifier TREC 7');

%     display('VB: uncertain worker labels')
%     
%     classifierMethod = 'VB-workerUncert';
%     crowdTrust = [20 80];
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP
%     
%     display('VB: uncertain worker labels')
%     
%     classifierMethod = 'VB-workerUncert';
%     crowdTrust = [10 90];
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP
%     
    display('VB: uncertain worker labels')
    
    classifierMethod = 'VB-workerUncert';
    crowdTrust = [4 16];

    featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
    fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';

    runSelP
    title('VB Classifier TREC 7');
   
% 
%     display('VB: uncertain labels, crowd as unit')
%     
%     classifierMethod = 'VB-uncertLabels';
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP
%     
%     display('VB: assume correct labels')
%     
% 
%     classifierMethod = 'VB';
% 
%     featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% 
%     runSelP
%     
%     display('1500 features');
% 
%     featureFile = '/homes/49/edwin/data/trec/trec8_sample/matrices_1500_topic.txt';
%     fileMapFile = '/homes/49/edwin/data/trec/trec8_sample/filemap.csv';
% 
%     runSelP
end