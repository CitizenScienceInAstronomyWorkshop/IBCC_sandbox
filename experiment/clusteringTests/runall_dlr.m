%Run this code from this directory
clear
import settings.*

%load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');

experimental_settings_x5_20_dlr

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
expSettings.clusterData = 'post';
expSettings.weightClusterInput = false;
expSettings.topSaveDir = '/homes/49/edwin/results/dlr_clustering_mean_mag_unweighted';
%reuseSubs = true;

runner = ExpRunner(expSettings);
runner.runWithGraphs(true, true, 52088, 'On');

close all

%WEIGHTED

%clear

% addpath('dlr');
% 
% %load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');
% %experimental_settings_x5_20_dlr
% 
% 
% graph_folder = 52011;
% 
% %cluster agents over a data sequence of 20. use posteriors as cluster input
% %data.
% topSaveDir = ('/homes/49/edwin/results/dlr_clustering_mean_mag_unweighted');
% 
% clusterData = 'post';
% weightClusterInput = true;
% 
% reuseSubs = true;
% 
% %generateData
% %trainAndTestRun
% cluster_agents_kmeans
% 
% graph_agent_performance
% graph_clusters
% 
% save('/homes/49/edwin/results/dlr_clustering_mean_mag_unweighted/data/52011.mat');
% 
% close all

% addpath('dlr');
% 
% %load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');
% %experimental_settings_x5_20_dlr
% 
% 
expSettings.graph_folder = 53011;

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
expSettings.topSaveDir = ('/homes/49/edwin/results/dlr_clustering_mean_mag_unweighted');

expSettings.clusterData = 'post';
expSettings.weightClusterInput = true;

expSettings.reuseSubs = true;
expSettings.seqLength = 30;

runner = ExpRunner(expSettings);
runner.runWithGraphs(false, false, expSettings.graph_folder, 'On');

close all