
%bcc multiple test settings: easy tests

exp_settings = settings.default_settings();

%IBCC
exp_settings.propKnown = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];
%exp_settings.propKnown = 0;
%exp_settings.propKnown = [0 0.1 0.4];
%exp_settings.propKnown = [0 0.1 0.25 0.4];

exp_settings.iPropKnown = 1;

exp_settings.lambdaSym = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
exp_settings.iLambdaSym = 1;

exp_settings.lambdaMag = [1 10 50 100];
exp_settings.iLambdaMag = 1;

%detailed runs use these settings
%propKnown = [0 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.2];
%lambdaSym = [0.5, 0.6, 0.7, 0.8, 0.9];
%lambdaMag = [0.1 0.5 1 5];

exp_settings.nDatasets = 10;
exp_settings.nRepeats = 3;

exp_settings.test_data_size = 100;

exp_settings.initCombinedPosteriors();

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
exp_settings.cluster_data = 'post';
exp_settings.weight_cluster_input = false;
exp_settings.top_save_dir = '/homes/49/edwin/matlab/combination/data/sandbox';
exp_settings.noOverwriteData = true;
%exp_settings.expLabel = 'dlrCombiners'; %experiment label
exp_settings.expLabel = 'comparison';

%other more boring things
exp_settings.save_imgs = true;
exp_settings.agent_type='dlr';

%data set sizes
exp_settings.training_set_size = 0;
exp_settings.length_history = 50;

%data generation
exp_settings.means = [-3 -2 -1; 3 2 1]; %[-3 -4 -5 -9 8; 3 4 5 9 -8];

%test with sensors all of same quality
exp_settings.deviations = [1 2 1 ; 1 1 2];

exp_settings.p_c1 = 0.5;

exp_settings.noninf_mean = 0;
exp_settings.noninf_dev = 5;

%sensors

%settings used: 
% inf =5 noninf = 11; max =10 min = 3
% inf = 3 noninf = 8
% inf = 2, noninf = 2; max = 3, min = 1
% inf = 1 noninf = 1

exp_settings.num_inf_sensors = 3;
exp_settings.num_noninf_sensors = 7;

exp_settings.p_switch = 0.0; % probability that a sensor will stop being useful
exp_settings.p_missing = 0.00;
exp_settings.p_corrupt = 0.0;
exp_settings.p_flip_label = 0.0;
exp_settings.p_flip_sensor = 0.0;

%agents and sensor allocation
exp_settings.num_agents = 10;

exp_settings.max_sensors = 6;
exp_settings.min_sensors = 2;

%number of clusters and clustering window length
exp_settings.sequence_length = 20;
exp_settings.K = 2;