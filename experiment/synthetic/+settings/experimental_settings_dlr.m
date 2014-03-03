exp_settings.save_imgs = true;
%top_save_dir = ('/homes/49/edwin/results/dlr_clustering/images');
%graph_folder = 1;

%data set sizes
exp_settings.training_set_size = 20;
exp_settings.test_data_size = 200;

exp_settings.length_history = 20;
exp_settings.agent_type = 'dlr';

%data generation
exp_settings.means = [-3 -4 -5 -9 8; 3 4 5 9 -8];

%test with sensors all of same quality
exp_settings.deviations = [1 2 1 3 1; 1 1 2 1 3];

exp_settings.p_c1 = 0.5;

exp_settings.noninf_means = 0;
exp_settings.noninf_dev = 5;

%missing labels, corrupted labels etc.

% Here we can overwrite defaults for each sensor. If the vector contains
% fewer elements than there are sensors, the default values will be used.

%p_missing_a = [0.07 0 0.3 0 0 0]; % probability of sensor switching to "I don't know" labels 
% for at least one agent. How do different agents deal with uncertainty
% after unknown labels?
exp_settings.mean_ml_a = [];
exp_settings.dev_ml_a = [];

%p_corrupt_a = [0 0.05 0.1 0 0 0 0 0 0];%0.05; %probability of a single label being flipped for one or more agents.

%flip labels for agents
exp_settings.p_flip_a = [];

%sensors
%p_flip_s = [0.015 0 0 0 0 0 0 0 0 0 0]; %class flips - may have similar effect to switching sensors?

%settings used: 
% inf =5 noninf = 11; max =10 min = 3
% inf = 3 noninf = 8
% inf = 2, noninf = 2; max = 3, min = 1
% inf = 1 noninf = 1

exp_settings.num_inf_sensors = 5;
exp_settings.num_noninf_sensors = 10;

exp_settings.p_switch = 0.01; % probability that a sensor will stop being useful

%agents and sensor allocation
exp_settings.num_agents = 100;

exp_settings.max_sensors = 4;
exp_settings.min_sensors = 2;

%cluster settings
exp_settings.K = 2;
exp_settings.sequence_length = 20;