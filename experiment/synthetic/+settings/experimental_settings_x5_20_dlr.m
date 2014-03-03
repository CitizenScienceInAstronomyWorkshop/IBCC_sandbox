
exp_settings = settings.default_settings();

exp_settings.save_imgs = true

exp_settings.agent_type='dlr';

%data set sizes
exp_settings.training_set_size = 50;
exp_settings.test_data_size = 300;

exp_settings.length_history = 20;

%data generation
exp_settings.means = [-3 -2 -1 3 2; 3 2 1 -3 -2];

%test with sensors all of same quality
exp_settings.deviations = [1 2 1 3 1; 1 2 1 1 3];

exp_settings.p_c1 = 0.5;

exp_settings.noninf_mean = 0;
exp_settings.noninf_dev = 5;

%sensors

%settings used: 
% inf =5 noninf = 11; max =10 min = 3
% inf = 3 noninf = 8
% inf = 2, noninf = 2; max = 3, min = 1
% inf = 1 noninf = 1

exp_settings.num_inf_sensors = 5;
exp_settings.num_noninf_sensors = 5;

exp_settings.p_switch = 0.0; % probability that a sensor will stop being useful
exp_settings.p_missing = 0.00;
exp_settings.p_corrupt = 0.0;
exp_settings.p_flip_label = 0.0;
exp_settings.p_flip_sensor = 0.0;

%agents and sensor allocation
exp_settings.num_agents = 10;

exp_settings.max_sensors = 8;
exp_settings.min_sensors = 3;

%number of clusters and clustering window length
exp_settings.sequence_length = 40;
exp_settings.K = 2;