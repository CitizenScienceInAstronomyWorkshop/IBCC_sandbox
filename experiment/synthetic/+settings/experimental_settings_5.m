save_imgs = true
%data set sizes
training_set_size = 20;
test_data_size = 500;

length_history = 20;

%data generation
means = [-3 -2 -1 3 2; 3 2 1 -3 -2];

%test with sensors all of same quality
deviations = [1 2 1 3 1; 1 1 2 1 3];

p_c1 = 0.5;

noninf_means = 0;
noninf_dev = 10;

%sensors

%settings used: 
% inf =5 noninf = 11; max =10 min = 3
% inf = 3 noninf = 8
% inf = 2, noninf = 2; max = 3, min = 1
% inf = 1 noninf = 1

num_inf_sensors = 1;
num_noninf_sensors = 1;
num_sensors = num_inf_sensors + num_noninf_sensors;

p_switch = 0.017; % probability that a sensor will stop being useful

%agents and sensor allocation
num_agents = 5;

max_sensors = 1;
min_sensors = 1;