%% This program retrieves time-domain data from the relevant case folders and generates waterfall plots in frequency domain
clc
clear

%% Creating save file locations
SplitResultsFolder = ("C:\Users\wshq\OneDrive - Loughborough University\PhD R1\7.0 Models\WSHQ_ANN_100_General_MATLAB_16_06_22\matlab\ANN_lubricated_CPU_0.001\Split Results");

%% Import data

for i = [5000:5000:25000]
    
load(['Simulation Files/' num2str(i) 'rpm_ANN.mat'])

%% Bearing 1

% Forces
dispy1_ANN = cell2mat({out.Bearing1Forces.signals(1).values}.');
dispz1_ANN = cell2mat({out.Bearing1Forces.signals(2).values}.');
forcey1_ANN = cell2mat({out.Bearing1Forces.signals(3).values}.');
forcez1_ANN = cell2mat({out.Bearing1Forces.signals(4).values}.');
angdispx1_ANN = cell2mat({out.Bearing1Forces.signals(5).values}.');
elementforcesz1_ANN = permute(cell2mat({out.Bearing1Forces.signals(6).values}.'),[3 2 1]);

%Film
delta1_ANN = permute(cell2mat({out.Bearing1Film.signals(1).values}.'),[3 2 1]); % cell2mat changes from cell array to matlab array, permute changes from 3d to 2d array
film1_ANN = permute(cell2mat({out.Bearing1Film.signals(2).values}.'),[3 2 1]);
W1_ANN = permute(cell2mat({out.Bearing1Film.signals(3).values}.'),[3 2 1]);
entrainment1_ANN = permute(cell2mat({out.Bearing1Film.signals(4).values}.'),[3 2 1]);
def1_ANN = permute(cell2mat({out.Bearing1Film.signals(5).values}.'),[3 2 1]);
angvel1_ANN = permute(cell2mat({out.Bearing1Film.signals(6).values}.'),[3 2 1]);

%Time
t1_ANN =  cell2mat({out.tout}');

%CPU Time
cpu_ANN = cell2mat({out.CPUtime.signals.values}');

%Save results
filename1 = sprintf('%drpm_bearing1_101_0_5um.mat', i);
save (fullfile(SplitResultsFolder, filename1), 'dispy1_ANN', 'dispz1_ANN', 'forcey1_ANN', 'forcez1_ANN', 'angdispx1_ANN', 'elementforcesz1_ANN', 'delta1_ANN', 'film1_ANN', 'W1_ANN', 'entrainment1_ANN', 'def1_ANN', 'angvel1_ANN', 't1_ANN', 'cpu_ANN')


%% Bearing 2

% Forces
dispy2_ANN = cell2mat({out.Bearing2Forces.signals(1).values}.');
dispz2_ANN = cell2mat({out.Bearing2Forces.signals(2).values}.');
forcey2_ANN = cell2mat({out.Bearing2Forces.signals(3).values}.');
forcez2_ANN = cell2mat({out.Bearing2Forces.signals(4).values}.');
angdispx2_ANN = cell2mat({out.Bearing2Forces.signals(5).values}.');
elementforcesz2_ANN = permute(cell2mat({out.Bearing2Forces.signals(6).values}.'),[3 2 1]);

%Film
delta2_ANN = permute(cell2mat({out.Bearing2Film.signals(1).values}.'),[3 2 1]); 
film2_ANN = permute(cell2mat({out.Bearing2Film.signals(2).values}.'),[3 2 1]);
W2_ANN = permute(cell2mat({out.Bearing2Film.signals(3).values}.'),[3 2 1]);
entrainment2_ANN = permute(cell2mat({out.Bearing2Film.signals(4).values}.'),[3 2 1]);
def2_ANN = permute(cell2mat({out.Bearing2Film.signals(5).values}.'),[3 2 1]);
angvel2_ANN = permute(cell2mat({out.Bearing2Film.signals(6).values}.'),[3 2 1]);

%Time
t2_ANN =  cell2mat({out.tout}');

%Save results
filename2 = sprintf('%drpm_bearing2_111_0_5um.mat', i);
save (fullfile(SplitResultsFolder, filename2), 'dispy2_ANN', 'dispz2_ANN', 'forcey2_ANN', 'forcez2_ANN', 'angdispx2_ANN', 'elementforcesz2_ANN', 'delta2_ANN', 'film2_ANN', 'W2_ANN', 'entrainment2_ANN', 'def2_ANN', 'angvel2_ANN', 't2_ANN', 'cpu_ANN')
end