% Results processing for MATLAB/Excite co-simulation of lateral degree of
% freedom Lubricated bearing model

clear
clc

for i = 20000
%% Choose speed for results

%% Load input data
load(['Split Results/' num2str(i) 'rpm_bearing1_101_0_5um.mat'])

%% Calculate cage angles and speeed
ang1_ANN = 0.4*angdispx1_ANN; %cage angular position
elementforcesy1_ANN = W1_ANN.*sin(ang1_ANN);
                                
speed = round((60/(2*pi))*((max(angdispx1_ANN)-min(angdispx1_ANN))/(max(t1_ANN)-min(t1_ANN)))); %calculate speed in RPM from angular displacement and time

%Creating upper and lower for plots based on one orbit of cage
incr = length(ang1_ANN)/max(ang1_ANN);
lower = 400000; %round((max(ang1)-2*pi)*incr); %lower limit of plot /rad
upper = 500001; %round((lower + 2*pi*incr));

%Stiffness Calculation
stiffness = W1_ANN./(0.5*delta1_ANN); %calculate linear stiffness

delta_r = dispz1_ANN.*cos(ang1_ANN)+dispy1_ANN.*sin(ang1_ANN);

if i == 5000
change = 57299;
% Copy all values from film1_ANN to film1_ANN_new
film1_ANN = film1_ANN(:,1);
film1_ANN_new = film1_ANN;
% Adjust values from change index where computation methodology switches
film1_ANN_new(change:end) = film1_ANN(change:end) .* 1.108;
%flip sinusoid around mean (load/film response)
film1_ANN_new_switch = film1_ANN_new;
mean = mean(film1_ANN_new(change:end));
diff = film1_ANN_new(change:end)- mean;
film1_ANN_new_switch(change:end) = film1_ANN_new(change:end)- 2*diff;

film1_ANN_new_switch(change-1:change+1) = NaN;
film1_ANN_new_switch(1:10) = NaN;

elseif i == 10000
%% 10000
change = 28651;
% Copy all values from film1_ANN to film1_ANN_new
film1_ANN = film1_ANN(:,1);
film1_ANN_new = film1_ANN;
% Adjust values from change index where computation methodology switches
film1_ANN_new(change:end) = film1_ANN(change:end) .* 1.095;
%flip sinusoid around mean (load/film response)
film1_ANN_new_switch = film1_ANN_new;
mean = mean(film1_ANN_new(change:end));
diff = film1_ANN_new(change:end)- mean;
film1_ANN_new_switch(change:end) = film1_ANN_new(change:end)- 2*diff;

film1_ANN_new_switch(change-1:change+1) = NaN;
film1_ANN_new_switch(1:10) = NaN;

elseif i == 15000
%% 15000
change = 19100;
% Copy all values from film1_ANN to film1_ANN_new
film1_ANN = film1_ANN(:,1);
film1_ANN_new = film1_ANN;
% Adjust values from change index where computation methodology switches
film1_ANN_new(change:end) = film1_ANN(change:end) .* 1.079;
%flip sinusoid around mean (load/film response)
film1_ANN_new_switch = film1_ANN_new;
mean = mean(film1_ANN_new(change:end));
diff = film1_ANN_new(change:end)- mean;
film1_ANN_new_switch(change:end) = film1_ANN_new(change:end)- 2*diff;

film1_ANN_new_switch(change-1:change+1) = NaN;
film1_ANN_new_switch(1:10) = NaN;

elseif i == 20000
%% 20000
change = 14325;
% Copy all values from film1_ANN to film1_ANN_new
film1_ANN = film1_ANN(:,1);
film1_ANN_new = film1_ANN;
% Adjust values from change index where computation methodology switches
film1_ANN_new(change:end) = film1_ANN(change:end) .* 1.055;
%flip sinusoid around mean (load/film response)
film1_ANN_new_switch = film1_ANN_new;
mean = mean(film1_ANN_new(change:end));
diff = film1_ANN_new(change:end)- mean;
film1_ANN_new_switch(change:end) = film1_ANN_new(change:end)- 2*diff;

film1_ANN_new_switch(change-1:change+1) = NaN;
film1_ANN_new_switch(1:10) = NaN;
end

plot(-angdispx1_ANN*(180/pi),film1_ANN_new_switch(:,1)*1e6)

hold on

end