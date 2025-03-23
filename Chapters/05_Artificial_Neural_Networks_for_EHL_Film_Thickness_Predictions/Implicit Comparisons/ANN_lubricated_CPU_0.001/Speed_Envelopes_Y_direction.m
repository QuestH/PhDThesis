% Results processing for MATLAB/Excite co-simulation of lateral degree of
% freedom Dry bearing model for accelerationof inner race

clear
clc

for i = 1:25

speed(i) = i*1000;

%% Inner Race Displacement Data
load(['Split Results/' num2str(speed(i)) 'rpm_bearing1_101.mat']);
load(['../DryNP175um\Split Results/' num2str(speed(i)) 'rpm_bearing1_101.mat']);

A = length(dispy1_lub);
B = length(dispy1_dry);
maxA = min(A,B);
minA = maxA-1e5;

dispdatalub = dispy1_lub(minA:maxA-1);
dispdatadry = dispy1_dry(minA:maxA-1);

maxdisplub(i) = max(dispdatalub);
mindisplub(i) = min(dispdatalub);

maxdispdry(i) = max(dispdatadry);
mindispdry(i) = min(dispdatadry);

%% Contact Stiffness Data

%Stiffness Calculation
stiffnesslub = W1_lub./(0.5*delta1_lub); %calculate linear stiffness
stiffnessdry = W1_dry./(0.5*delta1_dry); %calculate linear stiffness

stiffdatalub = stiffnesslub(minA:maxA-1);
stiffdatadry = stiffnessdry(10*minA:10*maxA-1);

maxstifflub(i) = max(stiffdatalub);
minstifflub(i) = min(stiffdatalub);

maxstiffdry(i) = max(stiffdatadry);
minstiffdry(i) = 0;


%% Contact Deformation Data
defdatalub = def1_lub(minA:maxA-1);
defdatadry = def1_dry(minA:maxA-1);

maxdeflub(i) = max(defdatalub);
mindeflub(i) = min(defdatalub);

maxdefdry(i) = max(defdatadry);
mindefdry(i) = min(defdatadry);


%% Contact Force Data
cforcedatalub = W1_lub(minA:maxA-1);
cforcedatadry = W1_dry(minA:maxA-1);

maxcforcelub(i) = max(cforcedatalub);
mincforcelub(i) = min(cforcedatalub);

maxcforcedry(i) = max(cforcedatadry);
mincforcedry(i) = min(cforcedatadry);


%% Bearing Centre Stiffness
forcedatalub = forcey1_lub(minA:maxA-1);
forcedatadry = forcey1_dry(minA:maxA-1);

centrestifflub = -forcedatalub./dispdatalub;
centrestiffdry = -forcedatadry./dispdatadry;

maxcentrestifflub(i) = max(centrestifflub);
mincentrestifflub(i) = min(centrestifflub);

maxcentrestiffdry(i) = max(centrestiffdry);
mincentrestiffdry(i) = min(centrestiffdry);

%% Inner Race Force Data

maxforceylub(i) = max(forcedatalub);
minforceylub(i) = min(forcedatalub);

maxforceydry(i) = max(forcedatadry);
minforceydry(i) = min(forcedatadry);

%% EHL Film Stiffness
film_lub = film1_lub(minA:maxA-1);

EHL_stiff(i) = cforcedatalub/film_lub;

%% EHL Film Thickness
film_thick_data = film1_lub(minA:maxA-1);
filmthick(i) = max(film_thick_data);

% %% Inner Race Force Envelopes
hold on
figure(6)
plot(speed,maxforceylub,'-ob',speed,maxforceydry,'-or',speed,minforceylub,'--ob',speed,minforceydry,'--or','LineWidth',1)
set(gca,'FontSize',24)
title('Inner Race Force Operating Envelope')
xlabel ('Rotational Speed (rpm)')
ylabel ('Force (N)')
legend ('Lubricated Model','Dry Model')



end