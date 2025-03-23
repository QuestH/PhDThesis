% Results processing for MATLAB/Excite co-simulation of lateral degree of
% freedom Dry bearing model for acceleration of inner race

clear
clc

for i = 1:25

speed(i) = i*1000;

%% Inner Race Displacement Data
load(['Split Results/' num2str(speed(i)) 'rpm_bearing1_101.mat']);
load(['../DryNP175um\Split Results/' num2str(speed(i)) 'rpm_bearing1_101.mat']);

A = length(dispz1_lub);
B = length(dispz1_dry);
maxA = min(A,B);
minA = maxA-1e5;

dispdatalub = dispz1_lub(minA:maxA-1);
dispdatadry = dispz1_dry(minA:maxA-1);

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
forcedatalub = forcez1_lub(minA:maxA-1);
forcedatadry = forcez1_dry(minA:maxA-1);

centrestifflub = -forcedatalub./dispdatalub;
centrestiffdry = -forcedatadry./dispdatadry;

maxcentrestifflub(i) = max(centrestifflub);
mincentrestifflub(i) = min(centrestifflub);

maxcentrestiffdry(i) = max(centrestiffdry);
mincentrestiffdry(i) = min(centrestiffdry);

%% Inner Race Force Data

maxforcezlub(i) = max(forcedatalub);
minforcezlub(i) = min(forcedatalub);

maxforcezdry(i) = max(forcedatadry);
minforcezdry(i) = min(forcedatadry);

%% EHL Film Stiffness
film_lub = film1_lub(minA:maxA-1);

EHL_stiff(i) = cforcedatalub/film_lub;

%% EHL Film Thickness
film_thick_data = film1_lub(minA:maxA-1);
filmthick(i) = max(film_thick_data);

% %% Inner Race Force Envelopes
hold on
figure(6)
plot(speed,maxforcezlub,'-ob',speed,maxforcezdry,'-or',speed,minforcezlub,'--ob',speed,minforcezdry,'--or','LineWidth',1)
set(gca,'FontSize',24)
title('Inner Race Force Operating Envelope')
xlabel ('Rotational Speed (rpm)')
ylabel ('Force (N)')
legend ('Lubricated Model','Dry Model')



end


%% Plot figures


% %% Inner Race Displacement Envelopes
% figure(1)
% plot(speed,maxdisplub,'-ob',speed,maxdispdry,'-or',speed,mindisplub,'--ob',speed,mindispdry,'--or','LineWidth',1)
% set(gca,'FontSize',24)
% title('Inner Race Displacement Operating Envelope')
% xlabel ('Rotational Speed (rpm)')
% ylabel ('Displacement (m)')
% legend ('Lubricated Model','Dry Model')
% ylim([0 10e-7])
% 
%% Contact Stiffness Envelopes
% figure(2)
% plot(speed,maxstifflub,'-ob',speed,maxstiffdry,'-or',speed,minstifflub,'--ob',speed,minstiffdry,'--or','LineWidth',1)
% set(gca,'FontSize',24)
% title('Contact Stiffness Operating Envelope')
% xlabel ('Rotational Speed (rpm)')
% ylabel ('Stiffness (N/m)')
% legend ('Lubricated Model','Dry Model')
% ylim([0 3e8])
% 
% %% Contact Deformation Envelopes
% figure(3)
% plot(speed,maxdeflub,'-ob',speed,maxdefdry,'-or',speed,mindeflub,'--ob',speed,mindefdry,'--or','LineWidth',1)
% set(gca,'FontSize',24)
% title('Contact Deformation Operating Envelope')
% xlabel ('Rotational Speed (rpm)')
% ylabel ('Deformation (m)')
% legend ('Lubricated Model','Dry Model')
% 
%% Contact Force Envelopes
% figure(4)
% plot(speed,maxcforcelub,'-ob',speed,maxcforcedry,'-or',speed,mincforcelub,'--ob',speed,mincforcedry,'--or','LineWidth',1)
% set(gca,'FontSize',24)
% title('Contact Force Operating Envelope')
% xlabel ('Rotational Speed (rpm)')
% ylabel ('Force (N)')
% legend ('Lubricated Model','Dry Model')
% 
% %% Centre Stiffness
% figure(5)
% plot(speed,maxcentrestifflub,'-ob',speed,maxcentrestiffdry,'-or',speed,mincentrestifflub,'--ob',speed,mincentrestiffdry,'--or','LineWidth',1)
% set(gca,'FontSize',24)
% title('Total Bearing Stiffness Operating Envelope')
% xlabel ('Rotational Speed (rpm)')
% ylabel ('Stiffness (N/m)')
% legend ('Lubricated Model','Dry Model')
% ylim([0 2e9])


%Individual Rolling Element Stiffness and Film Thickness

% figure(5)
% left_color = [0 0 0];
% right_color = [0 0 1];
% set(figure(5),'defaultAxesColorOrder',[left_color; right_color]);
% ax = gca;
% ax.FontSize = 24;
% 
% yyaxis left
% plot(speed,maxstifflub,'LineWidth',1)
% ylabel('Stiffness (N/m)','color','k')
% 
% yyaxis right
% plot(speed,filmthick/1e-6,'--b','LineWidth',1)
% ylabel('Film Thickness (\mum)','color', 'b')
% 
% title([num2str(speed),' rpm Lubricated - Contact Force and Corresponding Film Thickness'])
% ytickformat('%.1f')
% xlabel('Time (s)')
% % xlim ([1 1.2])
% ytickformat('%.1f')



