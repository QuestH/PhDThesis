% Results processing for MATLAB/Excite co-simulation of lateral degree of
% freedom Dry bearing model for accelerationof inner race

clear
clc

speedrange = [1000:1000:11000 11250:250:12750 13000:1000:25000];

for i = 1:1:length(speedrange)

speed(i) = speedrange(i);

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

velydry = diff(dispy1_dry)./diff(t1_dry);
velylub = diff(dispy1_lub)./diff(t1_lub);

velydry(end+1:numel(t1_dry))=0;
velylub(end+1:numel(t1_lub))=0;

accelydry = diff(velydry)./diff(t1_dry);
accelylub = diff(velylub)./diff(t1_lub);
 

%% Acceleration Data

acceldatalub = accelylub(minA:maxA-1);
acceldatadry = accelydry(minA:maxA-1);

maxaccellub(i) = max(acceldatalub);
minaccellub(i) = min(acceldatalub);

maxacceldry(i) = max(acceldatadry);
minacceldry(i) = min(acceldatadry);

percentagediff(i) = ((maxacceldry(i)-maxaccellub(i))/maxaccellub(i))*100;


% %% Inner Race Acceleration Envelopes

% figure(1)
% hold on
% plot(speed,maxaccellub,'-ob',speed,maxacceldry,'-or',speed,minaccellub,'--ob',speed,minacceldry,'--or','LineWidth',1)
% set(gca,'FontSize',24)
% title('Inner Race Acceleration Operating Envelope')
% xlabel ('Rotational Speed (rpm)')
% ylabel ('Acceleration (m/s^2)')
% legend ('Lubricated Model','Dry Model')

figure(2)
hold on
plot(speed,percentagediff,'-b','LineWidth',1)
set(gca,'FontSize',24)
title('Dry and Lubricated Percentage Difference')
xlabel ('Rotational Speed (rpm)')
ylabel ('Acceleration (m/s^2)')



end
