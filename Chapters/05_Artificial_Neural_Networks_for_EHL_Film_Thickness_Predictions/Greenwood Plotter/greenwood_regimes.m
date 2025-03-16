% clear
% clc

%% Load data

%Input data
%load range_test.mat
% [FILENAME, PATHNAME] = uigetfile('*.mat');
% fullfilepath = strcat(PATHNAME,'/',FILENAME);
% 
% load(fullfilepath);
%Boundaries

function []=greenwood_regimes(A,B)
load Region1s.txt
load Region2s.txt
load Region4s.txt
load Region5s.txt

x1 = Region1s(:,1);
y1 = Region1s(:,2);
x2 = Region2s(:,1);
y2 = Region2s(:,2);
x4 = Region4s(:,1);
y4 = Region4s(:,2);

%% Plotting

% 0.5 blue
% 1.5 red
% 5 green
% 10 magenta
% 15 black

h1 = figure(1);
loglog(x1,y1,'k',x2,y2,'k',x4,y4,'k','LineWidth',4)
%plot(x1,y1,'k',x2,y2,'k',x4,y4,'k','LineWidth',4)
hold on
loglog(B,A,'co','LineStyle', 'none')
title('Boundaries')
hold off

text(2,805,'\bf PE','FontSize',24);
text(0.11, 6.9,'\bf PR','FontSize',24);
text(0.165, 0.31,'\bf IR','FontSize',24);
text(90, 2.7,'\bf IE','FontSize',24);

xlim([10e-2 10e2])
ylim([10e-2 10e3])

set(gca,'FontSize',24)
title ('Greenwood Regimes Line Contact')
xlabel ('$B = (W_{s}^{2}/\eta _{0}uE_{r}R_{r})^\frac{1}{2}$','Interpreter','latex') %https://latex.codecogs.com/eqneditor/editor.php
ylabel ('$A = (\alpha ^{2}W_{s}^{3}/\eta_{0}uR_{r}^{2})^\frac{1}{2}$','Interpreter','latex')

x0=10;
y0=10;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
savefig(h1,'JOUNRAL_AI_Greenwood_Regime_v2')
end

