clear all
clc

%% Load data

load Greenwood_All.mat

% For this code, I had arrays with multiple values at different time steps
% for each variable saved in file Greenwood_All.mat 

% You will need to upload your own data for the variables in the next
% section of the code, or define them in the input variable section



%% Input variables

% RSTAR = Ry/Rx for rolling element
% neta0 = atmospheric viscosity
% U = speed of entraining motion
% Er = reduced modulus of elasticity of contacting bodies
% Rrout = reduced radius of outer race and roller contact
% alpha = pressure-viscosity coefficient
% Wt_out = outer race contact force
% Vn = Dimensionless speed parameter
% Gn = Dimensionless equivalent geometry
% Wn = Dimensionless load parameter

RSTAR =1; 
Vn = neta0.*U./(Er*Rrout); 
Gn = Er*alpha;
Wn = Wt_out./(Er*(Rrout^2)); 



%% Start of Code

FI=((1+(2/(3*RSTAR)))^(-1));
H=32*FI*(RSTAR^2)*((1.68*atan(RSTAR/2)+0.13)^2); %(7.21)

%Step 3 - Fixes point A
g2(1)=(H/(0.104*(1-exp(-0.7*(RSTAR^0.64)))))^(1/0.67);%(7.23) 
g3(1)=0;

%Step 5 - Fixes point C
g2(2)=g2(1);
g3(2)=(H/(2.419*(g2(2)^0.49)*(1-exp(-0.67*(RSTAR^0.66)))))^(1/0.17);%(7.20)

%Step 6 - Fixes point E
g2(4)=0;
g3(4)=(H/(8.7*(1-0.85*exp(-0.32*(RSTAR^0.64)))))^(1/0.67);%(7.24)

%Step 7 - Fixes point D
g3(3)=g3(4);
g2(3)=(H/(2.419*(g3(3)^0.17)*(1-exp(-0.67*(RSTAR^0.66)))))^(1/0.49);%(7.20)

H=3*H; %Increase H

gg2(1)=g2(2);
gg3(1)=g3(2);

gg2(2)=(H/(0.104*(1-exp(-0.7*(RSTAR^0.64)))))^(1/0.67);%(7.23)
gg3(2)=(H/(2.419*(gg2(2)^0.49)*(1-exp(-0.67*(RSTAR^0.66)))))^(1/0.17);%(7.20)

ggg2(1)=g2(3);
ggg3(1)=g3(3);

ggg3(2)=(H/(8.7*(1-0.85*exp(-0.32*(RSTAR^0.64)))))^(1/0.67);%(7.24)
ggg2(2)=(H/(2.419*(ggg3(2)^0.17)*(1-exp(-0.67*(RSTAR^0.66)))))^(1/0.49);%(7.20)


%% Overlaying results on axes
gggg2=4*(Gn.*(Wn.^3))./(Vn.^2); %point contact 
gggg3=0.85*(Wn.^(8/3))./(Vn.^2); %point contact

ge = Wn./(Vn.^0.5); %line contact 
gv = Wn.^(3/2)*Gn./(Vn.^0.5); %line contact 

x_flat = [g3(2) 1e-5];
y_flat = [g2(2) g2(2)];

x2_flat = [ggg3(1) ggg3(1)];
y2_flat = [ggg2(1) 1e0];



%% Plotting

figure(1)
loglog(g3,g2,'k',gg3,gg2,'k',ggg3,ggg2,'k',real(ge),real(gv),'b',x_flat,y_flat,'k',x2_flat,y2_flat,'k','LineWidth',4)
set(gca,'FontSize',24)
    title ('Line Contact Formulation')
    xlabel ('G_e')
    ylabel ('G_v')

xlim([1e-6 1e4])
ylim([1e-2 1e6])
