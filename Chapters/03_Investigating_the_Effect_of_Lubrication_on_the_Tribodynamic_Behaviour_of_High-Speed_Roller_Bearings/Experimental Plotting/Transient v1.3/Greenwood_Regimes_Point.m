clear all
clc

load Greenwood_All.mat

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

% figure(1)
% loglog(g3,g2,gg3,gg2,ggg3,ggg2,x_flat,y_flat,x2_flat,y2_flat)
% title('Boundaries')

% figure(2)
% loglog(g3,g2,gg3,gg2,ggg3,ggg2,real(gggg3),real(gggg2),x_flat,y_flat,x2_flat,y2_flat'.');
% title('Point Contact Formulation')
% xlim([1e-5 1e3])
% ylim([1e0 1e6])

figure(3)
loglog(g3,g2,'k',gg3,gg2,'k',ggg3,ggg2,'k',real(ge),real(gv),'b',x_flat,y_flat,'k',x2_flat,y2_flat,'k','LineWidth',4)
set(gca,'FontSize',24)
    title ('Line Contact Formulation')
    xlabel ('G_e')
    ylabel ('G_v')

xlim([1e-6 1e4])
ylim([1e-2 1e6])



% 
% figure(4)
% loglog(g3,g2,gg3,gg2,ggg3,ggg2,real(ge),real(gv),'.')
% hold on
% plot3(real(ge),real(gv),landa_out)

%% Contour plot

a = randi(9, 10, 3);
x = real(ggg3);
y = real(ggg2);
z = landa_out;
xv = linspace(min(x), max(x));
yv = linspace(min(y), max(y));
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);
figure (4)
contourf(X, Y, Z)
colorbar

% plot(gggg3,gggg2);