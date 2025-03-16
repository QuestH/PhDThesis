function [input_vars,A,B] = LHS_data(data_points)
%% LHS sampling - function to produce a sample of EHL input parameters

% latin hyper cube sampling function
no_var = 9;
% load [N]
min_force = 150;
max_force = 5000;
% speed - 1,2 [m/s]
min_speed = 0.6;
max_speed = 30;
% youngs [Pa]
min_youngs = 200e9;
max_youngs = 250e9;
% poissons ratio [-]
min_poisr = 0.30;
max_poisr = 0.35;
% reduced radius [-]
min_rreduced = 1.0e-04;
max_rreduced = 0.02;
% length [m]
min_length = 0.001;
max_length = 0.05;
% density [kg/m^3]
min_density = 7750;
max_density = 8050;
% viscosity [Pas]
min_visc = 0.0005;
max_visc = 0.1;
% pressure-viscosity term [1/Pa]
min_alpha_p = 1.0e-8;
max_alpha_p = 3.0e-8;

%% Max Marian - input
% https://link.springer.com/content/pdf/10.1007/s40544-022-0641-6.pdf
% load [N]
% min_force = 750;
% max_force = 1500;
% % speed [m/s]
% min_speed = 0.025;
% max_speed = 0.4;
% % youngs [Pa]
% min_youngs = 200e9;
% max_youngs = 440e9;
% % poissons ratio [-]
% min_poisr = 0.30;
% max_poisr = 0.35;
% % reduced radius [-]
% min_rreduced = 0.0075;
% max_rreduced = 0.02;
% % length [m]
% min_length = 0.00025;
% max_length = 0.01;
% % density [kg/m^3]
% min_density = 8500;
% max_density = 11000;
% % viscosity [Pas]
% min_visc = 0.005;
% max_visc = 0.05;
% % pressure-viscosity term [1/Pa]
% min_alpha_p = 1.25e-8;
% max_alpha_p = 2.5e-8;

load Region1s.txt

%load Region1s.txt
x1 = Region1s(:,1);
y1 = Region1s(:,2);
idx=0;

% convert to variables
X = lhsdesign(data_points,no_var);
force = min_force + (max_force-min_force) * X(:,1);
speed = min_speed + (max_speed-min_speed) * X(:,2);
youngs = min_youngs + (max_youngs-min_youngs)* X(:,3);
poisr = min_poisr + (max_poisr -min_poisr)* X(:,4);
rreduced = min_rreduced + (max_rreduced -min_rreduced)* X(:,5);
length = min_length + (max_length-min_length) * X(:,6);
density = min_density + (max_density-min_density) * X(:,7);
visc = min_visc + (max_visc-min_visc) * X(:,8);
alpha_p = min_alpha_p + (max_alpha_p-min_alpha_p) * X(:,9);

B = ((force./length).^2./(visc.*speed.*youngs.*rreduced)).^(1/2);
A = ((alpha_p.^2.*(force./length).^3)./(visc.*speed.*(rreduced.^2))).^(1/2);
p0 = (youngs.*force/pi./length./rreduced).^0.5;
idx =0;

% for i = 1:size(A,1)
%         check_1 = any(and(B(i)> x1, A(i)> y1));
%         check_2 = and(check_1,p0(i)<3.5e9);
%         check_3 = and(check_2,p0(i)>600e6);
%     while check_3 == 0
%         Xcorr = lhsdesign(1,no_var);
%         force(i) = min_force + (max_force-min_force)* Xcorr(1,1);
%         speed(i) = min_speed + (max_speed-min_speed) * Xcorr(1,2);
%         youngs(i) = min_youngs + (max_youngs-min_youngs)* Xcorr(1,3);
%         poisr(i) = min_poisr + (max_poisr-min_poisr)* Xcorr(1,4);
%         rreduced(i) = min_rreduced + (max_rreduced-min_rreduced)* Xcorr(1,5);
%         length(i) = min_length + (max_length-min_length)* Xcorr(1,6);
%         density(i) = min_density + (max_density-min_density)* Xcorr(1,7);
%         visc(i) = min_visc + (max_visc-min_visc)*Xcorr(1,8);
%         alpha_p(i) = min_alpha_p + (max_alpha_p-min_alpha_p)* Xcorr(1,9);
%         B(i) = ((force(i)./length(i)).^2./(visc(i).*speed(i).*youngs(i).*rreduced(i))).^(1/2);
%         A(i) = ((alpha_p(i).^2.*(force(i)./length(i)).^3)./(visc(i).*speed(i).*(rreduced(i).^2))).^(1/2);
%         p0(i) = (youngs(i).*force(i)./length(i)/pi./rreduced(i)).^0.5;
%         X(i,:) = Xcorr ;
%         idx = 1+idx;
%         
%         check_1 = any(and(B(i)> x1, A(i)> y1));
%         check_2 = and(check_1,p0(i)<3.5e9);
%         check_3 = and(check_2,p0(i)>600e6);
%          
%     end
% end
input_vars = [force,speed,rreduced,youngs,alpha_p,visc,poisr,density,length];
end