clear% clc
% clear all
tic
warning('off')

%% Precalculations
omegaINNERrpm = 1000:1000:25000;
omegaINNERrad = omegaINNERrpm .* ((2*pi)/60);

R11=8.8e-3/2; %roller radius
R21=51.2e-3/2; %race radius 60mm pitch radius
u_current1 =  abs((R21*omegaINNERrad/2)*((R21+2*R11)/(R21+R11))); %entrainment velocity
Rrin1 = 1/(1/R11+1/R21); %reduced radius inner race

% Dimensionalised Parameters
w1 = repmat(3000,size(omegaINNERrad));
u1 =  u_current1;
Rr1 = repmat(Rrin1,size(omegaINNERrad));
Er1 =   repmat(2.3187e11,size(omegaINNERrad));
alpha1 =   repmat(2.1e-8,size(omegaINNERrad));
neta1 =   repmat(0.08,size(omegaINNERrad));
length1 =   repmat(15e-3,size(omegaINNERrad)); % [m] roller length

% Dimensionless Parameters
Ue_in = (neta1.*u1)./(Er1.*Rr1);
Ge = Er1.*alpha1;
We_in = w1./(Er1.*Rr1.*length1);

data_points = length(omegaINNERrad);

%% EHL distribution X-domain 
nNodes = 200;
inlet = 4.5;
outlet = 1.5;
io_dis = inlet+outlet;
X=-inlet:io_dis/(nNodes-1):outlet; % dimensionless coordinate

%% Create the output parameters for storage
output_film = zeros(data_points,1);
tic

%% begin parfor loop to generate data
parfor i = 1 :data_points
    % select input parameters for loop
    w0 = w1(i); 
    u0 = u1(i);
    Rr0 = Rr1(i);
    Er0 = Er1(i);
    alpha0 = alpha1(i);
    neta0 = neta1(i);
    length0 = length1(i);
    Ue = Ue_in(i);
    Ge1 = Ge(i);
    We = We_in(i)
    
    % call function to compute 1D EHL
    [~,~,hmin,hc1,~,~,~,~,~,~,~]= OneD_EHL(u0,w0,Rr0,length0,alpha0, neta0, nNodes, X, Er0 );     
    % store output data
    output_film(i,:) = hc1;
    filmNEW(i) = Rr0*3.06*(Ue^0.69)*(Ge1^0.56)*(We^-0.1);
    
    percentage_error(i) = abs((filmNEW(i) - output_film(i,:)) / output_film(i,:)) * 100;
    
end
toc
%% save workspace
save output_file_central.mat

% Generate plots
figure;
yyaxis left % Left y-axis for film thickness
plot(omegaINNERrpm, filmNEW * 1e6, '--k', 'LineWidth', 1); % Analytical (dashed black)
hold on;
plot(omegaINNERrpm, output_film * 1e6, '-k', 'LineWidth', 1); % Numerical (solid black)
ylabel('Film Thickness / \mum', 'Color', 'k'); % Black y-axis label
ylim([0 5]); % Adjust as needed
set(gca, 'YColor', 'k'); % Set left y-axis color to black

% Formatting for left y-axis
set(gca, 'FontSize', 24);
xlabel('Rotational Speed / rpm');
title('EHL Film Thickness Analytical vs Numerical');

% Right y-axis for percentage error
yyaxis right
plot(omegaINNERrpm, percentage_error, 'Color', '#D95319', 'LineWidth', 1.5); % Percentage error in blue
ylabel('Percentage Difference / %', 'Color', '#D95319'); % Blue y-axis label
ylim([0 25]); % Adjust as needed
set(gca, 'YColor', '#D95319'); % Set right y-axis color to blue

% Add legend
legend('Analytical', 'Numerical', 'Percentage Difference');
