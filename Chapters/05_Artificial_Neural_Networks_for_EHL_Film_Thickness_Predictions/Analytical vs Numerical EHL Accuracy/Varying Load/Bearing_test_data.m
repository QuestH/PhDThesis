clear
clc

%% Spur Gear - ANN Application
%% EHL distribution X-domain
nNodes = 200;
inlet = 4.5;
outlet = 1.5;
io_dis = inlet+outlet;
X=-inlet:io_dis/(nNodes-1):outlet; % dimensionless coordinate
disp_update = 'test on application data...';
disp(disp_update);

% load test data
%load test_data.mat
load Bearing_Test_Data_FMBD_20m_s.mat
load PS1PS2.mat
% create dimensional data
deg = test_data(:,10);
[test_data] = mapminmax('apply',test_data(:,1:9)',PS1); %normalised test data
[test_dim] = mapminmax('reverse',test_data,PS1); %dimensionalised test data
test_dim = test_dim';

% assign individual parameters from test data
w_i = test_dim(:,1);
u_e_i =  test_dim(:,2);
Rr_i = test_dim(:,3);
Er2_i =   test_dim(:,4) ;
alpha1_i =   test_dim(:,5) ;
neta_i =   test_dim(:,6) ;
pois_i =   test_dim(:,7) ;
density_i =  test_dim(:,8) ;
length2_i =   test_dim(:,9) ;

% compute 1-D EHL distribution using numerical code
disp_update = 'time to complete a single numerical 1-D EHL distribution...';
disp(disp_update);
test_array_length = size(test_dim,1);
tic
for i = 1:test_array_length
    
    W1(i) =w_i(i)/(length2_i(i)*Er2_i(i)*Rr_i(i)); %dimensionless load
    PH(i)=Er2_i(i)*sqrt(W1(i)/(2*pi)); %dimensionless maximum Hertzian pressure?
    U(i)=u_e_i(i)*neta_i(i)./(Rr_i(i)*Er2_i(i)); %dimensionless speed of entraining motion
    G(i)=alpha1_i(i)*Er2_i(i);
    a(i)=Rr_i(i)*sqrt(8*W1(i)/pi); %dimension contact half width
    %-------------------------------------------------------------------------%
    %% central film thickness
    %-------------------------------------------------------------------------%
    Hc_ext(i) =3.09*(pi/(8*W1(i)))*G(i)^0.56*U(i)^0.69*W1(i)^(-0.1);
    Hc_ext_dim(i) = Hc_ext(i)*(a(i)^2)/Rr_i(i);
%     Hc_ext(i) =1.95*(alpha1(i)*neta(i)*u_e(i)/Rr(i))^(8/11) * (Er2*Rr(i)/(w(i)/length2))^(1/11);
%     Hc_ext_dim(i) = Hc_ext(i)*Rr(i);
    
    %-------------------------------------------------------------------------%
end
time_ext = toc;
tic
parfor i = 1:test_array_length
    w2 = w_i(i);
    u2 = u_e_i(i);
    Rr2 = Rr_i(i);
    alpha2 = alpha1_i(i);
    neta2  = neta_i(i);
    length2 = length2_i(i);
    Er2 = Er2_i(i);
    
    [Pht,ht,hmin,hc1,a,Pc1,Pmax,F_Ph1,Hmin,Etta, Ro]= OneD_EHL(u2,w2,Rr2,length2,alpha2, neta2, nNodes, X, Er2 );
    hc_num(i) = hc1;

end
time_num = toc;
disp_update = 'time to complete a single ANN 1-D EHL distribution...';
disp(disp_update);


tic
% compute 1-D EHL distribution points using ANN
for i = 1:test_array_length
    w0 = w_i(i);
    u0 = u_e_i(i);
    Rr0 = Rr_i(i);
    Er0 = Er2_i(i) ;
    alpha0 = alpha1_i(i);
    neta0  = neta_i(i);
    pois0 = pois_i(i);
    density0 = density_i(i);
    length0 = length2_i(i);
    
    in_vars = [w0 u0 Rr0 Er0 alpha0 neta0 pois0 density0 length0];
    [in_vars_norm] = mapminmax('apply',in_vars',PS1) ; %normalise input data by mapping between -1 and 1
        
    % compute 1-D EHL distribution points using ANN
    [h_c] = JOURNALNetworkFunction(in_vars_norm);
    h_ann = h_c;
    
    % dimensionalise 1-D EHL distribution points computed by ANN
    [h_ann_dim] = mapminmax('reverse',h_ann,PS2);
    
    hc_ann(i) = h_ann_dim;
    
    % dimensionalise 1-D EHL distribution points computed by ANN
end
    
    time_ann = toc;
    
    %% Plot comparisons
    disp_update = 'plotting Numerical versus ANN computed 1-D EHL distribution...';
    disp(disp_update);
    figure
    hold on
    plot(ang, hc_ann*1e6, 'b-') %plot(deg,hc_ann,'b-')
    plot(ang, hc_num*1e6, 'r-')
    plot(ang, Hc_ext_dim*1e6, 'g-')
    legend('ANN','Numerical','Analytical')
    xlabel('Deg [/degree]')
    ylabel('Film thickness [m]')
    hold off
    
    
    % %% Generate plots for EHL Film Thickness ANN vs Analytical vs Numerical
figure;
plot(deg, Hc_ext_dim*1e6, ':k', 'LineWidth', 4); % Analytical (dashed black)
hold on;
plot(deg, hc_num*1e6, 'Color', '#0072BD', 'LineWidth', 4); % Numerical (solid black)
hold on;
plot(deg, hc_ann*1e6, '--', 'Color', '#D95319', 'LineWidth', 4); % Numerical (solid black)
ylabel('Film Thickness / \mum', 'Color', 'k'); % Black y-axis label
ylim([2 2.5]); % Adjust as needed

% Formatting for left y-axis
set(gca, 'FontSize', 24);
xlabel('Cage angular displacement / $^{\circ}$');
title('EHL Film Thickness ANN vs Analytical vs Numerical');
legend('Analytical', 'Numerical', 'ANN');

%Calculate MSE:
ANNmse = immse(hc_ann*1e6,hc_num*1e6);
AnalytMSE = immse(Hc_ext_dim*1e6,hc_num*1e6);

