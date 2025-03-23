%% LHS sampling - function to produce a sample of EHL input parameters
load Bearing_Test_Data_FMBD_20m_s.mat
%load spur_gear_test_data.mat
% assign individual parameters from test data
w_i = test_data(:,1);
u_e_i =  test_data(:,2);
Rr_i = test_data(:,3);
Er2_i =   test_data(:,4) ;
alpha1_i =   test_data(:,5) ;
neta_i =   test_data(:,6) ;
pois_i =   test_data(:,7) ;
density_i =  test_data(:,8) ;
length2_i =   test_data(:,9) ;

% 
B = ((w_i./length2_i).^2./(neta_i.*u_e_i.*Er2_i.*Rr_i)).^(1/2);
A = ((alpha1_i.^2.*(w_i./length2_i).^3)./(neta_i.*u_e_i.*(Rr_i.^2))).^(1/2);

greenwood_regimes(A,B)