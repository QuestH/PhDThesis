%% Main Script for 1D EHL Data Generation and ANN Training and Function Generator

% Author: Jack Walker
% 4/12/2022
% Citation:
% Walker, J., Questa, H., Raman, A. et al.
% Application of Tribological Artificial Neural Networks in Machine Elements.
% Tribol Lett 71, 3 (2023).
% https://doi.org/10.1007/s11249-022-01673-5

%% Required functions:
% OneD_EHL.m
% defcoeff.m
% deformation.m
% film_thickness.m
% LHS_data.m
% greenwood_regimes.m
% COEFF.m

%% Data files required:
% Region1s.txt
% Region2s.txt
% Region4s.txt
% Region5s.txt
%% -----------------------------BEGIN SCRIPT ----------------------------- %%
% clear
% clear all
tic
warning('off')
% -------------------------------------------------------------------------- %
%% Generate data points based in Latin Hyper Cube Sampling
data_points = 5000;

disp_update = 'generate datapoints based on Latin Hypercube Sampling...';
disp(disp_update);

[input_vars,A,B] = LHS_data(data_points);

disp_update = 'plotting Greenwood regime...';
disp(disp_update);

greenwood_regimes(A,B);

disp_update = 'setting up  1-D EHL data generation...';
disp(disp_update);

w1 = input_vars(:,1);
u1 =  input_vars(:,2);
Rr1 = input_vars(:,3);
Er1 =   input_vars(:,4) ;
alpha1 =   input_vars(:,5) ;
neta1 =   input_vars(:,6) ;
pois1 =   input_vars(:,7) ;
density1 =   input_vars(:,8) ;
length1 =   input_vars(:,9) ;

%% EHL distribution X-domain
nNodes = 200;
inlet = 4.5;
outlet = 1.5;
io_dis = inlet+outlet;
X=-inlet:io_dis/(nNodes-1):outlet; % dimensionless coordinate

%% Create the output parameters for storage
output_film = zeros(data_points,1);
tic

% %% begin parfor loop to generate data
% disp_update = 'beginning 1-D EHL data generation...';
% disp(disp_update);
% for i = 1 :data_points
%     % select input parameters for loop
%     w0 = w1(i);
%     u0 = u1(i);
%     Rr0 = Rr1(i);
%     Er0 = Er1(i);
%     alpha0 = alpha1(i);
%     neta0 = neta1(i);
%     pois0 = pois1(i);
%     density0 = density1(i);
%     length0 = length1(i);
%     % call function to compute 1D EHL
%     [~,~,hmin,hc1,~,~,~,~,~,~,~]= OneD_EHL(u0,w0,Rr0,length0,alpha0, neta0, nNodes, X, Er0 );
%     % store output data
%     output_film(i,:) = hc1;
% end
% toc
% 
% disp_update = '1-D EHL data generation completed...';
% disp(disp_update);
% %% save workspace
% save output_file_JOURNAL_EDITION_central.mat
% 
% %% data formatting for neural network
% % append domain position to the input data
% 
% disp_update = 'format 1-D EHL data for ANN use...';
% disp(disp_update);
% in_vars = input_vars;
% out_vars = output_film;
% 
% %% Normalise data using mapminmax fnc
% disp_update = 'normalising 1-D EHL data using mapminmax function...';
% disp(disp_update);
% 
% [in_vars_norm,PS1] = mapminmax(in_vars') ; %normalise input data by mapping between -1 and 1
% [out_vars_norm,PS2] = mapminmax(out_vars') ; %normalise output data by mapping between -1 and 1
% in_vars_norm = in_vars_norm';
% out_vars_norm = out_vars_norm';
% 
% save 'PS1PS2_240622.mat' PS1 PS2; %save PS1 and PS2 values
% %%
% tic
% x = in_vars_norm';
% t = out_vars_norm';
% 
% % Choose a Training Function
% % For a list of all training functions type: help nntrain
% % 'trainlm' is usually fastest.
% % 'trainbr' takes longer but may be better for challenging problems.
% % 'trainscg' uses less memory. Suitable in low memory situations.
% trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.
% 
% % Create a Fitting Network
% disp_update = 'creating feedfoward network...';
% disp(disp_update);
% 
% hiddenSizes = [14,14,14]; % neurons per layer
% net = feedforwardnet(hiddenSizes,trainFcn);
% 
% % Choose Input and Output Pre/Post-Processing Functions
% % For a list of all processing functions type: help nnprocess
% net.input.processFcns = {'mapminmax'};
% net.output.processFcns = {'mapminmax'};
% 
% %     compet - Competitive transfer function.
% %     elliotsig - Elliot sigmoid transfer function.
% %     hardlim - Positive hard limit transfer function.
% %     hardlims - Symmetric hard limit transfer function.
% %     logsig - Logarithmic sigmoid transfer function.
% %     netinv - Inverse transfer function.
% %     poslin - Positive linear transfer function.
% %     purelin - Linear transfer function.
% %     radbas - Radial basis transfer function.
% %     radbasn - Radial basis normalized transfer function.
% %     satlin - Positive saturating linear transfer function.
% %     satlins - Symmetric saturating linear transfer function.
% %     softmax - Soft max transfer function.
% %     tansig - Symmetric sigmoid transfer function.
% %     tribas - Triangular basis transfer function.
% tf_fnc_name = 'tansig';
% net.layers{1}.transferFcn = tf_fnc_name; %'logsig'
% net.layers{2}.transferFcn = tf_fnc_name; %'logsig'
% net.layers{3}.transferFcn = tf_fnc_name; % 'logsig'
% % Setup Division of Data for Training, Validation, Testing
% % For a list of all data division functions type: help nndivision
% 
% disp_update = 'setting up training | validation | testing data sets...';
% disp(disp_update);
% 
% net.divideFcn = 'dividerand';  % Divide data randomly
% net.divideMode = 'sample';  % Divide up every sample
% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;
% 
% % Choose a Performance Function
% % For a list of all performance functions type: help nnperformance
% net.performFcn = 'mse';%'myperf';  % Mean Squared Error
% net.trainParam.max_fail = 6; % change the number of validation checks 6 is standard
% % Choose Plot Functions
% % For a list of all plot functions type: help nnplot
% net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
%     'plotregression', 'plotfit'};
% 
% % Train the Network
% 
% disp_update = 'beginning to train network...';
% disp(disp_update);
% [net,tr] = train(net,x,t);
% 
% disp_update = 'completed training network...';
% disp(disp_update);
% 
% % Test the Network
% disp_update = 'testing trained network...';
% disp(disp_update);
% 
% [ net, tr, y, e ] = train(net,x,t);
% time_2 = toc;
% 
% P = net(x);
% Error = t - P;
% MSE = mse(t,P);
% NMSE = (mean(Error.^2)) / (mean(var(t',1)));
% R2    = 1-NMSE  ; % coefficient of determination
% R     = sqrt(R2) ;
% 
% y = net(x);
% e = gsubtract(t,y);
% performance = perform(net,t,y)
% 
% % Recalculate Training, Validation and Test Performance
% trainTargets = t .* tr.trainMask{1};
% valTargets = t .* tr.valMask{1};
% testTargets = t .* tr.testMask{1};
% trainPerformance = perform(net,trainTargets,y)
% valPerformance = perform(net,valTargets,y)
% testPerformance = perform(net,testTargets,y)
% 
% % View the Network
% % view(net)
% 
% % Plots
% % Uncomment these lines to enable various plots.
% %figure, plotperform(tr)
% %figure, plottrainstate(tr)
% %figure, ploterrhist(e)
% figure, plotregression(t,y)
% %figure, plotfit(net,x,t)
% 
% % Deployment
% % Change the (false) values to (true) to enable the following code blocks.
% % See the help for each generation function for more information.
% if (true)
%     % Generate MATLAB function for neural network for application
%     % deployment in MATLAB scripts or with MATLAB Compiler and Builder
%     % tools, or simply to examine the calculations your trained neural
%     % network performs.
%     genFunction(net,'myNeuralNetworkFunction');
%     y = myNeuralNetworkFunction(x);
% end
% % if (false)
% disp_update = 'generating new myNeuralNetworkFunction...';
% disp(disp_update);
% 
% genFunction(net,'JOURNALNetworkFunction','MatrixOnly','yes');
% y = JOURNALNetworkFunction(x);
% % end

