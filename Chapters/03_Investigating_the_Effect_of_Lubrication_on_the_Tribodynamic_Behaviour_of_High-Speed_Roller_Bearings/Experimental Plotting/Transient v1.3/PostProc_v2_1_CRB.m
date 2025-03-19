clear all
clc

tic 

load 'inputDATA_CRB_Paper5.mat'  % Import the .mat file which contains data

steps=417000; % This is number of time steps  
lower = 1;
upper = 417000;
EHLrelax=1.3;% relaxation factor of EHL film calculations original 0.5

%% Properties of bearing

%Bearing geometry
    Nelem=13; % number of element
    c= 0; %clearance + for clearance and - for fitting (-1e-6) - m
    length=9e-3; % bearing length - m
    Rball=0.00375; %roller radius - m
    Rinner=0.01575; %inner race radius - m
    Router = 0.0232; %outer race radius - m
    Rrin = 1/(1/Rball+1/Rinner); %reduced radius inner race
    Rrout = 1/(1/Rball-1/Router); %reduced radius outer race
    Rc = 19.45e-3; %radius from bearing centre to centre of roller
    density_r = 7850; %density of steel used in roller - kg/m^3
    mass_r = pi*Rball^2*length*density_r; %kg

    %Lubricant properties
    ep=inf;
    neta0=0.08;
    alpha=21.1e-9;

    %Solid properties
    poas1=0.3;
    poas2=0.3;
    E1=211000000000;
    E2=211000000000;
    Er=1/(((1-(poas1)^2)/E1)+((1-(poas2)^2)/E2)); %note pi in reduced elastic modulus
    
    %Stiffness properties of bearing
    gama=10/9; % power of load deflection
    kn=3.366e9; % inner race, roller, outer race stiffness
    g=9.81; % gravity
    
    
for i=lower:1:upper
        
        %% Angular velocity and position
        
        omegaCAGE(i)=0.4*(2*pi*15012/60*i/steps); % fundamental train frequency
        if i==lower
            ang(i)=omegaCAGE(i)*time(i,1);
        else
            ang(i)=ang(i-1)+omegaCAGE(i)*(time(i)-time(i-1)); % instantaneouse position of each element
        end
                
        %Inner race velocity
        omegaINNER(i)=10+2*pi*15012/60*i/steps;
        
        %Surface velocity
        U(i)=abs(omegaCAGE(i)-omegaINNER(i))*Rinner/2; %speed of entraining motion definition
        
        %Centrifugal force
        cf(i) = 0; %mass_r*(omegaCAGE(i)^2)*Rc; %centrifugal force on roller 
        
          
%% Instantaneous deflection of each element
        
        film_in(i)=0; % instantaneouse film thickness at inner race contact
        film_out(i)=0; % instantaneous film thickness at outer race contact
        film(i) = 0;
        
        delta(i)=(film_in(i)+film_out(i)-c)+(yt(i)-ysli(i))*cos(ang(i))+(yt(i)-ysli(i))*sin(ang(i));%+y(i,3)*sin(ang(i)); % instantaneouse deflection of each element
        delta_in(i) = 0.5*delta(i);
        delta_out(i) = delta_in(i);
               
        if delta(i)==0   %0.1e-6 delta means film is near negligible and very low load
   
        Wt_in(i)=0; 
        Ht_in(i)=0;
        Wt_out(i)=0; 
        Ht_out(i)=0;
            
        elseif (delta(i)<0)
            
        delta_in(i) = 0.5*delta(i);
        delta_out(i) = delta_in(i);
        
        if delta_in(i)<0 && delta_in(i)>-0.025e-6 %if 0-2 micron separation, still assume 2 micron separation
                delta_in(i)=-0.025e-6;
        end
        
        if delta_out(i)<0 && delta_out(i)>-0.025e-6 %if 0-2 micron separation, still assume 2 micron separation
                delta_out(i)=-0.025e-6;
        end
            
        Ht_in(i) = abs(delta_in(i)); % film still plots correctly, but load not calculated hydrodynamically if between 0 and 4 micron as this is not valid
        Ht_out(i) = abs(delta_out(i));
                        
        W_in(i) = 2*(length/2)*U(i)*neta0*Rrin/abs(delta_in(i))-cf(i);
        W_out(i) = 2*(length/2)*U(i)*neta0*Rrout/abs(delta_out(i))+cf(i);
        
        if W_in(i)<0 %cannot have negative
            W_in(i)=0;
        end
        
        if W_out(i)<0 %cannot have negative
            W_out(i)=0;
        end
        
        Wt_in(i)= real(W_in(i));
        Wt_out(i) = real(W_out(i));

 
              else
         W_in(i) = kn*delta_in(i)^gama; % initial guess! instantaneouse reaction load of each element
        W_out(i) = kn*delta_out(i)^gama; % initial guess! instantaneouse reaction load of each element
        
 
        %Non-dimensional parameters
        Ue_in = (neta0*U(i))/(Er*Rrin);
        Ue_out = (neta0*U(i))/(Er*Rrout);
        Ge = Er*alpha;
        We_in = W_in(i)/(Er*Rrin*length);
        We_out = W_out(i)/(Er*Rrout*length);
        
        filmNEW_in(i) = Rrin*3.06*(Ue_in^0.69)*(Ge^0.56)*(We_in^-0.1);
        filmNEW_out(i) = Rrout*3.06*(Ue_out^0.69)*(Ge^0.56)*(We_out^-0.1);
        filmNEW(i) = filmNEW_in(i)+filmNEW_out(i);
        
        nit=1; %number of iterations
       
        while abs((filmNEW(i)-film(i))/film(i))>=0.01 % Convergence criteria - until error is 0.05
            
        nit=nit+1;
        NIT(nit)=nit; %NIT is vector
        errit(nit)= abs((filmNEW_in(i)-film_in(i))/film_in(i)); %Error
                
        %film
        film_in(i)=film_in(i)+EHLrelax*(filmNEW_in(i)-film_in(i));
        film_out(i)=film_out(i)+EHLrelax*(filmNEW_out(i)-film_out(i));
        film(i) = film_in(i)+film_out(i);
        
        %stiffness inner and outer
        k_in(i) = 3.366e9; %from curve fitting
        k_out(i) = 3.078e9; %from curve fitting
        k_total(i) = 1/(1/k_in(i)+1/k_out(i));
        
        %total delta        
        delta(i)=(film_in(i)+film_out(i)-c)+(yt(i)-ysli(i))*cos(ang(i))+(yt(i)-ysli(i))*sin(ang(i));%+y(i,3)*sin(ang(i)); % instantaneouse deflection of each element
        
        %delta inner and outer based on stiffness
        delta_out(i) = delta(i)/((k_out(i)/k_in(i))^(1/gama)+1);
        delta_in(i) = delta(i)-delta_out(i);
        
        %Load inner and outer
        W_in(i)=k_in(i)*delta_in(i)^gama;
        W_out(i)=k_out(i)*delta_out(i)^gama;
        
        Ue = (neta0*U(i))/(Er*Rball);
        Ge = Er*alpha;
        We_in = W_in(i)/(Er*Rball*length);
        We_out = W_out(i)/(Er*Rball*length);
        
        %new film calculation
        filmNEW_in(i)= real(Rrin*3.06*(Ue^0.69)*(Ge^0.56)*(We_in^-0.1));
        filmNEW_out(i)= real(Rrout*3.06*(Ue^0.69)*(Ge^0.56)*(We_out^-0.1));
        filmNEW(i) = filmNEW_in(i)+filmNEW_out(i);
        
        end
        
        Wt_in(i)=real(W_in(i));
        Ht_in(i)=filmNEW_in(i);
        Wt_out(i)=real(W_out(i));
        Ht_out(i)=filmNEW_out(i);
        
        end
              

end

%     
% %% Film stiffness
% for i = lower+1:upper-1
%     dh_dw(i) = (Ht_out(i+1)-Ht_out(i-1))/(Wt_out(i+1)-Wt_out(i-1));
% end
% 
% stiff = -1*(1./dh_dw);

%% Selecting nodes to pass to Waterfall

save Waterfall_load.mat Wt_out time yt delta_out ysli delta
save Greenwood_All.mat Wt_out alpha neta0 Rrout Er U


%% Plotting
    
% %Load vs Film Thickness
%     figure(1)
%     left_color = [0 0 0];
%     right_color = [0 0 0];
%     set(figure(1),'defaultAxesColorOrder',[left_color; right_color]);
%     ax = gca;
%     ax.FontSize = 20;
%     
%     yyaxis left
%     plot(time(1:upper)*3600,Wt_out,'k')
%     ylabel('Load (N)','color','k')
% 
%     yyaxis right
%     plot(time(1:upper)*3600,Ht_out/1e-6,'--k')
%     ylabel('Film Thickness (\mum)','color', 'k')
% 
%     xlabel('Rotational Speed (rpm)')
%     xlim([0 15000])
%     title('Outer Ring - Load Carried by Each Roller and Corresponding Film Thickness')
%   
%   
%     figure(6)
%     plot(speed,Wt_out,'k','LineWidth',0.75)
%     set(gca,'FontSize',18)
%     title ('Speed and Load Nodes for EHL Model')
%     xlabel ('Rotational Speed (rpm)')
%     ylabel ('Load (N)')
%     xlim ([8310 8420])
%     
%     figure(10)
%     plot(time(1:upper)*3600,Wt_out,'k','LineWidth',0.75)
%     set(gca,'FontSize',18)
%     xlabel('Rotational Speed (rpm)')
%     xlim([8340 8400])
%     ylabel('Load (N)','color','k')
%     ylim([0 2200])
    
    figure(11)
    left_color = [0 0 0];
    right_color = [0 0 0];
    set(figure(11),'defaultAxesColorOrder',[left_color; right_color]);
    ax = gca;
    ax.FontSize = 24;
    
    yyaxis left
    plot(time(1:upper)*3600,Wt_out,'k','LineWidth',1)
    ylabel('Load (N)','color','k')
    ylim([0 2200])

    yyaxis right
    plot(time(1:upper)*3600,Ht_out/1e-6,'--k','LineWidth',1)
    ylabel('Film Thickness (\mum)','color', 'k')
    ylim([0 4])

    xlabel('Rotational Speed (rpm)')
    xlim([8340 8400])
    title('Outer Ring - Load Carried by Each Roller and Corresponding Film Thickness')

toc