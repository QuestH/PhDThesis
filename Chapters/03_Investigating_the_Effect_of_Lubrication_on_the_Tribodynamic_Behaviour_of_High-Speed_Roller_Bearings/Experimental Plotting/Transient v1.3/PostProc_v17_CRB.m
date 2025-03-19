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
    kn=5e8; % inner race, roller, outer race stiffness
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
        
        if delta_in(i)<0 && delta_in(i)>-0.1e-6 %if 0-2 micron separation, still assume 2 micron separation
                delta_in(i)=-0.1e-6;
        end
        
        if delta_out(i)<0 && delta_out(i)>-0.1e-6 %if 0-2 micron separation, still assume 2 micron separation
                delta_out(i)=-0.1e-6;
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
        
        filmNEW_in(i) = Rrin*1.93*(Ue_in^0.69)*(Ge^0.56)*(We_in^-0.1);
        filmNEW_out(i) = Rrout*1.93*(Ue_out^0.69)*(Ge^0.56)*(We_out^-0.1);
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
        filmNEW_in(i)= real(Rrin*1.93*(Ue^0.69)*(Ge^0.56)*(We_in^-0.1));
        filmNEW_out(i)= real(Rrout*1.93*(Ue^0.69)*(Ge^0.56)*(We_out^-0.1));
        filmNEW(i) = filmNEW_in(i)+filmNEW_out(i);
        
        end
        
        Wt_in(i)=real(W_in(i));
        Ht_in(i)=filmNEW_in(i);
        Wt_out(i)=real(W_out(i));
        Ht_out(i)=filmNEW_out(i);
        
        end
              
%% Friction Calculations

tau0 = 2e6;                  %Eyring Stress set to 3MPa
zeta1 = 0.3;                 %Matrix of Shear Strength of Asperities
kappa = 1/(0.18e6);              %Matrix of Average Asperity Tip Radius of Curvature (1/m)
sigma = 0.197e-6;                %Matrix of Composite Surface Roughness (um)
zeta = 0.00116e12;            %Matrix of Aperity Density (1/m^2)
Lub_cond= 1600;                %lubricant conductivity (W/mK)
Sol_cond= 46;                %solid conductivity (W/mK)
Sol_dens= 7850;               %solid density (kq/m3)
Sol_SHeat= 470;      %solid specific heat (J/kgK)

% b_in(i) = ((4*Wt_out(i)*(((1-poas1^2)/E1)+((1-poas2^2)/E2)))/(pi*length*((1/Rrout)+(1/inf))))^(1/2); %half width of contact area
% b_out(i) = ((4*Wt_in(i)*(((1-poas1^2)/E1)+((1-poas2^2)/E2)))/(pi*length*((1/Rrin)+(1/inf))))^(1/2); %half width of contact area

%half-width of contact inner and outer
b_in(i) = sqrt(8*W_in(i)*Rrin/(pi*Er));
b_out(i) = sqrt(8*W_out(i)*Rrout/(pi*Er));
        
A_out(i) = 2*b_out(i)*length;                %Apparent Contact Area (m^2)
A_in(i) = 2*b_in(i)*length;

P_bar_out=Wt_out(i)/A_out(i); %average pressure at the contact point
P_bar_in=Wt_in(i)/A_in(i);

%% Viscous Friction

ixi_out(i)= (4/pi)*(Lub_cond/(Ht_out(i)/Rrout))*((P_bar_out/(Er*Rrout*Sol_cond*Sol_dens*Sol_SHeat*U(i)))^(1/2));
ixi_in(i)= (4/pi)*(Lub_cond/(Ht_in(i)/Rrin))*((P_bar_in/(Er*Rrin*Sol_cond*Sol_dens*Sol_SHeat*U(i)))^(1/2));

Mu_out(i)= 0.87*alpha*tau0+1.74*(tau0/P_bar_out)*log((1.2/(tau0*Ht_out(i)))*((2*Lub_cond*neta0)/(1+9.6*ixi_out(i)))^(1/2));
Mu_in(i)= 0.87*alpha*tau0+1.74*(tau0/P_bar_in)*log((1.2/(tau0*Ht_in(i)))*((2*Lub_cond*neta0)/(1+9.6*ixi_in(i)))^(1/2));

fv_out(i)= real(Mu_out(i)*Wt_out(i));
fv_in(i)= real(Mu_in(i)*Wt_in(i));

%% Boundary Friction

%landa outer race
landa_out(i) = Ht_out(i)/sigma;

if landa_out(i) <= 4                                                          %Integrals of the Gaussian Distribution

F2_out(i) = 3.5e-5.*(4.00001-landa_out(i)).^6.9; %excite

elseif landa_out(i) > 4
    F2_out(i) = 0;
else F2_out(i) = NaN;
end

if landa_out(i) <= 4

F5_2_out(i) = 4.4086e-5.*(4.00001-landa_out(i)).^6.804; %Excite

    
elseif landa_out(i) > 4
    F5_2_out(i) = 0;
else F5_2_out(i) = NaN;
end

%landa inner race
landa_in(i) = Ht_in(i)/sigma ;

if landa_in(i) <= 4                                                            %Integrals of the Gaussian Distribution
    
F2_in(i) = 3.5e-5.*(4.00001-landa_in(i)).^6.9; %excite

elseif landa_in(i) > 4
    F2_in(i) = 0;
else F2_in(i) = NaN;
end

if landa_in(i) <= 4

F5_2_in(i) = 4.4086e-5.*(4.00001-landa_in(i)).^6.804; %Excite

elseif landa_in(i) > 4
    F5_2_in(i) = 0;
else F5_2_in(i) = NaN;
end

% Aa_out(i) = (pi^2)*((zeta*kappa*sigma)^2)*(sqrt(sigma/kappa))*A_out(i)*F2_out(i); %Area occupied by the asperities
% Aa_in(i) = (pi^2)*((zeta*kappa*sigma)^2)*(sqrt(sigma/kappa))*A_in(i)*F2_in(i); %Area occupied by the asperities

Aa_out(i) = (pi^2)*((0.055)^2)*A_out(i)*F2_out(i); % Mohammadpour pg.247 thesis
Aa_in(i) = (pi^2)*((0.055)^2)*A_in(i)*F2_in(i); % Mohammadpour pg.247 thesis

% Wa_out(i) = ((16*sqrt(2))/15)*pi*((zeta*kappa*sigma)^2)*(sqrt(sigma/kappa))*Er*A_out(i)*F5_2_out(i); %Asperity load
% Wa_in(i) = ((16*sqrt(2))/15)*pi*((zeta*kappa*sigma)^2)*(sqrt(sigma/kappa))*Er*A_in(i)*F5_2_in(i); %Asperity load

Wa_out(i) = ((8*sqrt(2))/15)*pi*((0.055)^2)*(sqrt(1e-3))*Er*A_out(i)*F5_2_out(i); %Asperity load 0.055 and 10^-3 from mohammadpour paper https://doi.org/10.1177/1350650114537805
Wa_in(i) = ((8*sqrt(2))/15)*pi*((0.055)^2)*(sqrt(1e-3))*Er*A_in(i)*F5_2_in(i); %Asperity load


fb_out(i) = (tau0*Aa_out(i))+(zeta1*Wa_out(i)); 
fb_in(i) = (tau0*Aa_in(i))+(zeta1*Wa_in(i)); 


%% Total Friction
friction_out(i)= real(fv_out(i))+real(fb_out(i));
friction_in(i)= real(fv_in(i))+real(fb_in(i));


end
 coeff_f = friction_out./Wt_out;
%     
% %% Film stiffness
% for i = lower+1:upper-1
%     dh_dw(i) = (Ht_out(i+1)-Ht_out(i-1))/(Wt_out(i+1)-Wt_out(i-1));
% end
% 
% stiff = -1*(1./dh_dw);

%% Selecting nodes to pass to EHL

% selected_nodes = [231869,231928,231983,232040,232094,232151,232210,232267,232326,232380];
% Wt_out_EHL = Wt_out([1],selected_nodes);
% u_EHL = U([1],selected_nodes);
%                 
% save outputTribological.mat Wt_out_EHL length u_EHL Rrout

save Wt_out_lub.mat Wt_out speed

%% Plotting
    
%Load vs Film Thickness
    figure(1)
    left_color = [0 0 0];
    right_color = [0 0 0];
    set(figure(1),'defaultAxesColorOrder',[left_color; right_color]);
    ax = gca;
    ax.FontSize = 20;
    
    yyaxis left
    plot(time(1:upper)*3600,Wt_out,'k')
    ylabel('Load (N)','color','k')

    yyaxis right
    plot(time(1:upper)*3600,Ht_out/1e-6,'--k')
    ylabel('Film Thickness (\mum)','color', 'k')

    xlabel('Rotational Speed (rpm)')
    xlim([0 15000])
    title('Outer Ring - Load Carried by Each Roller and Corresponding Film Thickness')
  
 %Boundary Friction vs Lambda   
    figure(2)
    left_color = [0 0 0];
    right_color = [0 0 0];
    set(figure(2),'defaultAxesColorOrder',[left_color; right_color]);
    ax = gca;
    ax.FontSize = 20;
    
    yyaxis left
    plot(time(1:upper)*3600,(real(fb_out)),'k')
    ylabel('Boundary Friction/ N','color','k')

    yyaxis right
    plot(time(1:upper)*3600,landa_out,'--k')
    ylabel('Lambda Ratio','color', 'k')

    xlabel('Rotational Speed/ rpm')
    xlim([0 15000])
    title('Boundary Friction vs Lambda')

    figure(2)
    plot(speed,real(fb_out),'k','LineWidth',0.75)
    set(gca,'FontSize',18)
    title ('Boundary Friction')
    xlabel ('Speed/ rpm')
    ylabel ('Boundary Friction/ N')
    xlim ([0 15000])
    
    figure(3)
    plot(speed,real(fv_out),'k','LineWidth',0.75)
    set(gca,'FontSize',18)
    title ('Viscous Friction')
    xlabel ('Speed/ rpm')
    ylabel ('Viscous Friction/ N')
    xlim ([0 15000])
    
    figure(4)
    plot(speed,real(coeff_f),'k','LineWidth',0.75)
    set(gca,'FontSize',18)
    title ('Coefficient of Friction')
    xlabel ('Speed/ rpm')
    ylabel ('Friction Coefficient')
    xlim ([0 15000])
    

toc