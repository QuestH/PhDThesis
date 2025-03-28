function [Pht,ht,hmin,hc1,a,Pc1,Pmax,F_Ph1,Hmin,Etta, Ro]= OneD_EHL(u_current,w_current,Rrout,length,alpha, neta0, nNodes, X, Er )      
%% Input Parameters
 %u_current = u(i); % Surface velocity m/s
 %w_current = w(i); % Load N
%% Dimensionless parameters

%Not the same dimensionless parameters as MMP 

W1=w_current/(length*Er*Rrout); %dimensionless load
U=u_current*neta0/(Rrout*Er); %dimensionless speed of entraining motion
G=alpha*Er;
a=Rrout*sqrt(8*W1/pi); %dimension contact half width
PH=Er*sqrt(W1/(2*pi)); %dimensionless maximum Hertzian pressure?
Lambda=(12*Er*U/PH)*(pi/(8*W1))^1.5; %dimensionless constant
Z=alpha/(5.1*10^-9*(log(neta0)+9.67)); %z for Roelands
Domain= X.*a; %total contact domain
%-------------------------------------------------------------------------%
%% minimum film thickness
%-------------------------------------------------------------------------%
Hmin=3.09*(pi/(8*W1))*G^0.56*U^0.69*W1^(-0.1); 
%-------------------------------------------------------------------------%
%% Time parameters for transient
%-------------------------------------------------------------------------%
%Define Time Parameters
% dt=0.000625;
% DT=4.0393; %dimensonless dt 0.2578*dt/a
%-------------------------------------------------------------------------%
%% Influence Matrix for deformation
%-------------------------------------------------------------------------%
[AK]=defcoeff(nNodes); %influence matrix
%-------------------------------------------------------------------------%
%% Initial Matrices and Boundary Conditions
%-------------------------------------------------------------------------%
DX=X(2)-X(1);
C=log(DX); %for influence loop
A3=0.59/(PH*1e-9); % 
%-------------------------------------------------------------------------%
%%Initial Dimensionless-Herztian Pressure Matrix
%-------------------------------------------------------------------------%
P=zeros(1,nNodes);
for i=2:nNodes-1
    if abs(X(i))<1
        P(i)=sqrt(1-X(i)^2);
    else
        P(i)=0;
    end
end
Pa0=zeros(1,nNodes);
PN=P;
%-------------------------------------------------------------------------%
%% Intial Dimensionless film thickness parabola
%-------------------------------------------------------------------------%
% To compute film thicknes call thick which calls defor function
% (deformation)
[H] = film_thickness(X,nNodes,Hmin,PN,1,AK);% dimensionless  film thickness 
H_old = H;
%-------------------------------------------------------------------------%
%% Initial values for viscosity, density and epsilon
[Etta_old,Ro_old,EPSILON_old]=COEFF(PH,PN,H,nNodes,neta0,Z);
%-------------------------------------------------------------------------%
%% Start of iterative procedure
%-------------------------------------------------------------------------%
% Lambda=(12*Er*U/PH)*(pi/(8*W1))^1.5; % Constant removed from none dimesnionalising?
% Above computed twice? JMW 27/01/21
error=1;%initialise first error value
j0=1;  % ?
Itr=1; % iteration
ItrN=1;% iteration N?

%Begin pressure loop
while error>1e-5 %Pressure iteration criteria
     
    Itr = Itr+1; %Add 1 to iteration count each time
    
    %Function calculating viscosity, density and epsilon using COEFF function
    [Etta,Ro,EPSILON]=COEFF(PH,PN,H,nNodes,neta0,Z); %EPSILON(i)=Ro(i)*H(i)^3/(Etta(i));
    
    for i=2:nNodes-1
        EPSIL1=0.5*(EPSILON(i)+EPSILON(i-1));
        EPSIL2=0.5*(EPSILON(i)+EPSILON(i+1));
        F(i-1)=(EPSIL1*PN(i-1)-PN(i)*(EPSIL1+EPSIL2)+EPSIL2*PN(i+1))/(DX^2)-(Lambda)*(Ro(i)*H(i)-Ro(i-1)*H(i-1))/DX; %F term jacobian
    end
    
    %%  K is coefficient either 1 or zero depending on i related to j - zeroes terms in reynolds
    for i=2:nNodes-1
        EPSIL1=0.5*(EPSILON(i)+EPSILON(i-1));
        EPSIL2=0.5*(EPSILON(i)+EPSILON(i+1));
        for j=2:nNodes-1
            
            if (i-1)~=j
                K(i-1,j)=0;
            else
                K(i-1,j)=1;
            end
            if i~=j
                K(i,j)=0;
            else
                K(i,j)=1;
            end
            if (i+1)~=j
                K(i+1,j)=0;
            else
                K(i+1,j)=1;
            end
            
            IJ=abs(i-j)+1;
            IJ1=abs(i-1-j)+1;
            IJ2=abs(i+1-j)+1;
            
            %Film
            dH_dp(i,j)=-(1/pi)*(AK(IJ)+C)*DX; %influence matrix
            dH_dp(i-1,j)=-(1/pi)*(AK(IJ1)+C)*DX; %influence matrix
            dH_dp(i+1,j)=-(1/pi)*(AK(IJ2)+C)*DX; %influence matrix
            
            %Density
            dRo_dp1(i,j)=(K(i,j)*1.34)/(A3+PN(i));
            dRo_dp1(i-1,j)=(K(i-1,j)*1.34)/(A3+PN(i-1));
            dRo_dp1(i+1,j)=(K(i+1,j)*1.34)/(A3+PN(i+1));
            
            dRo_dp2(i,j)=(K(i,j)*(A3+1.34*PN(i)))/(A3+PN(i))^2;
            dRo_dp2(i-1,j)=(K(i-1,j)*(A3+1.34*PN(i-1)))/(A3+PN(i-1))^2;
            dRo_dp2(i+1,j)=(K(i+1,j)*(A3+1.34*PN(i+1)))/(A3+PN(i+1))^2;
                       
            dRo_dp(i,j)=dRo_dp1(i,j)-dRo_dp2(i,j);
            dRo_dp(i-1,j)=dRo_dp1(i-1,j)-dRo_dp2(i-1,j);
            dRo_dp(i+1,j)=dRo_dp1(i+1,j)-dRo_dp2(i+1,j);
            
            %Viscosity
            dEtta_dp(i,j)=(log(neta0)+9.67)*Z*K(i,j)*(PH/1.96e8)*(1+PN(i)*PH/(1.96e8))^(Z-1)*Etta(i);
            dEtta_dp(i-1,j)=(log(neta0)+9.67)*Z*K(i-1,j)*(PH/1.96e8)*(1+PN(i-1)*PH/(1.96e8))^(Z-1)*Etta(i-1);
            dEtta_dp(i+1,j)=(log(neta0)+9.67)*Z*K(i+1,j)*(PH/1.96e8)*(1+PN(i+1)*PH/(1.96e8))^(Z-1)*Etta(i+1);
            
            %Epsilon
            EPSILON_dp(i,j)=(Ro(i)*3*H(i)^2/(Etta(i)))*dH_dp(i,j)+(H(i)^3/(Etta(i)^2))*(dRo_dp(i,j)*Etta(i)-Ro(i)*dEtta_dp(i,j));
            EPSILON_dp(i-1,j)=(Ro(i-1)*3*H(i-1)^2/(Etta(i-1)))*dH_dp(i-1,j)+(H(i-1)^3/(Etta(i-1)^2))*(dRo_dp(i-1,j)*Etta(i-1)-Ro(i-1)*dEtta_dp(i-1,j));
            EPSILON_dp(i+1,j)=(Ro(i+1)*3*H(i+1)^2/(Etta(i+1)))*dH_dp(i+1,j)+(H(i+1)^3/(Etta(i+1)^2))*(dRo_dp(i+1,j)*Etta(i+1)-Ro(i+1)*dEtta_dp(i+1,j));
            
            
            dEPSIL1_dp(i,j)=0.5*(EPSILON_dp(i,j)+EPSILON_dp(i-1,j));
            dEPSIL2_dp(i,j)=0.5*(EPSILON_dp(i,j)+EPSILON_dp(i+1,j));
            
            
            dFdp(i-1,j-1)=(1/(DX^2))*(((dEPSIL1_dp(i,j)*PN(i-1))+(EPSIL1*K(i-1,j)))-((dEPSIL1_dp(i,j)+dEPSIL2_dp(i,j))*PN(i)+(EPSIL1+EPSIL2)*K(i,j))+...
                (dEPSIL2_dp(i,j)*PN(i+1)+EPSIL2*K(i+1,j)))-(Lambda/DX)*((dRo_dp(i,j)*H(i)+Ro(i)*dH_dp(i,j))-(dRo_dp(i-1,j)*H(i-1)+Ro(i-1)*dH_dp(i-1,j)));
        end
        
        EPSILON_dH00(i)=(Ro(i)*3*H(i)^2/(Etta(i)));
        EPSILON_dH00(i-1)=(Ro(i-1)*3*H(i-1)^2/(Etta(i-1)));
        EPSILON_dH00(i+1)=(Ro(i+1)*3*H(i+1)^2/(Etta(i+1)));


                
        dEPSIL1_dH00(i)=0.5*(EPSILON_dH00(i)+EPSILON_dH00(i-1));
        dEPSIL2_dH00(i)=0.5*(EPSILON_dH00(i)+EPSILON_dH00(i+1));
        
        dFdH00(i-1)=(1/DX)^2*(dEPSIL1_dH00(i)*PN(i-1)-( dEPSIL1_dH00(i)+ dEPSIL2_dH00(i))*PN(i)+dEPSIL2_dH00(i)*PN(i+1))-(Lambda/DX)*(Ro(i)-Ro(i-1));
    end
    
    
    Loadinteg=0;
    for i=1:nNodes-1
        Loadinteg=Loadinteg+DX*0.5*(PN(i)+PN(i+1)); %Integrating pressure across domain to find hydrodynamic load
    end
    
    %% Jacobians
    FLoad=-pi/2+Loadinteg;
    FLoad_t=FLoad;
    F1=[F';FLoad];
    dload_dp=[DX*ones(1,nNodes-2) 0];
    Jacob=[dFdp dFdH00';dload_dp];
    Pold=PN;
    
    for i=2:nNodes-1
        Pold1(i-1)=PN(i);
    end
    
    PP=[Pold1';Hmin];
    cc=0.1;
    dZ=Jacob\F1;
    XX=PP-cc*(dZ);
    PN(1)=0;
    
    for i=2:nNodes-1
        PN(i)=XX(i-1);
    end
    
    PN(nNodes)=0;
    Hmin=XX(nNodes-1);
    
    for i=1:nNodes
        if PN(i)<0
            PN(i)=0;
            
        end
    end
    
    %Calculating error
    s=0;
    for i=1:nNodes
        s=s+abs(PN(i)-Pold(i));
    end
    
    error=abs(s/sum(PN));
    [j0 error];
    %       [j0 Pa0]
    j0=j0+1;
    [H]=film_thickness(X,nNodes,Hmin,PN,0,AK);
    %         H=H*10;
    
    ItrN(Itr)= Itr; %Save iteration numbers in matrix 
    itr_error(Itr) = error; %value of error at each iteration
end

%% Results

[Etta,Ro,EPSILON]=COEFF(PH, PN,H,nNodes,neta0,Z); %Pressure and viscosity 
Pht=PN*PH; %Hydrodynamic pressure

%Load
LoadPh=0;

for i=1:nNodes-1
    LoadPh=LoadPh+DX*0.5*(PN(i)+PN(i+1)); % Differentiate pressure to find load
end

F_Ph=LoadPh; %Dimensionless force
F_Ph1=length*a*PH*F_Ph; %Dimensionalised force

%Pressure
Pt=PN;

%Film Thickness
Ht=H;
ht= Ht*(a^2)/Rrout;


H00t=Hmin; %Dimensionless minimum film thickness as calculated using extrapolated formulae
h00= H00t*(a^2)/Rrout; %Dimensionalised minimum film thickness as calculated using extrapolated formulae

MinHt=min(H); %Dimensionless minimum film thickness
hmin= MinHt*(a^2)/Rrout; %Dimensionalised minimum film thickness

Hc1=H(150); %Dimensionless central film thickness
hc1=Hc1*(a^2)/Rrout; %Dimensionalised central film thickness

Pc1 = Pht(150); %Central pressure
Pmax = max(Pht); %Max pressure

H_old=H;
Ro_old=Ro;
Itr=Itr+1;

end


% %% Plotting
% 
% figure(1)
% 
% left_color = [0 0 0];
% right_color = [0 0 0];
% set(figure(1),'defaultAxesColorOrder',[left_color; right_color]);
% ax = gca;
% ax.FontSize = 20;
% 
% 
% yyaxis left
% plot(Domain,Pt,'k')
% ylabel('Pressure (GPa)')
% 
% yyaxis right
% plot(Domain,ht/1e-6,'--k')
% ylabel('Film Thickness (\mum)')
% 
% xlabel('Contact Domain (mm)')
% %     xlim([4 4.17])
% title('Pressure and Film Thickness')
% 
