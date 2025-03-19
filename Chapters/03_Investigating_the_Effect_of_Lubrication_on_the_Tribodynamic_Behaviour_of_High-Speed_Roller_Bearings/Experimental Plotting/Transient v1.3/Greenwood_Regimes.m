clc
clear all

%% Inputs

load Greenwood_EHL.mat

Vn = neta0*U/(Er*Rrout);
Gn = Er*alpha;
Wn = Wt_out/(Er*(Rrout^2));


%% Subgroups

% g1 = (2*alpha^2*Wt_out.^3./(neta0*U.*Rrout^2)).^0.5;
% g3 = (2*pi*Wt_out.^2./(neta0*U.*Er*Rrout)).^0.5;

for i =1:length(Wt_out)
    
   g1(i) = (alpha^2*Wt_out(i)^3/(neta0*U(i)*Rrout))^0.5;
   g3(i) = (Wt_out(i)^2/(neta0*U(i)*Er*Rrout))^0.5;
   
end

ge = Wn./(Vn.^0.5); %line contact 
gv = Wn.^(3/2)*Gn./(Vn.^0.5); %line contact 

% figure(1)
% loglog(real(ge),real(gv))
% xlim([10^-1 10^3])
% ylim([10^-1 10^4])

figure(2)
loglog(real(g3),real(g1))
xlim([10^-1 10^3])
ylim([10^-1 10^4])
