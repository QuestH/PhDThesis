%% Pressure, Viscosity and Coefficient Calculation
function [Etta,Ro,EPSILON]=COEFF(P,H,N)

global Lambda PH neta0
A3=0.59/(PH*1e-9);
for i=1:N
    Ro(i)=(A3+1.34*P(i))/(A3+P(i));
    Etta(i)=exp((log(neta0)+9.67)*(-1+(1+P(i)*PH/(1.96e8))^0.68));
    EPSILON(i)=Ro(i)*H(i)^3/(Etta(i));

end    




