
clc
clear all

length=9e-3; % bearing length - m
Rball=0.00375; %roller radius - m
Rinner=0.01575; %inner race radius - m
Router = 0.0232; %outer race radius - m
Rrin = 1/(1/Rball+1/Rinner); %reduced radius inner race
Rrout = 1/(1/Rball-1/Router); %reduced radius outer race
poas1=0.3;
poas2=0.3;
E1=211000000000;
E2=211000000000;
Er=pi/(((1-(poas1)^2)/E1)+((1-(poas2)^2)/E2)); %note pi in reduced elastic modulus

for i = 1:20000
    
 F(i) = i;
 delta_out(i) = (F(i)/(pi*Er*length))*(log(4*pi*Er*Rrout*length/(F(i)))+1);
 delta_in(i) = (F(i)/(pi*Er*length))*(log(4*pi*Er*Rrin*length/(F(i)))+1);

end

