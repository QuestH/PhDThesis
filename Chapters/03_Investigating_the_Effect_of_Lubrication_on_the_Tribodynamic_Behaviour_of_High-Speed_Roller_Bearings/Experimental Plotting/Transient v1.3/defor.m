%% Deformation Calculation
function [ V ] = defor( N,DX,PN)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global AK 

PAI1=0.318309886;
C=log(DX);
V=zeros(1,N);
for  I=1:N
    V(I)=0.0;
    for  J=1:N
        IJ=abs(I-J)+1;      
        V(I)=V(I)+(AK(IJ)+C)*DX*(PN(J));
    end
end
for I=1:N
    V(I)=-PAI1*V(I);
end





