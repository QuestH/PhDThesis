%% Film Thickness Calculation
function [H]=thick(X,N,Hmin,P,K1)

DX=X(2)-X(1);
if K1~=1
    Loadinteg=0;
for i=1:N-1
    Loadinteg=Loadinteg+DX*0.5*(P(i)+P(i+1));
end
DLoad=pi/2-Loadinteg;

end

H= Hmin*ones(N,1);
[ deflect ] = defor( N,DX,P);

for i=1:N
    H(i)=H(i)+deflect(i)+0.5*X(i)^2;
end
end



