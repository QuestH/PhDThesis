%% Influence Matrix
function [AK]=defcoeff(N)
for  I=1:N
    for  J=1:N
        IJ=abs(I-J)+1;
        for  I=1:IJ
            AK(I)=(I-1+0.5)*(log(abs(I-1+0.5))-1.)-(I-1-0.5)*(log(abs(I-1-0.5))-1.);
        end
       
    end
end