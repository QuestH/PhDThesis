function I=Gauss_Legendre(f,a,b,N,varargin)
%Gauss_Legendre integration off over [a,b] with N gridpoints
%Never try N larger than 25
[t,w]=Gausslp(N);
x=((b-a)*t+a+b)/2;%Eq.(5.9.9)
fx=feval(f,x,varargin{:});
I=w*fx'*(b-a)/2;%Eq.(5.9.10)

function[t,w]=Gausslp(N)
if N<0, fprintf('\nGauss-Legendrepolynomialofnegativeorder??\n');
else
    t=roots(Lgndrp(N))';% make it arow vector
    A(1,:)=ones(1,N);b(1)=2;
    for n=2:N %Eq.(5.9.7)
        A(n,:)=A(n-1,:).*t;
        if mod(n,2)==0, b(n)=0;
        else b(n)=2/n; %Eq.(5.9.8)
        end
    end
    w=b/A';
end

function p=Lgndrp(N)%Legendrepolynomial
if N<=0
    p=1; %n*Ln(t)=(2n-1)tLn-1(t)-(n-1)Ln-2(t)Eq.(5.9.6b)
elseif N==1
    p=[1 0];
else
    p=((2*N-1)*[Lgndrp(N-1) 0]-(N-1)*[0 0 Lgndrp(N-2)])/N;
end

