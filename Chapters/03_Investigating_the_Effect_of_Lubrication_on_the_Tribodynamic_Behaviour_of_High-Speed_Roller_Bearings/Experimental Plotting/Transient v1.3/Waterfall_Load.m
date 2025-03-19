clear
clc

load Waterfall_Load.mat

Wt = Wt_out'*10^-9;
delta = delta';
delta_o = delta_out';
delta_original = (yt-ysli);

num=1;
for tt=2:0.1:14


rot_freq(num)=tt*1000/60;
one_period_steps(num)=floor((1/rot_freq(num))/1e-5);

begin(num)=tt*floor(417000/15)-6*one_period_steps(num);
endd(num)=tt*floor(417000/15)+6*one_period_steps(num);
Signal=Wt(begin(num):1:endd(num),1)/1e-3;
ttt=time(begin(num):1:endd(num),1);

timestep=1e-5;

L=length(Signal);

Fs=1/timestep;
T=1/Fs;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(Signal,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1); 

value=(2*abs(Y(1:NFFT/2+1)));

rot_speed=tt*1000+zeros(1,size(f,2));

plot3(rot_speed,f,value,'k');
    set(gca,'FontSize',18)
    title('FFT of Displacement Signal - X Direction')
    xlabel('Rotational Speed (RPM)')
    ylabel('Frequency (Hz)')
    zlabel('Displacement FFT Spectra (m)')
    
ylim([0 3000])
view([290 40])
hold on

% fullvalue=[fullvalue;value'];

% toExcite=[ttt Signal];
% 
% save (str1)

num=num+1;

end

