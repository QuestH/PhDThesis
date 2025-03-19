%% Importing and Integrating Data
clear
clc
load('inputDATA_CRB_Paper5.mat')

%Select data length
vsl = vsl;
ysli = ysli;

%Shaft Accelerometer
t = 1e-5:1e-5:4.17; %Time from ODE45

        
%% Fast-Fourier Transform 

%FFT original signal

        Fs = 100000;      % Sampling frequency                    
        T = 1/Fs;         % Sampling period       
        L = length(t);    % Length of signal
        t = (0:L-1)*T;

        %Compute the fourier transform of the signal
        YT = fft(ysli);
        VT = fft(vsl);
                
        % AT - Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
        YTP2 = abs(YT/L);
        YTP1 = YTP2(1:L/2+1);
        YTP1(2:end-1) = 2*YTP1(2:end-1);
        
        % YT - Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
        VTP2 = abs(VT/L);
        VTP1 = VTP2(1:L/2+1);
        VTP1(2:end-1) = 2*VTP1(2:end-1);
        
        
        %Define the frequency domain f and plot the single-sided amplitude spectrum P1       
        f = Fs*(0:(L/2))/L;    %fourier original
        



   %% Plotting
   
    figure(1)
    plot(t,ysli,'k')
    set(gca,'FontSize',18)
    title ('Shaft Displacement X Direction')
    xlabel ('Time (s)')
    ylabel ('Displacement (m)')

    figure(2)
    plot(f,VTP1,'k') 
    set(gca,'FontSize',18)
    title('FFT of Velocity Signal - X Direction')
    xlabel('Frequency (Hz)')
    ylabel('Veloctity FFT Spectra (m/s)')
    
    figure(3)
    plot(f,YTP1,'k') 
    set(gca,'FontSize',18)
    title('FFT of Displacement Signal - X Direction')
    xlabel('Frequency (Hz)')
    ylabel('Acceleration FFT Spectra (m)')
 
    
    figure(4)
    left_color = [0 0 0];
    right_color = [0 0 1];
    set(figure(4),'defaultAxesColorOrder',[left_color; right_color]);
    ax = gca;
    ax.FontSize = 18;

    yyaxis left
    plot(f,YTP1,'k')
    ylabel('Displacement FFT Spectra (m)','color','k')
    
    yyaxis right
    plot(f,VTP1,'b')
    ylabel('Velocity FFT Spectra (m)','color','b')
    
    xlabel('Frequency (Hz)')
    xlim ([0 10000])
    title('Displacement vs Velocity Spectra') 
    
    
    
    %limit cycles
    %plot(AT(2:360),YT(2:360))