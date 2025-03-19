%% Importing and Integrating Data
clear all
clc
load('Paper5-4.17V-4.17S-12345-NoShaker.mat')

%Select data length
dselect = 1:417000;

%Shaft Accelerometer
t = time(dselect,1);
    %extract shaft acceleration data
    as = data(dselect,2)/1e-2; %acceleration at load ball holder (shaft)

    %detrend shaft results to remove accelerometer drift.
    detrend_as = detrend(as); % Apply detrend, which performs a linear fit to as data and then removes the trend from it.

    %Integration
    
    vs = cumtrapz(t,detrend_as);    % Velocity
    cs = cumtrapz(t,vs);    % Displacement
    
    %Normalising Displacement about Zero
    mean(cs);
    cs_norm = cs - mean(cs);

%Test Bracket 12 o'clock Accelerometer
    %extract test bracket acceleration data
    at = data(dselect,3)/1e-2; %acceleration at test bracket

    %detrend test bracket results to remove accelerometer drift.
    detrend_at = detrend(at); % Apply detrend, which performs a linear fit to as data and then removes the trend from it.
   
    %Integration
    t = time(dselect,1);
    vt = cumtrapz(t, detrend_at);    % Velocity
    detrend_vt = detrend(vt);
    ct = cumtrapz(t, vt);    % Displacement
    detrend_ct = detrend(ct);
    
    %Normalising Displacement about Zero
    mean(ct);
    ct_norm = ct - mean(ct);

%Shaft Laser Vibrometer - Velocity
    %extract laser vibrometer velocity data
    vsl = data(:,5)*200e-4; % laser vibrometer acceleration measured at bearing - Vertical velocity 200mm/s/V

    %detrend results
    detrend_vs = detrend(vsl);

    %Integration
    t = time(:,1);
    csli = cumtrapz(t,detrend(vsl)); %displacement laser shaft integrated
    
    %Normalising Displacement about Zero
    mean(csli);
    csli_norm = csli - mean(csli);

 %% Butterworth Filter  

    %Buttord - find order  and  cutoff frequency
        Wp = [70 10000]/50000; %passband between 20 and 10000hz with 100kHz sampling rate
        Ws = [35 12500]/50000; %stopband set to 0.5 and 1.25 times either side of passband
        Rp = 3; % 3 dB of passband ripple
        Rs = 7; % 7 dB attenuation in the stopbands
        
        %Wp = [50 5000]/50000; %passband between 50 and 5000hz with 100kHz sampling rate
        %Ws = [25 6250]/50000; %stopband set to 0.5 and 1.25 times either side of passband
        %Rp = 3; %3 dB of passband ripple
        %Rs = 5; % 5 dB attenuation in the stopbands

        [n,Wn] = buttord(Wp,Ws,Rp,Rs);

    %Butter - returns the transfer function coefficients of an nth-order lowpass digital 
    %butterworth filter with normalized cutoff frequency Wn.

        [b,a] = butter(n,Wn);

    %tf2sos - Converts a transfer function representation of a given digital
    %filter to an equivalent second-order section representation

        [sos,g] = tf2sos(b,a);

    %apply filter to shaft accelerometer displacement
        ys = filtfilt(sos,g,cs);
        
    %apply filter to test bracket accelerometer displacement
        yt = filtfilt(sos,g,ct);
        
    %apply filter to shaft laser vibrometer integrated displacement
        ysli = filtfilt(sos,g,csli); %filtered shaft laser integrated
        
%% Fast-Fourier Transform 

%FFT original signal

        Fs = 100000;            % Sampling frequency                    
        T = 1/Fs;               % Sampling period       
        L = length(dselect);    % Length of signal

    %Shaft
        %Compute the fourier transform of the signal
        Ys = fft(cs_norm); 
        
        %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L
        P2s = abs(Ys/L);
        P1s = P2s(1:L/2+1);
        P1s(2:end-1) = 2*P1s(2:end-1);
        
        %Define the frequency domain f and plot the single-sided amplitude spectrum P1       
        fo = Fs*(0:(L/2))/L;    %fourier original

    %Test Bracket
        Yt = fft(ct_norm);

        P2t = abs(Yt/L);
        P1t = P2t(1:L/2+1);
        P1t(2:end-1) = 2*P1t(2:end-1);
        
        fo = Fs*(0:(L/2))/L;

%FFT filtered data

    %Shaft
        Yfs = fft(ys);

        P2sf = abs(Yfs/L); %p2 shaft filtered
        P1sf = P2sf(1:L/2+1);
        P1sf(2:end-1) = 2*P1sf(2:end-1);
        ff = Fs*(0:(L/2))/L;    %fourier filtered
        
        %Test Bracket
        Yft = fft(ys);

        P2tf = abs(Yft/L);
        P1tf = P2tf(1:L/2+1);
        P1tf(2:end-1) = 2*P1tf(2:end-1);
        ff = Fs*(0:(L/2))/L;
        
 %% Plotting
    
%Shaft Displacement
    figure(1)
    subplot (2,1,1)
    plot (time(dselect),-cs)
    title('Unfiltered CRB Shaft Displacement - Accelerometer')
    xlabel('Time (s)')
    ylabel('Displacement (m)')

    
    subplot (2,1,2)
    plot(time,csli)  %Vertical displacement 200um/s/V
    title('Unfiltered CRB Shaft Displacement - Laser Vibrometer Integrated Velocity')
    xlabel('Time (s)')
    ylabel('Displacement (m)')
    
%Shaft Displacement Filtered
    figure(2)
    subplot (2,1,1)
    plot (time(dselect),ys)
    title('Filtered CRB Shaft Displacement - Accelerometer')
    xlabel('Time (s)')
    ylabel('Displacement (m)')

 
    subplot (2,1,2)
    plot(time,ysli)  %Vertical displacement 200um/s/V
    title('Filtered CRB Shaft Displacement - Laser Vibrometer Displacement (Integrated Velocity)')
    xlabel('Time (s)')
    ylabel('Displacement (m)')
     
 %Relative Displacement
    figure(3)
    subplot (2,1,1)
    plot (time(dselect),(ys-yt))
    title('FIltered CRB Relative Displacement - Accelerometer')
    xlabel('Time (s)')
    ylabel('Displacement (m)')

    subplot (2,1,2)
    plot(time(dselect),(ysli-yt))  %Vertical displacement 200um/s/V
    title('Filtered CRB Relative Displacement - Laser Vibrometer Displacement (Integrated Velocity)')
    xlabel('Time (s)')
    ylabel('Displacement (m)')
    
    %Test Bracket Displacement
    figure(4)
    plot (time(dselect),yt)
    title('Filtered CRB Test Bracket Displacement - 12 o-clock Accelerometer')
    xlabel('Time (s)')
    ylabel('Displacement (m)')

    figure(5)
    plot(time,ysli)  %Vertical displacement 200um/s/V
    title('Filtered CRB Shaft Displacement - Laser Vibrometer Displacement (Integrated Velocity)')
    xlabel('Time (s)')
    ylabel('Displacement (m)')
    
    figure(6)
    plot (time(dselect),(ys))
    title('Filtered CRB Shaft Displacement - Accelerometer')
    xlabel('Time (s)')
    ylabel('Displacement (m)')
    
 yr = (ys-yt);
 
 speed = time*3600;
    
    %Export Data
%     inputDATA = [time,as,vs,ys,at,vt,yt,ysli];
%     save inputDATA_CRB_Paper5.mat
save inputDATA_CRB_Paper5.mat time yt ysli at speed vsl


%filter displacement (filter low frequency) high pass filter
%run 0.05 second then integration then filtering then relative then use
%for film thickness eqn.