clc
clear all

load 'Wt_out_EHL_dry.mat'
Wt_out_dry = Wt_out;
load 'Wt_out_EHL_lub.mat'
Wt_out_lub = Wt_out;

for i = 1:417000
    
if Wt_out_dry(i) == 0
 Wt_out_dry(i) = NaN;
end

if Wt_out_lub(i) == 0
 Wt_out_lub(i) = NaN;
end

end


Wt_diff = Wt_out_lub-Wt_out_dry;
    
Wt_diff_perc = (Wt_out_lub- Wt_out_dry)./Wt_out_dry;


percdiff = Wt_diff_perc*100;

total = size(Wt_out);

for i = 1:417000
if Wt_diff(i) == 0
 Wt_diff(i) = NaN;
end
end

kb = 100000;
kf = kb;

%moving averages

Wt_diff_fit = movmean(Wt_diff,[kb kf],'omitnan'); %omits NaN values

Wt_diff_perc_fit = movmean(Wt_diff_perc,[kb kf],'omitnan');

%% Plotting

    
figure(1)
plot(speed,Wt_diff_perc_fit,'k','LineWidth',1)
set(gca,'FontSize',24)
title ('Percentage Contact Load Difference Between Dry and Lubricated Model - Moving Average')
xlabel ('Rotational Speed/ rpm')
ylabel ('Percentage Difference')
xlim ([0 15000])

figure(2)
plot(speed,Wt_diff_perc,'k','LineWidth',1)
set(gca,'FontSize',24)
title ('Percentage Contact Load Difference Between Dry and Lubricated Model')
xlabel ('Rotational Speed/ rpm')
ylabel ('Percentage Difference')
xlim ([0 15000])

figure(3)
plot(speed,Wt_diff_fit,'k','LineWidth',1)
set(gca,'FontSize',24)
title ('Contact Load Difference Between Dry and Lubricated Model - Moving Average')
xlabel ('Rotational Speed/ rpm')
ylabel ('Load Difference/ N')
xlim ([0 15000])


figure(4)
plot(speed,Wt_diff,'k','LineWidth',1)
set(gca,'FontSize',24)
title ('Contact Load Difference Between Dry and Lubricated Model')
xlabel ('Rotational Speed/ rpm')
ylabel ('Load Difference/ N')
xlim ([0 15000])

figure(5)
plot(speed,Wt_out_lub,'--k','LineWidth',1)
hold on
plot(speed,Wt_out_dry,'k','LineWidth',1)
set(gca,'FontSize',24)
title ('Contact Load Difference Between Dry and Lubricated Model')
xlabel ('Speed/ rpm')
ylabel ('Load Difference/ N')
xlim ([0 15000])