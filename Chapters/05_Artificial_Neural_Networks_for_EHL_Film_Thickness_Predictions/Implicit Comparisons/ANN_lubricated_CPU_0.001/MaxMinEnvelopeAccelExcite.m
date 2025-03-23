drymax = [2 12 8 10 12 18 28 40 62 100 210 500 410 300 230 210 210 220 245 280 320 400 480 600 745];
drymin = -drymax;

lubmax = [2 12 8 9 12 18 26 40 60 110 210 470 420 285 240 220 215 220 240 270 308 370 450 560 695];
lubmin = -lubmax;


speed = [1000:1000:25000];


figure(1)
plot(speed,drymax,'r',speed,lubmax,'b',speed,drymin,'r',speed,lubmin,'b')
set(gca,'FontSize',24)
title('EXCITE Results - Dry vs Lubricated Speed Envelope DOF2 5um Preload')
xlabel ('Speed (rpm)')
ylabel ('Acceleration (m/s^2)')
legend ('Dry MATLAB','Lubricated MATLAB')

% figure(2)
% plot(speed,lubmax,'g',speed,lubmin,'g')
% set(gca,'FontSize',24)
% title('Lubricated MATLAB Speed Envelope')
% xlabel ('Speed (rpm)')
% ylabel ('Acceleration (m/s^2)')
% legend ('Lubricated MATLAB')