drymax = [890 1150 1250 1175 1050 880 775 700 625];
drymin = [-890 -1150  -1250 -1175 -1050 -880 -775 -700 -625];

lubmax = [650 850 1150 1400 1500 1450 1300 1150 950];
lubmin = [-650 -850 -1150 -1400 -1500 -1450 -1300 -1150 -950];

speed = [11000 11250 11500 11750 12000 12250 12500 12750 13000];


figure(1)
plot(speed,drymax,'r',speed,lubmax,'b',speed,drymin,'r',speed,lubmin,'b')
set(gca,'FontSize',24)
title('Dry vs Lubricated Acceleration Speed Envelope')
xlabel ('Speed (rpm)')
ylabel ('Acceleration (m/s^2)')
legend ('Dry MATLAB','Lubricated MATLAB')

