clc
clear all

load 'outputTribological'

global Rrout length alpha Er neta0

u = u_EHL';
w = Wt_out_EHL;

total = size(u);


% Landa = [1.065843621	1.056325172	1.046193094	1.035448626	1.024093968	1.012132289	0.999567726	0.986405387	0.972651348	0.958312646];
% w = [460 460 460 460 460 460 460 460 460 460];

for i = 1 :total
    w_current= w(i);
    u_current= u(i);
    
    [Pht,ht,a,X] = Main_Transient(u_current,w_current);

    Pressure_h(:,i) = Pht(:);
    Film_h(:,i) = ht(:);
    
   
end

Domain= X.*a;

save 'Pressure_h'
save 'Film_h'

for k = 1:total
    
    
    figure(k)
    
    left_color = [0 0 0];
    right_color = [0 0 0];
    set(figure(k),'defaultAxesColorOrder',[left_color; right_color]);
    ax = gca;
    ax.FontSize = 24;
    

    yyaxis left
    plot(Domain*1e3,Pressure_h(:,k),'k','LineWidth',2)
    ylabel('Pressure (GPa)')
    ylim([0 1.2])
    
    yyaxis right
    plot(Domain*1e3,Film_h(:,k)*1e6,'--k','LineWidth',2)
    ylabel('Film Thickness (\mum)','color', 'k')
    ylim([0 8])
    
    xlabel('Direction of Entraining Motion (mm)')
    title(sprintf('Node %g',(k)), 'fontSize', 24)
    xlim([-0.3 0.2])
    
    annotation(figure(k),'ellipse',...
    [0.582770833333332 0.20603813559322 0.0250416666666667 0.0524706123564786],...
    'LineWidth',2);
    
end



%%

% Pressure_h = Pht;
% save (['Pressure_h_',num2str(i)])



