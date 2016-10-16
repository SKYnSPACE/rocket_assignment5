%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of ideal nozzle (RSE Presentation #5)
% Seongheon Lee, AE Dept., 20165234
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

M           = 0.01:0.01:100;
k           = [1.2 1.3 1.4];

for i=1:length(k)
    Mach_rel(i,:)    = 1 + (k(i)-1)/2 * M.^2;                                      % Mach# Relationship
    T_over_T0(i,:)   = Mach_rel(i,:).^(-1);
    p_over_p0(i,:)   = Mach_rel(i,:).^(k(i)/(1-k(i)));
    A_over_At(i,:)   = 1./M .* sqrt( (2/(k(i)+1) .* Mach_rel(i,:)).^((k(i)+1)/(k(i)-1)) );
end

%% Plotting
for i=1:length(k)
    [ax, h1, h2] = plotyy(M, [T_over_T0(i,:); p_over_p0(i,:)],M, A_over_At(i,:),'plot');
    hold(ax(1),'on');
    hold(ax(2),'on');
    
    set(ax(1),'xscale','log','yscale','log',...
        'XLim',[M(1) M(length(M))],...
        'YLim',[0.01 10],'YTick',[0.01 0.1 1 10]);
    set(ax(2),'xscale','log','yscale','log',...
        'XLim',[M(1) M(length(M))],...
        'YLim',[1 1000], 'YTick',[1 10 100 1000]);
end
grid on;

legend('T/T_0', 'p/p_0', 'A/A_t', 'Location', 'northwest');

title('Solutions of ideal nozzle');

xlabel('Mach Number');
set(get(ax(1), 'Ylabel'), 'String', 'Pressure ratio & Temperature ratio');
set(get(ax(2), 'Ylabel'), 'String', 'Area ratio');