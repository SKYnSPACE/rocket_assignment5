clc;
clear all;
close all;

M           = 0.01:0.01:100;
k           = 1.2;

Mach_rel    = 1 + (k-1)/2 * M.^2;                                          % Mach# Relationship

T_over_T0   = Mach_rel.^(-1);
p_over_p0   = Mach_rel.^(k/(1-k));
A_over_At   = 1./M .* sqrt( (2/(k+1) .* Mach_rel).^((k+1)/(k-1)) );

%% Plotting
[ax, h1, h2] = plotyy(M, [T_over_T0; p_over_p0],M, A_over_At,'plot');
grid on;

set(ax(1),'xscale','log','yscale','log',...
    'YLim',[0.01 10],'YTick',[0.01 0.1 1 10]);
set(ax(2),'xscale','log','yscale','log',...
    'YTick',[1 10 100 1000 10000]);
legend('T/T_0', 'p/p_0', 'A/A_t', 'Location', 'northwest');

title('Solutions of ideal nozzle');

xlabel('Mach Number');
set(get(ax(1), 'Ylabel'), 'String', 'Pressure ratio & Temperature ratio');
set(get(ax(2), 'Ylabel'), 'String', 'Area ratio');

% legend();
% grid on;
% hold on;
% loglog(M, p_over_p0);
% loglog(M, A_over_At);
