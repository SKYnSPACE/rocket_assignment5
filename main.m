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
[ax, h1, h2] = plotyy(M, [T_over_T0(1,:); p_over_p0(1,:)],M, A_over_At(1,:),'plot');
hold(ax(1),'on');
hold(ax(2),'on');
set(ax(1),'xscale','log','yscale','log',...
    'XLim',[M(1) M(length(M))],...
    'YLim',[0.01 10],'YTick',[0.01 0.1 1 10]);
set(ax(2),'xscale','log','yscale','log',...
    'XLim',[M(1) M(length(M))],...
    'YLim',[1 1000], 'YTick',[1 10 100 1000]);
legend('T/T_0', 'p/p_0', 'A/A_t', 'Location', 'northwest');
    
for i=2:length(k)
    plot(ax(1),M, T_over_T0(i,:),'color',[0      0.4470 0.7410],'Linewidth',i);
    plot(ax(1),M, p_over_p0(i,:),'color',[0.8500 0.3250 0.0980],'Linewidth',i);
    plot(ax(2),M, A_over_At(i,:),'color',[0.9290 0.6940 0.1250],'Linewidth',i);
    set(ax(1),'xscale','log','yscale','log',...
        'XLim',[M(1) M(length(M))],...
        'YLim',[0.01 10],'YTick',[0.01 0.1 1 10]);
    set(ax(2),'xscale','log','yscale','log',...
        'XLim',[M(1) M(length(M))],...
        'YLim',[1 1000], 'YTick',[1 10 100 1000]);
end
grid on;

title('Solutions of ideal nozzle');

xlabel('Mach Number');
set(get(ax(1), 'Ylabel'), 'String', 'Pressure ratio & Temperature ratio');
set(get(ax(2), 'Ylabel'), 'String', 'Area ratio');

%% Default MATLAB Plotting Color Sets from R2014b
% 'Color', []
%          0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840