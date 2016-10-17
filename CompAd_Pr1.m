%% Compressible Aerodynamics Computer Project #1
clear all; clc;

% Let air with gamma=1.4, perfect gas flows
gamma=1.4;
a_cric=sqrt((gamma-1)/(gamma+1));

% Grid Setting
x_step=0.005;

% Nozzle geometry
P=polyfit([0, 1, 2], [0.75, 0.5, 0.75], 2);
x=0:x_step:2.5;
r=polyval(P,x);
A=pi*r.^2;

% To calculate, qbar(x=0) which reaches critical state at nozzle throat
% (i.e qbar(x=0) with chocked condition)
idx_throat=find(min(A)==A);
qbar_chocked=a_cric-0.00005;
for j=1:idx_throat-1
    qbar_chocked=Conv_Div_Nozzle_RK4(x(idx_throat-j+1),qbar_chocked,-x_step,gamma);
    if length(qbar_chocked)==2
        qbar_chocked=qbar_chocked(1);
    end
end

ini_qbar=[qbar_chocked/3 qbar_chocked/3*2 qbar_chocked];
sol_q=zeros(length(x),length(ini_qbar));
sol_q(1,:)=ini_qbar;
chocked_condition=0;

% 4th_Order_Runge_Kutta method and its solution
for i=1:length(ini_qbar)
    for j=2:length(x)
        if chocked_condition==0
            sol_s=Conv_Div_Nozzle_RK4(x(j-1),sol_q(j-1,i),x_step,gamma);
            if length(sol_s)==2
                sol_q(j,i+1)=sol_s(2);
                sol_q(j,i)=sol_s(1);
                sol_q(1:j-1,i+1)=sol_q(1:j-1,i);
                chocked_condition=1;
            else
                sol_q(j,i)=sol_s;
            end
        elseif chocked_condition==1
            sol_q(j,i)=Conv_Div_Nozzle_RK4(x(j-1),sol_q(j-1,i),x_step,gamma);
            sol_q(j,i+1)=Conv_Div_Nozzle_RK4(x(j-1),sol_q(j-1,i+1),x_step,gamma);
        end
    end
end


% Calculate P/P0 from q_bar
sol_p=zeros(size(sol_q));
for i=1:length(x)
    for j=1:min(size(sol_q))
        sol_p(i,j)=(1-sol_q(i,j)^2)^(gamma/(gamma-1));
    end
end

% Calculate Mach number from q_bar
sol_M=zeros(size(sol_q));
for i=1:length(x)
    for j=1:min(size(sol_q))
        sol_M(i,j)=sol_q(i,j)/sqrt(1-sol_q(i,j)^2)*sqrt(2/(gamma-1));
    end
end

x_after=2.5:0.005:2.6;

% Consider Normal Shock at outlet
idx=find(sol_M==max(max(sol_M)));
M2_NSE=sqrt(((gamma-1)*sol_M(idx)^2+2)/(2*gamma*sol_M(idx)^2-(gamma-1)));
P2_P1NE=(2*gamma*sol_M(idx)^2-(gamma-1))/(gamma+1);
P0_P1NE=(1+(gamma-1)/2*sol_M(idx)^2)^(gamma/(gamma-1));
P2_P0NSE=P2_P1NE/P0_P1NE;
qbar_NSE=M2_NSE/sqrt(2/(gamma-1))*sqrt(1-sol_q(idx)^2)*sqrt((2*gamma*sol_M(idx)^2-(gamma-1))*(2+(gamma-1)*sol_M(idx)^2)/sol_M(idx)^2/(gamma+1)^2);

% Consider Normal Shock inside the nozzle
% Assume normal shock exists inside nozzle at x=1.75
idx=find(x==1.75);
M2_NSI=sqrt(((gamma-1)*sol_M(idx,3)^2+2)/(2*gamma*sol_M(idx,3)^2-(gamma-1)));
P2_P1NI=(2*gamma*sol_M(idx,3)^2-(gamma-1))/(gamma+1);
P20_P10=((1+(gamma-1)/2*M2_NSI^2)/(1+(gamma-1)/2*sol_M(idx,3)^2))^(gamma/(gamma-1))*P2_P1NI;
qbar_NSI=M2_NSI/sqrt(2/(gamma-1))*sqrt(1-sol_q(idx,3)^2)*sqrt((2*gamma*sol_M(idx,3)^2-(gamma-1))*(2+(gamma-1)*sol_M(idx,3)^2)/sol_M(idx,3)^2/(gamma+1)^2);
sol_NSI=zeros(3,length(x)-idx+1);
sol_NSI(1,1)=qbar_NSI; sol_NSI(2,1)=P2_P1NI*sol_p(idx,3); sol_NSI(3,1)=M2_NSI;
for i=2:length(sol_NSI)
    sol_NSI(1,i)=Conv_Div_Nozzle_RK4(x(idx+i-1),sol_NSI(1,i-1),x_step,gamma);
end
for i=2:length(sol_NSI)
    sol_NSI(2,i)=(1-sol_NSI(1,i)^2)^(gamma/(gamma-1))*P20_P10;
    sol_NSI(3,i)=sol_NSI(1,i)/sqrt(1-sol_NSI(1,i)^2)*sqrt(2/(gamma-1));
end

    


% Q_Bar graph
figure(1), set(gcf,'color','w')
plot([x,x_after],[sol_q(:,[1,2,4]);ones(length(x_after),1)*sol_q(end,[1,2,4])], 'linewidth',1.5)
hold on
plot([x,x_after],[sol_q(:,3)',ones(1,length(x_after))*qbar_NSE],'c','linewidth',1.5)
plot([x(1:idx),x(idx:end),x_after],[sol_q(1:idx,3)',sol_NSI(1,:),ones(1,length(x_after))*sol_NSI(1,end)],'m--','linewidth',1.5)
axis([0 2.6 0 1]), grid on
legend('Subsonic(M=0.2 at throat)','Subsonic(M=0.4 at throat)','Chocked(Subsonic)', 'Chocked(Supersonic)+Normalshock(Exit)','Chocked(Supersonic)+NormalShock(Inside)','location','NW')
xlabel('x(m)','fontsize',13), ylabel('q/q_m_a_x','fontsize',13), title('X-q/q_m_a_x Graph','fontsize',15)

% P/P0 graph
figure(2), set(gcf,'color','w')
plot([x,x_after],[sol_p(:,[1,2,4]);ones(length(x_after),1)*sol_p(end,[1,2,4])], 'linewidth',1.5)
hold on
plot([x,x_after],[sol_p(:,3)',ones(1,length(x_after))*P2_P0NSE],'c','linewidth',1.5)
plot([x(1:idx),x(idx:end),x_after],[sol_p(1:idx,3)',sol_NSI(2,:),ones(1,length(x_after))*sol_NSI(2,end)],'m--','linewidth',1.5)
axis([0 2.6 0 1]), grid on
legend('Subsonic(M=0.2 at throat)','Subsonic(M=0.4 at throat)','Chocked(Subsonic)', 'Chocked(Supersonic)+Normalshock(Exit)','Chocked(Supersonic)+NormalShock(Inside)','location','SW')
xlabel('x(m)','fontsize',13), ylabel('P/P_0','fontsize',13), title('X-P/P_0 Graph','fontsize',15)


% M graph
figure(3), set(gcf,'color','w')
plot([x,x_after],[sol_M(:,[1,2,4]);ones(length(x_after),1)*sol_M(end,[1,2,4])], 'linewidth',1.5)
hold on
plot([x,x_after],[sol_M(:,3)',ones(1,length(x_after))*M2_NSE],'c','linewidth',1.5)
plot([x(1:idx),x(idx:end),x_after],[sol_M(1:idx,3)',sol_NSI(3,:),ones(1,length(x_after))*sol_NSI(3,end)],'m--','linewidth',1.5)
axis([0 2.6 0 3.5]), grid on
legend('Subsonic(M=0.2 at throat)','Subsonic(M=0.4 at throat)','Chocked(Subsonic)', 'Chocked(Supersonic)+Normalshock(Exit)','Chocked(Supersonic)+NormalShock(Inside)','location','NW')
xlabel('x(m)','fontsize',13), ylabel('Mach Num','fontsize',13), title('X-M Graph','fontsize',15)

% Nozzle configuration
figure(4), set(gcf,'color','w')
plot(x,r,'linewidth',1.5),hold on, grid on, plot(x,-r,'linewidth',1.5)
xlabel('x(m)','fontsize',13), ylabel('Nozzle Radius(m)','fontsize',13), title('Nozzle Configuration','fontsize',15)