function val=Conv_Div_Nozzle_RK4(x,dqdx,del_x,gamma)

% RK-4 coefficient
a1=1/6; a2=1/3; a3=1/3; a4=1/6;
% Critical aucoustic speed
a_cric=sqrt((gamma-1)/(gamma+1));

k1=dqbar_dx(dqdx,x,gamma);
k2=dqbar_dx(dqdx+k1*del_x/2,x+del_x/2,gamma);
k3=dqbar_dx(dqdx+k2*del_x/2,x+del_x/2,gamma);
k4=dqbar_dx(dqdx+k3*del_x,x+del_x,gamma);
val1=dqdx+(a1*k1+a2*k2+a3*k3+a4*k4)*del_x;

if abs(dqdx-a_cric)/a_cric<0.001
    k1=dqbar_dx_throat_critical(dqdx,x,gamma);
    k2=dqbar_dx_throat_critical(dqdx+k1*del_x/2,x+del_x/2,gamma);
    k3=dqbar_dx_throat_critical(dqdx+k2*del_x/2,x+del_x/2,gamma);
    k4=dqbar_dx_throat_critical(dqdx+k3*del_x,x+del_x,gamma);
    val2=dqdx+(a1*k1+a2*k2+a3*k3+a4*k4)*del_x;
    val=[val1; val2];
else
    val=val1;
end

end


function s=dqbar_dx(qbar,x,gamma)
a_cric=sqrt((gamma-1)/(gamma+1));
if abs(qbar-a_cric)/a_cric<0.001
    s=(a_cric^2/(1-a_cric^2)^2/2/An(x)*d2An_dx2(x))^0.5;
else
    s=(-1)*(1-qbar^2)*qbar/An(x)/(1-(gamma+1)/(gamma-1)*qbar^2)*dAn_dx(x);
    
end
end

function s=dqbar_dx_throat_critical(qbar,x,gamma)
a_cric=sqrt((gamma-1)/(gamma+1));
if abs(qbar-a_cric)/a_cric<0.001
    s=-(a_cric^2/(1-a_cric^2)^2/2/An(x)*d2An_dx2(x))^0.5;
else
    s=(-1)*(1-qbar^2)*qbar/An(x)/(1-(gamma+1)/(gamma-1)*qbar^2)*dAn_dx(x);
    
end
end

function y=An(x)
P=polyfit([0, 1, 2], [0.75, 0.5, 0.75], 2);
r=polyval(P,x);
y=pi*r.^2;
end

function y=dAn_dx(x)
P=polyfit([0, 1, 2], [0.75, 0.5, 0.75], 2);
Rsq=conv(P,P);
y=pi*polyval(polyder(Rsq),x);
end

function y=d2An_dx2(x)
P=polyfit([0, 1, 2], [0.75, 0.5, 0.75], 2);
Rsq=conv(P,P);
y=pi*polyval(polyder(polyder(Rsq)),x);
end
