%Gali_RegimeA liear zlb


var
y c g xi ni pi rho_hat mu m gm b tstar r;

varexo epsg epst epsz;

parameters
sigm nu bet thet alph ep lamb ph et rho bs chi delt psib
;

sigm=1.5;
nu=0;
bet=0.995;

thet=3/4;
alph=0.25;
ep=9;
lamb=(1-thet)*(1-bet*thet)*(1-alph)/thet/(1-alph+alph*ep);
ph=5;
et=7;
rho=1/bet-1;
bs=2.4;
chi=3;
delt=0;
psib=0.02;


model;
y=c+g;
xi=xi(+1)+(ni-pi(+1)-rho_hat);
xi=-sigm*c+nu*m;
pi=bet*pi(+1)-lamb*mu;
mu=xi-(alph+ph)/(1-alph)*y;
%m=c-et*ni;
ni=max(-0.05,(c-m)/et);
m(-1)=m+pi-gm;
b=(1+rho-psib)*b(-1)+bs*(1+rho)*(ni(-1)-pi)+g-tstar-chi*gm;
gm=1/chi*(g-tstar+bs*(1+rho)*(ni(-1)-pi));

g=delt*g(-1)+epsg;
tstar=delt*tstar(-1)+epst;
rho_hat=epsz;
r=ni-pi(+1);
end;

steady;
check;

shocks;
 var epsz;
 periods 1:5;
 values -2;

var epst;
periods 1:5;
values -1;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;


for i=1:length(M_.endo_names)
assignin('base', M_.endo_names{i}, oo_.endo_simul(i,:));
end

leg=15;
start=2;

figure;
subplot(3,3,1);
plot(ni(1:leg));
title('i');
xlim([start,leg]);

subplot(3,3,2);
plot(tstar(1:leg));
title('tstar');
xlim([start,leg]);

subplot(3,3,3);
plot(c(1:leg));
title('c');
xlim([start,leg]);

subplot(3,3,4);
plot(y(1:leg));
title('y');
xlim([start,leg]);

subplot(3,3,5);
plot(pi(1:leg));
title('\pi');
xlim([start,leg]);

subplot(3,3,6);
plot(rho_hat(1:leg));
title('rho_hat');
xlim([start,leg]);

subplot(3,3,7);
plot(gm(1:leg));
title('gm');
xlim([start,leg]);

subplot(3,3,8);
plot(r(1:leg));
title('r');
xlim([start,leg]);

subplot(3,3,9);
plot(m(1:leg));
title('m');
xlim([start,leg]);
