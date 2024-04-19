var
G b R pi T m Tstar gm um uc
Z un w A mc N pistar x1 x2 pitile
D Y C lamb r
;

varexo
epsA epsG epsT epsZ
;

parameters
bet psib delt alph ep xi iot chi xbar sigm ps ph G_Y b_Y m_Y
deltz thet
;

G_Y=0;
b_Y=1;
m_Y=3;
bet=0.99;
psib=0.1;
delt=0.5;
alph=0.35; 
ep=9; 
xi=0.75; 
iot=0; 
chi=0.018; 
xbar=0.1;  %recalibrate
sigm=2.5;
ps=0.1;    %recalibrate
ph=1.5;    %inverse frich 
deltz=0;
thet=3.5;

model;
%1
G + b(-1)*R(-1)/pi= T + b + m-m(-1)/pi;
%2
(T-steady_state(T))/steady_state(Y)=
    psib*(b(-1)-steady_state(b))/steady_state(Y)
    +(Tstar-steady_state(Tstar))/steady_state(Y);
%3
(G-steady_state(G))/steady_state(Y)=
   delt * (G(-1)-steady_state(G))/steady_state(Y) + epsG;
%4
(Tstar-steady_state(Tstar))/steady_state(Y)=
   delt * (Tstar(-1)-steady_state(Tstar))/steady_state(Y) + epsT;
%5a
m-m(-1)/pi=G+(R(-1)/pi-1)*b-T;
%6
m/m(-1)*pi=gm;
%7 ZLB
% [name = 'Nominal Interest Rate',mcp='R>1']
% um/uc=(R-1)/R;
R=max(1,1/(1-um/uc));
%8
lamb=bet*lamb(+1)*R/pi(+1)/Z;
%9
un/uc=w;
%10
w=mc*(1-alph)*Y/N;
%11
pistar=ep/(ep-1)*x1/x2;
%12
x1=lamb*mc*Y+bet*xi*(pitile(+1)/pi(+1))^(-ep)*x1(+1);
%13
x2=lamb*Y+bet*xi*(pitile(+1)/pi(+1))^(1-ep)*x2(+1);
%14
pitile=pi(-1)^iot*steady_state(pi)^(1-iot);
%15
1=xi*(pi/pitile)^(ep-1)+(1-xi)*pistar^(1-ep);
%16
D*Y=A*N^(1-alph);
%17
D=(1-xi)*pistar^(-ep)+xi*(pi/pitile)^ep*D(-1);
%18
Y=C+G;
%19
log(A/steady_state(A))=
  delt*log(A(-1)/steady_state(A))+epsA;
%20
log(Z/steady_state(Z))=
  deltz*log(Z(-1)/steady_state(Z))+epsZ;
%21
lamb=uc;
%22
un=ps*N^ph;
%23
uc=C^(-thet);
%24
um=chi*(m/C-xbar)^(-sigm)/C;
%25
r=R/pi(+1);
end;


steady;
check;

shocks;
 var epsZ;
 periods 1:5;
 values -0.45;

% var epsT;
% periods 1:5;
% values -1;
end;

perfect_foresight_setup(periods=1000);
perfect_foresight_solver(lmmcp);

for i=1:length(M_.endo_names)
assignin('base', M_.endo_names{i}, oo_.endo_simul(i,:));
end

leg=15;
start=2;

figure;
subplot(3,3,1);
plot(R(1:leg));
title('R');
xlim([start,leg]);

subplot(3,3,2);
plot(T(1:leg));
title('T');
xlim([start,leg]);

subplot(3,3,3);
plot(C(1:leg));
title('C');
xlim([start,leg]);

subplot(3,3,4);
plot(Y(1:leg));
title('Y');
xlim([start,leg]);

subplot(3,3,5);
plot(pi(1:leg));
title('\pi');
xlim([start,leg]);

subplot(3,3,6);
plot(Z(1:leg));
title('Z');
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
