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
;

G_Y=0.2;
b_Y=1;
m_Y=2;
bet=0.96;
psib=0.1; 
delt=0.9;
alph=0.35; 
ep=8; 
xi=0.5; 
iot=0; 
chi=0.1; 
xbar=0.1;  %recalibrate
sigm=1;
ps=0.1;    %recalibrate
ph=0.5;    %inverse frich 



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
   delt * (Tstar(-1)-steady_state(Tstar))/steady_state(Y) - epsT;
%5b
pi=steady_state(pi);
%6
m/m(-1)*pi=gm;
%7
um/uc=(R-1)/R;
%8
Z(+1)/Z*uc(+1)/uc=pi/bet/R;
%9
un/uc=w;
%10
w=mc*(1-alph)*Y/N;
%11
pistar=ep/(ep-1)*x1/x2;
%12
x1=lamb*mc*Y+bet*xi*(pitile/pi)^(-ep)*x1;
%13
x2=lamb*Y+bet*xi*(pitile/pi)^(1-ep)*x2;
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
  delt*log(Z(-1)/steady_state(Z))+epsZ;
%21
lamb=Z*uc;
%22
un=ps*N^ph;
%23
uc=1/C;
%24
um=chi*(m/C-xbar)^(-sigm)/C;
%25
r=R/pi;
end;

steady;
check;

shocks;
var epsT=1;
end;

stoch_simul(order=1,irf=16,periods=0)
T m b Y C R pi r 
;