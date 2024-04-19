%Gali_RegimeA liear


var
y c g xi i pi rho_hat mu m gm b tstar r;

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
delt=0.5;
psib=0.02;


model(linear);
y=c+g;
xi=xi(+1)+(i-pi(+1)-rho_hat);
xi=-sigm*c+nu*m;
pi=bet*pi(+1)-lamb*mu;
mu=xi-(alph+ph)/(1-alph)*y;
m=c-et*i;
m(-1)=m+pi-gm;
b=(1+rho-psib)*b(-1)+bs*(1+rho)*(i(-1)-pi)+g-tstar-chi*gm;
gm=1/chi*(g-tstar+bs*(1+rho)*(i(-1)-pi));
r=i-pi(+1);

g=delt*g(-1)+epsg;
tstar=delt*tstar(-1)-epst;
rho_hat=epsz;
end;

steady;
check;

shocks;
var epst=1;
end;

stoch_simul(order=1,irf=16,periods=0)
tstar m b y c pi i gm r
;
