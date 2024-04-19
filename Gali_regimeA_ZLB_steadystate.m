function [ys,params,check] = Gali_regimeA_ZLB_steadystate(ys,exo,M_,options_)
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] 0 if steady state computation worked and to
%                        1 of not (allows to impose restrictions on parameters)
check=0;
%% Step 1: read out parameters to access them with their name
for ii = 1:M_.param_nbr
  eval([ M_.param_names{ii} ' = M_.params(' int2str(ii) ');']);
end

%% Step 2: Enter model equations here
pi=1; gm=1; pitile=1;
D=1; pistar=1; 
N=1/3;
Y=N^(1-alph);
G=G_Y*Y;
C=Y-G;
R=1/bet;

uc=C^(-thet);
m=m_Y*Y;
um=(R-1)/R*uc;

funa=@(y) Fun_a(y,m,C,sigm,um,chi);
y0=[0.1];
options = optimoptions('fsolve', 'Display', 'none','TolX', 1e-8,'TolFun', 1e-8);
%options = optimoptions('fsolve', 'TolX', 1e-8,'TolFun', 1e-10);
y = fsolve(funa, y0, options);
xbar=y(1);

A=Y/ (N^(1-alph));
mc=(ep-1)/ep;
w=mc*(1-alph)*Y/N;
un=uc*w;
ps=un*N^(-ph);
Z=1;
lamb=Z*uc;
b=Y*b_Y;
T=G+(R-1)*b;
Tstar=G;
x1=lamb*mc*Y/(1-bet*xi);
x2=lamb*Y/(1-bet*xi);
r=R/pi;
%% Step 3: Update parameters and variables
params=NaN(M_.param_nbr,1);
for iter = 1:M_.param_nbr %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

for ii = 1:M_.orig_endo_nbr %auxiliary variables are set automatically
  eval(['ys(' int2str(ii) ') = ' M_.endo_names{ii} ';']);
end
end