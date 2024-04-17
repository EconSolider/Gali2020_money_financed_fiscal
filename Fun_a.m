function F=Fun_a(y,m,C,sigm,um,chi)
xbar=y(1);

F(1)=um*C-chi*(m/C-xbar)^(-sigm);
end