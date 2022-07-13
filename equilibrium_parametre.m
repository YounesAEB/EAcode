function [l] = equilibrium_parametre(L1,L2,M,Me,g)

syms x
lambda1= M/(4*(L1+L2))+3*M/(2*L2^2)*(L1-x);
lambda2= M/(4*(L1+L2));
q1= (0.8-0.2*cos(pi*x/L1));
q2= (1-(x-L1)/L2)*(1+(x-L1)/L2);

L=L1+L2;
TW=(Me+int(lambda1,x,0,L1)+int(lambda2,x,L1,L))*g;
TW=round(TW,10);
Lp=int(q1,x,0,L1)+int(q2,x,L1,L);
Lp=round(Lp,10);
l=TW/Lp;
end