function [u,R] = solveSys(vL,vR,uR,KG,Fext)

KLL=KG(vL,vL);
KLR=KG(vL,vR);
KRL=KG(vR,vL);
KRR=KG(vR,vR);

FL_ext=Fext(vL,1);
FR_ext=Fext(vR,1);

uL=KLL\(FL_ext-KLR*uR);
R=KRR*uR+KRL*uL-FR_ext;

u(vL,1)=uL;
u(vR,1)=uR;

end