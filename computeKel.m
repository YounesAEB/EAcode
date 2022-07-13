function [Kel] = computeKel(n_ne,n_i,n_el,x,Tnod,mat,Tmat)

Kel=zeros(n_ne*n_i, n_ne*n_i, n_el);
for e=1:n_el
    x1=x(Tnod(e,1),1);
    x2=x(Tnod(e,2),1);
    le=abs(x2-x1);
    Ize=mat(Tmat(e),3);
    Ee=mat(Tmat(e),1);
    
    Kel(:,:,e)=Ize*Ee/le^3*[
        12 6*le -12 6*le;
        6*le 4*le^2 -6*le 2*le^2;
        -12 -6*le 12 -6*le;
        6*le 2*le^2 -6*le 4*le^2];
end

end