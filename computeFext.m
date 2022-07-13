function [Fext] = computeFext(n_dof,n_el,n_ne,n_i,Td,Me,g,Fel,x,L1)
Fext_puntual=zeros(n_dof,1);
for i=1:size(x,1)
    if (x(i,1)==L1)
        Fext_puntual(2*i-1)=-Me*g;
    else 
        if(x(i,1)>L1) && (x(i-1,1)<L1)
        Fext_puntual(2*(i-1)-1)=-Me*g/2;
        Fext_puntual(2*i-1)=-Me*g/2;
        end
    end
end

Fext=Fext_puntual;
for e=1:n_el
    for i=1:(n_ne*n_i)
        I=Td(e,i);
        Fext(I)=Fext(I)+Fel(i,e);
    end
end

end