function [KG] = assemblyKG(n_dof,n_el,n_ne,n_i,Td,Kel)

KG=zeros(n_dof, n_dof);
for e=1:n_el
    for i=1:(n_ne*n_i)
        I=Td(e,i);
        for j=1:(n_ne*n_i)
            J=Td(e,j);
            KG(I,J)=KG(I,J)+Kel(i,j,e);
        end
    end

end
end
