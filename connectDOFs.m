function [Td] = connectDOFs(n_el,n_ne,n_i,Tnod)

Td=zeros(n_el,n_ne*n_i);
for i=1:n_el
    Td(i,1)=Tnod(i,1)*n_i-1;
    Td(i,2)=Tnod(i,1)*n_i;
    Td(i,3)=Tnod(i,2)*n_i-1;
    Td(i,4)=Tnod(i,2)*n_i;
end
