function [Tn] = nodal_connectivities(n_el,n_ne)
    Tn=zeros(n_el,n_ne);
    for i=1:size(Tn,1)
        Tn(i,1)=i;
        Tn(i,2)=i+1;
    end
end