function [x] = nodal_coordinates(n_nod,n_d,n_el,L)
    x=zeros(n_nod,n_d);
    x(1)=0;
    pos=x(1);
    for i=2:size(x,1)
        x(i)=pos+L/n_el;
        pos=x(i);
    end
end