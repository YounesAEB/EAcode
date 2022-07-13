function [pu,pt,Fy,Mz]=lastcomputations(n_el,n_i,n_ne,Tnod,Td,u,Kel,x)
    
    ue=zeros(n_i*n_ne,1);
    Fy=zeros(n_el,n_ne);
    Mz=zeros(n_el,n_ne);
    pu=zeros(n_el,4);
    pt=zeros(n_el,3);

for e=1:n_el
    x1e=x(Tnod(e,1),1);
    x2e=x(Tnod(e,2),1);
    l=abs(x2e-x1e);

    for j=1:(n_ne*n_i)
        I=Td(e,j);
        ue(j,1)=u(I);
    end
    
    %Internal forces, shear force and bending moment al element nodes
    Fe_int=Kel(:,:,e)*ue;
    Fy(e,1)=-Fe_int(1);
    Fy(e,2)=Fe_int(3);
    Mz(e,1)=-Fe_int(2);
    Mz(e,2)=Fe_int(4);

    %Third-order polynomial coefficients for element's deflection and
    %section rotations

    abcd=1/l^3*[2 l -2 l;
        -3*l -2*l^2 3*l -l^2;
        0 l^3 0 0;
        l^3 0 0 0;]*ue;

    pu(e,:)=abcd;
    pt(e,:)=[3*abcd(1),2*abcd(2),abcd(3)];


end

end