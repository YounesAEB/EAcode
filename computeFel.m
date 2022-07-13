function [Fel] = computeFel(n_ne, n_i, n_el, Tnod,x, M, L1, L2, g, l)

    syms t
lambda1= M/(4*(L1+L2))+3*M/(2*L2^2)*(L1-t);
lambda2= M/(4*(L1+L2));
q1= l*(0.8-0.2*cos(pi*t/L1));
q2= l*(1-(t-L1)/L2)*(1+(t-L1)/L2);

Fel=zeros(n_ne*n_i, n_el);
for e=1:n_el

    x1=x(Tnod(e,1),1);
    x2=x(Tnod(e,2),1);
    le=abs(x2-x1);

    if (x1<=L1) && (x2<=L1)
        qe=(int(q1,t,x1,x2)-int(lambda1,t,x1,x2)*g)/le;
    else 
        if (x1>=L1) && (x2>=L1)
            qe=(int(q2,t,x1,x2)-int(lambda2,t,x1,x2)*g)/le;
        else 
            if (x1<L1) && (x2>L1)
            qe=(int(q1,t,x1,L1)+int(q2,t,L1,x2)-int(lambda1,t,x1,L1)*g-int(lambda2,t,L1,x2)*g)/le;
            end
        end
    end
    Fel(:,e)=qe*le/2*[1; le/6; 1; -le/6;];
end

end