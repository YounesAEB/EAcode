function [vL,vR,uR] = applyCond(n_dof,fixNod)
vR=zeros(size(fixNod,1),1);
uR=ones(size(fixNod,1),1);

for i=1:size(fixNod,1)    
    if fixNod(i,2)==1
       vR(i)=2*fixNod(i,1)-1;
       uR(i)=fixNod(i,3);
    end

    if fixNod(i,2)==2
       vR(i)=2*fixNod(i,1);
       uR(i)=fixNod(i,3);
    end
end

vL=zeros(n_dof-size(fixNod,1),1);

suma=1;
for i=1:n_dof
    cont=0;
    for j=1:size(vR)
        if i==vR(j)
           cont=cont+1;
        end
    end
    if cont==0 
       vL(suma)=i;
       suma=suma+1;
    end
end
end