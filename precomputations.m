function [A, Iz, ycdg, zcdg] = precomputations(t1,t2,h1,h2,b)

% A  - Section area
d=sqrt(b^2+(h1-h2)^2/4);
A=t2*(h1+h2)+2*t1*d;
Ai=[t2*h1; t2*h2; t1*d; t1*d];

% Iz - Section inertia
y=[h1/2; h1/2; h1-(h1-h2)/4 ;(h1-h2)/4];
z=[0; b; b/2; b/2];

sumyA=0;
sumzA=0;
for i=1:size(Ai,1)
    sumyA=sumyA+Ai(i)*y(i);
    sumzA=sumzA+Ai(i)*z(i);
end

ycdg=sumyA/A;
zcdg=sumzA/A;

Izi=[
    1/12*t2*h1^3;
    1/12*t2*h2^3;
    1/12*t1*d^3*((h1-h2)/(2*d))^2+(ycdg-y(3))^2*Ai(3);
    1/12*t1*d^3*((h1-h2)/(2*d))^2+(ycdg-y(4))^2*Ai(4);
    ];

Iz=0;
for i=1:size(Izi,1)
    Iz=Iz+Izi(i);
end
end