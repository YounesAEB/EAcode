function [tau_max,sigma_max, qs0, tausigma_max] = VonMises(Syaux,Mzaux,Iz,b,h1,h2,t1,t2,zcdg)

    % Calculem la tensió normal on és màxima
    sigma_max = (Mzaux/Iz)*(h1/2);
    
    % Angle del perfil theta
    theta = atan((h1-h2)/(2*b));
    d = b/cos(theta);
    
    % Integrals del perfil tancat
    syms s
    
    % Primer tram 12
    q12 = (Syaux/Iz)*int(s*t2);
    q2_open = double(subs(q12,s,h1/2));
  
    % Segon tram 23
    q23 = q2_open +(Syaux/Iz)*int((h1/2-sin(theta)*s)*t1);
    q3_open = double(subs(q23,s,d));
    
    % Tercer tram 34
    q34 = q3_open + (Syaux/Iz)*int((h2/2-s)*t2);
    q4_open = double(subs(q34,s,h2));
    
    % Quart tram 45
    q45 = q4_open +(Syaux/Iz)*int((-h2/2-sin(theta)*s)*t1);
    q5_open = double(subs(q45,s,d));
    
    % Cinquè tram 56
    q56 = q5_open+(Syaux/Iz)*int((-h1/2+s)*t2);
    
    % Moments de cada tram
    rq1 = zcdg*int(q12,0,h1/2);
    rq2 = ((h1/2)-((h1-h2)/4))*cos(theta)*int(q23,0,d);
    rq3 = (b-zcdg)*int(q34,0,h2);
    rq4 = ((-h1/2)+((h1-h2)/4))*cos(theta)*int(q45,0,d);
    rq5 = zcdg*int(q56,0,h1/2);
    rq = double(rq1 + rq2 + rq3 + rq4 + rq5); 
    
    % qs0
    Ain = b*((h1+h2)/2);
    qs0 = rq/(2*Ain);
    
    %Distribucions finals
    
    q12closed = qs0 + q12;
    q23closed = qs0 + q23;
    q34closed = qs0 + q34;
    q45closed = qs0 + q45;
    q56closed = qs0 + q56;
    
    %Sabem com el màxim es troba a q34 o entre q12 i q56:
    qmax12 = double(subs(q12closed,s,0));
    qmax34 = double(subs(q34closed,s,h2/2));

    derivada23=diff(q23closed,s)==0;
    derivada45=diff(q45closed,s)==0;

    qmax23=double(solve(derivada23,s));
    qmax45=double(solve(derivada45,s));
    
    maxq=[qmax12,qmax23,qmax34,qmax45];
    qmax = max(maxq);
    
    if (qmax==qmax12) || (qmax==qmax34)
        tau_max = qmax/t2;
    else
        tau_max = qmax/t1;
    end

    %I trobem tau en el punt de sigma maxim
    tausigma_max = double(subs(q12closed,s,h1/2));
    
end


