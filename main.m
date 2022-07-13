%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date: 28/03/2022
% Author/s: Younes Akhazzan, Ricard Arbat
%

clear;
close all;
clc;

%% INPUT DATA

% Material properties
E = 85e9;

% Cross-section parameters
t1 = 1.5e-3;
t2 = 4e-3;
h1 = 500e-3;
h2 = 250e-3;
b = 775e-3;


% Other data
g = 9.81;
L1 = 5;
L2 = 10;
L = L1+L2;
Me = 2550;
M = 35000;

% Number of elements for each part
nel = [3,6,12,24,48,96];
Nel = 96; % Number of elements for the "exact" solution

%% PRECOMPUTATIONS

% Compute section: 
% A  - Section area
% Iz - Section inertia
[A, Iz, ycdg, zcdg] = precomputations(t1,t2,h1,h2,b);

% Compute parameter l:
% l - Equilibrium parameter
l=equilibrium_parametre(L1,L2,M,Me,g);

% Plot analytical solution
fig = plotBeamsInitialize(L1+L2);

%Constant dimensions
    n_d = 1;                      % Problem dimension.
    n_ne= 2;                      % Number of nodes in a beam.
    n_i = 2;                      % Number of DOFs for each node (deflection and section rotation)

% Boundary conditions
fixNod = [
            1, 1, 0;
            1, 2, 0;
];

% Function to study the convergence
er=zeros(1,length(nel)-1);

% Loop through each of the number of elements
for k = 1:length(nel)

    %% PREPROCESS
    
    % Variable dimensions
    n_el = nel(k);                % Total number of beams.
    n_nod = nel(k)+1;             % Total number of nodes.
    n_dof = n_nod*n_i;            % Total number of degrees of freedom

    % Nodal coordinates
    %  x(a,j) = coordinate of node a in the dimension j
    x= nodal_coordinates(n_nod,n_d,n_el,L);

    % Nodal connectivities  
    %  Tnod(e,a) = global nodal number associated to node a of element e
    Tnod= nodal_connectivities(n_el,n_ne);

    % Material properties matrix
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
    %  mat(m,3) = Section inertia of material m
    mat = [% Young M.        Section A.    Inertia 
              E,                A,         Iz;  % Material (1)
    ];

    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    Tmat=ones(n_el,1);   
    %% SOLVER
    
    % Computation of the DOFs connectivities
    Td = connectDOFs(n_el,n_ne,n_i,Tnod);

    % Computation of element stiffness matrices
    [Kel] = computeKel(n_ne,n_i,n_el,x,Tnod,mat,Tmat);

    % Global matrix assembly
    [KG] = assemblyKG(n_dof,n_el,n_ne,n_i,Td,Kel);
    
    % Computation of the element force vector
    [Fel] = computeFel(n_ne, n_i, n_el, Tnod, x, M, L1, L2, g, l);

    % Computation of the external forces vector
    [Fext] = computeFext(n_dof,n_el,n_ne,n_i,Td,Me,g,Fel,x,L1);

    % Apply conditions 
    [vL,vR,uR] = applyCond(n_dof,fixNod);

    % System resolution
    [u,R] = solveSys(vL,vR,uR,KG,Fext);


    % Compute:
    % u  - Displacements and rotations vector [ndof x 1]
    % pu - Polynomial coefficients for displacements for each element [nel x 4]
    % pt - Polynomial coefficients for rotations for each element [nel x 3]
    % Fy - Internal shear force at each elements's nodes [nel x nne]
    % Mz - Internal bending moment at each elements's nodes [nel x nne]

    [pu,pt,Fy,Mz]=lastcomputations(n_el,n_i,n_ne,Tnod,Td,u,Kel,x);

    
    %% POSTPROCESS
    
    % Number of subdivisions and plots
    nsub = Nel/nel(k);
    plotBeams1D(fig,x,Tnod,nsub,pu,pt,Fy,Mz)
    drawnow;

    % Study of the convergence
    load('uREF.mat'); %Deflection of the tip with nel=96;
    if (k<length(nel))
        er(k)=abs(u(length(u)-1)-uREF)/uREF;
    end
    
end

 %% VON MISES
   
  tau=zeros(size(Mz,1)+1,1);
  taumax=zeros(size(Mz,1)+1,1);
  sigma=zeros(size(Mz,1)+1,1);
  sVM=zeros(size(Mz,1)+1,1);
    
    for j=1:(size(Mz,1)+1)

        if (j<=size(Mz,1))
            Mzaux = Mz(j,1);
            Syaux = Fy(j,1);
        else
            Mzaux = Mz(j-1,2);
            Syaux = Fy(j-1,2);
        end
        
        [tau_max,sigma_max, qs0,tausigma_max] = VonMises(Syaux,Mzaux,Iz,b,h1,h2,t1,t2,zcdg);
        tausigma_max = tausigma_max*1e-6; %MPa
        sigma_max = sigma_max*1e-6; %MPa
        tau_max = tau_max*1e-6; %MPa
        tau(j,1)=tausigma_max;
        taumax(j,1) = tau_max;
        sigma(j,1)=sigma_max;
        sVM(j,1)=sqrt(sigma(j,1)^2+3*tau(j,1)^2);
    end
 
    %% PLOTS
%Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel))),'location','northeast','Interpreter','latex');

% Plot of the convergence
figure()
semilogx(nel(1:(length(nel)-1)),er(:)*100)
title('Study of the convergence','Interpreter','latex')
xlabel('Number of elements (logarithmic scale)','Interpreter','latex')
ylabel('Relative error (\%)','Interpreter','latex')
grid on
grid minor

% Plot of the normal and shear max stresses
figure()
position=0:15/96:15;
plot(position,taumax,'b',position,sigma,'r')
title('Maximum shear and bending stress','Interpreter','latex')
xlabel('Position of the section (m)','Interpreter','latex')
ylabel('Stress (MPa)','Interpreter','latex')
legend('Shear','Bending','Interpreter','latex')
grid on
grid minor

% Plot of the Von Mises stress
figure()
position=0:15/96:15;
plot(position,sVM,'g')
title('Von Mises stress','Interpreter','latex')
xlabel('Position of the section (m)','Interpreter','latex')
ylabel('Stress (MPa)','Interpreter','latex')
legend('Von Mises','Interpreter','latex')
grid on
grid minor

%% Coses varies
VonMisesMax = max(sVM);