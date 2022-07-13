function fig = plotBeamsInitialize(L,x_an,u_an,theta_an,Fy_an,Mz_an)
% PLOTBEAMSINITIALIZE - Initialize plots for beams
% If arguments are given, then it plots the analytical solution
% Inputs:
%  L         Beam's total length
%  x_an      Vector with the beam's coordinates
%  u_an      Vector with the beam's displacements (analytical)
%  theta_an  Vector with the beam's section rotation (analytical)
%  Fy_an     Vector with the beam's shear force (analytical)
%  Mz_an     Vector with the beam's bending moment (analytical)

% Open figure window
fig = figure('units','centimeters','Position',[2,2,16,16]);
% Plot beam deflection
subplot(2,2,1)
if nargin > 1
    plot(x_an,u_an(:),'-k','linewidth',1.5);
end
hold on
box on
xlim([0,L]);
xlabel('x (m)','Interpreter','latex');
ylabel('$U_y$ (m)','Interpreter','latex');
title('Deflection','Interpreter','latex');
% Plot beam section rotation
subplot(2,2,2)
if nargin > 1
    plot(x_an,theta_an,'-k','linewidth',1.5);
end
hold on
box on
xlim([0,L]);
xlabel('x (m)','Interpreter','latex');
ylabel('$\theta_z$ (rad)','Interpreter','latex');
title('Rotation','Interpreter','latex');
% Plot beam internal shear force
subplot(2,2,3)
if nargin > 1
    plot(x_an,Fy_an,'-k','linewidth',1.5);
end
hold on
box on
xlim([0,L]);
xlabel('x (m)','Interpreter','latex');
ylabel('$F_y$ (N)','Interpreter','latex');
title('Shear force','Interpreter','latex');
% Plot beam internal bending moment
subplot(2,2,4)
if nargin > 1
    plot(x_an,Mz_an,'-k','linewidth',1.5);
end
hold on
box on
xlim([0,L]);
xlabel('x (m)','Interpreter','latex');
ylabel('$M_z$ (Nm)','Interpreter','latex');
title('Bending moment','Interpreter','latex');

end