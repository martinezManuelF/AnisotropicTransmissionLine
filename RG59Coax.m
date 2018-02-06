%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILENAME:         RG59Coax.m
% COURSE:           EE5322--21st Century Electromagnetics
% INSTRUCTOR:       Raymond C. Rumpf
% NAME:             Manuel F. Martinez
% SEMESTER:         Spring 2018
% DUE DATE:         02/06/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESTORE MATLAB STATE
clear all;
close all;
clc;

% UNITS
meters      = 1;
seconds     = 1;
degrees     = pi/180;
F           = 1;
H           = 1;

% CONSTANTS
e0 = 8.85418782e-12 * F/meters;
u0 = 1.25663706e-6 * H/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

% OPEN FIGURE WINDOW
figure('Color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TRANSMISSION LINE PARAMETERS
erin  = 2.3 * eye(3,3);     % Permittivity tensor of inside of Coax cable
erout = 1.0 * eye(3,3);     % Permittivity tensot of outside of Coax cable
r1    = 0.35;               % Radius of inner conductor
r2    = 2.5;                % Radius of dielectric core
r3    = 2.7;                % Thickness of cladding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SPACER REGIONS
BUFF  = 3*r3;
Sx    = BUFF;
Sy    = BUFF;

% GRID SIZE
Nx = 512;
Ny = 512;

% INITIAL GUESS AT RESOLUTION
dx = Sx/Nx;
dy = Sy/Ny;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(r1/dx);
dx = r1/nx;
ny = ceil(r1/dy);
dy = r1/ny;

% COMPUTE 2X GRID
Nx2 = 2*Nx;
dx2 = dx/2;
Ny2 = 2*Ny;
dy2 = dy/2;

% GRID AXES
xa = [0:Nx-1]*dx; xa = xa - mean(xa);
ya = [0:Ny-1]*dy; ya = ya - mean(ya);

% 2x GRID AXES
xa2 = [0:Nx2-1]*dx2; xa2 = xa2 - mean(xa2);
ya2 = [0:Ny2-1]*dy2; ya2 = ya2 - mean(ya2);

% CREATE MESH
[Y,X] = meshgrid(ya,xa);
RSQ = (X.^2 + Y.^2);

% CREATE 2X MESH
[Y2,X2] = meshgrid(ya2,xa2);
RSQ2 = (X2.^2 + Y2.^2);
CIN = RSQ2 < ((r3 + r2)/2)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE SIGNALS
SIG.V = [0 1 0];
SIG.GND = zeros(Nx,Ny);
SIG.SIG1 = SIG.GND;
SIG.SIG2 = SIG.GND;

% FORCE CONDUCTORS
SIG.GND(:,[1 Ny]) = 1;
SIG.GND([1 Nx],:) = 1;
SIG.SIG1 = RSQ <= r1^2;
SIG.SIG2 = RSQ >= r2^2 & RSQ < r3^2;

% BUILD PERMITTIVITIES ON 2X GRID
% Outside dielectrics
A = RSQ2 <= r1^2;
ER2xx = erout(1,1) * (1 - CIN);
ER2xy = erout(1,2) * (1 - CIN);
ER2yx = erout(2,1) * (1 - CIN);
ER2yy = erout(2,2) * (1 - CIN);
% Inside dielectrics
ER2xx = ER2xx + erin(1,1) * CIN - 1.3*erout(1,1) * A;
ER2xy = ER2xy + erin(1,2) * CIN - 1.3*erout(1,2) * A;
ER2yx = ER2yx + erin(2,1) * CIN - 1.3*erout(2,1) * A;
ER2yy = ER2yy + erin(2,2) * CIN - 1.3*erout(2,2) * A;

% PARSE TO 1x GRID
DEV.ERxx = ER2xx(2:2:Nx2,1:2:Ny2);
DEV.ERxy = ER2xy(1:2:Nx2,2:2:Ny2);
DEV.ERyx = ER2yx(2:2:Nx2,1:2:Ny2);
DEV.ERyy = ER2yy(1:2:Nx2,2:2:Ny2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATE DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALL anisotropicTL.m
RES = [dx dy];
TL = anisotropicTL(RES,DEV,SIG);

% CALCULATE TOTAL FIELD
E = sqrt(abs(TL.Ex).^2 + abs(TL.Ey).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SET SPACING FOR QUIVER
sp = 18;

% SHOW NUMERICAL PARAMETERS ON CONSOLE
disp(['C    = ' num2str(TL.C/1e-12,'%3.5f') ' pF/m']);
disp(['L    = ' num2str(TL.L/1e-09,'%3.5f') ' nH/m']);
disp(['Z0   = ' num2str(TL.Z0) ' Ohms']);
disp(['nEff = ' num2str(TL.nEff)]);

% VISUALIZE POTENTIAL AND FIELDS
imagesc(xa,ya,TL.V');
colormap(hot);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar;
axis equal tight;
title('Electric Potential V','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);

figure('Color','w');
imagesc(xa,ya,E');
caxis([0 0.75]);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar
colormap(hot);
axis equal tight;
title('|E|','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
hold on;

% ADD QUIVER
[Y,X] = meshgrid(ya,xa);
quiver(X(1:sp:Nx,1:sp:Ny),Y(1:sp:Nx,1:sp:Ny),TL.Ex(1:sp:Nx,1:sp:Ny),...
      TL.Ey(1:sp:Nx,1:sp:Ny),'Color','w');
hold off;

% VISUALIZE DIELECTRIC TENSORS
figure('Color','w');
a = subplot(2,2,1);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2xx');
title('$\varepsilon_{xx}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
%caxis([1.5 +1.5*max(ER2xx(:))]);
colorbar;
colormap(hot);

a = subplot(2,2,2);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2xy');
title('$\varepsilon_{xy}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
%caxis([-2*max(ER2xy(:)) +2*max(ER2xy(:))]);
colorbar;
colormap(hot);

a = subplot(2,2,3);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2yx');
title('$\varepsilon_{yx}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
%caxis([-2*max(ER2yx(:)) +2*max(ER2yx(:))]);
colorbar;
colormap(hot);

a = subplot(2,2,4);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2yy');
title('$\varepsilon_{yy}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
%caxis([1.5 +1.5*max(ER2yy(:))]);
colorbar;
colormap(hot);

% PLOT CONDUCTORS
figure('Color','w');
imagesc(xa,ya,SIG.GND');
title('GND','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
colormap(gray);
colorbar;
caxis([0 1]);

figure('Color','w');
imagesc(xa,ya,SIG.SIG1' | SIG.SIG2');
axis equal tight;
title('SIG1 | SIG2','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
colormap(gray);
colorbar;
caxis([0 1]);