%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILENAME:         AnisoMicroStrip.m
% COURSE:           EE5322--21st Century Electromagnetics
% INSTRUCTOR:       Raymond C. Rumpf
% NAME:             Manuel F. Martinez
% SEMESTER:         Spring 2018
% DUE DATE:         02/06/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZE MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESTORE STATE
clear all;
close all;
clc;

% UNITS
millimeters = 1;
meters      = 1e3*millimeters;
centimeters = 1e2*millimeters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;
degrees     = pi/180;
F           = 1;
H           = 1;

% CONSTANTS
e0 = 8.85418782e-12 * F/meters;
u0 = 1.25663706e-6 * H/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TRANSMISSION LINE PARAMETERS
ersub = [2 0 0 ; 0 9 0 ; 0 0 6];
h = 3; % Height of substrate
w = 4; % Width of TL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID PARAMETERS
BUFF = 3*w;
Sx = BUFF + w + BUFF;
Sy = h + BUFF;
Nx = 512;
Ny = 512;

% INITIAL GUESS AT GRID RESOLUTION
dx = Sx/Nx;
dy = Sy/Ny;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(w/dx);
dx = w/nx;
ny = ceil(h/dy);
dy = h/ny;

% COMPUTE 2x GRID
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE SIGNALS
SIG.V       = [0 1]; % [GND,SIG1,SIG2,...,SIGn] Forced Potentials
SIG.GND     = zeros(Nx,Ny);
SIG.SIG1    = SIG.GND;

% FORCE CONDUCTORS
SIG.GND(:,1)     = 1;
SIG.GND(:,Ny)    = 1;
SIG.GND(1,:)     = 1;
SIG.GND(Nx,:)    = 1;
nx1 = 1 + floor((Nx-nx)/2);
nx2 = nx1 + nx - 1;
ny  = Ny  - 1 - (h/dy);
SIG.SIG1(nx1:nx2,ny) = 1;


% BUILD PERMITTIVITIES IN 2x GRID
ER2xx = zeros(Nx2,Ny2);
ER2xy = ER2xx;
ER2yx = ER2xx;
ER2yy = ER2xx;

for ny = 1 : Ny2
  
  % CALCULATE FUNCTION TO LINEARLY ROTATE TENSOR
  y = (ny)*dx2;
  
  % CALCULATE ANGLE OF ROTATION
  PHI = (y/h)*(10*degrees);
  
  % ROTATION MATRIX
  R = RotMat(PHI,'z');
  
  % ROTATE TENSOR
  er = R*ersub*R^-1;
  
  % FORM PERMITIVITY ARRAYS
  ER2xx(:,ny) = er(1,1) * ones(Nx2,1);
  ER2xy(:,ny) = er(1,2) * ones(Nx2,1);
  ER2xz(:,ny) = er(1,3) * ones(Nx2,1);
  ER2yx(:,ny) = er(2,1) * ones(Nx2,1);
  ER2yy(:,ny) = er(2,2) * ones(Nx2,1);
  ER2yz(:,ny) = er(2,3) * ones(Nx2,1);
  ER2zx(:,ny) = er(3,1) * ones(Nx2,1);
  ER2zy(:,ny) = er(3,2) * ones(Nx2,1);
  ER2zz(:,ny) = er(3,3) * ones(Nx2,1);
  
end

ER2xx = fliplr(ER2xx);
ER2xy = fliplr(ER2xy);
ER2xz = fliplr(ER2xz);
ER2yx = fliplr(ER2yx);
ER2yy = fliplr(ER2yy);
ER2yz = fliplr(ER2yz);
ER2zx = fliplr(ER2zx);
ER2zy = fliplr(ER2zy);
ER2zz = fliplr(ER2zz);

% PARSE TO 1x GRID
DEV.ERxx = ER2xx(2:2:Nx2,1:2:Ny2);
DEV.ERxy = ER2xy(1:2:Nx2,2:2:Ny2);
DEV.ERyx = ER2yx(2:2:Nx2,1:2:Ny2);
DEV.ERyy = ER2yy(1:2:Nx2,2:2:Ny2); 

% CALL anisotropicTL.m
RES =[dx dy];
TL = anisotropicTL(RES,DEV,SIG);

% Calculate total field
E = sqrt(abs(TL.Ex).^2 + abs(TL.Ey).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW NUMERICAL PARAMETERS ON CONSOLE
disp(['C    = ' num2str(TL.C/1e-12,'%3.5f') ' pF/m']);
disp(['L    = ' num2str(TL.L/1e-09,'%3.5f') ' nH/m']);
disp(['Z0   = ' num2str(TL.Z0) ' Ohms']);
disp(['nEff = ' num2str(TL.nEff)]);

% VISUALIZE POTENTIAL
figure('Color','w');
imagesc(xa,ya,TL.V');
colormap(hot);
colorbar;
axis equal tight;
set(gca,'FontSize',12,'FontWeight','bold');
title('Electric Potential V','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);

% PLOT ELECTRIC FIELD
figure('Color','w');
imagesc(xa,ya,E');
caxis([min(E(:)) max(E(:))/10]);
colorbar
colormap(hot);
axis equal tight;
set(gca,'FontSize',12,'FontWeight','bold');
title('|E|','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
hold on;

% ADD QUIVER
[Y,X] = meshgrid(ya,xa);
quiver(X(1:16:Nx,1:16:Ny),Y(1:16:Nx,1:16:Ny),TL.Ex(1:16:Nx,1:16:Ny),...
      TL.Ey(1:16:Nx,1:16:Ny),'Color','w');
hold off;

% VISUALIZE TENSORS
figure('Color','w');
a = subplot(3,3,1);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2xx');
title('$\varepsilon_{xx}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2xx(:)) max(ER2xx(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,2);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2xy');
title('$\varepsilon_{xy}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2xy(:)) max(ER2xy(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,3);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2xz');
title('$\varepsilon_{xz}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2xz(:)) max(ER2xz(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,4);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2yx');
title('$\varepsilon_{yx}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2yx(:)) max(ER2yx(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,5);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2yy');
title('$\varepsilon_{yy}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2yy(:)) max(ER2yy(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,6);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2yz');
title('$\varepsilon_{yz}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2yz(:)) max(ER2yz(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,7);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2zx');
title('$\varepsilon_{zx}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2zx(:)) max(ER2zx(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,8);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2zy');
title('$\varepsilon_{zy}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([min(ER2zy(:)) max(ER2zy(:))]);
colorbar;
colormap(hot);

a = subplot(3,3,9);
set(a,'FontSize',12);
imagesc(xa2,ya2,ER2zz');
title('$\varepsilon_{zz}$','FontSize',14,'Interpreter','LaTex');
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
caxis([5 max(ER2zz(:))+1]);
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
imagesc(xa,ya,SIG.SIG1');
title('SIG1','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
colormap(gray);
colorbar;
caxis([0 1]);


