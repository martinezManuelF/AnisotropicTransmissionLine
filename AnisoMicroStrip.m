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

% OPEN FIGURE WINDOW
figure('Color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TRANSMISSION LINE PARAMETERS
ersup = [3 3.4641 0;     % Permitivitty tensor of Superstrate
         3.4641 7 0;
         0 0 1];
ersub = [3 3.4641 0;
         3.4641 7 0; % Permittivity tensor of Substrate
         0 0 1];
h = 3*millimeters; % Height of substrate
w = 4*millimeters; % Width of TL

% MISC. PARAMETERS
phi = -30*degrees;;        % Angle of tensor rotation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rot = [cos(phi) -sin(phi) 0;
       sin(phi) cos(phi) 0;
       0  0  1];
       
ersup = Rot*ersup*inv(Rot);
ersub = Rot*ersub*inv(Rot);


% GRID PARAMETERS
BUFF = 3*w;
Sx = BUFF + w + BUFF;
Sy = h + BUFF;
Nx = 200;
Ny = 150;

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
ER2xx = ersup(1,1)*ones(Nx2,Ny2);       % Fill with superstrate
ER2xy = ersup(1,2)*ones(Nx2,Ny2);
ER2yx = ersup(2,1)*ones(Nx2,Ny2);
ER2yy = ersup(2,2)*ones(Nx2,Ny2);
ER2xx(:,Ny2-h/dy2-1:Ny2) = ersub(1,1);  % Fill with substrate
ER2xy(:,Ny2-h/dy2-1:Ny2) = ersub(1,2);
ER2yx(:,Ny2-h/dy2-1:Ny2) = ersub(2,1);
ER2yy(:,Ny2-h/dy2-1:Ny2) = ersub(2,2);

% PARSE TO 1x GRID
DEV.ERxx = ER2xx(2:2:Nx2,1:2:Ny2);
DEV.ERxy = ER2xy(1:2:Nx2,2:2:Ny2);
DEV.ERyx = ER2yx(2:2:Nx2,1:2:Ny2);
DEV.ERyy = ER2yy(1:2:Nx2,2:2:Ny2); 

TL = anisotropicTL([dx dy],DEV,SIG);

% Calculate total field
E = sqrt(abs(TL.Ex).^2 + abs(TL.Ey).^2);

imagesc(xa,ya,TL.V');
colormap(jet);
colorbar;

figure('Color','w');
imagesc(xa,ya,E');
colorbar
caxis([0 0.2]);
colormap(jet);
