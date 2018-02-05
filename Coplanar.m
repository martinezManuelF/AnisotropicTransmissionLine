%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILENAME:         Coplanar.m
% COURSE:           EE5322--21st Century Electromagnetics
% INSTRUCTOR:       Raymond C. Rumpf
% NAME:             Manuel F. Martinez
% SEMESTER:         Spring 2018
% DUE DATE:         02/06/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESTORE MATLAB STATE
clear all;
close all;
clear;

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
w = 2.5;                % Width of trace
s = 0.1;                % Spacing between traces
ersup = 1.0 * eye(3,3); % Superstrate Tensor
ersub = 2.5 * eye(3,3); % Substrate Tensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID PARAMETERS
BUFF  = 3*w;
Sx    = BUFF + s + w + BUFF;
Sy    = 2*BUFF + 1;
Nx    = 512;
Ny    = 512;

% FIRST GUESS AT RESOLUTION
dx = Sx/Nx;
dy = Sy/Ny;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(s/dx);
dx = s/nx;

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
SIG.V     = [0 1];
SIG.GND   = zeros(Nx,Ny);
SIG.SIG1  = SIG.GND;

% FORCE CONDUCTORS
SIG.GND(:,[1 Ny]) = 1;
SIG.GND([1 Nx],:) = 1;
nx1 = round((w + 2*s)/dx);
nx2 = 1 + floor((Nx-nx1)/2);
nx3 = nx2 + nx1 - 1;
ny = 1 + floor(Ny/2);
SIG.GND(1:nx2,ny) = 1;
SIG.GND(nx3:Nx,ny) = 1;
nx4 = w/dx;
nx5 = 1 + floor((Nx-nx4)/2);
nx6 = nx5 + nx4 - 1;
SIG.SIG1(nx5:nx6,ny) = 1;


