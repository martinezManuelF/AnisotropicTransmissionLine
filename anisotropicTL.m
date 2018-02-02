function TL = anisotropicTL(RES,DEV,SIG)
% anisotropicTL: Anisotropic transmission line analysis.
%
% INPUT ARGUMENTS
% ==============================================================================
% RES         Grid resolution
%   RES(1)    Resolution along x
%   RES(2)    Resolution along y
%
% DEV         Device structure
%   .ERxx     Permittivity tensor
%   .ERxy     Permittivity tensor
%   .ERyx     Permittivity tensor
%   .ERyy     Permittivity tensor
%
% SIG         Conductor structure
%   .GND      Ground plane matrix
%   .SIG1     Signal 1 matrix
%   .SIGn     nth Signal matrix
%   .V        Array containing forced potential values for each conductor. ...
%             nSIG + GND = lenght(SIG.V)
%
% OUTPUT ARGUMENTS
% ==============================================================================
% TL          Transmission line analysis results
%   .V        Electric Potential
%   .Ex       Electric field intensity along x
%   .Ey       Electric field intensity along y
%   .C        Distributed capacitance of transmission line
%   .L        Distributed inductance of transmission line
%   .Z0       Characteristic impedance of transmission line
%   .nEff     Effective refractive index of transmission line

% Units and constants
mm  = 1;
m   = 1e3 * mm;
e0  = 8.8541878176e-12 * 1/m;
c0  = 299792458 * m;

% Extract Parameters
dx        = RES(1);
dy        = RES(2);
[Nx, Ny]  = size(SIG.GND);
M         = Nx*Ny;
I         = speye(M,M);

% Determine number of signal traces
pot     = SIG.V;
SIG     = rmfield(SIG,'V');
signals = fieldnames(SIG);
nSig    = length(signals);

% Construct derivative matrix operators
NS = [Nx Ny];
BC = [0 0];
[DVX,DVY,DEX,DEY] = yeeder(NS,RES,BC);

% Form interpolation matrix from derivative operators
RXP = (dx/2)*abs(DVX);
RYM = (dy/2)*abs(DEY);
R   = RXP*RYM;

% Diagonalize material tensors
ERxx = diag(sparse(DEV.ERxx(:)));
ERxy = diag(sparse(DEV.ERxy(:)));
ERyx = diag(sparse(DEV.ERyx(:)));
ERyy = diag(sparse(DEV.ERyy(:)));

% Construct composite tensor
ER = [ERxx , R*ERxy ; R'*ERyx , ERyy];

% Build inhomogeneous case
L = [DEX DEY] * ER * [DVX ; DVY];

% Build homogeneous case
Lh = [DEX DEY] * [DVX ; DVY];

% Force known potentials
F   = zeros(Nx,Ny);
vf  = F;
for i = 1 : nSig
  F   = F | SIG.(char(signals(i)));
  vf  = vf + pot(i)*SIG.(char(signals(i)));
end 

F = diag(sparse(F(:)));

L   = (I - F)*L + F;
Lh  = (I - F)*Lh + F;
b   = F*vf(:);

% Compute potentials
v   = L\b;
vh  = Lh\b;

% Compute E Fields
e   = -[DVX ; DVY]*v;
eh  = -[DVX ; DVY]*vh;

% Compute D fields
d   = ER*e;
dh  = eh;

% Compute TL parameters
TL.C  = (e0*dx*dy) * d' * e;    % Distributed capacitance
Ch    = (e0*dx*dy) * dh' * eh;
TL.L  = 1/(c0^2*Ch);            % Distributed inductance
TL.Z0 = sqrt(TL.L/TL.C);        % Characteristic impedance

% Obtain fields
TL.Ex = e(1:M);
TL.Ey = e(M+1:2*M);

% Reshape to 2D Grid
TL.V  = reshape(v,Nx,Ny);
TL.Ex = reshape(TL.Ex,Nx,Ny);
TL.Ey = reshape(TL.Ey,Nx,Ny);
 
% Calculate effective refractive index
TL.nEff = c0*sqrt(TL.L/TL.C);
 
end
