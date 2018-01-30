function [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc)
% YEEDER Construct Yee Grid Derivative Operators on a 2D Grid
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc);
%
% Note for normalized grid, use this function as follows:
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0);
%
% Input Arguments
% =================
% NGRID         [Nx Ny] grid size
% RES           [dx dy] grid resolution of the 1X grid
% BC            [xbc ybc] boundary conditions
%                   -2: periodic (requires kinc)
%                   0: Dirichlet
% kinc          [kx ky] incident wave vector
%               This argument is only needed for periodic boundaries.

% DECLARE VARIABLES
Nx  = NGRID(1);
Ny  = NGRID(2);
dx  = RES(1);
dy  = RES(2);
xbc = BC(1);
ybc = BC(2);

lamx = Nx * dx;
lamy = Ny * dy;

% DETERMINE MATRIX SIZE
M = Nx * Ny;

% INITIALIZE MATRICES
DEX = sparse(M,M);
DEY = sparse(M,M);
I   = speye(M,M);

% DEX
if(Nx == 1)
    DEX = 1i * kinc(1) * I;
    DHX = DEX;
else
    % PLACE MAIN DIAGONALS
    DEX = spdiags(-ones(M,1),0,DEX);
    DEX = spdiags(+ones(M,1),1,DEX);
    
    % CORRECT OFF-CENTER MISTAKES (DEFAULT TO DIRICHLET)
    for ny = 1 : Ny-1
        m = Nx * (ny - 1) + Nx;
        DEX(m,m+1) = 0;
    end
    
    % ENFORCE PERIODIC BC
    if(xbc == -2)
        for ny = 1 : Ny
            m = Nx * (ny - 1) + Nx;
            DEX(m,m-(Nx-1)) = exp(1i * (kinc(1) * lamx));
        end
    end
    DEX = DEX/dx;
    DHX = -DEX';
end

% DEY
if(Ny == 1)
    DEY = 1i * kinc(2) * I;
    DHY = DEY;
else
    % PLACE MAIN DIAGONALS
    DEY = spdiags(-ones(M,1),0,DEY);
    DEY = spdiags(+ones(M,1),Nx,DEY);
    
    % ENFORCE PERIODIC BC
    if(ybc == -2)
        DEY = spdiags(+exp(1i*(kinc(2)*lamy))*ones(M,1),-M+Nx,DEY);
    end
    
    DEY = DEY/dy;
    DHY = -DEY';
end
end
