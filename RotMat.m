function R = RotMat(PHI,AXIS)
% RotMax(THETA,AXIS) : Tensor rotation matrices
%
% INPUT ARGUMENTS
% ==============================================================================
% * PHI         Angle of rotation.
% * AXIS        Axis of rotation
%   - x         Rotate along x
%   - y         Rotate along y
%   - z         Rotate along z
%
% OUTPUT ARGUMENTS
% ==============================================================================
% * R           3-by-3 rotation matrix

% Handle axis cases
switch AXIS
  case 'x'
    R = [ 1 0 0 ; 0 cos(PHI) -sin(PHI) ; 0 sin(PHI) cos(PHI) ]; 
  case 'y'
    R = [ cos(PHI) 0 sin(PHI) ; 0 1 0 ; -sin(PHI) 0 cos(PHI) ];
  case 'z'
    R = [ cos(PHI) -sin(PHI) 0 ; sin(PHI) cos(PHI) 0 ; 0 0 1 ];
end

end  