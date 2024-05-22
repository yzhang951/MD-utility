N = 20;
total_segments = 2*N+1;

% ---------- property of HEA ---------------
DeltaU = 0.0846; % Peierls barrier (eV)
a0 = 3.220; % lattice parameter (A)
Gamma = 2.5;%4.22;%2.86;%3.7475*0.7; % line tension (eV/A)

std = 0.0327;
% ---------- property of Fe ---------------
% DeltaU = 0.03487; % Peierls barrier (eV)
% a0 = 2.86652; % lattice parameter (A)
% Gamma = 3.52; % line tension (eV/A)

% ---------- property of Ta ---------------
% DeltaU = 0.03641;
% a0 = 3.3013;
% Gamma = 4.22;

% ---------- property of Nb ---------------
% DeltaU = 0.03500;
% a0 = 3.3004;
% Gamma = 4.02;

Lp = a0*sqrt(2/3);
b = a0*sqrt(3)/2;

converter = 1.602e5; % convert MPa to ev/A^3
period = 1;
%tau_xy = 40; % MPa
