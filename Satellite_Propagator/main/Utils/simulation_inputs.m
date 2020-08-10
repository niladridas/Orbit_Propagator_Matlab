no_of_repeat = 10;
samples = 500; % SAMPLE SIZE
R = 5; % Measurement noise sigma
radius_limit = 1;
options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
% OT inputs
cost = @(x) distance_matrix(x);
OT_constantshdl = @(x) OT_constants(x);
