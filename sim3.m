clear; clc;
% pf = 1.225; %kg/m^3 (fluid density)
d = 2.5 * 10 ^ (-6); %m (diameter)
pp = 2000; %kg/m^3 (PM2.5 Density)
% g = 9.81; %m/s^2 (accel due to gravity)
% mu = 1.81 * 10 ^ (-5); %kg/(ms) (dynamic viscosity of air)


q = 10 ^ -11; %Coulombs
vpm = 10 ^ 5; %v/m (electric field strength)
H = 0.1524 / 2; %m (distance to next collector plate)
e0 = 8.845 * 10 ^ (-12); %F/m (electric constant)
sig =  e0 * vpm; %C/m^2 (Collector charge density)

fallTime = 4 / 0.000444; %s (time to fall)

% ihat = [1; 0; 0];
% jhat = [0; 1; 0];
% khat = [0; 0; 1];

m1 = (pp * 4 * pi * (d/2) ^ 2) / 3; %kg (mass of one particle)
avgConc = 2.05 * 10 ^ (-8); %kg/m^3 (PM10, PM2.5 average concentration Hong Kong)
c = avgConc / m1; %n/m^3 (Concentration up to H away from collector)


% rng("default");
n = 1000;
% Dx = zeros(1, n);
Dx = H * rand(1, n);
q = q * rand(1, n);

vt = zeros(1, n);

dt = 0.01;
count = zeros(1, n);
% totalCount = 0;
while (any(Dx > 0))
    a = - (q * sig / (2 * e0)) + (q .* q * c / (2 * e0)) .* (2 * Dx - H);
    vt = vt + a * dt;
    Dx = Dx + vt * dt;
    count(Dx > 0) = count(Dx > 0) + 1;
    % totalCount = totalCount + 1;
end
% disp(Dx);
disp("Avg time to plate:");
disp(mean(count * dt));
