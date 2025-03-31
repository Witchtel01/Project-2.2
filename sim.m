clear; clc;
pf = 1.225; %kg/m^3 (fluid density)
d = 2.5 * 10 ^ (-6); %m (diameter)
pp = 2000; %kg/m^3 (PM2.5 Density)
g = 9.81; %m/s^2 (accel due to gravity)
mu = 1.81 * 10 ^ (-5); %kg/(ms) (dynamic viscosity of air)
q = 1.602 * 10 ^ -15; %Coulombs
sig = 10 ^ -6; %C/m^2 (Collector charge density)
H = 0.1524; %m (distance to next collector plate)
e0 = 8.845 * 10 ^ (-12); %F/m (electric constant)

fallTime = 4 / 0.000444; %m/s (terminal particle velocity)

ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = [0; 0; 1];

m1 = (pp * 4 * pi * (d/2) ^ 2) / 3; %kg (mass of one particle)
avgConc = 2.05 * 10 ^ (-8); %kg/m^3 (PM10, PM2.5 average concentration Hong Kong)
c = avgConc / m1; %n/m^3 (Concentration up to H away from collector)


rng("default");
n = 5;
Dx = H * rand(1, n);

vt = zeros(3, n);

y = zeros(1, n);
z = zeros(1, n);

D = [Dx; y; z];

Cd = @() (24 * mu) ./ (pf * d * vecnorm(vt));

cdVal = Cd();
cdVal(isinf(cdVal) | isnan(cdVal)) = 0;

a = @() ((pi/6) * (pf-pp) * g * (d ^ 3) * khat + ...
    (1/4) * pf * cdVal * ((pi/4) * d^2) .* (vt .^ 2) - ...
    ((q * sig) / (2 * e0)) * ihat + ...
    (((q * q * c) / (2 * e0)) * (2 * D(1, :) - H)) .* ihat) / ...
    ((pi / 6) * pp * (d ^ 3));

steps = 10000;
ending = 1;
dt = 0.001;
for i = 0:dt:ending
    output = a();
    vt = vt + output .* dt;
    D = D + vt .* dt;
    cdVal = Cd();
    cdVal(isinf(cdVal) | isnan(cdVal)) = 0;
end