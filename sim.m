clear; clc;
pf = 1.225; %kg/m^3 (fluid density)
d = 2.5 * 10 ^ (-6); %m (diameter)
pp = 2000; %kg/m^3 (PM2.5 Density)
g = 9.81; %m/s^2 (accel due to gravity)
mu = 1.81 * 10 ^ (-5); %kg/(ms) (dynamic viscosity of air)
q = 1; %Coulombs
sig = 1; %Collector charge density
H = 0.1524; %m (distance to next collector plate)
e0 = 8.845 * 10 ^ (-12); %F/m (electric constant)

ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = [0; 0; 1];

m1 = (pp * 4 * pi * (d/2) ^ 2) / 3; %kg (mass of one particle)
aqi = 13 * 10 ^ (-9); %kg/m^3 (AQI Hong Kong)
c = aqi / m1; %n/m^3 (Concentration up to H away from collector)

vt = [0; 0; 0];

Cd = @() (24 * mu) / (pf * d * norm(vt));
a = @() ((pi/6) * (pf-pp) * g * (d ^ 3) * khat + ...
    (1/4) * pf * Cd() * ((pi/4) * d^2) * (vt .^ 2) - ...
    ((q * sig) / (2 * e0)) * ihat + ...
    (((q * q * c) / (2 * e0)) * (2 * Dx - H)) * ihat);

rng("default");
n = 6000;
Dx = H * rand(1, n);
