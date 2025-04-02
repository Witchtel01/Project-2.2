clear; clc;
d = 2.5 * 10 ^ (-6); %m (diameter)
pp = 2000; %kg/m^3 (PM2.5 Density)


q = 10 ^ -18; %Coulombs
vpm = 20000/4; %v/m (electric field strength)
H = 0.1524 / 2; %m (distance to next collector plate)
e0 = 8.845 * 10 ^ (-12); %F/m (electric constant)
sig =  e0 * vpm; %C/m^2 (Collector charge density)

fallTime = 1000; %s (time to fall)

m1 = (pp * 4 * pi * (d/2) ^ 2) / 3; %kg (mass of one particle)
avgConc = 2.05 * 10 ^ (-8); %kg/m^3 (PM10, PM2.5 average concentration Hong Kong)
c = avgConc / m1; %n/m^3 (Concentration up to H away from collector)

n = 10000;
Dx = H * rand(1, n);
dxi = max(Dx);
q = q * rand(1, n);

vt = zeros(1, n);

dt = 0.01;
count = zeros(1, n);

f = waitbar(0.0, "Calculating...");
total = 0;
while (any(Dx > 0))
    a = (- (q * sig / (2 * e0)) + (q .* q * c / (2 * e0)) .* (2 * Dx - H)) / ((pi / 6) * pp * (d ^ 3));
    vt = vt + a * dt;
    Dx = Dx + vt * dt;
    count(Dx > 0) = count(Dx > 0) + 1;
    if (mod(total, 1000) == 0)
        waitbar(1 - max(Dx) / dxi, f);
    end
    total = total + 1;
end
delete(f);
disp("Avg time to plate:");
disp(mean(count * dt));
disp(mean(count * dt) < fallTime);
