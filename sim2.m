clear;clc;
collectorLength = 4; %m
freeFallVelocity = 0.000444; %m/s
fallTime = collectorLength / freeFallVelocity; %s (time to fall at terminal velocity)
volts = 12; %v (plate voltage difference)
plateDistance = 0.1525; %m
eField = volts / plateDistance; % v/m (electric field between each plate)

fieldForce = q * eField;

