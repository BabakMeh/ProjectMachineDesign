%% Project Machine Design 2
% Babak Mehdizadeh & Nima mohammadi
% 995241050

clc;
clear;

% Initial Setup
initialTemp = 60;  % Ti in °C
loadVal = 3000;    % W in N
deltaTemp = input('Guess the temperature difference (°C): ');
revPerSec = 8;     % N
pressureMPa = 0.469;
pressurePa = pressureMPa * 1e6;
diameter = 80;     % mm
lag = input('Input the radial clearance (mm): ');

radius = diameter / 2;
midTempC = initialTemp + deltaTemp / 2;
midTempF = 32 + (midTempC * 1.8);
fprintf('Temperature (Tf): %.2f °F\n', midTempF);

viscInput1 = input(sprintf('Using chart 12.14, enter the absolute viscosity at %.2f°F: ', midTempF));
fprintf('Viscosity: %.3f reyn\n', viscInput1);
viscPaS1 = viscInput1 * 6.89 * 1e-3;

% Sommerfeld Number
Sval = ((radius / lag)^2) * viscPaS1 * revPerSec / pressurePa;
fprintf('Sommerfeld Number: %.6f\n', Sval);

% Coefficients (assuming l/d = 1)
coeffU = 0.349109 + Sval * 6.0094 + (Sval^2) * 0.047467;
deltaC = 1e-6 * pressurePa * coeffU / 0.12;
deltaTempUpdated = (deltaTemp + deltaC) / 2;
fprintf('Updated ΔT: %.3f °C\n', deltaTempUpdated);

isClose = input('Is updated ΔT close enough? (Yes = 1 / No = 0): ');

% Iteration if needed
while isClose == 0
    avgTempC = initialTemp + deltaTempUpdated / 2;
    avgTempF = 32 + (avgTempC * 9 / 5);

    viscInputLoop = input(sprintf('Using chart 12.14, enter absolute viscosity at %.2f°F: ', avgTempF));
    viscPaSLoop = viscInputLoop * 6.89 * 1e-3;

    Sval = ((radius / lag)^2) * viscPaSLoop * revPerSec / pressurePa;
    fprintf('Sommerfeld Number: %.6f\n', Sval);

    coeffU = 0.349109 + Sval * 6.0094 + (Sval^2) * 0.047467;
    deltaC = 1e-6 * pressurePa * coeffU / 0.12;
    fprintf('ΔC: %.5f\n', deltaC);

    deltaTempUpdated = (deltaTempUpdated + deltaC) / 2;
    fprintf('Updated ΔT: %.3f °C\n', deltaTempUpdated);

    isClose = input('Is updated ΔT now close enough? (Yes = 1 / No = 0): ');
end

fprintf('Final ΔT: %.3f °C\n', deltaTempUpdated);

% Final calculations
finalMidTempC = initialTemp + deltaTempUpdated / 2;
finalMidTempF = 32 + finalMidTempC * 9 / 5;
finalTemp = initialTemp + deltaTempUpdated;

finalVisc = input(sprintf('Final viscosity value at %.2f°F: ', finalMidTempF));
viscFinalPaS = finalVisc * 6.89 * 1e-3;
Sval = ((radius / lag)^2) * viscFinalPaS * revPerSec / pressurePa;
fprintf('Definitive Sommerfeld Number: %.6f\n', Sval);

hRatio = input('Enter h0/c value: ');
h0 = hRatio * lag;

rfRatio = input('Enter rf/c value: ');
fFactor = rfRatio * lag / radius;

Qratio = input('Input Q/(rcNl): ');
Qactual = Qratio * radius * lag * revPerSec * diameter;

QsRatio = input('Enter Qs/Q value: ');
Qs = Qactual * QsRatio;

pressureRatio = input('Enter P/Pmax value: ');
Pmax = (pressurePa / pressureRatio) * 1e-6;

% Power loss calculation
powerLoss = fFactor * loadVal * diameter * pi * revPerSec / 1000;
%% Results
fprintf('Final Sommerfeld Number (S) : %.6f\n', Sval);
fprintf('Mid-film Temperature (Tf)   : %.2f °C\n', finalMidTempC);
fprintf('Outlet Temperature (To)     : %.2f °C\n', finalTemp);
fprintf('Minimum Film Thickness (h₀) : %.4f m\n', h0);
fprintf('Radial Factor (f)           : %.4f\n', fFactor);
fprintf('Heat Loss (H_loss)          : %.2f W\n', powerLoss);
fprintf('Volumetric Flow Rate (Q)    : %.2f mm³/s\n', Qactual);
fprintf('Side Flow (Qs)              : %.2f mm³/s\n', Qs);
fprintf('Maximum Pressure (Pₘₐₓ)      : %.4f MPa\n', Pmax);
