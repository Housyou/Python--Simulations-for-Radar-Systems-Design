% This helper function plots Pd vs SNR using a non-linear scale along the
% y-axis based on the normal distribution.

%   Copyright 2020 The MathWorks, Inc.

function h = helperRadarPdVsSNRPlot(SNRdB, Pd, ylims)
    yaxisScale = @(y) erfinv(2*y - 1); 
    pdScaled = yaxisScale(Pd);
    h = plot(SNRdB, pdScaled, 'LineWidth', 2);

    xlabel('Effective Detectability Factor (dB)');
    ylabel('Probability of Detection');

    yaxisTicks = [0.0001 0.0005 0.001 0.01 0.02 0.05 0.1:0.1:0.9 ...
        0.95 0.98 0.99 0.995 0.999 0.9995 0.9999 0.99999 0.999999];
    set(gca(), 'ytick', yaxisScale(yaxisTicks), 'yticklabel', num2cell(yaxisTicks));

    ylim(yaxisScale([ylims(1) ylims(end)]));
    grid on;
    grid minor;
end