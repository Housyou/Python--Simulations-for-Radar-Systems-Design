function helperTargetRCSReturnHistogramPlot(rcsval,mu_rcs)
% This function helperTargetRCSReturnHistogramPlot is only in support of
% TargetRCSExample. It may change in a future release.

% Copyright 2015 The MathWorks, Inc.

clf
histogram(rcsval,50,'Normalization','pdf')
hold on;
rcs_plot = linspace(0,max(rcsval),50);
pdf_plot = 1/mu_rcs*exp(-rcs_plot/mu_rcs);
env_plot = plot(rcs_plot,pdf_plot);
xlabel('RCS (m^2)');
ylabel('Probability Density');
title('Histogram of Swerling 1 Target Return');
legend(env_plot,'Theory');
