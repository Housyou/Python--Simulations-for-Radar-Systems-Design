% This helper function creates a waterfall plot for the components of the
% detectability factor

%   Copyright 2020 The MathWorks, Inc.
        

function helperDetectabilityWaterfallPlot(data, names)
    BarWidth = 0.8;
    BarLabelFontSize = 8;
%     XTickLabelFontSize = 9;
    RedColor = [1 0.75 0.75];
    GreenColor = [0.75 1 0.75];
    BlueColor = [0, 0.4470, 0.7410];
    axes = newplot;

    ypos = zeros(numel(data) + 1, 2);
    k = 0;
    for i = 1 : numel(data)
        ypos(i, 1) = k;
        k = k + data(i);
        ypos(i, 2) = k;
    end

    data = [data sum(data)];
    ypos(end, 2) = data(end);

    strvalues = arrayfun(@(s)num2str(s, '%.2f'), data, 'UniformOutput', false);
    
    x = [-BarWidth/2 BarWidth/2 BarWidth/2 -BarWidth/2];

    numBars = size(ypos, 1);

    for i = 1 : numBars - 1               
        if ypos(i, 1) > ypos(i, 2)
            c = GreenColor;
        else
            c = RedColor;
        end

        patch(axes, x + i, [ypos(i, 1) ypos(i, 1) ypos(i, 2) ypos(i, 2)], c);

        if abs(ypos(i, 1) - ypos(i, 2)) > 0.5
            text(axes, i, (ypos(i, 1) + ypos(i, 2))/2, strvalues{i},...
                'Color', [0 0 0], 'HorizontalAlignment', 'center',...
                'FontSize', BarLabelFontSize, 'FontWeight', 'bold');
        else
            text(axes, i, max(ypos(i, 1), ypos(i, 2)), strvalues{i},...
                'Color', [0 0 0], 'HorizontalAlignment', 'center',...
                'VerticalAlignment', 'bottom', 'FontSize', BarLabelFontSize, 'FontWeight', 'bold');
        end
    end

    for i = 1 : numBars - 2
         line(axes, [i + BarWidth/2 i+1 - BarWidth/2], [ypos(i, 2) ypos(i+1, 1)], 'LineStyle', '-', 'Color', 'k');
    end

    h = line(axes, [0 numBars], [ypos(end, 2) ypos(end, 2)], 'LineStyle', '--', 'Color', BlueColor, 'LineWidth', 2);
    uistack(h, 'bottom');

    text(axes, numBars, ypos(end, 2), strvalues{end}, 'Color', BlueColor, 'HorizontalAlignment', 'right',...
        'VerticalAlignment', 'bottom', 'FontSize', BarLabelFontSize, 'FontWeight', 'bold');

    axes = gca;
%     axes.XAxis.FontSize = XTickLabelFontSize;

    axes.XTick = 1 : size(ypos, 1)-1;
    axes.XTickLabel = names;
    ylabel('SNR (dB)');
    title('Detectability Factor');
    grid on;
end