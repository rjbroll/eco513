function plot_band(x, y, lower, upper, plot_title, plot_ylim)

    % Plot IRF/FVD and error bands

    figure('Units', 'normalize', 'Position', [0.1 0.1 0.8 0.8]);
    plot(x, y, '-kx', 'LineWidth', 2);
    hold on;
    plot(x, [lower; upper], 'LineStyle', '--', 'Color', 'k');
    line([min(x) max(x)], [0 0], 'Color', 'k');
    hold off;
    
    xlim([min(x) max(x)]);
    ylim(plot_ylim);
    set(gca, 'FontSize', 14);
    title(plot_title, 'FontSize', 18);

end