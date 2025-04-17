function [fig] = CASPErFigSettings(fig)
%CASPERFIGSETTINGS Formats the figure to be "publication quality" 
%   intent is that everyone in-lab uses this so we have consistent
%   formatting for all of the figures being created. Requires passing
%   the specified figure in as a variable.
%   Andrew Winter August 8, 2023
    
    fig.Color = [1, 1, 1]; % Change figure background color

    % Use handle ("fig") of the figure to access all properties
    axesProperties = findobj(fig, 'type', 'axe'); % also known as "gca"
    xlabelString = get(get(axesProperties(1), 'xlabel'), 'string');
    ylabelString = get(get(axesProperties(1), 'ylabel'), 'string');
    set(axesProperties, 'Fontweight', 'normal', 'FontSize', 14); % Make fonts readable
    
    % Check if axes are labeled:
    if isempty(xlabelString)
        set(get(axesProperties(1), 'xlabel'), 'string', 'LABEL ME!'); % User should see this and correct the labels
    end
    if isempty(ylabelString)
        set(get(axesProperties(1), 'Ylabel'), 'string', 'LABEL ME!'); % User should see this and correct the labels
    end
    
    box on;
    % Maybe make a default figure size

    return;
end

