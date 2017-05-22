function setupAxisCircular(varargin)
    % Changes the settings of a given axis (x,y,z) to visualize a circular
    % quantitiy on [0,2pi)
    %
    % Parameters:
    %   axisName (string)
    %       either 'x' or 'y' or 'z'
    for axisName = varargin
        if strcmp(axisName,'x')
            xlim([0,2*pi]);
            set(gca,'XTick', [0,pi,2*pi]);
            set(gca,'XTickLabel', {'0', '\pi', '2\pi'}); %this does not work in older matlab versions (prior to 2014b?)
        elseif strcmp(axisName,'y')
            ylim([0,2*pi]);
            set(gca,'YTick', [0,pi,2*pi]);
            set(gca,'YTickLabel', {'0', '\pi', '2\pi'}); %this does not work in older matlab versions (prior to 2014b?)
        elseif strcmp(axisName,'z')
            zlim([0,2*pi]);
            set(gca,'ZTick', [0,pi,2*pi]);
            set(gca,'ZTickLabel', {'0', '\pi', '2\pi'}); %this does not work in older matlab versions (prior to 2014b?)
        else
            error('invalid axis')
        end    
    end
end