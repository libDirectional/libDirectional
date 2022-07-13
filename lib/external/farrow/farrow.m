function a = farrow(x0, y0, z0, x1, y1, z1, color, width)
        
    detail = 50;                                % Number of cone lines
    arrowcolor = [.8 .2 .2];                    % Arrow color
    arrowwidth = 3;                             % Arrow width
    
    if nargin > 6                               % If no 7th argument:
        arrowcolor = color;                     % -- default color
    end
    
    if nargin > 7                               % If no 8th argument:
        arrowwidth = width;                     % -- default width
        
        if width < .1                           % If width too small              
            arrowwidth = .01;                   % -- set to minimum width
            detail = 500;                       % -- set high detail
        end
    end
                                                % Start buffer data 
    X = [x0, x1, nan];                          % (tail --> head)
    Y = [y0, y1, nan];                          % 'nan' separates lines
    Z = [z0, z1, nan];
    
    u = [x1 - x0; y1 - y0; z1 - z0];            % Directional vector
    uN = norm(u);                               % Length
    
    for i = 1:detail                            % Loop: Cone buffer
        random_vec = ones(3,1)-2*rand(3,1);     % - Generate random vector
        rad = cross(u, random_vec);             % - Cross product
        rad = rad/norm(rad)*uN/30;              % - Normalize radius
        X = [X, x1, x1 - (x1 - x0)/7 + rad(1), nan]; % - Draw cone-line
        Y = [Y, y1, y1 - (y1 - y0)/7 + rad(2), nan];
        Z = [Z, z1, z1 - (z1 - z0)/7 + rad(3), nan];
    end
    
    a = line(X, Y, Z, ...                       % Output 'farrow'
        'linewidth', arrowwidth, ...
        'color', arrowcolor); 
end
