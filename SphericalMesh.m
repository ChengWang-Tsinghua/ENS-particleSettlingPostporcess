function [r,theta,rgrid,tgrid,X,Y] = SphericalMesh(radius,nbinr,nbint, theta)
    % Default values if not provided
    
    if nargin < 2
        nbinr = 20;
        nbint = 20;
    end
    
    if nargin < 4
        thetaMin = 0;
        thetaMax = pi;
    else
        thetaMin = theta(1);
        thetaMax = theta(2);
    end
    
    Rmin = radius(1);
    Rmax = radius(2);
    % Define spherical grid parameters
    r = linspace(Rmin, Rmax, nbinr); % radial coordinates
    theta = linspace(thetaMin, thetaMax, nbint); % azimuthal coordinates

    % Create the mesh grid in spherical coordinates
    [rgrid, tgrid] = meshgrid(r, theta);

    % Convert spherical coordinates to Cartesian coordinates
    X = rgrid .* cos(tgrid);
    Y = rgrid .* sin(tgrid);

%     % Plot the 2D mesh with radial and azimuthal grid lines
%     figure;
%     plot(X, Y, '.');
%     hold on;
% 
%     % Connect points along the azimuthal angle to show azimuthal grid lines
%     for i = 1:size(X, 2)
%         plot(X(:, i), Y(:, i), '-k');
%     end
% 
%     % Connect points along the radial direction to show radial grid lines
%     for i = 1:size(X, 1)
%         plot(X(i, :), Y(i, :), '-k');
%     end
% 
%     hold off;
% 
%     title('Spherical Mesh Grid (2D) with Grid Lines');
%     xlabel('X');
%     ylabel('Y');
%     axis equal;
%     grid on;
end
