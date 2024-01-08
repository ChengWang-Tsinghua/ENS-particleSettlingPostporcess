function [rho, theta, z] = cartToCylind(X1, X2, RM)

    % cartsianToCylind - Convert Cartesian coordinates to cylindrical coordinates.
    %
    %   [r, theta, z] = cartianToCylind(X1, X2) converts the Cartesian coordinates
    %   of point X2 to cylindrical coordinates (r, theta, z) with respect to the origin at X1.

    % Calculate differences in coordinates

    R = eye(3);
    if nargin>2
        R = RM;
    end
    delta_r = X2-X1;
%     delta_r(3)
    
    rotated_delta_r = R*delta_r';
%     rotated_delta_r(3)

    % Calculate cylindrical coordinates
    rho = sqrt(rotated_delta_r(1)^2 + rotated_delta_r(2)^2);
    theta = atan2(rotated_delta_r(2), rotated_delta_r(1));
    z = rotated_delta_r(3);

%     % Convert theta to degrees if needed
%     theta_deg = rad2deg(theta);

%     % Display the results
%     fprintf('Cylindrical Coordinates:\n');
%     fprintf('r = %.4f\n', r);
%     fprintf('theta = %.4f degrees\n', theta_deg);
%     fprintf('z = %.4f\n', z);
end
