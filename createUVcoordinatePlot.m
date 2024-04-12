% Define the vector in spherical coordinates
azimuth = 45; % degrees
elevation = 30; % degrees
radius = 1; % length of the vector

% Convert degrees to radians for MATLAB functions
azimuth_rad = deg2rad(azimuth);
elevation_rad = deg2rad(elevation);

% Convert spherical to Cartesian coordinates
[x, y, z] = sph2cart(azimuth_rad, elevation_rad, radius);

% Create a new figure
figure;

% Plot the vector
plot3([0 x], [0 y], [0 z], 'b', 'LineWidth', 2); % Plot vector from origin to point
hold on;

% Mark the origin
plot3(0, 0, 0, 'ko'); % 'ko' plots a black circle at the origin

% Draw explicit X, Y, Z axes
axis_limit = 1.5; % Set limit for how far the axes should extend
plot3([-axis_limit axis_limit], [0 0], [0 0], 'k--', 'LineWidth', 1); % X-axis
plot3([0 0], [-axis_limit axis_limit], [0 0], 'k--', 'LineWidth', 1); % Y-axis
plot3([0 0], [0 0], [-axis_limit axis_limit], 'k--', 'LineWidth', 1); % Z-axis

% Enhance plot appearance
grid on;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Vector Plot with Azimuth and Elevation');

% Annotate the azimuth and elevation
text(x, y, z, sprintf('Azimuth: %d°\nElevation: %d°', azimuth, elevation));

% Label the ends of the axes
% text(axis_limit, 0, 0, 'X', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% text(0, axis_limit, 0, 'Y', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% text(0, 0, axis_limit, 'Z', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Set view to 3D for better visualization
view(3);

