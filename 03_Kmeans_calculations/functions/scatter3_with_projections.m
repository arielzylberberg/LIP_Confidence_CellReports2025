function p = scatter3_with_projections(data, color)

n = size(data,1);
x = data(:,1);
y = data(:,2);
z = data(:,3);


% % Generate random data points for demonstration
% n = 100; % Number of data points
% x = rand(n, 1) * 10; % Random x-coordinates
% y = rand(n, 1) * 10; % Random y-coordinates
% z = rand(n, 1) * 10; % Random z-coordinates

% Define the condition for color assignment (for demonstration)
% condition = randi([1, 3], n, 1);

% Define color map for the conditions
% color_map = [1, 0, 0; 0, 1, 0; 0, 0, 1]; % Red, Green, Blue

% Initialize figure
p = publish_plot(1,1);
hold on;
% view(116.7000, 33.8158);
view(131.2127, 39.3530);
xlim([-2,2]);
ylim([-2,2]);
zlim([-2,2]);

if 1
    % Plot spheres for 3D data points
    radius = 0.06; % Define the radius of spheres
    for i = 1:n
        [xs, ys, zs] = sphere;
        xs = xs * radius + x(i);
        ys = ys * radius + y(i);
        zs = zs * radius + z(i);
        surf(xs, ys, zs, 'FaceColor', color(i,:), 'EdgeColor', 'none');
    end
else
    scatter3(x, y, z, 60, color,'filled');
    hold all
end
xli = xlim;
yli = ylim;
zli = zlim;
% 2D projections with the same color as the 3D points
scatter3(x, y, zli(1)*ones(size(z)), 20, color,'MarkerFaceColor','w');
scatter3(xli(1)*ones(size(x)), y, z, 20, color,'MarkerFaceColor','w');
scatter3(x, yli(1)*ones(size(y)), z, 20, color,'MarkerFaceColor','w');

% Labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
% title('3D Scatter Plot with 2D Projections');
grid on;

hold off;


