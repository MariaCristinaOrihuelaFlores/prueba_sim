% Simulations In Three Dimensions Example
%
% This example provides a simple demonstration of using k-Wave for the
% simulation and detection of the pressure field generated by an initial
% pressure distribution within a three-dimensional heterogeneous
% propagation medium. It builds on the Homogeneous Propagation Medium and
% Heterogeneous Propagation Medium examples.    
%
% author: Bradley Treeby
% date: 1st July 2009
% last update: 13th April 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clearvars;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 64;            % number of grid points in the x direction
Ny = 64;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
dz = 0.1e-3;        % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500 * ones(Nx, Ny, Nz);	% [m/s]
medium.sound_speed(1:Nx/2, :, :) = 1800;        % [m/s]
medium.density = 1000 * ones(Nx, Ny, Nz);       % [kg/m^3]
medium.density(:, Ny/4:end, :) = 1200;          % [kg/m^3]

% create initial pressure distribution using makeBall
ball_magnitude = 10;    % [Pa]
ball_x_pos = 38;        % [grid points]
ball_y_pos = 32;        % [grid points]
ball_z_pos = 32;        % [grid points]
ball_radius = 5;        % [grid points]
ball_1 = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

ball_magnitude = 10;    % [Pa]
ball_x_pos = 20;        % [grid points]
ball_y_pos = 20;        % [grid points]
ball_z_pos = 20;        % [grid points]
ball_radius = 3;        % [grid points]
ball_2 = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

source.p0 = ball_1 + ball_2;

% define a series of Cartesian points to collect the data
x = (-22:2:22) * dx;            % [m]
y = 22 * dy * ones(size(x));    % [m]
z = (-22:2:22) * dz;            % [m]
sensor.mask = [x; y; z];

% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false, ...
    'DataCast', 'single', 'CartInterp', 'nearest'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the simulation layout using voxelplot
voxelPlot(double(source.p0 | cart2grid(kgrid, sensor.mask)));
view([50, 20]);

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;