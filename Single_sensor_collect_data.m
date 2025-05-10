% Recording Particle Velocity Example
%
% This example demonstrates how to record the particle velocity using a
% Cartesian or binary sensor mask. It builds on the Homogeneous Propagation
% Medium and Heterogeneous Propagation Medium examples.  
%
% author: Bradley Treeby
% date: 1st November 2010
% last update: 3rd May 2017
% 
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby

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
close all
% =========================================================================
% SIMULATION
% =========================================================================
load('Acoustic_pressure.mat')
% imagesc(Acoustic_pressure)
% create the computational grid
PML_size = 10;              % size of the PML in grid points
Nx = 140 - 2 * PML_size;    % number of grid points in the x direction
Ny = 140 - 2 * PML_size;    % number of grid points in the y direction
dx = 2.5e-3;        % grid point spacing in the x direction [m]
dy = 2.5e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium——水
% medium.sound_speed = 1500;  % [m/s]
% medium.density = 1000;      % [kg/m^3]
% medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
% medium.alpha_power = 1.5; 
%% define the properties of the propagation medium

% 材料密度、声速、比热容、体积热膨胀系数、格律乃森参数定义
Water_density = 1000; % 水的密度 1000 kg /m3
Water_sound = 1500;    % 水的传播声速为 1500 m /s
Water_thermal_expansion_coefficient = 210 * 1e-6; % 水的体积热膨胀系数为 210 * 1e-6
Water_specific_heat_capacity = 4181; % 水的比热容 4181 J/( kg·K)

% Bone
Bone_density = 1450; % 密度 kg /m3
Bone_sound = 4080; % 传播声速为 m /s
Bone_thermal_expansion_coefficient = 50 * 1e-6;% 体积热膨胀系数K-1
Bone_specific_heat_capacity = 1760; % 比热容  J/(kg·K)

PMMA_density = 1190; % 密度 kg/m3
PMMA_sound = 2500; % 传播声速为 m/s
PMMA_thermal_expansion_coefficient = 75 * 1e-6;% 体积热膨胀系数K-1
PMMA_specific_heat_capacity = 1500; % 比热容  J/( kg·K)
% medium.density参数定义
m = 120; n = 120; % 网格范围
x0 = 60.5; y0 = 60.5; % 中心点坐标，16个(（60,60）只有15个网格)
r1 = 20; % Pixel换算
r2 = 8;
medium.density =ones(m, n)*Water_density; % 定义一个维度为120*120的全1矩阵，表示水的密度
[x, y] = meshgrid(1:n, 1:m); 
circleMask = (x - x0).^2 + (y - y0).^2 < r1^2; % logical变量
medium.density(circleMask) = Bone_density; % 内环材料密度

ringMask = (x - x0).^2 + (y - y0).^2 < r1^2 & (x - x0).^2 + (y - y0).^2 >=r2^2;
medium.density(ringMask) = PMMA_density; % 外圆环材料密度

% medium.sound_speed参数定义
medium.sound_speed =ones(m, n)*Water_sound; % 定义一个维度为120*120的全1矩阵，表示水的声速
[x, y] = meshgrid(1:n, 1:m);
medium.sound_speed(circleMask) = Bone_sound; % 内环材料声速
medium.sound_speed(ringMask) = PMMA_sound; % 外圆环材料声速
%% 
% set the input arguements: force the PML to be outside the computational
% grid; switch on p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size,'Smooth', false, 'PlotPML', true,'RecordMovie', true, 'MovieName', 'single_sensor_measurement','PlotLayout', true};
% input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false,'PlotLayout', true};
% create time array
t_end = 300e-6;       % [s]  % 声波时间设置——根据声波情况来
kgrid.makeTime(medium.sound_speed, [], t_end);

% create initial pressure distribution using makeDisc

source.p0 = Acoustic_pressure;

% define four sensor points centered about source.p0
sensor_radius = 50; % [grid points]  12.5cm
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2 + sensor_radius, Ny/2) = 1; % 下
% sensor.mask(Nx/2 - sensor_radius, Ny/2) = 1; % 上
% sensor.mask(Nx/2, Ny/2 + sensor_radius) = 1; % 右
% sensor.mask(Nx/2, Ny/2 - sensor_radius) = 1; % 左

% set the acoustic variables that are recorded
sensor.record = {'p', 'u'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, source.p0 + sensor.mask*5); % 5表示sensor的位置, show the sensor
colormap Parula;
% colormap(getColorMap);
ylabel('X-position [mm]');
xlabel('Y-position [mm]');
axis image;

% get a suitable scaling factor for the time array
[t, t_sc, t_prefix] = scaleSI(kgrid.t_array(end));

% set y-axis limits
p_lim = 4;  % 声压坐标
u_lim = 12e-7;  % 速度坐标

% plot the simulated sensor data
figure;

% plot the pressure
    plot(t_sc * kgrid.t_array, sensor_data.p(1, :), 'k-');
    set(gca, 'YLim', [-p_lim, p_lim+1], 'XLim', [0, t_end * t_sc]);
    xlabel(['Time [' t_prefix 's]']);
    ylabel('XA Signal Pressure (Pa)');

% for sensor_num = 1:4
%     
%     plot the pressure
%     subplot(4, 3, 3 * sensor_num - 2);
%     plot(t_sc * kgrid.t_array, sensor_data.p(sensor_num, :), 'k-');
%     set(gca, 'YLim', [-p_lim, p_lim], 'XLim', [0, t_end * t_sc]);
%     xlabel(['Time [' t_prefix 's]']);
%     ylabel('p');    
% 
%     plot the particle velocity ux
%     subplot(4, 3, 3 * sensor_num - 1);
%     plot(t_sc * kgrid.t_array, sensor_data.ux(sensor_num, :), 'k-');
%     set(gca, 'YLim', [-u_lim, u_lim], 'XLim', [0, t_end * t_sc]);
%     xlabel(['Time [' t_prefix 's]']);
%     ylabel('ux'); 
% 
%     plot the particle velocity uz
%     subplot(4, 3, 3 * sensor_num);
%     plot(t_sc * kgrid.t_array, sensor_data.uy(sensor_num, :), 'k-');
%     set(gca, 'YLim', [-u_lim, u_lim], 'XLim', [0, t_end * t_sc]);
%     xlabel(['Time [' t_prefix 's]']);
%     ylabel('uy'); 
%     
% end
% scaleFig(1, 1.5);