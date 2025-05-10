% 2D Time Reversal Reconstruction For A Line Sensor Example
%
% This example demonstrates the use of k-Wave for the time-reversal
% reconstruction of a two-dimensional photoacoustic wave-field recorded
% over a linear array of sensor elements. The sensor data is simulated and
% then time-reversed using kspaceFirstOrder2D. It builds on the 2D FFT 
% Reconstruction For A Line Sensor Example. 
%
% author: Bradley Treeby
% date: 6th July 2009
% last update: 25th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby

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

clearvars; % data delate
% Firstly, run the procedure: plot_radiation_distribution
close all; % close the plots
% =========================================================================
% SIMULATION
% =========================================================================
load('Acoustic_pressure.mat')
% create the computational grid
PML_size = 10;              % size of the PML in grid points
Nx = 140 - 2 * PML_size;    % number of grid points in the x direction
Ny = 140 - 2 * PML_size;    % number of grid points in the y direction
dx = 2.5e-3;                 % grid point spacing in the x direction [m]
dy = 2.5e-3;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% define the properties of the propagation medium,and the propagation medium only need to define the density and speed
% medium.sound_speed = 1500;	% [m/s]

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


% assign to the source structure
source.p0 = Acoustic_pressure;

% define a binary line sensor
%% Sensor define
% Linear array sensors
% sensor.mask = zeros(Nx, Ny);
% sensor.mask(2, 6:115) = 1;  % 传感器数量 216个
% sensor.mask(118, 1:120) = 1;  % 传感器数量 216个

% sensor.mask(6:115, 2) = 1;  % 传感器数量 216个
% sensor.mask(5:110, 118) = 1;  % 传感器数量 216个

% Circle array sensors
% Parameters for the circular array
N = 128;           % Number of sensors
radius = 50;     % Radius of the circular array 50 grids * 0.25 = 12.5 cm
% Create a mask for the sensors
sensor.mask = zeros(Ny, Nx); % Initialize the mask
% Center of the circle (can adjust as needed)
centerX = Nx / 2;
centerY = Ny / 2;
% Calculate sensor positions
for k = 1:N
angle = 2 * pi * (k - 1) / N; % Angle for each sensor
x = round(centerX + radius * cos(angle)); % x-coordinate
y = round(centerY + radius * sin(angle)); % y-coordinate
    
% Set the sensor position in the mask
if x > 0 && x <= Nx && y > 0 && y <= Ny % Check bounds
sensor.mask(y, x) = 1; % Mark sensor position
end
end
%% Simulation process

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false,'RecordMovie', true, 'MovieName', 'example_movie1','PlotLayout', true};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false,'RecordMovie', true, 'MovieName', 'example_movie2'};
% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add first order compensation for only recording over a half plane
% 初始线性传感器仅设置在半平面，因此需要一阶步长，如果全平面，则不需要
% p0_recon = 2 * p0_recon;

% % repeat the FFT reconstruction for comparison
% p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, medium.sound_speed, ...
%     'PosCond', true, 'Interp', '*linear');
% 
% % define a second k-space grid using the dimensions of p_xy
% [Nx_recon, Ny_recon] = size(p_xy);
% kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy);
% 
% % resample p_xy to be the same size as source.p0
% p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure;
% imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, Acoustic_pressure + sensor.mask,[-1, 1]); % normalization
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, Acoustic_pressure + sensor.mask); % normal 
colormap(parula);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% plot the reconstructed initial pressure 
figure;
% imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon,[-1, 1]); % normalization
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon); % normal
colormap(parula);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% apply a positivity condition
p0_recon(p0_recon < 0) = 0;

% plot the reconstructed initial pressure with positivity condition
figure;
% imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon,[-1, 1]); % normalization
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon); % normal
colormap(parula);
colormap parula;
ylabel('X (mm)');
xlabel('Y (mm)');
h=colorbar;
set(get(h, 'Title'), 'string', 'Pa/Pulse')% 单位
axis image;
% scaleFig(1, 0.65);

% plot the initial pressure
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, Acoustic_pressure); % normal
colormap(parula);
colormap parula;
ylabel('X (mm)');
xlabel('Y (mm)');
h=colorbar;
set(get(h, 'Title'), 'string', 'Pa/Pulse')% 单位
axis image;

%% 重建剂量 VS MC仿真
% plot the relative error
figure;
reconstruction_dose_distribution = p0_recon ./ Gruneisen_coefficient ./ Density_matrix;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3,reconstruction_dose_distribution);
% surfc(reconstruction_dose_distribution)；
colormap(parula);  % parula、hot
colorbar;  % 色阶
ylabel('X (mm)');
xlabel('Y (mm)');
h=colorbar;
set(get(h, 'Title'), 'string', 'Gy/Pulse')% 单位
axis image;

figure;
error=(abs(reconstruction_dose_distribution-Dose)./ Dose)*100; % 乘以100%,误差分布图
error_radiation_field = error(51:70,51:70); % 照射野内部相对误差（百分比）
% error_radiation_field = error; % 整体相对误差（百分比）
imagesc(error_radiation_field); % 照射野内部相对误差（百分比），比较这个相对合理
% imagesc(error); % 水箱内部相对误差（百分比）
colormap parula;
colorbar;  % 色阶
ylabel('X (mm)');
xlabel('Y (mm)');
set(gca,'xtick', 1:5:20)
set(gca,'ytick', 1:5:20)
set(gca,'XTickLabel',{'-10','-5','0','5','10'}); % 
set(gca,'YTickLabel',{'-10','-5','0','5','10'}); %

h=colorbar;
t=get(h,'YTickLabel');
t=strcat(t,'%');
set(h,'YTickLabel',t); %单位
axis image;
% scaleFig(1, 0.65);

%% 计算剂量重建误差误差
% average_relative_error = mean(error_radiation_field(:));
% max_relative_error = max(error_radiation_field(:));
Reconstruction_dose = reshape(reconstruction_dose_distribution(51:70,51:70), [], 1); % 转换为列向量
Ground_truth = reshape(Dose(51:70,51:70), [], 1); % 转换为列向量
[mae1,mse1,rmse1,mape1,error1,errorPercent1,mre1]=calc_error(Ground_truth,Reconstruction_dose);