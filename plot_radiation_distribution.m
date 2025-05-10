clear all;
close all
clc;

% Energy-deposition distribution
A=load('F:\博士课题研究\热声成像-人体三维剂量建模\X光声成像数值仿真\MC Simulation\The phantom simulation based on the Geant4\2D中心切平面\Bone_PMMA material\Bone_PMMA material_2.5mm\certain_plane.txt');
A=A*100; % Geant4仿真时粒子数不够，减小了100倍
% imagesc(A)
surfc(A)

%% 材料密度、声速、比热容、体积热膨胀系数、格律乃森参数定义
Water_density = 1000; % 水的密度 1000 kg /m3
Water_sound = 1500;    % 水的传播声速为 1500 m /s
Water_thermal_expansion_coefficient = 210 * 1e-6; % 水的体积热膨胀系数为 210 * 1e-6
Water_specific_heat_capacity = 4181; % 水的比热容 4181 J/( kg·K)

Bone_density = 1450; % 密度 kg /m3
Bone_sound = 4080; % 传播声速为 m /s
Bone_thermal_expansion_coefficient = 50 * 1e-6;% 体积热膨胀系数K-1
Bone_specific_heat_capacity = 1760; % 比热容  J/(kg·K)

PMMA_density = 1190; % 密度 kg/m3
PMMA_sound = 2500; % 传播声速为 m/s
PMMA_thermal_expansion_coefficient = 75 * 1e-6;% 体积热膨胀系数K-1
PMMA_specific_heat_capacity = 1500; % 比热容  J/( kg·K)

Water_Gruneisen_coefficient = Water_sound*Water_sound*Water_thermal_expansion_coefficient/Water_specific_heat_capacity; % 格律乃森常数
Bone_Gruneisen_coefficient = Bone_sound*Bone_sound*Bone_thermal_expansion_coefficient/Bone_specific_heat_capacity; % 格律乃森常数
PMMA_Gruneisen_coefficient = PMMA_sound*PMMA_sound*PMMA_thermal_expansion_coefficient/PMMA_specific_heat_capacity; % 格律乃森常数
%%

%% Dose distribution
% 环形圆柱体参数定义
m = 120; n = 120; % 网格范围
x0 = 60.5; y0 = 60.5; % 中心点坐标，16个(（60,60）只有15个网格)
r1 = 20; % Pixel换算
r2 = 8;
Density_matrix=ones(m, n)*Water_density; % 定义一个维度为120*120的全1矩阵，表示水的密度
[x, y] = meshgrid(1:n, 1:m); 
circleMask = (x - x0).^2 + (y - y0).^2 < r1^2; % logical变量
Density_matrix(circleMask) = Bone_density; % 内环材料密度

ringMask = (x - x0).^2 + (y - y0).^2 <= r1^2 & (x - x0).^2 + (y - y0).^2 >=r2^2;
Density_matrix(ringMask) = PMMA_density; % 外圆环材料密度

Dose =(A*1.6021892*1e-13)./(Density_matrix*0.0025*0.0025*0.0025);% 区域剂量分布情况(Gy) 
%% Dose-XA signal relation
Gruneisen_coefficient = ones(120, 120)*Water_Gruneisen_coefficient; % 定义一个维度为120*120的全1矩阵，格律乃森矩阵
Gruneisen_coefficient(circleMask) = Bone_Gruneisen_coefficient;  % 
Gruneisen_coefficient(ringMask) = PMMA_Gruneisen_coefficient; % 
Acoustic_pressure = Gruneisen_coefficient.* Density_matrix.* Dose; % 声压分布情况
%% 绘图
figure % Energy-deposition distribution
imagesc(A);
colormap(parula);  % parula、hot
colorbar;  % 色阶
set(gca,'xtick',0:10:120)
set(gca,'ytick',0:10:120)
set(gca,'XTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); % FontSize=25
set(gca,'YTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); %

xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','MeV'); % 单位
% A(A > 4*1e6)= 0; % 检查中心部分
figure % 3D Dose-deposition distribution
surfc(A)
xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10),zlabel('Z (Dose deposition)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','MeV/Pulse'); % 单位

figure % 2D Dose-deposition distribution
imagesc(Dose);
colormap(parula);  % parula、hot
colorbar;  % 色阶
set(gca,'xtick',0:10:120)
set(gca,'ytick',0:10:120)
set(gca,'XTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); % FontSize=25
set(gca,'YTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); %

xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','Gy/Pulse'); % 单位

figure % 3D Dose-deposition distribution
surfc(Dose)
xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10),zlabel('Z (Dose deposition)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','Gy/Pulse'); % 单位

figure % 2D Pression distribution
imagesc(Acoustic_pressure);
colormap(parula);  % parula、hot
colorbar;  %色阶
set(gca,'xtick',0:10:120)
set(gca,'ytick',0:10:120)

set(gca,'XTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); % FontSize=25
set(gca,'YTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); %25

xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10);
% caxis([0 max(max(C))]) 
h=colorbar;
set(get(h,'Title'),'string','Pa/Pulse'); % 单位

figure % 3D Pression distribution
surfc(Acoustic_pressure)
xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10),zlabel('Z (XA signal pressure)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','Pa/Pulse'); % 单位

save('Acoustic_pressure.mat','Acoustic_pressure')