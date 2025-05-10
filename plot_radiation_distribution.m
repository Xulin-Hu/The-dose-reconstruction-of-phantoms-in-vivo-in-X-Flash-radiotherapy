clear all;
close all
clc;

% Energy-deposition distribution
A=load('F:\��ʿ�����о�\��������-������ά������ģ\X����������ֵ����\MC Simulation\The phantom simulation based on the Geant4\2D������ƽ��\Bone_PMMA material\Bone_PMMA material_2.5mm\certain_plane.txt');
A=A*100; % Geant4����ʱ��������������С��100��
% imagesc(A)
surfc(A)

%% �����ܶȡ����١������ݡ����������ϵ����������ɭ��������
Water_density = 1000; % ˮ���ܶ� 1000 kg /m3
Water_sound = 1500;    % ˮ�Ĵ�������Ϊ 1500 m /s
Water_thermal_expansion_coefficient = 210 * 1e-6; % ˮ�����������ϵ��Ϊ 210 * 1e-6
Water_specific_heat_capacity = 4181; % ˮ�ı����� 4181 J/( kg��K)

Bone_density = 1450; % �ܶ� kg /m3
Bone_sound = 4080; % ��������Ϊ m /s
Bone_thermal_expansion_coefficient = 50 * 1e-6;% ���������ϵ��K-1
Bone_specific_heat_capacity = 1760; % ������  J/(kg��K)

PMMA_density = 1190; % �ܶ� kg/m3
PMMA_sound = 2500; % ��������Ϊ m/s
PMMA_thermal_expansion_coefficient = 75 * 1e-6;% ���������ϵ��K-1
PMMA_specific_heat_capacity = 1500; % ������  J/( kg��K)

Water_Gruneisen_coefficient = Water_sound*Water_sound*Water_thermal_expansion_coefficient/Water_specific_heat_capacity; % ������ɭ����
Bone_Gruneisen_coefficient = Bone_sound*Bone_sound*Bone_thermal_expansion_coefficient/Bone_specific_heat_capacity; % ������ɭ����
PMMA_Gruneisen_coefficient = PMMA_sound*PMMA_sound*PMMA_thermal_expansion_coefficient/PMMA_specific_heat_capacity; % ������ɭ����
%%

%% Dose distribution
% ����Բ�����������
m = 120; n = 120; % ����Χ
x0 = 60.5; y0 = 60.5; % ���ĵ����꣬16��(��60,60��ֻ��15������)
r1 = 20; % Pixel����
r2 = 8;
Density_matrix=ones(m, n)*Water_density; % ����һ��ά��Ϊ120*120��ȫ1���󣬱�ʾˮ���ܶ�
[x, y] = meshgrid(1:n, 1:m); 
circleMask = (x - x0).^2 + (y - y0).^2 < r1^2; % logical����
Density_matrix(circleMask) = Bone_density; % �ڻ������ܶ�

ringMask = (x - x0).^2 + (y - y0).^2 <= r1^2 & (x - x0).^2 + (y - y0).^2 >=r2^2;
Density_matrix(ringMask) = PMMA_density; % ��Բ�������ܶ�

Dose =(A*1.6021892*1e-13)./(Density_matrix*0.0025*0.0025*0.0025);% ��������ֲ����(Gy) 
%% Dose-XA signal relation
Gruneisen_coefficient = ones(120, 120)*Water_Gruneisen_coefficient; % ����һ��ά��Ϊ120*120��ȫ1���󣬸�����ɭ����
Gruneisen_coefficient(circleMask) = Bone_Gruneisen_coefficient;  % 
Gruneisen_coefficient(ringMask) = PMMA_Gruneisen_coefficient; % 
Acoustic_pressure = Gruneisen_coefficient.* Density_matrix.* Dose; % ��ѹ�ֲ����
%% ��ͼ
figure % Energy-deposition distribution
imagesc(A);
colormap(parula);  % parula��hot
colorbar;  % ɫ��
set(gca,'xtick',0:10:120)
set(gca,'ytick',0:10:120)
set(gca,'XTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); % FontSize=25
set(gca,'YTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); %

xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','MeV'); % ��λ
% A(A > 4*1e6)= 0; % ������Ĳ���
figure % 3D Dose-deposition distribution
surfc(A)
xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10),zlabel('Z (Dose deposition)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','MeV/Pulse'); % ��λ

figure % 2D Dose-deposition distribution
imagesc(Dose);
colormap(parula);  % parula��hot
colorbar;  % ɫ��
set(gca,'xtick',0:10:120)
set(gca,'ytick',0:10:120)
set(gca,'XTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); % FontSize=25
set(gca,'YTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); %

xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','Gy/Pulse'); % ��λ

figure % 3D Dose-deposition distribution
surfc(Dose)
xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10),zlabel('Z (Dose deposition)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','Gy/Pulse'); % ��λ

figure % 2D Pression distribution
imagesc(Acoustic_pressure);
colormap(parula);  % parula��hot
colorbar;  %ɫ��
set(gca,'xtick',0:10:120)
set(gca,'ytick',0:10:120)

set(gca,'XTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); % FontSize=25
set(gca,'YTickLabel',{'120','10','20','30','40','50','60','70','80','90','100','110'},'FontSize',10); %25

xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10);
% caxis([0 max(max(C))]) 
h=colorbar;
set(get(h,'Title'),'string','Pa/Pulse'); % ��λ

figure % 3D Pression distribution
surfc(Acoustic_pressure)
xlabel('X (Grid)','FontSize',10), ylabel('Y (Grid)','FontSize',10),zlabel('Z (XA signal pressure)','FontSize',10);
h=colorbar;
set(get(h,'Title'),'string','Pa/Pulse'); % ��λ

save('Acoustic_pressure.mat','Acoustic_pressure')