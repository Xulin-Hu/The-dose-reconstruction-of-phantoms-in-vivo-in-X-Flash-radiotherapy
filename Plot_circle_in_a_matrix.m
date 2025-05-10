% ������Բ����������
m = 120; n = 120; 
x0 = 60; y0 = 60; 
r1 = 20; 
r2 = 8;
matrix = ones(m, n); 
[x, y] = meshgrid(1:n, 1:m); 
circleMask1 = (x - x0).^2 + (y - y0).^2 < r1^2; % logical����
matrix(circleMask1) = 20;  % ������������Բ����Ĥ�����������ֵ��Ϊ `20`������λ�ñ���Ϊ `1`

ringMask = (x - x0).^2 + (y - y0).^2 < r1^2 & (x - x0).^2 + (y - y0).^2 >=r2^2;
matrix(ringMask) = 10;

imagesc(matrix);
% surfc(matrix)
axis equal; colormap(gray);