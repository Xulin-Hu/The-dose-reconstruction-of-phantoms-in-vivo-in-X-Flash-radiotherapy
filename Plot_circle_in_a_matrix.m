% 矩阵中圆环参数设置
m = 120; n = 120; 
x0 = 60; y0 = 60; 
r1 = 20; 
r2 = 8;
matrix = ones(m, n); 
[x, y] = meshgrid(1:n, 1:m); 
circleMask1 = (x - x0).^2 + (y - y0).^2 < r1^2; % logical变量
matrix(circleMask1) = 20;  % 将矩阵中满足圆形掩膜条件的区域的值设为 `20`，其他位置保持为 `1`

ringMask = (x - x0).^2 + (y - y0).^2 < r1^2 & (x - x0).^2 + (y - y0).^2 >=r2^2;
matrix(ringMask) = 10;

imagesc(matrix);
% surfc(matrix)
axis equal; colormap(gray);