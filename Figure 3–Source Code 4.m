%% 已知离散数据点（x,y）坐标，用最小二乘法拟合圆圈半径和圆心坐标，利用生成的角度数组可以画出拟合的圆圈。
function [center_x,center_y,radius] = Fit_circle(x,y)
% 假设已知圆圈的一些离散数据点坐标为 (x, y)
% 线性回归
A = [x(:), y(:), ones(size(x(:)))]; % 构建设计矩阵
b = -(x(:).^2 + y(:).^2); % 构建响应向量
coefficients = (A' * A) \ (A' * b); % 最小二乘法求解拟合系数

% 获取圆心和半径
center_x = -0.5 * coefficients(1);
center_y = -0.5 * coefficients(2);
radius = sqrt(center_x^2 + center_y^2 - coefficients(3));

% 绘制拟合圆圈
theta = linspace(0, 2*pi, 100); % 生成角度数组
circle_x = center_x + radius * cos(theta); % 根据圆的生成坐标数组
circle_y = center_y + radius * sin(theta);
plot(circle_x, circle_y, 'r');
hold on;
scatter(x, y, 'b'); % 原始数据点
end

