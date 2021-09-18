function [ab] = coefficient(xData, yData, zData, k)
%求每一列直线的方向向量

%试试第一列
x = xData(:,k)';
y = yData(:,k)';
z = zData(:,k)';

%去除有效数据
x_ = [];
j = 1;
for i = 1:length(x)
    if(~isnan(x(i))) 
        x_(j) = x(i);
        j = j+1;
    end
end

y_ = [];
j = 1;
for i = 1:length(y)
    if(~isnan(y(i)))
    y_(j) = y(i);
    j = j+1;
    end
end

z_ = [];
j = 1;
for i = 1:length(z)
    if(~isnan(z(i)))
    z_(j) = z(i);
    j = j+1;
    end
end

z1_min = min(z_);
z2_max = max(z_);

%求解直线方程
L = length(x_);
zm = [z_;ones(1,L)];
mm = zm*zm';
xm = x_*zm';
ym = y_*zm';
a = (xm/mm);
b = (ym/mm);

% A = a(1);
% B = b(1);

%求椎角
% z1 = z1_min-5:z2_max+5;
% z1=0:200;
% x1=a(1)*z1+a(2);%空间直线子方程1：x=a*z+b
% y1=b(1)*z1+b(2);%空间直线子方程2：y=c*z+d
% plot3(x1,y1,z1,'r',x_,y_,z_,'.');

%求方向向量
ab = [a(1) a(2) b(1) b(2)];