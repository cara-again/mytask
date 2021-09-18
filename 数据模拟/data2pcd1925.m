%将无效值剔除并计算椎角
%其中偏差值为86.2068//86.1226\86.2910

%读取测试数据
clear
clc
close all

[Filename, Pathname]=uigetfile({'*.mat'},'请选择文件');  

load([Pathname,Filename]);

%原始数据不足，需要反推
R = {};
theta = {};
for i = 1 : length(X)
    [theta{i},R{i}] = cart2pol(X{i}, Y{i});
end

%最小值
R1 = R;

%重置空间
M1 = R;
for i=1:length(R)
    M1{i}(1:end) = nan;
end


for i=1:length(R)
    min1 = min(R{i});
    for j = 1:length(R{i})
        if(abs(R{i}(j) - min1)<0.5)
            M1{i}(j) = R{i}(j);
        end
    end
end

% figure(1)
% for i =1:length(M1)
%     plot([1:length(M1{i}-1)],[M1{i}{2:end}],'.');  //要么全是nan，要么全部是
%     hold on
% end


 %上面的部分去掉
 for i=1                               %根据data1811
    for j = 1:length(R{i})
            M1{i}(j) = nan;
    end
 end
 
 for i = 1 : length(M1)
     M1{i} = M1{i} + 86.2068;
 end

%坐标转换
X1 = {};
Y1 = {};
for i=1:length(M1)
    for j=1:length(M1{1,i})
        X1{i}(j) = M1{i}(j).*cos(theta{i}(j));  %直接用pol2cart转换就好了
        Y1{i}(j) = M1{i}(j).*sin(theta{i}(j));
    end
end

Z1 = Z;
for i=1:length(M1)
    for j=1:length(M1{1,i})
        if(isnan(X1{i}(j)) || isnan(Y1{i}(j)))  %isnan
            Z1{i}(j) = NaN;
        end
    end
end

% 取有效数据
xData = [];
yData = [];
zData = [];
for i =1:length(X1)
    for j=2:length(X1{i})
        xData(i,j-1) = X1{i}(j);
        yData(i,j-1) = Y1{i}(j);
        zData(i,j-1) = Z1{i}(j);
    end
end

%取有效值
%这里似乎可以用all(~isnan)
X2 = {};
Y2 = {};
Z2 = {};

for k = 1:length(xData)
    x_ = youxiaozhix(xData, k);
    X2{k} = x_;
    y_ = youxiaozhiy(yData, k);
    Y2{k} = y_;
    z_ = youxiaozhiz(zData, k);
    Z2{k} = z_;
end

%转化为行向量

X2_ = [];
Y2_ = [];
Z2_ = [];

k = 1;
for i = 1:length(X2)
    for j = 1:length(X2{i})
        X2_(k) = X2{i}(j);
        Y2_(k) = Y2{i}(j);
        Z2_(k) = Z2{i}(j);
        k = k+1;
    end
end 

X3 = X2_;
Y3 = Y2_;
Z3 = Z2_;


%转化为pcd文件格式
Data = [X3;Y3;Z3]';%确保为n*3格式

Data = single(Data);
ptCloud = pointCloud(Data(:,1:3));
pcwrite(ptCloud, 'test.pcd', 'Encoding', 'ascii'); %将程序中的xyz数据写入pcd文件中
pc = pcread('test.pcd');

figure(2)
pcshow(pc); %显示点云