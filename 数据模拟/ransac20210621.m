%最终的轴向量为[0 0 -1]  以dir2 =====》错误
%迭代10次
%最终的轴向量[0.030942986,0.0093595507,-1]
%锥角的余角为acos([-0.19566616])===》锥角为 11.2871°
%最优的为：acos([-0.17274714]) ====》锥角为 9.9463°

%迭代100次  
%最优的 锥角余角为acos([-0.16980059])  ====> 锥角为 9.7744°

%求表面法矢
%最小二乘法求平面
clear all;
close all;
clc;
warning off;

pc = pcread('test.pcd')
% pc=pcdownsample(pc,'random',0.3);       %0.3倍降采样 %已经经过去燥处理，不要再降噪了
pcshow(pc);

pc_point = pc.Location';                %得到点云数据
kdtree = vl_kdtreebuild(pc_point);      %使用vlfeat建立kdtree

normE = [];
dir2 = [];
for i=1:length(pc_point)
    
    p_cur = pc_point(:,i);
    [index, distance] = vl_kdtreequery(kdtree, pc_point, p_cur, 'NumNeighbors', 15);    %寻找当前点最近的15个点
    p_neighbour = pc_point(:,index)';
    p_cent = mean(p_neighbour);     %得到局部点云平均值，便于计算法向量长度和方向
    
    %最小二乘估计平面  （不太准确）
    X=p_neighbour(:,1);
    Y=p_neighbour(:,2);
    Z=p_neighbour(:,3);
    XX=[X Y ones(length(index),1)];
    YY=Z;
    %得到平面法向量
    C = pinv(XX)*YY;  %半锥角为：9.8517 （+86.2068）   否则为10.1252（k=10的时候）

    %局部平面指向局部质心的向量
    dir1 = p_cent-p_cur';
    %局部平面法向量.
    dir2(i,:)=[C(1) C(2) -1];  %抽取出来
    dir3(i,:)=[C(1) C(2) -1]/norm([C(1) C(2) -1]);  %对dir2做归一化处理
    
    %计算两个向量的夹角
    ang = dir1.*dir2 / (sqrt(dir1(1)^2 +dir2(1)^2) + sqrt(dir1(2)^2 +dir2(2)^2)+sqrt(dir1(3)^2 +dir2(3)^2) );
    
    %根据夹角判断法向量正确的指向   需要使用
    flag = acos(ang);
    dis = norm(dir1);
    if flag<0
        dis = -dis;
    end
    
    %画出当前点的表面法向量
    t = (0:0.1:2)';
    x = p_cur(1) + C(1)*t;
    y = p_cur(2) + C(2)*t;
    z = p_cur(3) + (-1)*t;
    
    normE =[normE;x y z];
    i;
end
pcshowpair(pc,pointCloud(normE));

%构造矩阵S和A，求轴线向量和半顶角
%%数据矩阵
data = dir3';
%迭代次数
iter = 100; 
 
%%% 绘制数据点
 figure (3);plot3(data(1,:),data(2,:),data(3,:),'o');hold on; % 显示数据点
 number = size(data,2); % 总点数
 bestParameter1=0; bestParameter2=0; bestParameter3=0; % 最佳匹配的参数
 sigma = 1;
 pretotal=0;     %符合拟合模型的数据的个数
 
for i=1:iter
 %%% 随机选择三个点
     idx = randperm(number,3); 
     sample = data(:,idx); 
 
     %%%拟合锥角方程 z=ax+by+c
     plane = zeros(1,3);
     x = sample(1,:);
     y = sample(2,:);
     z = sample(3,:);
 
     a = ((z(1)-z(2))*(y(1)-y(3)) - (z(1)-z(3))*(y(1)-y(2)))/((x(1)-x(2))*(y(1)-y(3)) - (x(1)-x(3))*(y(1)-y(2)));
     b = ((z(1) - z(3)) - a * (x(1) - x(3)))/(y(1)-y(3));
     c = z(1) - a * x(1) - b * y(1);
     plane = [a b -1 c];
 
     mask=abs(plane*[data; ones(1,size(data,2))]);    %求每个数据到拟合平面的距离
     total=sum(mask<sigma);              %计算数据距离平面小于一定阈值的数据的个数
 
     if total>pretotal            %找到符合拟合平面数据最多的拟合平面
         pretotal=total;
         bestplane=plane;          %找到最好的拟合平面
    end  
 end
 %显示符合最佳拟合的数据
mask=abs(bestplane*[data; ones(1,size(data,2))])<sigma;    
hold on;
k = 1;
for i=1:length(mask)
    if mask(i)
        inliers(1,k) = data(1,i);
        inliers(2,k) = data(2,i);
        plot3(data(1,i),data(2,i),data(3,i),'r+');
        k = k+1;
    end
end
 
 %%% 绘制最佳匹配平面
 bestParameter1 = bestplane(1);
 bestParameter2 = bestplane(2);
 bestParameter3 = bestplane(4);
 
 