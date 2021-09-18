%读取测试数据
clear
clc
close all

[Filename, Pathname]=uigetfile({'*.mat'},'请选择文件');  

load([Pathname,Filename]);

figure(1)
for i =1:length(data)
    plot([1:length(data{i})-1],[data{i}{2:end}],'.');
    hold on
end

figure(2)
for i =1:length(data)
    scatter3(X{i},Y{i},Z{i},'.','b');
    axis equal
    hold on
end

figure(3)
for i =1:length(data)
    scatter([1:length(data{i})-1],[data{i}{2:end}],'o','b','filled');
    hold on
end

xlabel('i');
ylabel('R');

figure(4)
for i =1:length(data)
    scatter3(X{i},Y{i},Z{i},'o','b','filled');
    hold on
end
xlabel('x');
ylabel('y');
zlabel('z');