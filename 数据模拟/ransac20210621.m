%���յ�������Ϊ[0 0 -1]  ��dir2 =====������
%����10��
%���յ�������[0.030942986,0.0093595507,-1]
%׶�ǵ����Ϊacos([-0.19566616])===��׶��Ϊ 11.2871��
%���ŵ�Ϊ��acos([-0.17274714]) ====��׶��Ϊ 9.9463��

%����100��  
%���ŵ� ׶�����Ϊacos([-0.16980059])  ====> ׶��Ϊ 9.7744��

%����淨ʸ
%��С���˷���ƽ��
clear all;
close all;
clc;
warning off;

pc = pcread('test.pcd')
% pc=pcdownsample(pc,'random',0.3);       %0.3�������� %�Ѿ�����ȥ�ﴦ����Ҫ�ٽ�����
pcshow(pc);

pc_point = pc.Location';                %�õ���������
kdtree = vl_kdtreebuild(pc_point);      %ʹ��vlfeat����kdtree

normE = [];
dir2 = [];
for i=1:length(pc_point)
    
    p_cur = pc_point(:,i);
    [index, distance] = vl_kdtreequery(kdtree, pc_point, p_cur, 'NumNeighbors', 15);    %Ѱ�ҵ�ǰ�������15����
    p_neighbour = pc_point(:,index)';
    p_cent = mean(p_neighbour);     %�õ��ֲ�����ƽ��ֵ�����ڼ��㷨�������Ⱥͷ���
    
    %��С���˹���ƽ��  ����̫׼ȷ��
    X=p_neighbour(:,1);
    Y=p_neighbour(:,2);
    Z=p_neighbour(:,3);
    XX=[X Y ones(length(index),1)];
    YY=Z;
    %�õ�ƽ�淨����
    C = pinv(XX)*YY;  %��׶��Ϊ��9.8517 ��+86.2068��   ����Ϊ10.1252��k=10��ʱ��

    %�ֲ�ƽ��ָ��ֲ����ĵ�����
    dir1 = p_cent-p_cur';
    %�ֲ�ƽ�淨����.
    dir2(i,:)=[C(1) C(2) -1];  %��ȡ����
    dir3(i,:)=[C(1) C(2) -1]/norm([C(1) C(2) -1]);  %��dir2����һ������
    
    %�������������ļн�
    ang = dir1.*dir2 / (sqrt(dir1(1)^2 +dir2(1)^2) + sqrt(dir1(2)^2 +dir2(2)^2)+sqrt(dir1(3)^2 +dir2(3)^2) );
    
    %���ݼн��жϷ�������ȷ��ָ��   ��Ҫʹ��
    flag = acos(ang);
    dis = norm(dir1);
    if flag<0
        dis = -dis;
    end
    
    %������ǰ��ı��淨����
    t = (0:0.1:2)';
    x = p_cur(1) + C(1)*t;
    y = p_cur(2) + C(2)*t;
    z = p_cur(3) + (-1)*t;
    
    normE =[normE;x y z];
    i;
end
pcshowpair(pc,pointCloud(normE));

%�������S��A�������������Ͱ붥��
%%���ݾ���
data = dir3';
%��������
iter = 100; 
 
%%% �������ݵ�
 figure (3);plot3(data(1,:),data(2,:),data(3,:),'o');hold on; % ��ʾ���ݵ�
 number = size(data,2); % �ܵ���
 bestParameter1=0; bestParameter2=0; bestParameter3=0; % ���ƥ��Ĳ���
 sigma = 1;
 pretotal=0;     %�������ģ�͵����ݵĸ���
 
for i=1:iter
 %%% ���ѡ��������
     idx = randperm(number,3); 
     sample = data(:,idx); 
 
     %%%���׶�Ƿ��� z=ax+by+c
     plane = zeros(1,3);
     x = sample(1,:);
     y = sample(2,:);
     z = sample(3,:);
 
     a = ((z(1)-z(2))*(y(1)-y(3)) - (z(1)-z(3))*(y(1)-y(2)))/((x(1)-x(2))*(y(1)-y(3)) - (x(1)-x(3))*(y(1)-y(2)));
     b = ((z(1) - z(3)) - a * (x(1) - x(3)))/(y(1)-y(3));
     c = z(1) - a * x(1) - b * y(1);
     plane = [a b -1 c];
 
     mask=abs(plane*[data; ones(1,size(data,2))]);    %��ÿ�����ݵ����ƽ��ľ���
     total=sum(mask<sigma);              %�������ݾ���ƽ��С��һ����ֵ�����ݵĸ���
 
     if total>pretotal            %�ҵ��������ƽ�������������ƽ��
         pretotal=total;
         bestplane=plane;          %�ҵ���õ����ƽ��
    end  
 end
 %��ʾ���������ϵ�����
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
 
 %%% �������ƥ��ƽ��
 bestParameter1 = bestplane(1);
 bestParameter2 = bestplane(2);
 bestParameter3 = bestplane(4);
 
 