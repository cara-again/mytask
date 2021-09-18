%ʹ������˲����������ܷ���ĵ�
%��������
% for i = 1 : length(R)
%      R{i} = R{i} + 86.2068;
% end

%����ת��
X = {};
Y = {};
for i=1:length(R)
    for j=1:length(R{1,i})
        X{i}(j) = R{i}(j).*cos(Theta{i}(j));  %ֱ����pol2cartת���ͺ���
        Y{i}(j) = R{i}(j).*sin(Theta{i}(j));
    end
end

% ���ݺϲ�
Data = {};
a = [];
a1 = [];
a2 = [];
a3 = [];

for i = 1:length(data)
    a1 = R{1,i};
    a2 = Theta{1,i};
    a3 = Z{1,i};
    a = [a1; a2; a3];
    Data{1,i} = a;
end

%���Ⱦ���
R1 = {};
for i = 1:length(R)
    for j = 1:length(R{i})-1
        R1{i}(j) = abs(R{i}(j) - R{i}(j+1));
    end
end

Data_ = Data;

Data2 = Data;

% for i = 1:length(Data2)
%     Data2{i}(:,:) = 0;
% end

% ������Χ�ڵĵ㣨ֻ����һ�Σ�
for i = 1:length(Data)
    for j = 1:length(Data{i})-1
        if(abs(Data{i}(1,j)-Data{i}(1,j+1)) <= 0.05)
            Data2{i}(:,j+1) = Data{i}(:,j+1);
        end
    end
end

%���������������ϵ�����ڵĵ�
Data3 = Data2;
% Data3 = Data;
% % for i = 1:length(Data3)
% %     Data3{i}(:,:) = 0;
% % end

% for i = 1:length(Data)
%     for j = 2:length(Data{i})-2
%         if(abs(Data{i}(1,j-1)-Data{i}(1,j)) > 0.05 && abs(Data{i}(1,j)-Data{i}(1,j+1)) > 0.05 && abs(Data{i}(1,j+1)-Data{i}(1,j+2)) > 0.05)
%             Data3{i}(:,j-1) = 0;
%             Data3{i}(:,j) = 0;
%             Data3{i}(:,j+1) = 0;
%         end
%     end
% end

%������0��
Data4 = Data3;
for i = 1:length(Data4)
    Data4{i}(:,all(Data4{i} == 0,1)) = [];%����all
end

%���������ϵ����ݣ���ֵ��
Data5 = Data4;
for i = 1:length(Data4)
    for j = 1:length(Data4{i})
        if Data4{i}(1,j) > (min(Data4{i}(1,:)) + (max(Data4{i}(1,:)) - min(Data4{i}(1,:)))/5 )  %%
            Data5{i}(:,j) = 0;
        end
    end
end

for i = 1:length(Data5)
    Data5{i}(:,all(Data5{i} == 0,1)) = [];%����all
end

%����ֱ������
Data6 = Data5;
for i = 1:size(Data5,2)
    for j = 1:size(Data5{i},2)   
        Data6{i}(1,j) = Data5{i}(1,j).*cos(Data5{i}(2,j));
        Data6{i}(2,j) = Data5{i}(1,j).*sin(Data5{i}(2,j));
        Data6{i}(3,j) = Data5{i}(3,j);
    end
end

%��ͼ
figure(11)
for i =10:size(Data6,2)
    scatter3(Data6{i}(1,:),Data6{i}(2,:),Data6{i}(3,:),'b','o','filled');
    hold on
end

% X1 = {};
% Y1 = {};
% Z1 = {};

% k = 1;
% for i = 8:size(Data6,2)
%     X1{k} = Data6{i}(1,:);
%     Y1{k} = Data6{i}(2,:);
%     Z1{k} = Data6{i}(3,:);
%     k = k+1;
% end
% 
% X2 = cell2mat(X1);
% Y2 = cell2mat(Y1);
% Z2 = cell2mat(Z1);
% 
% 
% 
% figure(1007)
% [Z3] = meshgrid(Z2);
% surf(X2,Y2,Z3)


