% 1901
% +86.2068 数据前处理环节
data_ = [];
k = 1;
for i = 1 : length(data)
  for j = 1 : length(data{i})
      data_(k,1) = j;
      data_(k,2) = cell2mat(data{i}(j));
      k = k + 1;
  end
end
data_(:,2) = data_(:,2) + 86.2068;
data_1 = mapminmax(data_');
data_1 = data_1';

% distmat = zeros(10790,10790);   %先限制下矩阵范围
for i=1:size(data_1,1)   %这里时间有点长啊（果然是n^2的时间复杂度啊）
    for j=i:size(data_1,1)
        %distmat(i,j)=abs(data(i,2)-data(j,2) );
        distmat(i,j)=sqrt((data_1(i,1:2)-data_1(j,1:2)) * (data_1(i,1:2)-data_1(j,1:2))');
    end
end

for i=1:size(data_1,1)
    for j=i:size(data_1,1)
        distmat(j,i)=distmat(i,j);
    end
end

Eps=0.03;
MinPts=4;
Clust = DBSCAN(distmat,Eps,MinPts);

x=data_1(:,1);
y=data_1(:,2);
%plot different clusters
figure
for i=1:max(Clust)
    hold on
    scatter(x(Clust==i),y(Clust==i),'filled');
end

%plot the noise
hold on
scatter(x(Clust==0),y(Clust==0),'*')