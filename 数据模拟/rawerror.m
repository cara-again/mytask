% （重要）承接data2pcd.m

% x = [];
% y = [];
% z = [];
% [x,y,z] = pol2cart(Theta{1,20},R{1,20}+86.2068,Z{1,20});
% scatter(x,y,z,'.','b');
% axis([-200 200 -200 200]);

% abc1 = [-0.00867827897544223,-0.00567428356203975,0.999946243545163];
% alp = [10.012912];
% X0 = [-9.5819063;-6.2926908;1101.0614];

%1901
% abc1 = [-0.00498661576121135,-0.00432139212413837,0.999978229379700];
% alp1 = [10.036967];
% X0 = [-5.6275134,-4.7783532,1098.4368];

%1910
% abc1 = [-0.00567535661242441,-0.00456568261927486,0.999973472082906];
% alp1 = [10.034208];
% X0 =[-6.4212599;-5.0545988;1098.7336];
% 
%        a =     0.02591  (0.01464, 0.03717)
%        b =   -0.002124  (-0.005625, 0.001377)
%        h =        1117  (1114, 1120)
%        k =       5.746  (5.73, 5.762)
% abc1 = [0 0 1];
% alp1 = 9.8725;
% X0 = [-0.02591, 0.002124, 1117];
%1910 法矢
% abc1 = [-0.00567535661242441,-0.00456568261927486,0.999973472082906];
% alp1 = [10.034208];
% X0 = [-6.4212599,-5.0545988,1098.7336];

% 1910 
% abc1 = [0 0 1]
% alp1 = 9.8742
% X0 = [-0.02545, 0, 1116]

%1901
X0 = [-5.4282565;-4.7781353;1098.3915]';
abc1 = [-0.00498661576121135;-0.00432139212413837;0.999978229379700]';
alp1 = [10.036967];

%1901
% X0 = [-2.591;2.124;1116]';
% abc1 = [0 0 1]';
% alp1 = [10.036967];


%密封面各点误差求解
vi = [];

% 初始值
% abc1 = [-0.00867827897544223;-0.00567428356203975;0.999946243545163];
%  abc1 = [-0.00238966254262207;-0.00458918673965534;0.999986614349413];
a = abc1(1); b = abc1(2); c = abc1(3);

% X0 = [-9.8746052;-6.3026247;1101.1158];
% X0 = [-9.7708988,-6.2991443,1101.0964];
%  X0 = [-2.4907751,-4.9241958,1097.7065];
x0 = X0(1); y0 = X0(2); h0 = X0(3);
% 
% alp1 = 10.012912;
% alp1 = 10.043347;
alpha = alp1/180*pi;


%模拟误差
%测量点的坐标为
% xi = pc_point(1,:); yi = pc_point(2,:); hi = pc_point(3,:);
xi = xData; yi = yData; hi = zData;
% xi = newData(1,:); yi = newData(2,:); hi = newData(3,:); 
%将行向量转化为列向量：
% xi = xi'; yi = yi'; hi = hi';

%交点P（轴线上）的坐标迭代方式为
xp=x0+a*a*(xi-x0)+a*b*(yi-y0)+a*c*(hi-h0);
yp=y0+a*b*(xi-x0)+b*b*(yi-y0)+b*c*(hi-h0);
hp=h0+a*c*(xi-x0)+b*c*(yi-y0)+c*c*(hi-h0);

%距离的差为
ppi = sqrt((xp-xi).^2+(yp-yi).^2+(hp-hi).^2);            
pp0 = sqrt((xp-x0).^2+(yp-y0).^2+(hp-h0).^2); 
vi = ppi - pp0 * tan(alpha);
    
%作图
% X4 = repmat(X3, length(X3), 1);
% Y4 = repmat(Y3, length(Y3), 1);
% Z4 = repmat(Z3, length(Z3), 1);
% vi = vi'; 

% surf(X3,Y3,vi)

% figure (3)
% plot3(pc_point(1,:),pc_point(2,:),vi);



% colorbar
% scatter3(X3,Y3,vi,'.')

% x4 =X3; 
% y4 = Y3; 
% z4 = -5.6638*sqrt((x4'+2.4907751).^2+(y4+4.9241958).^2)+1097.7065; %矩阵
% surf(x4,y4,z4,vi)

% figure(3)
% for i =1:length(vi)
%     plot(i, vi(i), '.-');
%     hold on
% end

% figure(4)
% for i =1:length(vi)
%     plot(X3(i), vi(i), '.-');
%     hold on
% end

% figure(5)
% for i =1:length(vi)
%     plot(rho(i), vi(i), '.-');
%     hold on
% end

% vi_positive = [];
% vi_negative = [];
% kn = 1;
% for i = 1 : length(vi)
%     if vi(i) < 0
%        vi_negative(kn) = vi(i); 
%        kn = kn + 1;   
%     end 
% end
% 
% kp = 1;
% for j = 1 : length(vi)
%     if vi(j) > 0 || vi(j) == 0
%        vi_positive(kp) = vi(j); 
%        kp = kp + 1;   
%     end   
% end

% analysis_vi = [max(vi);min(vi);mean(vi);std(vi);mean(abs(vi))];
vi2 = vi;  % 备份

% 插值估计要弄个一个新的function,但是首先要计算点数
%对于超出x 范围的xi 的分量，使用方法’nearest’、’linear’、’v5cubic’的插值算法，相应地将返回NaN。对其他的方法，interp1 将对超出的分量执行外插值算法。
%不能直接使用笛卡尔坐标中的来进行插值。会变形。需要先对R进行插值。
% for i = 1 : length(xData)
%     index = find(~isnan(xData(:,i)));
%     len_x(i) = length(index);
% end
% l = max(len_x);
% 
% xf = zeros(l,length(xData));
% yf = zeros(l,length(xData));
% zf = zeros(l,length(xData));
% vif = zeros(l,length(xData));
% for k = 1 : length(xData)
%     %无法执行赋值，因为左侧的大小为 15×1，右侧的大小为 16×1。需要先对矩阵规定内存
%     xf1 = interruptxData(xData,k);
%     yf1 = interruptyData(yData,k);
%     zf1 = interruptzData(zData,k);
%     vif1 = interruptvi2Data(vi2,k);
%     
%     if length(xf1) == l
%         xf(:,k) = xf1;
%         yf(:,k) = yf1;
%         zf(:,k) = zf1;
%         vif(:,k) = vif1;
%     else
%         xf1 = interp1(linspace(1,length(xf1),length(xf1)),xf1, linspace(1,l,l),'ppval');
%         yf1 = interp1(linspace(1,length(yf1),length(yf1)),yf1, linspace(1,l,l),'ppval');
%         zf1 = interp1(linspace(1,length(zf1),length(zf1)),zf1, linspace(1,l,l),'ppval');
%         vif1 = interp1(linspace(1,length(vif1),length(vif1)),vif1, linspace(1,l,l),'ppval');
%         xf(:,k) = xf1;
%         yf(:,k) = yf1;
%         zf(:,k) = zf1;
%         vif(:,k) = vif1;        
%         
%     end
% end

% figure (1)
% surf(xf,yf,zf,vif)
% colorbar
    
figure (1)
surf(xData,yData,zData,vi)
colorbar
    
figure (6)
scatter(X3,Z3,'k','.')
    
    
    
