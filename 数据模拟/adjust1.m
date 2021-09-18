%����취����һ��������ȱ�ݣ����ǵ������㲻�Ǿ��ȷֲ���Բ׶���ϵ�ʱ��ģ�ͽ���ʧЧ

% %��ʼ����ֵ
% a = abc(1); b = abc(2); c = abc(3);
% x0 = midxyz(1); y0 = midxyz(2); h0 = midxyz(3);
% alpha = 10/180*pi;

function delta_ = adjust1(a,b,c,x0,y0,h0,alp,X3,Y3,Z3)

%�����������Ϊ
xi = X3; yi = Y3; hi = Z3;
%��������ת��Ϊ��������
xi = xi'; yi = yi'; hi = hi';

%����P�������ϣ������������ʽΪ
xp=x0+a*a*(xi-x0)+a*b*(yi-y0)+a*c*(hi-h0);
yp=y0+a*b*(xi-x0)+b*b*(yi-y0)+b*c*(hi-h0);
hp=h0+a*c*(xi-x0)+b*c*(yi-y0)+c*c*(hi-h0);

%����Ĳ�Ϊ
vi=sqrt((xp-xi).^2+(yp-yi).^2+(h0-hi).^2)-sqrt((xp-x0).^2+(yp-y0).^2+(hp-h0).^2)*tan(alp);
p1=sqrt((xp-xi).^2+(yp-yi).^2+(h0-hi).^2);               %�����ƺ������⣬��Ϊÿһ��xp����仯����liͬʱҲ��仯
p2=sqrt((xp-x0).^2+(yp-y0).^2+(hp-h0).^2)/tan(alp);      %�����ƺ�������

% ���ε��������ֵ��
li = p2.*tan(alp).*tan(alp) - p1;

%��ƫ������
dxpda = 2*a*(xi-x0)+b*(yi-y0)+c*(hi-h0);
dypda = b*(xi-x0);
dhpda = c*(xi-x0);

dxpdb = a*(yi-y0);
dypdb = a*(xi-x0)+2*b*(yi-y0)+c*(hi-h0);
dhpdb = c*(yi-y0);

dxpdc = a*(hi-h0);
dypdc = b*(hi-h0);
dhpdc = a*(xi-x0)+b*(yi-y0)+2*c*(hi-h0);

dxpdx0 = 1-a*a;
dypdx0 = -a*b;
dhpdx0 = -a*c;

dxpdy0 = -a*b;
dypdy0 = 1-b*b;
dhpdy0 = -b*c;

dxpdh0 = -a*c;
dypdh0 = -b*c;
dhpdh0 = 1-c*c;

xpi = xp-xi;
ypi = yp-yi;
hpi = hp-hi;
xp0 = xp-x0;
yp0 = yp-y0;
hp0 = hp-h0;

xpip1 = xpi./p1;
ypip1 = ypi./p1;
hpip1 = hpi./p1;
xp0p2 = xp0./p2;
yp0p2 = yp0./p2;
hp0p2 = hp0./p2;

%�ۺ�
dvida = xpip1.*(dxpda)+ypip1.*(dypda)+hpip1.*(dhpda)-(xp0p2.*(dxpda)+yp0p2.*(dypda)+hp0p2.*(dhpda));
dvidb = xpip1.*(dxpdb)+ypip1.*(dypdb)+hpip1.*(dhpdb)-(xp0p2.*(dxpdb)+yp0p2.*(dypdb)+hp0p2.*(dhpdb));
dvidc = xpip1.*(dxpdc)+ypip1.*(dypdc)+hpip1.*(dhpdc)-(xp0p2.*(dxpdc)+yp0p2.*(dypdc)+hp0p2.*(dhpdc));
dvidx0 = xpip1.*(dxpdx0)+ypip1.*(dypdx0)+hpip1.*(dhpdx0)-(xp0p2.*(dxpdx0 - 1)+yp0p2.*(dypdx0)+hp0p2.*(dhpdx0));
dvidy0 = xpip1.*(dxpdy0)+ypip1.*(dypdy0)+hpip1.*(dhpdy0)-(xp0p2.*(dxpdy0)+yp0p2.*(dypdy0 - 1)+hp0p2.*(dhpdy0));
dvidh0 = xpip1.*(dxpdh0)+ypip1.*(dypdh0)+hpip1.*(dhpdh0)-(xp0p2.*(dxpdh0)+yp0p2.*(dypdh0)+hp0p2.*(dhpdh0 - 1));
dvidalp = -p2.*tan(alp).*sec(alp).*sec(alp);

%���Ի�
A = [dvida dvidb dvidc dvidx0 dvidy0 dvidh0 dvidalp];
% delta1 = [a b c x0 y0 h0 alp];
L = li;

%��������
G=[a b c 0 0 0 0];
W = 0;

%�������Է�����
M = [A'*A G'; G 0]; N = [A'*L; W];
% M * delta = N;

%���
% [U_,S_,V_] = svd(M);
% T = S_;
% T(find(S_~=0)) = 1./S_(find(S_~=0));
% svdInvM = V_ * T' * U_';
% delta_svd=svdInvM*N;

delta_svd = inv(M)*N;
delta_ = [delta_svd(1) delta_svd(2) delta_svd(3) delta_svd(4) delta_svd(5) delta_svd(6) delta_svd(7)];

%��ʼ����
% delat = [a b c x0 y0 z0 alp k];
% a = a + delta_svd(1);
% b = b + delta_svd(2);
% c = c + delta_svd(3);
% x0 = x0 + delta_svd(4);
% y0 = y0 + delta_svd(5);
% z0 = z0 + delta_svd(6);
% alp = alp + delta_svd(7);
% 
% a = a/sqrt(a*a+b*b+c*c);
% b = b/sqrt(a*a+b*b+c*c);
% c = c/sqrt(a*a+b*b+c*c);

% ������ֹ��������
% if delta_'*delta < 1e-3
%     %����
% end