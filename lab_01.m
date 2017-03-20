clc;
clear variables;
close all force;

%������� 2
%�������������� �������
A = randn(6,6);   %������� +
B = randn(6,6);    %�������+

%�������������� �������
v = randn(6,1);    %�������+
w = randn(6,1);   %�������+

%������������� ����������
ATemp = A^(-1);
wTemp = w';      

%�������� ���������
calc = (A*B - B*ATemp)*v*wTemp; %�������+

%�������� ���� � ��� �������� � � �
maxAElem = max(A);
maxA = max(maxAElem); %�������+

minAElem = min(A); 
minA = min(minAElem); %�������+

maxBElem = max(B);
maxB = max(maxBElem);%�������+

minBElem = min(B);
minB = min(minBElem);%�������+

meanAElem = mean(A);
meanA = mean(meanAElem); %�������+

meanBElem = mean(B);
meanB = mean(meanBElem); %�������+

%����������� 
ordA = length(A);
ordB = length(B);
ordv = length(v);
ordw = length(w);
ordCalc = length(calc);

%����� �����
% fprintf('����������� � ����� = %2d', ordA );
% fprintf('\n');
% fprintf('A = \n');
% disp(A);
% fprintf('\n');
% fprintf('����������� B ����� = %2d', ordB );
% fprintf('\n');
% fprintf('B = \n');
% disp(B);
% fprintf('\n');
% fprintf('����������� v ����� = %2d ', ordv);
% fprintf('\n');
% fprintf('v = \n');
% disp(v);
% fprintf('\n');
% fprintf('����������� w ����� = %2d',ordw);
% fprintf('\n');
% fprintf('w = \n');
% disp(w);
% fprintf('\n');
% fprintf('����������� calc ����� = %2d',ordCalc);
% fprintf('\n');
% fprintf('calc = \n');
% disp(calc);
% fprintf('\n');
% fprintf('max A =   %2f', maxA);
% fprintf('\n');
% fprintf('min A =  %2f', minA);
% fprintf('\n');
% fprintf('mean A =  %2f', meanA);
% fprintf('\n');
% fprintf('max B =   %2f', maxB);
% fprintf('\n');
% fprintf('min B =  %2f', minB);
% fprintf('\n');
% fprintf('mean B = %2f', meanB);
% fprintf('\n');

%������� 3
%c������ ������
t = (0.0 : 0.001 : 1.0)';
%������� 2 ��
f = 2;
%����������� �������
omega = 2*pi*f;
%������� ����� �
cosWT = cos(omega*t);
y = cosWT;
ySz = size(y);
% z ��� �� ����������� , ��� � y
z = randn(ySz);

%���������� ����� y
yLen = length(y);
%���� ������ ���� - ��
noisePow = 2;
%��������� ����� ��� 
q = wgn(yLen, 1, noisePow);

%������� 4

%����� ���������� ��� randn
y1 = y+z;
%����� ���������� ��� wgn
y2 = y+q;
%����� ���������� ��� awgn
y3 = awgn(y, noisePow);

%�������� �������������� ������������ t,y,y1,y2,y3 � �������
Y = horzcat(t, horzcat(y,horzcat(y1, horzcat(y2, y3))));
%�������� � ���� .mat
% ����� save
save('Yfile.mat','Y');
% ����� dlmwrite
dlmwrite('YfileASCII.csv', Y , 'delimiter',';');

% �������� ,���� �����
%type('YfileASCII.csv');

% ����� ����
%������� ���� �� ������
fileID = fopen('matrix.txt', 'w');

%������� ������ � �������
Ysz = size(Y);
rowsY = Ysz(1);
columnsY = Ysz(2);

%���� �� ������� � �� �������� � ������� � ����
% �� �������
for i = 1 : 1 : rowsY
%     �� ��������
    for j = 1 : 1 : columnsY
%         ������ � ���� � ������� 2f ����������� �������� ������� Y
        fprintf(fileID, '%2f', Y(i,j) );
%         ����� ������������ ����� ;
        fprintf(fileID, ';');
    end   
%     ������� �� ��������� ������
    fprintf(fileID, '\n');
end

%������ �������� - ������� ����
fclose(fileID);

%��� ����������� ������� ��� ���� � �������
% type('matrix.txt');

% ������� 5 - ���������
%������� 6
% ������� � �������������� ������� 6 �� 2 �� ���������� �������� v � w
 barData = [v(1), w(1); 
            v(2), w(2);
            v(3), w(3);
            v(4), w(4);
            v(5), w(5);
            v(6), w(6)];

%  �� ��������� ������� ������� ����������� 
    b = bar(barData);
    
    %��������� �� gca 
%     current axes
    ax = gca;
    
%     ������ �������� ������� ������
    ax.Title.String = 'v & w bar graph '; 
    ax.Title.FontSize = 20; 
    ax.Title.FontName = 'Times New Roman Cyr ';
    ax.Title.Color = 'g'; 
    ax.XLabel.String = 'i-Position';
        ax.YLabel.String = 'i-Value';
    
%     ������� ����� 1� � 2� �������� �����������
    set(b(1),'FaceColor','r');
    set(b(2),'FaceColor','g');

%     ������� 7

figure(2);
clf;
subplot(3,1,1)
% ������� ������
plot(t, y);
ax1 = gca;
subplot(3,1,2)
% ���������� ��� ��� �������������
stairs(t,y);
ax2 = gca;
subplot(3,1,3)
% ������ ��� ��������� (stem - �������)
stem(t, y);
ax3 = gca;

ax1.Title.String = 'Y(t) - Plot '; 
    ax1.Title.FontSize = 20; 
    ax1.Title.FontName = 'Times New Roman Cyr ';
    ax1.Title.Color = 'r '; 
    ax1.XLabel.String = 't';
        ax1.YLabel.String = 'Y';

    
ax2.Title.String = 'Y(t) - Stairs '; 
    ax2.Title.FontSize = 20; 
    ax2.Title.FontName = 'Times New Roman Cyr ';
    ax2.Title.Color = 'g ';     
     ax2.XLabel.String = 't';
        ax2.YLabel.String = 'Y';
    
ax3.Title.String = ' Y(t) - Stem '; 
    ax3.Title.FontSize = 20; 
    ax3.Title.FontName = 'Times New Roman Cyr ';
    ax3.Title.Color = 'b ';   
     ax3.XLabel.String = 't';
        ax3.YLabel.String = 'Y';

%     ������� 8
%������� ��� ��� ���� �� ������� Y
YLength = length(Y);
Y_ = randn(YLength, 1);
Y1 = randn(YLength, 1);
Y2 = randn(YLength, 1);
Y3 = randn(YLength, 1);

for i = 1 : 1 : YLength
    for j = 1 : 1 : 5
        if j == 2
        Y_(i,1) = Y(i,j);
        end
         if j == 3
        Y1(i,1) = Y(i,j);
         end
         if j == 4
        Y2(i,1) = Y(i,j);
         end
         if j == 5
        Y3(i,1) = Y(i,j);
        end
    end
end

%�������� �������
figure(3);
clf;
fx = subplot(3,1,1);
% ������� ������
plot(t, Y1);
fx.XLabel.String = 't';
fx.YLabel.String = 'Y';
fx.Title.String = 'Y1(t) and Y(t)';
hold on;
plot(t, Y_);
hold off;

fx1 = subplot(3,1,2);
% ������� ������
plot(t, Y2);
fx1.XLabel.String = 't';
fx1.YLabel.String = 'Y';
fx1.Title.String = 'Y2(t) and Y(t)';
hold on;
plot(t, Y_);
hold off;

fx2 = subplot(3,1,3);
% ������� ������
plot(t, Y3);
fx2.XLabel.String = 't';
fx2.YLabel.String = 'Y';
fx2.Title.String = 'Y3(t) and Y(t)';
hold on;
plot(t, Y_);

hold off

    % ������� 9
figure(4);
clf;
[X,Y] = meshgrid(-4:0.1:4);
R = sqrt(X.^2 + Y.^2) + X.^2 + Y.^2;
Z = (sin(R))./R;

ax1 = subplot(2,2,1);
mesh(X,Y,Z);
colormap(ax1, winter);
view(ax1,240,25);
ax1.XLabel.String = 'X';
ax1.YLabel.String = 'Y';
ax1.ZLabel.String = 'Z';
ax1.Title.String = 'Sin(R)/R mesh colormap winter';
hold on;

ax2 = subplot(2,2,3);
surf(X,Y,Z);
colormap(ax2,summer);
view(ax2,240,25);
ax2.XLabel.String = 'X';
ax2.YLabel.String = 'Y';
ax2.ZLabel.String = 'Z';
ax2.Title.String = 'Sin(R)/R surf colormap summer';
hold on;

ax3 = subplot(2,2,[2,4]);
surf(X,Y,Z);
colormap(ax3,parula);
shading flat;
c = colorbar('southoutside');
c.Label.String = 'Function values';
view(ax3,240,15);
ax3.XLabel.String = 'X';
ax3.YLabel.String = 'Y';
ax3.ZLabel.String = 'Z';
ax3.Title.String = 'Sin(R)/R surf no grid, flat shading, colormap parula ';

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
clf;
[X,Y] = meshgrid(-4:0.1:4);
F = sqrt(X.^2 + Y.^2) + X.^2 + Y.^2;
C = (sin(F))./F;
surf(X,Y,C);
colormap cool
shading flat;
view(100,50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















