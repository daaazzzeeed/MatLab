clc;
clear variables;
close all force;

%задание 2
%инициализируем матрицы
A = randn(6,6);   %вывести +
B = randn(6,6);    %вывести+

%инициализируем векторы
v = randn(6,1);    %вывести+
w = randn(6,1);   %вывести+

%промежуточные вычисления
ATemp = A^(-1);
wTemp = w';      

%вычислим выражение
calc = (A*B - B*ATemp)*v*wTemp; %вывести+

%вычислим макс и мин значения А и В
maxAElem = max(A);
maxA = max(maxAElem); %вывести+

minAElem = min(A); 
minA = min(minAElem); %вывести+

maxBElem = max(B);
maxB = max(maxBElem);%вывести+

minBElem = min(B);
minB = min(minBElem);%вывести+

meanAElem = mean(A);
meanA = mean(meanAElem); %вывести+

meanBElem = mean(B);
meanB = mean(meanBElem); %вывести+

%размерности 
ordA = length(A);
ordB = length(B);
ordv = length(v);
ordw = length(w);
ordCalc = length(calc);

%вывод всего
% fprintf('Размерность А равна = %2d', ordA );
% fprintf('\n');
% fprintf('A = \n');
% disp(A);
% fprintf('\n');
% fprintf('Размерность B равна = %2d', ordB );
% fprintf('\n');
% fprintf('B = \n');
% disp(B);
% fprintf('\n');
% fprintf('Размерность v равна = %2d ', ordv);
% fprintf('\n');
% fprintf('v = \n');
% disp(v);
% fprintf('\n');
% fprintf('Размерность w равна = %2d',ordw);
% fprintf('\n');
% fprintf('w = \n');
% disp(w);
% fprintf('\n');
% fprintf('Размерность calc равна = %2d',ordCalc);
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

%задание 3
%cоздаем сигнал
t = (0.0 : 0.001 : 1.0)';
%частота 2 Гц
f = 2;
%циклическая частота
omega = 2*pi*f;
%косинус омега т
cosWT = cos(omega*t);
y = cosWT;
ySz = size(y);
% z той же размерности , что и y
z = randn(ySz);

%количество строк y
yLen = length(y);
%сила белого шума - дб
noisePow = 2;
%формируем белый шум 
q = wgn(yLen, 1, noisePow);

%задание 4

%белый аддитивный шум randn
y1 = y+z;
%белый аддитивный шум wgn
y2 = y+q;
%белый аддитивный шум awgn
y3 = awgn(y, noisePow);

%проведем горизонтальную конкатенацию t,y,y1,y2,y3 в матрицу
Y = horzcat(t, horzcat(y,horzcat(y1, horzcat(y2, y3))));
%сохраним в файл .mat
% через save
save('Yfile.mat','Y');
% через dlmwrite
dlmwrite('YfileASCII.csv', Y , 'delimiter',';');

% просмотр ,если нужен
%type('YfileASCII.csv');

% через цикл
%откроем файл на запись
fileID = fopen('matrix.txt', 'w');

%объявим строки и столбцы
Ysz = size(Y);
rowsY = Ysz(1);
columnsY = Ysz(2);

%цикл по строкам и по столбцам с записью в файл
% по строкам
for i = 1 : 1 : rowsY
%     по столбцам
    for j = 1 : 1 : columnsY
%         запись в файл в формате 2f конкретного элемента матрицы Y
        fprintf(fileID, '%2f', Y(i,j) );
%         пусть разделителем будет ;
        fprintf(fileID, ';');
    end   
%     прыгаем на следующую строку
    fprintf(fileID, '\n');
end

%запись окончена - закроем файл
fclose(fileID);

%для уверенности выведем наш файл в консоль
% type('matrix.txt');

% задание 5 - прочитано
%задание 6
% объявим и инициализируем матрицу 6 на 2 со значениями векторов v и w
 barData = [v(1), w(1); 
            v(2), w(2);
            v(3), w(3);
            v(4), w(4);
            v(5), w(5);
            v(6), w(6)];

%  по значениям матрицы создаем гистограмму 
    b = bar(barData);
    
    %указатель на gca 
%     current axes
    ax = gca;
    
%     всякие красивые свойств тайтла
    ax.Title.String = 'v & w bar graph '; 
    ax.Title.FontSize = 20; 
    ax.Title.FontName = 'Times New Roman Cyr ';
    ax.Title.Color = 'g'; 
    ax.XLabel.String = 'i-Position';
        ax.YLabel.String = 'i-Value';
    
%     зададим цвета 1х и 2х столбцов гистограммы
    set(b(1),'FaceColor','r');
    set(b(2),'FaceColor','g');

%     задание 7

figure(2);
clf;
subplot(3,1,1)
% обычный график
plot(t, y);
ax1 = gca;
subplot(3,1,2)
% ступенечки как без антиалиазинга
stairs(t,y);
ax2 = gca;
subplot(3,1,3)
% внутри все закрашено (stem - стебель)
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

%     задание 8
%возьмем все что надо из матрицы Y
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

%построим графики
figure(3);
clf;
fx = subplot(3,1,1);
% обычный график
plot(t, Y1);
fx.XLabel.String = 't';
fx.YLabel.String = 'Y';
fx.Title.String = 'Y1(t) and Y(t)';
hold on;
plot(t, Y_);
hold off;

fx1 = subplot(3,1,2);
% обычный график
plot(t, Y2);
fx1.XLabel.String = 't';
fx1.YLabel.String = 'Y';
fx1.Title.String = 'Y2(t) and Y(t)';
hold on;
plot(t, Y_);
hold off;

fx2 = subplot(3,1,3);
% обычный график
plot(t, Y3);
fx2.XLabel.String = 't';
fx2.YLabel.String = 'Y';
fx2.Title.String = 'Y3(t) and Y(t)';
hold on;
plot(t, Y_);

hold off

    % задание 9
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
















