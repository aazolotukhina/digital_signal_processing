close all;
clear;
clc;
%%                              Условия
% Дано 
s=-200;
r=50;
n1=1; n2=1.75;
b1=[0; -0.1];
%%                          Матричный метод
tStart1=tic; 
% Матрица переноса 
D0=[1 s; 0 1]; % создаем матрицу 2х2, где в первой строке два элемента - 1 s, 
% и во второй строке также два элемента - 0 1 

% Матрица преломения
F=(1/r)*(1-n1/n2); 
N=n1/n2; % делим n1 на n2
R=[1 0; F N];

b2=R*D0*b1;
s2=b2(1)/b2(2); % делим первый элемент вектора b2 на второй элемент вектора b2
tEnd1= toc(tStart1);
%%                        Построение графика
% Задаем уравнения в аналитическом виде
syms z k
y=sqrt(2*r*z-z.^2); % уравнение сферы 
y1=-b1(2)*(z-s); % уравнение падающего луча
y2=-b2(2)*z+b2(1); %уравнение преломленного луча

figure(1)
subplot(2,1,1)
fplot(y, [0 6],'Linewidth',2.5); hold on;
fplot(y1, [s 0],'Linewidth',1.5);
fplot(y2, [0 s2],'Linewidth',1.5); hold off;

axis([-200 200 0 30]); 
title({'Расчет хода луча матричным методом'}); % название графика
legend({'y', 'y_1', 'y_2'});
time_text=strcat('Время вычислений = ', num2str(tEnd1), ' c');
text(-175, 20, time_text); % вывод строки, определенной на предыдущей строке

% Помещаем ось Y в ноль
ax=gca;
ax.YAxisLocation = 'origin';

%% Расчет хода реального луча с использованием закона преломления
tStart2=tic; 
% %Определение координат точки падения 
equation=y==y1;
zp=double(vpasolve(equation,z,[3 5]));

y_zp_analytical=subs(y1,z,zp); % подставляем в уравнение y1 значение zp, 
y_zp=double(y_zp_analytical); 

% переводим результат из аналитического вида  в числовое значение
y_p=double(subs(y,z,zp));
sigma1=atan(b1(2));

% Угол наклона нормали 
f=atan(1/(diff(y,z)));
f_zp=double(subs(f,z,zp));

% Уравнение нормали 
yn=-tan(subs(f,z,zp))*(z-zp)+y_p;

% Углы падения и преломления
e1=f_zp-sigma1;  
e2=-asin(n1/n2*sin(e1));
sigma2=f_zp+e2;

% Уравнение преломленного луча
y3=-tan(sigma2)*(z-zp)+y_p; 

% Задний отрезок
s3=zp+(y_p/tan(sigma2)); % добавляем в название цифру 3 к букве s

% Продольная аберрация
aberration=s3 - s2;
tEnd2=toc(tStart2);

subplot(2,1,2)
fplot(y, [0 zp],'Linewidth', 2.5); hold on;
fplot(y1, [s zp],'Linewidth', 1.5);
fplot(y3, [zp s3],'Linewidth', 1.5);
fplot(yn, [zp s3], '-.m','Linewidth',1.5);

axis([-200 200 0 30]); 
title('Расчет хода реального луча с использованием закона преломления'); % название графика
legend('y', 'y_1', 'y_2');
time_text=strcat('Время вычислений = ', num2str(tEnd2), ' c');
text(-175, 20, time_text);

% Помещаем ось Y в ноль
ax = gca;
ax.YAxisLocation = 'origin';