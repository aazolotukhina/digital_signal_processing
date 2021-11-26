close all;
clear;
clc;
%% Моделирование регистрации спектральной характеристики 

% Моделирование спектра из функций Гаусса
lam=[580, 670, 710, 760, 830]; 
sko=[40.5, 4.8, 10, 10.1, 20.1]; 
A=[0.3, 0.8,0.5, 0.2, 1.0];
lambda=500:880;
t=gaussfunction(lambda, lam, sko, A);
Spectr=sum(t); % формирование спектра объекта из отдельных гауссоид

PSF=gaussfunction(530:670, 600, 20, 0.46); %создание аппаратной функции спектрометра

% Создание матрицы оператора интегрального выражения, описывающего работы спектрометра, из аппаратной функции спектрометра

for i=1:381
    OperatorFull(i,i:i+141-1)=PSF;  
end

Operator=OperatorFull(:,71:end-70);
Spectrum=Spectr';
Signal=Operator*Spectrum; % моделирование сигнала, получаемого в результате работы
% спектрометра

%Введение погрешности в правую часть. Помехи. Снижение ОСШ
PogrTeor=0.01000; % 1% от максимального значения неискаженного сигнала
ErrorSignal = (randi([-10 10], size(Signal)))/10;
Signal=Signal+ErrorSignal*(max(Signal)*PogrTeor);

%% Решение обратной задачи
mu=cond(Operator); % определение числа обусловленности матрица оператора

% МНК и Метод Гаусса
L_Gauss=Operator\Signal;

% МНК и Метод обратной матрицы
Inverse_matrix=inv(Operator);
L_Inverse=Inverse_matrix*Signal; 

% LU разложение
[L,UU]=lu(Operator);
k=L\Signal;
L_LU=UU\k;

% %QR разложение
[Q,R] = qr(Operator);
L_QR=R\(Q'*Signal);

% Метод псевдообратной матрицы
PsevdoInverse_matrix=pinv(Operator);
L_PsevdoInverse=PsevdoInverse_matrix*Signal; 
 
% figure(2)
% plot(lambda, PsevdoInverse_matrix); 

% Метод Тихонова c SVD
[U,s,V] = svd(Operator); 
s=diag(s);
regul=0.56;
L_TikhonovSVD=tikhonovSVD(U,s,V, Signal, regul);

% figure(2);
% plot(lambda, f);

% Метод Тихонова
L_Tikhonov =((regul*eye(381)+Operator'*Operator)\(Operator'*Signal))';

% Метод псевдообратной матрицы с усеченным SVD
u=10;
L_TSVD=V(:,1:u)*inv(diag(s(1:u)))*U(:,1:u)'*Signal;

% Аппроксимация 
p_deg=15;
v8=polyfit(lambda, Signal/20, p_deg); 
L_deg8=polyval(v8,lambda); 

figure (1)

plot(lambda, Spectrum,'--k', 'Linewidth', 2); hold on;
plot(lambda, Signal/20,'k'); hold on;
% plot(lambda, L_LU, '--m'); hold on;
% plot(lambda, L_Gauss, '--g'); hold on;
% plot(lambda, L_QR, '--c'); hold on;
% plot(lambda, L_Inverse, '--y'); hold on;
plot(lambda, L_PsevdoInverse, 'r'); hold on; 
plot(lambda, L_TikhonovSVD,'-.g');
plot(lambda, L_Tikhonov,'--b'); 
plot(lambda, L_TSVD,'m');
plot(lambda, L_deg8,'c');
axis([500 880 0 1]);
legend ({'Spectrum','Signal','PsevdoInverse','TikhonovSVD','Tikhonov', 'TSVD', 'deg8'});

%% Оценка результатов
MethodName={'PsevdoInverse'; 'TikhonovSVD'; 'Tikhonov'; 'TSVD'; 'deg8'};
Results=[L_PsevdoInverse'; L_TikhonovSVD'; L_Tikhonov; L_TSVD'; L_deg8;];

for q=1:5
    scoreTeor= corrcoef(Results(q,:),Spectrum');
    Corrcoef(q,:)= scoreTeor(1,2);
end
Error=table(MethodName, Corrcoef);