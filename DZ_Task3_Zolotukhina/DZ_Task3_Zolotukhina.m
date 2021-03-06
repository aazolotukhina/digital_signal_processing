clear;
clc;
close all;
% Решение прямой и обратной задачи спектроскопии
%% Моделирование регистрации спектральным прибором (Этап 1)
% Исходные данные:
Procent1 = 0.01;      % 1% от максимального значения неискаженного сигнала
Procent2 = 0.15;      % 15% от максимального значения неискаженного сигнала
load('PSF_3.mat');
load('Spectr_4.mat');

% Моделирование спектра из функций Гаусса
t=gaussfunction(lambda, lam, sko, A);
Spectr=sum(t); % формирование спектра объекта из отдельных гауссоид

% Создание матрицы оператора интегрального выражения, описывающего работу
% спектрометра, из аппаратной функции спектрометра
for i=1:size(lambda,2)
    ind(i) =  i+151-1;
    OperatorFull(i,i:i+151-1)=PSF_3;  
end

Operator=OperatorFull(:,76:end-75);
mu=cond(Operator);    % определение числа обусловленности матрица оператора
fprintf('Число обусловленности матрицы оператора - %d \n',mu)
Spectrum=Spectr';
% моделирование сигнала, получаемого в результате работы спектрометра
Signal=Operator*Spectrum; 

% Cлучаи наличия случайной равномерно распределенной помехи
Error = (randi([-10 10], size(Signal)))/10;  
SignalError1=Signal + Error*(max(Signal)*Procent1);
Error = (randi([-10 10], size(Signal)))/10;
SignalError2=Signal + Error*(max(Signal)*Procent2);
%%  Решение обратной задачи (Этап 2)
% Метод псевдообратной матрицы
PsevdoInverse_matrix=pinv(Operator);
L_start = tic;
L_PsevdoInverse0=PsevdoInverse_matrix*Signal;
L_end_PI0 = toc(L_start);
L_start = tic;
L_PsevdoInverse1=PsevdoInverse_matrix*SignalError1;
L_end_PI1 = toc(L_start);
L_start = tic;
L_PsevdoInverse2=PsevdoInverse_matrix*SignalError2;
L_end_PI2 = toc(L_start);
% QR разложение
[Q,R] = qr(Operator);
L_start = tic;
L_QR0=R\(Q'*Signal);
L_end_QR0 = toc(L_start);
L_start = tic;
L_QR1=R\(Q'*SignalError1);
L_end_QR1 = toc(L_start);
L_start = tic;
L_QR2=R\(Q'*SignalError2);
L_end_QR2 = toc(L_start);
% Метод Тихонова c SVD
[U,s,V] = svd(Operator); 
s=diag(s);
regul=0.5;
L_start = tic;
L_TikhonovSVD0=tikhonovSVD(U,s,V, Signal, regul);
L_end_TikhSVD0 = toc(L_start);
L_start = tic;
L_TikhonovSVD1=tikhonovSVD(U,s,V, SignalError1, regul);
L_end_TikhSVD1 = toc(L_start);
L_start = tic;
L_TikhonovSVD2=tikhonovSVD(U,s,V, SignalError2, regul);
L_end_TikhSVD2 = toc(L_start);
% Метод Тихонова
L_start = tic;
L_Tikhonov0 =((regul*eye(381)+Operator'*Operator)\(Operator'*Signal))';
L_end_Tikh0 = toc(L_start);
L_start = tic;
L_Tikhonov1 =((regul*eye(381)+Operator'*Operator)\(Operator'*SignalError1))';
L_end_Tikh1 = toc(L_start);
L_start = tic;
L_Tikhonov2 =((regul*eye(381)+Operator'*Operator)\(Operator'*SignalError2))';
L_end_Tikh2 = toc(L_start);
% Метод псевдообратной матрицы с усеченным SVD
u=10;
L_start = tic;
L_TSVD0=V(:,1:u)*inv(diag(s(1:u)))*U(:,1:u)'*Signal;
L_end_TSVD0 = toc(L_start);
L_start = tic;
L_TSVD1=V(:,1:u)*inv(diag(s(1:u)))*U(:,1:u)'*SignalError1;
L_end_TSVD1 = toc(L_start);
L_start = tic;
L_TSVD2=V(:,1:u)*inv(diag(s(1:u)))*U(:,1:u)'*SignalError2;
L_end_TSVD2 = toc(L_start);
% Глобальная интерполяция
p_deg=15; 
v0=polyfit(lambda, Signal, p_deg);
L_start = tic;
L_deg0=polyval(v0,lambda);
L_end_deg0 = toc(L_start);
v1=polyfit(lambda, SignalError1, p_deg);
L_start = tic;
L_deg1=polyval(v1,lambda);
L_end_deg1 = toc(L_start);
v2=polyfit(lambda, SignalError2, p_deg);
L_start = tic;
L_deg2=polyval(v2,lambda);
L_end_deg2 = toc(L_start);
%%  Количественная оценка полученных решений (Этап 3)
% для случая без погрешности
MethodName0={'PsevdoInverse';'QR'; 'TikhonovSVD'; 'Tikhonov'; 'TSVD'; 'deg'};
Results0=[L_PsevdoInverse0'; L_QR0'; L_TikhonovSVD0'; L_Tikhonov0; L_TSVD0'; L_deg0];

for q=1:6
    scoreTeor= corrcoef(Results0(q,:),Spectrum');
    Corrcoef0(q,:)= scoreTeor(1,2);
end
% для случая с погрешностью 1%
Results1=[L_PsevdoInverse1'; L_QR1'; L_TikhonovSVD1'; L_Tikhonov1; L_TSVD1'; L_deg1];

for q=1:6
    scoreTeor= corrcoef(Results1(q,:),Spectrum');
    Corrcoef1(q,:)= scoreTeor(1,2);
end
% для случая с погрешностью 15%
Results2=[L_PsevdoInverse2'; L_QR2'; L_TikhonovSVD2'; L_Tikhonov2; L_TSVD2'; L_deg2];

for q=1:6
    scoreTeor= corrcoef(Results2(q,:),Spectrum');
    Corrcoef2(q,:)= scoreTeor(1,2);
end
Error=table(MethodName0, Corrcoef0, Corrcoef1, Corrcoef2,'VariableNames',...
    {'Название метода', 'Коэфф корреляции (сигнал без погрешн)',...
    'Коэфф корреляции (сигнал с погрешн 1%)','Коэфф корреляции (сигнал с погрешн 15%)'});
writetable(Error, 'Error.xlsx'); %запись структуры типа таблица
disp(Error);
%%  Поиск 4 лучших решений (Этап 3)
[Coeff_sorted,Ind] = sort(Corrcoef0);
% Сортировка значений по сумме коэффициентов
Method_Time0 = [L_end_PI0,L_end_QR0,L_end_TikhSVD0,...           % Массив времен выполнения вычислений
L_end_Tikh0,L_end_TSVD0,L_end_deg0];
Methods_sorted = [];
for k = 1:4
    Methods_sorted(k,:) = Results0(Ind(k),:);
    MethodName_sorted{:,k} = MethodName0{Ind(k),1};
    Method_Time_sorted(k,:) = Method_Time0(1,Ind(k)); 
end
best_plotting(Methods_sorted, Method_Time_sorted,MethodName_sorted,lambda,Signal,1)
%%  Построение графиков
figure();
% plot(lambda, Spectrum, 'Linewidth', 2); hold on;
% plot(lambda,Signal/20,'LineWidth',1.5); hold on;
% plot(lambda,SignalError1/20); hold on;
plot(lambda,SignalError2/20); hold on;
plot(lambda, L_PsevdoInverse2, 'r'); 
plot(lambda, L_QR2, 'k'); 
plot(lambda, L_TikhonovSVD2,'-.g');
plot(lambda, L_Tikhonov2,'--b'); 
plot(lambda, L_TSVD2,'m');
plot(lambda, L_deg2/20,'c'); grid on;
axis([500 880 0 1]);
% legend ({'Spectrum','Signal','SignalError1','SignalError2','PsevdoInverse','QR','TikhonovSVD','Tikhonov', 'TSVD', 'deg8'});
legend({'Signal','PsevdoInverse','QR','TikhonovSVD','Tikhonov', 'TSVD', 'deg'});