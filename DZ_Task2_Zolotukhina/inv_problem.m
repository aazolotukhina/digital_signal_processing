function [L_deg3, L_deg7, L_nearest, L_spline, L_PsevdoInverse, L_QR,...
    L_TikhonovSVD, L_Tikhonov, L_TSVD, L_deg3_end,L_deg7_end,L_nearest_end,...
    L_spline_end,L_PsevdoInverse_end,L_QR_end,L_TikhonovSVD_end,...
    L_Tikhonov_end,L_TSVD_end] = inv_problem(signal_teor,signal_teor0,...
    lambdaMulti,lambda,ChannelMatrix,SpectrMatrix,regul,u)
% Функция inv_problem решает обратную задачу следующими методами:
% глобальная интерполяция; аппроксимация; локальная интерполяция; 
% метод псевдообратной матрицы Мура-Пенроуза; метод QR-разложения; 
% метод регуляризации Тихонова с SVD; метод регуляризации Тихонова; 
% метод псевдообратной матрицы c усеченным SVD.
% А также подсчитывает время получения решения каждым из рассматриваемых 
% методов.

% Входные данные:
% signal_teor - сигнал ПИ для объекта съемки
% signal_teor0 - сигнал ПИ для эталона
% lambdaMulti - длины волн для интерполяции и аппроксимации
% lambda - весь диапазон длин волн
% ChannelMatrix - чувствительность каналов ПИ
% SpectrMatrix - спектры объектов
% regul - параметр регуляризации для метода Тиханова
% u - кол-во используемых сингулярных чисел

% Выходные данные: 
% L_deg3, L_deg7 - результат аппроксимации и глобальной интерполяции
% L_nearest, L_spline - рез-т локально интерполяции методом ближайшего 
%                       соседа и методом сплайнов
% L_PsevdoInverse - рез-т решения методом псевдообратной матрицы Мура-Пенроуза
% L_QR - рез-т решения методом QR-разложения
% L_TikhonovSVD - рез-т решения методом регуляризации Тихонова с SVD
% L_Tikhonov - рез-т решения методом регуляризации Тихонова
% L_TSVD - рез-т решения методом псевдообратной матрицы c усеченным SVD
% L_deg3_end,L_deg7_end,L_nearest_end,L_spline_end,L_PsevdoInverse_end,...
% L_QR_end,L_TikhonovSVD_end, L_Tikhonov_end,L_TSVD_end - время получения 
% спектральной характеристики неизвестного текущего объекта

%%  Решение обратной задачи
% Аппроксимация и глобальная интерполяция
p_deg=3;    % аппроксимация

m=size(signal_teor,1);
p_deg2=m-1; % интерполяция

v3=polyfit(lambdaMulti, signal_teor, p_deg); % определение по табличным данным 
% функции (из SignalTeor) и аргумента (из lambdaMulti)
% коэффициентов аппроксимирующего их полинома p_deg-ной степени.
% При p_deg=m-1 - интерполяция, при p_deg<m-1 - аппроксимация
L_deg3_start = tic;
L_deg3=polyval(v3,lambda); % получение вектора значений полинома
L_deg3_end = toc(L_deg3_start);
% с коэффициентами, найденными на предыдущем этапе
v7=polyfit(lambdaMulti, signal_teor, p_deg2);
L_deg7_start = tic;
L_deg7=polyval(v7,lambda);
L_deg7_end = toc(L_deg7_start);

% Кусочная интерполяция 
L_nearest_start = tic;
L_nearest = interp1(lambdaMulti, signal_teor, lambda,'nearest'); 
L_nearest_end = toc(L_nearest_start);
L_spline_start = tic;
L_spline = interp1(lambdaMulti, signal_teor, lambda,'spline');
L_spline_end = toc(L_spline_start);
% интерполяция табличных данных функции (из SignalTeor1) и аргумента (из
% lambdaMulti). получение результатов в точках из lambda
% 'nearest' - метод ближайшего соседа
% 'spline' - метод сплайнов
signal_teor=signal_teor*max(signal_teor0); % действие, обратное калибровке сигнала

% Метод псевдообратной матрицы
PsevdoInverse_matrix=pinv(ChannelMatrix);
L_PsevdoInverse_start = tic;
L_PsevdoInverse=PsevdoInverse_matrix*signal_teor; 
% калибровка решения
L_PsevdoInverse1=PsevdoInverse_matrix*signal_teor0;
L_PsevdoInverse=L_PsevdoInverse.*(SpectrMatrix(13,:)'./L_PsevdoInverse1);
L_PsevdoInverse_end = toc(L_PsevdoInverse_start);

% осуществление svd разложения 
[U,s,V] = svd(ChannelMatrix); 
s=diag(s);

% QR разложение
[Q,R] = qr(ChannelMatrix);
[Q1,R1] = qr(ChannelMatrix);
L_QR_start = tic;
L_QR=R'*(Q'*signal_teor);
L_QR1=R1'*(Q1'*signal_teor0);
L_QR=L_QR.*(SpectrMatrix(13,:)'./L_QR1);
L_QR_end = toc(L_QR_start);

%Метод регуляризации Тихонова c SVD
L_TikhonovSVD_start = tic;
L_TikhonovSVD = tikhonovSVD(U,s,V, signal_teor, regul);
L_TikhonovSVD1 = tikhonovSVD(U,s,V, signal_teor0, regul);
L_TikhonovSVD = L_TikhonovSVD.*(SpectrMatrix(13,:)'./L_TikhonovSVD1);
L_TikhonovSVD_end = toc(L_TikhonovSVD_start);

%Метод регуляризации Тихонова
L_Tikhonov_start = tic;
L_Tikhonov =((regul*eye(381)+ChannelMatrix'*ChannelMatrix)\(ChannelMatrix'*signal_teor))';
L_Tikhonov1=((regul*eye(381)+ChannelMatrix'*ChannelMatrix)\(ChannelMatrix'*signal_teor0))';
L_Tikhonov = L_Tikhonov.*(SpectrMatrix(13,:)./L_Tikhonov1);
L_Tikhonov_end = toc(L_Tikhonov_start);

% Метод псевдообратной матрицы с усеченным SVD, u<m, u=m - псевдобратная матрица с SVD 
L_TSVD_start = tic;
L_TSVD=V(:,1:u)*inv(diag(s(1:u)))*U(:,1:u)'*signal_teor;
L_TSVD1=V(:,1:u)*inv(diag(s(1:u)))*U(:,1:u)'*signal_teor0;
L_TSVD=L_TSVD.*(SpectrMatrix(13,:)'./L_TSVD1);
L_TSVD_end = toc(L_TSVD_start);
end

