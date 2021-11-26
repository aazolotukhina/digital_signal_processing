function [t1] = gaussfunction(x,lam,sko, A)
% Функция создает массив функций гаусса с заданными параметрами
% x = область определения
% lambda=мат. ожидание
% sko=СКО
if size(lam)==size(sko) & size(sko)==size(A)
    for i=1:size(lam,2)
        t1(i,:)=A(i)*exp(-((x-lam(i)).^2)/(2*sko(i).^2));
    end
else
    disp ('error size of vectors')
end
end

