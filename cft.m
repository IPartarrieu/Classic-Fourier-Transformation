function [A0,Aq,Bq] = cft(x)
% x es la serie de tiempo estacionaria
% dt en este caso es constante al paso de tiempo de la serie
N = length(x);
n = (1:N);
%% A_0
A0 = mean(x);
%
%%
% Calculamos theta de dimensiones qxn
for q=1:floor(N/2)-1
    theta(q,:)=((2*pi)/N).*q.*n;
end

%sin y cos de theta
cos_t = cos(theta);
sen_t = sin(theta);
%
%% Calculamos Aq y Bq
for q=1:floor(N/2)-1
    Aq(q,:) = (2/N)*sum(x.*cos_t(q,:));
    Bq(q,:) = (2/N)*sum(x.*sen_t(q,:));
end

% Cabe destacar que A_{N/2} tiene el siguiente valor
Aq(end) = (sum(x.*cos(n.*pi)))/N;
%
end