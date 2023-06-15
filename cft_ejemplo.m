clear all
%% EJERCICIO 2

% Generamos una serie de armónicos con distintos periodos
dt = 1;             % Delta-t de la serie
N = 1000;           % Largo de la serie
t = (0:N-1)'*dt;    % Vector de tiempos
fr = (1./t)';
% Formamos 5 señales con periodos entre 10 y 100
% (frecuencias angulares)
w1 = 2*pi/10*dt;    % 0.6283
w2 = 2*pi/30*dt;    % 0.2094
w3 = 2*pi/50*dt;    % 0.1257
w4 = 2*pi/70*dt;    % 0.0898
w5 = 2*pi/100*dt;   % 0.0628

% Series de amplitudes 0.1, 0.3, 0.5, 0.7 y 1

ss(1,:) = 0.1*sin(w1*t);
ss(2,:) = 0.3*sin(w2*t);
ss(3,:) = 0.5*sin(w3*t);
ss(4,:) = 0.7*sin(w4*t);
ss(5,:) = 1*sin(w5*t);

% Serie sintética de amplitudes 0.1, 0.3, 0.5, 0.7 y 1
S = ss(1,:) + ss(2,:) + ss(3,:) + ss(4,:) + ss(5,:);
%
clear w1 w2 w3 w4 w5
%
%%
% Agregamos ruido distribuido normalmente
SR = S' + 2*randn(N,1);
SR = SR';

%% Ploteamos ambas series (S y SR)
close all
figure()
subplot(2,1,1)
plot(t,S,LineWidth=1.5)
xlabel('Tiempo [s]','FontSize',12)
ylabel('Amplitud [m]','FontSize',12)
title('Serie Sintética','FontSize',12)
axis([t(1) t(end) -3 3])
grid minor
subplot(2,1,2)
plot(t,SR,LineWidth=1.5)
xlabel('Tiempo [s]','FontSize',12)
ylabel('Amplitud [m]','FontSize',12)
title('Serie Sintética con Ruido','FontSize',12)
axis([t(1) t(end) -10 10])
grid minor
%
%% Ploteamos cada armónico 
close all
figure()
%
for i=1:5
subplot(5,1,i)
plot(t,ss(i,:),LineWidth=1.5)
xlabel('Tiempo [s]','FontSize',12)
ylabel('Amplitud [m]','FontSize',12)
title(['Armónico ',num2str(i)],'FontSize',12)
grid minor
axis([t(1) t(end) -1.5 1.5])
end
%
%% Tabla de varianzas
var_table(1,1) = round(var(S),2);
for i=1:5
var_table(i+1,1) = round(var(ss(i,:)),2);
end
var_table(7,1) = round(var(SR),2);
%
%% Varianza mediante transformadas de fourier
% C = sqrt(Aq.^2+Bq.^2);
clear EV
% Serie sintetica
[~,Aq,Bq] = cft(S);
EV(1,:) = (Aq.^2+Bq.^2)/2;
var_table(1,2) = round(sum(EV(1,:)),2);
%
% Series componente
for i=1:5
    clear Aq Bq
    [~,Aq,Bq] = cft(ss(i,:));
    EV(i+1,:) = (Aq.^2+Bq.^2)/2;
    var_table(i+1,2) = round(sum(EV(i+1,:)),2);
end
%
% Serie sintética ruidosa
clear Aq Bq
[~,Aq,Bq] = cft(SR);
EV(7,:) = (Aq.^2+Bq.^2)/2;
var_table(7,2) = round(sum(EV(7,:)),2);
%
%% Graficamos para ambas series
close all
figure()
subplot(2,1,1)
plot(fr(1:499),EV(1,:),LineWidth=1.5)
xlabel('Frecuencia [s^-^1]','FontSize',12)
ylabel('Varianza [c_k^2/2]','FontSize',12)
title('Serie Sintética','FontSize',12)
grid minor
%
subplot(2,1,2)
plot(fr(1:499),EV(7,:),LineWidth=1.5)
xlabel('Frecuencia [s^-^1]','FontSize',12)
ylabel('Varianza [c_k^2/2]','FontSize',12)
title('Serie Sintética Ruidosa','FontSize',12)
grid minor
sgtitle('Espectro de Varianza','FontSize',16)
%
%% Graficamos para cada armónico
close all
figure()
%
for i=1:5
subplot(5,1,i)
plot(fr(1:499),EV(i+1,:),LineWidth=1.5)
xlabel('Frecuencia [s^-^1]','FontSize',12)
ylabel('Varianza [c_k^2/2]','FontSize',12)
title(['Armónico ',num2str(i)],'FontSize',12)
grid minor
axis([0 0.5 min(EV(i+1,:)) max(EV(i+1,:))])
end
sgtitle('Espectro de Varianza','FontSize',16)
%