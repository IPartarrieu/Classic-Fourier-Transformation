%%
clear all
% Generamos una serie de armónicos con distintos periodos
dt = 1;             % Delta-t de la serie
N = 1000;           % Largo de la serie
t = (0:N-1)'*dt;    % Vector de tiempos
fr = (1./t)';

% Formamos 5 señales con periodos entre 10 y 100
% (frecuencias angulares)
w1 = 2*pi/10*dt;
w2 = 2*pi/30*dt;
w3 = 2*pi/50*dt;
w4 = 2*pi/70*dt;
w5 = 2*pi/100*dt;

% Series de amplitudes 1, 3, 5, 7, 10

ss(1,:) = 0.1*sin(w1*t);
ss(2,:) = 0.3*sin(w2*t);
ss(3,:) = 0.5*sin(w3*t);
ss(4,:) = 0.7*sin(w4*t);
ss(5,:) = 1*sin(w5*t);

% Serie sintética de amplitudes 1, 3, 5, 7, 10
S = ss(1,:) + ss(2,:) + ss(3,:) + ss(4,:) + ss(5,:);
%
clear w1 w2 w3 w4 w5
%
%% EJERCICIO 3

% Dividimos en 5 segmentos
seg = reshape(S,200,5);
seg = seg';
t_s = reshape(t,200,5);
t_s = t_s';
% Aplicamos una ventana de Hanning a cada segmento
for i=1:5
    seg_h(i,:) = seg(i,:)'.*hann(200);
end
%
%% Ploteamos resultado
close all
cont=1;
figure()
for i=1:5
subplot(5,2,cont)
plot(t_s(i,:),seg(i,:),LineWidth=1.5)
xlabel('Tiempo [s]','FontSize',12)
ylabel('Amplitud [m]','FontSize',12)
title(['Segmento Original: ', num2str(i)],'FontSize',12)
ylim([-2.5 2.5])
grid minor
cont = cont+1;
%
subplot(5,2,cont)
plot(t_s(i,:),seg_h(i,:),LineWidth=1.5)
xlabel('Tiempo [s]','FontSize',12)
ylabel('Amplitud [m]','FontSize',12)
title(['Segmento con Ventana: ', num2str(i)],'FontSize',12)
ylim([-2.5 2.5])
grid minor
cont = cont+1;
end


%% Calculamos varianza de segmentos
clear seg_var EV_seg EV_seg_h Aq_sh Bq_sh Aq_s Bq_s
% Calculamos varianza con la función
for i=1:5
    seg_var(i,1) = round(var(seg(i,:)),2);
    seg_var(i,3) = round(var(seg_h(i,:)),2);
end

% Calculamos los coeficientes de los segmentos originales
for i=1:5
    [~,Aq_s(:,i),Bq_s(:,i)] = cft(seg(i,:));
    EV_seg(i,:) = (Aq_s(:,i).^2+Bq_s(:,i).^2)/2;
    seg_var(i,2) = round(sum(EV_seg(i,:)),2);
    %
    [~,Aq_sh(:,i),Bq_sh(:,i)] = cft(seg_h(i,:));
    EV_seg_h(i,:) = (Aq_sh(:,i).^2+Bq_sh(:,i).^2)/2;
    seg_var(i,4) = round(sum(EV_seg_h(i,:)),2);
end
Aq_s = Aq_s';
Bq_s = Bq_s';
Aq_sh = Aq_sh';
Bq_sh = Bq_sh';
%
%% Multiplicamos por el factor
% De las clases, ocupamos el ejemplo de sqrt(8/2)
F = sqrt(8/3);
Aq_sh_aj = Aq_sh.*F;
Bq_sh_aj = Bq_sh.*F;
%
% Clculamos varianza con las componentes ajustadas

for i=1:5
    EV_seg_aj(i,:) = (Aq_sh_aj(i,:).^2+Bq_sh_aj(i,:).^2)/2;
    seg_var(i,5) = round(sum(EV_seg_aj(i,:)),2);
end
%
%% Calculamos espectro promedio
fr_x = (1:99)'/200;
A_x = mean(Aq_sh_aj);
B_x = mean(Bq_sh_aj);
EV_x = (A_x.^2+B_x.^2)/200;
PSD = EV_x*200;

close all
figure()
subplot(2,1,1)
plot(fr_x,EV_x,LineWidth=1.5)
xlabel('Frecuencia [Hz]','FontSize',12)
ylabel('Amplitud','FontSize',12)
title('Método de Welch','FontSize',12)
grid minor
%
subplot(2,1,2)
plot(fr_x,PSD,LineWidth=1.5)
xlabel('Frecuencia [Hz]','FontSize',12)
ylabel('Potencia [Hz^-^1]','FontSize',12)
title('Densidad Espectral','FontSize',12)
grid minor
%