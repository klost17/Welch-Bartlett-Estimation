%% ACTIVITAT 5.1

h=[1 -0.8 0.64]; %resposta impulsional h[n]
P1=1; P2=0.5; omega1=pi/4; omega2=pi/2; sigma2=1; %paràmetres del procés

Nf=1024;%número de msotres frequencials dels estimadors, és una potencia de dos
eje_freq=2*pi*(0:(Nf-1))/Nf;%definim l'eix frequencial per a les posteriors representacions

H=fft(h,Nf);%calculem la FFT de la resposta impulsional

figure(1)

subplot(2,1,1)
plot(eje_freq,(abs(H)).^2); grid; xlim([0 2*pi])
title('DEP exacta (escala lineal)'); xlabel('\Omega')

subplot(2,1,2)
plot(eje_freq,20*log10(abs(H))); grid; xlim([0 2*pi])
title('DEP exacta (dB)'); xlabel('\Omega')

%% ACTIVITAT 5.2

N=250;
mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);

% Periodograma convencional
pConvencional=periodograma_modificat(Nf,mostres_x,ones(1,N));

% Periodograma modificat finestra Hanning
pModHanning=periodograma_modificat(Nf,mostres_x,hanning(N)');

% Estimador de Bartlett
estBartlett=bartlett_welch(Nf,mostres_x,ones(1,50),0);

% Estimador de Welch
estWelch=bartlett_welch(Nf,mostres_x,hamming(49)',50);

figure(2)

subplot(2,2,1)
plot(eje_freq,pConvencional); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('Periodograma Convencional'); xlabel('\Omega')

subplot(2,2,2)
plot(eje_freq,pModHanning); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('Periodograma Modificat'); xlabel('\Omega')

subplot(2,2,3)
plot(eje_freq,estBartlett); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('Estimador de Bartlett'); xlabel('\Omega')

subplot(2,2,4)
plot(eje_freq,estWelch); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('Estimador de Welch'); xlabel('\Omega')

%% ACTIVITAT 5.3

preAverage1=zeros(1000,Nf);
preAverage2=zeros(1000,Nf);
preAverage3=zeros(1000,Nf);
preAverage4=zeros(1000,Nf);

for j=1:1000
    
    mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);
    
    % Periodograma convencional
    pConvencional=periodograma_modificat(Nf,mostres_x,ones(1,N));
    preAverage1(j,:)=pConvencional;
    
    % Periodograma modificat finestra Hanning
    pModHanning=periodograma_modificat(Nf,mostres_x,hanning(N)');
    preAverage2(j,:)=pModHanning;
    
    % Estimador de Bartlett // amb 5 segments es fa la mitjana
    estBartlett=bartlett_welch(Nf,mostres_x,ones(1,50),0);
    preAverage3(j,:)=estBartlett;
    
    % Estimador de Welch
    estWelch=bartlett_welch(Nf,mostres_x,hamming(49)',50);
    preAverage4(j,:)=estWelch;
end

% 1000 estimacions solapades en escala lineal
figure(3)
subplot(2,4,1)
plot(eje_freq,preAverage1); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Periodograma Convencional'); xlabel('\Omega')
subplot(2,4,2)
plot(eje_freq,preAverage2); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Periodograma Modificat'); xlabel('\Omega')
subplot(2,4,3)
plot(eje_freq,preAverage3); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Estimador de Bartlett'); xlabel('\Omega')
subplot(2,4,4)
plot(eje_freq,preAverage4); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Estimador de Welch'); xlabel('\Omega')
% 1000 estimacions solapades en dB
subplot(2,4,5)
plot(eje_freq,10*log10(preAverage1)); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Periodograma Convencional (dB)');xlabel('\Omega')
subplot(2,4,6)
plot(eje_freq,10*log10(preAverage2)); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Periodograma Modificat (dB)');xlabel('\Omega')
subplot(2,4,7)
plot(eje_freq,10*log(preAverage3)); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Estimador de Bartlett (dB)');xlabel('\Omega')
subplot(2,4,8)
plot(eje_freq,10*log10(preAverage4)); xlim([0 2*pi]); grid; title('N=250, Nf=1024')
ylabel('1000 Estimador de Welch (dB)');xlabel('\Omega')

% mitjana en escala lineal
figure(4)
subplot(2,4,1)
plot(eje_freq,mean(preAverage1)); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Periodograma Convencional');xlabel('\Omega')
subplot(2,4,2)
plot(eje_freq,mean(preAverage2)); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Periodograma Modificat');xlabel('\Omega')
subplot(2,4,3)
plot(eje_freq,mean(preAverage3)); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Estimador de Bartlett');xlabel('\Omega')
subplot(2,4,4)
plot(eje_freq,mean(preAverage4)); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Estimador de Welch');xlabel('\Omega')
% mitjana en dB
subplot(2,4,5)
plot(eje_freq,10*log10(mean(preAverage1))); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Periodograma Convencional (dB)');xlabel('\Omega')
subplot(2,4,6)
plot(eje_freq,10*log10(mean(preAverage2))); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Periodograma Modificat (dB)');xlabel('\Omega')
subplot(2,4,7)
plot(eje_freq,10*log10(mean(preAverage3))); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Estimador de Bartlett (dB)');xlabel('\Omega')
subplot(2,4,8)
plot(eje_freq,10*log10(mean(preAverage4))); grid; xlim([0 2*pi]); title('N=250, Nf=1024')
ylabel('Mitjana 1000 Estimador de Welch (dB)');xlabel('\Omega')

% mitjana +- desviació típica
figure(5)
subplot(2,2,1)
plot(eje_freq,mean(preAverage1),eje_freq,mean(preAverage1)+...
    std(preAverage1),eje_freq,mean(preAverage1)-std(preAverage1));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Periodograma Convencional');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,2)
plot(eje_freq,mean(preAverage2),eje_freq,mean(preAverage2)+...
    std(preAverage2),eje_freq,mean(preAverage2)-std(preAverage2));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Periodograma Modificat');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,3)
plot(eje_freq,mean(preAverage3),eje_freq,mean(preAverage3)+...
    std(preAverage3),eje_freq,mean(preAverage3)-std(preAverage3));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Estimador de Bartlett');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,4)
plot(eje_freq,mean(preAverage4),eje_freq,mean(preAverage4)+...
    std(preAverage4),eje_freq,mean(preAverage4)-std(preAverage4));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Estimador de Welch');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')

%% ACTIVITAT 5.4

N=150;
preAverage1=zeros(1000,Nf);
preAverage2=zeros(1000,Nf);
preAverage3=zeros(1000,Nf);
preAverage4=zeros(1000,Nf);

for j=1:1000
    
    mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);
    
    % Periodograma convencional
    pConvencional=periodograma_modificat(Nf,mostres_x,ones(1,N));
    preAverage1(j,:)=pConvencional;
    
    % Periodograma modificat finestra Hanning
    pModHanning=periodograma_modificat(Nf,mostres_x,hanning(N)');
    preAverage2(j,:)=pModHanning;
    
    % Estimador de Bartlett // amb 5 segments es fa la mitjana
    estBartlett=bartlett_welch(Nf,mostres_x,ones(1,50),0);
    preAverage3(j,:)=estBartlett;
    
    % Estimador de Welch
    estWelch=bartlett_welch(Nf,mostres_x,hamming(49)',50);
    preAverage4(j,:)=estWelch;
end

% 1000 estimacions solapades en escala lineal
figure(6)
subplot(2,4,1)
plot(eje_freq,preAverage1); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Periodograma Convencional'); xlabel('\Omega')
subplot(2,4,2)
plot(eje_freq,preAverage2); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Periodograma Modificat'); xlabel('\Omega')
subplot(2,4,3)
plot(eje_freq,preAverage3); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Estimador de Bartlett'); xlabel('\Omega')
subplot(2,4,4)
plot(eje_freq,preAverage4); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Estimador de Welch'); xlabel('\Omega')
% 1000 estimacions solapades en dB
subplot(2,4,5)
plot(eje_freq,10*log10(preAverage1)); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Periodograma Convencional (dB)');xlabel('\Omega')
subplot(2,4,6)
plot(eje_freq,10*log10(preAverage2)); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Periodograma Modificat (dB)');xlabel('\Omega')
subplot(2,4,7)
plot(eje_freq,10*log(preAverage3)); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Estimador de Bartlett (dB)');xlabel('\Omega')
subplot(2,4,8)
plot(eje_freq,10*log10(preAverage4)); xlim([0 2*pi]); grid; title('N=150, Nf=1024')
ylabel('1000 Estimador de Welch (dB)');xlabel('\Omega')

% mitjana en escala lineal
figure(7)
subplot(2,4,1)
plot(eje_freq,mean(preAverage1)); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Periodograma Convencional');xlabel('\Omega')
subplot(2,4,2)
plot(eje_freq,mean(preAverage2)); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Periodograma Modificat');xlabel('\Omega')
subplot(2,4,3)
plot(eje_freq,mean(preAverage3)); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Estimador de Bartlett');xlabel('\Omega')
subplot(2,4,4)
plot(eje_freq,mean(preAverage4)); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Estimador de Welch');xlabel('\Omega')
% mitjana en dB
subplot(2,4,5)
plot(eje_freq,10*log10(mean(preAverage1))); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Periodograma Convencional (dB)');xlabel('\Omega')
subplot(2,4,6)
plot(eje_freq,10*log10(mean(preAverage2))); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Periodograma Modificat (dB)');xlabel('\Omega')
subplot(2,4,7)
plot(eje_freq,10*log10(mean(preAverage3))); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Estimador de Bartlett (dB)');xlabel('\Omega')
subplot(2,4,8)
plot(eje_freq,10*log10(mean(preAverage4))); grid; xlim([0 2*pi]); title('N=150, Nf=1024')
ylabel('Mitjana 1000 Estimador de Welch (dB)');xlabel('\Omega')

% mitjana +- desviació típica
figure(8)
subplot(2,2,1)
plot(eje_freq,mean(preAverage1),eje_freq,mean(preAverage1)+...
    std(preAverage1),eje_freq,mean(preAverage1)-std(preAverage1));
xlim([0 2*pi]);grid;title('N=150, Nf=1024')
ylabel('Periodograma Convencional');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,2)
plot(eje_freq,mean(preAverage2),eje_freq,mean(preAverage2)+...
    std(preAverage2),eje_freq,mean(preAverage2)-std(preAverage2));
xlim([0 2*pi]);grid;title('N=150, Nf=1024')
ylabel('Periodograma Modificat');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,3)
plot(eje_freq,mean(preAverage3),eje_freq,mean(preAverage3)+...
    std(preAverage3),eje_freq,mean(preAverage3)-std(preAverage3));
xlim([0 2*pi]);grid;title('N=150, Nf=1024')
ylabel('Estimador de Bartlett');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,4)
plot(eje_freq,mean(preAverage4),eje_freq,mean(preAverage4)+...
    std(preAverage4),eje_freq,mean(preAverage4)-std(preAverage4));
xlim([0 2*pi]);grid;title('N=150, Nf=1024')
ylabel('Estimador de Welch');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')



N=1000;
preAverage1=zeros(1000,Nf);
preAverage2=zeros(1000,Nf);
preAverage3=zeros(1000,Nf);
preAverage4=zeros(1000,Nf);

for j=1:1000
    
    mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);
    
    % Periodograma convencional
    pConvencional=periodograma_modificat(Nf,mostres_x,ones(1,N));
    preAverage1(j,:)=pConvencional;
    
    % Periodograma modificat finestra Hanning
    pModHanning=periodograma_modificat(Nf,mostres_x,hanning(N)');
    preAverage2(j,:)=pModHanning;
    
    % Estimador de Bartlett // amb 5 segments es fa la mitjana
    estBartlett=bartlett_welch(Nf,mostres_x,ones(1,50),0);
    preAverage3(j,:)=estBartlett;
    
    % Estimador de Welch
    estWelch=bartlett_welch(Nf,mostres_x,hamming(49)',50);
    preAverage4(j,:)=estWelch;
end

% 1000 estimacions solapades en escala lineal
figure(9)
subplot(2,4,1)
plot(eje_freq,preAverage1); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Periodograma Convencional'); xlabel('\Omega')
subplot(2,4,2)
plot(eje_freq,preAverage2); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Periodograma Modificat'); xlabel('\Omega')
subplot(2,4,3)
plot(eje_freq,preAverage3); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Estimador de Bartlett'); xlabel('\Omega')
subplot(2,4,4)
plot(eje_freq,preAverage4); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Estimador de Welch'); xlabel('\Omega')
% 1000 estimacions solapades en dB
subplot(2,4,5)
plot(eje_freq,10*log10(preAverage1)); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Periodograma Convencional (dB)');xlabel('\Omega')
subplot(2,4,6)
plot(eje_freq,10*log10(preAverage2)); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Periodograma Modificat (dB)');xlabel('\Omega')
subplot(2,4,7)
plot(eje_freq,10*log(preAverage3)); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Estimador de Bartlett (dB)');xlabel('\Omega')
subplot(2,4,8)
plot(eje_freq,10*log10(preAverage4)); xlim([0 2*pi]); grid; title('N=1000, Nf=1024')
ylabel('1000 Estimador de Welch (dB)');xlabel('\Omega')

% mitjana en escala lineal
figure(10)
subplot(2,4,1)
plot(eje_freq,mean(preAverage1)); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Periodograma Convencional');xlabel('\Omega')
subplot(2,4,2)
plot(eje_freq,mean(preAverage2)); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Periodograma Modificat');xlabel('\Omega')
subplot(2,4,3)
plot(eje_freq,mean(preAverage3)); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Estimador de Bartlett');xlabel('\Omega')
subplot(2,4,4)
plot(eje_freq,mean(preAverage4)); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Estimador de Welch');xlabel('\Omega')
% mitjana en dB
subplot(2,4,5)
plot(eje_freq,10*log10(mean(preAverage1))); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Periodograma Convencional (dB)');xlabel('\Omega')
subplot(2,4,6)
plot(eje_freq,10*log10(mean(preAverage2))); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Periodograma Modificat (dB)');xlabel('\Omega')
subplot(2,4,7)
plot(eje_freq,10*log10(mean(preAverage3))); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Estimador de Bartlett (dB)');xlabel('\Omega')
subplot(2,4,8)
plot(eje_freq,10*log10(mean(preAverage4))); grid; xlim([0 2*pi]); title('N=1000, Nf=1024')
ylabel('Mitjana 1000 Estimador de Welch (dB)');xlabel('\Omega')

% mitjana +- desviació típica
figure(11)
subplot(2,2,1)
plot(eje_freq,mean(preAverage1),eje_freq,mean(preAverage1)+...
    std(preAverage1),eje_freq,mean(preAverage1)-std(preAverage1));
xlim([0 2*pi]);grid;title('N=1000, Nf=1024')
ylabel('Periodograma Convencional');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,2)
plot(eje_freq,mean(preAverage2),eje_freq,mean(preAverage2)+...
    std(preAverage2),eje_freq,mean(preAverage2)-std(preAverage2));
xlim([0 2*pi]);grid;title('N=1000, Nf=1024')
ylabel('Periodograma Modificat');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,3)
plot(eje_freq,mean(preAverage3),eje_freq,mean(preAverage3)+...
    std(preAverage3),eje_freq,mean(preAverage3)-std(preAverage3));
xlim([0 2*pi]);grid;title('N=1000, Nf=1024')
ylabel('Estimador de Bartlett');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,4)
plot(eje_freq,mean(preAverage4),eje_freq,mean(preAverage4)+...
    std(preAverage4),eje_freq,mean(preAverage4)-std(preAverage4));
xlim([0 2*pi]);grid;title('N=1000, Nf=1024')
ylabel('Estimador de Welch');xlabel('\Omega')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')

%% ACTIVITAT 5.5

Nf=512; N=500;
eje_freq=2*pi*(0:(Nf-1))/Nf;

Average10=zeros(1000,Nf);
for j=1:1000
    
    mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);
    
    % Estimador de Bartlett 
    estBartlett=bartlett_welch(Nf,mostres_x,ones(1,10),0);% Q=10
    Average10(j,:)=estBartlett;
end
mean10=mean(Average10);
std10=std(Average10);

Average50=zeros(1000,Nf);
for j=1:1000
    
    mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);
    
    % Estimador de Bartlett 
    estBartlett=bartlett_welch(Nf,mostres_x,ones(1,50),0);% Q=50
    Average50(j,:)=estBartlett;
end
mean50=mean(Average50);
std50=std(Average50);

Average250=zeros(1000,Nf);
for j=1:1000
    
    mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);
    
    % Estimador de Bartlett 
    estBartlett=bartlett_welch(Nf,mostres_x,ones(1,250),0);% Q=250
    Average250(j,:)=estBartlett;
end
mean250=mean(Average250);
std250=std(Average250);


figure(12)
subplot(3,1,1)
plot(eje_freq,mean10,eje_freq,mean10+std10,eje_freq,mean10-std10);grid; xlim([0 2*pi])
xlabel('\Omega'); title('N=500, Nf=512, Q=10')
ylabel('Estimador de Bartlett')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica')
subplot(3,1,2)
plot(eje_freq,mean50,eje_freq,mean50+std50,eje_freq,mean50-std50);grid; xlim([0 2*pi])
xlabel('\Omega'); title('N=500, Nf=512, Q=50')
ylabel('Estimador de Bartlett')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(3,1,3)
plot(eje_freq,mean250,eje_freq,mean250+std250,eje_freq,mean250-std250);grid; xlim([0 2*pi])
xlabel('\Omega'); title('N=500, Nf=512, Q=250')
ylabel('Estimador de Bartlett')
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')

%% ACTIVITAT 5.6

% Deixem omega2=pi/2 fixada i fem que omega1 prengui valors entre pi/4 i
% pi/2. Fem un drawnow per veure com evoluciona el gràfic a mesura que
% ambdós valors s'apropen l'un a l'altre.

N=250; Nf=1024; eje_freq=2*pi*(0:(Nf-1))/Nf;
for omega1=linspace(pi/4,pi/2,30); omega2=pi/2;

preAverage1=zeros(1000,Nf);
preAverage2=zeros(1000,Nf);
preAverage3=zeros(1000,Nf);
preAverage4=zeros(1000,Nf);

for j=1:1000
    
    mostres_x=mostres_proces(N,h,sigma2,P1,P2,omega1,omega2);
    
    % Periodograma convencional
    pConvencional=periodograma_modificat(Nf,mostres_x,ones(1,N));
    preAverage1(j,:)=pConvencional;
    
    % Periodograma modificat finestra Hanning
    pModHanning=periodograma_modificat(Nf,mostres_x,hanning(N)');
    preAverage2(j,:)=pModHanning;
    
    % Estimador de Bartlett // amb 5 segments es fa la mitjana
    estBartlett=bartlett_welch(Nf,mostres_x,ones(1,50),0);
    preAverage3(j,:)=estBartlett;
    
    % Estimador de Welch
    estWelch=bartlett_welch(Nf,mostres_x,hamming(49)',50);
    preAverage4(j,:)=estWelch;
end

% mitjana +- desviació típica
figure(13)
subplot(2,2,1)
plot(eje_freq,mean(preAverage1),eje_freq,mean(preAverage1)+...
    std(preAverage1),eje_freq,mean(preAverage1)-std(preAverage1));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Periodograma Convencional');xlabel(['\Omega_1 =' num2str(omega1) '; \Omega_2 =' num2str(omega2)])
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,2)
plot(eje_freq,mean(preAverage2),eje_freq,mean(preAverage2)+...
    std(preAverage2),eje_freq,mean(preAverage2)-std(preAverage2));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Periodograma Modificat');xlabel(['\Omega_1 =' num2str(omega1) '; \Omega_2 =' num2str(omega2)])
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,3)
plot(eje_freq,mean(preAverage3),eje_freq,mean(preAverage3)+...
    std(preAverage3),eje_freq,mean(preAverage3)-std(preAverage3));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Estimador de Bartlett');xlabel(['\Omega_1 =' num2str(omega1) '; \Omega_2 =' num2str(omega2)])
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')
subplot(2,2,4)
plot(eje_freq,mean(preAverage4),eje_freq,mean(preAverage4)+...
    std(preAverage4),eje_freq,mean(preAverage4)-std(preAverage4));
xlim([0 2*pi]);grid;title('N=250, Nf=1024')
ylabel('Estimador de Welch');xlabel(['\Omega_1 =' num2str(omega1) '; \Omega_2 =' num2str(omega2)])
legend('Mitjana','Mitjana + Desviació típica','Mitjana - Desviació típica','Location','north')

drawnow
end