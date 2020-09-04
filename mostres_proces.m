function mostres_x=mostres_proces(N,h,sigma2,P1,P2,Omeg1,Omeg2)
L=length(h);
Lw=N+L-1;
w=normrnd(0,sqrt(sigma2),1,Lw);% generem el soroll w[n] com una variable gaussiana de mitja 0 i vari�ncia sigma al quadrat
y=conv(w,h); %calculem la convoluci� de w[n] amb h[h]
y=y(L:end-L+1);% traiem les L-1 primeres i �ltimes mostres que corresponen a transitorios
s=sqrt(2*P1)*cos(Omeg1*[0:N-1]+2*pi*rand)+sqrt(2*P2)*cos(Omeg2*[0:N-1]+2*pi*rand);%generem el senyal s[n] amb fases les quals son variables aleat�ries uniformes[0,2p�]
mostres_x=y+s;% creem la sortida x[n] sumant y[n] i s[n]
end