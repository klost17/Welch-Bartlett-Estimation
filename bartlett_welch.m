function mostres_bartlett_welch=bartlett_welch(Nf,x,v,solap)
Q=length(v); norma=1/Q*sum(v.^2); v=v/sqrt(norma); % normalitzem la finestra
D=round(Q*(1-solap/100));%calculem el nombre de mostres entre segments consecutius
N=length(x);
K=floor(N/Q);%calculem el nombre de segments necessaris, hi haurà N-K*Q mostres que no seran utilitzades si N/Q no es enter
preAverage=zeros(1,Nf);
for k=0:K-1
    kSegment = x((1+k*D):(Q+k*D));%agafem el segment de Q mostres
    x_v_k = kSegment.*v;%multipliquem el segment per la finestra
    r_x_k = 1/Q*conv(x_v_k,conj(flip(x_v_k)));%estimem l'autocorrelació de les dades del segment
    period_mod=abs(fft(r_x_k,Nf));%calculem el periodograma modificat del segment: FFT de l'autocorrelació estimada amb el segment
    preAverage=preAverage+period_mod;%sumem els periodogrames modificats de tots els segments
end
mostres_bartlett_welch=1/K*preAverage;%fem la mitjana dels periodogrames modificats de tots els segments
end