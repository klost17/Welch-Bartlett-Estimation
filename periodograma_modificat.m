function mostres_period_mod=periodograma_modificat(Nf,x,v)
N = length(v);
norma = 1/N*sum(v.^2);
v = v/sqrt(norma);%normalitzem la finestra
x_v = x.*v;%multipliquem el vector que conté les mostres de la realització del proces per la finestra 
r_x = 1/N*conv(x_v,conj(flip(x_v)));%estimem l'autocorrelació 
mostres_period_mod = abs(fft(r_x,Nf));%calculem DEP estimada: FFT de l'autocorrelació estimada
end