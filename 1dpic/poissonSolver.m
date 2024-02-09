function phi=poissonSolver(rho,filter)
len=length(rho)-1;
kk=[1:len/2.,floor(-len/2):-1];
phi=fft(rho);
phi=-ifft(filter.*[0.,phi(2:end)./kk.^2.*(sin(kk/(2.*len))./(kk/(2.*len))).^2]);
end