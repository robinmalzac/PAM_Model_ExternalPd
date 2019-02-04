function [Fn,Yn,fn,Ze] = impedance_cyl(l,r)
    
    c = 340;
    rho = 1.125;
    w = linspace(2*pi*10,2*pi*3000,1000);
    k = w/c;
    S = pi*r^2;
    
    Zc = rho*c/S; % imp�dance caract�ristique
    Zs = Zc*(1j*k*0.6*r + 0.25*(k*r^2)); % Impedance de rayonnement
    GAMMA = 1j*w/c +(1+1j)*3e-5*sqrt(w/(2*pi))/r; % amortissement et dissipation
    Ze = Zc*tanh(GAMMA*l + atanh(Zs/Zc)); % Impedance d'entr�e
    [ZMn,ind] = findpeaks(abs(Ze));
    Yn = 1./ZMn;
    fn = w(ind)/(2*pi);
    Fn = ones(length(Yn),1)*2*c/l;
end