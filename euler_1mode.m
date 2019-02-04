function p1 = euler_1mode(gamma,zeta,time)
    fe = 44100;
    Te = 1/fe;
    
    l = 0.660;
    r = 0.007;
    
    [Fn,Yn,fn,~] = impedance_cyl(l,r);

    F0 = zeta*(1-gamma)*sqrt(gamma);
    A = zeta*(3*gamma-1)/(2*sqrt(gamma));
    B = -zeta*(3*gamma+1)/(8*gamma.^(3/2));
    C = -zeta*(gamma+1)/(16*gamma.^(5/2));
    
    p1 = zeros(time*fe,1); % p(t)
    p2 = p1; % p'(t)
    p1(1) = F0/(1-A);
    p2(1) = gamma/Te;

    for i = 1:time*fe-1
        p1(i+1) = p1(i) + Te*p2(i);
        p2(i+1) = p2(i) - Te*(Fn(1)*((Yn(1)-A) ...
                              -2*B*p1(i)-3*C*p1(i)^2)*p2(i) ...
                              +(2*pi*fn(1))^2*p1(i));
    end
end


