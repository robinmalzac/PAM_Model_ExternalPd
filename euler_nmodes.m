function p = euler_nmodes(gamma,zeta,time,n)
    fe = 44100;
    Te = 1/fe;
    
    l = 0.660;
    r = 0.007;
    
    [Fn,Yn,fn,~] = impedance_cyl(l,r);
    Fn = Fn(1:n);
    Yn = Yn(1:n);
    fn = fn(1:n);

    F0 = zeta*(1-gamma)*sqrt(gamma);
    A = zeta*(3*gamma-1)/(2*sqrt(gamma));
    B = -zeta*(3*gamma+1)/(8*gamma.^(3/2));
    C = -zeta*(gamma+1)/(16*gamma.^(5/2));
    
    pn1 = zeros(time*fe,n); % pn(t)
    pn2 = pn1; % pn'(t)
    pn1(1,:) = 0.0001; % initialization
    
    p1tot = zeros(time*fe,1); % p(t) total
    p2tot = p1tot; % p'(t) total
    
    for i = 1:time*fe-1
        p1tot(i) = sum(pn1(i,:));
        p2tot(i) = sum(pn2(i,:));
        for j = 1:n
            pn1(i+1,j) = pn1(i,j) + Te*pn2(i,j);
            pn2(i+1,j) = pn2(i,j) - Te*(Fn(j)*(Yn(j)*pn2(i,j) ...
                                         -p2tot(i)*(A+2*B*p1tot(i)+3*C*p1tot(i)^2)) ...
                                         +(2*pi*fn(j))^2*pn1(i,j));
        end
    end
    
    p = p1tot/n;
end
