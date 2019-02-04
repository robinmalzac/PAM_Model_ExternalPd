function p = euler_nmodes(gamma,ksi,time,n)
    fe = 44100;
    Te = 1/fe;
    
    l = 0.660;
    r = 0.015;
    
    [Fn,Yn,fn,~] = impedance_cyl(l,r);

    F0 = ksi*(1-gamma)*sqrt(gamma);
    A = ksi*(3*gamma-1)/(2*sqrt(gamma));
    B = -ksi*(3*gamma+1)/(8*gamma.^(3/2));
    C = -ksi*(gamma+1)/(16*gamma.^(5/2));
    
    p_0 = F0/(1-A);
    dp_0 = 0;
    p_1 = p_0 + Te*dp_0;
    p_mod = zeros(time*fe,n);
    p_tot = zeros(time*fe,2);

    p_mod(1,:) = p_0/n;
    p_mod(2,:) = p_1/n;

    for i=3:time*fe
        p_tot(i-2,1) = sum(p_mod(i-2,:));
        p_tot(i-2,2) = 1/Te*sum(p_mod(i-1,:)-p_mod(i-2,:));
        for j = 1:n
            p_mod(i,j) = 2*p_mod(i-1,j) - p_mod(i-2,j) ...
                - Te*Fn(j)*Yn(j)*(p_mod(i-1,j)-p_mod(i-2,j)) ...
                + Fn(j)*Te^2*p_tot(i-2,2)*(A+2*B*p_tot(i-2,1) ...
                + 3*C*p_tot(i-2,1)^2)-(Te*2*pi*fn(j))^2*p_mod(i-2,j);
        end
    end
    
    p = p_tot(:,1);
end