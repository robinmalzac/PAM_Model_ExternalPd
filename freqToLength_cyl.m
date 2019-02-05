% r = 0.007
% l = 0.660 max

r = 0.007;
c = 340;
% length_cyl_max= 660/1000;
fn = zeros(660,1);

for i = 100:660
    l = i/1000;
    [Fn,Yn,fni,Ze] = impedance_cyl(l,r);
    fn(i) = fni(1);
end


l = (100:660)/1000;

figure(1)
plot(l(1:1:length(l)), fn(100:1:660))
grid on
hold on
plot(l(1:1:length(l)), 2./(l(1:1:length(l))))
