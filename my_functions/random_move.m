% Function to generate random move direction
function p = random_move(k)
a1 = 2*pi*rand;
r1 = k*sqrt(rand);
a2 = 2*pi*rand;
r2 = k*sqrt(rand);
a3 = 2*pi*rand;
r3 = k*sqrt(rand);
a4 = 2*pi*rand;
r4 = k*sqrt(rand);
a5 = 2*pi*rand;
r5 = k*sqrt(rand);
a6 = 2*pi*rand;
r6 = k*sqrt(rand);
a7 = 2*pi*rand;
r7 = k*sqrt(rand);

p = [r1*sin(a1); r1*cos(a1); r2*sin(a2); r2*cos(a2); r3*sin(a3); r3*cos(a3); r4*sin(a4); r4*cos(a4); ...
    r5*sin(a5); r5*cos(a5); r6*sin(a6); r6*cos(a6); r7*sin(a7); r7*cos(a7);];
end