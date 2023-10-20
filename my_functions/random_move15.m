% Function to generate random move direction
function p = random_move15(k)
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
a8 = 2*pi*rand;
r8 = k*sqrt(rand);
a9 = 2*pi*rand;
r9 = k*sqrt(rand);
a10 = 2*pi*rand;
r10 = k*sqrt(rand);
a11 = 2*pi*rand;
r11 = k*sqrt(rand);
a12 = 2*pi*rand;
r12 = k*sqrt(rand);
a13 = 2*pi*rand;
r13 = k*sqrt(rand);
a14 = 2*pi*rand;
r14 = k*sqrt(rand);
a15 = 2*pi*rand;
r15 = k*sqrt(rand);

p = [r1*sin(a1); r1*cos(a1); r2*sin(a2); r2*cos(a2); r3*sin(a3); r3*cos(a3); r4*sin(a4); r4*cos(a4); ...
    r5*sin(a5); r5*cos(a5); r6*sin(a6); r6*cos(a6); r7*sin(a7); r7*cos(a7); r8*sin(a8); r8*cos(a8); ...
    r9*sin(a9); r9*cos(a9); r10*sin(a10); r10*cos(a10); r11*sin(a11); r11*cos(a11); r12*sin(a12); r12*cos(a12); ...
    r13*sin(a13); r13*cos(a13); r14*sin(a14); r14*cos(a14); r15*sin(a15); r15*cos(a15)];
end