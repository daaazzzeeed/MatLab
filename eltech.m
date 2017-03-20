I = zeros(2,1);
R = zeros(2,2);
E = zeros(2,1);
R1 = 100;
R5 = 100;
R2 = 50;
R3 = 50;
R4 = 50;
E1 = 400;
E4 = 200;
E3 = 400;
Y = 6;
I33 = -Y;

R(1,1) = R3 + R4 + R5;
R(1,2) = - R5;
R(2,1) = - R5;
R(2,2) = R1 + R5;

E(1,1) = E3 + E4 + I33*R3;
E(2,1) = E1 + I33*R1;

invR = R^(-1);

I = invR*E;

fprintf('I11 = %2f',I(1,1));
fprintf('\n');
fprintf('I22 = %2f',I(2,1));
fprintf('\n');
fprintf('I33 = %2f',I33);
fprintf('\n');

Fi = zeros(2,1);
G  = zeros(2,2);
YY  = zeros(2,1);

G11 = 1/R1 + 1/R3 + 1/R5;
G22 =1/R2 + 1/R3 + 1/R4;

G12 = -1/R3;
G23 = - 1/R4;
G21 = G12;
G13 = -1/R5;

Y11 = -E3/R3;
Y22 = -Y + E3/R3 - E4/R4;

YY(1,1) = Y11;
YY(2,1) = Y22;

G(1,1) = G11;
G(1,2) = G12;
G(2,1) = G21;
G(2,2) = G22;


invG = (G)^(-1);

Fi = invG*YY;

fprintf('fi1 = %2f',Fi(1,1));
fprintf('\n');
fprintf('fi2 = %2f',Fi(2,1));
fprintf('\n');
fprintf('fi3 = %2f',E1);
fprintf('\n');





