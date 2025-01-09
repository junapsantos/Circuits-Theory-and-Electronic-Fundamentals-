close all
clear all

pkg load symbolic

[a] = textread("../data.txt", "%f", 'delimiter', '= ', "endofline", "\n");
#a
R1 = a(28);
R2 = a(30);
R3 = a(32);
R4 = a(34);
R5 = a(36);
R6 = a(38);
R7 = a(40);
Vs = a(42);
C = a(44);
kb = a(46);
kd = a(48);
#printf("%12.11f", R1);

filename = "datatab1.tex";
f = fopen (filename, "w");
fprintf(f, "R1 (kOhm)& R2 (kOhm)& R3 (kOhm)& R4 (kOhm)& R5 (kOhm) \\\\ \n%e & %e & %e & %e & %e \\\\ \\hline \n", R1, R2, R3, R4, R5);
fclose (f);
filename = "datatab2.tex";
f = fopen (filename, "w");
fprintf(f, "R6 (kOhm)& R7 (kOhm)& Vs (V)& C ($\\mu$ F)& Kb (mS)& Kd (kOhm)\\\\ \n%e & %e & %e & %e & %e & %e  \\\\ \\hline \n", R6, R7, Vs, C, kb, kd);
fclose (f);

#Printing data to a .cir file

printf("printing table with data...\n")

filename = "octsim1.cir";
f = fopen (filename, "w");
fprintf(f, "R1 1 2 %12.11fK \nR2 2 3 %12.11fK \nR3 2 5 %12.11fK \nR4 5 0 %12.11fK \nR5 6 5 %12.11fK \nR6 7 0 %12.11fK \nR7 8 9 %12.11fK \n", R1, R2, R3, R4, R5, R6, R7);
fprintf(f, "\nVs 1 0 dc %12.11f \nC1 6 8 %12.11fU \nGb 6 3 2 5 %12.11fM \nV2 7 9 dc 0 \nHd 5 8 v2 %12.11fk \n", Vs, C, kb, kd);
fclose (f);

#Node Analysis
#1)t<0
printf("Node Analysis\n\n")

G1 = 1/(R1*1000);
G2 = 1/(R2*1000);
G3 = 1/(R3*1000);
G4 = 1/(R4*1000);
G5 = 1/(R5*1000);
G6 = 1/(R6*1000);
G7 = 1/(R7*1000);

C = C/10^6;
Kb = kb/1000;
Kd = kd*1000;

l1 = [1,0,0,0,0,0,0,0];
l2 = [-G1,G1+G2+G3,-G2,0,-G3,0,0,0];
l3 = [0,-G2-Kb, G2,0,Kb,0,0,0];
l4 = [0,0,0,1,0,0,0,0];
l5 = [0,0,0,-Kd*G6,1,0,Kd*G6,-1];
l6 = [0,Kb,0,0,-Kb-G5,G5,0,0];
l7 = [0,0,0,-G6,0,0, G6+G7,-G7];
l8 = [0,-G3,0,-G4,G3+G4+G5,-G5,-G7,G7];
A = [l1;l2;l3;l4;l5;l6;l7;l8];
b = [Vs;0;0;0;0;0;0;0];
V = A\b

Ia = (V(2) - V(1))*G1;
Ib = Kb*(V(2) - V(5));
Id = (V(5)-V(8))/Kd;


#tabel with the results 1

printf("\nprinting table with results...\n")

filename = "node1.tex";
f = fopen (filename, "w");
fprintf(f, "Is & %e \\\\ \\hline \nIb & %e \\\\ \\hline \nIc & 0.000000e+00 \\\\ \\hline \nId & %e \\\\ \\hline \n", Ia, Ib, Id);
fprintf(f, "V1 & %e \\\\ \\hline \nV2 & %e \\\\ \\hline \nV3 & %e \\\\ \\hline \nV4 & %e \\\\ \\hline \nV5 & %e \\\\ \\hline \nV6 & %e \\\\ \\hline \nV7 & %e \\\\ \\hline \nV8 & %e \\\\ \\hline \n", V)
fclose (f);

#2) 
Vx = V(6)-V(8)

filename = "octsim2.cir";
f = fopen (filename, "w");
fprintf(f, "R1 1 2 %12.11fK \nR2 2 3 %12.11fK \nR3 2 5 %12.11fK \nR4 5 0 %12.11fK \nR5 6 5 %12.11fK \nR6 7 0 %12.11fK \nR7 8 9 %12.11fK \n", R1, R2, R3, R4, R5, R6, R7);
fprintf(f, "\nVs 1 0 dc 0 \nVx 6 8 %12.11f \nGb 6 3 2 5 %12.11fM \nV2 7 9 dc 0 \nHd 5 8 v2 %12.11fk \n", Vx, kb, kd);
fclose (f);

l1 = [1,0,0,0,0,0,0,0];
l2 = [-G1,G1+G2+G3,-G2,0,-G3,0,0,0];
l3 = [0,-G2-Kb, G2,0,Kb,0,0,0];
l4 = [0,0,0,1,0,0,0,0];
l5 = [0,0,0,-Kd*G6,1,0,Kd*G6,-1];
l6 = [0,0,0,0,0,1,0,-1];
l7 = [0,0,0,-G6,0,0, G6+G7,-G7];
l8 = [0,Kb-G3,0,-G4,G3+G4-Kb,0,-G7,G7];
A = [l1;l2;l3;l4;l5;l6;l7;l8];
b = [0;0;0;0;0;Vx;0;0];

V1 = A\b

#Determinar Ix
Ix = Kb*(V1(2)-V1(5)) + G5*(V1(6)-V1(5))

Req = Vx/Ix
Tau = Req*C

filename = "node2.tex";
f = fopen (filename, "w");
#fprintf(f, "Is & %e \\\\ \\hline \nIb & %e \\\\ \\hline \nIc & 0.000000e+00 \\\\ \\hline \nId & %e \\\\ \\hline \n", Ia, Ib, Id);
fprintf(f, "V1 & %e \\\\ \\hline \nV2 & %e \\\\ \\hline \nV3 & %e \\\\ \\hline \nV4 & %e \\\\ \\hline \nV5 & %e \\\\ \\hline \nV6 & %e \\\\ \\hline \nV7 & %e \\\\ \\hline  \nV8 & %e \\\\ \\hline \n", V1)
fprintf(f, "Ix & %e \\\\ \\hline \n$R_{eq}$ & %e \\\\ \\hline \nRC & %e \\\\ \\hline \n", Ix, Req, Tau);
fclose (f);

#3 Natural Solution
filename = "octsim3.cir";
f = fopen (filename, "w");
fprintf(f, "R1 1 2 %12.11fK \nR2 2 3 %12.11fK \nR3 2 5 %12.11fK \nR4 5 0 %12.11fK \nR5 6 5 %12.11fK \nR6 7 0 %12.11fK \nR7 8 9 %12.11fK \n", R1, R2, R3, R4, R5, R6, R7);
fprintf(f, "\nVs 1 0 dc 0 \nC1 6 8 %12.11fU ic = %12.11f \nGb 6 3 2 5 %12.11fM \nV2 7 9 dc 0 \nHd 5 8 v2 %12.11fk \n", C*10^6, Vx, kb, kd);
fprintf(f, "\n.ic v(6)=%12.11f v(8)=0.0\n", V1(6));
fclose (f);

t = 0:0.0001:0.02;
hf = figure ();
hold on;
plot (t, V1(6)*exp(-t/Tau), "color", "r", "linewidth", 3);
axis ([0, 0.02, 0, 9]);
xlabel ("t(s)");
ylabel ("V6(V)");
title ("V6(t) - Natural Solution");
print -djpg V6n.jpg

#4 Forced Solution

filename = "octsim4.cir";
f = fopen (filename, "w");
fprintf(f, "R1 1 2 %12.11fK \nR2 2 3 %12.11fK \nR3 2 5 %12.11fK \nR4 5 0 %12.11fK \nR5 6 5 %12.11fK \nR6 7 0 %12.11fK \nR7 8 9 %12.11fK \n", R1, R2, R3, R4, R5, R6, R7);
fprintf(f, "\nVs 1 0 0.0 ac 1.0 sin (0.0 1.0 1k) \nC1 6 8 %12.11fU ic = %12.11f \nGb 6 3 2 5 %12.11fM \nV2 7 9 dc 0 \nHd 5 8 v2 %12.11fk \n", C*10^6, Vx, kb, kd);
fprintf(f, "\n.ic v(6)=%12.11f v(8)=0.0\n", V1(6));
fclose (f);

f = 1000;
w = 2*pi*f

l1 = [1,0,0,0,0,0,0,0];
l2 = [-G1,G1+G2+G3,-G2,0,-G3,0,0,0];
l3 = [0,-G2-Kb, G2,0,Kb,0,0,0];
l4 = [0,0,0,1,0,0,0,0];
l5 = [0,0,0,-Kd*G6,1,0,Kd*G6,-1];
l6 = [0,Kb,0,0,-Kb-G5,G5+C*w*j,0,-C*w*j];
l7 = [0,0,0,-G6,0,0, G6+G7,-G7];
l8 = [0,-G3,0,-G4,G3+G4+G5,-G5-C*w*j,-G7,G7+C*w*j];
A = [l1;l2;l3;l4;l5;l6;l7;l8];
b = [1;0;0;0;0;0;0;0];

V2 = A\b

A = abs(V2(6))
fase = arg(V2(6))
printf("Amplitude:\n");

Amplitudes = [
abs(V2(1));
abs(V2(2));
abs(V2(3));
abs(V2(4));
abs(V2(5));
abs(V2(6));
abs(V2(7));
abs(V2(8));
]

printf("Fases:\n");
Fases = [
arg(V2(1));
arg(V2(2));
arg(V2(3));
arg(V2(4));
arg(V2(5));
arg(V2(6));
arg(V2(7));
arg(V2(8));
]

printf("\nprinting table with results...\n")

filename = "amplitudes.tex";
f = fopen (filename, "w");
fprintf(f, "A1 & %e \\\\ \\hline \nA2 & %e \\\\ \\hline \nA3 & %e \\\\ \\hline \nA4 & %e \\\\ \\hline \nA5 & %e \\\\ \\hline \nA6 & %e \\\\ \\hline \nA7 & %e \\\\ \\hline \nA8 & %e \\\\ \\hline \n", Amplitudes);
fclose (f);

printf("\nprinting table with results...\n")

filename = "fases.tex";
f = fopen (filename, "w");
fprintf(f, "$\\phi$1 & %e \\\\ \\hline \n$\\phi$2 & %e \\\\ \\hline \n$\\phi$3 & %e \\\\ \\hline \n$\\phi$4 & %e \\\\ \\hline \n$\\phi$5 & %e \\\\ \\hline \n$\\phi$6 & %e \\\\ \\hline \n$\\phi$7 & %e \\\\ \\hline \n$\\phi$8 & %e \\\\ \\hline \n", Fases);
fclose (f);


#5 Final Solution

filename = "octsim5.cir";
f = fopen (filename, "w");
fprintf(f, "R1 1 2 %12.11fK \nR2 2 3 %12.11fK \nR3 2 5 %12.11fK \nR4 5 0 %12.11fK \nR5 6 5 %12.11fK \nR6 7 0 %12.11fK \nR7 8 9 %12.11fK \n", R1, R2, R3, R4, R5, R6, R7);
fprintf(f, "\nVs 1 0 0.0 ac 1.0 sin (0.0 1.0 1k) \nC1 6 8 %12.11fU ic = %12.11f \nGb 6 3 2 5 %12.11fM \nV2 7 9 dc 0 \nHd 5 8 v2 %12.11fk \n", C*10^6, Vx, kb, kd);
fprintf(f, "\n.ic v(6)=%12.11f v(8)=0.0\n", V1(6));
fclose (f);

t = 0:0.0001:0.02;
hf = figure ();
plot (t,V1(6)*exp(-t/Tau) + A*sin(w*t + fase), "color", "r", "linewidth", 3);

hold on;

t = -0.005:0.0001:0;
plot (t,V(6) + 0*t, "color", "r", "linewidth", 3);

hold on;

t = 0:0.0001:0.02;
plot (t,sin(w*t), "color", "b", "linewidth", 3);

hold on;

t = -0.005:0.0001:0;
plot (t,V(1) + 0*t, "color", "b", "linewidth", 3);
axis ([-0.005, 0.02, -1, 12]);
xlabel ("t(s)");
ylabel ("V(V)");
title ("Vs and V6 - Final Solution");
print -djpg VsV6f.jpg

#6 
f = logspace(-1, 6, 7*15);

for i = 1 : 7*15
w = 2*pi*f(i);

l1 = [1,0,0,0,0,0,0,0];
l2 = [-G1,G1+G2+G3,-G2,0,-G3,0,0,0];
l3 = [0,-G2-Kb, G2,0,Kb,0,0,0];
l4 = [0,0,0,1,0,0,0,0];
l5 = [0,0,0,-Kd*G6,1,0,Kd*G6,-1];
l6 = [0,Kb,0,0,-Kb-G5,G5+C*w*j,0,-C*w*j];
l7 = [0,0,0,-G6,0,0, G6+G7,-G7];
l8 = [0,-G3,0,-G4,G3+G4+G5,-G5-C*w*j,-G7,G7+C*w*j];
A = [l1;l2;l3;l4;l5;l6;l7;l8];
b = [1;0;0;0;0;0;0;0];

V3 = A\b;

Vc = V3(6) - V3(8);

amplitude6(i) = abs(V3(6));
amplitudeC(i) = abs(Vc);

fase6(i) = arg(V3(6));
faseC(i) = arg(Vc);

endfor

n = 1:1:7*15;
hf = figure ();
plot (log10(f(n)),20*log10(amplitude6(n)), "color", "r", "linewidth", 3);
hold on;
plot (log10(f(n)),20*log10(amplitudeC(n)), "color", [1, 0.49, 0.314], "linewidth", 3);
hold on;
plot (log10(f(n)), 1 + 0*n , "color", "b", "linewidth", 3);

#axis ([-0.005, 0.02, -1, 12]);
xlabel ("log(f(Hz))");
ylabel ("Magnitude (dB)");
title ("Vs,  V6,  Vc - Magnitude");
print -djpg V6Amplitude1.jpg

n = 1:1:7*15;
hf = figure ();
plot (log10(f(n)),180/pi * fase6(n), "color", "r", "linewidth", 3);
hold on;
plot (log10(f(n)),180/pi * faseC(n), "color", [1, 0.49, 0.314], "linewidth", 3);
hold on;
plot (log10(f(n)), 0*n , "color", "b", "linewidth", 3);
xlabel ("log(f(Hz))");
ylabel ("Phase (Degrees)");
title ("Vs,  V6,  Vc - Phase");
axis ([-1, 6, -200, 10]);
#axis ([-0.005, 0.02, -1, 12]);
print -djpg V6fASE1.jpg
