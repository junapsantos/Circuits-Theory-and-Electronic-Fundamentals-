close all
clear all

pkg load symbolic

#Mesh Analysis
printf("Mesh Analysis\n\n")

syms Id;

R1 = 1.01072152255*1000; 
R2 = 2.0468229318*1000; 
R3 = 3.05917032304*1000; 
R4 = 4.1965026363*1000;
R5 = 3.07879619547*1000;
R6 = 2.08425994597*1000;
R7 = 1.01419841564*1000;

Va = 5.03865114929;
Id = 1.01596988715/1000;
Kb = 7.14516521228/1000;
Kc = 8.13498446601*1000;

#Matrix construction

l1 = [R1+R3+R4,R3,R4];
l2 = [Kb*R3,Kb*R3-1,0];
#l2 = [R3,R2+R3+R5-(1/Kb),0];
l3 = [R4,0,R4+R6+R7-Kc];

A = [l1;l2;l3]
b = [Va;0;0]

#Computing linear sistem

syms Is;
Is = A\b


#Node Analysis
printf("Node Analysis\n\n")

syms Vs;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G5 = 1/R5;
G6 = 1/R6;
G7 = 1/R7;

Id = 1.01596988715/1000;

Va = 5.03865114929;
kb = 7.14516521228/1000;
kc = 8.13498446601*1000;

l1 = [1,0,0,0,0,0,0];
l2 = [-G1,G1+G2+G3,-G2,0,-G3,0,0];
l3 = [0,-G2-kb, G2,0,kb,0,0];
l4 = [0,kb,0,G5,-kb-G5,0,0];
l5 = [0,0,0,0,1,kc*G6,-1];
l6 = [0,0,0,0,0,G6+G7,-G7];
l7 = [0,-G3,0,-G5,G3+G4+G5,-G7,G7];

A = [l1;l2;l3;l4;l5;l6;l7]
b = [Va;0;0; Id; 0;0 ;-Id]

Vs = A\b



printf("printing table with results...")

filename = "octvalues.tex";
f = fopen (filename, "w");
fprintf(f, "Ia & %e \\\\ \\hline \nIb & %e \\\\ \\hline \nIc & %e \\\\ \\hline \nId & %e \\\\ \\hline \n", Is, Id);
fprintf(f, "V1 & %e \\\\ \\hline \nV2 & %e \\\\ \\hline \nV3 & %e \\\\ \\hline \nV4 & %e \\\\ \\hline \nV5 & %e \\\\ \\hline \nV6 & 0.000000e+00 \\\\ \\hline \nV7 & %e \\\\ \\hline \nV8 & %e \\\\ \\hline \n", Vs)
fclose (f);

