%%%%%%%%% Analise Teorica %%%%%%%%%%%%%%%
close all
clear all

pkg load symbolic
format long

V = 12.;
Vin = 230.;
V1 = 13.0855;
n_espiras = Vin/V1;
nr_diodos = 26; 
R_vr = 13600; %resistor of voltage regulator
C = 8*10^-6; 


f= 50;
T = 1/f;
w = 2*pi*f;

t=linspace(0,1/(2*f),10000);

%%%%%%%%% Bridge Rectifier %%%%%%%%%%%%%%
Vs = V1*cos(w*t);
VS = abs(Vs);



%%%%%%%%% Voltage Rectifier - Parte 1 %%%%%%%%%%%%%%
eta = 1;
V_D = V/nr_diodos;

I_S = 1e-14;
V_T = 0.02335; 
r_d = eta*V_T/(I_S*exp(V_D/(eta*V_T)));
R_d = nr_diodos*r_d;
R_aux = R_d/(R_d + R_vr)
R_total = R_d + R_vr

%V0 = R_aux*V1*cos(w*t);


%%%%%%%%% Envelope Detector %%%%%%%%%%%%%%
vOhr = zeros(1, length(t));
vO = zeros(1, length(t));

tOFF = 1/w * atan(1/w/R_total/C);
vOnexp = V1*cos(w*tOFF)*exp(-(t-tOFF)/R_total/C);


figure
for i=1:length(t)
  if (Vs(i) > 0)
    vOhr(i) = Vs(i);
  else
    vOhr(i) = -Vs(i);
  endif
endfor

plot(t*1000, vOhr)
hold

for i=1:length(t)
  if t(i) < tOFF
    vO(i) = abs(Vs(i));
  elseif vOnexp(i) > vOhr(i)
    vO(i) = vOnexp(i);
  else 
    vO(i) = abs(Vs(i));
  endif
endfor

plot(t*1000, vO, "linewidth", 1)
title("Envelope Voltage")
xlabel ("t[ms]")
%legend("rectified","envelope")
%print ("venvelope.eps", "-depsc");
print -djpg venvelope.jpg

%%%%%%%%% Ripple %%%%%%%%%%%%%%
t0=linspace(0,1/(2*f),20);

Vs = V1*cos(w*t0);
VS = abs(Vs);

vOhr = zeros(1, length(t0));
vO = zeros(1, length(t0));

tOFF = 1/w * atan(1/w/R_total/C);
vOnexp = V1*cos(w*tOFF)*exp(-(t0-tOFF)/R_total/C);

for i=1:length(t0)
  if (Vs(i) > 0)
    vOhr(i) = Vs(i);
  else
    vOhr(i) = -Vs(i);
  endif
endfor

plot(t0*1000, vOhr)
hold

for i=1:length(t0)
  if t0(i) < tOFF
    vO(i) = abs(Vs(i));
  elseif vOnexp(i) > vOhr(i)
    vO(i) = vOnexp(i);
  else 
    vO(i) = abs(Vs(i));
  endif
endfor
num = 5;

t1=linspace(0,num*T/2,100);

hold off;
count = 1;

for i=1:num
  for j=1:length(t0)
  vripple(count) = (vO(j) - V1);
  count = count+1;
  endfor
endfor

plot(t1*1000, vripple)
title("Ripple voltage v(t)")
xlabel ("t[ms]")
%print ("vripple.eps", "-depsc");
print -djpg vripple.jpg

%ripple 
ripple = abs(max(vripple) - min(vripple))

%%%%%%%%% Voltage Rectifier - Parte 2 %%%%%%%%%%%%%%
count = 1;
for i=1:num
  for j=1:length(t0)
  vC(count) = (vO(j))*R_aux;
  count = count+1;
  endfor
endfor

v_0 = mean(vC)

plot(t1*1000, vC)
title("Voltage Rectifier v(t)")
xlabel ("t[ms]")
%print ("voltage_rectifier.eps", "-depsc");
print -djpg voltage_rectifier.jpg

%%%%%%%%%% Print Table with data %%%%%%%%%

filename = "data.tex";
f = fopen (filename, "w");
fprintf(f, "{\\bf Quantaties} & {\\bf Volts} \\\\ \\hline \naverage & %e \\\\ \\hline \nripple & %e \\\\ \\hline \n deviation & %e \\\\ \\hline \n", v_0, ripple, v_0-12);
fclose (f);

