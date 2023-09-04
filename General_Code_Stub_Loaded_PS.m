%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stub-Loaded Planar Wide-band Phase Shifter                 %
% General MATLAB Code for Deriving the ABCD and S-parameters %
% Programmed by Dr. Falih M. Alnahwi                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables:
%%%%%%%%%%%%
% Zo: The characteristic impedance of the reference line and the main line.
% Zp: The chacteristic impedance of the stubs.
% the1: The electrical length of the line sections (Theta 1)
% the2: The electrical length of the stub (Theta 2)
% N: Number of Stubs
% f_norm: Normalized frequency (f/fo) w.r.t the center frequency (fo)
% l_line_norm : Normalized length of the line w.r.t center wanelength.
% l_stub__norm : Normalized length of the stub w.r.t center wanelength.
% l_ref_norm : Normalized length of the reference line w.r.t center wanelength.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PART 1 %%
%%%%%%%%%%%%
% Generation of ABCD Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
N = input('Enter the number of stubs = ');
if N<2 | N~=round(N)
    error('Error!! N should be an integer larger than 1')
end

j = sqrt(-1);
syms Zo Zp the1 the2
Yo = 1/Zo;
Yp = 1/Zp;
T_line = [cos(the1) j*Zo*sin(the1);j*Yo*sin(the1) cos(the1)]; % T-matrix of the line
T_stub = [1 0;-j*Yp*cot(the2) 1];                             % T-matrix of the stub

T = T_stub;
for i = 2:N
    T = T*T_line*T_stub;
end
ABCD = expand(T);           % For multiplying the brackets
ABCD = simplify(ABCD)       % To simplify some trigonometric identities
A = ABCD(1,1); B = ABCD(1,2); C = ABCD(2,1); D = ABCD(2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S-Parameters of the Stub-Loaded Phase Sifter %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s11 = ((A-D)+(B/Zo-C*Zo))/((A+D)+B/Zo+C*Zo);
s21 = 2/((A+D)+B/Zo+C*Zo);
S = [s11 s21;s21 s11]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PART 2 %%
%%%%%%%%%%%%
% Substitution of Zo, Zp, the1, and the2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line_imp = input('The value of the line characteristic impedance is = ');
stub_imp = input('The value of the stub characteristic impedance is = ');

f_norm = 0:0.01:2;
l_line_norm = 0.25;
l_stub_norm = 0.25;

for i=1:length(f_norm)
    S11(i) = subs(s11,{Zo,Zp,the1,the2},{line_imp,stub_imp,2*pi*l_line_norm*f_norm(i),2*pi*l_stub_norm*f_norm(i)});
    S21(i) = subs(s21,{Zo,Zp,the1,the2},{line_imp,stub_imp,2*pi*l_line_norm*f_norm(i),2*pi*l_stub_norm*f_norm(i)});
end

S11_mag = 20*log10(abs(S11));
S21_mag = 20*log10(abs(S21));
S21_ang = phase(S21)*180/pi;   % Please use the command "phase" not "angle"

l_ref = (N-1)*l_line_norm+0.25;   % 0.25 to provide 90 deg. phase shift at fo.
the_ref = 2*pi*f_norm*l_ref;
ref = phase(exp(-j*the_ref))*180/pi;  % Phase of S21 of the reference line.

Phase = S21_ang-ref;    % Differential phase shift.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Demostration of the phase shifter parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(f_norm,S11_mag,'b','linewidth',2); grid on;
xlabel('Normalized frequency (f/fo)','fontsize',16);
ylabel('|S11| (dB)','fontsize',16);
set(gca,'fontsize',16);

figure(2)
plot(f_norm,S21_mag,'b','linewidth',2); grid on;
xlabel('Normalized frequency (f/f_o)','fontsize',16);
ylabel('|S21| (dB)','fontsize',16);
set(gca,'fontsize',16);

figure(3)
plot(f_norm,S21_ang,'b','linewidth',2);
hold on;
plot(f_norm,ref,'r--','linewidth',2);
hold off; grid on;
l = legend('Phase of main line','Phase of ref line')
l.FontSize = 14;
xlabel('Normalized frequency (f/fo)','fontsize',16);
ylabel('Phase (deg.)','fontsize',16);
set(gca,'fontsize',16);

figure(4)
plot(f_norm,Phase,'b','linewidth',2); grid on;
xlabel('Normalized frequency (f/fo)','fontsize',16);
ylabel('Differential Phase (deg.)','fontsize',16);
set(gca,'fontsize',16);
axis([0.4 1.6 85 95]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
