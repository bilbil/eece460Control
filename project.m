%% EECE 460 Project 2013 April 10
% Bill Liu #63461081

%% Plant and system model creation setup

% Plant Model represented by 2nd order damped mass-spring system
% constants: m = 7*10^-4 kg, b= 2*10^-2 Ns/m, k = 32.2 N/m
Gnumerator = [1];
Gdenominator = [7*10^-4, 2*10^-2, 32.3];
G = tf(Gnumerator, Gdenominator);
G.InputName = 'u';
G.OutputName = 'p';
bode(G);

% Controller yet to be implemented
Cnumerator = [1];
Cdenominator = [1];
C_unimplemented = tf(Cnumerator, Cdenominator);
C_unimplemented.InputName = 'e2';
C_unimplemented.OutputName = 'u';

% r disturbance signal generation
s = tf('s');
% modified disturb generation gain to match 3 micrometer peak amplitude
r_disturb_generate = 0.986*3.5*10^-2/(s^2 + 75.4*s + 1.58*10^4);
r_disturb_generate.InputName = 'r-disturb-step';
r_disturb_generate.OutputName = 'r-disturb';

% rt step change filter for known rt signal unimplemented
rt_filter = tf(1,1);
rt_filter.InputName = 'rt-input';
rt_filter.OutputName = 'rt';
bode(rt_filter)

% n disturbance signal generation
n_disturb_generate = tf([3.18*10^-12, 0],[3.18*10^-5, 1]);
n_disturb_generate.InputName = 'n-disturb-step';
n_disturb_generate.OutputName = 'n';

% Feedback loop: e1 = r - p
feedback = sumblk('e1', 'r', 'p', '+-');

% Noise Adder: e2 = e1 + n
sumNoise = sumblk('e2', 'e1', 'n', '++');

% r disturbance adder: r = rt + r_disturb
sum_R = sumblk('r', 'rt', 'r-disturb', '++');

%connect the blocks together and see unimplemented closed loop
T_unimplemented = connect(C_unimplemented,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt'},{'e1'});
bode(T_unimplemented);

%% Specification 1 and 2 
%?DC zero steady state error and limit low frequency r-disturb gain to < -18dB for f < 200Hz 

% Implement affine parameterization for low frequency disturbance Trial 1
% natural frequency of closed loop
W_cl = 200*2*pi;
% damping of closed loop
Damp_cl = 2;

% prefilter
Fq = tf(W_cl^2,[1, 2*W_cl*Damp_cl, W_cl^2]);

% open plant model inverse
G_inverse = tf(Gdenominator, Gnumerator);

% affine Q
Q = Fq * G_inverse;
Q = minreal(Q);

% Convert Q to C
C = Q/(1 - Q*G);
C = minreal(C);
C.InputName = 'e2';
C.OutputName = 'u';

% transfer function from rt to error
T_rt_to_e1 = connect(C,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt'},{'e1'});
T_rt_to_e1 = minreal(T_rt_to_e1);
% bode(T_rt_to_e1);
[mag_C, phase_C] = bode(T_rt_to_e1,[200*2*pi])
%gain = 1.0308dB at 200Hz > -18dB: does not satisfy low frequency
%disturbance requirement

% Implement affine parameterization for low frequency disturbance Trial 2
% open plant model inverse
G_inverse = tf(Gdenominator, Gnumerator);

% prefilter trial 2
% shift natural frequency to higher value
W_cl_trial2 = 2500*2*pi;
Damp_cl_trial2 = 0.7;
s = tf('s');
Fq2 = W_cl_trial2^2/(s^2 + 2*W_cl_trial2*Damp_cl_trial2*s + W_cl_trial2^2);
Fq2 = minreal(Fq2);

% affine Q trial 2
Q2 = Fq2 * G_inverse;
Q2 = minreal(Q2);

% Convert Q to C
C2 = Q2/(1 - Q2*G);
C2 = minreal(C2);
C2.InputName = 'e2';
C2.OutputName = 'u';

% transfer function from rt to error
T_rt_to_e1_trial2 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'r'},{'e1'});
T_rt_to_e1_trial2 = minreal(T_rt_to_e1_trial2);

ltiview('bode',T_rt_to_e1_trial2);

% get rt to e1 gain at DC
[mag_C2,phase_C2]=bode(T_rt_to_e1_trial2,[0*2*pi])
% gain of 0 satisfy zero steady state tracking error at DC

% get rt to e1 gain at 200Hz
[mag_C2,phase_C2]=bode(T_rt_to_e1_trial2,[200*2*pi])
% gain of 0.1122 = 10*log(0.1122^2) dB = -19dB < -18dB: satisfy low
% frequency disturbance requirement

%% Specification 3 
% limit high frequency input r-disturb to <11dB at f>3000Hz

% see what the gain of transfer function from n to e1 is like
T_n_to_e1 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'n'},{'e1'});
T_n_to_e1 = minreal(T_n_to_e1);
ltiview('bode',T_n_to_e1);
[mag_C3,phase_C3]=bode(T_n_to_e1,[3000*2*pi])
% gain is 0.5758 = 10*log(0.5758^2)dB = -4.7946dB < 11dB : satisfy high
% frequency measurement noise specification

%% Specification 4 
% Check step response from rt-input to u in order to ensure maximum gain
% of 2/0.001 = 66.02dB, assuming maximum of 1mm step change is allowed for rt-input.

% transfer from unfiltered r to u
T_r_to_u_trial2 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'r'},{'u'});
bode(T_r_to_u_trial2);
% from the bode plot, we see maximum gain of 105dB which exceeds
% specification
ltiview('step',0.001*T_r_to_u_trial2);
% from the step response, we see maximum gain of 104.76dB, which agrees with
% bode plot since step input is high in frequency

% add a filter of 40Hz cutoff to known rt-input signal
s = tf('s');
rt_filter = (40*2*pi)^2/(s + 0.7071*40*2*pi*s + (40*2*pi)^2);
rt_filter.InputName = 'rt-input';
rt_filter.OutputName = 'rt';
bode(rt_filter)

%test gain from rt-input to u with the rt-input filter
T_rt_filtered_to_u_trial2 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt-input'},{'u'});
ltiview('step',0.001*T_rt_filtered_to_u_trial2);
% from the step response plot, we see a gain of 1.75/0.001 = 62.43dB :
% satisfy specification of <66.02dB

%% Specification 5 
% Check step response from rt to e1 or p to ensure < 5% overshoot.

% overshoot check from rt to e
T_rt_to_p_trial2 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt'},{'p'});
T_rt_to_p_trial2 = minreal(T_rt_to_e1_trial2);
% step change in rt
ltiview('step',T_rt_to_p_trial2);
% overshoot without rt filter = 4.6% < 5% : satisfy step input rt tracking requirement

% overshoot check from filtered rt to e
T_rt_filt_to_e1_trial2 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt-input'},{'p'});
T_rt_filt_to_e1_trial2 = minreal(T_rt_filt_to_e1_trial2);
% step change in rt
ltiview('step',T_rt_filt_to_e1_trial2);
% overshoot with rt filter = 2.43% < 5%: satisfy step input rt tracking requirement

%% Specification 6 
% Check response of system to rt-input step change of 1mm and disturbances
% by generation of disturbance signal from r_disturb_generate and
% n_disturb_generate to ensure error tracking within 370nm

[u_rt_step,t_rt_step] = gensig('square',0.5,10,0.00001);        %2Hz   r(t) step changes
% test system with rt step changes
T_system_test_with_rt_change = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt-input'},{'e1'});
lsim(T_system_test_with_rt_change,[0.001.*u_rt_step],[t_rt_step]);
% 32.1 micrometer maximum tracking error when performing 1mm step changes from filtered rt

% generate a low frequency input disturbance
T_r_disturb_to_e1_trial2 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'r-disturb-step'},{'e1'});
ltiview('step',T_r_disturb_to_e1_trial2);
% 16.4nm maximum tracking error

% generate a high frequency measurement disturbance
T_n_disturb_to_e1 = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'n-disturb-step'},{'e1'});
ltiview('step',T_n_disturb_to_e1);
% 20.3nm maximum tracking error

%generate periodic input rt signal and two disturbance signals rt and n as well
[u_r_disturb,t_r_disturb] = gensig('square',1/200,10,0.00001);   %200Hz r disturbance and apply to disturbance generators
[u_n_disturb,t_n_disturb] = gensig('square',1/3000,10,0.00001);  %3kHz  n disturbance and apply to disturbance generators
[u_rt_step,t_rt_step] = gensig('square',0.5,10,0.00001);        %2Hz   r(t) step changes

% test system with r and n disturbances
T_system_test = connect(C2,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt-input','r-disturb-step','n-disturb-step'},{'e1'});
lsim(T_system_test,[0.001.*u_rt_step,u_r_disturb,u_n_disturb],[t_rt_step]);
% displays adequate error tracking of 19.3nm in presences of both disturbances

%% Specification 7 - design for zero steady state tracking error in presence of 15Hz disturbances

%(e1(15*2?))/(r(15*2?))= |S_o (15*2?)|=0
% plot response to <r> = 10^-3 * sin(15*2?)
% Tracking error < 370nm also in presence of disturbance generated by 
% <r_disturb_generate> and <n_disturb_generate>.

% write transfer function using 's'
s = tf('s');

% get the magnitude and phase at 15Hz for G inverse
[mag_test,phase_test]=bode(G_inverse,[15*2*pi])
% 26.1502 magnitude, 4.1336 degrees phase

% add complex zeros at 15Hz to original Q design
% also add poles at 15Hz to cancel effect of complex zeros at
% other frequencies

newZero = (s^2 + (15*2*pi)^2)/((15*2*pi)^2);
bode(newZero)
hold on
newPole = ((15*2*pi)^2)/(s^2 + 0.000001*15*2*pi*s + (15*2*pi)^2);
bode(newPole)
newZeroAndPole = newZero*newPole;
bode(newZeroAndPole);
hold off
[mag_test,phase_test]=bode(newZeroAndPole,[15*2*pi])

Qa = newZeroAndPole*Fq2*G_inverse;
bode(Qa,'y');
hold on;
bode(Q2,'g');
hold off;
[mag_test,phase_test]=bode(Qa,[15*2*pi])

% make a bandpass filter at 15Hz
% bandpass filter needs to have relative degree of 2 due to G_inverse
bp = ((0.1*s)/(s^2 + 0.1*s + (15*2*pi)^2))^2;
bode(bp);
[mag_test,phase_test]=bode(bp*G_inverse*G,[15*2*pi])

% create Qb to get perfect G inverse at 15Hz and apply bandpass filter
Qb = G_inverse*bp;
bode(Qb);
hold on;
bode(G_inverse)
hold on;
hold off;
[mag_test,phase_test]=bode(Qb,[15*2*pi])
% from plot wee se both curves cross intersect at 15Hz for gain and phase

% final Q design as combination of Qa and Qb
testQ3_filt = Qa + Qb;
bode(testQ3_filt,'r');
[mag_test,phase_test]=bode(testQ3_filt,[15*2*pi])
hold on
bode(Q2)

% Convert Q to C
testC3 = testQ3_filt/(1 - testQ3_filt*G);
testC3 = minreal(testC3);
testC3.InputName = 'e2';
testC3.OutputName = 'u';

% transfer function from r to error
T_test = connect(testC3,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'r'},{'e1'});
T_test = minreal(T_test);
bode(T_test, 'y')
hold on;
bode(T_rt_to_e1_trial2,'g');

% check gain requirements
[mag_test,phase_test]=bode(T_test,[200*2*pi])
% gain = 0.1122, same as Part A, satisfy requirement for f < 200Hz
[mag_test,phase_test]=bode(T_test,[15*2*pi])
% 2.919*10^-4 gain at 15Hz
[mag_test,phase_test]=bode(T_test,[0*2*pi])
% 0 gain and 90 degrees phase at DC

% transfer function from noise to error
T_test_noise = connect(testC3,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'n'},{'e1'});
T_test_noise=minreal(T_test_noise);
bode(T_test_noise, 'y')
[mag_test,phase_test]=bode(T_test_noise,[3000*2*pi])
% low pass filter and gain = 0.5758 and 75.3236 degrees at 3000Hz, satisfy < 3.7 gain
% requirement for f > 3000Hz

%generate periodic input rt signal and two disturbance signals rt and n as well
[u_r_disturb,t_r_disturb] = gensig('sin',1/15,8,0.00001);   %15Hz r disturbance and apply to disturbance generators
[u_n_disturb,t_n_disturb] = gensig('square',1/3000,8,0.00001);  %3kHz  n disturbance and apply to disturbance generators
[u_rt_step,t_rt_step] = gensig('square',1/3,8,0.00001);        %3Hz   r(t) step changes
[u_n_sine,t_n_sine] = gensig('sin',1/3000,8,0.00001);  %3kHz  n disturbance and apply to disturbance generators

% test system with 15Hz eccentric r disturbance
T_system_test_with_r_disturb_trial3 = connect(testC3,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'r'},{'e1'});
lsim(T_system_test_with_r_disturb_trial3,[u_r_disturb],[t_r_disturb]);

% test system with noise sine wave
T_system_test_with_n_sine_trial3 = connect(testC3,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt','n'},{'e1'});
lsim(T_system_test_with_n_sine_trial3,[0.0001.*u_rt_step,u_n_sine],[t_rt_step]);

% test system with n disturbance
T_system_test_with_n_disturb_trial3 = connect(testC3,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt','n-disturb-step'},{'e1'});
lsim(T_system_test_with_n_disturb_trial3,[0.0001.*u_rt_step,u_n_disturb],[t_rt_step]);
ltiview(T_system_test_with_n_disturb_trial3)

% test system with r and n disturbances
T_system_test_trial3 = connect(testC3,G,feedback,sumNoise,sum_R,r_disturb_generate,n_disturb_generate,rt_filter,...
{'rt-input','r','n-disturb-step'},{'e1'});
lsim(T_system_test_trial3,[0.001.*u_rt_step,0*u_r_disturb,u_n_disturb],[t_rt_step]);
