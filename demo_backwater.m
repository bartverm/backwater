restoredefaultpath
clearvars
close all

B=Backwater;

%% Make A2 curve
figure
subplot(3,1,2)
B.x_end=-4e4;
B.So=-1e-4;
B.a0=B.a_critical+1;
plot(B)

%% Make A3 curve
subplot(3,1,3)
B.a0=B.a_critical-1;
B.x_end=B.x_target;
B.bed_offset=1;
plot(B)


%% Make H2 curve
figure
subplot(3,1,2)
B.So=0;
B.a0=4;
B.x_end=B.x_target;
B.bed_offset=1;
plot(B)

%% Make H3 curve
subplot(3,1,3)
B.a0=1;
B.x_end=B.x_target;
plot(B)


%% Make M1 curve
figure
subplot(3,1,1)
B.So=1e-4;
B.a0=B.a_equilibrium+2;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)

%% Make M2 curve
subplot(3,1,2)
B.So=1e-4;
B.a0=B.a_equilibrium-2;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)

%% Make M3 curve
subplot(3,1,3)
B.So=1e-4;
B.a0=B.a_critical-1;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)


%% Make C1 curve
figure
subplot(3,1,1)
B.So=B.Sc;
B.a0=B.a_equilibrium+1;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)

%% Make C2 curve
subplot(3,1,2)
B.So=B.Sc;
B.a0=B.a_equilibrium;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)
text(300,0.5,'unstable flow','HorizontalAlignment','center')

%% Make C3 curve
subplot(3,1,3)
B.So=B.Sc;
B.a0=B.a_equilibrium-1;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)




%% Make S1 curve
figure
subplot(3,1,1)
B.So=1e-2;
B.a0=B.a_critical+1;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)

%% Make S2 curve
subplot(3,1,2)
B.So=1e-2;
B.a0=(B.a_critical+B.a_equilibrium)/2;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)

%% Make S3 curve
subplot(3,1,3)
B.So=1e-2;
B.a0=1;
B.bed_offset=0;
B.x_end=B.x_target;
plot(B)




