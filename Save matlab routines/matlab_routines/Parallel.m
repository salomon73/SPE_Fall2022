clear all
N=[1,2,4,8];
t=[469.8, 294.1,253.2, 158.5];

Speedup=t(1)./t;
Efficiency=Speedup./N*100;
figure
plot(N,Speedup)
xlabel('N')
ylabel('Speed up')

figure
plot(N,Efficiency)
xlabel('N')
ylabel('Efficiency %')

N=[    1,   2,     4,     8,      12,      16,    24];
t=[617.1, 389.4, 287.2, 211.9,   182.7, 142.4, 132.7];
Speedup=t(1)./t;
Efficiency=Speedup./N*100;
figure
plot(N,Speedup)
hold on
plot(N,N)
xlabel('N')
ylabel('Speed up')

figure
plot(N,Efficiency)
xlabel('N')
ylabel('Efficiency %')