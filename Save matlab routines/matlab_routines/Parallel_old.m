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

N=[1,2,4];
t=[182.6, 113.2,64.4];
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