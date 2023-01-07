

x = linspace(-300,300,1000);
y = linspace(1,10,100);
[X,Y] = meshgrid(x,y);
Dens = NaN*X+NaN*Y;

Slice_1 = find((-300 <= x) & (x<=-100));
Slice_2 = find((-100 <= x) & (x<=100));
Slice_3 = find(( 100 <= x) & (x<=300));

Slice1y = find((2.9 <=y)& (y<=3.1));
Slice2y = find((4.9 <=y)& (y<=5.1));
Slice3y = find((7.9 <=y)& (y<=8.1));

Dens(Slice1y, Slice_1) = 1000/length(Slice_1);
Dens(Slice2y, Slice_2) = 1000/length(Slice_2);
Dens(Slice3y, Slice_3) = 1000/length(Slice_3);
%%
figure
hold on 
pcolor(X,Y,Dens)
colorbar
shading interp
xlimit([-320 320])