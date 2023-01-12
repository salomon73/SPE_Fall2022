%% add electrode to slice config  plot

R_up_grid = linspace(10,10.5,2);
Z_up_grid = linspace(-320,320,2);

[R,Z] = meshgrid(R_up_grid, Z_up_grid);
p = surf(Z,R,ones(2,2));
    set(p, 'FaceAlpha', 0.4)
    set(p, 'FaceColor', 'k')

R_down_grid = linspace(.5,1,2);
Z_down_grid = linspace(-320,320,2);

[R,Z] = meshgrid(R_down_grid, Z_down_grid);
p = surf(Z,R,ones(2,2));
    set(p, 'FaceAlpha', 0.6)
    set(p, 'FaceColor', 'k')
fig.Children(2).YLim = [.5 , 10.5]
