function out = Copper()
E0 = [.001 .00125 .00150 .00175 .002 .00225 .00250 .00275 .003 .00350 ...
      .004 .00450 .005 .00550 .006 .00650 .007 .00750 .008 .00850  ...
      .009 .00950];

E1 = [.01 .0125 .015 .0175 .02 .0225 .0250 .0275 .03 .0350 .04 .0450 .05 ...
      .0550 .06 .0650 .07 .0750 .08 .0850  .09 .095 ];

E2 = [.1 .125 .150 .175 .2 .3 .4 .5 .6 .7 .8 .9 ];

E3 = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5];

E4 = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95];


Eloss0 = [281.23 314.31 344.21 371.71 397.31 421.35 444.08 465.71 486.37 525.26 ...
          561.46 595.46 627.61 662.76 696.43 728.76 759.86 789.83 818.74 846.67 ...
          873.66 899.76 ];

Eloss1 = [925.03 1039.9 1138.2 1222.0 1293.1 1347.6 1402.2 1447.8 1488.8  1560.2 1620.3 1671.5...
          1715.5 1753.6 1786.6 1815.2 1840.0 1861.5 1880.1 1896.1 1909.7 1921.3 ];

Eloss2 = [1931.0 1957.4 1957.6 1941.5 1915.1 1767.5 1612.3 1474.5 1357.4 1259.2 1180.1 1114.5];

Eloss3 = [1058.9 855.87 719.47 626.46 557.63 504.76 461.60 426.43 396.88 ...
          371.66 349.82 330.31 313.86 298.80 285.27 273.08 262.10 252.06 ];

Eloss4 = [242.84 180.32 145.38 122.81 106.95 95.132 85.965 78.632 72.614  ...
          67.596 63.335 59.610 56.481 53.640 51.200 48.988 47.001 45.207 ];

    out.E       = cat(2, E0, E1, E2, E3, E4);
    out.Eloss   = cat(2, Eloss0, Eloss1, Eloss2, Eloss3, Eloss4);
    out.element = 'Cu';
    out.EFermi  = 7.0; %eV

end