%%
timesteppart=length(M.tpart);

fig=figure;
ax1=gca;

timestepN=find(M.tpart(end)==M.t2d,1);
Ndistrib=M.N(:,:,timestepN);
h=contourf(ax1,M.zgrid,M.rgrid,Ndistrib);
%set(h, 'EdgeColor', 'none');
hold on
[r,z]=find(Ndistrib~=0);
xlim(ax1,[M.zgrid(min(z)) M.zgrid(max(z))])
ylim(ax1,[M.rgrid(min(r)) M.rgrid(max(r))])
xlabel(ax1,'Z [m]')
ylabel(ax1,'R [m]')
c = colorbar(ax1);
c.Label.String= 'n [m^{-3}]';
view(ax1,2)
%set(ax1,'colorscale','log')


[x,y]=ginput(1);
Zindex=find(x>M.zgrid,1,'last');
Rindex=find(y>M.rgrid,1,'last');
%  Rindex=155;
%  Zindex=344;

% Rindex=123;
% Zindex=658;
gcs=false;

f=figure();
ax1=subplot(1,2,1);
ax2=subplot(1,2,2);
%[p,maxnb,c]=M.displayPhaseSpace(1:2,Rindex,Zindex,'',sprintf('t=%1.3g [s]',M.tpart(1)),ax1);
% ax1=subplot(1,2,1);
% hold(ax1,'off')
c=cell(2,1);
c{1}=linspace(-2e7,2e7,21);
c{2}=linspace(-0e6,3e7,51);


[p,maxnb,c]=M.displayPhaseSpace('parper',2,Rindex,Zindex,'',sprintf('t=%1.3g [s]',M.tpart(2)),ax1,-1,c,gcs);
M.displayPhaseSpace('parper',length(M.tpart)+(0),Rindex,Zindex,'',sprintf('t=%1.3g [s]',M.tpart(end)),ax2,-1,c,gcs);
xlimits=xlim(ax1);
xlimits=[min([xlimits,xlim(ax2)]) max([xlimits,xlim(ax2)])];
ylimits=ylim(ax1);
ylimits=[min([ylimits,ylim(ax2)]) max([ylimits,ylim(ax2)])];
xlim([ax1,ax2],xlimits)
ylim([ax1,ax2],ylimits)
climitsmax=max([caxis(ax1),caxis(ax2)]);
caxis(ax1,[0 climitsmax]);
caxis(ax2,[0 climitsmax]);
xtickformat(ax1,'%.3g')
ytickformat(ax1,'%.3g')
ax1.YAxis.Exponent = 0;
ax1.XAxis.Exponent = 0;
legend(ax2,'off');
ax1.Children(5).LineStyle='none';

xtickformat(ax2,'%.3g')
ytickformat(ax2,'%.3g')
ax2.YAxis.Exponent = 0;
ax2.XAxis.Exponent = 0;

sgtitle(sprintf('r=%1.2f [mm] z=%1.2f [mm] \\Delta\\phi=%1.1f[kV] R=%1.1f',M.rgrid(Rindex)*1e3,M.zgrid(Zindex)*1e3,(M.potout-M.potinn)*M.phinorm/1e3,M.Rcurv))
M.savegraph(f,sprintf('%s/%s_phasespaceR%dZ%dbegendnogcs',M.folder,M.name,Rindex,Zindex),[16,10]);