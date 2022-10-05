function createPartFile(filename,datas)
fileid=fopen(filename,'w');

fprintf(fileid,'&partsload\n');
fprintf(fileid,'mass=%19.14e,\n',datas.m);
fprintf(fileid,'charge=%19.14e,\n',datas.q);
fprintf(fileid,'weight=%19.14e,\n',datas.weight);
fprintf(fileid,'H0=%19.14e,\n',datas.H0);
fprintf(fileid,'P0=%19.14e,\n',datas.P0);
fprintf(fileid,'nblock=%d,\n',length(datas.ra));
fprintf(fileid,'npartsalloc=%d,\n',sum(datas.npartsslice));
fprintf(fileid,'velocitytype=%d,\n',datas.velocitytype);
fprintf(fileid,'radialtype=%d,\n',datas.radialtype);
fprintf(fileid,'temperature=%19.14e,\n',datas.temperature);

fprintf(fileid,'/\n\n');
fprintf(fileid,'//slices\n');
for i=1:length(datas.ra)
    fprintf(fileid,'%19.14e %19.14e %19.14e %d\n',datas.z(i),datas.ra(i),datas.rb(i),datas.npartsslice(i));
end
fprintf(fileid,'%19.14e \n',datas.z(end));
fclose(fileid);
end