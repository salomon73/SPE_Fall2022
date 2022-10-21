function model=potentialwellmodel(obj,timestep,calcfull,mag,nblv)
% Computes the potential well at the given timestep and return the model to be able to
% interpolate either using grid coordinates or magnetic field line coordinates
% the potential well is calculated along each magnetic field
% line by taking the difference between a local potential well maximum
% and the maximum local minimum along the line on
% each side of the maximum

% model.z and model.r are the axial and radial positions where
% the potential depth has been calculated and is used with a
% scatteredinterpolant to calculate the well depth at the
% desired position
% model.pot is the potential well depth

% model.rathet is the magnetic field vector value at the
% corresponding position multiplied by the radial position
if iscell(timestep)
    timestep=cell2mat(timestep);
end
% Defines if only the well region is kept or if the hills are also kept
if nargin<3
    calcfull=false;
end
if nargin<4
    mag=[];
end
if nargin<5
    nblv=400;
end
% if one of the timesteps is 0 we take the external potential
id=find(timestep==0);
if(~isempty(id))
    timestep(id)=[];
end
Phi=-obj.pot(:,:,timestep);
if(~isempty(id))
    phiext=-obj.potxt(:,:,1);
    if isempty(obj.potxt)
        phiext=-obj.pot(:,:,1);
    end
    Phi=cat(3,Phi(:,:,1:id-1),phiext,Phi(:,:,id:end));
end

% We get the magnetic field lines rA_theta values
%lvls=sort(unique([obj.rAthet(:,1)',obj.rAthet(:,end)', obj.rAthet(1,:),obj.rAthet(end,:)]));
if isempty(mag)
Blines=obj.rAthet;
lvls=logspace(log10(min(Blines(obj.geomweight(:,:,1)>0))),log10(max(Blines(obj.geomweight(:,:,1)>0))),nblv);
%lvls=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),200);
Blines(obj.geomweight(:,:,1)<=0)=NaN;
lvls(isnan(lvls))=[];

% We obtain the magnetic field lines coordinates for each
% level in r and z
contpoints=contourc(obj.zgrid,obj.rgrid,Blines,lvls(1:end));
else
[z,r]=ndgrid(obj.zgrid,obj.rgrid);
locw=griddedInterpolant(z,r,obj.geomweight(:,:,1)','linear','none');
idrmin=find(obj.rgrid(1)<mag.r,1,'first');
idrmax=find(obj.rgrid(end)>mag.r,1,'last');
idzmin=find(obj.zgrid(1)<mag.z,1,'first');
idzmax=find(obj.zgrid(end)>mag.z,1,'last');

[z,r]=ndgrid(mag.z(idzmin:idzmax),mag.r(idrmin:idrmax));
geomwB=locw(z,r)'; 
geomwB(isnan(geomwB))=-1;
Blines=mag.Athet.*mag.r;
Blines=Blines((idrmin:idrmax),(idzmin:idzmax));
lvls=logspace(log10(min(Blines(geomwB>0))),log10(max(Blines(geomwB>0))),nblv);
Blines(geomwB<0)=NaN;
lvls(isnan(lvls))=[];

% We obtain the magnetic field lines coordinates for each
% level in r and z
contpoints=contourc(mag.z(idzmin:idzmax),mag.r(idrmin:idrmax),Blines,lvls(1:end));
    
end
[x,y,zcont]=C2xyz(contpoints);

% memory allocation for spped

potfinal=zeros(numel(cell2mat(x)),length(timestep));
Pot=cell(1,size(x,2)-1);
rathet=cell(1,size(x,2)-1);
[z,r]=ndgrid(obj.zgrid,obj.rgrid);

% for each timestep we calculate the well
for i=1:length(timestep)+~isempty(id)
    locPhi=Phi(:,:,i);
    % We delete points outside of the domain
    locPhi(obj.geomweight(:,:,1)<0)=NaN;
    
    
    locPhi=griddedInterpolant(z,r,locPhi','makima');
    % For each field line we remove the lowest maximum
    for j=1:size(x,2)
        %for i=1:size(pot,1)
        xloc=x{j};
        yloc=y{j};
        
        rathet{j}=zcont(j)*ones(1,length(xloc));
        
        if(length(xloc)>=3)% We neglect if the field length is too small
            
            Pot{j}=locPhi(xloc,yloc);
            Pot{j}=movmean(Pot{j},1);
            locpot=Pot{j};
            %locpot=locpot-min(locpot);
            % along the given field line j we calculate the
            % minimum on the left and right side of the
            % position k and calculate the local well depth
            for k=2:length(locpot)-1
                left=max(locpot(1:k-1));
                right=max(locpot(k+1:end));
                Pot{j}(k)=locpot(k)-min(left,right);
            end
            Pot{j}(1)=0;
            Pot{j}(end)=0;
        else
            Pot{j}=zeros(size(xloc));
        end
    end
    potfinal(:,i)=cell2mat(Pot);
end
%potfinal=cell2mat(potfinal);
if ~calcfull
    potfinal(potfinal>=0)=NaN;
end
model.z=cell2mat(x);
model.r=cell2mat(y);
model.pot=-potfinal;
model.rathet=cell2mat(rathet);
end