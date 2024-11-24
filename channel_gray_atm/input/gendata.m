% This is a matlab script that generates the input data
kwr=0; prec='real*4';

% Dimensions of grid
%nx=80;
nx=300;
ny=182;
nz=29;
% Nominal depth of model (meters)
% Size of domain
% Lx=200.e3;
% dy = Lx/nx;

%-- Close the domain with a wall @ Northern edge:
PSL=1.e+5; %- sea-level reference pressure
ps=ones(nx,ny)*PSL;
ps(:,ny)=0;

%-- make SST field:
c2K=273.15;
sst0=12+c2K; sst1=6; sst2=2;

 xx=([1:nx]-0.5)/nx*2*pi; 
 yy=([1:ny]-1)/(ny-1)*pi;
 ssty=sst1*ones(nx,1)*cos(yy);
 sstx=sst2*sin(xx')*ones(1,ny);
 sst=sst0+ssty+sstx;
 mnV=min(sst(:))-c2K; MxV=max(sst(:))-c2K;
fprintf(' SST min,max = %7.2f , %7.2f (oC)\n',mnV,MxV);

ps0=0; ps1=6; ps2=2; %- Surf. Pres Anom (in mb)
 psy=ps1*ones(nx,1)*cos(yy);
 psx=ps2*sin(xx')*ones(1,ny);
 %- convert to Pa and change sign:
 psIni=-100*(psy+psx);
%-------------------------------------------------------------------------
%-- Write binary files:

if kwr > 1,
 fNam='refSP.bin';
 fprintf('writing file: %s ...',fNam);
 fid=fopen(fNam,'w','b'); fwrite(fid,ps,prec); fclose(fid);
 fprintf(' done\n');
end

if kwr > 1,
 fNam='SST_ini.bin';
 fprintf('writing file: %s ...',fNam);
 fid=fopen(fNam,'w','b'); fwrite(fid,sst,prec); fclose(fid);
 fprintf(' done\n');
end

if kwr > 2,
 fNam='Eta_ini.bin';
 fprintf('writing file: %s ...',fNam);
 fid=fopen(fNam,'w','b'); fwrite(fid,psIni,prec); fclose(fid);
 fprintf(' done\n');
end


%-- Do some plots to check:
%kmFac=1.e-3; xax=x*kmFac;
% xBnd=[0 Lx]*kmFac; 

%return
figure(1);clf;
 xax=[1:nx]-0.5; yax=[1:ny]-0.5;
subplot(211);
 CI=12; clrb=[1 0 0]; thick=1; ccB=[0 0];
%CI=[15:0.5:20]*0.1; ccB=[1.5 2];
 var=ps;
 imagesc(xax,yax,var'); set(gca,'YDir','normal');
%contourf(xax,z,var',CI);
 if ccB(2)> ccB(1); caxis(ccB); end
 colorbar
%axis([xBnd yBnd]);
 grid

subplot(212);
 var=ps(1,:); yBnd=[0 ny];
 plot(yax,var,'r-');
 AA=axis; axis([yBnd AA(3:4)]);
 grid 

figure(2);clf;
%var=sst; ccB=[0 0]; titv='SST [K]';
 var=psIni; ccB=[0 0]; titv='PS-anom [Pa]';
 imagesc(xax,yax,var'); set(gca,'YDir','normal');
%contourf(xax,z,var',CI);
 if ccB(2)> ccB(1); caxis(ccB); end
 colorbar
%axis([xBnd yBnd]);
 grid
 title(titv);
