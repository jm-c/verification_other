
dz=25; dx=1000;
nx=60; nz=10; nDead=5;
nr=nz+nDead;

Lx=nx*dx;
x=[1:nx]*dx; x=x-mean(x); %x=x+Lx/2;

%- initial free-surface: single mode
Ampli=0.25;
 var=2*x*pi/Lx; et1=Ampli*(cos(var)+1)*0.5;
%et(1)=0; et(nx)=0;

%- initial free-surface: exponential shape (width: 2*nyD*dy) height: Ampli)
 nyD=5;
 var=exp(-(x/dx/nyD).^2);   % Exponential width: 2*nyD*dy
 et2=Ampli*var;

%- plot to check:

figure(1);clf;
%subplot(211);
 var=et1;
 plot(x/dx,var,'b-');
 hold on;
 var=et2;
 plot(x/dx,var,'r-');
 hold off;
 grid
%axis([9.4 11 -300 0]);
%title('t\_ini');

hBot=nz*dz;
zTop=nDead*dz;
zDeep=nr*dz;
fprintf('dz= %3.0f , hBot= %4.0f , zTop= %4.0f , zDeep= %4.0f , zDeep-zTop= %4.0f\n', ...
         dz,hBot, zTop, zDeep, zDeep-zTop )

%return

var=-hBot*ones(nx,1);
fid=fopen('flat_bottom.bin','w','b'); fwrite(fid,var,'real*8'); fclose(fid);

var=-zTop*ones(nx,1);
fid=fopen('flat_top.bin','w','b'); fwrite(fid,var,'real*8'); fclose(fid);

var=-zDeep*ones(nx,1);
fid=fopen('flat_deep.bin','w','b'); fwrite(fid,var,'real*8'); fclose(fid);

fid=fopen('Eta_ini.bin','w','b'); fwrite(fid,et2,'real*8'); fclose(fid);

return
