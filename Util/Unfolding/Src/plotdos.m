% plotdos.m
% Plots a DOS file written by unfold.F90
% SGM & JMS, Oct.2018   (mod. May 2019)

clear all
wdir = '../path_to_file/';
syslabel = 'SystemLabel';
fname = [wdir,syslabel,'.refoldedBands'];  % .path, .spin    

file = fopen(fname);                       % Read DOS files
dat = fscanf(file,'%d %d %f %f %f',5);
 nq = dat(1);            
 ne = dat(2);            
 emin = dat(3);
 emax = dat(4);
 Fermi = dat(5);         

q = zeros(3,nq);
dos = zeros(nq,ne);
label = string(zeros(nq,1));    
for iq = 1:nq
    q(:,iq) = fscanf(file,'%f %f %f',3);
    iline(iq) = fscanf(file,'%i',1);
    label(iq) = fscanf(file,'%10c',1);
    dos(iq,:) = fscanf(file,'%f',ne);
end
fclose(file);

% adjust labels to MATLAB style
for iq=1:nq               
   if (label(iq)=='  Gamma   ')
     label(iq)='\Gamma';
   elseif (label(iq)=='  '' ''     ')
     label(iq)=' ';
   else
     label(iq)=strtrim(label(iq));  
   end
end


for ix=1:3
  dq(ix,:) = diff(q(ix,:));
end
dq = sqrt(sum(dq.^2));
qpath = [0,cumsum(dq)];

emin = emin - Fermi;
emax = emax - Fermi;
de = (emax-emin)/(ne-1);
e = emin:de:emax;
[qpath,e] = ndgrid(qpath,e);

    
% convolution with a gaussian (smooth)
a = 0;
b = 0.03;
x = -5:0.1:5;
f = exp(-((x-a).^2)/(2*(b^2))) / (sqrt(2*pi*(b^2)));
f = f/trapz(f);
dos = conv2(f,f,dos,'same');
  
% --------- Plot of the density ---------     figure size:
figure(1)                                 %  [x0,   y0,   w,   h]  
set(gcf,'Units','Normalized','OuterPosition',[0.40, 0.05, 0.6, 0.7]);

PlotLog = 0;
zlim = max(dos(:));
zsat = zlim*1.0e-3;      % saturation for log scale

% emax = 1.5;            
% emin = -1.5;


if PlotLog == 0                            % normal      
  s = surf(qpath',e',dos');  
  s.EdgeColor = 'none';             
  axis([qpath(1),qpath(end),emin,emax,0,zlim])     
  caxis([0,zlim]);         
else                                        % log
  s = surf(qpath',e',max(dos,zsat)'); 
  set(gca,'colorscale','log');
  axis([qpath(1),qpath(end),emin,emax,zsat,zlim]) 
  caxis([zsat,zlim]);
end

s.EdgeColor = 'none';
set(gca,'FontSize',14);
ylabel('Energy  (eV)');
zlabel('LDOS');

view(2);
colormap(flipud(bone));  
h = colorbar; 
h.TickLength = 0.03;
set(get(h,'label'),'string',' LDOS (eV^{-1})');

% line edges and fermi level
hold on
cpl = [0.8,0.7,0.75];   
cfl = [0.92,0.59,0.59];              

lineEnd = find(iline(1:nq-1)<iline(2:nq));      
for jl = 1:numel(lineEnd)
    il = lineEnd(jl);
    plot3(qpath(il)*[1,1],[emin,emax],[zlim,zlim],'Color',cpl,'LineStyle',':')
end
plot3(qpath(end)*[1,1],[emin,emax],[zlim,zlim],'Color',cpl,'LineStyle',':')

plot3([qpath(1),qpath(end)],0*[1,1],[zlim,zlim],'Color',cfl,'LineStyle',':','LineWidth',0.5)

xticks( [qpath((lineEnd(1:1:numel(lineEnd)))),qpath(end)] );
xticklabels( [ label(lineEnd(1:1:numel(lineEnd))); label(nq) ] );
hold off  

% To specify labels:
% xticklabels({'\Gamma\color{gray}(000)','X','L'});

% ----------------------------------------------------------

