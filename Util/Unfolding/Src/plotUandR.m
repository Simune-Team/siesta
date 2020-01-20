% plotUandR.m
% Plots fully-unfolded and refolded qLDOS written by unfold.F90
% SGM & JMS, Oct.2018   (mod. May 2019)

clear all
wdir = './path_to_file/';   
syslabel = 'systemLabel';           
fnameA = [wdir,syslabel,'.unfoldedBands'];    % .path, .spin    
fnameB = [wdir,syslabel,'.refoldedBands'];     
 
fileA = fopen(fnameA);                            % Read DOS
datA = fscanf(fileA,'%d %d %f %f %f',5);
 nq = datA(1);
 ne = datA(2);
 emin = datA(3);
 emax = datA(4);
 efermi = datA(5);
fileB = fopen(fnameB);
datB = fscanf(fileB,'%d %d %f %f %f',5);

q = zeros(3,nq);
dosA = zeros(nq,ne);
dosB = zeros(nq,ne);
label = string(zeros(nq,1));
for iq = 1:nq
    q(:,iq) = fscanf(fileA,'%f %f %f',3);
    iline(iq) = fscanf(fileA,'%i',1);
    label(iq) = fscanf(fileA,'%10c',1);
    dosA(iq,:) = fscanf(fileA,'%f',ne);
    fscanf(fileB,'%f %f %f',3);
    fscanf(fileB,'%i',1);
    fscanf(fileB,'%10c',1);
    dosB(iq,:) = fscanf(fileB,'%f',ne);
end
fclose(fileA);  
fclose(fileB);

% set zero at Fermi level
emin = emin - efermi;
emax = emax - efermi;

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

if nq==1
    qpath = 0;
else
    for ix=1:3
        dq(ix,:) = diff(q(ix,:));
    end
    dq = sqrt(sum(dq.^2));
    qpath = [0,cumsum(dq)];
end
de = (emax-emin)/(ne-1);
e = emin:de:emax;
[qpath,e] = ndgrid(qpath,e);

% --------- Plot of the density ---------     figure size:
%                                             [x0,   y0,   w,    h]  
 set(gcf,'Units','Normalized','OuterPosition',[0.55, 0.05, 0.45, 0.9]);

if nq==1            % gamma point only
  subplot(2,1,1)
  plot(e,dosA)
  xlabel('Energy'); ylabel('DOS');

  subplot(2,1,2)
  plot(e,dosB)
  xlabel('Energy');  ylabel('DOS');
  grid on

else
      
%     NORMAL scale:  PlotLog == 0
%        LOG scale:  PlotLog == 1
  
  PlotLogA = 1;
  PlotLogB = 0;
  
  % --- Full unfolding ---
  subplot(2,1,1)
  
  zlim = max(dosA(:));
  zsat = zlim*1.0e-3;    % saturation for log scale
  
  %emax =  10.0;           
  %emin = -15.0;
  
  if PlotLogA == 0                            % normal      
      sA = surf(qpath',e',dosA');  
      axis([qpath(1),qpath(end),emin,emax,0,zlim])     
      caxis([0,zlim]);         
  else                                          % log
      sA = surf(qpath',e',max(dosA,zsat)'); 
      set(gca,'colorscale','log');
      axis([qpath(1),qpath(end),emin,emax,zsat,zlim]) 
      caxis([zsat,zlim]);   %caxis('auto');
  end
    
  sA.EdgeColor = 'none';               
  set(gca,'FontSize',13);
  ylabel('Energy  (eV)');
  zlabel('LDOS');

  view(2);
  colormap(flipud(bone));
  h = colorbar; 
  h.TickLength = 0.05;
  set(get(h,'label'),'string',' LDOS (eV^{-1})');
                                                 
  % line edges and fermi level
  hold on
  cpl = [0.8,0.7,0.75];  cfl = [0.92,0.59,0.59];     
  
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
  

  % --- Refolding ---
  subplot(2,1,2)
  
  zlim = max(dosB(:));
  zsat = zlim*1.0e-2;      % saturation for log scale
  
  %emax =  10.0;           
  %emin = -15.0;
    
  if PlotLogB == 0                           % normal       
      sB = surf(qpath',e',dosB'); 
      axis([qpath(1),qpath(end),emin,emax,0,zlim])     
      caxis([0,zlim]);      
  else                                         % log
      sB = surf(qpath',e',max(dosB,zsat)'); 
      set(gca,'colorscale','log');
      axis([qpath(1),qpath(end),emin,emax,zsat,zlim]) 
      caxis([zsat,zlim]);  %caxis('auto');
  end
    
  sB.EdgeColor = 'none';
  set(gca,'FontSize',13);
  ylabel('Energy  (eV)');
  zlabel('LDOS');

  view(2);
  colormap(flipud(bone));
  h = colorbar; 
  h.TickLength = 0.05;
  set(get(h,'label'),'string',' LDOS (eV^{-1})');
                                                  
% line edges and fermi level
  hold on
  cpl = [0.8,0.7,0.75];  cfl = [0.92,0.59,0.59];     
  
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
  
end

% --------------------------------------------
