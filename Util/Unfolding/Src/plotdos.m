% plotdos.m
% Plots unfolded and refolded qLDOS written by unfold.F90.
% SGM & JMS, Oct.2018

clear
wdir = './path_to_file/';   
syslabel = 'systemLabel';           
fnameA = [wdir,syslabel,'.unfoldedBands'];  % .path, .spin    
fnameB = [wdir,syslabel,'.refoldedBands'];     
 
fileA = fopen(fnameA);                          % Read DOS
datA = fscanf(fileA,'%d %d %f %f %f',5);
nq = datA(1);
ne = datA(2);
emin = datA(3);
emax = datA(4);
efermi = dat(5);
fileB = fopen(fnameB);
datB = fscanf(fileB,'%d %d %f %f %f',5);
q = zeros(3,nq);
dosA = zeros(nq,ne);
dosB = zeros(nq,ne);
ST = zeros(nq,1);
label = string(ST);
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

% --------- Plot of the density ---------
% if only Gamma point
if nq==1       
  subplot(2,1,1)
  plot(e,dosA)
  xlabel('Energy')

  subplot(2,1,2)
  plot(e,dosB)
  xlabel('Energy')
  ylabel('DOS')
  grid on
  
else
      
  % --- PLOTS --- set for normal and log scales. If you wish to
  % use NORMAL scale set parameter PlotLog == 0. If you wish to
  % use LOG scale set PlotLog == 1.
  
  % To specify labels different from those on systemLabel.fdf: 
  % xticklabels({'\Gamma\color{gray}(000)','X','L'});
    
  PlotLogA = 1;
  PlotLogB = 1;
  
 % -----------------------------------------------------------
  
 % size of figure window:  [x0,y0,w,h]  
 set(gcf,'Units','Normalized','OuterPosition',[0.535, -0.002, 0.47, 1]);

  % --- Subplot number 1 (A) ---
  subplot(2,1,1)
  
  if PlotLogA == 0                         % normal scale        
      sA = surf(qpath',e',dosA');  
      sA.EdgeColor = 'none';
      zlimit = max(max(dosA));               
      axis([0,qpath(end),emin,emax,0,zlimit])     
      caxis([0,zlimit]);         
  else                                     % logaritmic scale 
      sA = surf(qpath',e',max(dosA,0.001)'); 
      sA.EdgeColor = 'none';
      zlimit = max(max(dosA));
      set(gca,'colorscale','log');
      axis([0,qpath(end),emin,emax,0,zlimit]) 
      caxis([0.001,zlimit]);  %caxis('auto');
  end
    
  set(gca,'FontSize',12);
  ylabel('Energy  (eV)');
  zlabel('LDOS');

  view(2);
  colormap(flipud(bone));
  colorbar;
  h = colorbar; 
  h.TickLength = 0.05;
  set(get(h,'label'),'string',' LDOS (eV^{-1})'); % 'title': over cbar
                                                  % 'label': side of cb

  lineEnd = find(iline(1:nq-1)<iline(2:nq));      % Line edges  
  hold on
  cpl = [0.635,0.078,0.184];            % burgundy: [0.635,0.078,0.184];
                                        % orange:   [0.850,0.325,0.098];
  cfl = [0.92,0.59,0.59];               % salmon:   [0.92,0.59,0.56];
                                        % gold:     [0.98,0.86,0.41];
  for jl = 1:numel(lineEnd)
      il = lineEnd(jl);
      plot3(qpath(il)*[1,1],[emin,emax],[zlimit,zlimit],'Color',cpl,'LineStyle',':')
  end
plot3(qpath(end)*[1,1],[emin,emax],[zlimit,zlimit],'Color',cpl,'LineStyle',':')
plot3([qpath(1),qpath(end)],0*[1,1],[zlimit,zlimit],'Color',cfl,'LineStyle','-.','LineWidth',0.5)

  vtick = [qpath((lineEnd(1:1:numel(lineEnd)))),qpath(end)];
  xticks(vtick);
  listoflabels = [ label(lineEnd(1:1:numel(lineEnd))); label(nq) ];
  xticklabels( listoflabels );
  hold off  
  
  
  % --- Subplot number 2 (B) ---
  subplot(2,1,2)
  
  if PlotLogB == 0                         % normal scale        
      sB = surf(qpath',e',dosB'); 
      sB.EdgeColor = 'none';
      zlimit = max(max(dosB));               
      axis([0,qpath(end),emin,emax,0,zlimit])     
      caxis([0,zlimit]);      
  else                                     % log scale 
      sB = surf(qpath',e',max(dosB,0.001)'); 
      sB.EdgeColor = 'none';
      zlimit = max(max(dosB));
      set(gca,'colorscale','log');
      axis([0,qpath(end),emin,emax,0,zlimit]) 
      caxis([0.001,zlimit]);  %caxis('auto');
  end
    
  set(gca,'FontSize',12);
  ylabel('Energy  (eV)');
  zlabel('LDOS');

  view(2);
  colormap(flipud(bone));
  colorbar;
  h = colorbar; 
  h.TickLength = 0.05;
  set(get(h,'label'),'string',' LDOS (eV^{-1})'); % 'title': over cbar
                                                  % 'label': side of cb
                                                  
  lineEnd = find(iline(1:nq-1)<iline(2:nq));      % Line edges  
  hold on
  cpl = [0.635,0.078,0.184];            % burgundy: [0.635,0.078,0.184];
                                        % orange:   [0.850,0.325,0.098];
  cfl = [0.92,0.59,0.59];               % salmon:   [0.92,0.59,0.56];
                                        % gold:     [0.98,0.86,0.41];
  for jl = 1:numel(lineEnd)       
      il = lineEnd(jl);
      plot3(qpath(il)*[1,1],[emin,emax],[zlimit,zlimit],'Color',cpl,'LineStyle',':')
  end
  plot3(qpath(end)*[1,1],[emin,emax],[zlimit,zlimit],'Color',cpl,'LineStyle',':') % last sym point
  plot3([qpath(1),qpath(end)],0*[1,1],[zlimit,zlimit],'Color',cfl,'LineStyle','-.','LineWidth',0.5)

  vtick = [qpath((lineEnd(1:1:numel(lineEnd)))),qpath(end)];
  xticks(vtick);
  listoflabels = [ label(lineEnd(1:1:numel(lineEnd))); label(nq) ];
  xticklabels( listoflabels ); 
  hold on

end

hold off

% -------------------------------------------------------------------------
