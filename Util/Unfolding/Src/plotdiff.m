% plotdiff.m
% Plots the difference between qDOS of two systems written by unfold 
% in order to highlight the differences between them.
% It can be useful for obtaining changes in the bands of a defective 
% system with respect to the crystaline system.
% In the plot, new states arising from the defects are depicted in
% red, while states that vanish appear in blue.
% Both systems must have identical 'UnfoldedBandLines' blocks in 
% their fdf files (that is, same sampling energies and q points). 
% SGM & JMS, Oct.2018

clear
wdir = './WU_OUTPUTS/';    
syslabelA = 'systemA';              % System A : non-defective
syslabelB = 'systemB';              % System B : defective
fnameA = [wdir,syslabelA,'.refoldedBands'];          % '.refoldedBands' or
fnameB = [wdir,syslabelB,'.refoldedBands']; %.path1  %   '.unfoldedBands'.

% Read DOS
fileA = fopen(fnameA);
datA = fscanf(fileA,'%d %d %f %f %f',5);
nq = datA(1);
ne = datA(2);
emin = datA(3);
emax = datA(4);
efermiA = datA(5);
fileB = fopen(fnameB);
datB = fscanf(fileB,'%d %d %f %f %f',5); % nq,ne,emin,emax are read from fileA
efermiB = datB(5);
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

% set zero at efermiA (non-defective)
emin = emin - efermiA;
emax = emax - efermiA;

for iq=1:nq                % Adjust labels to MATLAB style
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

% Align spectra if necessary (disp of an integer number of pixels)

plim = 15;
dmin = 333333;
for ip = -plim:plim
   ipmin = 1+plim;  ipmax = ne-plim;
   pdiff = sum(sum( abs((dosA(:,ipmin:ipmax)-dosB(:,ipmin+ip:ipmax+ip))) ));
   if pdiff < dmin
       dmin = pdiff;
       ipopt = ip;
   end
end
if (ipopt>0)
    dosB(:,1+ipopt:ne-ipopt)=dosB(:,1+ipopt+ipopt:ne);
else
    dosB(:,1:ne+ipopt)=dosB(:,1-ipopt:ne);
end

% Convolution with a gaussian (smooths minor changes, highlights relevants)

a = 0;
b = 0.15;
b2 = b^2;
x = -5:0.1:5;
f = exp(-((x-a).^2)/(2*b2)) / (sqrt(2*pi*b2));
f = f/trapz(f);

dos_dif = dosB - dosA;

if nq==1
  dos_dif = conv(f,dos_dif,'same');
else
  dos_dif = conv2(f,f,dos_dif,'same');
end

% --------- Plot of the density ---------

if nq==1            % if Gamma point only
  plot(e,dos_dif)
  xlabel('Energy')
  ylabel('DOS')
  grid on

else
  s = surf(qpath',e',dos_dif');    
  s.EdgeColor = 'none';
  set(gca,'FontSize',12);
  ylabel('Energy  (eV)')
  zlabel('LDOS')
  zlimit = max(max(max(dos_dif)),abs(min(min(dos_dif))))*1.01;
  axis([0,qpath(end),emin,emax,-zlimit,zlimit])
  caxis([-zlimit,zlimit])

  view(2);
  colormap jet;
  colorbar

  % Plot line edges  

  lineEnd = find(iline(1:nq-1)<iline(2:nq));      % Line edges  
  hold on
  clr = [0.85,0.325,0.098];            % burgundy: [0.635,0.078,0.184];
  for jl = 1:numel(lineEnd)            %   orange: [0.850,0.325,0.098];
      il = lineEnd(jl);
      plot3(qpath(il)*[1,1],[emin,emax],[zlimit,zlimit],'Color',clr,'LineStyle',':')
  end
  plot3(qpath(end)*[1,1],[emin,emax],[zlimit,zlimit],'Color',clr,'LineStyle',':') % last sym point

  plot3([qpath(1),qpath(end)],0*[1,1],[zlimit,zlimit],'Color','k') % Fermi A
  plot3([qpath(1),qpath(end)],(efermiB-efermiA)*[1,1],[zlimit,zlimit],'Color','r') % Fermi B (defective)

  vtick = [qpath((lineEnd(1:1:numel(lineEnd)))),qpath(end)];
  xticks(vtick);
  listoflabels = [ label(lineEnd(1:1:numel(lineEnd))); label(nq) ];
  xticklabels( listoflabels ); 
  hold off
end


% -------------------------------------------------------------------------
