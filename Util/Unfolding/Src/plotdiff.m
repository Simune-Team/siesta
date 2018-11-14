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
wdir = './path_to_files/';
syslabelA = 'si_supercell';           % System A : non-defective
syslabelB = 'si_vacancy';           % System B : defective
fnameA = [wdir,syslabelA,'.refoldedBands'];          % '.refoldedBands' or
fnameB = [wdir,syslabelB,'.refoldedBands'];          %   '.unfoldedBands'.

% Read DOS
fileA = fopen(fnameA);
datA = fscanf(fileA,'%d %d %f %f',4);
nq = datA(1);
ne = datA(2);
emin = datA(3);
emax = datA(4);
fileB = fopen(fnameB);
datB = fscanf(fileB,'%d %d %f %f',4);
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

% Convolution with a gaussian (smooths minor changes, highlights relevants)
a = 0;
b = 0.13;
b2 = b^2;
x = -5:0.1:5;
f = exp(-((x-a).^2)/(2*b2)) / (sqrt(2*pi*b2));
f = f/trapz(f);

dos_dif = dosB - dosA;
dos_dif = conv2(f,f,dos_dif,'same');

% -------------- Plot of the density --------------

if nq==1
    plot(e,dosA)
    xlabel('energy')
    ylabel('dos')
    grid on
else
    s = surf(qpath',e',dos_dif');    
    s.EdgeColor = 'none';
    xlabel('qpath')
    ylabel('energy')
    zlabel('dos')
    zlimit = max(max(max(dos_dif)),abs(min(min(dos_dif))))*1.01;
    axis([0,qpath(end),emin,emax,-zlimit,zlimit])
    caxis([-zlimit,zlimit])
end

% Set viewpoint and colormap
view(2);
colormap jet(21);
colorbar

% Plot line edges  
if nq ~= 1
    lineEnd = find(iline(1:nq-1)<iline(2:nq));
    hold on
    clr=[0.85,0.325,0.098]; % burgundy:[0.635,0.078,0.184];orange:[0.85,0.325,0.098]
    for jl = 1:numel(lineEnd)
        il = lineEnd(jl);
        plot3(qpath(il)*[1,1],[emin,emax],[zlimit,zlimit],'Color',clr,'LineStyle',':')
    end
    plot3(qpath(end)*[1,1],[emin,emax],[zlimit,zlimit],'Color',clr,'LineStyle',':') % last sym point

    ntick = numel(lineEnd)+1;
    vtick = [qpath((lineEnd(1:1:ntick-1))),qpath(end)];
    xticks(vtick);
    listoflabels = [ label(lineEnd(1:1:ntick-1)); label(nq) ];
    xticklabels( listoflabels );
    
   % xticklabels({'\Gamma\color{gray}(000)','X','L'});    % to specify labels
    
    hold off
end


% -------------------------------------------------------------------------
