% plotdos.m
% Plot qDOS written by unfold
% SGM & JMS, Oct.2018

clear
wdir = './path_to_file/';
syslabel = 'si_fcc';
fname = [wdir,syslabel,'.refoldedBands'];

% Read DOS
file = fopen(fname);
dat = fscanf(file,'%d %d %f %f',4);
nq = dat(1);
ne = dat(2);
emin = dat(3);
emax = dat(4);
q = zeros(3,nq);
dos = zeros(nq,ne);
ST = zeros(nq,1);
  label = string(ST);
for iq = 1:nq
    q(:,iq) = fscanf(file,'%f %f %f',3);
    iline(iq) = fscanf(file,'%i',1);
    label(iq) = fscanf(file,'%10c',1);
    dos(iq,:) = fscanf(file,'%f',ne);
end
fclose(file);

% Plot qDOS
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

if nq==1
    plot(e,dos)              % normal
%    plot(e,cumtrapz(dos))     % integrated
    xlabel('Energy (eV)')
    ylabel('DOS')
    set(gca,'FontSize',12.5)
    grid on

else
    s = surf(qpath',e',dos');                      % normal (refolding)
%    s = surf(qpath',e',cumtrapz(dos,2)');          % integrated
%    s = surf(qpath',e',max(log10(dos),-2)');       % log scale (unfolding)
    s.EdgeColor = 'none';
    
    xlabel('qpath')
    ylabel('Energy (eV)')
    zlabel('DOS')
    set(gca,'FontSize',11.5)
    
    zlimit = max(max(dos))*1.01;        % 'auto' saturates at max dos
    axis([0,qpath(end),emin,emax,0,zlimit])        % normal
    caxis([0,zlimit]); caxis('auto'); 
%    axis([0,qpath(end),emin,emax,0,26])           % integrated
%    axis([0,qpath(end),emin,emax,-2,0.25])        % log scale
    caxis('auto');
end

% Set viewpoint and colormap
view(2);
 colormap(flipud(bone))    % black on white
% colormap(bone);          % white on black
% colormap default
colorbar

% Plot line edges  
if nq ~= 1
    lineEnd = find(iline(1:nq-1)<iline(2:nq));
    hold on
    clr=[0.85,0.325,0.098];    % burgundy:[0.635,0.078,0.184];orange:[0.85,0.325,0.098]
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