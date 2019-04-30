% plotmesh.m 
% Plot a surface of qDOS writen by unfold in a 2D mesh of
% q points. Reads a colection of .pathN files written by 
% unfold. The UnfoldedBandLines block can be generated with 
% the program meshgen.
% SGM & JMS, April 2019

clear all
ry = 13.60580;  % one ry in eV
wmin = 1.e-4;

% Change this:
  npaths = 22;                  % number of files to read
  wdir = '../divacancy/mesh02/';       % path to files
  syslabel = 'Divacancy-11X11';        % Siesta label
  flabel = '.refoldedBands.path';      % unfolded or refolded

for ipath=1:npaths              
  spath = num2str(ipath);         % turns number into string
  fname = [wdir,syslabel,flabel,spath];

  file = fopen(fname);
  dat = fscanf(file,'%d %d %f %f %f',5);
  nq = dat(1);       ne = dat(2);
  emin = dat(3);     emax = dat(4);   
  Fermi = dat(5);
  
  if ipath==1
    q = zeros(3,nq,npaths);
    dos = zeros(nq,npaths,ne);  % label = string(zeros(nq,1,npaths));
  end

  for iq = 1:nq
    q(:,iq,ipath) = fscanf(file,'%f %f %f',3);
    fscanf(file,'%i',1);        %   iline(iq,ipath) = 
    fscanf(file,'%10c',1);      %   label(iq) = 
    dos(iq,ipath,:) = fscanf(file,'%f',ne);
  end
  fclose(file);
end
  
emin = emin - Fermi;
emax = emax - Fermi;
de = (emax-emin)/(ne-1);
emesh = emin:de:emax;
  
qx(:) = q(1,1,:);
qy(:) = q(2,:,1);
[qx,qy,emesh] = ndgrid(qx,qy,emesh);

dos = permute(dos,[2,1,3]);

dos = dos+1e-30;
ie = 2:ne-1;
dos2(:,:, 1,2) = dos(:,:,1);         % unfold splits dos in two ie 
dos2(:,:,ne,1) = dos(:,:,ne);
dos2(:,:,ie,1) = dos(:,:,ie).*dos(:,:,ie-1)./(dos(:,:,ie-1)+dos(:,:,ie+1));
dos2(:,:,ie,2) = dos(:,:,ie).*dos(:,:,ie+1)./(dos(:,:,ie-1)+dos(:,:,ie+1));
ie = 1:ne-1;
w(:,:,ie) = (dos2(:,:,ie,2)+dos2(:,:,ie+1,1));
e(:,:,ie) = (emesh(:,:,ie).*dos2(:,:,ie,2)+...
             emesh(:,:,ie+1).*dos2(:,:,ie+1,1))./w(:,:,ie);
ne = ne-1;
e(:,:,ie+1)=e(:,:,ie);
w(:,:,ie+1)=w(:,:,ie);

% Select points with some significant weight
iqe = find(w(:)>wmin);
qx = qx(iqe);
qy = qy(iqe);
e = e(iqe);
w = w(iqe);
w = w/max(w);

[~,iw] = sort(w);
w = w(iw);
qx = qx(iw);
qy = qy(iw);
e = e(iw);
nqe = numel(e);

ptsize = 75;    % size of dots
scatter3(qx,qy,e,ptsize*w.^1.1,w,'filled');

colormap(flipud(bone));
xlabel('qy')
xlabel('qx')
zlabel('energy')
 
