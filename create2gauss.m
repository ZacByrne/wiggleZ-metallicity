function [] = create2gauss(nocol)
%   function takes filename and then creates spectra file
%   
%   input; filename of the starlight output file
%   output; labels of data
%           x = lambda
%           y = matrix of [flux, mass, weighting]

%Read in Error bars
filestr = fileread(['/home/uqzbyrne/WiggleZ-Metallicity/Spectra/wig02bin' num2str(nocol) '.txt']);
c =textscan(filestr, '%f %f %f %f');

z = [c{1}, c{3}];

er4363 = z((find(z == 4300)):(find(z==4500)),:);

erhb = z((find(z == 4800)):(find(z==5100)),:);

er3728 = z(1:(find(z==4000)),:);


% Find Cont level
filename = ['/home/uqzbyrne/WiggleZ-Metallicity/wiggleZ-metallicity-master/metal/sl_fits/wig02hb' num2str(nocol) '.BN'];

header = findspectra(filename);

[labels, x, y] = readColData(filename, 4, header-1, 2);

average1 = (mean(y(find(x == 4800)):(find(x==4835))));
average2 = (mean(y(find(x == 4885)):(find(x==4925))));

contlevel = (average1+average2)/2;


linenum = 15;

filename = ['/home/uqzbyrne/WiggleZ-Metallicity/wiggleZ-metallicity-master/metal/sl_fits/wig02fit' num2str(nocol) '.BN'];
norma = textscan(filename,' %f %s',1, 'delimiter', '\n','headerlines', linenum-1);
normfull = norma(1);

filename = ['/home/uqzbyrne/WiggleZ-Metallicity/wiggleZ-metallicity-master/metal/sl_fits/' num2str(nocol) '4363.BN'];
norma = textscan(filename,' %f %s',1, 'delimiter', '\n','headerlines', linenum-1);
norm4363 = norma(1);

filename = ['/home/uqzbyrne/WiggleZ-Metallicity/wiggleZ-metallicity-master/metal/sl_fits/wig02hb' num2str(nocol) '.BN'];
norma = textscan(filename,' %f %s',1, 'delimiter', '\n','headerlines', linenum-1);
normhb = norma(1);

fclose('all')


%create outfile
filename = ['/home/uqzbyrne/WiggleZ-Metallicity/wiggleZ-metallicity-master/metal/sl_fits/wig02fit' num2str(nocol) '.BN'];

header = findspectra(filename);

[labels, x, y] = readColData(filename, 4, header-1, 2);

index = find(x == 4000);
spectrafile = fopen(['fitwig02' num2str(nocol) '.txt'], 'w');
%fprintf(spectrafile, ['Wavelength   flux  \n']);
fprintf(spectrafile, '%f %f %f \n', [x(1:index), (y(1:index,1)-y(1:index,2))*normfull+contlevel, er3728(:,2)].');
fclose(spectrafile);

filename = ['/home/uqzbyrne/WiggleZ-Metallicity/wiggleZ-metallicity-master/metal/sl_fits/' num2str(nocol) '4363.BN'];

header = findspectra(filename);

[labels, x, y] = readColData(filename, 4, header-1, 2);


spectrafile = fopen(['fitwig02' num2str(nocol) '.txt'], 'a+');
%fprintf(spectrafile, ['Wavelength   flux  \n']);
fprintf(spectrafile, '%f %f %f \n', [x, (y(:,1)-y(:,2))*norm4363+contlevel, er4363(:,2)].');
fclose(spectrafile);

filename = ['/home/uqzbyrne/WiggleZ-Metallicity/wiggleZ-metallicity-master/metal/sl_fits/wig02hb' num2str(nocol) '.BN'];

header = findspectra(filename);

[labels, x, y] = readColData(filename, 4, header-1, 2);


spectrafile = fopen(['fitwig02' num2str(nocol) '.txt'], 'a+');
%fprintf(spectrafile, ['Wavelength   flux  \n']);
fprintf(spectrafile, '%f %f %f \n', [x, (y(:,1)-y(:,2))*normhb + contlevel, erhb(:,2)].');
fclose(spectrafile);
