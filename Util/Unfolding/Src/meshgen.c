// Generates the UnfoldedBandLines block for a squared 2D 
// q-point mesh. Reads the input from the 'meshgen.dat' 
// file, with format:
//    ne   emin   emax
//    nx   x0   xend
//    ny   y0   yend
// and writes the formated output 'meshgen.fdf'.
// After running unfold, plot the output files with 
// 'plotmesh.m'.
//
// April 2019 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.14159265

int main() {

FILE *dat;
dat = fopen("meshgen.dat","r");

FILE *fdf;
fdf = fopen("meshgen.fdf","w");

int ne, nx, ny;
double emin, emax;
double x0, xend, y0, yend;

int i;
double xi, dx;

fscanf(dat, "%d %lf %lf", &ne, &emin, &emax);
fscanf(dat, "%d %lf %lf", &nx, &x0, &xend);
fscanf(dat, "%d %lf %lf", &ny, &y0, &yend);

xi = x0;
dx = (xend - x0)/(nx - 1.0);
ny = ny - 1;

fprintf(fdf, "%%block UnfoldedBandLines\n");
fprintf(fdf, "   %d\t%lf\t%lf\t eV\n", ne, emin, emax);

i = 0;
for ( i=0; i<nx; i++ ) {
  xi = xi + dx; 
  fprintf(fdf, " 1\t%lf\t%lf\t0.0\n", xi, y0);
  fprintf(fdf, " %d\t%lf\t%lf\t0.0\n", ny, xi, yend);
}
fprintf(fdf, "%%endblock UnfoldedBandLines\n");
printf("'meshgen.fdf' written\n");

fclose(dat);
fclose(fdf);

return 0;
}

