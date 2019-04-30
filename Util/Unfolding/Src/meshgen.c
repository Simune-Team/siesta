// Generates the UnfoldedBandLines block for a 2D qpoint mesh
// Reads the input from a 'meshgen.dat' file with the format:
//    ne   emin   emax
//    nx   ny
//    x0   xend
//    y0   yend
// and writes the formated output 'meshgen.fdf'
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
fscanf(dat, "%d %d", &nx, &ny);
fscanf(dat, "%lf %lf", &x0, &xend);
fscanf(dat, "%lf %lf", &y0, &yend);

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

