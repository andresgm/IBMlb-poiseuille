/*
// Librerias estandar
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fluid.h"
#include "fronteras.h"

// Weight for the equilibrium distribution according with Hecht
double w[19] = {(2./36.),(2./36.),(2./36.),(2./36.),(2./36.),(2./36.),
      (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
      (1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
      (12./36.)};

double e_x[19] = {1.,-1.,0.,0.,0.,0.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,0.,0.,0.,0.,0.};
double e_y[19] = {0.,0.,1.,-1.,0.,0.,1.,-1.,0.,0.,1.,-1.,0.,0.,1.,1.,-1.,-1.,0.};
double e_z[19] = {0.,0.,0.,0.,1.,-1.,0.,0.,1.,-1.,0.,0.,1.,-1.,1.,-1.,1.,-1.,0.};
int dfInv[19] = {1,0,3,2,5,4,11,10,13,12,7,6,9,8,17,16,15,14,18};

// different types of cells
const int FLUIDO  = 3, TOP = 1, BOTTOM  = 2, NOSLIP = 0;
// velocity for the moving wall
const double V = 0.00, W = 0.00, cs=1.0/sqrt(3.0);

// Declare some constants and globals ...
const double   omega = 2.0;     // viscosity of the fluid, 0..2
int current = 0, other = 1; // which grid is the current and the old one?


 * 1. Construye las estructuras de información macro
 * 2. Construye las estructuras de informacion micro
 * 3. Inicializa el valor de las distribuciones en la malla
 * @param double x
 * @param double y
 * @param double z
 * @return void

void fluid::inicializar(int x, int y, int z,double r)
{
	    X = x;
		Y = y;
		Z = z;

		// Construye las estructuras de información macroscopicas codeblocks, ubuntu
		vel = new double***[x];
		fuerza = new double***[x];
		rho = new double**[x];
		flags = new double**[x];
		for(int i=0;i<x;i++)
		{
			vel[i] = new double**[y];
			fuerza[i] = new double**[y];
			rho[i] = new double*[y];
			flags[i] = new double*[y];
			for(int j = 0; j<y;j++)
			{
				vel[i][j] = new double*[z];
				fuerza[i][j] = new double*[z];
				rho[i][j] = new double[z];
				flags[i][j] = new double[z];
				for(int k=0; k<z ; k++)
				{
					vel[i][j][k] = new double[3];
					fuerza[i][j][k] = new double[3];
					rho[i][j][k] = 0;
					flags[i][j][k]=0;
				}
			}
		}

		// Construye las estructuras de información microscopicas
		cells = new double****[2];
		for(int s=0;s<2;s++)
		{
			cells[s]=new double***[x];
			for(int i=0;i<x;i++)
			{
				cells[s][i]=new double**[y];
				for(int j = 0; j<y;j++)
				{
					cells[s][i][j]=new double*[z];
					for(int k=0; k<z ; k++)
					{
						cells[s][i][j][k]=new double[19];
						for(int l=0;l<19;l++)
						{
							cells[s][i][j][k][l] = w[l];
						}
					}
				}
			}
		}

		for (int i=0;i<X;i++)
				for (int j=0;j<Y;j++)
					for (int k=0;k<Z;k++) {
						for (int l=0;l<19;l++)  {
							cells[0][i][j][k][l] = cells[1][i][j][k][l] = w[l];
						}
						if(k==0){flags[i][j][k]==BOTTOM;}
						else if(k==Z-1){flags[i][j][k]==TOP;}
						else{flags[i][j][k]==FLUIDO;}
					}
}

void fluid::stream()
{
	for (int i=0;i<X;i++)
				for (int j=0;j<Y;j++)
					for (int k=1;k<Z-1;k++)
						for (int l=0;l<19;l++) {
							int inv = dfInv[l];
							int a = i + e_x[inv];
						    int b = j + e_y[inv];
						    int c = k + e_z[inv];

						        // Periodico en x
						        if(a<0){a=X-1;}
						        if(a>(X-1)){a=0;}

						        // Periodico en y
						        if(b<0){b=Y-1;}
						        if(b>(Y-1)){b=0;}
						        if(flags[a][b][c]==TOP){
								// Bounce - back
								cells[current][i][j][k][l] = cells[other][i][j][k][inv];}
							else{
								// Streaming - normal
								cells[current][i][j][k][l] = cells[other][a][b][c][l];}
						}//Stream

			for (int i=0;i<X;i++)
				for (int j=0;j<Y;j++){
				velNodoInferior(cells[current][i][j][0],cells[other][i][j][0], U, V, W);
				velNodoSuperior(cells[current][i][j][Z-1],cells[other][i][j][Z-1], U, V, W);
				}
}

void fluid::collide()
{
	// collision step
			for (int i=0;i<X;i++)
				for (int j=0;j<Y;j++)
					for (int k=0;k<Z;k++) {

						// standard collide without turbulence model
						double rho = 0.0, u_x=0.0, u_y=0.0, u_z=0.0;
						// normal fluid cell
						for (int l=0;l<19;l++) {
							const double fi = cells[current][i][j][k][l];
							rho += fi;
							u_x += e_x[l]*fi;
							u_y += e_y[l]*fi;
							u_z += e_z[l]*fi;
						}
						if(k==0){u_x=-U;}
						if(k==Z-1){u_x=U;}
						for (int l=0;l<19;l++) {
							const double tmp = (e_x[l]*u_x + e_y[l]*u_y + e_z[l]*u_z);
							double feq = w[l] * ( rho -
								(3.0/2.0 * (u_x*u_x + u_y*u_y + u_z*u_z)) +
								(3.0 *     tmp) +
								(9.0/2.0 * tmp*tmp ) );
							cells[current][i][j][k][l] =
								(1.0-omega) * cells[current][i][j][k][l] +
								omega * feq;
						}
					} // ijk
			// We're done for one time step, switch the grid...
			other = current;
			current = (current+1)%2;
}

void fluid::calcularMacro()
{
	for(int i = 0 ;i<X;i++)
			for(int j = 0 ;j<Y;j++)
				for(int k = 0 ;k<Z;k++){
					double rho=0.0;
					double u_x=0.0;
					double u_y=0.0;
					double u_z=0.0;
					for(int l = 0 ;l<19;l++){
						const double fi = cells[current][i][j][k][l];
						rho+= fi;
						u_x+=fi*e_x[l];
						u_y+=fi*e_y[l];
						u_z+=fi*e_z[l];
					}
					u_x=u_x/rho;
					u_y=u_y/rho;
					u_z=u_z/rho;
					vel[i][j][k][0] = (u_x)/rho;
					vel[i][j][k][1] = (u_y)/rho;
					vel[i][j][k][2] = (u_z)/rho;
					fuerza[i][j][k][0]=0.0;
					fuerza[i][j][k][1]=0.0;
					fuerza[i][j][k][2]=0.0;
				}
}

// Save fluid in structured grid format .vts
int fluid::guardar(int s)
{
	FILE *archivo;El manejador de archivo
	char ruta[80];
	char numero[10];
	strcpy(ruta, "temp/fluido-");
	sprintf(numero,"%d",s);
	strcat(ruta,numero);
	strcat(ruta,".vtk");
        archivo=fopen(ruta, "w");
        if(archivo==NULL){Si no lo logramos abrir, salimos
		printf("No se puede guardar archivo");
		return 1;}
	else{
	// Escribir datos al archivo
	// 1. Escribir cabecera.
	fprintf(archivo, "# vtk DataFile Version 3.0\n");
	fprintf(archivo, "vtk output\n");
	fprintf(archivo, "ASCII\n");
	fprintf(archivo, "DATASET STRUCTURED_GRID\n");
	fprintf(archivo, "DIMENSIONS %i %i %i\n", X, Y, Z);
	fprintf(archivo, "POINTS %i double\n",X*Y*Z);
	// 2. Escribir coordenadas de cada punto
	for(int i = 0 ;i<X;i++)
		for(int j = 0 ;j<Y;j++)
			for(int k = 0 ;k<Z;k++){
				fprintf(archivo,"%d %d %d\n", i,j,k);
			}
	// 3. Escribir datos sobre puntos

	// Densidad
	fprintf(archivo, "POINT_DATA %d\n", X*Y*Z);
	fprintf(archivo, "SCALARS Densidad double\n");
	fprintf(archivo, "LOOKUP_TABLE default\n");
	for(int i = 0 ;i<X;i++)
		for(int j = 0 ;j<Y;j++)
			for(int k = 0 ;k<Z;k++){
				if(flags[i][j][k] == FLUIDO)
				{
					fprintf(archivo, "%f \n", darDensidad(i,j,k));
				}
				else{
					fprintf(archivo, "%f \n", 0.0);
				}
			}

	// Velocidad
	fprintf(archivo, "\nVECTORS Velocidad double\n");
	for(int i = 0 ;i<X;i++)
				for(int j = 0 ;j<Y;j++)
					for(int k = 0 ;k<Z;k++){
						fprintf(archivo, "%f %f %f\n", darVelocidad(i,j,k,0),darVelocidad(i,j,k,1),darVelocidad(i,j,k,2));
					}
	// Escribir vectores fuerza en el algoritmo
	fprintf(archivo, "VECTORS fuerza double\n");
		for(int i = 0 ;i<X;i++)
			for(int j = 0 ;j<Y;j++)
				for(int k = 0 ;k<Z;k++)
				{
					fprintf(archivo, "%f %f %f\n", fuerza[i][j][k][0], fuerza[i][j][k][0], fuerza[i][j][k][0]);
				}
        fclose(archivo);Cerramos el archivo
        return 0;
	}
}


double fluid::darVelocidad(int x, int y, int z, int f)
{
	return vel[x][y][z][f];
}


double fluid::darDensidad(int x, int y, int z)
{
	double rho = 0.0;
	for(int l = 0; l < 19 ; l++)
	{
		rho += cells[current][x][y][z][l];
	}
	return rho;
}

void fluid::setFuerza(int x, int y, int z, double f[3])
{
	fuerza[x][y][z][0] = f[0];
	fuerza[x][y][z][1] = f[1];
	fuerza[x][y][z][2] = f[2];
}

void fluid::addFuerza(int x, int y, int z, double f[3])
{
	fuerza[x][y][z][0] += f[0];
	fuerza[x][y][z][1] += f[1];
	fuerza[x][y][z][2] += f[2];
}

void fluid::setVelocidad(double u)
{
	U = u;
}
*/
