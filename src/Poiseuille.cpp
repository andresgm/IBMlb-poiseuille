#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ibm.h"
#include "fluid.h"
#include "mesh.h"

using namespace std;

int main(int argc, char *argv[])
{
	int N=10;
	mesh membrana[N];
	mesh referencia[N];
	fluid fluido;
	double dt = 1.0;
	double dx = 1.0;
	int X = 40;
	int Y = 100;
	int Z = 40;
	int VTK = 10;

	// Parametros adimensionales
	double rho = 1.0;
	double nu = 1./6.;
	double Re = 50.0;
	double G = 0.1;
	double R = 10.0;
	double gamma_dot = (Re*nu)/(rho*pow(R,2));
	double ks = (gamma_dot*nu*R)/(G);
	double kb = ks*1.0e-6;
	double kp = (gamma_dot)/(G);
	double STEPS = 100000/kp;
	printf("A completar %f iteraciones\n", STEPS);

	// Membrana
	for(int i = 0; i < N; i++)
	{
		membrana[i].setID(i);
		membrana[i].mesh_refine_tri4();
		membrana[i].mesh_refine_tri4();
		membrana[i].mesh_refine_tri4();
		membrana[i].proyectarEsfera(R);
		membrana[i].proyectarRBC(R);
		membrana[i].rotarEstructura(3.1416/2.0,0.0,0.0);
		membrana[i].moverCentro((X-1)/2.0,10.0*i, (Z-1)/2.0);
		membrana[i].iniciarGeometria();

		referencia[i].setID(i);
		referencia[i].mesh_refine_tri4();
		referencia[i].mesh_refine_tri4();
		referencia[i].mesh_refine_tri4();
		referencia[i].proyectarEsfera(R);
		referencia[i].proyectarRBC(R);
		referencia[i].rotarEstructura(3.1416/2.0,0.0,0.0);
		referencia[i].moverCentro((X-1)/2.0,10.0*i, (Z-1)/2.0);
		referencia[i].iniciarGeometria();
		referencia[i].actualizarGeometria();
	}

	// Fluido
	fluido.inicializar(X,Y,Z);
	fluido.setVelocidad(gamma_dot);

	for(int ts = 0 ; ts < STEPS ; ts++)
		{
		// -----------------------------------------------------------------------//
		// 1. Interpolation
		// -----------------------------------------------------------------------//
		for(int i = 0; i<N; i++)
		{
			interpolation(fluido, membrana[i], X, Y, Z);
		}

		// -----------------------------------------------------------------------//
		// 2. Encontrar nuevas posiciones de la membrana
		// -----------------------------------------------------------------------//
		for(int i = 0; i < N; i++)
		{
			membrana[i].moverNodos(dt, dx, Y);
		}

		// -----------------------------------------------------------------------//
		// 3. Calcular fuerzas en los nodos de la membrana
		// -----------------------------------------------------------------------//
		for(int i = 0; i < N; i++)
		{
			membrana[i].calcularFuerzasHelfrich(kb);
			membrana[i].calcularFuerzasFEM(referencia[i], ks);
		}

		// -----------------------------------------------------------------------//
		// 4. Propagar la densidad de fuerza hacia el fluido
		// -----------------------------------------------------------------------//
		for(int i = 0; i < N; i++)
		{
			spread(fluido, membrana[i], X, Y, Z);
		}

		// -----------------------------------------------------------------------//
		// 5. Solucionar la dinámica del fluido
		// -----------------------------------------------------------------------//
		fluido.collide();
		fluido.stream();

		// -----------------------------------------------------------------------//
		// 6. Calcular propiedades macro del fluido
		// -----------------------------------------------------------------------//
		fluido.calcularMacro();

		// -----------------------------------------------------------------------//
		// 7. Calcular propiedades macro de la membrana
		// -----------------------------------------------------------------------//
		for(int i = 0; i<N;i++)
		{
			membrana[i].calcularCambioArea(referencia[i]);
			membrana[i].actualizarGeometria();
		}

		// -----------------------------------------------------------------------//
		// 9. Visualización
		// -----------------------------------------------------------------------//
		if(ts%VTK==0)
		{
			fluido.guardar(ts);
			for(int i=0;i<N;i++)
			{
				membrana[i].guardarVTU(ts);
			}
			printf("%d\n",ts);
		}
	}//Ciclo principal

	return 0;
}
