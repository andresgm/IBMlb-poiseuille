#ifndef ELEMENTO_H
	#define ELEMENTO_H
	using namespace std;

	/*
	 * elemento.cpp
	 *
	 *  Created on: Mar 23, 2011
	 *      Author: oscar
	 */

	#include "elemento.h"
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <string.h>

	using namespace std;

	/**
	*	Función para calcular las fuerzas en un elemento triangular
	*   @param double referencia[3][3], Vertices iniciales del elemento
	*	@param double deformado[3][3], Vertices del elemento deformado en el plano
	*	@return double** fuerzas, Apuntador hacia las fuerzas calculadas en cada nodo
	*/

	void fuerzas(double referencia[3][3], double deformado[3][3], double fuerzas[3][3], double ks)
	{
			// Coordenadas de los tres vertices iniciales
			double xi=referencia[0][0];
			double yi=referencia[0][1];
			//double zi=referencia[0][2];

			double xj=referencia[1][0];
			double yj=referencia[1][1];
			//double zj=referencia[1][2];

			double xk=referencia[2][0];
			double yk=referencia[2][1];
			//double zk=referencia[2][2];

			// Coordenadas de los tres vertices del elemento deformado
			double Xi=deformado[0][0];
			double Yi=deformado[0][1];
			//double Zi=deformado[0][2];

			double Xj=deformado[1][0];
			double Yj=deformado[1][1];
			//double Zj=deformado[1][2];

			double Xk=deformado[2][0];
			double Yk=deformado[2][1];
			//double Zk=deformado[2][2];

			// Desplazamientos de cada nodo
			double ui=Xi-xi;
			double vi=Yi-yi;

			double uj=Xj-xj;
			double vj=Yj-yj;

			double uk=Xk-xk;
			double vk=Yk-yk;

			// Coeficientes de funciones de forma
			double ai = yj-yk;
			double bi = xk-xj;
			double ci = xj*yk - xk*yj;
			double Li = ai*xi + bi*yi + ci;

			double aj = yk - yi;
			double bj = xi - xk;
			double cj = xk*yi - xi*yk;
			double Lj = aj*xj + bj*yj + cj;

			double ak = yi - yj;
			double bk = xj - xi;
			double ck = xi*yj - xj*yi;
			double Lk = ak*xk + bk*yk + ck;

			// Derivadas parciales para calcular el vector [G] revision Sep/14/2010
			double dudx = ui*ai/Li + uj*aj/Lj + uk*ak/Lk;
			double dudy = ui*bi/Li + uj*bj/Lj + uk*bk/Lk;
			double dvdx = vi*ai/Li + vj*aj/Lj + vk*ak/Lk;
			double dvdy = vi*bi/Li + vj*bj/Lj + vk*bk/Lk;

			// Componentes del vector [G] revision Sep/14/2010
			double g11 = (1.0+dudx)*(1.0+dudx) + (dvdx)*(dvdx);
			double g12 = (1.0+dudx)*(dudy) + (1.0+dvdy)*(dvdx);
			//double g21 = g12;
			double g22 = (1.0+dvdy)*(1.0+dvdy) + (dudy)*(dudy);
			//double G[2][2] = {{g11,g12},{g21,g22}};

			// Calculo de lambda1 y lambda2 revision Sep/14/2010
			double l1 = sqrt((g11 + g22 + sqrt((g11-g22)*(g11-g22) + 4.*g12*g12))/2.);
			double l2 = sqrt((g11 + g22 - sqrt((g11-g22)*(g11-g22) + 4.*g12*g12))/2.);

			// Derivadas de la funcion Strain Energy respecto a lambda 1 y lambda 2
			// Modelo de energia Skalak 1973 revision Sep/14/2010
			double I1 = (l1*l1) + (l2*l2) - 2.0;
			double I2 = (l1*l1)*(l2*l2) - 1.0;
			double dI1dl1 = 2.0*l1;
			double dI1dl2 = 2.0*l2;
			double dI2dl1 = 2.0*l1*(l2*l2);
			double dI2dl2 = 2.0*l2*(l1*l1);

			//    dwdl1 = (B/4)*(I1*dI1dl1 + dI1dl1 - dI2dl1) + (C/4)*(I2)*(dI2dl1)
			//    dwdl2 = (B/4)*(I1*dI1dl2 + dI1dl2 - dI2dl2) + (C/4)*(I2)*(dI2dl2)
			double dwdl1 = (ks/12.)*(2.*I1*dI1dl1 + 2.0*dI1dl1 - 2.0*dI2dl1) + (ks/6.0)*I2*dI2dl1;
			double dwdl2 = (ks/12.)*(2.*I1*dI1dl2 + 2.0*dI1dl2 - 2.0*dI2dl2) + (ks/6.0)*I2*dI2dl2;

			// Calculo de diferenciales sobre l1 y l2 respecto a desplazamientos de nodos
			// 1. Derivadas de [G] respecto a desplazamiento de nodos revision Sep/14/2010
			double dg11dui = 2.0*(1.0+ dudx)*(ai/Li);
			double dg11duj = 2.0*(1.0+ dudx)*(aj/Lj);
			double dg11duk = 2.0*(1.0+ dudx)*(ak/Lk);
			double dg11dvi = 2.0*dvdx*(ai/Li);
			double dg11dvj = 2.0*dvdx*(aj/Lj);
			double dg11dvk = 2.0*dvdx*(ak/Lk);

			double dg12dui = (1.0+dudx)*(bi/Li) + (ai/Li)*(dudy);
			double dg12duj = (1.0+dudx)*(bj/Lj) + (aj/Lj)*(dudy);
			double dg12duk = (1.0+dudx)*(bk/Lk) + (ak/Lk)*(dudy);
			double dg12dvi = (1.0+dvdy)*(ai/Li) + (bi/Li)*(dvdx);
			double dg12dvj = (1.0+dvdy)*(aj/Lj) + (bj/Lj)*(dvdx);
			double dg12dvk = (1.0+dvdy)*(ak/Lk) + (bk/Lk)*(dvdx);

			double dg22dui = 2.0*dudy*(bi/Li);
			double dg22duj = 2.0*dudy*(bj/Lj);
			double dg22duk = 2.0*dudy*(bk/Lk);

			double dg22dvi = 2.0*(1.0+dvdy)*(bi/Li);
			double dg22dvj = 2.0*(1.0+dvdy)*(bj/Lj);
			double dg22dvk = 2.0*(1.0+dvdy)*(bk/Lk);

			// 2. Calculo de las derivadas de lambda 1 y lambda 2 respecto desplazamientos
			// nodales revision Sep/14/2010
			double t0 = sqrt((g11-g22)*(g11-g22) + 4.*g12*g12);

			// Derivadas de lambda1  y lambda 2

			double dt0dui;
			if(t0>10e-3){
				dt0dui = ((((g11-g22)*(dg11dui-dg22dui))+(4.*g12*dg12dui))/(t0));}
			else{
				dt0dui = 0.0;}
			double dl1dui = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dui + dg22dui + dt0dui);
			double dl2dui = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dui + dg22dui - dt0dui);

			double dt0duj;
			if(t0>10e-3){
				dt0duj = ((((g11-g22)*(dg11duj-dg22duj))+(4.*g12*dg12duj))/(t0));}
			else{
				dt0duj = 0.0;}
			double dl1duj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11duj + dg22duj + dt0duj);
			double dl2duj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11duj + dg22duj - dt0duj);

			double dt0duk;
			if(t0>10e-3){
				dt0duk = ((((g11-g22)*(dg11duk-dg22duk))+(4.*g12*dg12duk))/(t0));}
			else{
				dt0duk = 0.0;}
			double dl1duk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11duk + dg22duk + dt0duk);
			double dl2duk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11duk + dg22duk - dt0duk);

			double dt0dvi;
			if(t0>10e-3){
				dt0dvi = ((((g11-g22)*(dg11dvi-dg22dvi))+(4.*g12*dg12dvi))/(t0));}
			else{
				dt0dvi = 0.0;}
			double dl1dvi = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dvi + dg22dvi + dt0dvi);
			double dl2dvi = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dvi + dg22dvi - dt0dvi);

			double dt0dvj;
			if(t0>10e-3){
				dt0dvj = ((((g11-g22)*(dg11dvj-dg22dvj))+(4.*g12*dg12dvj))/(t0));}
			else{
				dt0dvj = 0.0;}
			double dl1dvj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dvj + dg22dvj + dt0dvj);
			double dl2dvj = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dvj + dg22dvj - dt0dvj);

			double dt0dvk;
			if(t0>10e-3){
				dt0dvk = ((((g11-g22)*(dg11dvk-dg22dvk))+(4.*g12*dg12dvk))/(t0));}
			else{
				dt0dvk = 0.0;}
			double dl1dvk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22+t0)))*(dg11dvk + dg22dvk + dt0dvk);
			double dl2dvk = ((sqrt(0.5)*0.5)/(sqrt(g11+g22-t0)))*(dg11dvk + dg22dvk - dt0dvk);

			// 3. Calculo de las derivadas de w respecto a los desplazamientos nodales
			// revision Sep/14/2010
			double dwdui = dwdl1*dl1dui + dwdl2*dl2dui;
			double dwdvi = dwdl1*dl1dvi + dwdl2*dl2dvi;
			double dwduj = dwdl1*dl1duj + dwdl2*dl2duj;
			double dwdvj = dwdl1*dl1dvj + dwdl2*dl2dvj;
			double dwduk = dwdl1*dl1duk + dwdl2*dl2duk;
			double dwdvk = dwdl1*dl1dvk + dwdl2*dl2dvk;

			// 4. Volumen del elemento revision Sep/14/2010
			double a0 = ((xj-xi)*(yk-yi) + (xk-xi)*(yj-yi))/2.0;
			double espesor = 0.5;//0.0025;

			// 5. Calculo de las componentes de fuerza revision Sep/14/2010
			double fxi = dwdui*a0*espesor;
			double fyi = dwdvi*a0*espesor;
			double fzi = 0.0;
			double fxj = dwduj*a0*espesor;
			double fyj = dwdvj*a0*espesor;
			double fzj = 0.0;
			double fxk = dwduk*a0*espesor;
			double fyk = dwdvk*a0*espesor;
			double fzk = 0.0;

			fuerzas[0][0] = fxi;
			fuerzas[0][1] = fyi;
			fuerzas[0][2] = fzi;
			fuerzas[1][0] = fxj;
			fuerzas[1][1] = fyj;
			fuerzas[1][2] = fzj;
			fuerzas[2][0] = fxk;
			fuerzas[2][1] = fyk;
			fuerzas[2][2] = fzk;
	}


	/**
	*	Función para realizar cross product entre dos vectores
	*   @param double a[3], Vector a
	*	@param double b[3], Vector b
	*	@param double &resultado[3], Paso por parametros
	*/
	void cross(double a[3], double b[3], double &x, double &y, double &z)
	{
		x = a[1]*b[2] - a[2]*b[1];
		y = a[2]*b[0] - a[0]*b[2];
		z = a[0]*b[1] - a[1]*b[0];
	}

	/**
	*	Función para calcular la norma de un vector
	*   @param double a[3], Vector a
	*	@return double norma, Norma del vector
	*/
	double norm(double a[3])
	{
		return sqrt( (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]));
	}

	/**
	*	Función para calcular el producto Matriz dot Vector
	*   @param M[3][3], Matriz de tamaño 3x3
	*   @param V[3], Vector de tamaño 3
	*   @param V[3], Vector de tamaño 3
	*/
	void MdotV(double M[3][3], double V[3], double &x, double &y, double &z)
	{
		x = M[0][0]*V[0] + M[0][1]*V[1] + M[0][2]*V[2];
		y = M[1][0]*V[0] + M[1][1]*V[1] + M[1][2]*V[2];
		z = M[2][0]*V[0] + M[2][1]*V[1] + M[2][2]*V[2];
	}

	/**
	*	Función para llevar los elementos de referencia y deformado a una misma base
	*   @param double referencia[3][3], Vertices iniciales del elemento
	*	@param double deformado[3][3], Vertices del elemento deformado en 3D
	*	@return double** vertices, Apuntador hacia los vertices en un plano común
	*/

	void rotacion(double referencia[3][3], double deformado[3][3], double nfuerzas[3][3], double ks)
	{
		// Coordenadas iniciales de los tres nodos respecto al sistema
	    // global de coordenadas i1, i2, i3
	    double xi = referencia[0][0];
	    double yi = referencia[0][1];
	    double zi = referencia[0][2];

	    double xj = referencia[1][0];
	    double yj = referencia[1][1];
	    double zj = referencia[1][2];

	    double xk = referencia[2][0];
	    double yk = referencia[2][1];
	    double zk = referencia[2][2];

	    // Coordenadas del elemento deformado de los tres nodos
	    // respecto al sistema global de coordenadas i1, i2, i3
	    double Xi = deformado[0][0];
	    double Yi = deformado[0][1];
	    double Zi = deformado[0][2];

	    double Xj = deformado[1][0];
	    double Yj = deformado[1][1];
	    double Zj = deformado[1][2];

	    double Xk = deformado[2][0];
	    double Yk = deformado[2][1];
	    double Zk = deformado[2][2];

	    // Unit vectors Local undeformed local coordinate axis
	    double m1 = sqrt((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi));
		double ve1[3] = {(xj-xi)/m1, (yj-yi)/m1, (zj-zi)/m1};
	    double m2 = sqrt((xk-xi)*(xk-xi) + (yk-yi)*(yk-yi) + (zk-zi)*(zk-zi));
		double ve4[3] = {(xk-xi)/m2, (yk-yi)/m2, (zk-zi)/m2};
	    double ve3[3];
		cross(ve1, ve4, ve3[0], ve3[1], ve3[2]);
		double m3 = norm(ve3);

		ve3[0] = ve3[0]/m3;
		ve3[1] = ve3[1]/m3;
		ve3[2] = ve3[2]/m3;

	    double ve2[3];
		cross(ve3,ve1,ve2[0], ve2[1], ve2[2]);

		// Formar la base e
		double e[3][3];
		e[0][0] = ve1[0];
		e[0][1] = ve1[1];
		e[0][2] = ve1[2];

		e[1][0] = ve2[0];
		e[1][1] = ve2[1];
		e[1][2] = ve2[2];

		e[2][0] = ve3[0];
		e[2][1] = ve3[1];
		e[2][2] = ve3[2];

		// Matriz de rotacion para la configuracion NO deformada
		// [r] hacia las coordenadas globales
	    double d1 = (xj-xi)/m1;
	    double d2 = (yj-yi)/m1;
	    double d3 = (zj-zi)/m1;
	    double e1 = (xk-xi)/m2;
	    double e2 = (yk-yi)/m2;
	    double e3 = (zk-zi)/m2;
	    double f1 = (d2*e3-d3*e2)/m3;
	    double f2 = (d3*e1-d1*e3)/m3;
	    double f3 = (d1*e2-d2*e1)/m3;
	    double g1 = f2*d3 - f3*d2;
	    double g2 = f3*d1 - f1*d3;
	    double g3 = f1*d2 - f2*d1;
	    double r[3][3] = {{d1, d2, d3},{g1, g2, g3},{f1, f2, f3}};

	    // Unit vectors Local deformed coordinate axis
	    double M1 = sqrt((Xj-Xi)*(Xj-Xi) + (Yj-Yi)*(Yj-Yi) + (Zj-Zi)*(Zj-Zi));
	    double E1[3] = {(Xj-Xi)/M1, (Yj-Yi)/M1, (Zj-Zi)/M1};
	    double M2 = sqrt((Xk-Xi)*(Xk-Xi) + (Yk-Yi)*(Yk-Yi) + (Zk-Zi)*(Zk-Zi));
	    double E4[4] = {(Xk-Xi)/M2, (Yk-Yi)/M2, (Zk-Zi)/M2};
		double E3[3];
		cross(E1,E4,E3[0], E3[1], E3[2]);
	    double M3 = norm(E3);
	    E3[0] = E3[0]/M3;
	    E3[1] = E3[1]/M3;
	    E3[2] = E3[2]/M3;
	    double E2[3];
		cross(E3,E1,E2[0],E2[1],E2[2]);

		// Formar la base E
		double E[3][3];
		E[0][0] = E1[0];
		E[0][1] = E1[1];
		E[0][2] = E1[2];

		E[1][0] = E2[0];
		E[1][1] = E2[1];
		E[1][2] = E2[2];

		E[2][0] = E3[0];
		E[2][1] = E3[1];
		E[2][2] = E3[2];

		// Matriz de rotacion para la configuracion deformada
		// [R] hacia la base global i
	    d1 = (Xj-Xi)/M1;
	    d2 = (Yj-Yi)/M1;
	    d3 = (Zj-Zi)/M1;
	    e1 = (Xk-Xi)/M2;
	    e2 = (Yk-Yi)/M2;
	    e3 = (Zk-Zi)/M2;
	    f1 = (d2*e3-d3*e2)/M3;
	    f2 = (d3*e1-d1*e3)/M3;
	    f3 = (d1*e2-d2*e1)/M3;
	    g1 = f2*d3 - f3*d2;
	    g2 = f3*d1 - f1*d3;
	    g3 = f1*d2 - f2*d1;
	    double R[3][3] = {{d1, d2, d3},{g1, g2, g3},{f1, f2, f3}};

		// Vectores posicion de los nodos en coordenadas locales estado deformado
	    // tiene origen en el nodo i deformado
	    double Xil = 0.0;
	    double Yil = 0.0;
	    double Zil = 0.0;
	    double Xjl = Xj - Xi;
	    double Yjl = Yj - Yi;
	    double Zjl = Zj - Zi;
	    double Xkl = Xk - Xi;
	    double Ykl = Yk - Yi;
	    double Zkl = Zk - Zi;


		// Transformar cada coordenada local del elemento deformado
	    // a la base comun i,j,k
		double Pi[3];
		double v1[3] = {Xil, Yil, Zil};
		MdotV(R, v1, Pi[0], Pi[1],Pi[2]);

		double Pj[3];
		double v2[3]={Xjl, Yjl, Zjl};
		MdotV(R,v2,Pj[0],Pj[1],Pj[2]);

		double Pk[3];
		double v3[3]={Xkl, Ykl, Zkl};
		MdotV(R,v3,Pk[0],Pk[1],Pk[2]);

		// Vectores posicion de los nodos en coordenadas locales estado NO deformado
	    // tiene origen en el nodo i NO deformado
	    double xil = 0.0;
	    double yil = 0.0;
	    double zil = 0.0;
	    double xjl = xj - xi;
	    double yjl = yj - yi ;
	    double zjl = zj - zi;
	    double xkl = xk - xi;
	    double ykl = yk - yi;
	    double zkl = zk - zi;

	    // Transformar cada coordenada local del elemento sin deformar a la base i,j,k
		double V1[3] = {xil, yil, zil};
		double pi[3];
		MdotV(r,V1,pi[0],pi[1],pi[2]);

		double V2[3] = {xjl, yjl, zjl};
		double pj[3];
		MdotV(r,V2,pj[0],pj[1],pj[2]);

		double V3[3] = {xkl, ykl, zkl};
		double pk[3];
		MdotV(r,V3,pk[0],pk[1],pk[2]);

		// Vectores de desplazamientos en el plano
	    double di[3];
		di[0] = Pi[0] - pi[0];
		di[1] = Pi[1] - pi[1];
		di[2] = Pi[2] - pi[2];

		double dj[3];
		dj[0] = Pj[0] - pj[0];
		dj[1] = Pj[1] - pj[1];
		dj[2] = Pj[2] - pj[2];

		double dk[3];
		dk[0] = Pk[0] - pk[0];
		dk[1] = Pk[1] - pk[1];
		dk[2] = Pk[2] - pk[2];

		// Llama rutina para calcular las fuerzas en el elemento deformado
		double ref[3][3], def[3][3];
		ref[0][0] = pi[0];
		ref[0][1] = pi[1];
		ref[0][2] = pi[2];

		ref[1][0] = pj[0];
		ref[1][1] = pj[1];
		ref[1][2] = pj[2];

		ref[2][0] = pk[0];
		ref[2][1] = pk[1];
		ref[2][2] = pk[2];

		def[0][0] = Pi[0];
		def[0][1] = Pi[1];
		def[0][2] = Pi[2];

		def[1][0] = Pj[0];
		def[1][1] = Pj[1];
		def[1][2] = Pj[2];

		def[2][0] = Pk[0];
		def[2][1] = Pk[1];
		def[2][2] = Pk[2];
		fuerzas(ref, def, nfuerzas, ks);

		// Transformar las fuerzas calculadas a la base original
		double tR[3][3];
		for(int i = 0; i<3; i++ )
			for(int j = 0; j<3; j++ )
			{
				tR[i][j] = R[j][i];
			}
		double nfn1[3], fn1[3] = {nfuerzas[0][0], nfuerzas[0][1], nfuerzas[0][2]};
		double nfn2[3], fn2[3] = {nfuerzas[1][0], nfuerzas[1][1], nfuerzas[1][2]};
		double nfn3[3], fn3[3] = {nfuerzas[2][0], nfuerzas[2][1], nfuerzas[2][2]};
		MdotV(tR, fn1, nfn1[0], nfn1[1], nfn1[2]);
		MdotV(tR, fn2, nfn2[0], nfn2[1], nfn2[2]);
		MdotV(tR, fn3, nfn3[0], nfn3[1], nfn3[2]);

		nfuerzas[0][0] = nfn1[0];
		nfuerzas[0][1] = nfn1[1];
		nfuerzas[0][2] = nfn1[2];

		nfuerzas[1][0] = nfn2[0];
		nfuerzas[1][1] = nfn2[1];
		nfuerzas[1][2] = nfn2[2];

		nfuerzas[2][0] = nfn3[0];
		nfuerzas[2][1] = nfn3[1];
		nfuerzas[2][2] = nfn3[2];
	}

#endif
