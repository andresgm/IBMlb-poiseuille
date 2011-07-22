#ifndef MESH_H
	#define MESH_H
	using namespace std;

	class mesh{

		private:
		    int id;
		    int nNodos;
		    int nCeldas;
			double cX, cY, cZ;
			double **vertex;
			double **velocidad;
			double **velocidad2;
			double **fuerza;
			double *area;
			int **faces;
			double **normalesPorNodo;
			double **normalesPorCara;
			double **carasPorNodo;
			double **angulosPorNodo;
			double **vecinosPorNodo;
			double *laplaceKg;
			double *laplaceKh;
			int *nodosProblema;
			double areaS;
			double volumenE;

		public:
			void setID(int ID){id=ID;}
			int getID(){return id;}
			double** darNodos(){return vertex;}
			int** darCeldas(){return faces;}
			int darNumeroNodos(){return nNodos;}
			int darNumeroCeldas(){return nCeldas;}
			int posicionNodo(double x, double y, double z);
			int guardarVTU(int t);
			int agregarNodo(double x, double y, double z);
			int agregarCelda(int a, int b, int c);
			int existeNodo(double x, double y, double z);
			void mesh_refine_tri4();
			void proyectarEsfera(double r);
			void proyectarElipsoide(double a, double b, double c);
			void proyectarRBC(double r);
			void moverCentro(double x, double y, double z);
			void rotarEstructura(double alpha, double phi, double theta);
			double darKgPromedioPorNodo(int nodo);
			double darKhPromedioPorNodo(int nodo);
			double darLaplaceKgPromedioPorNodo(int nodo);
			double darLaplaceKhPromedioPorNodo(int nodo);
			double darKgPorNodo(int nodo);
			double darKhPorNodo(int nodo);
			double darK1PorNodo(int nodo);
			double darK2PorNodo(int nodo);
			void darMeanCurvatureOperator(int nodo, double dmco[3]);
			void iniciarGeometria();
			void actualizarGeometria();
			bool contieneNodo(int nodo, int cara);
			void darNormalPromedio(int nodo, double normal[3]);
			void calcularCurvaturaGaussiana();
			double calcularAngulosPorNodo(int nodo, double angulos[7]);
			void darCarasPorNodo(int nodo, int caras[7]);
			void darNodosPorElemento(int cara, int nodos[3]);
			double darAreaAlrededorPorNodo(int nodo);
			double darAreaVoronoiPorNodo(int nodo);
			double darAreaVoronoiParcial(int nodoA, int nodoB);
			double darAreaBaricentricaPorNodo(int nodo);
			void darNodosVecinos(int nodo, int vecinos[7]);
			void darCarasSegunDosNodos(int nodoA, int nodoB, int caras[2]);
			double darLaplaceKg(int nodo);
			double darLaplaceKh(int nodo);
			void darNormalCara(int i, double normal[3]);
			void darPosNodo(int n, double pos[3]);
			void darFuerzaNodo(int n, double f[3]);
			void setVelocidad(int n, double ux, double uy, double uz);
			void moverNodos(double dt, double dx, int Y);
			void calcularFuerzasFEM(mesh referencia, double ks);
			void calcularFuerzasHelfrich(double kb);
			void actualizarNodos(double **);
			void calcularCambioArea(mesh ref);
			double calcularAreaSuperficial();
			double darAreaElemento(int i);
			double darVolumenElemento(int i);
			double calcularVolumen();
			void calcularFuerzasVolumen(double v0, double ka);
			void calcularFuerzaNeta(double fNeta[3]);
			double calcularEnergia();
			void calcularMomentoNeto(double fMomento[3]);
			void encontrarNodosProblema();
			bool esNodoProblema(int nodo);
			mesh();
			~mesh();
};
#endif
