#include<ilcplex/ilocplex.h>
#include<string>
#include <fcntl.h>
#include <JV/lap.cpp>
#include <JV/lap.h>
#include <deque>
#include <vector>
#define NUMTRIPS  500
#define NUMDEPOTS 4
#define NUM_STABILIZED 10
using namespace std;

ILOSTLBEGIN

struct Edge {
    int source;
    int grau;
    long int dest[NUMTRIPS];
    long long int weight[NUMTRIPS];
};
//
// Constantes
//
const unsigned int ITERLIMIT = 100000;                        // Limite de iterações
const unsigned int ZMIN = 0;                                  // Valor mínimo de redução
const unsigned short int ITERMIN = 20;                       // Número mínimo de iterações
const float OMEGAMIN = 0.7;                                  // Valor minimo de omega para fixação da coluna
const long int BIGM = 1e6;                                    // big-m inicial
const unsigned long int INFEASIBLEVALUE = 1e9;                // Valor marcado como infeasible
const float DEPOTSJVCOEF = 1.23;                              // Coeficiente do JV (número mínimo de veículos)
/*
 * Leitura do arquivo, transformando os dados nas variáveis correspondentes
 */
static void readdata(const char* filename, IloArray<IloNumArray>& costMatrix, unsigned short int& depots, unsigned int& nodes, IloNumArray& maxVehiclesPerDepot)
{
    ifstream in(filename);
//cout << "Teste" << endl;
    if (in)
    {
        in >> depots;
        in >> nodes;
        in >> maxVehiclesPerDepot;
        in >> costMatrix;
    }
    else
    {
        cerr << "Arquivo não encontrado..." << filename << endl;
    }
}

/*
 * Redução 1 do rohde, mantém apenas as n variáveis por linha e coluna com valores ótimos, em relação ao modelo maior
 */
void reducao1Rohde(IloArray<IloNumArray> costMatrix, int depots, int nodes, int matrixSize, vector< vector<bool> >& assignMatrix)
{
    //cout << "Aplicando redução 1 - Menores caminhos" << endl;
    unsigned int i, j, x;
    vector< vector<int> > reduceItens(matrixSize, vector<int>(matrixSize) );
    unsigned int reduction = ceil(sqrt(nodes)); // TODO: Revisar histograma de redução
    vector<unsigned int> reducei(reduction);
    vector<unsigned int> reducej(reduction);
    vector<int> indexj(reduction);
    vector<int> indexi(reduction);
    int posChangei;
    int posChangej;

    // Cria array de redução
    for (i=0; i<matrixSize; i++)
        for (j=0; j<matrixSize; j++)
            reduceItens[i][j] = 0;

    //
    // Inicia procedimento de redução
    //
    for (i=depots; i<matrixSize; i++)
    {
        // Reinicia variáveis
        for (x=0; x<reduction; x++)
        {
            reducej[x] = 1e6;
            reducei[x] = 1e6;
            indexj[x] = -1;
            indexi[x] = -1;
            posChangej = -1;
            posChangei = -1;
        }
        for (j=depots; j<matrixSize; j++)
        {
            posChangei = -1;
            posChangej = -1;

            // Redução das linhas
            if (costMatrix[i][j] != -1)
            {
                for (x=0; x<reduction; x++)
                {
                    if (costMatrix[i][j]<=reducej[x])
                    {
                        posChangej = x;
                        break;
                    }
                }
                if (posChangej!=-1)
                {
                    for (x=reduction-1; x>posChangej; x--)
                    {
                        reducej[x] = reducej[x-1];
                        indexj[x] = indexj[x-1];
                    }
                    reducej[posChangej] = costMatrix[i][j];
                    indexj[posChangej] = j;
                }
            }

            // Redução das colunas
            if (costMatrix[j][i] != -1)
            {
                for (x=0; x<reduction; x++)
                {
                    if (costMatrix[j][i]<=reducei[x])
                    {
                        posChangei = x;
                        break;
                    }
                }
                if (posChangei != -1)
                {
                    for (x=reduction-1; x>posChangei; x--)
                    {
                        reducei[x] = reducei[x-1];
                        indexi[x] = indexi[x-1];
                    }
                    reducei[posChangei] = costMatrix[j][i];
                    indexi[posChangei] = i;
                }
            }

        }
        for (x=0; x<reduction; x++)
        {
            if (indexj[x]!=-1)
                reduceItens[i][indexj[x]] = 1;
            if (indexi[x]!=-1)
                reduceItens[indexi[x]][i] = 1;
        }

    }
    for (i=depots; i<matrixSize; i++)
        for (j=depots; j<matrixSize; j++)
            if (reduceItens[i][j] == 1)
            {
                assignMatrix[i][j] = true;
            }
}

/*
 *
 */
void reducao2Rohde(IloArray<IloNumArray> costMatrix, int depots, int nodes, int matrixSize, vector< vector<bool> >& assignMatrix)
{
    //
    // JONKER VOLGENANT
    //
    #define COSTRANGE 1000.0
    #define PRINTCOST 0
    int bestlap, bestdim, depot;
    int dim, startdim, enddim;
    cost **assigncost, *u, *v, lapcost;
    cost **origAssignCost;
    row i, *colsol;
    col j, *rowsol;
    double runtime;
    unsigned int k = 0;
    unsigned int i_, j_;

    startdim = nodes;
    enddim = nodes+nodes;

    assigncost = new cost*[enddim];
    origAssignCost = new cost*[enddim];

    for (i = 0; i < enddim; i++)
    {
        assigncost[i] = new cost[enddim];
        origAssignCost[i] = new cost[enddim];
    }
    rowsol = new col[enddim];
    colsol = new row[enddim];
    u = new cost[enddim];
    v = new cost[enddim];
  
    for (depot = 0; depot<depots; depot++)
    {
        for (i = 0; i <nodes; i++)
        {
            for (j = 0; j <nodes; j++)
            {
                origAssignCost[i][j] = (cost)(costMatrix[i+depots][j+depots]==-1 ? INFEASIBLEVALUE : costMatrix[i+depots][j+depots]);
                origAssignCost[i+nodes][j] = (cost)(costMatrix[depot][j+depots]==-1 ? INFEASIBLEVALUE : costMatrix[depot][j+depots]);
                origAssignCost[j][i+nodes] = (cost)(costMatrix[j+depots][depot]==-1 ? INFEASIBLEVALUE : costMatrix[j+depots][depot]);
            }
        }
        
        bestlap = -1;
        dim = ceil(startdim*DEPOTSJVCOEF);
        int iterations = 0;
        bool dimrevert = false;
        bool best = false;
        for (;;)
        {
            bestlap = lapcost;  
            assigncost = origAssignCost;
            lapcost = lap(dim, assigncost, rowsol, colsol, u, v);
            if ((!dimrevert) && (bestlap == lapcost))
            {
                if (iterations<2)
                {
                    dimrevert = true;
                }
                else
                {
                    bestlap = lap(dim-1, assigncost, rowsol, colsol, u, v);
                    bestdim = dim-1;
                    break;
                }
            }
            if (dimrevert)
            {
                if (bestlap != lapcost)
                {
                    bestlap = lap(dim+1, assigncost, rowsol, colsol, u, v);
                    bestdim = dim+1;
                    break;
                }
                dim--;
            }
            else
            {
                dim++;
            }
            iterations++;
        }
        for (dim=0; dim<bestdim; dim++)
        {
            
            assignMatrix[dim>=nodes ? depot : dim+depots][rowsol[dim]>=nodes ? depot : rowsol[dim]+depots] = true;    
        }
    }
    delete[] assigncost;
    delete[] rowsol;
    delete[] colsol;
    delete[] u;
    delete[] v;
}

bool reducoesRohde(IloArray<IloNumArray>& costMatrix, int depots, int nodes, int matrixSize, bool reducao1 = true, bool reducao2 = true)
{
    unsigned int i, j;
    vector< vector<bool> >assignMatrix(matrixSize, vector<bool>(matrixSize));

    for (i=0; i<matrixSize; i++)
        for (j=0; j<matrixSize; j++)
            assignMatrix[i][j] = false;

    if (!reducao1 && !reducao2)
        return true;

    if (reducao1 && !reducao2)
    {
        reducao1Rohde(costMatrix, depots, nodes, matrixSize, assignMatrix);
        for (i=depots; i<matrixSize; i++)
            for (j=depots; j<matrixSize; j++)
                if (assignMatrix[i][j] == false)
                {
                    costMatrix[i][j] = -1;
                }
    }
    if (reducao2 && !reducao1)
    {
        reducao2Rohde(costMatrix, depots, nodes, matrixSize, assignMatrix);
        for (i=0; i<matrixSize; i++)
            for (j=0; j<matrixSize; j++)
                if (assignMatrix[i][j] == false)
                {
                    costMatrix[i][j] = -1;
                }
    }
    if (reducao1 && reducao2)
    {
        reducao1Rohde(costMatrix, depots, nodes, matrixSize, assignMatrix);
        reducao2Rohde(costMatrix, depots, nodes, matrixSize, assignMatrix);
        for (i=0; i<matrixSize; i++)
            for (j=0; j<matrixSize; j++)
                if (assignMatrix[i][j] == false)
                {
                    costMatrix[i][j] = -1;
                }
    }

    return true;
}


bool find(deque<int> &List, int k) 
{
    for(deque<int>::iterator it = List.begin(); it != List.end(); it++) 
        if(*it == k) return true;
    return false;
}



int SLF(vector<int> &pred, IloArray<IloNumArray>& costMatrix, long int nodes, short int depots, unsigned short int depot, IloNumArray& pi, IloNumArray& sigma) {

    long long int distances[nodes+2];
    distances[nodes] = 0;
    pred[nodes] = -1;
    distances[nodes+1] = INFEASIBLEVALUE;
    pred[nodes+1] = -1;
    vector<Edge> Arestas(nodes+1);

    for (int i=0; i<nodes; i++)
    {
        // Alimenta as viagens entre os nodos
        for (int j=0; j<nodes; j++)
        {
            if (costMatrix[i+depots][j+depots] == -1)
                continue;    

            Arestas[i].dest[Arestas[i].grau] = j;
            Arestas[i].weight[Arestas[i].grau] = costMatrix[i+depots][j+NUMDEPOTS]-pi[i];
            Arestas[i].grau++;
        }

        // Adicionando a viagem de pull-out
        if (costMatrix[depot][i+depots] != -1)
        {
            Arestas[nodes].dest[Arestas[nodes].grau] = i;
            Arestas[nodes].weight[Arestas[nodes].grau] = costMatrix[depot][i+depots]-sigma[depot];
            Arestas[nodes].grau++;
        }

        // Adicionando viagem de pull-in
        if (costMatrix[i+depots][depot] != -1)
        {
            Arestas[i].dest[Arestas[i].grau] = nodes+1;
            Arestas[i].weight[Arestas[i].grau] = costMatrix[i+depots][depot]-pi[i]; 
            Arestas[i].grau++;
        }
        distances[i] = INFEASIBLEVALUE;
        pred[i] = -1;
    }
     deque<int> List;
     List.push_front(nodes);
     int i;
     while(!List.empty()) 
     {
         i = List.front();
         List.pop_front();
         for(int k = 0; k < Arestas[i].grau; k++) 
         {
            if(distances[Arestas[i].dest[k]] > distances[i] + Arestas[i].weight[k]) 
            {
                distances[Arestas[i].dest[k]] = distances[i] + Arestas[i].weight[k];
                pred[Arestas[i].dest[k]] = i;
                if(!find(List, Arestas[i].dest[k])) 
                    List.push_front(Arestas[i].dest[k]);
            }
         }
     }
     return distances[nodes+1];
}

/*
 * Feasibility test 2 - Para todos os omegas, testa se foi encontrada a melhor solução
 * se foi encontrada, retorna true, caso contrário, false.
 */
static bool feasibility2Test(int nodes, IloNumArray masterData)
{
    unsigned int j;
    for (j=0; j<nodes; j++)
        if (masterData[j] != 0)
            return false;
    return true;
}

/*
 * Feasibility test 6 - Para todos os deltas, analisa se foi encontrada a melhor solução
 * se foi encontrado, retorna true, caso contrário, false
 */ 
static bool feasibility6Test(int nodes, IloNumArray delta)
{
    unsigned int j;
    for (j=0; j<nodes; j++)
    {
        if (delta[j] == 0)
        {
//            cout << "Infeasible solution, trying again" << endl;
            return false;
        }
    }
//    cout << "Feasible" << endl;
    return true;
}

/*
 * Função principal
 */
int main(int argc, char **argv)
{
    IloEnv env;
    try
    {
        //
        // STEP 0: Initialization
        //
        bool lb = false;
        unsigned short int depots;                              // Número de garagens
        unsigned int nodes;                                     // Número de nodos
        IloNumArray maxVehiclesPerDepot(env);                   // Matriz com o Número máximo de carros por garagem
        int masterValues[ITERLIMIT];                                    
        IloArray<IloNumArray> costMatrix(env);                              // Matriz de custos entre garagens e os pontos "i" e "j"
        if (argc > 1)
            readdata(argv[1], costMatrix, depots, nodes, maxVehiclesPerDepot);
        else
            readdata("instancias/pepin/m4n500s0.inp", costMatrix, depots, nodes, maxVehiclesPerDepot);
        unsigned short int matrixSize = depots+nodes;           // Tamanho da matriz
        unsigned long int i;                                    // Contador "i" -> Saída nodo
        unsigned long int j;                                    // Contador "j" -> Chegada nodo
        unsigned  int k;                                        // Contador "k" -> Garagens

        // Reduções do Rohde 
        // - penultimo parâmetro -> True para reducao 1 on e false para reducao 1 off
        // - ultimo parâmetro -> True para reducai 2 on e false para reducao 2 off
        //reducoesRohde(costMatrix, depots, nodes, matrixSize, true, true);

        //
        // Variáveis de auxílio
        //
        //
		// Vari�veis de aux�lio
		//
		unsigned short int p = 0;                                   // Contador "p" -> Caminhos
		unsigned short int p_ = 0;                                  // Contador "p_" -> Controle extra dos caminhos
		unsigned int p__ = 0;                                       // 
		unsigned short int pAnt = 0;                                // Contador "pAnt" -> Controle dos caminhos adicionados ao problema
		bool solution = false;                                          // Vari�vel de controle para teste se houve ou n�o solu��o
		vector<bool> checkp(nodes);                                     // Vari�vel que "marca" os caminhos utilizados no RMP
		unsigned int pos = 0;                                       // posi��o para processamento do predecessor para a[j][p]
		vector<int> pred(matrixSize);                           // Valores predecessores
		bool testFeasibility6 = false;                                  // Booleana de testes do STEP 6
		bool checkOmega = false;                                        // Verifica��o se o Omega foi ou n�o alterado
		float lastOmega = 0;                                    // �ltimo valor de omega entre os maiores
		unsigned short int omegasSelected[ITERLIM];            // Omegas selecionados
		IloIntArray omega(env, ITERLIM);                       // Valores de omega com depot
		bool continueGc = true;                                // Booleana de controle para gera��o de colunas
		unsigned int iteracoes = 1;                            // N�mero de itera��es
		unsigned int nO = 0;                                   // "�ltima itera��o com rounding"
		unsigned int lastp = 0;                                // �ltimo p da rodada 
		bool feasibility2 = false;

		//
		// Cria��o do problema mestre 
		//
		IloNumArray c(env, ITERLIM);
		IloNumVarArray y(env, ITERLIM);
		IloArray<IloNumArray> a(env, ITERLIM);
		IloNumVarArray delta(env, matrixSize);
		IloNumArray masterData(env);
		IloModel masterModel(env);
		IloObjective pathsSelect = IloAdd(masterModel, IloMinimize(env));
		IloRangeArray pi = IloAdd(masterModel, IloRangeArray(env, nodes, 1, IloInfinity));
		IloRangeArray sigma = IloAdd(masterModel, IloRangeArray(env, 0, maxVehiclesPerDepot));
		IloNumVarArray deltapi(env, nodes, 0, IloInfinity, ILOFLOAT);
		IloNumVarArray deltasigma(env, depots, 0, IloInfinity, ILOFLOAT);

		for (int j=0; j<nodes; j++)
		{
			deltapi[j] = IloNumVar(pathsSelect(bigM)+pi[j](1));
		}

		for (int k=0; k<depots; k++)
		{
			deltasigma[k] = IloNumVar(sigma[k](1));
		}

		IloCplex master(masterModel);
		//master.setParam(IloCplex::RootAlg, CPX_ALG_BARRIER);
		//master.setParam(IloCplex::NodeAlg, CPX_ALG_BARRIER);
		//master.setParam(IloCplex::MIPEmphasis, 2);
		//master.setParam(IloCplex::WorkMem ,8000.0);
		master.setOut(env.getNullStream());
		IloNumArray piDuals(env);
		IloNumArray sigmaDuals(env);
		p = 0;

		// var para tomada de tempo

		time(&before);
		int stabilized = 0;
		masterValues[0] = NUMTRIPS * bigM;
		while (continueGc)
		{
			pAnt = p;			
			testFeasibility6 = false;
			continueGc = false;
			master.solve(); 
			cout << iteracoes << " - " << master.getObjValue() << endl;
			if(masterValues[iteracoes-1] ==  master.getObjValue()) {
				stabilized++; 
			} else {
				stabilized = 0;
			}
			masterValues[iteracoes] = master.getObjValue();
			iteracoes++;
			if(stabilized > NUM_STABILIZED) continueGc = false;
			else {
				master.getValues(deltapi, masterData);
				// Envio ao CPLEX e obten��o da solu��o da formula��o
				// Se foi encontrada a solu��o, verifica se foi rodado um m�nimo 
				// de itera��es antes entre o in�cio ou o �ltimo arredondamento
				// para habilitar novamente o processo de arredondamento
				if (feasibility2Test(nodes, masterData))
				{
					if ((iteracoes-nO>ITER_MIN) && (masterValues[iteracoes-2] - masterValues[iteracoes-1] <= ZMIN))
						testFeasibility6=true;
				}
				// Se n�o for arredondar, ent�o, passa pelo processo de gera��o de novas colunas
				if (!testFeasibility6)
				{
					// Cria as viagens entre os nodos (sempre ser�o iguais para todas as garagens) - guarda valor de arestas
					master.getDuals(piDuals, pi);
					master.getDuals(sigmaDuals, sigma);
					for (int k=0; k<depots; k++)
					{
						// Roda o Bellman-ford
						if (SLF(pred, costMatrix, nodes, depots, k, piDuals, sigmaDuals)<0)
						{
							a[p] = IloNumArray(env, nodes);
							// - Ajusta os valores para verificar se o valor do shortest path 
							//   pode adicionar um caminho de redu��o ao problema mestre
							// - Se o custo total for menor que zero (retorno do shortest path),
							//   significa que o valor pode ser utilizado para reduzir a fun��o objetivo do RMP
							c[p] = 0;
							pos = pred[nodes+1];
							c[p]+= costMatrix[pos+depots][k];
							omega[p] = k;
							while (pos != nodes)
							{
								if (pred[pos] == nodes)
								{
									c[p]+= costMatrix[k][pos+depots];
									a[p][pos] = 1;
								}
								else
								{
									c[p]+= costMatrix[pred[pos]+depots][pos+depots];
									a[p][pos] = 1;
								}
								pos = pred[pos];
							}
							y[p] = IloNumVar(pathsSelect(c[p])+pi(a[p])+sigma[k](1));
							y[p].setBounds(0, 1);
							omegasSelected[p] = false;
							p++;
							
							continueGc = true;
							master.solve(); 
							master.getDuals(piDuals, pi);
							master.getDuals(sigmaDuals, sigma);
							
							//cout << " - " << master.getObjValue() << endl;
							/*if(masterValues[iteracoes-1] ==  master.getObjValue()) {
								stabilized++; 
							} else {
								stabilized = 0;
							}*/
							//masterValues[iteracoes] = master.getObjValue();
							//iteracoes++;
							/*if(stabilized > NUM_STABILIZED) {
								continueGc = false;
								break;
							}*/
						}
					}
				}
			}
			if (!continueGc)
			{
				cout << "Stabilized" << endl;
				testFeasibility6 = true;
			}
			if (testFeasibility6)
			{
				if (!feasibility6Test(nodes, masterData))
				{
					cout << "No solution found, sorry :/" << endl;
					env.end();
					time(&after);
					seconds = difftime(after,before);
					cout << "tempo total = " << seconds << " segundos" << endl;
					return 1; 
				}
				if (!lb)
				{
					for (int j=0; j<nodes; j++)
						pi[j].setBounds(1,1);
					nO = iteracoes;
					lb = true;
					//p__ = p;
				}
				else
				{
					if (iteracoes>ITER_MIN)
					{
						solution = true;
						//
						// Verifica se todos os resultados do problema s�o inteiros
						//
						cout << "Integrability" << endl;
						for (p_=0; (p_<pAnt); p_++)
						{
							if ((master.getValue(y[p_]) != 0) && (master.getValue(y[p_]) != 1))
							{
								cout << "No best solution found, yet..." << endl;
								// Se n�o forem, indica que n�o houve solu��o
								solution = false;
								break;
							}
						}


						// Se a integrabilidade est� correta, mostra as informa��es e encerra o processo
						if (solution)
						{
							// Conta o n�mero de ve�culos
							int nVehicles = 0;
							//for (p_=p__; ((p_<p) && (solution)); p_++)
							for (p_=0; ((p_<p) && (solution)); p_++)
							{
								if (master.getValue(y[p_]) == 0)
									continue;
								nVehicles++;
								//cout << "Path " << p_ << endl;
								//cout << "Depot " << omega[p_] << endl;
								//for (int i=0; i<nodes; i++)
								//{
								//	if (a[p_][i] == 1)
								//		cout << i << " - ";
								//}

//								cout << endl;
							}
							// Mostra a solu��o
							cout << "Feasible solution found" << endl;
							cout << "Number of vehicles: " << nVehicles << endl;
							cout << "Objective value: " << master.getObjValue() << endl;
							env.end();
							time(&after);
							seconds = difftime(after,before);
							cout << "tempo total = " << seconds << " segundos" << endl;
							return 0;    
						}
					}
					//
					// ============================= STEP 7: END INTEGRABILITY TEST ======================================
					//
					//
					// ==================================== STEP 8: ROUNDING =============================================
					//
					//
					// Caso n�o tenha ocorrido a integrabilidade
					// ent�o precisamos de mais rodadas, com perturba��o
					//
					cout << "Rounding" << endl;
					checkOmega = false;
					lastOmega = 0;
					for (p_=0; p_<pAnt; p_++)
					//for (p_=p__; p_<pAnt; p_++)
					{
						//
						// Se o "p" j� foi selecionado anteriormente, pula
						//
						if (omegasSelected[p_])
							continue;

						//
						// Se for maior que o "omegaMin", ent�o marca como um caminho 
						// a ser utilizado
						//
						if (master.getValue(y[p_])>= OMEGA_MIN && master.getValue(y[p_]) != 1)
						{
							y[p_].setLB(1);
							omegasSelected[p_] = true;
							master.solveFixed();
							checkOmega = true;
						}

						if (lastOmega<master.getValue(y[p_]) && master.getValue(y[p_]) != 1)
						{
							lastOmega = master.getValue(y[p_]);
							lastp = p_;
						}
					}
					//
					// Se n�o achou, pega o maior da �ltima rodada
					//
					if (!checkOmega)
					{
						omegasSelected[lastp] = true;
						masterModel.add(y[lastp] == 1);
					}

					//rounding++;
					nO = iteracoes;
					feasibility2 = false;
					//
					// ============================== STEP 8 END: ROUNDING =================================
					//
				}

			}
			continueGc = true;
			// Guarda o valor da �ltima itera��o
		}
		// Checa limite de itera��es
		if (p >= ITERLIM)
		{
			cout << "Limite de itera��es excedido, resultado obtido: " << masterValues[iteracoes-1] << endl;
		} 
	}
    }
    //
    // Exceções (que nunca acontecem)...
    //
    catch (IloException& ex)
    {
        cerr << "Error: " << ex << endl;
    }
    catch (...)
    {
        cerr << "Error: " << endl;
    }
    env.end();
    return 0;
}
