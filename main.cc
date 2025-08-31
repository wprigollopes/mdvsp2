#include<ilcplex/ilocplex.h>
#include<string>
#include <fcntl.h>
#include <JV/lap.cpp>
#include <JV/lap.h>
#include <deque>
#include <vector>
#define NUMTRIPS 1500
#define NUMDEPOTS 8
#define NUM_STABILIZED 1e6
#define bigM 1e6
#define INFEASIBLEVALUE 1e8
using namespace std;

ILOSTLBEGIN
// teste denis
struct Edge {
    int source;
    int grau;
    long int dest[NUMTRIPS];
    long long int weight[NUMTRIPS];
};
//
// Constantes
//
const unsigned int ITERLIM = 100000;      // Limite de iterações
const unsigned int ZMIN = 0;              // Valor mínimo de redução
const unsigned short int ITER_MIN = 30;  // Número mínimo de iterações
const float OMEGA_MIN = 0.9;              // Valor minimo de omega para fixação da coluna
const long int BIGM = 1e6;                // big-m inicial
const int CGTYPE = 1;                     // 1: Geração com relaxamento; 2: Geração sem relaxamento.

//
// Leitura do arquivo, transformando os dados nas variáveis correspondentes
//
static void readdata(const char* filename, IloArray<IloNumArray>& costMatrix, unsigned short int& depots, unsigned int& nodes, IloNumArray& maxVehiclesPerDepot)
{
    ifstream in(filename);
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

 void reducao1Denis(IloArray<IloNumArray> costMatrix, int depots, int nodes, int matrixSize, vector< vector<bool> >& assignMatrix)
{
    unsigned int i, j, x;
    vector< vector<int> > reduceItens(matrixSize, vector<int>(matrixSize) );
    unsigned int reduction = ceil(sqrt((float)nodes));
    vector<unsigned int> reducei(reduction);
    vector<unsigned int> reducej(reduction);
    vector<int> indexj(reduction);
    vector<int> indexi(reduction);
    int posChangei;
    int posChangej;

    // Cria array de reducao
    for (i=0; i<matrixSize; i++)
    {
        for (j=0; j<matrixSize; j++)
        {
            reduceItens[i][j] = 0;
        }
    }

    // Inicia procedimento de reducao
    for (i=depots; i<matrixSize; i++)
    {
        // Reinicia variaveis
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

            // Reducao das linhas
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

            // Reducao das colunas
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
    {
        for (j=depots; j<matrixSize; j++)
        {
            if (reduceItens[i][j] == 1)
            {
                assignMatrix[i][j] = true;
            }
        }
    }
}

void reducao2Denis(IloArray<IloNumArray> costMatrix, int depots, int nodes, int matrixSize, vector<vector<bool> >& assignMatrix, vector<  vector<vector<bool> > >& assignSImatrix)
{
    // JONKER VOLGENANT
    #define COSTRANGE 1000.0
    #define PRINTCOST 0
    int depot;
    int dim, enddim;
    cost *u, *v;
    int lapcostvalue = 0;
    cost **origAssignCost;
    col j, *rowsol;
    row i, *colsol;

    //col *colsol1;
    //row *rowsol1;
    //cost *u1, *v1;
    enddim = nodes+nodes;

    origAssignCost = new cost*[enddim];
    
    for (unsigned int i = 0; i < enddim; i++)
    {
        origAssignCost[i] = new cost[enddim];
    }
    
     for (int i = 0; i < enddim; i++)
    {
        origAssignCost[i] = new cost[enddim];
    }
    
    rowsol = new col[enddim];
    colsol = new row[enddim];
    u = new cost[enddim];
    v = new cost[enddim];

		 
    for (depot = 0; depot<depots; depot++)
    {
        for (int i = 0; i <nodes; i++)
        {
            for (j = 0; j <nodes; j++)
            {
                origAssignCost[i][j] = (cost)(costMatrix[i+depots][j+depots]==-1 ? INFEASIBLEVALUE : costMatrix[i+depots][j+depots]);
                origAssignCost[i+nodes][j] = (cost)(costMatrix[depot][j+depots]==-1 ? INFEASIBLEVALUE : costMatrix[depot][j+depots]);
                // origAssignCost[i][j+nodes] = (cost)(costMatrix[j+depots][depot]==-1 ? INFEASIBLEVALUE : costMatrix[j+depots][depot]);
		origAssignCost[j][i+nodes] = (cost)(costMatrix[j+depots][depot]==-1 ? INFEASIBLEVALUE : costMatrix[j+depots][depot]);
		origAssignCost[i+nodes][j+nodes] = 0;
            }
        }
        
	 /*	
	 for(int i = 0; i < enddim; i++) {
		for(int j = 0; j < enddim; j++) {
			cout << origAssignCost[i][j] << " ";
			// printf{"%d ", origAssignCost[i][j]);
		}
			cout << endl;
		// printf("\n");
	}
	*/	
        	
	lapcostvalue = lap(enddim, origAssignCost, rowsol, colsol, u, v);
		

	/*for (int i = 0; i < enddim ; i++)
	{
		cout << "linha = " << i << " assigned to column - " << rowsol[i] << endl;
		// cout << "column = " << i << "assigned to column - " << colsol[i] << endl; 
				
	}
	cout << endl;
	*/
	for (dim=0; dim<enddim; dim++)
        {
             assignSImatrix[depot][dim>=nodes ? depot : dim+depots][rowsol[dim]>=nodes ? depot : rowsol[dim]+depots] = true; 
             assignMatrix[dim>=nodes ? depot : dim+depots][rowsol[dim]>=nodes ? depot : rowsol[dim]+depots] = true;   
        }

	assignSImatrix[depot][depot][depot] = false;
	assignMatrix[depot][depot] = false;

	/*
	for(int i = 0; i < nodes+depots; i++) {
		for(int j = 0; j < nodes+depots; j++) {
			cout << assignMatrix[i][j] << " ";
			// printf{"%d ", origAssignCost[i][j]);
		}
		cout << endl;
		// printf("\n");
	}
	*/
    }
    // delete[] assigncost;
    delete[] origAssignCost;
    delete[] rowsol;
    delete[] colsol;
    delete[] u;
    delete[] v;         
}

void reducao3Denis(IloArray<IloNumArray> costMatrix, int depots, int nodes, int matrixSize, vector<vector<bool> >& assignMatrix, IloNumArray maxVehiclesPerDepot)
{
        IloEnv env;
        IloModel relax_MDVSP(env);
        IloCplex reduction3(relax_MDVSP);
        IloRangeArray pi = IloAdd(relax_MDVSP, IloRangeArray(env, nodes, 1, 1));
        IloRangeArray sigma1 = IloAdd(relax_MDVSP, IloRangeArray(env, 0, maxVehiclesPerDepot));
        IloRangeArray sigma2 = IloAdd(relax_MDVSP, IloRangeArray(env, 0, maxVehiclesPerDepot));
        
        IloArray<IloNumVarArray> x(env, matrixSize);

        for (unsigned int k=0; k< matrixSize; k++)
        {
            x[k] = IloNumVarArray(env,matrixSize);
        }

        for (unsigned int i=0; i<matrixSize; i++)
        {
            for (unsigned int j=0; j<matrixSize; j++)
            {
                if(i < depots && j >= depots) 
                {
                    x[i][j] = IloNumVar(pi[j-depots](1)+sigma1[i](1));
                }
                else if(i >= depots && j < depots) {
                    x[i][j] = IloNumVar(sigma2[j](1));
                    
                } else if(i >= depots && j >= depots) {
                    x[i][j] = IloNumVar(pi[j-depots](1));
                } else {
                    x[i][j] = IloNumVar(env);
                }
            }
        }
        
        for (int i = depots; i < matrixSize; i++)
        {
            IloExpr fluxo1(env);
            IloExpr fluxo2(env);
            for (int j = 0; j < matrixSize; j++)
            {
                    fluxo1 += x[j][i];
                    fluxo2 += x[i][j];
            }
            relax_MDVSP.add(fluxo1-fluxo2 == 0);
            fluxo1.end();
            fluxo2.end();
        }
        
        IloExpr terms_fo(env);

        for (unsigned int i = 0; i < matrixSize; i++)
        {
            for (unsigned int j = 0; j < matrixSize; j++)
            {
                    
                    if (costMatrix[i][j] == -1)
                    {
                        terms_fo += INFEASIBLEVALUE * x[i][j];
                    }
                    else terms_fo += costMatrix[i][j] * x[i][j];
            }
        }
        
            
    relax_MDVSP.add(IloMinimize(env,terms_fo));    

    if (reduction3.solve())
    {
        for (unsigned int i=0; i<matrixSize; i++)
        {
            for (unsigned int j=0; j<matrixSize; j++)
            {
                
                if(reduction3.getValue(x[i][j]) > 0) 
                {
                    assignMatrix[i][j] = true;
		    // int denis = getchar();
		    // cout << "i=" << i << "j=" << j << endl;
                }
            }
        }
    }
    reduction3.end();
    relax_MDVSP.end();
    env.end();
}


bool reducoesDenis(IloArray<IloNumArray>& costMatrix, int depots, int nodes, int matrixSize, IloNumArray maxVehiclesPerDepot, vector<vector<int> >& predSI, vector<  vector< vector<bool> > >& assignSImatrix, bool reducao1 = true, bool reducao2 = true, bool reducao3 = true, bool SI = true)
{
    cout << "Redução 1 - Critérios estatísticos: " << (reducao1 ? "Sim" : "Não") << endl;
    cout << "Redução 2 - Jonker Volgenant: " << (reducao2 ? "Sim" : "Não") << endl;
    cout << "Redução 3 - Min Cost Max Flow: " << (reducao3 ? "Sim" : "Não") << endl;
    //cout << "==============================================" << endl;
    
    if (!reducao1 && !reducao2 && !reducao3 && !SI) return true;
    
     unsigned int i, j;
      
    vector< vector<bool> >assignMatrix(matrixSize);

    for (i=0; i<matrixSize; i++)
     	assignMatrix[i] = vector<bool>(matrixSize);
     	
    for (i=0; i<matrixSize; i++)
        for (j=0; j<matrixSize; j++)
            assignMatrix[i][j] = false;

      
    if (SI){
    	   
    	   reducao2Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, assignSImatrix);
    	   for (int k=0; k<depots; k++) {
		  for(int j = 0; j < matrixSize; j++) {
		    for(int i = 0; i < matrixSize; i++) {
		        if(assignSImatrix[k][j][i])
		            predSI[k][i] = j;
		    }
	          }	
  		}
    
    }

    if (reducao1 && reducao2 && !reducao3)
    {
        reducao1Denis(costMatrix, depots, nodes, matrixSize, assignMatrix);
	if (!SI) reducao2Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, assignSImatrix);
    }
    if (reducao2 && !reducao1 && !reducao3)
    {
        if (!SI) reducao2Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, assignSImatrix);
    }
    
    if (reducao1 && reducao2 && reducao3)
    {
        reducao1Denis(costMatrix, depots, nodes, matrixSize, assignMatrix);
        if (!SI) reducao2Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, assignSImatrix);
        reducao3Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, maxVehiclesPerDepot);
    }

    if (reducao2 && reducao3 && !reducao1)
    {
    	if (!SI) reducao2Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, assignSImatrix);
        reducao3Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, maxVehiclesPerDepot);
    }

    if (!reducao2 && reducao3 && reducao1)
    {
        reducao1Denis(costMatrix, depots, nodes, matrixSize, assignMatrix);
        reducao3Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, maxVehiclesPerDepot);
    }

    if (!reducao2 && reducao3 && !reducao1)
    {
        reducao3Denis(costMatrix, depots, nodes, matrixSize, assignMatrix, maxVehiclesPerDepot);
    }

    if (reducao1 || reducao2 || reducao3){
    	for (i=depots; i<matrixSize; i++)
        {
            for (j=depots; j<matrixSize; j++)
            {
                if (assignMatrix[i][j] == false)
                {
                    costMatrix[i][j] = -1;
                }
            }
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

    long int distances[nodes+2];
    distances[nodes] = 0;
    pred[nodes] = -1;
    distances[nodes+1] = INFEASIBLEVALUE;
    pred[nodes+1] = -1;
    vector<Edge> Arestas(nodes+1);

    for (unsigned int i=0; i<nodes; i++)
    {
        // Alimenta as viagens entre os nodos
        for (unsigned int j=0; j<nodes; j++)
        {
            if (costMatrix[i+depots][j+depots] == -1)
                continue;    

            Arestas[i].dest[Arestas[i].grau] = j;
            Arestas[i].weight[Arestas[i].grau] = costMatrix[i+depots][j+depots]-pi[i];
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
     unsigned int i;
     while(!List.empty()) 
     {
         i = List.front();
         List.pop_front();
         for(unsigned int k = 0; k < Arestas[i].grau; k++) 
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

//
// Feasibility test 2 - Para todos os omegas, testa se foi encontrada a melhor solução
// se foi encontrada, retorna true, caso contrário, false.
//
bool feasibility2Test(int nodes, IloNumArray masterData)
{
    unsigned int j;
    for (j=0; j<nodes; j++)
        if (masterData[j] != 0) {
            return false;
        }
    return true;
}

//
// Feasibility test 6 - Para todos os deltas, verifica se foi encontrado algum com valor diferente
// de zero, se foi encontrado, retorna false, caso contrário, retorna true.
//
bool feasibility6Test(int nodes, IloNumArray delta)
{
    unsigned int j;
    for (j=0; j<nodes; j++)
    {
        if (delta[j] != 0)
        {
            cout << "Infeasible solution, trying again" << endl;
            return false;
        }
    }
    return true;
}

// Tirei IloIntArray omegadepois de costMatrixSI

bool constroiPath(int i, int k, int p, vector< vector<int> > pred, IloArray<IloNumArray> costMatrixSI, int depots, IloArray<IloNumArray>& a, int& cost) {
    cost += costMatrixSI[i][k];
    //cout << i-depots << " - ";
    int pos = i;
    // omega[p] = k;
    
    while (pos != -1)
    {  
	//cout << pred[k][pos] << endl;
	if (pred[k][pos] < depots && pred[k][pos] != -1)
        {
            if (costMatrixSI[k][pos] == -1)
            {
		//cout << "to retornando false" << endl;                
		return false;
            } else {
                cost += costMatrixSI[k][pos];
                a[p][pos-depots] = 1;
		//cout << "depot " << k << endl;
                return true;
            }
        } else if (pred[k][pos] == -1) {
		cout << "Caiu no caso de erro" << endl;
		return false;
        }
        else
        {
            if (costMatrixSI[pred[k][pos]][pos] == -1)
            {
                return false;
            } else {
                cost += costMatrixSI[pred[k][pos]][pos];
                a[p][pos-depots] = 1;
		pos = pred[k][pos];
		//cout << pos-depots << " - ";
            }
        }
    }
    return true;
}

//
// Função principal
//
int main(int argc, char **argv)
{
    IloEnv env;
    try
    {
        //
        // STEP 0: Initialization
        //
        bool lb = (CGTYPE == 1) ? false : true;
        unsigned short int depots;                              // Número de garagens
        unsigned int nodes;                                     // Número de nodos
        IloNumArray maxVehiclesPerDepot(env);                   // Matriz com o Número máximo de carros por garagem
        int masterValues[ITERLIM];                              // Controle de valores mestre      
        IloArray<IloNumArray> costMatrix(env);                  // Matriz de custos entre garagens e os pontos "i" e "j"
        
        if (argc > 1)
        {
            readdata(argv[1], costMatrix, depots, nodes, maxVehiclesPerDepot);
            cout << argv[1] << endl;
        }
        else
        {
            readdata("instancias/pepin/m4n500s0.inp", costMatrix, depots, nodes, maxVehiclesPerDepot);
            cout << "instancias/pepin/m4n500s0.inp" << endl;
        }
        
        unsigned short int matrixSize = depots+nodes;           // Tamanho da matriz
        
            	
        	
	cout << endl << endl << " ========================================= " << endl;
	cout << (argc > 1 ? argv[1] : "---") << endl <<
         "Depots: " << depots << endl << 
         "Nodes: " << nodes << endl << 
         "OmegaMin: " << OMEGA_MIN << endl << 
         "IterMin: " <<  ITER_MIN << endl <<
         "Z min: " << ZMIN << endl << 
         "Forma de geração: " << ((CGTYPE == 1) ? "Com relaxamento" : "Sem relaxamento") << endl;
	
	
	// Variaveis de auxilio
        unsigned short int p = 0;                             // Contador "p" -> Caminhos
        unsigned short int p_ = 0;                            // Contador "p_" -> Controle extra dos caminhos
        unsigned int p__ = 0;                                 // 
        unsigned short int pAnt = 0;                          // Contador "pAnt" -> Controle dos caminhos adicionados ao problema
        bool solution = false;                                // Variavel de controle para teste se houve ou nao solucao
        vector<bool> checkp(nodes);                           // Variavel que "marca" os caminhos utilizados no RMP
        unsigned int pos = 0;                                 // Posicao para processamento do predecessor para a[j][p]
        vector<int> pred(matrixSize);                         // Valores predecessores
        vector<int> actPred(matrixSize);                      // Valores predecessores - melhor coluna
        bool testFeasibility6 = false;                        // Booleana de testes do STEP 6
        bool checkOmega = false;                              // Verificacao se o Omega foi ou nao alterado
        float lastOmega = 0;                                  // Ultimo valor de omega entre os maiores
        unsigned short int omegasSelected[ITERLIM];           // Omegas selecionados
        IloIntArray omega(env, ITERLIM);                      // Valores de omega com depot
        bool continueGc = true;                               // Booleana de controle para geracao de colunas
        unsigned int iteracoes = 1;                           // Numero de iteracoes
        unsigned int nO = 0;                                  // Ultima iteracao com rounding
        unsigned int lastp = 0;                               // Ultimo p da rodada 
        bool feasibility2 = false;                            //

        // Criacao do problema mestre 
        IloNumArray c(env, ITERLIM);
        IloNumVarArray y(env, ITERLIM);
        IloArray<IloNumArray> a(env, ITERLIM);
       // IloNumVarArray delta(env, matrixSize);
        IloNumArray masterData(env);
        IloModel masterModel(env);
        IloObjective pathsSelect = IloAdd(masterModel, IloMinimize(env));
        IloRangeArray pi = IloAdd(masterModel, IloRangeArray(env, nodes, 1, (CGTYPE == 1) ? IloInfinity : 1));
        IloRangeArray sigma = IloAdd(masterModel, IloRangeArray(env, 0, maxVehiclesPerDepot));
	
	
	IloNumVarArray deltapi(env, nodes, 0, IloInfinity, ILOFLOAT);
	// Mudei aqui!!!! Tirei a negativa!!!
	     
	// Inserção dos Big-M
	
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
       // master.setParam(IloCplex::Probe, 3);
        //master.setParam(IloCplex::MIPEmphasis, 1);
        // master.setParam(IloCplex::DisjCuts, 2);
        // master.setParam(IloCplex::MCFCuts, 2);
        // master.setParam(IloCplex::FlowCovers, 2);
        // master.setParam(IloCplex::DisjCuts, 2);
        // master.setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_BALANCED);
        //master.setParam(IloCplex::EpGap, 0.02);
        //master.setParam(IloCplex::NodeFileInd, 2);
        //master.setParam(IloCplex::PreLinear, 0);
        // master.setParam(IloCplex::RootAlg, CPX_ALG_BARRIER);
        // master.setParam(IloCplex::NodeAlg, CPX_ALG_BARRIER);
        master.setOut(env.getNullStream());
        IloNumArray piDuals(env);
        IloNumArray sigmaDuals(env);
        p = 0;

                int stabilized = 0;
		masterValues[0] = NUMTRIPS * bigM;

// ------------------------------------------------------------------------------------------------------------------------------------------------
    // solução inicial
    
    	vector< vector< vector<bool> > > assignSImatrix(depots);
 		
    	for (int k=0; k < depots; k++)
        	assignSImatrix[k] =  vector< vector<bool> >(matrixSize);

    	for (int k=0; k < depots; k++)
		for (int i=0; i<matrixSize; i++)
        		assignSImatrix[k][i] =  vector<bool>(matrixSize);

    	for (int k=0; k < depots; k++)
        	for (int i=0; i<matrixSize; i++)
            		for (int j=0; j<matrixSize; j++)
                		assignSImatrix[k][i][j] = false;
    
    	// Roda JV e prepara para solução inicial
        bool JV_SI = true;
        cout << "Com Solução Inicial: " << (JV_SI ? "Sim" : "Não") << endl;
    	cout << "==============================================" << endl;

        // Reducoes do Denis
        vector<vector<int> > predSI(depots);
	for (int k=0; k<depots; k++) predSI[k] = vector<int>(matrixSize);
	
        
        reducoesDenis(costMatrix, depots, nodes, matrixSize, maxVehiclesPerDepot, predSI, assignSImatrix, false, false, false, JV_SI);
	if(JV_SI) {
		/*
		IloArray<IloNumArray> costMatrixSI(env, matrixSize);                  // Matriz de custos entre garagens e os pontos "i" e "j"
        	for(int i = 0; i < matrixSize; i++)
        		costMatrixSI[i] = IloNumArray(env, matrixSize);
        
        	for(int i = 0; i < matrixSize; i++)
        		for(int j = 0; j < matrixSize; j++) 
        			costMatrixSI[i][j] = costMatrix[i][j];
        	*/
        			
          // constroi os caminhos - solução inicial
	  for(int k = 0; k < depots; k++) {
            for(int i = 0; i < matrixSize; i++) {
                if(assignSImatrix[k][i][k]) {
		   // cout << "depot Ini" << k << " - ";
                    a[p] = IloNumArray(env, nodes);
                    int cost = 0;
                    // tirei um omega depois de costMatrixSI
                    if(constroiPath(i, k, p, predSI, costMatrix, depots, a, cost)) {
                        y[p] = IloNumVar(pathsSelect(cost)+pi(a[p])+sigma[k](1));
                        y[p].setBounds(0, 1);
                        omegasSelected[p] = false;
                        p++;
                    }
                }
            }
        }
 
}
// ------------------------------------------------------------------------------------------------------------------------------------------------
    		//master.exportModel("model.lp");
		while (continueGc)
		{
			pAnt = p;			
			testFeasibility6 = false;
			continueGc = false;
			
			master.solve(); 
			//master.exportModel("model.lp");
			// cout << iteracoes << " - " << master.getObjValue() << endl;
			if(masterValues[iteracoes-1] ==  master.getObjValue()) {
				stabilized++; 
			} else {
				stabilized = 0;
			}
			masterValues[iteracoes] = master.getObjValue();
			iteracoes++;
			if(stabilized > NUM_STABILIZED) continueGc = false;
			else {
				if(!JV_SI) master.getValues(deltapi, masterData);
				// Envio ao CPLEX e obtenção da solução da formulação
				// Se foi encontrada a solução, verifica se foi rodado um mínimo 
				// de iterações antes entre o início ou o último arredondamento
				// para habilitar novamente o processo de arredondamento
				
				if(JV_SI) {
					if ((iteracoes-nO>ITER_MIN) && (masterValues[iteracoes-2] - masterValues[iteracoes-1] <= ZMIN))
						testFeasibility6=true;
				} else if (feasibility2Test(nodes, masterData))
				{
					if ((iteracoes-nO>ITER_MIN) && (masterValues[iteracoes-2] - masterValues[iteracoes-1] <= ZMIN))
						testFeasibility6=true;
				}
				
				// Se não for arredondar, então, passa pelo processo de geração de novas colunas
				if (!testFeasibility6)
				{
					// Cria as viagens entre os nodos (sempre serão iguais para todas as garagens) - guarda valor de arestas
					master.getDuals(piDuals, pi);
					master.getDuals(sigmaDuals, sigma);
					for (int k=0; k<depots; k++)
					{
						// Roda o Bellman-ford
						if (SLF(pred, costMatrix, nodes, depots, k, piDuals, sigmaDuals)<0)
						{
							a[p] = IloNumArray(env, nodes);
							// - Ajusta os valores para verificar se o valor do shortest path 
							//   pode adicionar um caminho de redução ao problema mestre
							// - Se o custo total for menor que zero (retorno do shortest path),
							//   significa que o valor pode ser utilizado para reduzir a função objetivo do RMP
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
				// cout << "Stabilized" << endl;
				testFeasibility6 = true;
			}
			if (testFeasibility6)
			{
				if(!JV_SI) {
					if (!feasibility6Test(nodes, masterData))
					{
						cout << "No solution found, sorry :/" << endl;
						env.end();
						return 1; 
					}
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
						// Verifica se todos os resultados do problema são inteiros
						//
						//6 cout << "Integrability" << endl;
						for (p_=0; (p_<pAnt); p_++)
						{
							if ((master.getValue(y[p_]) != 0) && (master.getValue(y[p_]) != 1))
							{
								// cout << "No best solution found, yet..." << endl;
								// Se não forem, indica que não houve solução
								solution = false;
								break;
							}
						}


						// Se a integrabilidade está correta, mostra as informações e encerra o processo
						if (solution)
						{
							// Conta o número de veículos
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
							// Mostra a solução
							// cout << "Feasible solution found" << endl;
							cout << "Number of vehicles: " << nVehicles << endl;
							cout << "Objective value: " << master.getObjValue() << endl;
							env.end();
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
					// Caso não tenha ocorrido a integrabilidade
					// então precisamos de mais rodadas, com perturbação
					//
					// cout << "Rounding" << endl;
					checkOmega = false;
					lastOmega = 0;
					for (p_=0; p_<pAnt; p_++)
					//for (p_=p__; p_<pAnt; p_++)
					{
						//
						// Se o "p" já foi selecionado anteriormente, pula
						//
						if (omegasSelected[p_])
							continue;

						//
						// Se for maior que o "omegaMin", então marca como um caminho 
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
					// Se não achou, pega o maior da última rodada
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
			// Guarda o valor da última iteração
		}
		// Checa limite de iterações
		if (p >= ITERLIM)
		{
			cout << "Limite de iterações excedido, resultado obtido: " << masterValues[iteracoes-1] << endl;
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
