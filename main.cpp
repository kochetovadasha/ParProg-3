#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;

int* BFSM(int**mat, int size, int n , int point);
int BalanceVertexM(int** mat, int size, int* v1, int v1size, int* v2, int v2size);
int BalanceEdgeM(int** mat, int size, int* v1, int v1size, int* v2, int v2size);
int BalanceGraphM(int** mat, int size, int* v1, int v1size, int* v2, int v2size);
int ReturnPartitionBalance(int**mat, int*v, int size, int n);
int** ReturnSubMatrix(int**mat, int size, int*v, int n);

class Vertex
{
public:
	friend class Graph;
	Vertex() {}
	~Vertex() {}
};

class Graph
{
	Vertex** vert;
	int nofv;
public:
	int** adj;
	Graph()
	{
		vert = NULL;
		adj = NULL;
	}
	Graph(int n)
	{
		nofv = n;
		vert = new Vertex*[nofv];
		adj = new int*[nofv];
		for (unsigned i = 0; i < nofv; i++)
		{
			vert[i] = new Vertex();
			adj[i] = new int[nofv];
		}
	}
	Graph(int**mat, int n)
	{
		nofv = n;
		vert = new Vertex*[nofv];
		adj = new int*[nofv];
		for (unsigned i = 0; i < nofv; i++)
		{
			vert[i] = new Vertex();
		}
		adj = mat;
	}
	Graph(Graph& g1)
	{
		nofv = g1.nofv;
		vert = new Vertex*[nofv];
		adj = new int*[nofv];
		for (unsigned i = 0; i < nofv; i++)
		{
			vert[i] = g1.vert[i];
			if (i < (nofv+1)/2)
				for (unsigned j = i; j < nofv; j++)
					adj[j][i] = adj[i][j] = g1.adj[i][j];
		}
	}
	~Graph()
	{
		for (unsigned i = 0; i < nofv; i++)
		{
			if (vert[i])
				delete[]vert[i];
			if (adj[i])
				delete[]adj[i];
		}
		delete []vert;
		delete []adj;
	}

	void PrintMatrix()
	{
		for (unsigned i = 0; i < nofv; i++)
		{
			for (unsigned j = 0; j < nofv; j++)
				cout << adj[i][j] << ' ';
			cout << endl;
		}
	}

	void RandomGraph()
	{
		for (unsigned i = 0; i < nofv; i++)
		{
			for (unsigned j = i; j < nofv; j++)
			{
				if (i != j)
					adj[j][i] = adj[i][j] = (rand() % 100 < 45 ? (rand() % 9 + 1) : 0);
				else
					adj[i][j] = rand() % 9 + 1;
			}
		}
	}

	int* BFS(int n, int point)
	{
		if (point < 0 || point >= nofv || n <= 0 || n > nofv)
			return 0;
		int i;
		bool* visited = new bool[nofv];
		for (i = 0; i < nofv; i++)
			visited[i] = false;
		int *queue = new int[n];
		int count = 0, head = 0;
		queue[count++] = point;
		visited[point] = true;
		while (head < count)
		{
			if (count == n)
				break;
			point = queue[head++];
			for (i = 0; i < nofv; i++)
				if (adj[point][i] && !visited[i])
				{
					queue[count++] = i;
					visited[i] = true;
					if (count == n)
						break;
				}
		}
		return queue;
	}

	int BalanceVertex(Graph& g1)
	{
		int s1 = 0, s2 = 0;
		for (unsigned i = 0; i < nofv; i++)
		{
			s1 += adj[i][i];
		}
		for (unsigned i = 0; i < g1.nofv; i++)
		{
			s2 += g1.adj[i][i];
		}
		return abs(s1 - s2);
	}

	int BalanceEdge(Graph& g1)
	{
		int s1 = 0, s2 = 0;
		for (unsigned i = 0; i < nofv; i++)
		{
			for (unsigned j = i; j < nofv; j++)
				s1 += adj[i][j];
		}
		for (unsigned i = 0; i < g1.nofv; i++)
		{
			for (unsigned j = i; j < g1.nofv; j++)
				s2 += g1.adj[i][j];
		}
		return abs(s1 + s2);
	}

	int BalanceGraph(Graph& g1)
	{
		return BalanceVertex(g1) + BalanceEdge(g1);
	}

	Graph* ReturnGraph(int* v, int n)
	{
		Graph* res = new Graph(n);
		for (unsigned i = 0; i < n; i++)
		{
			if (v[i] > -1 && v[i] < nofv)
			{
				res->vert[i] = vert[v[i]];
				for (unsigned j = 0; j < n; j++)
				{
					if (v[j] > -1 && v[j] < nofv)
					res->adj[i][j] = adj[v[i]][v[j]];
				}
			}
			else
				return 0;
		}
		return res;
	}

	Graph* CreateBFSGraph(int n, int point)
	{
		int* v = new int[n];
		v = BFS(n, point);
		return ReturnGraph(v, n);
	}

	Graph** ReturnPartition(int*v, int n)
	{
		if (n <= 0 || n > nofv)
			return 0;
		int* arr = new int[nofv - n];
		int k = 0;
		for (unsigned i = 0; i < nofv; i++)
		{
			if (k == nofv - n)
				break;
			for (unsigned j = 0; j < n; j++)
			{
				if (v[j] >= nofv || v[j] < 0)
					return 0;
				if (v[j] == i)
				{
					break;
				}
				else
				if (j == n - 1)
				{
					arr[k] = i;
					k++;
				}
			}
		}
		
		Graph** gres = new Graph*[2];
		gres[0] = new Graph(n);
		gres[1] = new Graph(nofv - n);
		gres[0] = ReturnGraph(v, n);
		gres[1] = ReturnGraph(arr, nofv - n);
		return gres;
	}

	Graph** CreateBFSPartition(int n, int point)
	{
		int* v = new int[n];
		v = BFS(n, point);
		return ReturnPartition(v, n);
	}
};

int main(int argc, char* argv[])
{
	int procnum, procrank;
	double startwtime = 0.0;
	double endwtime;
	MPI_Status status;
	srand(time(NULL));	
	
	//int numofv = atoi(argv[1]);  //Число вершин в графе
	int numofv = 1000;  //Число вершин в графе
	int mid = numofv / 2;  //Число вершин в одной из частей разделения
	int randpoint;  // Случайная вершина в графе
	int** adjm = new int*[numofv];  //Матрица смежности
	int numofparts = 2;  // Количество частей на которые нужно разделить

	int nofusedprocesses;  //Число задействованных процессов
	int remainder = 0;     //Остаток операций, идущий 0 процессу
	Graph* g = NULL;

	unsigned i;  //Счётчики в циклах

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procrank);

	int searchcount = numofv / procnum;   //Количество поисков, проводимх каждым процессом
	if (searchcount > 1000)                //Ограничиваем количество поисков для огромных графов
		searchcount = 1000;

	if (procrank == 0)
	{
		startwtime = MPI_Wtime();
		g = new Graph(numofv);
		g->RandomGraph();
		//g->PrintMatrix();
		adjm = g->adj;

		for (i = 0; i < numofv; i++)
			MPI_Bcast(g->adj[i], numofv, MPI_INT, 0, MPI_COMM_WORLD);   // Рассылка остальным процессам матрицы

		if (searchcount == 0)
		{
			searchcount = 1;
			nofusedprocesses = numofv;
		}
		else
		{
			nofusedprocesses = procnum;
			remainder = numofv - procnum * searchcount;
			searchcount = searchcount + remainder;
		}
	}
	else
	{
		if (searchcount != 0 || procrank < numofv)  //Если searchcount==0, то задействованы только первые numofv процессов
		{
			if (searchcount == 0)  //Если попали в тело с условием, что searchcount==0, то процесс должен совершить 1 поиск
				searchcount = 1;
			
			for (i = 0; i < numofv; i++)  // Получение матрицы смежности исходного графа, созданного 0 процессом
			{
				adjm[i] = new int[numofv];
				MPI_Bcast(adjm[i], numofv, MPI_INT, 0, MPI_COMM_WORLD);
			}
		}
	}

	int bestbalance = -1, bestbalanceindex = -1;
	int* bestbalancearr = new int[mid];  //Массив вершин из получаемых BFS подграфов
	if (procrank == 0 || searchcount != 0 || procrank < numofv)
	{
		
		int* arr = new int[mid];  //Массив для последовательности вершин в BFS подграфе
		for (i = 0; i < searchcount; i++)
		{
			randpoint = rand() % numofv;
			arr = BFSM(adjm, numofv, mid, randpoint);
			int b = ReturnPartitionBalance(adjm, arr, numofv, mid);		// Создание разделения на основе BFS дерева
			if (bestbalance == -1 || bestbalance > b) // Определение, а не лучший ли это баланс
			{
				bestbalanceindex = i;
				bestbalance = b;
				bestbalancearr = arr;
			}
		}
	}
	
	if (procrank == 0)
	{
		int *bmas = new int[nofusedprocesses];
		bmas[0] = bestbalance;
		for (i = 1; i < nofusedprocesses; i++)
		{
			MPI_Recv(&bmas[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);  // Получение от процессов лучшего баланса
		}
		int min = bmas[0], mini = 0;
		for (i = 1; i < nofusedprocesses; i++)  // Поиск лучшего баланса из полученных
		{
			if (bmas[i] < min)
			{
				min = bmas[i];
				mini = i;
			}
		}
		int* arr = new int[mid];
		if (mini != 0)
		{
			MPI_Recv(arr, mid, MPI_INT, mini, 0, MPI_COMM_WORLD, &status);  // Получение результата с лучшим балансом
		}
		else
			arr = bestbalancearr;
		Graph** g2 = new Graph*[2];
		g2[0] = new Graph(mid);
		g2[1] = new Graph(numofv - mid);
		g2 = g->ReturnPartition(arr, mid);  // Создание двух графов разделением исходного
		
		/*for (i = 0; i < 2; i++)
		{
			cout << endl;
			g2[i]->PrintMatrix();
		}*/
		//cout << "BESTBALANCE = "<<bestbalance << endl;

		endwtime = MPI_Wtime();
		printf("Time = %f \n", (endwtime - startwtime) * 1000);
	}
	else
	if (searchcount != 0 || procrank < numofv)
	{
		MPI_Send(&bestbalance, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);  // Отправка лучшего баланса для определения, стоит ли рассматривать это разделение
		MPI_Send(bestbalancearr, mid, MPI_INT, 0, 0, MPI_COMM_WORLD); // Отправка разделения с лучшим балансом (будет получено 0 только, если тот согласится, что баланс лучший)
	}

	MPI_Finalize();
	return 0;
}

int* BFSM(int** mat, int size, int n, int point)
{
	if (point < 0 || point >= size || n <= 0 || n > size)
		return 0;
	int i;
	bool* visited = new bool[size];
	for (i = 0; i < size; i++)
		visited[i] = false;
	int *queue = new int[n];
	int count = 0, head = 0;
	queue[count++] = point;
	visited[point] = true;
	while (head < count)
	{
		if (count == n)
			break;
		point = queue[head++];
		for (i = 0; i < size; i++)
		if (mat[point][i] && !visited[i])
		{
			queue[count++] = i;
			visited[i] = true;
			if (count == n)
				break;
		}
	}
	return queue;
}

int BalanceVertexM(int** mat, int size, int* v1, int v1size, int* v2, int v2size)
{
	unsigned i;
	if (v1size > size || v2size > size)
		return 0;
	int s1 = 0, s2 = 0;
	for (i = 0; i < v1size; i++)
	{
		s1 += mat[v1[i]][v1[i]];
	}
	for (i = 0; i < v2size; i++)
	{
		s2 += mat[v2[i]][v2[i]];
	}
	return abs(s1 - s2);
}

int BalanceEdgeM(int** mat, int size, int* v1, int v1size, int* v2, int v2size)
{
	unsigned i, j;
	if (v1size > size || v2size > size)
		return 0;
	int s1 = 0, s2 = 0;
	for (i = 0; i < v1size; i++)
	{
		for (j = i + 1; j < v1size; j++)
			s1 += mat[v1[i]][v1[j]];
	}
	for (i = 0; i < v2size; i++)
	{
		for (j = i + 1; j < v2size; j++)
			s2 += mat[v2[i]][v2[j]];
	}
	return abs(s1 + s2);
}

int BalanceGraphM(int** mat, int size, int* v1, int v1size, int* v2, int v2size)
{
	return BalanceVertexM(mat, size, v1, v1size, v2, v2size) + BalanceEdgeM(mat, size, v1, v1size, v2, v2size);
}

int ReturnPartitionBalance(int**mat, int*v, int size, int n)
{
	unsigned i, j;
	if (n <= 0 || n > size)
		return 0;
	int* arr = new int[size - n];
	int k = 0;
	for (i = 0; i < size; i++)
	{
		if (k == size - n)
			break;
		for (j = 0; j < n; j++)
		{
			if (v[j] >= size || v[j] < 0)
				return 0;
			if (v[j] == i)
			{
				break;
			}
			else
			if (j == n - 1)
			{
				arr[k] = i;
				k++;
			}
		}
	}
	return BalanceGraphM(mat, size, v, n, arr, size - n);
}

int** ReturnSubMatrix(int**mat, int size, int* v, int n)
{
	unsigned i, j;
	int** res = new int*[n];
	for (i = 0; i < n; i++)
	{
		res[i] = new int[n];
		if (v[i] > -1 && v[i] < size)
			return 0;			
		for (j = 0; j < n; j++)
		{
			if (v[j] > -1 && v[j] < size)
				return 0;
			res[i][j] = mat[v[i]][v[j]];				
		}
	}
	return res;
}