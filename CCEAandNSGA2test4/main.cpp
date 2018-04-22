#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include <algorithm>
#define etax 20
#define etam 20
#define EPS 1e-8
using namespace std;
const double LOWBOUND = 0;
const double UPPBOUND = 1;
const int F1_NVAR = 1;
const int F2_NVAR = 29;
const int IN_POP_SIZE = 100;
const double P_CROSS = 0.5;
const double P_MUTAT = 0.15; //0.02
const int IN_GENE_TIMES = 10;
const int OUT_GENE_TIMES = 10;
struct F1_INDIVIDUAL
{
	double m_var[F1_NVAR];
	double m_obj;
	int m_fit;
	void calObj();
};
struct F1_POPULATION
{
	vector<F1_INDIVIDUAL> m_population;
	void initPopulation();
	void geneticOperation();
	void select(int &p1, int &p2);
	void crossover(const F1_INDIVIDUAL &p1, const F1_INDIVIDUAL &p2, F1_INDIVIDUAL &c);
	void mutation(F1_INDIVIDUAL &ind, double rate);
};
struct F2_INDIVIDUAL
{
	double m_var[F2_NVAR];
	double m_obj;
	int m_fit;
	F1_INDIVIDUAL f1;
	void calObj();
};
struct F2_POPULATION
{
	vector<F2_INDIVIDUAL> m_population;

	void initPopulation();
	void geneticOperation();
	void select(int &p1, int &p2);
	void crossover(const F2_INDIVIDUAL &p1, const F2_INDIVIDUAL &p2, F2_INDIVIDUAL &c);
	void mutation(F2_INDIVIDUAL &ind, double rate);
};
struct F_INDIVIDUAL
{
	//double m_f1_var[F1_NVAR];
	//double m_f2_var[F2_NVAR];
	double m_f1_obj;
	double m_f2_obj;
	bool operator>(const F_INDIVIDUAL &ind) const;
	bool operator<(const F_INDIVIDUAL &ind) const;
};
double createZeroToOneRand();
void rankSort(int *rank, int *id, int size);  //升序
void distanceSort(double *dis, int *id, int size);  //降序
bool fIndividualEqual(const F_INDIVIDUAL &ind1, const F_INDIVIDUAL &ind2); //ind1 == ind2 ？
bool fIndividualCmp(const F_INDIVIDUAL &ind1, const F_INDIVIDUAL &ind2); //ind1支配ind2 ？
void nonDomiSort(const vector<F_INDIVIDUAL> &f, int *rank, int size); //rank越小表示越优
//从f中选择较优的个体集合给r，size_of_f是f的大小，size_of_r是r的大小
void selectOptInds(const vector<F_INDIVIDUAL> &f, int size_of_f, vector<F_INDIVIDUAL> &r, int size_of_r);
void run();
F1_POPULATION P1;
F2_POPULATION P2;
vector<F_INDIVIDUAL> pareto_space;
vector<F_INDIVIDUAL> res_space;

int main()
{
	srand(time(NULL));
	run();
	ofstream out("outdata.txt");
	for (int i = 0; i < res_space.size(); ++i)
	{
		out << res_space[i].m_f1_obj << "\t" << res_space[i].m_f2_obj << endl;
	}
	return 0;
}

void run()
{
	P1.initPopulation();
	P2.initPopulation();
	for (int out_gene_times = 0; out_gene_times < OUT_GENE_TIMES; out_gene_times++)
	{
		vector<F_INDIVIDUAL> res_space_tmp1;
		vector<F_INDIVIDUAL> res_space_tmp2;
		for (int in_gene_times = 0; in_gene_times < IN_GENE_TIMES; in_gene_times++)
		{
			pareto_space.clear();
			for (int i = 0; i < P1.m_population.size(); i++)
			{
				F_INDIVIDUAL f;
				P1.m_population[i].calObj();
				f.m_f1_obj = P1.m_population[i].m_obj;
				for (int j = 0; j < 6; j++)
				{
					P2.m_population[j].f1 = P1.m_population[i];
					P2.m_population[j].calObj();
					f.m_f2_obj = P2.m_population[j].m_obj;
					pareto_space.push_back(f);
				}
			}
			int rank[6 * IN_POP_SIZE];
			nonDomiSort(pareto_space, rank, 6 * IN_POP_SIZE);  //pareto_space.size()==6*IN_POP_SIZE
			for (int i = 0; i < IN_POP_SIZE; i++)
			{
				int fit = 0;
				for (int j = 0; j < 6; j++)
				{
					fit += rank[6 * i + j];
				}
				P1.m_population[i].m_fit = fit;  //fit值越小越优
			}
			P1.geneticOperation();
		}
		selectOptInds(pareto_space, pareto_space.size(), res_space_tmp1, IN_POP_SIZE);
		for (int i = 0; i < res_space_tmp1.size(); i++) res_space_tmp2.push_back(res_space_tmp1[i]);
		pareto_space.clear();

		for (int in_gene_times = 0; in_gene_times < IN_GENE_TIMES; in_gene_times++)
		{
			pareto_space.clear();
			for (int i = 0; i < P2.m_population.size(); i++)
			{
				F_INDIVIDUAL f;
				for (int j = 0; j < 6; j++)
				{
					P1.m_population[j].calObj();
					P2.m_population[i].f1 = P1.m_population[j];
					P2.m_population[i].calObj();
					f.m_f1_obj = P1.m_population[j].m_obj;
					f.m_f2_obj = P2.m_population[i].m_obj;
					pareto_space.push_back(f);
				}
			}
			int rank[6 * IN_POP_SIZE];
			nonDomiSort(pareto_space, rank, 6 * IN_POP_SIZE);  //pareto_space.size()==6*IN_POP_SIZE
			for (int i = 0; i < IN_POP_SIZE; i++)
			{
				int fit = 0;
				for (int j = 0; j < 6; j++)
				{
					fit += rank[6 * i + j];
				}
				P2.m_population[i].m_fit = fit;  //fit值越小越优
			}
			P2.geneticOperation();
		}
		selectOptInds(pareto_space, pareto_space.size(), res_space_tmp1, IN_POP_SIZE);
		for (int i = 0; i < res_space_tmp1.size(); i++) res_space_tmp2.push_back(res_space_tmp1[i]);
		selectOptInds(res_space_tmp2, res_space_tmp2.size(), res_space, IN_POP_SIZE);
		pareto_space.clear(); 
		res_space_tmp1.clear();
		res_space_tmp2.clear();
	}
}

double createZeroToOneRand()
{
	return rand() * 1.0 / RAND_MAX;
}

void F1_POPULATION::initPopulation()
{
	F1_INDIVIDUAL ind;
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		for (int j = 0; j < F1_NVAR; ++j)
		{
			ind.m_var[j] = LOWBOUND + createZeroToOneRand()*(UPPBOUND - LOWBOUND);
		}
		m_population.push_back(ind);
	}
}
void F2_POPULATION::initPopulation()
{
	F2_INDIVIDUAL ind;
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		for (int j = 0; j < F2_NVAR; ++j)
		{
			ind.m_var[j] = LOWBOUND + createZeroToOneRand()*(UPPBOUND - LOWBOUND);
		}
		m_population.push_back(ind);
	}
}
void F1_INDIVIDUAL::calObj()
{
	for (int i = 0; i < F1_NVAR; ++i)
	{
		m_obj = m_var[i];
	}
}
void F2_INDIVIDUAL::calObj()
{
	double g = 0;
	for (int i = 0; i < F2_NVAR; ++i)
	{
		g = g + m_var[i];
	}
	g = g / F2_NVAR;
	g = 1 + 9 * g;
	f1.calObj();
	double h = 1 - sqrt(f1.m_obj / g);
	m_obj = g * h;
}
void F1_POPULATION::geneticOperation()
{
	vector<F1_INDIVIDUAL> p_tmp;
	int is_chosed[6 * IN_POP_SIZE];
	for (int i = 0; i < 6 * IN_POP_SIZE; i++) is_chosed[i] = 0;
	for (int i = 0; i < 6; i++)   //从m_population中选出最优的6个个体，加入到集合p_tmp中
	{
		int minfit = INT16_MAX;
		int id;
		for (int j = 0; j < m_population.size(); j++)  //m_population.size() == IN_POP_SIZE
		{
			if (is_chosed[j] == 0 && minfit > m_population[j].m_fit)
			{
				minfit = m_population[j].m_fit;
				id = j;
			}
		}
		is_chosed[id] = 1;
		p_tmp.push_back(m_population[id]);
	}
	for (int i = 6; i < IN_POP_SIZE; i++)
	{
		int index1, index2;
		select(index1, index2);
		F1_INDIVIDUAL child;
		crossover(m_population[index1], m_population[index2], child);
		mutation(child, P_MUTAT);
		p_tmp.push_back(child);
	}
	m_population.clear();
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		m_population.push_back(p_tmp[i]);
	}
}
void F1_POPULATION::select(int &p1, int &p2)   //规模为2的锦标赛选择算法
{
	int t1 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	int t2 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	if (m_population[t1].m_fit < m_population[t2].m_fit)
		p1 = t1;
	else
		p1 = t2;
	t1 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	t2 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	if (m_population[t1].m_fit < m_population[t2].m_fit)
		p2 = t1;
	else
		p2 = t2;
}
void F1_POPULATION::crossover(const F1_INDIVIDUAL &p1, const F1_INDIVIDUAL &p2, F1_INDIVIDUAL &c)
{
	double rand;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double eta_c = etax;
	if (createZeroToOneRand() <= 1.0)
	{
		for (int i = 0; i<F1_NVAR; i++)
		{
			if (createZeroToOneRand() <= P_CROSS)
			{
				if (fabs(p1.m_var[i] - p2.m_var[i]) > EPS)
				{
					if (p1.m_var[i] < p2.m_var[i])
					{
						y1 = p1.m_var[i];
						y2 = p2.m_var[i];
					}
					else
					{
						y1 = p2.m_var[i];
						y2 = p1.m_var[i];
					}
					yl = LOWBOUND;
					yu = UPPBOUND;
					rand = createZeroToOneRand();
					beta = 1.0 + (2.0*(y1 - yl) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
					beta = 1.0 + (2.0*(yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));
					if (c1<yl)
						c1 = yl;
					if (c2<yl)
						c2 = yl;
					if (c1>yu)
						c1 = yu;
					if (c2>yu)
						c2 = yu;
					if (createZeroToOneRand() <= 0.5)
					{
						c.m_var[i] = c2;
					}
					else
					{
						c.m_var[i] = c1;
					}
				}
				else
				{
					c.m_var[i] = p1.m_var[i];
				}
			}
			else
			{
				c.m_var[i] = p1.m_var[i];
			}
		}
	}
	else
	{
		for (int i = 0; i<F1_NVAR; i++)
		{
			c.m_var[i] = p1.m_var[i];
		}
	}
	return;
}
void F1_POPULATION::mutation(F1_INDIVIDUAL &ind, double rate)
{
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	double eta_m = etam;
	int id_rnd = int(createZeroToOneRand()*F1_NVAR);
	for (int j = 0; j < F1_NVAR; j++)
	{
		if (createZeroToOneRand() <= rate)
		{
			y = ind.m_var[j];
			yl = LOWBOUND;
			yu = UPPBOUND;
			delta1 = (y - yl) / (yu - yl);
			delta2 = (yu - y) / (yu - yl);
			rnd = createZeroToOneRand();
			mut_pow = 1.0 / (eta_m + 1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0 - delta1;
				val = 2.0*rnd + (1.0 - 2.0*rnd)*(pow(xy, (eta_m + 1.0)));
				deltaq = pow(val, mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0 - delta2;
				val = 2.0*(1.0 - rnd) + 2.0*(rnd - 0.5)*(pow(xy, (eta_m + 1.0)));
				deltaq = 1.0 - (pow(val, mut_pow));
			}
			y = y + deltaq*(yu - yl);
			if (y<yl)
				y = yl;
			if (y>yu)
				y = yu;
			ind.m_var[j] = y;
		}
	}
	return;
}
void F2_POPULATION::geneticOperation()
{
	vector<F2_INDIVIDUAL> p_tmp;
	int is_chosed[6 * IN_POP_SIZE];
	for (int i = 0; i < 6 * IN_POP_SIZE; i++) is_chosed[i] = 0;
	for (int i = 0; i < 6; i++)   //从m_population中选出最优的6个个体，加入到集合p_tmp中
	{
		int minfit = INT16_MAX;
		int id;
		for (int j = 0; j < m_population.size(); j++)  //m_population.size() == IN_POP_SIZE
		{
			if (is_chosed[j] == 0 && minfit > m_population[j].m_fit)
			{
				minfit = m_population[j].m_fit;
				id = j;
			}
		}
		is_chosed[id] = 1;
		p_tmp.push_back(m_population[id]);
	}
	for (int i = 6; i < IN_POP_SIZE; i++)
	{
		int index1, index2;
		select(index1, index2);
		F2_INDIVIDUAL child;
		crossover(m_population[index1], m_population[index2], child);
		mutation(child, P_MUTAT);
		p_tmp.push_back(child);
	}
	m_population.clear();
	for (int i = 0; i < IN_POP_SIZE; ++i)
	{
		m_population.push_back(p_tmp[i]);
	}
}
void F2_POPULATION::select(int &p1, int &p2)
{
	int t1 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	int t2 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	if (m_population[t1].m_fit > m_population[t2].m_fit)
		p1 = t1;
	else
		p1 = t2;
	t1 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	t2 = int(createZeroToOneRand() * IN_POP_SIZE) % IN_POP_SIZE;
	if (m_population[t1].m_fit > m_population[t2].m_fit)
		p2 = t1;
	else
		p2 = t2;
}
void F2_POPULATION::crossover(const F2_INDIVIDUAL &p1, const F2_INDIVIDUAL &p2, F2_INDIVIDUAL &c)
{
	double rand;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double eta_c = etax;
	if (createZeroToOneRand() <= 1.0)
	{
		for (int i = 0; i<F2_NVAR; i++)
		{
			if (createZeroToOneRand() <= P_CROSS)
			{
				if (fabs(p1.m_var[i] - p2.m_var[i]) > EPS)
				{
					if (p1.m_var[i] < p2.m_var[i])
					{
						y1 = p1.m_var[i];
						y2 = p2.m_var[i];
					}
					else
					{
						y1 = p2.m_var[i];
						y2 = p1.m_var[i];
					}
					yl = LOWBOUND;
					yu = UPPBOUND;
					rand = createZeroToOneRand();
					beta = 1.0 + (2.0*(y1 - yl) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
					beta = 1.0 + (2.0*(yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha))
					{
						betaq = pow((rand*alpha), (1.0 / (eta_c + 1.0)));
					}
					else
					{
						betaq = pow((1.0 / (2.0 - rand*alpha)), (1.0 / (eta_c + 1.0)));
					}
					c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));
					if (c1<yl)
						c1 = yl;
					if (c2<yl)
						c2 = yl;
					if (c1>yu)
						c1 = yu;
					if (c2>yu)
						c2 = yu;
					if (createZeroToOneRand() <= 0.5)
					{
						c.m_var[i] = c2;
					}
					else
					{
						c.m_var[i] = c1;
					}
				}
				else
				{
					c.m_var[i] = p1.m_var[i];
				}
			}
			else
			{
				c.m_var[i] = p1.m_var[i];
			}
		}
	}
	else
	{
		for (int i = 0; i<F1_NVAR; i++)
		{
			c.m_var[i] = p1.m_var[i];
		}
	}
	return;
}
void F2_POPULATION::mutation(F2_INDIVIDUAL &ind, double rate)
{
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	double eta_m = etam;
	int id_rnd = int(createZeroToOneRand()*F1_NVAR) % F1_NVAR;
	for (int j = 0; j < F1_NVAR; j++)
	{
		if (createZeroToOneRand() <= rate)
		{
			y = ind.m_var[j];
			yl = LOWBOUND;
			yu = UPPBOUND;
			delta1 = (y - yl) / (yu - yl);
			delta2 = (yu - y) / (yu - yl);
			rnd = createZeroToOneRand();
			mut_pow = 1.0 / (eta_m + 1.0);
			if (rnd <= 0.5)
			{
				xy = 1.0 - delta1;
				val = 2.0*rnd + (1.0 - 2.0*rnd)*(pow(xy, (eta_m + 1.0)));
				deltaq = pow(val, mut_pow) - 1.0;
			}
			else
			{
				xy = 1.0 - delta2;
				val = 2.0*(1.0 - rnd) + 2.0*(rnd - 0.5)*(pow(xy, (eta_m + 1.0)));
				deltaq = 1.0 - (pow(val, mut_pow));
			}
			y = y + deltaq*(yu - yl);
			if (y<yl)
				y = yl;
			if (y>yu)
				y = yu;
			ind.m_var[j] = y;
		}
	}
	return;
}


void nonDomiSort(const vector<F_INDIVIDUAL> &f, int *rank, int size)
{ 
	int *np = new int[size];  //被支配的数目
	for (int i = 0; i < size; ++i) np[i] = 0;  //初始化
	vector<int> *sp = new vector<int>[size];  //支配的集合
	for (int i = 0; i < size; ++i)       //耗时太长
	{
		for (int j = i + 1; j < size; ++j)
		{
			if (f[i] > f[j])
			{
				np[j]++;
				sp[i].push_back(j);
			}
			else if (f[j] > f[i])
			{
				np[i]++;
				sp[j].push_back(i);
			}
		}
	}
	int curr_rank = 0;  //当前的rank值
	int *is_done = new int[size];  //是否已被赋过rank值，0为否，1为是
	for (int i = 0; i < size; ++i) is_done[i] = 0;
	while (1)
	{
		int *np_tmp = new int[size];
		for (int i = 0; i < size; ++i) np_tmp[i] = np[i];
		bool flag = false;
		for (int i = 0; i < size; ++i)
		{
			if (0 == np[i] && (0 == is_done[i]))   //剩下的里边没有组合能够支配它
			{
				flag = true;
				rank[i] = curr_rank;
				is_done[i] = 1;
				for (int j = 0; j < sp[i].size(); ++j)
				{
					np_tmp[sp[i][j]]--;
				}
			}
		}
		for (int i = 0; i < size; ++i) np[i] = np_tmp[i];
		delete[] np_tmp;
		curr_rank++;
		if (flag == false) break;
	}
	delete[] np;
	for (int i = 0; i < size; i++) sp[i].clear();
	delete[] sp;
	delete[] is_done;
}
bool fIndividualCmp(const F_INDIVIDUAL &ind1, const F_INDIVIDUAL &ind2)
{
	if (ind1.m_f1_obj > ind2.m_f1_obj || ind1.m_f2_obj > ind2.m_f2_obj)
	{
		return false;
	}
	if (fIndividualEqual(ind1, ind2) == true)
	{
		return false;
	}
	return true;
}
bool fIndividualEqual(const F_INDIVIDUAL &ind1, const F_INDIVIDUAL &ind2)
{
	if (ind1.m_f1_obj == ind2.m_f1_obj && ind1.m_f2_obj == ind2.m_f2_obj)
	{
		return false;
	}
	return false;
}
void selectOptInds(const vector<F_INDIVIDUAL> &f, int size_of_f, vector<F_INDIVIDUAL> &r, int size_of_r)
{
	int *rank = new int[size_of_f];
	nonDomiSort(f, rank, size_of_f);
	int *id = new int[size_of_f];
	for (int i = 0; i < size_of_f; i++)
	{
		id[i] = i;
	}
	rankSort(rank, id, size_of_f);
	r.clear();
	int curr_rank = 0;
	int b = 0, e = 0;
	while (1)
	{
		while (rank[e] == curr_rank) e++;
		if (r.size() + e - b > size_of_r) break;
		for (int i = b; i < e; i++) r.push_back(f[id[i]]);
		b = e;
		curr_rank++;
		if (r.size() >= size_of_r) break;
	}
	if (r.size() < size_of_r)
	{
		double *distance = new double[e - b];
		for (int i = 0; i < e - b; i++) distance[i] = 0;
		int *id2 = new int[e - b];
		for (int i = 0; i < e - b; i++) id2[i] = id[b + i];
		double *obj1 = new double[e - b];
		double *obj2 = new double[e - b];
		for (int i = 0; i < e - b; i++)
		{
			obj1[i] = f[id[b + i]].m_f1_obj;
		}
		distanceSort(obj1, id2, e - b);
		for (int i = 0; i < e - b; i++)
		{
			obj2[i] = f[id2[i]].m_f2_obj;
		}
		distance[0] = 1.0e+30;
		distance[e - b - 1] = 1.0e+30;
		for (int k = 1; k < e - b - 1; k++)
		{
			distance[k] = distance[k] + abs(obj1[k+1] - obj1[k - 1]) + abs(obj2[k + 1] - obj2[k - 1]);
		}
		distanceSort(distance, id2, e - b);
		int mm = size_of_r - r.size();
		for (int i = 0; i < mm; i++)
		{
			r.push_back(f[id2[i]]);
		}
		delete[] distance;
		delete[] id2;
		delete[] obj1;
		delete[] obj2;
	}
	delete[] rank;
	delete[] id;
}
void rankSort(int *rank, int *id, int size)   //升序
{
	for (int i = 0; i < size; ++i)
	{
		for (int j = i + 1; j < size; ++j)
		{
			if (rank[i] > rank[j])
			{
				int tmp = rank[i];
				rank[i] = rank[j];
				rank[j] = tmp;
				tmp = id[i];
				id[i] = id[j];
				id[j] = tmp;
			}
		}
	}
}

void distanceSort(double *dis, int *id, int size)  //降序
{
	for (int i = 0; i < size; ++i)
	{
		for (int j = i + 1; j < size; ++j)
		{
			if (dis[i] < dis[j])
			{
				double tmp1 = dis[i];
				dis[i] = dis[j];
				dis[j] = tmp1;
				int tmp2 = id[i];
				id[i] = id[j];
				id[j] = tmp2;
			}
		}
	}
}

bool F_INDIVIDUAL::operator>(const F_INDIVIDUAL &ind) const
{
	if (m_f1_obj > ind.m_f1_obj || m_f2_obj > ind.m_f2_obj)
	{
		return false;
	}
	if (m_f1_obj == ind.m_f1_obj && m_f2_obj == ind.m_f2_obj)
	{
		return false;
	}
	return true;
}
bool F_INDIVIDUAL::operator<(const F_INDIVIDUAL &ind) const
{
	if (m_f1_obj < ind.m_f1_obj || m_f2_obj < ind.m_f2_obj)
	{
		return false;
	}
	if (m_f1_obj == ind.m_f1_obj && m_f2_obj == ind.m_f2_obj)
	{
		return false;
	}
	return true;
}