#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <climits>

#define MAX 100001
#define MAXCTC 10001
#define dijMax 50001
#define MAXbf  1001
#define MAXapm 2001
#define MAXmtrx  101

using namespace std;

ifstream inbfs("bfs.in");
ofstream outbfs("bfs.out");

ifstream indfs("dfs.in");
ofstream outdfs("dfs.out");

ifstream inhakimi("hakimi.in");
ofstream outhakimi("hakimi.out");

ifstream insortaret("sortaret.in");
ofstream outsortaret("sortaret.out");

ifstream inctc("ctc.in");
ofstream outctc("ctc.out");

ifstream indijkstra("dijkstra.in");
ofstream outdijkstra("dijkstra.out");

ifstream inbellmanford("bellmanford.in");
ofstream outbellmanford("bellmanford.out");

ifstream inapm("apm.in");
ofstream outapm("apm.out");

ifstream indisjoint("disjoint.in");
ofstream outdisjoint("disjoint.out");

ifstream inroyfloyd("royfloyd.in");
ofstream outroyfloyd("royfloyd.out");

ifstream indarb("darb.in");
ofstream outdarb("darb.out");

ifstream inmaxflow("maxflow.in");
ofstream outmaxflow("maxflow.out");

ifstream ineuler("ciclueuler.in");
ofstream outeuler("ciclueuler.out");

class Graf
{
    int NrNoduri;
    vector<int> Adiacenta[MAX];

public:

    Graf(int NrNoduri);

    void AdaugaMuchie(int nod, int nodConectat);
    void AdaugaMUchieNeorientat(int nod, int nodConectat);

    //BFS
    vector<int> Bfs(int nod);

    //DFS
    void Dfs(int nod, bool Vizitat[MAX]);
    int DfsCompCnx (int nod, bool Vizitat[MAX], stack<int> &Stiva, int &total);

    //HAVEL HAKIMI
    int HavelHakimi(vector<int> Hakimi, int L);

    //SORTARE TOPOLOGICA

    void DfsTop(int nod, bool Vizitat[MAX], stack<int> &Stiva, int &total);
    void Sortare_Topologica(bool Vizitat[MAX], stack<int> &Stiva, int &total);

    //CTC
    void UmpleStiva(int nod, bool VizitatT[MAX], stack <int> &Stiva);
    void Transpus(vector <int> AdiacentaT[MAXCTC]);
    void Parcurgere(bool VizitatT[MAX], stack <int> &Stiva);
    void DfsT(int nod, vector <int> AdiacentaT[MAXCTC], vector <int> componente[MAXCTC], int total, bool Vizitat[MAX]);
    void Ctc(vector <int> AdiacentaT[MAX], vector <int> componente[MAXCTC], int &total, bool Vizitat[MAX], stack <int> &Stiva);

    //DIJKSTRA
    vector<int> Dijkstra(int nod, vector <pair<int, int>>G[dijMax]);

    //BELMANFORD
    vector<int> BellmanFord(int nod, vector<pair<int, int>>G[MAXbf]);

    //APM
    vector<pair<int, int>> Apm(int nod, vector <pair<int, int>> G[MAXapm], int &Tcost);

    //DISJOINT
    int gaseste(int nod, int v[MAX]);
    void uneste(int nod1, int nod2, int v[MAX], int inaltime[MAX]);

    //ROYFLOYD
    void RoyFloyd(int G[MAXmtrx][MAXmtrx]);

    //DARB
    void Bfsdarb(int nod, int *nodinitial, int *nodfinal, int distanta[]);

    //MAXFLOW
    bool bfs(int rGraf[MAXmtrx][MAXmtrx], int sursa, int destinatie, int parinte[]);
    int MaxFlow(int Graf[MAXmtrx][MAXmtrx], int sursa, int destinatie);

    //EULER
    vector<int> Euler(int nod, bool Vizitat[MAXCTC], vector<pair<int,int>> AdiacentaEuler[MAX]);
};

//ADAUGARE MUCHII
Graf::Graf(int NrNoduri)
{
    this->NrNoduri = NrNoduri;
}
void Graf::AdaugaMuchie(int nod, int nodConectat)
{
    Adiacenta[nod].push_back(nodConectat);
}

void Graf::AdaugaMUchieNeorientat(int nod, int nodConectat)
{
    Adiacenta[nod].push_back(nodConectat);
    Adiacenta[nodConectat].push_back(nod);
}

/*
Voi folosi Bfsdarb(implementeaza atat Bfs, cat si Darb.
Pt. Bfs ma intereseaza doar vectorul distanta. Il afisez
cu distanta[i] - 1.
///BFS
vector<int> Graf::Bfs(int nod)
{
    vector<int> distanta(NrNoduri);
    distanta[nod] = 1;

    queue<int> Q;
    bool Vizitat[NrNoduri] = {0};

    Q.push(nod);
    Vizitat[nod] = 1;

    while(!Q.empty())
    {
        nod = Q.front();
        Q.pop();

        for(auto i: Adiacenta[nod])
        {
            if(!Vizitat[i])
            {
                Q.push(i);
                Vizitat[i] = 1;
                distanta[i] = distanta[nod] + 1;
            }
        }
    }

    for(int i = 1; i <= NrNoduri; i++)
        distanta[i]--;

    return distanta;
}
*/
///DFS
/*
Voi folosi DfsTop(face Dfs, dar si sortare topologica. Nu ma intereseaza
                  Stiva si total cand vreau sa fac dfs)
void Graf::Dfs(int nod, bool Vizitat[MAX])
{
    Vizitat[nod] = 1;

    for (auto i: Adiacenta[nod])
        if (!Vizitat[i])
            Dfs(i, Vizitat);
}
*/
int Graf::DfsCompCnx (int nod, bool Vizitat[MAX], stack<int> &Stiva, int &total)
{

    int tot = 0;
    for(int i = 0; i < NrNoduri; i++)
        if(!Vizitat[i])
        {
            DfsTop(i, Vizitat, Stiva, total);
            tot++;
        }

    return tot;
}


///HAVEL HAKIMI
int Graf::HavelHakimi(vector<int> Hakimi, int L)
{
    while(1)
    {
        sort(Hakimi.begin(), Hakimi.end(), greater<>());

        int primul = Hakimi[0];

        if(primul == 0)
            return 1;

        if(primul > Hakimi.size())
            return 0;

        Hakimi.erase(Hakimi.begin());

        for(int i = 0; i < primul; i++)
        {
            Hakimi[i]--;

            if(Hakimi[i] < 0)
                return 0;
        }
    }
}

///SORTARE TOPOLOGICA

void Graf::DfsTop(int nod,bool Vizitat[MAX], stack<int> &Stiva, int &total)
{
    Vizitat[nod] = 1;

    for (auto i: Adiacenta[nod])
        if (!Vizitat[i])
            DfsTop(i, Vizitat, Stiva, total);

    Stiva.push(nod);
}

void Graf::Sortare_Topologica(bool Vizitat[MAX], stack<int> &Stiva, int &total)
{
    for(int i = 1; i <= NrNoduri; i++)
        if(!Vizitat[i])
            DfsTop(i, Vizitat, Stiva, total);
}

///CTC
void Graf :: Transpus(vector <int> AdiacentaT[MAXCTC])
{
    for(int i = 1; i <= NrNoduri; ++i)
        for(auto j : Adiacenta[i])
            AdiacentaT[j].push_back(i);
}

void Graf :: UmpleStiva(int nod, bool VizitatT[MAX], stack <int> &Stiva)
{
    VizitatT[nod] = 1;
    for(auto i : Adiacenta[nod])
        if(!VizitatT[i])
            UmpleStiva(i, VizitatT, Stiva);

    Stiva.push(nod);
}

void Graf :: Parcurgere(bool VizitatT[MAX], stack <int> &Stiva)
{
    for(int i = 1; i <= NrNoduri; ++i)
        if(!VizitatT[i])
            UmpleStiva(i, VizitatT, Stiva);
}

void Graf :: DfsT(int nod, vector <int> AdiacentaT[MAXCTC], vector <int> componente[MAXCTC], int total, bool Vizitat[MAX])
{
    Vizitat[nod] = 1;
    componente[total].push_back(nod);

    for(auto i: AdiacentaT[nod])
        if(!Vizitat[i])
            DfsT(i, AdiacentaT, componente, total, Vizitat);
}

void Graf :: Ctc(vector <int> AdiacentaT[MAXCTC], vector <int> componente[MAXCTC], int &total, bool Vizitat[MAX], stack <int> &Stiva)
{
    while(!Stiva.empty())
    {
        int nod = Stiva.top();
        Stiva.pop();

        if(!Vizitat[nod])
        {
            total++;
            DfsT(nod, AdiacentaT, componente, total, Vizitat);
        }
    }

    /*
    Voi folosi asta pentru a afisa in main
    outctc << total << "\n";

    for(int i = 1; i <= total; i++)
    {
        for(auto j: componente[i])
            outctc << j << " ";

        outctc<<"\n";
    }
    */
}


///DIJKSTRA
vector<int> Graf::Dijkstra(int nod, vector <pair<int, int>>G[dijMax])
{
    vector<int> distanta(dijMax);
    priority_queue <pair<int, int>, vector <pair<int,int>>, greater<pair<int,int>>> Coada;

    bool Vizitat[dijMax] = {0};

    for(int i = 1; i <= NrNoduri; i++)
        distanta[i] = INT_MAX;

    distanta[nod] = 0;
    Coada.push(make_pair(distanta[nod], nod));
    Vizitat[nod] = 1;

    while(!Coada.empty())
    {

        int nodCurent = Coada.top().second;
        Coada.pop();
        Vizitat[nodCurent] = 0;

        for(int i = 0; i < G[nodCurent].size(); i++)
        {
            int Vecin = G[nodCurent][i].first;
            int Cost = G[nodCurent][i].second;

            if(distanta[nodCurent] + Cost < distanta[Vecin])
            {
                distanta[Vecin] = distanta[nodCurent] + Cost;

                if(!Vizitat[Vecin])
                {
                    Coada.push(make_pair(distanta[Vecin], Vecin));
                    Vizitat[Vecin] = 1;
                }
            }
        }
    }

    return distanta;
    /*
    Voi folosi asta pentru a afisa in main
    for(int i = 2; i <= NrNoduri; i++)

        if(distanta[i] != INT_MAX)
            outdijkstra << distanta[i] << " ";

        else
            outdijkstra << 0 <<" ";
            */
}

///BELLMANFORD
vector<int> Graf :: BellmanFord(int nod, vector<pair<int, int>>G[MAXbf])
{

    queue <int> coada;
    int Vizitat[dijMax] = {0};
    vector<int> distanta(dijMax);
    bool InCoada[dijMax];

    for(int i = 1; i <= NrNoduri; i++)
        distanta[i] = INT_MAX;

    distanta[nod] = 0;
    coada.push(nod);
    InCoada[nod] = 1;

    while(!coada.empty())
    {
        int nodcurent = coada.front();
        Vizitat[nodcurent]++;

        if(Vizitat[nodcurent] >= NrNoduri)
        {
            distanta[0] = -100; //Verific in main daca asta e 0 si daca e, afisez "Ciclu negativ".
            return distanta;
        }

        coada.pop();
        InCoada[nodcurent] = 0;

        for(int i = 0; i < G[nodcurent].size(); i++)
        {

            int vecin = G[nodcurent][i].first;
            int cost = G[nodcurent][i].second;

            if(distanta[nodcurent] + cost < distanta[vecin])
            {
                distanta[vecin] = distanta[nodcurent] + cost;

                if(!InCoada[vecin])
                    coada.push(vecin);
            }
        }
    }

    /*
    Voi folosi asta pentru a afisa in main
    for(int i = 2; i <= NrNoduri; i++)
        outbellmanford << distanta[i] << " ";

    */

    return distanta;
}

///APM
vector<pair<int, int>> Graf::Apm(int nod, vector <pair<int, int>> G[MAXapm], int &Tcost)
{
    priority_queue <pair<int, int>, vector<pair <int, int>>, greater<pair<int,int>>> Coada;
    bool Vizitat[MAX] = {0};
    int distanta[MAX];
    vector <int> tata(MAX);
    vector <pair<int, int>> arbore;

    for(int i = 1; i <= NrNoduri; i++)
        distanta[i] = INT_MAX;



    distanta[nod] = 0;
    Coada.push(make_pair(0, nod));

    while(!Coada.empty())
    {
        int nodcurent = Coada.top().second;
        Coada.pop();

        if(!Vizitat[nodcurent])
        {
            Vizitat[nodcurent] = 1;
            Tcost += distanta[nodcurent];

            if(tata[nodcurent] > 0)
                arbore.push_back(make_pair(nodcurent, tata[nodcurent]));

            for(int i = 0; i < G[nodcurent].size(); i++)
            {
                int Vecin = G[nodcurent][i].first;
                int Cost = G[nodcurent][i].second;

                if(!Vizitat[Vecin] && distanta[Vecin] > Cost)
                {
                    distanta[Vecin] = Cost;
                    tata[Vecin] = nodcurent;
                    Coada.push(make_pair(distanta[Vecin], Vecin));
                }
            }
        }
    }

    /*
    Voi folosi asta pentru a afisa in main
    outapm << Tcost << "\n" << arbore.size() << "\n";

    for(int i = 0; i < arbore.size(); i++)
        outapm << arbore[i].first << " " << arbore[i].second << "\n";
    */

    return arbore;
}

///DISJOINT
int Graf::gaseste(int nod, int v[MAX])
{
    int i = nod;

    while(v[i] != i)
        i = v[i];

    while(v[nod] != nod)
    {
        int aux = v[nod];
        v[nod] = i;
        nod = aux;
    }

    return i;
}

void Graf::uneste(int nod1, int nod2, int v[MAX], int inaltime[MAX])
{
    if(inaltime[nod1] > inaltime[nod2])
        v[nod2] = nod1;

    else
        v[nod1] = nod2;

    if(inaltime[nod1] == inaltime[nod2])
        inaltime[nod2]++;

}

///ROYFLOYD
void Graf :: RoyFloyd(int G[MAXmtrx][MAXmtrx])
{
    for(int k = 1; k <= NrNoduri; k++)
        for(int i = 1; i <= NrNoduri; i++)
            for(int j = 1; j <= NrNoduri; j++)
                if ((G[i][j] > G[i][k] + G[k][j] || !G[i][j]) && i != j && G[i][k] && G[k][j])
                    G[i][j] = G[i][k] + G[k][j];



}

///DARB + BFS(distanta[i] - 1 la afisarea bfs)
void Graf::Bfsdarb(int nod, int* nodinit, int *dist, int distanta[])
{
    distanta[nod] = 1;

    queue<int> Q;
    bool Vizitat[NrNoduri] = {0};

    Q.push(nod);
    Vizitat[nod] = 1;

    while(!Q.empty())
    {
        nod = Q.front();
        Q.pop();

        for(auto i: Adiacenta[nod])
        {
            if(!Vizitat[i])
            {
                Q.push(i);
                Vizitat[i] = 1;
                distanta[i] = distanta[nod] + 1;
            }
        }
    }

    int distMax = -1;
    int nodinitial;

    for(int i = 1; i < NrNoduri; i++)
        if(distanta[i] > distMax)
            {
                distMax = distanta[i];
                nodinitial = i;
            }

    *nodinit = nodinitial;
    *dist = distMax;
}

bool Graf :: bfs(int rGraf[MAXmtrx][MAXmtrx], int sursa, int destinatie, int parinte[])
{
    bool Vizitat[MAX] = {0};
    queue<int> Q;

    Q.push(sursa);
    Vizitat[sursa] = 1;
    parinte[sursa] = -1;

    while(!Q.empty())
    {
        int nodcurent = Q.front();
        Q.pop();

        for (int nod = 1; nod <= NrNoduri; nod++)
        {
            if (!Vizitat[nod] && rGraf[nodcurent][nod] > 0)
            {
                Q.push(nod);
                parinte[nod] = nodcurent;
                Vizitat[nod] = 1;

                if (nod == destinatie)
                    return 1;
            }
        }
    }

    return 0;
}

int Graf :: MaxFlow(int Graf[MAXmtrx][MAXmtrx], int sursa, int destinatie)
{
    int rGraf[MAXmtrx][MAXmtrx];

    for (int i = 0; i < MAXmtrx; i++)
        for (int j = 0; j < MAXmtrx; j++)
            rGraf[i][j] = Graf[i][j];

    int totalFlow = 0;
    int parinte[MAX];

    while (bfs(rGraf, sursa, destinatie, parinte))
    {
        int drumFlow = INT_MAX;
        int j = destinatie;

        while (j != sursa)
        {
            int i = parinte[j];
            drumFlow = min(drumFlow, rGraf[i][j]);
            j = parinte[j];
        }

        j = destinatie;
        while (j != sursa)
        {
            int i = parinte[j];

            rGraf[i][j] -= drumFlow;
            rGraf[j][i] += drumFlow;

            j = parinte[j];
        }

        totalFlow += drumFlow;
    }

    return totalFlow;
}

vector<int> Graf::Euler(int nod, bool Vizitat[MAXCTC], vector<pair<int,int>> AdiacentaEuler[MAX])
{
    vector<int> Raspuns;

    stack<int> Stiva;
    Stiva.push(nod);

    while(!Stiva.empty())
    {
        nod = Stiva.top();

        if(!AdiacentaEuler[nod].empty())
        {
            int ordine = AdiacentaEuler[nod].back().first;
            int vecin = AdiacentaEuler[nod].back().second;

            AdiacentaEuler[nod].pop_back();

            if(!Vizitat[ordine])
            {
                Vizitat[ordine] = 1;
                Stiva.push(vecin);
            }
        }

        else
        {
            Raspuns.push_back(nod);
            Stiva.pop();
        }
    }

    return Raspuns;
}
int main()
{

    return 0;;
}

/*
    Solutiile sunt pe infoarena pe contul: https://www.infoarena.ro/utilizator/andreinova

    BFS (90 pct)                    - https://www.infoarena.ro/job_detail/2792577?action=view-source
    DFS (100 pct)                   - https://www.infoarena.ro/job_detail/2792927?action=view-source
    SORTARE TOPOLOGICA (100 pct)    - https://www.infoarena.ro/job_detail/2793653?action=view-source
    CTC (100 pct)                   - https://www.infoarena.ro/job_detail/2795568?action=view-source
    DIJKSTRA (100 pct)              - https://www.infoarena.ro/job_detail/2803039?action=view-source
    BELLMANFORD (100 pct)           - https://www.infoarena.ro/job_detail/2803151?action=view-source
    APM (100 pct)                   - https://www.infoarena.ro/job_detail/2803885?action=view-source
    DISJOINT (100 pct)              - https://www.infoarena.ro/job_detail/2805433?action=view-source
    ROYFLOYD (100 pct)              - https://www.infoarena.ro/job_detail/2811879?action=view-source
    DIAMETRUl ARBORELUI (100 pct)   - https://www.infoarena.ro/job_detail/2812848?action=view-source
    MAX FLOW (60 pct)               - https://www.infoarena.ro/job_detail/2814815?action=view-source
    EULER (100 pct)                 - https://www.infoarena.ro/job_detail/2820805?action=view-source

*/
