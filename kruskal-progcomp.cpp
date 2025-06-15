#include <bits/stdc++.h>
using namespace std;

typedef pair<int,int> ii;
typedef pair<int,ii> iii;
// Toma una lista de adyacencia con pares (vecino,peso)
// Retorna una lista de adyacencia con las aristas del MST
vector< vector<ii> > kruskal(vector< vector<ii> > &gr){
	int n = gr.size();
	vector< vector<ii> > ans(n);
	// Creamos una lista de aristas y agregamos todas las aristas del grafo
	vector<iii> edges;
	for (int i=0;i<n;i++){
		for (int j=0;j<gr[i].size();j++){
			edges.emplace_back(gr[i][j].second,i,gr[i][j].first);
		}
	}
	// Ordenamos las aristas por peso de menor a mayor
	sort(edges.begin(),edges.end());
	// Creamos nuestro Union Find
	ufset uf(n);
	for (int i=0;i<edges.size();i++){
		int repa = uf.findp(edges[i].second.first);
		int repb = uf.findp(edges[i].second.second);
		// Si la arista conecta dos nodos que están en conjuntos distintos
		if (repa != repb){
			// Agregamos la arista a la respuesta
			ans[edges[i].second.first].emplace_back(edges[i].second.second,edges[i].first);
			ans[edges[i].second.second].emplace_back(edges[i].second.first,edges[i].first);
			// Y unimos los dos conjuntos
			uf.uni(repa,repb);
		}
	}
	return ans;
}
