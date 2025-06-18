#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <queue>
#include <fstream>
#undef UNICODE

using namespace std;
using namespace std::chrono;

// -------------------- ESTRUCTURAS --------------------

// Estructura: Node
// Propósito: Representa un nodo con coordenadas 2D (x, y) e identificador único.
// Atributos:
//   - id: identificador del nodo
//   - x, y: coordenadas del nodo en el plano
struct Node {
    int id;
    double x, y;
};

// Estructura: Edge
// Propósito: Representa una arista entre dos nodos con un peso (distancia).
// Atributos:
//   - u, v: nodos conectados por la arista
//   - weight: peso de la arista
// Métodos:
//   - operator<: permite ordenar aristas de menor a mayor peso
struct Edge {
    int u, v;
    double weight;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

// Estructura: CompareEdges
// Propósito: Comparador de aristas para usar en una cola de prioridad (min-heap).
// Métodos:
//   - operator(): devuelve true si 'a' tiene mayor peso que 'b', para ordenar ascendentemente.
struct CompareEdges {
    bool operator()(const Edge& a, const Edge& b) {
        return a.weight > b.weight;
    }
};

// -------------------- UNION-FIND --------------------

// Estructura: UnionFind
// Propósito: Implementa la estructura Union-Find (o Disjoint Set Union) para gestionar componentes conexas.
// Métodos:
//   - UnionFind(int n): inicializa 'n' elementos, cada uno en su conjunto.
//   - find_no_opt(int u): retorna la raíz del conjunto de 'u' sin compresión de caminos.
//   - find(int u): retorna la raíz del conjunto de 'u' con compresión de caminos (optimizado).
//   - unite(int u, int v, bool optimize): une los conjuntos que contienen a 'u' y 'v'.
//   - connected(int u, int v, bool optimize): indica si 'u' y 'v' están en el mismo conjunto.
struct UnionFind {
    vector<int> parent, size;

    UnionFind(int n) : parent(n), size(n, 1) {
        for (int i = 0; i < n; i++) parent[i] = i;
    }

    int find_no_opt(int u) {
        while (u != parent[u]) u = parent[u];
        return u;
    }

    int find(int u) {
        if (u != parent[u]) parent[u] = find(parent[u]);
        return parent[u];
    }

    void unite(int u, int v, bool optimize) {
        int root_u = optimize ? find(u) : find_no_opt(u);
        int root_v = optimize ? find(v) : find_no_opt(v);
        if (root_u == root_v) return;
        if (size[root_u] < size[root_v]) {
            parent[root_u] = root_v;
            size[root_v] += size[root_u];
        } else {
            parent[root_v] = root_u;
            size[root_u] += size[root_v];
        }
    }

    bool connected(int u, int v, bool optimize) {
        return optimize ? find(u) == find(v) : find_no_opt(u) == find_no_opt(v);
    }
};

// -------------------- FUNCIONES AUXILIARES --------------------

// Función: euclidean_distance_squared
// Propósito: Calcula la distancia euclidiana al cuadrado entre dos nodos.
// Entrada: a, b - nodos a comparar
// Salida: distancia al cuadrado entre a y b
double euclidean_distance_squared(const Node& a, const Node& b) {
    return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
}

// Función: generate_nodes
// Propósito: Genera un conjunto de N nodos con coordenadas aleatorias entre 0 y 1.
// Entrada: N - cantidad de nodos a generar
// Salida: vector con nodos generados
vector<Node> generate_nodes(int N) {
    vector<Node> nodes(N);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < N; ++i) {
        nodes[i] = {i, dis(gen), dis(gen)};
    }
    return nodes;
}

// Función: generate_all_edges
// Propósito: Genera todas las aristas posibles entre nodos, con sus respectivos pesos (distancias).
// Entrada: vector de nodos
// Salida: vector con todas las aristas generadas
vector<Edge> generate_all_edges(const vector<Node>& nodes) {
    vector<Edge> edges;
    int N = nodes.size();
    for (int i = 0; i < N; ++i) {
        for (int j = i+1; j < N; ++j) {
            double w = euclidean_distance_squared(nodes[i], nodes[j]);
            edges.push_back({i, j, w});
        }
    }
    return edges;
}

// -------------------- KRUSKAL --------------------

// Función: kruskal_sorted
// Propósito: Ejecuta el algoritmo de Kruskal usando sort en las aristas.
// Entrada: edges - lista de aristas; N - cantidad de nodos; optimize - si se usa compresión de caminos
// Salida: vector con las aristas del árbol de expansión mínima (MST)
vector<Edge> kruskal_sorted(vector<Edge> edges, int N, bool optimize) {
    sort(edges.begin(), edges.end());
    UnionFind uf(N);
    vector<Edge> mst;

    for (Edge& e : edges) {
        if (!uf.connected(e.u, e.v, optimize)) {
            uf.unite(e.u, e.v, optimize);
            mst.push_back(e);
        }
        if ((int)mst.size() == N - 1) break;
    }
    return mst;
}

// Función: kruskal_heap
// Propósito: Ejecuta Kruskal usando una min-heap para seleccionar la siguiente arista mínima.
// Entrada: edges - lista de aristas; N - cantidad de nodos; optimize - si se usa compresión de caminos
// Salida: vector con las aristas del árbol de expansión mínima (MST)
vector<Edge> kruskal_heap(vector<Edge> edges, int N, bool optimize) {
    priority_queue<Edge, vector<Edge>, CompareEdges> heap(edges.begin(), edges.end());
    UnionFind uf(N);
    vector<Edge> mst;

    while (!heap.empty() && (int)mst.size() < N - 1) {
        Edge e = heap.top(); heap.pop();
        if (!uf.connected(e.u, e.v, optimize)) {
            uf.unite(e.u, e.v, optimize);
            mst.push_back(e);
        }
    }
    return mst;
}

// -------------------- EXPERIMENTO --------------------

// Función: run_experiments
// Propósito: Ejecuta múltiples experimentos comparando variantes del algoritmo de Kruskal,
//            generando grafos completos aleatorios para diferentes valores de N.
// Entrada: ninguna
// Salida: imprime resultados por consola y guarda detalles en archivos .txt
void run_experiments() {
    cout << "N;Sorted+Opt;Sorted+NoOpt;Heap+Opt;Heap+NoOpt;Grafo(ms)\n";
    
    for (int exp = 5; exp <= 12; ++exp) {
        int N = 1 << exp;
        double sum_time[5] = {0, 0, 0, 0, 0};  // [0..3] Kruskal, [4] grafo
        double sum_weight[4] = {0, 0, 0, 0};   // pesos totales MST por método

        string filename = "grafo_N=" + to_string(N) + ".txt";
        ofstream out(filename);
        out << "==== GRAFO PARA N = " << N << " ====\n\n";

        mt19937 global_gen(chrono::steady_clock::now().time_since_epoch().count());

        for (int rep = 0; rep < 5; ++rep) {
            int seed = global_gen();  // semilla distinta por repetición
            mt19937 gen(seed);
            uniform_real_distribution<> dis(0.0, 1.0);

            auto gen_start = high_resolution_clock::now();

            // Crear nodos
            vector<Node> nodes(N);
            for (int i = 0; i < N; ++i)
                nodes[i] = {i, dis(gen), dis(gen)};
            
            // Crear aristas
            vector<Edge> edges;
            for (int i = 0; i < N; ++i) {
                for (int j = i + 1; j < N; ++j) {
                    double w = euclidean_distance_squared(nodes[i], nodes[j]);
                    edges.push_back({i, j, w});
                }
            }

            auto gen_end = high_resolution_clock::now();
            double gen_time = duration_cast<milliseconds>(gen_end - gen_start).count();
            sum_time[4] += gen_time;

            out << "Repetición #" << rep + 1 << " (semilla " << seed << ")\n";
            out << "Tiempo de generación: " << gen_time << " ms\n";

            out << "Nodos:\n";
            for (const auto& node : nodes)
                out << "  " << node.id << ": (" << node.x << ", " << node.y << ")\n";

            out << "\nAristas (" << edges.size() << "):\n";
            for (const auto& edge : edges)
                out << "  " << edge.u << " - " << edge.v << " (peso " << edge.weight << ")\n";

            auto test = [&](auto func) {
                auto start = high_resolution_clock::now();
                auto mst = func(edges, N);
                auto end = high_resolution_clock::now();
                double duration = duration_cast<std::chrono::duration<double, std::milli>>(end - start).count();

                double total_weight = 0.0;
                for (const auto& e : mst) total_weight += e.weight;

                return make_pair(duration, total_weight);
            };

            auto [t1, w1] = test([&](auto e, auto n){ return kruskal_sorted(e, n, true); });
            auto [t2, w2] = test([&](auto e, auto n){ return kruskal_sorted(e, n, false); });
            auto [t3, w3] = test([&](auto e, auto n){ return kruskal_heap(e, n, true); });
            auto [t4, w4] = test([&](auto e, auto n){ return kruskal_heap(e, n, false); });

            sum_time[0] += t1; sum_weight[0] += w1;
            sum_time[1] += t2; sum_weight[1] += w2;
            sum_time[2] += t3; sum_weight[2] += w3;
            sum_time[3] += t4; sum_weight[3] += w4;
            sum_time[4] += gen_time;

            out << "\nTiempos de ejecución (ms):\n";
            out << "  Sorted + Optimized: " << t1 << " ms (peso total: " << w1 << ")\n";
            out << "  Sorted + No Opt:    " << t2 << " ms (peso total: " << w2 << ")\n";
            out << "  Heap + Optimized:   " << t3 << " ms (peso total: " << w3 << ")\n";
            out << "  Heap + No Opt:      " << t4 << " ms (peso total: " << w4 << ")\n";
            out << "--------------------------------------------------\n\n";
        }

        out.close();

        cout << N;
        for (int i = 0; i < 5; ++i)
            cout << ";" << (sum_time[i] / 5.0);
        cout << "\n";
    }
}

// Función: main
// Propósito: Función principal. Ejecuta los experimentos de evaluación de Kruskal.
// Entrada: ninguna
// Salida: 0 al finalizar correctamente
int main() {
    run_experiments();
    return 0;
}
