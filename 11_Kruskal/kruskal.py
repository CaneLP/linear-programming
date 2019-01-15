# Kruskal's algorithm for MST


class Graph:

    def __init__(self, vertices):
        self.num_vertex = vertices
        self.graph = []
        self.adj = [[] for _ in range(vertices)]

    def add_edge(self, u, v, w):
        if v in self.adj[u] or u in self.adj[v]:
            return
        else:
            self.graph.append([u, v, w])
            self.adj[u].append(v)
            self.adj[v].append(u)

    def find(self, parent, i):
        if parent[i] == i:
            return i
        return self.find(parent, parent[i])

    def union(self, parent, rank, x, y):
        xroot = self.find(parent, x)
        yroot = self.find(parent, y)

        if rank[xroot] < rank[yroot]:
            parent[xroot] = yroot
        elif rank[xroot] > rank[yroot]:
            parent[yroot] = xroot

        else:
            parent[yroot] = xroot
            rank[xroot] += 1

    def kruskal_mst(self):

        result = []

        i = 0
        e = 0

        # Sortirati grane u neopadajucem redosledu tezina
        self.graph = sorted(self.graph, key=lambda item: item[2])

        parent = []
        rank = []

        # Napravi onoliko podskupova koliko ima cvorova
        for node in range(self.num_vertex):
            parent.append(node)
            rank.append(0)

        # U MST ulaze svi cvorovi
        while e < self.num_vertex - 1:
            # Izaberi sledecu granu (koje su sortirane) i uvecaj taj brojac
            u, v, w = self.graph[i]
            i = i + 1
            x = self.find(parent, u)
            y = self.find(parent, v)

            # Provera ciklusa... Ako dodavanjem grane ne dobijamo ciklus, ukljucujemo tu granu u rezultat MST-a
            if x != y:
                e = e + 1
                result.append([u, v, w])
                self.union(parent, rank, x, y)

        print("-----------------------------")
        print("MST for the current component")
        print("-----------------------------")
        res = 0
        for u, v, weight in result:
            print("%d - %d = %d" % (u, v, weight))
            res += weight
        print("Value: %s" % res)

    def dfs_traversal(self, temp, v, visited):
        visited[v] = True
        temp.append(v)

        for i in self.adj[v]:
            if visited[i] is False:
                temp = self.dfs_traversal(temp, i, visited)
        return temp

    def get_components(self):
        visited = []
        comps = []
        for i in range(self.num_vertex):
            visited.append(False)
        for v in range(self.num_vertex):
            if visited[v] is False:
                temp = []
                comps.append(self.dfs_traversal(temp, v, visited))
        return comps


def main():
    # same_graph = 'n'
    # num_vertices = int(input("Number of vertices: "))
    # g = Graph(num_vertices)
    # while True:
    #     while True and same_graph == 'n':
    #         add_edge = input("Add edge to graph? (y/n) ")
    #         if add_edge == 'n':
    #             break
    #         v = int(input("Starting node: "))
    #         u = int(input("End node: "))
    #         weight = int(input("Edge weight: "))
    #         g.add_edge(v, u, weight)
    #         print("Edge (%d, %d) with weight %s added." % (v, u, weight))
    #
    #     print("***********************")
    #     print("* Kruskal's algorithm *")
    #     print("***********************\n")
    #     print("Given graph:")
    #     print(g.graph)
    #     print()
    #
    #     print("Graph components:")
    #     components = g.get_components()
    #     print(components)
    #     print()
    #     for i in range(len(components)):
    #         print("#####################")
    #         print("##Current component##")
    #         print("#####################")
    #         print(components[i])
    #         g1 = Graph(len(components[i]))
    #         for u, v, weight in g.graph:
    #             if u in components[i]:
    #                 # print(u, v, weight)
    #                 g1.add_edge(u % len(components[i]), v % len(components[i]), weight)
    #         g1.kruskal_mst()
    #
    #     ans = input("Continue testing? (y/n) ")
    #     if ans == 'n':
    #         break
    #     same_graph = input("Do you want to use the same graph? (y/n) ")
    #     if same_graph == 'n':
    #         num_vertices = int(input("Number of vertices: "))
    #         g = Graph(num_vertices)

    g = Graph(4)
    g.add_edge(0, 1, 10)
    g.add_edge(0, 2, 6)
    g.add_edge(0, 3, 5)
    g.add_edge(1, 3, 15)
    g.add_edge(2, 3, 4)
    g.add_edge(3, 2, 4)
    g.add_edge(3, 1, 15)

    g.kruskal_mst()

    g = Graph(4)
    g.add_edge(0, 1, 10)
    g.add_edge(0, 2, 6)
    g.add_edge(0, 3, 5)
    g.add_edge(1, 3, 15)
    g.add_edge(2, 3, 4)
    g.add_edge(3, 2, 4)
    g.add_edge(3, 1, 15)

    g.kruskal_mst()

    # g = Graph(7)
    # g.add_edge(0, 1, 10)
    # g.add_edge(0, 2, 6)
    # g.add_edge(1, 3, 5)
    # g.add_edge(4, 5, 15)
    # g.add_edge(5, 6, 4)


if __name__ == '__main__':
    main()
