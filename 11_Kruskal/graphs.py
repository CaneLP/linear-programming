# Oriented BFS and DFS
from collections import defaultdict
from collections import deque
import copy


class Graph:

    def __init__(self, vertices):
        self.num_vertices = vertices
        self.graph = defaultdict(list)
        self.marked = {}

    def print_graph(self):
        print("Graph:")
        for key in self.graph.keys():
            print("Node %s is connected with nodes: " % key, end="")
            for value in self.graph[key]:
                print(str(value), end="   ")
            print()

    # Dodavanje jednosmerne grane u graf
    def add_edge(self, u, v):
        self.graph[u].append(v)

    # Brisanje jednosmerne grane iz grafa
    def remove_edge(self, u, v):
        for index, key in enumerate(self.graph[u]):
            if key == v:
                self.graph[u].pop(index)

    def reset_marks(self):
        self.marked = {}

    # Nalazenje puta izmedju cvorova start i stop pomocu pretrage u sirinu
    def BFS(self, start, stop):

        # Ponistavamo tekuce oznake cvorova
        self.reset_marks()

        # Polazni cvor (start) se oznacava i dodaje u red
        self.marked[start] = True
        queue = deque([start])

        # Mapa parents cuva roditelje cvorova
        parents = {}
        for v in self.graph:
            parents[v] = None

        # Dok red nije prazan uklanjamo cvor sa pocetka reda i dodajemo u red njegove neoznacene susede
        while len(queue) > 0:
            v = queue.popleft()

            # Ako je na vrhu reda ciljni cvor (stop), rekonstruise se putanja do njega
            if v == stop:

                path = deque([])

                while parents[v] is not None:
                    path.appendleft(v)
                    v = parents[v]

                # Dodajemo polazni cvor
                path.appendleft(v)

                return list(path)

            # U suprotnom, oznacavamo i dodajemo u red sve neoznacene susede cvora sa pocetka reda
            for w in self.graph[v]:
                if w not in self.marked:
                    self.marked[w] = True
                    parents[w] = v
                    queue.append(w)

    # Nalazenje puta izmedju cvorova start i stop pomocu pretrage u dubinu
    def DFS(self, start, stop):

        # Ponistavamo tekuce oznake cvorova
        self.reset_marks()

        # Polazni cvor se oznacava i dodaje na putanju
        self.marked[start] = True
        path = [start]

        # Dok putanja nije prazna uzima se cvor sa vrha steka, koji predstavlja tekucu putanju
        # i dodaje njegov prvi neoznaceni sused na stek
        while len(path) > 0:
            v = path[-1]

            # Ukoliko je na vrhu steka ciljni cvor (stop), vraca se sadrzaj steka koji predstavlja
            # put od cvora start do cvora stop
            if v == stop:
                return path

            has_unvisited = False

            for w in self.graph[v]:
                if w not in self.marked:
                    path.append(w)
                    self.marked[w] = True
                    has_unvisited = True
                    break

            # Ako cvor nema vise neoznacenih suseda uklanja se sa steka
            if has_unvisited is False:
                path.pop()


def main():
    first_test = True
    g = None
    while True:
        input_from_file = int(input("Input from command line or from file? (1: command line; 2: file) "))
        if input_from_file == 2:
            file_name = input("Filename? ")
            file = open(file_name, "r")

            line_num = 0
            for line in file:
                if line_num == 0:
                    num_vertices = int(line.split()[0])
                    g = Graph(num_vertices)
                if line:
                    data = line.split()
                    curr_vertex = int(data[0])
                    data.reverse()
                    data.pop()
                    for v in data:
                        g.add_edge(curr_vertex, int(v))
                line_num += 1
        else:
            if not first_test:
                same_graph = input("Do you want to use the same graph? (y/n) ")
                if same_graph == 'n':
                    num_vertices = int(input("Number of vertices: "))
                    g = Graph(num_vertices)
            else:
                num_vertices = int(input("Number of vertices: "))
                g = Graph(num_vertices)
            while True:
                add_edge = input("Add more edges to graph? (y/n) ")
                if add_edge == 'n':
                    break
                v = int(input("Edge start node: "))
                u = int(input("Edge end node: "))
                g.add_edge(v, u)
                print("Edge (%d, %d) added." % (v, u))
                two_way = input("Add edge with other direction %d -> %d to graph? (y/n) " % (u, v))
                if two_way == 'y':
                    g.add_edge(u, v)

        test_algorithm = int(input("Test: 1) BFS; 2) DFS? "))
        if test_algorithm == 1:
            start = int(input("Start BFS from node: "))
            end = int(input("End BFS in node: "))
            g_temp = copy.deepcopy(g)
            print("********************")
            solution = g_temp.BFS(start, end)
            if solution is None:
                print("These two vertices are not connected")
            else:
                for i in range(len(solution)):
                    if i == len(solution) - 1:
                        print(solution[i])
                    else:
                        print("%s -> " % solution[i], end="")
            print("********************\n")
        elif test_algorithm == 2:
            start = int(input("Start DFS from node: "))
            end = int(input("End DFS in node: "))
            g_temp = copy.deepcopy(g)
            print("********************")
            solution = g_temp.BFS(start, end)
            if solution is None:
                print("These two vertices are not connected")
            else:
                for i in range(len(solution)):
                    if i == len(solution) - 1:
                        print(solution[i])
                    else:
                        print("%s -> " % solution[i], end="")
            print("********************\n")

        ans = input("Continue testing? (y/n) ")
        if ans == 'n':
            break

        first_test = False


if __name__ == '__main__':
    main()
