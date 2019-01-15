# Fleury's algorithm, BFS and DFS
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

    # Dodavanje grane u graf
    def add_edge(self, u, v):
        if v in self.graph[u] or u in self.graph[v]:
            return
        else:
            self.graph[u].append(v)
            self.graph[v].append(u)

    # Brisanje grane iz grafa
    def remove_edge(self, u, v):
        for index, key in enumerate(self.graph[u]):
            if key == v:
                self.graph[u].pop(index)
        for index, key in enumerate(self.graph[v]):
            if key == u:
                self.graph[v].pop(index)

    # DFS obilazak za proveru povezanosti cvorova nenula stepena
    def DFS_traversal(self, v, t, visited, print_nodes):
        visited[v] = True
        if print_nodes == 'y':
            print(v, end=" ")
        if v == t:
            return True
        for i in self.graph[v]:
            if visited[i] is False:
                if print_nodes == 'y':
                    self.DFS_traversal(i, t, visited, 'y')
                else:
                    self.DFS_traversal(i, t, visited, 'n')
        return False

    # Prebrojavanje cvorova dostupnih iz nekog cvora v
    def count_reachable(self, v, visited):
        count = 1
        visited[v] = True
        for i in self.graph[v]:
            if visited[i] is False:
                count = count + self.count_reachable(i, visited)
        return count

    # Funkcija koja proverava da li su svi cvorovi nenula stepena povezani
    def is_connected(self):
        # Resetuj oznake o posecenosti cvorova
        visited = [False] * self.num_vertices

        # Cvor nenula stepena
        for i in range(self.num_vertices):
            if len(self.graph[i]) > 1:
                break

        # Ako nema grana u grafu
        if i == self.num_vertices - 1:
            return True

        # Obidji graf iz cvora nenula stepena
        self.DFS_traversal(i, -1, visited, print_nodes='n')
        # Ako su svi cvorovi nenula stepena poseceni vrati True
        for i in range(self.num_vertices):
            if visited[i] is False and len(self.graph[i]) > 0:
                return False
        return True

    def is_eulerian(self):
        # Provera da li su svi cvorovi nenula stepena povezani
        if not self.is_connected():
            return 0
        else:
            # Prebroji cvorove neparnog stepena
            odd = 0
            for i in range(self.num_vertices):
                if len(self.graph[i]) % 2 != 0:
                    odd += 1

            # Ako je broj cvorova: 2 - postoji Ojlerov put, 0 - postoji Ojlerov cikl, > 2 - ni put ni cikl
            if odd == 0:
                return 2
            elif odd == 2:
                return 1
            elif odd > 2:
                return 0

    # Funkcija za testiranje
    def solve_fleury(self):
        res = self.is_eulerian()
        if res == 0:
            print("**************************************************************")
            print("Graph is not Eulerian - Euler path and Euler cycle don't exist")
            print("**************************************************************")
            self.print_graph()
            print("###############\n")
        elif res == 1:
            print("***********************")
            print("Graph has a Euler path:")
            print("***********************")
            self.print_graph()
            print("###############")
            print("Euler path is:")
            print("###############")
            self.print_euler()
            print("###############\n")
        else:
            print("***********************")
            print("Graph has a Euler cycle")
            print("***********************")
            self.print_graph()
            print("###############")
            print("Euler cycle is:")
            self.print_euler()
            print("###############\n")

    # Funkcija koja proverava da li grana moze uci u Ojlerov put/cikl
    def is_valid_edge(self, u, v):
        # Ako je v jedini sused cvora u - most
        if len(self.graph[u]) == 1:
            return True
        else:
            # Ako ima vise suseda
            # Broj cvorova do kojih se moze stici iz u
            visited = [False] * self.num_vertices
            num_reachable1 = self.count_reachable(u, visited)

            # Ukloni granu (u, v), pa ponovo izbroji cvorove do kojih se moze stici iz u
            self.remove_edge(u, v)
            visited = [False] * self.num_vertices
            num_reachable2 = self.count_reachable(u, visited)

            # Vrati uklonjenu granu u graf
            self.add_edge(u, v)

            # Ako je if broj dohvatljivih cvorova veci u prvom brojanju znaci da je grana (u, v) most
            if num_reachable1 > num_reachable2:
                return False
            else:
                return True

    def print_euler_from_node(self, u):
        # Obidji sve cvorove susedne cvoru u
        for v in self.graph[u]:
            # Ako je grana v validna, odstampaj, izbaci je iz grafa i obidji susede cvora v
            if self.is_valid_edge(u, v):
                print("%d-%d " % (u, v)),
                self.remove_edge(u, v)
                self.print_euler_from_node(v)

    def print_euler(self):
        # Pronadji cvor neparnog stepena
        u = 0
        for i in range(self.num_vertices):
            if len(self.graph[i]) % 2 != 0:
                u = i
                break
        # Odstampaj polazeci od pronadjenog cvora u
        self.print_euler_from_node(u)

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

        test_algorithm = int(input("Test: 1) Fleury's algorithm; 2) BFS; 3) DFS? "))
        if test_algorithm == 1:
            g_temp = copy.deepcopy(g)
            g_temp.solve_fleury()
        elif test_algorithm == 2:
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
        elif test_algorithm == 3:
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
