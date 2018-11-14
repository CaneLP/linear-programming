# Transportation Problem

# Revised Simplex Method with Eta factorization
import numpy as np
import sys
import copy


class FileLogger(object):
    def __init__(self, filename, mode):
        self.terminal = sys.stdout
        self.logfile = open(filename, mode)

    def write(self, message):
        self.terminal.write(message)
        self.logfile.write(message)

    def flush(self):
        pass


def start_logging_to_file(filename, mode):
    sys.stdout = FileLogger(filename, mode)


def stop_logging_to_file():
    sys.stdout.logfile.close()
    sys.stdout = sys.stdout.terminal


class TP:
    def __init__(self, num_storages, num_shops, a, b, matrix_c):
        self.m = num_storages
        self.n = num_shops
        self.a = a
        self.b = b
        self.c = matrix_c
        self.x = self.initial_solution()

    # Getting initial solution using minimum-cost method
    # Trazenje pocetnog resenja metodom minimalnih cena
    def initial_solution(self):
        # If sum a_i is equal to sum b_i - closed transportation problem, otherwise introduce fictitious storage row
        sum_a = np.sum(self.a)
        sum_b = np.sum(self.b)
        if sum_a < sum_b:
            row_to_add = np.full((1, self.c.shape[1]), int(np.mean(self.c)) * 10)
            self.c = np.r_[self.c, row_to_add]
            self.a = np.append(self.a, [sum_b - sum_a])
        elif sum_a > sum_b:
            col_to_add = np.full((1, self.c.shape[0]), int(np.mean(self.c)) * 10)
            self.c = np.c_[self.c, col_to_add.transpose()]
            self.b = np.append(self.b, [sum_a - sum_b])

        # Indicator that shows whether the row/column should be considered when searching for min c_i,j
        indicator = np.zeros(shape=(self.c.shape[0], self.c.shape[1]))
        x = np.zeros(shape=(self.c.shape[0], self.c.shape[1]))

        a = copy.deepcopy(self.a)
        b = copy.deepcopy(self.b)

        while True:
            # Find row (p) and column (q) of the current available minimum element and assign new values to x, a, b
            p = 0
            q = 0
            minimum_c = np.inf
            for i in range(self.c.shape[0]):
                for j in range(self.c.shape[1]):
                    if self.c[i][j] < minimum_c and indicator[i][j] == 0:
                        p = i
                        q = j
                        minimum_c = self.c[p][q]

            min_ab = min(a[p], b[q])
            x[p][q] = min_ab
            a[p] -= min_ab
            b[q] -= min_ab

            if a[p] == 0:
                indicator[p] = -1
            else:
                indicator[:, q] = -1

            if a[p] == b[q]:
                break

        return x

    # Positions in the matrix are encoded with a number whose digits represent positions, a = ij (4: i = 0, j = 4)
    def calculate_adjacency_list(self):
        adj_list = {}
        for i in range(self.x.shape[0]):
            for j in range(self.x.shape[1]):
                if self.x[i][j] != 0:
                    row_neighbours = [(i * 10 + k) for k in range(self.x.shape[1]) if self.x[i][k] != 0 and k != j]
                    col_neighbours = [(k * 10 + j) for k in range(self.x.shape[0]) if self.x[k][j] != 0 and k != i]
                    adj_list[i * 10 + j] = row_neighbours + col_neighbours

        return adj_list

    # Returns encoded values of positions where theta or theta negative should be added. Using DFS for cycle searching
    @staticmethod
    def find_theta_cycle(adj_list, marked, start):
        path = [start]
        marked[start] = True
        while len(path) > 0:
            v = path[-1]

            if v == start and len(path) > 3:
                return path

            has_unvisited = False

            for w in adj_list[v]:
                available = True
                # Disable stepping over already selected basic variable.
                wi = int(w / 10)
                wj = w % 10
                vi = int(v / 10)
                vj = v % 10
                for p in path:
                    if p != v and p != w:
                        pi = int(p / 10)
                        pj = p % 10
                        if ((vi == pi == wi and (vj < pj < wj or wj < pj < vj)) or
                           (vj == pj == wj and (wi < pi < vi or vi < pi < wi))) and w not in marked:
                            # print("v")
                            # print(v)
                            # print("p")
                            # print(p)
                            # print("w")
                            # print(w)
                            # print()
                            available = False
                if len(path) > 1:
                    p = path[-2]
                    pi = int(p / 10)
                    pj = p % 10
                    if ((vi == pi == wi and (vj < wj < pj or pj < wj < vj)) or
                       (vj == pj == wj and (pi < wi < vi or vi < wi < pi))) and w not in marked:
                        available = False

                if available and ((w == start and len(path) > 3 and len(path) % 2 == 0) or w not in marked):
                    path.append(w)
                    marked[w] = True
                    has_unvisited = True
                    break

            if not has_unvisited:
                path.pop()

        return []

    def method_of_potentials(self):
        # Determining the row/column with maximum number of basic variables
        i = 0
        max_basic_i = 0
        for k in range(self.x.shape[0]):
            num_basic = np.count_nonzero(self.x[k] > 0)
            if num_basic > max_basic_i:
                max_basic_i = num_basic
                i = k
        j = 0
        max_basic_j = 0
        for k in range(self.x.shape[1]):
            num_basic = np.count_nonzero(self.x[:, k] > 0)
            if num_basic > max_basic_j:
                max_basic_j = num_basic
                j = k

        u = np.full((1, len(self.a)), np.inf)[0]
        v = np.full((1, len(self.b)), np.inf)[0]

        # Setting that row/column of vector u/v to zero and calculating rest of u and v elements
        if max_basic_i >= max_basic_j:
            u[i] = 0
        else:
            v[j] = 0

        while np.any(u == np.inf) or np.any(v == np.inf):
            for i in range(len(u)):
                if u[i] != np.inf:
                    for j in range(len(v)):
                        if self.x[i][j] != 0 and v[j] == np.inf:
                            v[j] = self.c[i][j] - u[i]
            for j in range(len(v)):
                if v[j] != np.inf:
                    for i in range(len(u)):
                        if self.x[i][j] != 0 and u[i] == np.inf:
                            u[i] = self.c[i][j] - v[j]

        # Finding smallest value less than zero: c(non-basic)_i,j - (u_i + v_j).
        # If all greater than or equal zero, optimal solution is found
        neg_values = []
        for i in range(self.c.shape[0]):
            for j in range(self.c.shape[1]):
                if self.x[i][j] == 0:
                    curr_diff = self.c[i][j] - (u[i] + v[j])
                    if curr_diff < 0:
                        neg_values.append((curr_diff, [i, j]))

        # Optimal solution found
        if len(neg_values) == 0:
            print("Optimal solution found!")
            print(self.x)

        # Positioning theta on the smallest negative value. Creating cycle - start from theta position, go
        # through basic variables and return at the starting position.
        theta = np.zeros(shape=(self.c.shape[0], self.c.shape[1]))
        pos_i = sorted(neg_values)[0][1][0]
        pos_j = sorted(neg_values)[0][1][1]
        # theta[pos_i][pos_j] = 1
        self.x[pos_i][pos_j] = 1
        adjacency_list = self.calculate_adjacency_list()
        self.x[pos_i][pos_j] = 0

        print(adjacency_list)
        marked = {}
        start = pos_i * 10 + pos_j
        positions = TP.find_theta_cycle(adjacency_list, marked, start)[:-1]
        print(positions)
        sign = -1
        # Decode values and fill theta matrix
        for pos in positions:
            sign *= -1
            i = int(pos / 10)
            j = pos % 10
            theta[i][j] = sign
        theta[pos_i][pos_j] = 1
        print(pos_i, pos_j)
        # print(sorted(neg_values))
        print(theta)

    def solve_problem(self):
        print("Initial solution found using minimum-cost method: \n%s" % self.x)
        # print("a")
        # print(self.a)
        # print("b")
        # print(self.b)
        # print("c")
        # print(self.c)

        self.method_of_potentials()


def main():
    output_file = "output.txt"
    mode = "a"
    sys.stdout = FileLogger(output_file, mode)
    ans = 'y'
    while ans == 'y':
        stop_logging_to_file()
        # num_storages = int(input("Enter number of storages: "))
        # num_shops = int(input("Enter number of shops: "))
        #
        # a = []
        # print("Enter the a vector of length %s which contains quantity of merchandise in i-th storage: " % num_storages)
        # a.append(list(map(float, input().rstrip().split())))
        #
        # b = []
        # print("Enter the b vector of length %s which contains demand of j-th shop: " % num_shops)
        # b.append(list(map(float, input().rstrip().split())))
        #
        # matrix_c = []
        # print("Enter the %s x %s matrix C which contains transportation cost "
        #       "from i-th storage to j-th shop: " % (num_storages, num_shops))
        # for i in range(num_storages):
        #     matrix_c.append(list(map(float, input().rstrip().split())))
        # matrix_c = np.array(matrix_c)

        # print(num_storages)
        # print(num_shops)
        # print(a)
        # print(b)
        # print(matrix_c)
        # print(matrix_x)

        # num_storages = 3
        # num_shops = 4
        # a = np.array([100, 200, 100])
        # b = np.array([80, 120, 150, 50])
        # matrix_c = np.array([[5, 4, 8, 3], [4, 7, 4, 5], [5, 3, 6, 1]])

        num_storages = 4
        num_shops = 5
        a = np.array([28, 13, 19, 18])
        b = np.array([24, 16, 10, 20, 22])
        # matrix_c = np.array([[1, 0, 0, 0, 0], [4, 2, 0, 1, 4], [4, 5, 6, 2, 0], [0, 0, 0, 1, 0]])
        # matrix_c = np.array([[1, 0, 5, 6, 9], [4, 2, 7, 9, 4], [14, 25, 16, 21, 5], [3, 4, 2, 6, 1]])
        # matrix_c = np.array([[1, 2, 3, 4, 5], [14, 12, 10, 11, 24], [34, 15, 6, 2, 7], [8, 8, 9, 1, 2]])
        matrix_c = np.array([[3, 9, 8, 10, 4], [6, 10, 3, 2, 3], [3, 2, 7, 10, 3], [3, 2, 3, 2, 8]])

        problem = TP(num_storages, num_shops, a, b, matrix_c)
        start_logging_to_file(output_file, mode)
        problem.solve_problem()
        stop_logging_to_file()
        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
