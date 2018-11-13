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
        

    def solve_problem(self):
        pass


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

        num_storages = 3
        num_shops = 4
        a = np.array([100, 200, 100])
        b = np.array([80, 120, 150, 50])
        matrix_c = np.array([[5, 4, 8, 3], [4, 7, 4, 5], [5, 3, 6, 1]])

        problem = TP(num_storages, num_shops, a, b, matrix_c)
        start_logging_to_file(output_file, mode)
        problem.solve_problem()
        stop_logging_to_file()
        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
