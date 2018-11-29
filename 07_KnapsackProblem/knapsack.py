# Knapsack Problem
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


def mprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="   ")
        print("")


class KP:
    def __init__(self, num_items, knapsack_size, weights, values, recurrence_type):
        self.num_items = num_items
        self.knapsack_size = knapsack_size
        self.weights = weights
        self.values = values
        self.recurrence_type = recurrence_type
        self.table_f = None
        self.table_i = None
        self.num_nz_items = None

    def solve_in_front(self, table_f, table_i, curr):
        if curr < self.num_nz_items:
            for y in range(self.knapsack_size + 1):
                if y - self.weights[curr] < 0:
                    second_par = float('-inf')
                else:
                    second_par = table_f[curr - 1][y - self.weights[curr]] + self.values[curr]
                table_f[curr][y] = max(table_f[curr - 1][y], second_par)
                if table_f[curr - 1][y] > second_par:
                    table_i[curr][y] = table_i[curr - 1][y]
                else:
                    table_i[curr][y] = curr + 1
            self.solve_in_front(table_f, table_i, curr + 1)
        else:
            self.table_f = table_f
            self.table_i = table_i
            return

    # def solve_backwards(self, ):
    #     if curr < self.num_nz_items:
    #         for y in range(self.knapsack_size + 1):
    #             if y - self.weights[curr] < 0:
    #                 second_par = float('-inf')
    #             else:
    #                 second_par = table_f[curr - 1][y - self.weights[curr]] + self.values[curr]
    #             table_f[curr][y] = max(table_f[curr - 1][y], second_par)
    #             if table_f[curr - 1][y] > second_par:
    #                 table_i[curr][y] = table_i[curr - 1][y]
    #             else:
    #                 table_i[curr][y] = curr + 1
    #         self.solve_in_front(table_f, table_i, curr + 1)
    #     else:
    #         self.table_f = table_f
    #         self.table_i = table_i
    #         return

    def solve_problem(self):
        if self.recurrence_type == 1:
            self.num_nz_items = np.count_nonzero(self.weights != 0)
            self.table_f = np.zeros(shape=(self.num_nz_items, self.knapsack_size + 1))
            self.table_i = np.zeros(shape=(self.num_nz_items, self.knapsack_size + 1))
            iterable = (self.values[0] * np.floor(i / self.weights[0]) for i in range(self.knapsack_size + 1))
            self.table_f[0] = np.fromiter(iterable, int)
            for i in range(self.table_f.shape[1]):
                if self.table_f[0][i] != 0:
                    self.table_i[0][i] = 1
            self.solve_in_front(self.table_f, self.table_i, 1)
            print("Table F (k\y):")
            mprint(self.table_f)
            print("Table i (i\y):")
            mprint(self.table_i)

            f_max = np.max(self.table_f[-1])
            print("f_max = %s" % f_max)
            solution_x = np.zeros(shape=(1, self.num_items))
            i, j = np.where(self.table_f == f_max)
            i = i[0]
            j = j[0]
            space_left = self.knapsack_size
            while space_left != 0:
                pos = int(self.table_i[i][j]) - 1
                solution_x[0][pos] = 1
                space_left -= self.weights[pos]
                j = space_left
            print("x = %s" % solution_x[0])

        else:
            self.solve_backwards()


def main():
    output_file = "output.txt"
    mode = "a"
    sys.stdout = FileLogger(output_file, mode)
    ans = 'y'
    while ans == 'y':
        stop_logging_to_file()
        # knapsack_size = int(input("Enter knapsack size: "))
        # num_items = int(input("Enter number of items: "))
        #
        # values = []
        # print("Enter values of given %s items: " % num_items)
        # values.append(list(map(float, input().rstrip().split())))
        #
        # weights = []
        # print("Enter weights of given %s items: " % num_items)
        # weights.append(list(map(float, input().rstrip().split())))
        #
        # recurrence_type = int(input("Recurrence formula type (1: in front/ 2: backwards)? "))


        # knapsack_size = 4
        # num_items = 4
        # values = np.array([[10, 2, -1, 3]])
        # weights = np.array([[3, 1, 3, 0]])
        # recurrence_type = 1

        # knapsack_size = 8
        # num_items = 5
        # values = np.array([[3, 1, 7, 2, 5]])
        # weights = np.array([[4, 1, 2, 3, 6]])
        # recurrence_type = 1

        knapsack_size = 9
        num_items = 4
        values = np.array([[3, 4, 5, 2]])
        weights = np.array([[2, 3, 4, 5]])
        recurrence_type = 2

        problem = KP(num_items, knapsack_size, weights[0], values[0], recurrence_type)
        start_logging_to_file(output_file, mode)
        problem.solve_problem()
        stop_logging_to_file()
        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
