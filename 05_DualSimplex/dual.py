# The Two Phase Simplex Method
import numpy as np
import sys
import enum
import copy


class RelationSymbols(enum.Enum):
    greater = 1
    less = -1
    equals = 0


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


class Dual:
    def __init__(self, a, b, c, inequalities, var_num, f_opt_sign):
        self.a = a
        self.initial_a = copy.deepcopy(a)
        self.b = b
        self.c = c
        self.initial_var_num = var_num
        self.f_opt_sign = f_opt_sign
        self.inequalities = inequalities
        self.artificial_var_num = 0

    def solve_problem(self):
        # Adding slack variables
        # Dodavanje izravnavajucih promenljivih
        print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
        print()
        num_of_slack_vars = 0
        col_to_add = np.zeros(shape=(self.a.shape[0], 1))
        for i in range(self.a.shape[0]):
            if self.inequalities[i] == RelationSymbols.greater:
                col_to_add[i][0] = -1
            elif self.inequalities[i] == RelationSymbols.less:
                col_to_add[i][0] = 1
            else:
                continue
            self.a = np.c_[self.a, col_to_add]
            num_of_slack_vars += 1
            col_to_add[i][0] = 0

        for i in range(self.a.shape[0]):
            if self.inequalities[i] == RelationSymbols.greater:
                self.a[i] = np.negative(self.a[i])
                self.b[i] = -self.b[i]

        self.a = np.r_[self.a, np.c_[[self.c], np.zeros(shape=(1, self.a.shape[1] - len(self.c)))]]
        self.b.append(0)

        print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
        print()


        while np.count_nonzero(np.array([self.b]) < 0) != 0:
            start = 0
            # Blendovo pravilo, uvek se bira prvi po redu
            for i in range(len(self.b) - 1):
                if self.b[i] < 0:
                    break

            if np.count_nonzero(np.array([self.b[:-1]]) < 0) == 0:
                # Searching for an identity matrix for optimal solution
                # Potraga za jedinicnom matricom za optimalno resenje
                eye_matrices_pos = []
                for j in range(self.a.shape[1]):
                    for e in range(self.a.shape[0]):
                        if np.all(self.a[:, j] == np.eye(self.a.shape[0])[:, e]):
                            eye_matrices_pos.append([(np.where(self.a[:, j] == 1)[0][0]), j])

                x0 = np.zeros(shape=(1, self.a.shape[1]))
                for eye_pos in eye_matrices_pos:
                    x0[0][eye_pos[1]] = self.b[eye_pos[0]]

                print("Optimal solution is: ")
                print(x0)
                print("min_f: %s" % -self.b[-1])
                return

            if np.count_nonzero(self.a[i] >= 0) == self.a.shape[1]:
                print("Infeasible.")
                return

            pivot = 0

            print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()

            max_col = 0
            maximum = float('-inf')
            for j in range(len(self.c)):
                if self.a[i][j] < 0 and (self.a[-1][j] / self.a[i][j]) > maximum:
                    maximum = self.a[-1][j] / self.a[i][j]
                    max_col = j

            pivot = self.a[i][max_col]
            # print(maximum)
            # print(pivot)
            #
            # print(i)
            # print(max_col)

            for k in range(self.a.shape[0]):
                if k != i:
                    # if pivot * self.a[k][max_col] < 0:
                    temp = self.a[k][max_col]
                    self.a[k] = self.a[k] + (temp / pivot) * -self.a[i]
                    self.b[k] = self.b[k] + (temp / pivot) * -self.b[i]
                        # print(self.a[k] + (self.a[k][max_col] / pivot) * -self.a[i])
                    # else:
                    #     self.a[k] = self.a[k] + (self.a[k][max_col] / pivot) * self.a[i]
                    #     self.b[k] = self.b[k] + (self.b[k] / pivot) * self.b[i]
                        # print(self.a[k] + (self.a[k][max_col] / pivot) * self.a[i])
            self.a[i] = np.divide(self.a[i], pivot)
            self.b[i] = self.b[i] / pivot

            print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()

    def tableau_method(self):
        # Izvrsava se sve dok postoji rj < 0
        while np.count_nonzero(self.b < 0) > 0:
            # Choosing pivot
            # Odredjivanje pivota
            minimum_col = 0
            for j in range(self.a.shape[1]):
                if self.a[-1][j] < 0:
                    minimum_col = j
                    break
            minimum = float('inf')
            minimum_row = 0
            for i in range(len(self.b)):
                if self.a[i][minimum_col] > 0 and self.b[i]/self.a[i][minimum_col] < minimum:
                    minimum = self.b[i]/self.a[i][minimum_col]
                    minimum_row = i
            pivot = self.a[minimum_row][minimum_col]
            # Izracunavanje sledece simpleks tablice
            for k in range(self.a.shape[0]):
                for m in range(self.a.shape[1]):
                    if k != minimum_row and m != minimum_col:
                        self.a[k][m] -= self.a[minimum_row][m] * self.a[k][minimum_col] / pivot
            for i in range(len(self.b)):
                if i != minimum_row:
                    self.b[i] -= self.b[minimum_row] * self.a[i][minimum_col] / pivot
            self.b[minimum_row] = self.b[minimum_row] / pivot
            self.a[minimum_row, :] = np.divide(self.a[minimum_row, :], pivot)
            for k in range(self.a.shape[0]):
                if k != minimum_row:
                    self.a[k][minimum_col] = 0
            print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()


def main():
    output_file = "output.txt"
    mode = "a"
    sys.stdout = FileLogger(output_file, mode)
    ans = 'y'
    while ans == 'y':
        np.set_printoptions(suppress=True)
        stop_logging_to_file()
        # todo begin
        # rows = int(input("Enter number of rows in the matrix A: "))
        # columns = int(input("Enter number of columns in the matrix A: "))
        #
        # matrix_a = []
        # inequalities = []
        # print("Enter the %s x %s matrix A one row by one and (in)equality signs when asked: " % (rows, columns))
        # for i in range(rows):
        #     print("Row %d: " % (i + 1))
        #     matrix_a.append(list(map(float, input().rstrip().split())))
        #     c = input("Choose (in)equality type (<=, >=, =): ")
        #     if c == "<=":
        #         inequalities.append(RelationSymbols.less)
        #     elif c == ">=":
        #         inequalities.append(RelationSymbols.greater)
        #     elif c == "=":
        #         inequalities.append(RelationSymbols.equals)
        #
        # b = []
        # print("Enter the b vector of length %s: " % rows)
        # b.append(list(map(float, input().rstrip().split())))
        #
        # matrix_a = np.array(matrix_a)
        #
        # c = []
        # print("Enter vector c, which will be optimized: ")
        # c.append(list(map(float, input().rstrip().split())))
        # opt = input("Optimization (min/max)? ")
        # f_opt_sign = 1
        # if opt == "max":
        #     f_opt_sign = -1
        #     c = np.negative(c)
        # todo end

        # print(matrix_a)
        # print(inequalities)
        # print(b[0])
        # print(c)

        # rows = 3
        # columns = 4
        # matrix_a = np.array([[0, -1, -1, 1], [2, 0, 2, 4], [1, 1, 2, 1]])
        # b = [[3, 12, 3]]
        # c = [2, 0, 3, 1]
        # inequalities = [0, 0, 0]
        # f_opt_sign = 1

        # rows = 4
        # columns = 3
        # matrix_a = np.array([[1, 3, 1], [3, 1, -1], [3, 1, 3], [1, 0, 0]])
        # b = [[10, 2, 6, 1]]
        # c = [-3, -1, -4]
        # inequalities = [RelationSymbols.less, RelationSymbols.greater, RelationSymbols.less, RelationSymbols.less]
        # f_opt_sign = -1

        # rows = 3
        # columns = 3
        # matrix_a = np.array([[2, -1, -1], [1, 2, 1], [-1, 1, -2]])
        # b = [[0, 3, -4]]
        # c = [[7, 4, 1]]
        # inequalities = [RelationSymbols.greater, RelationSymbols.greater, RelationSymbols.greater]
        # f_opt_sign = 1

        # rows = 5
        # columns = 9
        # matrix_a = np.array([[1, 1, 0, 1, 4, 0, 0, 9, 0],
        #                      [0, 4, 0, 3, 0, 0, 1, 7, 0],
        #                      [0, 5, 0, 7, 6, 0, 0, 3, 1],
        #                      [0, 9, 0, 1, 9, 1, 0, 2, 0],
        #                      [0, 4, 1, 3, 0, 0, 0, 1, 0]])
        # b = [[1, 2, 3, 4, 5]]
        # c = [[1, 2, 3, 4, 5, 6, 7, 8, 9]]
        # inequalities = [RelationSymbols.equals, RelationSymbols.equals, RelationSymbols.equals,
        #                 RelationSymbols.equals, RelationSymbols.equals]
        # f_opt_sign = 1

        # rows = 2
        # columns = 3
        # matrix_a = np.array([[1, -1, -1], [2, 1, 0]])
        # b = [[-2, 1]]
        # c = [[2, 4, 3]]
        # inequalities = [RelationSymbols.less, RelationSymbols.greater]
        # f_opt_sign = 1

        # rows = 3
        # columns = 2
        # matrix_a = np.array([[3, 1], [4, 3], [1, 2]])
        # b = [[3, 6, 3]]
        # c = [[2, 1]]
        # inequalities = [RelationSymbols.greater, RelationSymbols.greater, RelationSymbols.less]
        # f_opt_sign = 1

        rows = 3
        columns = 3
        matrix_a = np.array([[3, 2, 1], [2, 3, 3], [1, 1, -1]])
        b = [[10, 15, 4]]
        c = [[2, 3, 4]]
        inequalities = [RelationSymbols.less, RelationSymbols.less, RelationSymbols.greater]
        f_opt_sign = -1

        # rows = 2
        # columns = 3
        # matrix_a = np.array([[3, -1, 2], [4, 2, -1]])
        # b = [[-1, 5]]
        # c = [[9, 1, 1]]
        # inequalities = [RelationSymbols.greater, RelationSymbols.greater]
        # f_opt_sign = 1

        problem = Dual(matrix_a, b[0], c[0], inequalities, columns, f_opt_sign)
        start_logging_to_file(output_file, mode)
        problem.solve_problem()
        stop_logging_to_file()

        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
