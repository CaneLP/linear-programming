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
    def __init__(self, a, b, c, inequalities, var_num, f_opt_sign, blend):
        self.a = a
        self.initial_a = copy.deepcopy(a)
        self.b = b
        self.c = c
        self.initial_var_num = var_num
        self.f_opt_sign = f_opt_sign
        self.inequalities = inequalities
        self.artificial_var_num = 0
        self.blend = blend

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

            # Blendovo pravilo, uvek se bira prvi po redu, koji ispunjava uslov
            # Na primer, izaberimo poslednji koji ispunjava uslov, ukoliko ova opcija nije izabrana
            # Biramo poslednje negativno b
            if self.blend != "y":
                options = []
                for i in range(len(self.b) - 1):
                    if self.b[i] < 0:
                        options.append(i)
                i = options[0]
            else:
                for i in range(len(self.b) - 1):
                    if self.b[i] < 0:
                        break

            if np.count_nonzero(self.a[i] >= 0) == self.a.shape[1]:
                print("Infeasible.")
                return

            print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()

            # I ovde koristimo Blendovo pravilo, kod biranja pivota iz prethodno odabranog reda.
            # Biramo poslednji maksimalni element
            max_col = 0
            maximum = float('-inf')
            if self.blend != "y":
                for j in range(len(self.c)):
                    if self.a[i][j] < 0 and (self.a[-1][j] / self.a[i][j]) >= maximum:
                        maximum = self.a[-1][j] / self.a[i][j]
                        max_col = j
            else:
                for j in range(len(self.c)):
                    if self.a[i][j] < 0 and (self.a[-1][j] / self.a[i][j]) > maximum:
                        maximum = self.a[-1][j] / self.a[i][j]
                        max_col = j
            pivot = self.a[i][max_col]

            for k in range(self.a.shape[0]):
                if k != i:
                    temp = self.a[k][max_col]
                    self.a[k] = self.a[k] + (temp / pivot) * -self.a[i]
                    self.b[k] = self.b[k] + (temp / pivot) * -self.b[i]
            self.a[i] = np.divide(self.a[i], pivot)
            self.b[i] = self.b[i] / pivot

            print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()


def main():
    output_file = "output.txt"
    mode = "a"
    sys.stdout = FileLogger(output_file, mode)
    ans = 'y'
    while ans == 'y':
        # Iskljucivanje stampanje eksponenta - samo zbog lepseg ispisa
        np.set_printoptions(suppress=True)
        stop_logging_to_file()

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
        # blend = input("Use Blends rule (y/n): ")

        rows = 3
        columns = 3
        matrix_a = np.array([[2, -1, -1], [1, 2, 1], [-1, 1, -2]])
        b = [[0, 3, -4]]
        c = [[7, 4, 1]]
        inequalities = [RelationSymbols.greater, RelationSymbols.greater, RelationSymbols.greater]
        f_opt_sign = 1
        blend = "n"

        problem = Dual(matrix_a, b[0], c[0], inequalities, columns, f_opt_sign, blend)
        start_logging_to_file(output_file, mode)
        problem.solve_problem()
        stop_logging_to_file()

        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
