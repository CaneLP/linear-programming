# Dual Simplex Method
import numpy as np
import sys
import enum
import copy


def mprint(mat, fmt="g", new_line=True):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="   ")
        if new_line:
            print("")


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
    def __init__(self, a, b, c, inequalities, var_num, f_opt_sign, first, blend="y"):
        self.a = a
        self.initial_a = copy.deepcopy(a)
        self.b = b
        self.initial_b = copy.deepcopy(self.b)
        self.c = c
        self.initial_var_num = var_num
        self.f_opt_sign = f_opt_sign
        self.inequalities = inequalities
        self.artificial_var_num = 0
        self.blend = blend
        self.solution = None
        self.curr_table = None
        self.first = first

    def dual(self):
        # Adding slack variables
        # Dodavanje izravnavajucih promenljivih
        print("Procedure of solving the problem: ")
        print()
        mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
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

        if self.first:
            self.a = np.r_[self.a, np.c_[[self.c], np.zeros(shape=(1, self.a.shape[1] - len(self.c)))]]
            self.b.append(0)

        # print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
        # print()

        # Izvrsava se dok ima negativnih b
        while True:
            # Ukoliko su svi b pozitivni, pronadjeno je optimalno resenje
            if np.count_nonzero(np.array([self.b]) < 0) == 0:
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
                self.curr_table = np.c_[self.a, np.array(self.b).transpose()]
                self.solution = x0

                print("******************************")
                print("Input:")
                if self.f_opt_sign == -1:
                    print("(max) c = %s" % np.negative(self.c))
                else:
                    print("(min) c = %s" % self.c)
                print("subject to: ")
                for i in range(self.initial_a.shape[0]):
                    print(self.initial_a[i], end="")
                    if self.inequalities[i] == RelationSymbols.greater:
                        print(" >= ", end="")
                    elif self.inequalities[i] == RelationSymbols.less:
                        print(" <= ", end="")
                    elif self.inequalities[i] == RelationSymbols.equals:
                        print(" = ", end="")
                    print(" = ", end="")
                    print(self.initial_b[i])
                if self.blend != "y":
                    print("Blend's rule was not used!")
                else:
                    print("Blend's rule was used!")
                print("******************************")
                print("Optimal solution is: ")
                print(x0)
                if self.f_opt_sign == -1:
                    print("max_f: %s" % self.b[-1])
                else:
                    print("min_f: %s" % str(-self.b[-1]))
                print("******************************\n\n")
                return

            # Trazimo red sa negativnim b
            for i in range(len(self.b)):
                if self.b[i] < 0:
                    break
            # Ukoliko su svi clanovi matrice a tog reda nenegativni, ne postoji resenje, koje zadovoljava
            # sva ogranicenja
            if np.count_nonzero(self.a[i] >= 0) == self.a.shape[1]:
                print("Infeasible.")
                return

            # print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            # print()

            # I ovde koristimo Blendovo pravilo, kod biranja pivota iz prethodno odabranog reda.
            # Biramo poslednji maksimalni element
            max_col = 0
            maximum = float('-inf')
            if self.blend != "y":
                for j in range(self.initial_var_num):
                    if self.a[i][j] < 0 and (self.a[-1][j] / self.a[i][j]) >= maximum:
                        maximum = self.a[-1][j] / self.a[i][j]
                        max_col = j
            else:
                for j in range(self.initial_var_num):
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

            # Handle bug with very small negative values
            epsilon = 1e-10
            for rj in range(len(self.b)):
                if abs(self.b[rj]) < epsilon:
                    self.b[rj] = 0
            self.curr_table = np.c_[self.a, np.array(self.b).transpose()]
            mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()
        # print("\n\n")

