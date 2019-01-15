# Gomory Method
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


class Gomory:
    def __init__(self, a, b, c, inequalities, var_num, f_opt_sign):
        self.a = a
        self.initial_a = copy.deepcopy(a)
        self.b = b
        self.initial_b = copy.deepcopy(b)
        self.c = c
        self.initial_c = copy.deepcopy(c)
        self.var_num = var_num
        self.initial_var_num = copy.copy(var_num)
        self.f_opt_sign = f_opt_sign
        self.inequalities = inequalities
        self.artificial_var_num = 0
        self.solution = None
        self.curr_table = None
        self.first = None
        self.blend = "y"
        self.slack_var_num = None

    # def __init__(self, a, b, c, inequalities, var_num, f_opt_sign, first):
    #     self.a = a
    #     self.initial_a = copy.deepcopy(a)
    #     self.b = b
    #     self.initial_b = copy.deepcopy(b)
    #     self.c = c
    #     self.initial_var_num = var_num
    #     self.f_opt_sign = f_opt_sign
    #     self.inequalities = inequalities
    #     self.artificial_var_num = 0
    #     self.slack_var_num = 0
    #     self.curr_table = None
    #     self.solution = None
    #     self.first = first
    #     self.f_pos = None

    def dual(self):
        # Adding slack variables
        # Dodavanje izravnavajucih promenljivih

        print("*****************************************")
        print("**************DUAL SIMPLEX***************")
        print("*****************************************")
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

                # print("******************************")
                # print("Input:")
                # if self.f_opt_sign == -1:
                #     print("(max) c = %s" % np.negative(self.c))
                # else:
                #     print("(min) c = %s" % self.c)
                # print("subject to: ")
                # for i in range(self.initial_a.shape[0]):
                #     print(self.initial_a[i], end="")
                #     if self.inequalities[i] == RelationSymbols.greater:
                #         print(" >= ", end="")
                #     elif self.inequalities[i] == RelationSymbols.less:
                #         print(" <= ", end="")
                #     elif self.inequalities[i] == RelationSymbols.equals:
                #         print(" = ", end="")
                #     print(" = ", end="")
                #     print(self.initial_b[i])
                # if self.blend != "y":
                #     print("Blend's rule was not used!")
                # else:
                #     print("Blend's rule was used!")
                print("******************************")
                print("Optimal solution is: ")
                print("x0 = %s" % x0)
                # if self.f_opt_sign == -1:
                #     print("max_f: %s" % self.b[-1])
                # else:
                #     print("min_f: %s" % str(-self.b[-1]))
                print("******************************")
                return

            # Trazimo red sa negativnim b
            for i in range(len(self.b)):
                if np.round(self.b[i], decimals=10) < 0:
                    break
            # print(self.b)
            # print(i)

            # Ukoliko su svi clanovi matrice a tog reda nenegativni, ne postoji resenje, koje zadovoljava
            # sva ogranicenja
            if np.count_nonzero(self.a[i] >= 0) == self.a.shape[1]:
                self.solution = None
                print("Infeasible.")
                return

            # print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            # print()
            # print('pre petlje', self.a)
            # print(self.b)

            # I ovde koristimo Blendovo pravilo, kod biranja pivota iz prethodno odabranog reda.
            # Biramo poslednji maksimalni element
            max_col = 0
            maximum = float('-inf')
            # if self.blend != "y":
            for j in range(self.a.shape[1]):
                if np.round(self.a[i][j], decimals=10) < 0 and (self.a[0][j] / self.a[i][j]) >= maximum:
                    maximum = self.a[0][j] / self.a[i][j]
                    max_col = j
            # print('posle petlje', self.a)
            # print(self.b)
            # else:
            #     for j in range(self.var_num):
            #         if np.round(self.a[i][j], decimals=10) < 0 and (self.a[0][j] / self.a[i][j]) > maximum:
            #             maximum = self.a[0][j] / self.a[i][j]
            #             max_col = j
            pivot = self.a[i][max_col]
            # print(pivot)
            for k in range(self.a.shape[0]):
                if k != i:
                    temp = self.a[k][max_col]
                    self.a[k] = self.a[k] + (temp / pivot) * -self.a[i]
                    self.b[k] = self.b[k] + (temp / pivot) * -self.b[i]
            self.a[i] = np.divide(self.a[i], pivot)
            self.b[i] = self.b[i] / pivot

            self.a = np.round(self.a, decimals=10)
            self.b = np.round(self.b, decimals=10)
            # print("round")
            # print(self.a)
            # print(self.b)

            # Handle bug with very small negative values
            epsilon = 1e-10
            for rj in range(len(self.b)):
                if abs(self.b[rj]) < epsilon:
                    self.b[rj] = 0
            self.curr_table = np.c_[self.a, np.array(self.b).transpose()]
            mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()
            print(np.around(self.solution, decimals=10))
        # print("\n\n")

    def two_phase(self):
        if not self.phase_one_simplex():
            print("The original problem is infeasible!\n\n")
            self.solution = None
            return
        self.phase_two_simplex()

    def tableau_method(self):
        epsilon = 1e-10
        # Izvrsava se sve dok postoji rj < 0
        while np.count_nonzero(self.a[0] < 0) > 0:
            # Choosing pivot
            # Odredjivanje pivota
            minimum_col = 0
            for j in range(self.a.shape[1]):
                if self.a[0][j] < 0:
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

            # Handle bug with very small negative values
            for rj in range(self.a.shape[1]):
                if abs(self.a[0][rj]) < epsilon:
                    self.a[0][rj] = 0

            self.curr_table = np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=40)
            mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=40))
            print()

    def phase_one_simplex(self):
        print("*****************************************")
        print("***********TWO PHASE SIMPLEX*************")
        print("*****************************************")
        print("--------------- PHASE ONE ---------------")
        print("*****************************************")
        # Adding slack variables
        # Dodavanje izravnavajucih promenljivih
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
        self.slack_var_num = self.a.shape[1] - self.initial_var_num
        # Adding artificial variables
        # Dodavanje vestackih promenljivih
        artificial_var_rows = []
        self.artificial_var_num = 0
        for i in range(self.a.shape[0]):
            curr_slack_position = np.where(self.a[i,
                                           self.initial_var_num: self.initial_var_num + num_of_slack_vars] != 0)
            if len(curr_slack_position[0]) == 0 \
                    or self.a[i][self.initial_var_num + curr_slack_position[0][0]] * self.b[i] < 0:
                self.artificial_var_num += 1
                artificial_var_rows.append(i)
                col_to_add[i][0] = np.sign(self.b[i])
                self.a = np.c_[self.a, col_to_add]
                col_to_add[i][0] = 0
        # print(self.a)
        # Add new row for minimization of new function - sum of artificial variables
        # Dodavanje novog reda za minimizaciju nove funkcije - sume vestackih promenljivih
        self.a = np.r_[self.a, np.c_[np.zeros(shape=(1, self.a.shape[1] - self.artificial_var_num)),
                                     np.ones(shape=(1, self.artificial_var_num))]]

        mprint(np.c_[self.a, np.c_[[self.b], np.zeros(shape=(1, 1))].transpose()])
        print()
        # Calculate starting matrix for the tableau method
        # Izracunavanje pocetne matrice za tablicni metod, tj. dodavanje poslednjeg reda
        row_to_add = np.zeros(shape=(1, self.a.shape[1]))
        b_to_add = 0
        for row in artificial_var_rows:
            b_to_add += self.b[row]
            row_to_add += self.a[row]
        row_to_add[:, self.a.shape[1] - self.artificial_var_num:] = 0
        row_to_add = np.negative(row_to_add)
        # print(row_to_add)
        # np.insert(self.a, 0, row_to_add, 0)
        self.a[-1, :] = row_to_add
        # self.b.insert(0, -b_to_add)
        self.b.append(-b_to_add)
        self.a[[0, -1]] = self.a[[-1, 0]]
        self.b[0], self.b[-1] = self.b[-1], self.b[0]
        mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
        print()

        self.tableau_method()

        if np.around(self.b[0], decimals=5) != 0:
            # print(self.b[0])
            return False
        else:
            return True

    def phase_two_simplex(self):
        print("*****************************************")
        print("--------------- PHASE TWO ---------------")
        print("*****************************************")
        # print(self.a)
        self.a = np.delete(self.a, 0, 0)
        # print(self.a)

        # Remove nonbasic artificial columns
        # Izbacivanje vestackih nebazisnih kolona
        for j in range(self.a.shape[1] - self.artificial_var_num, self.a.shape[1]):
            for e in range(self.a.shape[0]):
                if np.all(self.a[:, j] != np.eye(self.a.shape[0])[:, e]):
                    self.a[:, j] = np.zeros(shape=(self.a.shape[0], 1))[0]

        self.a = self.a[:, ~np.all(self.a == 0, axis=0)]

        if self.first:
            # Searching for an identity matrix before adding row for the function f
            # Potraga za jedinicnom matricom, pre dodavanja poslednje vrste, funkcije f
            eye_matrices_pos = []
            for j in range(self.a.shape[1]):
                for e in range(self.a.shape[0]):
                    if np.all(self.a[:, j] == np.eye(self.a.shape[0])[:, e]):
                        eye_matrices_pos.append([(np.where(self.a[:, j] == 1)[0][0]), j])

            # Add row for initial f
            # Dodavanje reda funkcije ciji optimum trazimo
            # print(self.a)
            # print(self.c)
            self.a = np.r_[np.c_[[self.c], np.zeros(shape=(1, self.a.shape[1] - len(self.c)))], self.a]

            mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
            print()

            # Cancel columns of f on positions which correspond to basic columns
            # Ponistavanje kolona od f na pozicijama koje odgovaraju bazisnim kolonama
            for eye_pos in eye_matrices_pos:
                self.b[0] -= (self.b[eye_pos[0]] * self.a[0][eye_pos[1]])
                self.a[0] -= (self.a[eye_pos[0]] * self.a[0][eye_pos[1]])

        # Izbaci sve vestacke
        while self.a.shape[1] > self.initial_var_num + self.slack_var_num:
            if np.where(self.a[:, -1] == 1)[0]:
                row_one = np.where(self.a[:, -1] == 1)[0][0]
                # Kada su u vrsti sve nule
                if np.count_nonzero(self.a[row_one] == 0) == (self.a.shape[1] - 1):
                    self.a = np.delete(self.a, self.a.shape[1] - 1, 1)
                    self.a = np.delete(self.a, row_one, 0)
                    self.b = np.delete(self.b, row_one, 0)
            # Kada nisu sve nule
            else:
                for i in range(self.a.shape[0]):
                    for j in range(self.initial_var_num + self.slack_var_num, self.a.shape[1]):
                        if self.a[i][j] == 1:
                            column_r = np.where(self.a[i] != 0)
                            pivot = self.a[i][column_r[0][0]]
                            col = column_r[0][0]
                            for k in range(self.a.shape[0]):
                                for m in range(self.a.shape[1]):
                                    if k != i and m != col:
                                        self.a[k][m] -= self.a[i][m] * self.a[k][col] / pivot
                            for o in range(len(self.b)):
                                if o != i:
                                    self.b[o] -= self.b[i] * self.a[o][col] / pivot
                            self.b[i] = self.b[i] / pivot
                            self.a[i, :] = np.divide(self.a[i, :], pivot)
                            for k in range(self.a.shape[0]):
                                if k != i:
                                    self.a[k][col] = 0
                        self.a = np.delete(self.a, j, 1)
                        # print(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
                        # print()
                        break

        # print("+++++++++++++++++++++++++++")
        # print(self.a)
        # print(np.array(self.b).transpose())
        self.curr_table = np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=40)
        mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
        print()

        self.tableau_method()

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

        self.solution = x0

        # print("******************************")
        # print("Input:")
        # if self.f_opt_sign == -1:
        #     print("(max) c = %s" % np.negative(self.c))
        # else:
        #     print("(min) c = %s" % self.c)
        # print("subject to: ")
        # for i in range(self.initial_a.shape[0]):
        #     mprint(np.array([self.initial_a[i]]), new_line=False)
        #     if self.inequalities[i] == RelationSymbols.greater:
        #         print(" >= ", end="")
        #     elif self.inequalities[i] == RelationSymbols.less:
        #         print(" <= ", end="")
        #     elif self.inequalities[i] == RelationSymbols.equals:
        #         print(" = ", end="")
        #     print(self.initial_b[i])
        print("******************************")
        print("Optimal solution is: ")
        print("x0 = %s" % x0[0])
        # if self.f_opt_sign == -1:
        #     print("max_f: %s" % self.b[-1])
        # else:
        #     print("min_f: %s" % str(-self.b[-1]))
        print("******************************")

    def is_integer_solution(self):
        eps = 1e-6

        for x in self.solution:
            if np.modf(x)[0] > eps:
                return False
        return True

    def use_tps(self):
        cond1 = np.count_nonzero(np.array([self.b]) >= 0) == len(self.b)
        slack_signs = np.array([i.value for i in self.inequalities])
        b_signs = np.array([np.sign(i) for i in self.b])
        for i in range(self.a.shape[0]):
            if b_signs[i] < 0 and self.inequalities[i] == RelationSymbols.equals.value:
                return False
        return True

    def solve_problem(self):
        num_iter = 0
        self.first = True
        print("###########################################")
        print("############## Gomory method ##############")
        print("###########################################")
        print()
        while True:
            if self.use_tps() and self.first:
                # print("###########################################")
                # print("In this step Two-Phase Simplex will be used")
                # print("###########################################")
                # problem = tps.TPSM(self.a, self.b, self.c, self.inequalities, self.var_num, self.f_opt_sign)
                # problem.two_phase()
                self.two_phase()
                if self.solution is None:
                    # self.dual()
                    # if self.solution is None:
                    return
                # self.curr_table = problem.curr_table
                # self.solution = problem.solution[0][:problem.solution.shape[1] -
                #                                     (problem.solution.shape[1] - self.initial_var_num)]
            else:
                self.dual()
                if self.solution is None:
                    return
                # print("######################################")
                # print("In this step Dual Simplex will be used")
                # print("######################################")
                # problem = du.Dual(self.a, self.b, self.c, self.inequalities, self.var_num, self.f_opt_sign, first)
                # problem.dual()
                # if problem.solution is None:
                #     return
                # self.solution = problem.solution[0][:problem.solution.shape[1] -
                #                                     (problem.solution.shape[1] - self.initial_var_num)]
                # self.curr_table = problem.curr_table
            i = 0
            self.b = self.curr_table[:, -1]
            for x in range(len(self.b)):
                if np.around(np.modf(self.b[x])[0], decimals=10) != 0:
                    i = x
                    break

            self.a = self.curr_table[:, :-1]
            # Just to be consistent with notebook examples
            # if first:
            #     self.a[[-1, 0]] = self.a[[0, -1]]
            #     self.b[0], self.b[-1] = self.b[-1], self.b[0]

            self.b = self.b.tolist()
            row_to_add = -np.modf(self.curr_table[i])[0]
            self.b.append(row_to_add[-1])
            row_to_add = row_to_add[:-1]
            self.a = np.r_[self.a, np.array([row_to_add])]
            col_to_add = np.zeros(shape=(self.a.shape[0], 1))
            col_to_add[-1][0] = 1
            self.a = np.c_[self.a, col_to_add]
            if self.first:
                self.inequalities = np.array(len(self.inequalities)*[RelationSymbols.equals]).tolist()
                self.inequalities.append(RelationSymbols.equals)
            self.inequalities.append(RelationSymbols.equals)
            # print("ineq", len(self.inequalities))
            self.var_num += 1

            # print("+++++++++++++++++++++++++++")
            # print(self.a)
            # print(np.array(self.b).transpose())
            # print("000000000000000000000000000")

            print()
            print("#########################################################")
            print("Current state in Gomory method (current number of cuts: %s)" % num_iter)
            print("#########################################################")
            print("Table:")
            mprint(self.curr_table)
            print("Current solution:")
            print(self.solution)
            self.solution = self.solution[0][:self.solution.shape[1] - (self.solution.shape[1] - self.initial_var_num)]
            print("x0 = %s" % self.solution)
            print()
            if self.is_integer_solution():
                break
            num_iter += 1
            self.first = False
            print("########################################")

        print("######################################")
        print("Solved using Gomory method with %s cuts" % num_iter)
        print("######################################")
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
            print(self.initial_b[i])
        # if self.f_opt_sign == -1:
        #     print("max_f: %s" % self.b[-1])
        # else:
        #     print("min_f: %s" % str(-self.b[-1]))
        print("******************************")
        print("Optimal integer solution is: ")
        print("x0 = %s" % self.solution)
        print("f_opt = %s" % self.calculate_f())
        print("******************************\n\n")

    def calculate_f(self):
        return -np.dot(self.solution, self.initial_c)


def main():
    output_file = "output.txt"
    mode = "a"
    sys.stdout = FileLogger(output_file, mode)
    ans = 'y'
    while ans == 'y':
        np.set_printoptions(suppress=True)
        stop_logging_to_file()

        rows = int(input("Enter number of rows in the matrix A: "))
        columns = int(input("Enter number of columns in the matrix A: "))

        matrix_a = []
        inequalities = []
        print("Enter the %s x %s matrix A one row by one and (in)equality signs when asked: " % (rows, columns))
        for i in range(rows):
            print("Row %d: " % (i + 1))
            matrix_a.append(list(map(float, input().rstrip().split())))
            c = input("Choose (in)equality type (<=, >=, =): ")
            if c == "<=":
                inequalities.append(RelationSymbols.less)
            elif c == ">=":
                inequalities.append(RelationSymbols.greater)
            elif c == "=":
                inequalities.append(RelationSymbols.equals)

        b = []
        print("Enter the b vector of length %s: " % rows)
        b.append(list(map(float, input().rstrip().split())))

        matrix_a = np.array(matrix_a)

        c = []
        print("Enter vector c, which will be optimized: ")
        c.append(list(map(float, input().rstrip().split())))
        opt = input("Optimization (min/max)? ")

        # rows = 2
        # columns = 4
        # matrix_a = np.array([[2, -2, -3, 2], [0, 3, 3, -1]])
        # c = [[-1, -1, 0, 0]]
        # b = [[5, 3]]
        # inequalities = [RelationSymbols.equals, RelationSymbols.equals]
        # opt = "min"

        # rows = 2
        # columns = 2
        # matrix_a = np.array([[-1, 3], [-1, 4]])
        # c = [[1, 1]]
        # b = [[5, 2]]
        # inequalities = [RelationSymbols.greater, RelationSymbols.less]
        # opt = "max"

        # rows = 2
        # columns = 2
        # matrix_a = np.array([[-4, 6], [1, 1]])
        # c = [[1, -2]]
        # b = [[9, 4]]
        # inequalities = [RelationSymbols.less, RelationSymbols.less]
        # opt = "max"

        # rows = 3
        # columns = 2
        # matrix_a = np.array([[2, 11], [1, 1], [4, -5]])
        # c = [[1, 1]]
        # b = [[38, 7, 5]]
        # inequalities = [RelationSymbols.less, RelationSymbols.less, RelationSymbols.less]
        # opt = "max"

        # rows = 3
        # columns = 2
        # matrix_a = np.array([[8, -7], [1, 6], [7, -3]])
        # c = [[1, 1]]
        # b = [[-14, 60, 16]]
        # inequalities = [RelationSymbols.greater, RelationSymbols.less, RelationSymbols.less]
        # opt = "max"

        # rows = 3
        # columns = 2
        # matrix_a = np.array([[5, -3], [9, -8], [2, 1]])
        # c = [[1, 1]]
        # b = [[10, -8, 24]]
        # inequalities = [RelationSymbols.less, RelationSymbols.greater, RelationSymbols.less]
        # opt = "max"

        # rows = 3
        # columns = 2
        # matrix_a = np.array([[2, 11], [1, 1], [4, -5]])
        # c = [[1, 1]]
        # b = [[38, 7, 5]]
        # inequalities = [RelationSymbols.less, RelationSymbols.less, RelationSymbols.less]
        # opt = "max"

        f_opt_sign = 1
        if opt == "max":
            f_opt_sign = -1
            c = np.negative(c)

        problem = Gomory(matrix_a, b[0], c[0], inequalities, columns, f_opt_sign)
        problem.solve_problem()

        # start_logging_to_file(output_file, mode)
        # stop_logging_to_file()

        ans = input("Do you want to continue testing? (y/n) ")
        # start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
