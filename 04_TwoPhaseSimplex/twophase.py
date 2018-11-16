# The Two Phase Simplex Method
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


class TPSM:
    def __init__(self, a, b, c, inequalities, var_num, f_opt_sign):
        self.a = a
        self.initial_a = copy.deepcopy(a)
        self.b = b
        self.initial_b = copy.deepcopy(b)
        self.c = c
        self.initial_var_num = var_num
        self.f_opt_sign = f_opt_sign
        self.inequalities = inequalities
        self.artificial_var_num = 0
        self.slack_var_num = 0

    def solve_problem(self):
        if not self.phase_one_simplex():
            print("The original problem is infeasible!\n\n")
            return
        self.phase_two_simplex()

    def tableau_method(self):
        epsilon = 1e-10
        # Izvrsava se sve dok postoji rj < 0
        while np.count_nonzero(self.a[-1] < 0) > 0:
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

            for rj in range(self.a.shape[1]):
                if abs(self.a[-1][rj]) < epsilon:
                    self.a[-1][rj] = 0

            mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=40))
            print()

    def phase_one_simplex(self):
        print("*****************************************")
        print("*************** PHASE ONE ***************")
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
        self.a[-1, :] = row_to_add
        self.b.append(-b_to_add)
        mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
        print()

        self.tableau_method()

        if np.around(self.b[-1], decimals=5) != 0:
            print(self.b[-1])
            return False
        else:
            return True

    def phase_two_simplex(self):
        print("*****************************************")
        print("*************** PHASE TWO ***************")
        print("*****************************************")

        self.a = np.delete(self.a, self.a.shape[0] - 1, 0)

        # Remove nonbasic artificial columns
        # Izbacivanje vestackih nebazisnih kolona
        for j in range(self.a.shape[1] - self.artificial_var_num, self.a.shape[1]):
            for e in range(self.a.shape[0]):
                if np.all(self.a[:, j] != np.eye(self.a.shape[0])[:, e]):
                    self.a[:, j] = np.zeros(shape=(self.a.shape[0], 1))[0]

        self.a = self.a[:, ~np.all(self.a == 0, axis=0)]

        # Searching for an identity matrix before adding row for the function f
        # Potraga za jedinicnom matricom, pre dodavanja poslednje vrste, funkcije f
        eye_matrices_pos = []
        for j in range(self.a.shape[1]):
            for e in range(self.a.shape[0]):
                if np.all(self.a[:, j] == np.eye(self.a.shape[0])[:, e]):
                    eye_matrices_pos.append([(np.where(self.a[:, j] == 1)[0][0]), j])

        # Add row for initial f
        # Dodavanje reda funkcije ciji optimum trazimo
        self.a = np.r_[self.a, np.c_[[self.c], np.zeros(shape=(1, self.a.shape[1] - len(self.c)))]]

        mprint(np.around(np.c_[self.a, np.array(self.b).transpose()], decimals=4))
        print()

        # Cancel columns of f on positions which correspond to basic columns
        # Ponistavanje kolona od f na pozicijama koje odgovaraju bazisnim kolonama
        for eye_pos in eye_matrices_pos:
            self.b[-1] -= (self.b[eye_pos[0]] * self.a[-1][eye_pos[1]])
            self.a[-1] -= (self.a[eye_pos[0]] * self.a[-1][eye_pos[1]])

        # Izbaci sve vestacke
        while self.a.shape[1] > self.initial_var_num + self.slack_var_num:
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

        print("******************************")
        print("Input:")
        if self.f_opt_sign == -1:
            print("(max) c = %s" % np.negative(self.c))
        else:
            print("(min) c = %s" % self.c)
        print("subject to: ")
        for i in range(self.initial_a.shape[0]):
            mprint(np.array([self.initial_a[i]]), new_line=False)
            if self.inequalities[i] == RelationSymbols.greater:
                print(" >= ", end="")
            elif self.inequalities[i] == RelationSymbols.less:
                print(" <= ", end="")
            elif self.inequalities[i] == RelationSymbols.equals:
                print(" = ", end="")
            print(self.initial_b[i])
        print("******************************")
        print("Optimal solution is: ")
        mprint(x0)
        if self.f_opt_sign == -1:
            print("max_f: %s" % self.b[-1])
        else:
            print("min_f: %s" % str(-self.b[-1]))
        print("******************************")
        print("\n\n")


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
        f_opt_sign = 1
        if opt == "max":
            f_opt_sign = -1
            c = np.negative(c)

        rows = 3
        columns = 4
        matrix_a = np.array([[1, -1, 3, -1],
                             [4, 0, -1, 3],
                             [1, 3, 0, -3]])
        b = [[9, 1, 3]]
        c = [[1, 1, 0, 5]]
        inequalities = [RelationSymbols.greater, RelationSymbols.greater, RelationSymbols.greater]
        f_opt_sign = 1

        problem = TPSM(matrix_a, b[0], c[0], inequalities, columns, f_opt_sign)
        start_logging_to_file(output_file, mode)
        problem.solve_problem()
        stop_logging_to_file()

        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
