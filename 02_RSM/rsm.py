# Revised Simplex Method
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


class RSM:
    def __init__(self, a, b, c, var_num, f_opt_sign, system_form):
        self.a = a
        self.initial_a = copy.deepcopy(a)
        self.b = b
        self.initial_b = copy.deepcopy(b)
        self.c = c
        self.initial_var_num = var_num
        self.f_opt_sign = f_opt_sign
        self.system_form = system_form

    def find_basic_feasible_solution(self):
        x0 = np.array(np.zeros(shape=(1, self.initial_var_num)))[0]
        # If it's in canonical form, Ax=b, find identity submatrix and assign them values of vector b respectively
        if self.system_form == 1:
            eye_matrices = 0
            for j in range(self.a.shape[1]):
                for e in range(self.a.shape[0]):
                    if np.all(self.a[:, j] == np.eye(self.a.shape[0])[:, e]):
                        x0[j] = self.b[e]
                        eye_matrices += 1
            if eye_matrices != self.a.shape[0]:
                return np.array([-1])

        # When in the form Ax<=b or Ax>=b, BFS is easily found by just adding vector b columns
        if self.system_form == 2 or self.system_form == 3:
            x0 = np.concatenate((x0, self.b), axis=None)
        return x0

    def solve_problem(self):
        # Ax<=b, x>=0
        if self.system_form == 2:
            # Adding positive slack variables
            self.a = np.c_[self.a, np.eye(self.a.shape[0])]
        # Ax>=b, x>=0
        elif self.system_form == 3:
            # Adding negative slack variables and multiplying the system by -1
            self.a = np.c_[np.negative(self.a), np.eye(self.a.shape[0])]
            self.b = np.negative(self.b)

        x0 = self.find_basic_feasible_solution()

        # Check if the initial solution is feasible
        if np.count_nonzero(x0 < 0) > 0:
            print("System can't be solved using Revised Simplex Method.")
            return
        # Setting initial P and Q
        p = [i for i in range(x0.shape[0]) if x0[i] > 0]
        q = [i for i in range(x0.shape[0]) if x0[i] == 0]
        self.c = np.c_[self.c, np.zeros(shape=(1, self.a.shape[1]-self.initial_var_num))]

        print("Input:\nA = \n%s,\n b = %s," % (self.initial_a, self.initial_b))
        if self.system_form == 1:
            print("Constraints are given in the form Ax=b")
        if self.system_form == 2:
            print("Constraints are given in the form Ax<=b")
        if self.system_form == 3:
            print("Constraints are given in the form Ax>=b")
        if self.f_opt_sign == -1:
            print("(max) c = %s" % np.negative(self.c[0, :-self.initial_var_num + 1]))
        else:
            print("(min) c = %s" % self.c[0, :-self.initial_var_num + 1])

        print("*********************************")
        print("Output:")
        print("*********************************")
        print("Initial x0 is: %s" % str(x0))
        print("Initial P is: %s" % str(p))
        print("Initial Q is: %s" % str(q))

        iter_num = 1
        while True:
            cb = np.array(self.c[0, p])
            cn = np.array(self.c[0, q])
            B = self.a[:, p]
            N = self.a[:, q]
            u = np.matmul(cb, np.linalg.inv(B))
            cn_p = cn - np.matmul(u, N)

            # Optimal solution is found if there are no negative rj values
            if np.count_nonzero(cn_p < 0) == 0:
                print("*********************************")
                print("Optimal solution is found after %d iterations:\nx0: %s" % (iter_num, str(x0)))
                print("f_opt = %s" % (u.dot(self.b)*self.f_opt_sign))
                break

            # Finding j for rj < 0
            j = 0
            for j in range(len(cn_p)):
                if cn_p[j] < 0:
                    break

            # Determining vectors y and parametric t vector
            y = np.zeros(shape=(1, self.a.shape[1]))
            y[:, p] = np.linalg.solve(B, self.a[:, q[j]])
            t = np.zeros(shape=(1, self.a.shape[1]))[0]

            for elem in range(len(t)):
                if elem == q[j]:
                    t[elem] = 1
                elif elem in p:
                    t[elem] = -y[0][elem]
                elif elem in q:
                    t[elem] = 0

            # If all yp <= 0, function is unlimited
            if np.count_nonzero(y[0] > 0) == 0:
                print("Problem is unlimited, f->-inf.")
                return

            # Calculating optimal t, t_kappa
            t_mins = []
            for i in range(len(t)):
                if i in p and y[0][i] > 0 and x0[i]*y[0][i] > 0:
                    t_mins.append(x0[i]/y[0][i])
            t_min = min(t_mins)

            prod = t_min*y

            # Choose l so that the xl - tkapa*yl > 0
            l = 0
            for l in range(len(x0)):
                if x0[l] != 0:
                    epsilon = 1e-10
                    if abs(x0[l] - prod[0][l]) <= epsilon and y[0][l] > 0:
                        break
            x0 = x0+t*t_min
            x0[l] = 0

            p.remove(l)
            p.append(q[j])
            q.remove(q[j])
            q.append(l)
            p.sort()
            q.sort()

            print("***********ITERATION %s***********" % iter_num)
            iter_num += 1
            print("Current x0 is: %s" % str(x0))
            print("Current P is: %s" % str(p))
            print("Current Q is: %s" % str(q))
            print("f = %s" % (u.dot(self.b)*self.f_opt_sign))

        print("\n\n")


def main():
    output_file = "output.txt"
    mode = "a"
    sys.stdout = FileLogger(output_file, mode)
    ans = 'y'
    while ans == 'y':
        stop_logging_to_file()
        system_form = int(input(
            "System form (1: Ax=b; 2: Ax<=b; 3: Ax>=b; 4: Combination of <= and >= inequalities)? "))
        if system_form == 4:
            print("System can't be solved using Revised Simplex Method.")
        else:
            rows = int(input("Enter number of rows in the matrix A: "))
            columns = int(input("Enter number of columns in the matrix A: "))
            matrix_a = []
            print("Enter the %s x %s matrix A: " % (rows, columns))
            for i in range(rows):
                matrix_a.append(list(map(float, input().rstrip().split())))

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

            problem = RSM(matrix_a, b[0], c, columns, f_opt_sign, system_form)
            start_logging_to_file(output_file, mode)
            problem.solve_problem()
            stop_logging_to_file()

        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
