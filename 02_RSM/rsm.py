import numpy as np


class RSM:
    def __init__(self, a, b, c, var_num, f_opt_sign):
        self.a = a
        self.b = b
        self.c = c
        self.initial_var_num = var_num
        self.f_opt_sign = f_opt_sign

    def solve_problem(self):
        x0 = np.array(np.zeros(shape=(1, self.initial_var_num)))
        x0 = np.concatenate((x0, self.b), axis=None)
        p = [i for i in range(self.initial_var_num, self.a.shape[1])]
        q = [i for i in range(self.initial_var_num)]
        self.c = np.c_[self.c, np.zeros(shape=(1, self.a.shape[1]-self.initial_var_num))]

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

            if np.count_nonzero(cn_p < 0) == 0:
                print("*********************************")
                print("Optimal solution is found after %d iterations:\nx0: %s" % (iter_num, str(x0)))
                print("f_opt = %s" % (u.dot(self.b)*self.f_opt_sign))
                break

            j = 0
            for j in range(len(cn_p)):
                if cn_p[j] < 0:
                    break

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

            if np.count_nonzero(y[0] > 0) == 0:
                print("Problem is unlimited, f->-inf.")
                return

            t_mins = []
            for i in range(len(t)):
                if i in p and y[0][i] > 0 and x0[i]*y[0][i] > 0:
                    t_mins.append(x0[i]/y[0][i])
            t_min = min(t_mins)

            prod = t_min*y

            l = 0
            for l in range(len(x0)):
                if x0[l] != 0:
                    if x0[l] - prod[0][l] == 0 and y[0][l] > 0:
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


def main():
    ans = 'y'
    while ans == 'y':
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
        matrix_a = np.c_[matrix_a, np.eye(matrix_a.shape[0])]

        c = []
        print("Enter vector c, which will be optimized: ")
        c.append(list(map(float, input().rstrip().split())))
        opt = input("Optimization (min/max)? ")
        f_opt_sign = 1
        if opt == "max":
            f_opt_sign = -1
            c = np.negative(c)

        problem = RSM(matrix_a, b[0], c, columns, f_opt_sign)
        problem.solve_problem()

        ans = input("Do you want to continue testing? (y/n) ")


if __name__ == '__main__':
    main()
