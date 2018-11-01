# Revised Simplex Method with Eta factorization
import numpy as np
import sys


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


class RSM_EF:
    def __init__(self, a, b, c, var_num, f_opt_sign):
        self.a = a
        self.b = b
        self.c = c
        self.initial_var_num = var_num
        self.f_opt_sign = f_opt_sign

    def solve_problem(self):
        x0 = np.array(np.zeros(shape=(1, self.initial_var_num)))
        x0 = np.concatenate((x0, self.b), axis=None)
        # Check if the initial solution is feasible
        # Ukoliko se u pocetnom bazisnom resenju x0 nalazi broj, koji je manji od nule, nije ispunjen uslov
        # pozitivnosti resenja, pa vracamo odgovarajucu poruku i izlazimo iz algoritma.
        if np.count_nonzero(x0 < 0) > 0:
            print("System can't be solved using Revised Simplex Method with Eta .")
            return
        # Odredjujemo inicijalne P i Q
        p = [i for i in range(self.initial_var_num, self.a.shape[1])]
        q = [i for i in range(self.initial_var_num)]
        self.c = np.c_[self.c, np.zeros(shape=(1, self.a.shape[1]-self.initial_var_num))]
        B = self.a[:, p]
        E = np.eye(len(p))
        eta_matrices = []

        # Stampamo i input radi boljeg formatiranja output fajla
        print("Input:\nA = \n%s,\n b = %s," % (self.a[:, :-self.initial_var_num], self.b))
        if self.f_opt_sign == -1:
            print("(max) c = %s" % self.c[0])
        else:
            print("(min) c = %s" % self.c[0])
        print("*********************************")
        print("Initial x0 is: %s" % str(x0))
        print("Initial P is: %s" % str(p))
        print("Initial Q is: %s" % str(q))

        iter_num = 1
        while True:
            cb = np.array(self.c[0, p])
            cn = np.array(self.c[0, q])
            # Izracunavanje nove matrice B
            B = np.matmul(B, E)
            v = []
            # Odredjivanje vektora v, w, ... iterativno, korisnjem informacija o prethodnim eta matricama
            for i in range(len(eta_matrices) - 1):
                if i == 0:
                    v.append(np.matmul(cb, np.linalg.inv(eta_matrices[-1])))
                else:
                    v.append(np.matmul(v[-i], np.linalg.inv(eta_matrices[-1-i])))
            # Odredjivanje vektora u
            u = np.matmul(cb, np.linalg.inv(B))

            N = self.a[:, q]
            cn_p = cn - np.matmul(u, N)

            # Optimal solution is found if there are no negative rj values
            # Ukoliko su svi elementi rj > 0, optimalno resenje je pronadjeno
            if np.count_nonzero(cn_p < 0) == 0:
                print("*********************************")
                print("Optimal solution is found after %d iterations:\nx0: %s" % (iter_num, str(x0)))
                for w in range(len(v)):
                    print("v%d = %s" % (w, v[w]))
                print("u = %s" % u)
                print("f_opt = %s" % (u.dot(self.b)*self.f_opt_sign))
                break

            # Pronalazimo odgovarajuce rj < 0
            j = 0
            for j in range(len(cn_p)):
                if cn_p[j] < 0:
                    break

            # Odredjujemo y, pa pravimo x(t) - cuvamo u vektoru t
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

            # Ukoliko su svi yp, p pripada P, nepozitivni, funkcija cilja neograniceno opada ka -inf
            if np.count_nonzero(y[0] > 0) == 0:
                print("Problem is unlimited, f->-inf.")
                return

            # Odredjujemo optimalno t, t kapa
            t_mins = []
            for i in range(len(t)):
                if i in p and y[0][i] > 0 and x0[i]*y[0][i] > 0:
                    t_mins.append(x0[i]/y[0][i])
            t_min = min(t_mins)

            prod = t_min*y

            # Biramo l takvo da je xl - tkapa*yl > 0
            l = 0
            for l in range(len(x0)):
                if x0[l] != 0:
                    if x0[l] - prod[0][l] == 0 and y[0][l] > 0:
                        break
            x0 = x0+t*t_min
            x0[l] = 0

            print("***********ITERATION %s***********" % iter_num)
            print("Current x0 is: %s" % str(x0))
            print("Current P is: %s" % str(p))
            print("Current Q is: %s" % str(q))
            print("f = %s" % (u.dot(self.b)*self.f_opt_sign))
            # Kreiranje nove eta matrice, koja ce se koristiti u narednoj iteraciji i smestanje trenutne u niz
            # eta matrica, koji cemo koristiti za iterativno racunanje vektora v, w, ...
            E = np.eye(len(p))
            E[:, p.index(l)] = y[0, p]
            eta_matrices.append(E)
            # Azuriramo skupove P i Q, koji sadrze bazisne, odnosno nebazisne kolone
            p[p.index(l)] = q[j]
            q.remove(q[j])
            q.append(l)
            print("E%d = " % iter_num)
            print(E)
            for w in range(len(v)):
                print("v%d = %s" % (w, v[w]))
            print("u = %s" % u)
            iter_num += 1
        print("\n\n")


def main():
    output_file = "output.txt"
    mode = "a"
    sys.stdout = FileLogger(output_file, mode)
    ans = 'y'
    while ans == 'y':
        stop_logging_to_file()
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

        problem = RSM_EF(matrix_a, b[0], c, columns, f_opt_sign)
        start_logging_to_file(output_file, mode)
        problem.solve_problem()
        stop_logging_to_file()
        ans = input("Do you want to continue testing? (y/n) ")
        start_logging_to_file(output_file, mode)


if __name__ == '__main__':
    main()
