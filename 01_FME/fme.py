# Fourier-Motzkin Elimination
import numpy as np
import copy


class InequalityInterval:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __str__(self):
        return "{" + str(self.a) + ", " + str(self.b) + "}"


class SystemOfInequalities:
    def __init__(self, matrix_a):
        self.initial_matrix_a = matrix_a
        self.matrix_a = matrix_a
        self.var_num = matrix_a.shape[1]-1

    def set_matrix_a(self, new_matrix):
        self.matrix_a = new_matrix

    def set_initial_matrix_a(self, new_matrix):
        self.initial_matrix_a = new_matrix

    def set_var_num(self, new_var_num):
        self.var_num = new_var_num

    # Checking if the system contains inaccurate inequality
    @staticmethod
    def bad_system(matrix):
        for i in range(matrix.shape[0]):
            if np.count_nonzero(matrix[i, :-1] == 0) == matrix.shape[1] - 1:
                if matrix[i][-1] > 0:
                    return True
        return False

    @staticmethod
    def change_b_sign(inequalities):
        inequalities[:, -1] = -inequalities[:, -1]
        return inequalities

    def solve_system_intervals(self):
        intervals = []
        matrix_a_curr = copy.deepcopy(self.matrix_a)

        # Loop for each variable x_1, x_2, ..., x_n
        for j in range(self.var_num):
            # Get cardinality of sets I, J, K
            i_plus_inequalities = np.zeros(shape=(np.count_nonzero(matrix_a_curr[:, j] > 0), self.var_num+1))
            i_minus_inequalities = np.zeros(shape=(np.count_nonzero(matrix_a_curr[:, j] < 0), self.var_num+1))
            i_zero_inequalities = np.zeros(shape=(np.count_nonzero(matrix_a_curr[:, j] == 0), self.var_num+1))

            row_plus = 0
            row_minus = 0
            row_zero = 0

            # Loop for each row in the current system and determine interval for x_1, x_2, ..., x_n
            for i in range(matrix_a_curr.shape[0]):

                # Getting I, J, K inequalities from the current system
                if matrix_a_curr[i][j] != 0:
                    if matrix_a_curr[i][j] > 0:
                        i_plus_inequalities[row_plus] = np.divide(np.negative(matrix_a_curr[i]), matrix_a_curr[i][j])
                        i_plus_inequalities[row_plus][j] = 0
                        row_plus += 1
                    if matrix_a_curr[i][j] < 0:
                        i_minus_inequalities[row_minus] = np.divide(np.negative(matrix_a_curr[i]), matrix_a_curr[i][j])
                        i_minus_inequalities[row_minus][j] = 0
                        row_minus += 1
                else:
                    i_zero_inequalities[row_zero] = (matrix_a_curr[i])
                    row_zero += 1

            # Determining the dimension of next system using the information about sets I, J, K
            next_system_dimension = i_plus_inequalities.shape[0] * i_minus_inequalities.shape[0] + \
                i_zero_inequalities.shape[0]
            matrix_a_temp = np.zeros(shape=(next_system_dimension, self.var_num+1))
            # Calculating all possible values for left and the right edges of the interval
            if next_system_dimension > 0:
                row_num = 0
                for k in range(i_plus_inequalities.shape[0]):
                    for l in range(i_minus_inequalities.shape[0]):
                        matrix_a_temp[row_num] = i_minus_inequalities[l]-i_plus_inequalities[k]
                        row_num += 1
                for k in range(i_zero_inequalities.shape[0]):
                    matrix_a_temp[row_num] = i_zero_inequalities[k]
                    row_num += 1
                # If there are no elements for some edge of the interval, set that edge to appropriate inf
                left_interval = SystemOfInequalities.change_b_sign(i_plus_inequalities)
                right_interval = SystemOfInequalities.change_b_sign(i_minus_inequalities)
                if len(left_interval) == 0:
                    left_interval = np.array([[float('-inf')]])
                if len(right_interval) == 0:
                    right_interval = np.array([[float('inf')]])
                intervals.append(InequalityInterval(left_interval, right_interval))
            else:
                if len(i_plus_inequalities) + len(i_minus_inequalities) == 0:
                    intervals.append(InequalityInterval(np.array([[float('-inf')]]), np.array([[float('inf')]])))
                elif len(i_plus_inequalities) == 0:
                    intervals.append(
                        InequalityInterval(np.array([[float('-inf')]]),
                                           np.array(SystemOfInequalities.change_b_sign(i_minus_inequalities))))
                elif len(i_minus_inequalities) == 0:
                    intervals.append(
                        InequalityInterval(np.array(SystemOfInequalities.change_b_sign(i_plus_inequalities)),
                                           np.array([[float('inf')]])))

            if SystemOfInequalities.bad_system(matrix_a_temp):
                print("System does not have a solution!")
                return []

            matrix_a_curr = matrix_a_temp

        return intervals

    def import_values(self, positions, values, matrix_state):
        """
        Importing variable values into the system of inequalities
        :param positions: Vector of the positions of the variable to be changed
        :param values: Vector of the values of the variable to be changed
        :param matrix_state: Initial/Current
        """
        curr_matrix = []
        if matrix_state == "cur":
            curr_matrix = self.matrix_a
        elif matrix_state == "initial":
            curr_matrix = self.initial_matrix_a
        for i in range(self.matrix_a.shape[0]):
            for j in positions:
                curr_matrix[i][-1] = curr_matrix[i][-1]-curr_matrix[i][j]*values[j]
                curr_matrix[i][j] = 0
        self.matrix_a = curr_matrix


def main():
    ans = 'y'
    while ans == 'y':
        fm_format = int(input("System of inequalities format? (1: Ax <= b; 2: Ax >=b) "))
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
        matrix_a = np.c_[matrix_a, np.transpose(np.array(b))]

        if fm_format == 1:
            matrix_a = np.negative(matrix_a)

        system = SystemOfInequalities(matrix_a)

        test_type = int(input("Enter the type of test:\n"
                              "(1) Determine variable intervals\n"
                              "(2) Determine whether the point is in the interval or not\n"
                              "(3) Solve LP problem using Fourier-Motzkin elimination\n"))
        if test_type == 1:
            var_values = []
            positions = []
            # Determining interval for first variable from the left in the system and asking the user
            # to choose a value from the interval. Import chosen value and repeat the procedure for
            # each variable in the system
            for i in range(matrix_a.shape[1]-1):
                if i == 0:
                    if len(system.solve_system_intervals()) == 0:
                        break
                # Using elementary transformations of matrices, put the current variable
                # to the last column of each row
                permutation = [j for j in range(matrix_a.shape[1]-1)]
                permutation.append(matrix_a.shape[1]-1)
                temp = permutation[-2]
                permutation[-2] = permutation[i]
                permutation[i] = temp
                m = np.argsort(permutation)
                system.set_matrix_a(system.matrix_a[:, m])
                intervals = system.solve_system_intervals()
                if len(intervals) == 0:
                    break
                print("{ %s, %s }" % (str(max(intervals[-1].a[:, -1:])[0]), str(min(intervals[-1].b[:, -1:])[0])))
                if i == matrix_a.shape[1]-2:
                    break
                var_values.append(float(input("Choose value for the variable from the range above: ")))
                positions.append(i)
                system.import_values(positions, var_values, "initial")
        elif test_type == 2:
            point = []
            print("Enter the point (x y z ...): ")
            point.append(list(map(float, input().rstrip().split())))
            system.import_values([i for i in range(matrix_a.shape[0]-1)], point[0], "initial")
            if SystemOfInequalities.bad_system(system.matrix_a):
                print("Point (%s, %s, %s) IS NOT in the interval" % (point[0][0], point[0][1], point[0][2]))
            else:
                print("Point (%s, %s, %s) IS in the interval" % (point[0][0], point[0][1], point[0][2]))
        elif test_type == 3:
            c = []
            print("Enter vector c, which will be optimized: ")
            c.append(list(map(float, input().rstrip().split())))
            optimization_problem = input("Optimization (min/max)? ")
            matrix_a = np.r_[matrix_a, np.eye(matrix_a.shape[1]-1, matrix_a.shape[1])]
            matrix_a = np.c_[matrix_a[:, 0:-1], np.zeros(shape=(matrix_a.shape[0], 1)), matrix_a[:, -1:]]
            if optimization_problem == "min":
                c = np.c_[c, np.ones(shape=(1, 1)), np.zeros(shape=(1, 1))]
                c = np.negative(c)
                c[0][-2] = -c[0][-2]
                matrix_a = np.r_[c, matrix_a]
                system.set_matrix_a(matrix_a)
                system.set_initial_matrix_a(matrix_a)
                system.set_var_num(matrix_a.shape[1]-1)
                intervals = system.solve_system_intervals()
                if len(intervals) == 0:
                    break
                print(max(intervals[-1].a[:, -1:])[0])
            elif optimization_problem == "max":
                c = np.c_[c, np.ones(shape=(1, 1)), np.zeros(shape=(1, 1))]
                c[0][-2] = -c[0][-2]
                matrix_a = np.r_[c, matrix_a]
                system.set_matrix_a(matrix_a)
                system.set_initial_matrix_a(matrix_a)
                system.set_var_num(matrix_a.shape[1]-1)
                intervals = system.solve_system_intervals()
                if len(intervals) == 0:
                    break
                print(min(intervals[-1].b[:, -1:])[0])

        ans = input("Do you want to continue testing? (y/n) ")


if __name__ == "__main__":
    main()
