import numpy as np
from tabulate import tabulate
def printMatrixE(a):
   print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in a]))

def spectral(matrix):
    matrix = np.matrix(matrix)
    return np.linalg.norm(matrix)*np.linalg.norm(matrix.I)
def ortega(matrix):
    num_rows, num_cols = np.shape(matrix)
    result = 1
    sum_of_col = 0
    for n in range(num_rows):
        for m in range(num_cols):
            sum_of_col += matrix[n,m] * matrix[n,m]
        sum_of_col = np.sqrt(sum_of_col)
        result *= sum_of_col
        sum_of_col = 0

    return result/np.linalg.det(matrix)
def angle(matrix):
    matrix_lengts = []
    inv_matrix_lenghts = []
    inv_matrix = np.linalg.inv(matrix)
    num_rows, num_cols = np.shape(matrix)
    num_rows_inv, num_cols_inv = np.shape(inv_matrix)
    cond_list = []
    matrix_rows = []
    inv_matrix_columns = []
    inv_matrix_T = inv_matrix.T
    sum_of_col = 0
    for n in range(num_rows):
        matrix_rows.extend([matrix[:][n]])
    for m in range(num_cols_inv):
        inv_matrix_columns.extend([inv_matrix_T[m][:]])
    for z in range(num_rows):
        for y in range(num_cols):
            sum_of_col += matrix_rows[z][y] * matrix_rows[z][y]
        sum_of_col = np.sqrt(sum_of_col)
        matrix_lengts.append(sum_of_col)
        sum_of_col = 0
    for z in range(num_rows):
        for y in range(num_cols):
            sum_of_col += inv_matrix_columns[z][y] * inv_matrix_columns[z][y]
        sum_of_col = np.sqrt(sum_of_col)
        inv_matrix_lenghts.append(sum_of_col)
        sum_of_col = 0
    for k in range(num_rows):
        cond_list.append(matrix_lengts[k]*inv_matrix_lenghts[k])
    return max(cond_list)


def main():
    A = np.array([[1,1/2,1/3],
         [1/2,1/3,1/4],
         [1/3,1/4,1/5]])
    x = np.array([[2.17],
        [-6.14],
        [0.15]])

    B = np.dot(A,x)
    print('Матрица A:')
    printMatrixE(A)
    print()
    print(tabulate([[0,A,B,np.linalg.solve(A,B),0,spectral(A),ortega(A),angle(A)], ['10^(-2)',A+0.01,B+0.01,np.linalg.solve(A+0.01,B+0.01),abs(x-np.linalg.solve(A+0.01, B+0.01)),spectral(A+0.01),ortega(A+0.01),angle(A+0.01)],['10^(-5)',A+10**(-5),B+10**(-5),np.linalg.solve(A+10**(-5),B+10**(-5)),abs(x-np.linalg.solve(A+10**(-5),B+10**(-5))),spectral(A+10**(-5)),ortega(A+10**(-5)),angle(A+10**(-5))],['10^(-8)',A+10**(-8),B+10**(-8),np.linalg.solve(A+10**(-8),B+10**(-8)),abs(x-np.linalg.solve(A+10**(-8),B+10**(-8))),spectral(A+10**(-8)),ortega(A+10**(-8)),angle(A+10**(-8))]], headers=['Погрешеность','A','B','X','|x - x~|','Спектральный', 'Ортега','Угловой'], tablefmt='orgtbl'))


if __name__ == "__main__":
    main()

