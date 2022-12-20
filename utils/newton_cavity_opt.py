import numpy as np


def evaluate_func(pn):
    return np.array([5, 5])

if __name__ == '__main__':

    # xn+1 = xn - J(xn)f(xn)
    pn = np.array([1, 2, 3])  # starting vector, geometric variables
    # calculate geometry

    # get objective function values
    fn = np.array([1, 2])

    # small offset from point
    delta_p = np.array([0.1, 0.2, 0.2])

    # compute Jacobian
    J = np.zeros((len(fn), len(pn)))
    pn_1 = np.zeros((len(pn)))
    for j in range(len(pn)):
        # evaluate objective functions
        for i in range(len(fn)):
            # evaluate funciton with variables pn[i] = pn[i] + delta_p[i]
            pn_1[i] = pn[i] + delta_p[i]
            fn_1 = evaluate_func(pn)

            # calculate jacobian
            J[i][j] = (fn_1[i]-fn[i])/(pn_1[j]-pn[j])

            # solve linear equation
            #
    # pn_1 = np.linalg.solve(J, J@pn - fn)

    print(J)