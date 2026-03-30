""" This is taken from cobrapy and modified to work with Pyomo. Cobrapy is licensed under GNU GPL v2.0.
## Licensing

This package, MyPackage, includes code from OriginalProject, which is licensed under the GNU General Public License v2. The following files are derived from OriginalProject:

- `src/modified_function1.py`: Modified to integrate with MyPackage's framework.
- `src/modified_function2.py`: Adapted for additional functionality in MyPackage.

MyPackage as a whole is licensed under the GNU General Public License v2. See the `COPYING` file for more details.

add copying file

"""

import scipy
import pyomo.environ as pe
from typing import TYPE_CHECKING, Dict, Optional, Union

import time
import numpy as np
from optlang.symbolics import Zero

import cobra
from cobra.core import get_solution
from cobra.util import create_stoichiometric_matrix, nullspace
from cobra.flux_analysis.helpers import normalize_cutoff

if TYPE_CHECKING:
    from cobra import Model, Reaction, Solution
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, hstack, vstack
from scipy.sparse.linalg import svds, splu
from scipy.linalg import lu
from scipy.sparse import csc_matrix
from copy import deepcopy

import scipy
import numpy as np
from scipy.sparse.linalg import svds
from scipy.linalg import lu


def sparseNull(S, tol=1e-9):
    # S = csr_matrix(S)
    SpLeft, SpRight = spspaces(S, 2, tol)
    N = SpRight[0][:, SpRight[2]]
    rankS = len(SpRight[1])
    N[np.abs(N) < tol] = 0
    return N, rankS


def spspaces(A, opt, tol):
    L, U, Q = luq(A, 0, tol)
    if opt == 1:
        SpLeft = spspaces_side(L, U, tol)
        SpRight = {}
    elif opt == 2:
        SpLeft = {}
        SpRight = spspaces_side(Q, U.T, tol)
    else:
        SpLeft = spspaces_side(L, U, tol)
        SpRight = spspaces_side(Q.T, U.T, tol)
    return SpLeft, SpRight


def spspaces_side(L, U, tol):
    if L.shape[0] == 0 or L.shape[1] == 0:
        QQ = L
    else:
        QQ = np.linalg.inv(L)
    S = np.max(np.abs(U), axis=1).T
    S_dense = S.toarray()
    I = np.where(S_dense > tol)[1]
    # check if S is empty:
    if S.shape[0] == 0 or S.shape[1] == 0:
        J = np.arange(S.shape[1])
    else:
        J = np.where(S_dense <= tol)[1]
    return QQ, I, J, L


def luq(A, do_pivot, tol):

    if type(A) is not csc_matrix:
        A = csc_matrix(A)
    (n, m) = A.shape

    if A.shape[0] == 0:
        L = np.eye(n)
        U = A
        Q = np.eye(m)
        return L, U, Q
    if A.shape[1] == 0:
        L = np.eye(n)
        U = A
        Q = np.eye(m)
        return L, U, Q

    if do_pivot:
        ### Could not find a way to reproduce Q from matlab's LU funciton
        # PL, U = lu(A.toarray(), permute_l=True)
        # P, L, U = lu(A.toarray())
        # _, n_ = A.shape
        # Q = np.eye(n_)[:, np.argsort(np.argmax(U != 0, axis=0))].T
        pass

    else:
        P, L, U = lu(A.toarray())
        Q = np.eye(m)
    small_p = n - L.shape[1]
    print(small_p)
    if small_p == 0:
        # LL = []
        L = L @ P
        U = U
    else:
        L = np.vstack((L.T @ P, P[n - small_p : n, :]))
        U = np.vstack((U, np.zeros((small_p, m))))

    if U.shape[0] == 1 or U.shape[1] == 1:
        S = U[0, 0]
    else:
        S = np.diag(U)

    I = np.where(np.abs(S) > tol)
    Jl = np.setdiff1d(np.arange(n), I)
    Jq = np.setdiff1d(np.arange(m), I)

    I_ori = deepcopy(I)
    I = np.ravel(I)
    Ubar1 = U[np.ix_(I, I)]
    Ubar2 = U[np.ix_(Jl, Jq)]
    Qbar1 = Q[I, :]
    Lbar1 = L[:, I]

    if len(I) > 0:
        Utmp = U[I, :][:, Jq]
        X = np.linalg.solve(Ubar1.T, U[Jl, :][:, I].T)
        Ubar2 = Ubar2 - X.T @ Utmp
        Lbar1 = Lbar1 + L[:, Jl] @ X.T

        X = np.linalg.solve(Ubar1, Utmp)
        Qbar1 = Qbar1 + X @ Q[Jq, :]
        Utmp = []
        X = []

    I2 = np.where(np.max(np.abs(Ubar2), axis=1) > tol)
    I5 = np.where(np.max(np.abs(Ubar2), axis=0) > tol)

    I3 = Jl[I2]
    I4 = Jq[I5]
    Jq = np.delete(Jq, I5)  # Remove elements at indices I5 from Jq1
    Jl = np.delete(Jl, I2)  # Remove elements at indices I2 from Jl
    U = []

    I2_ravel = np.ravel(I2)
    I5_ravel = np.ravel(I5)

    A = Ubar2[np.ix_(I2_ravel, I5_ravel)]
    L1, U1, Q1 = luq(A, do_pivot, tol)

    if len(L1) == 0:
        Lbar2 = []
        L = np.hstack((Lbar1, L[:, Jl]))
    else:
        Lbar2 = L[:, I3] @ L1
        L = np.hstack((Lbar1, Lbar2, L[:, Jl]))

    if len(Q1) == 0:
        Qbar2 = []
        Q3 = np.vstack((Qbar1, Q[Jq, :]))
    else:
        Qbar2 = Q1 @ Q[I4, :]
        Q = np.vstack((Qbar1, Qbar2, Q[Jq, :]))

    n1 = len(I)
    n2 = len(I3)
    m2 = len(I4)

    sparse1 = csr_matrix((n1, m - n1))
    sparse2 = csr_matrix((n2, n1))
    sparse3 = csr_matrix((n2, m - n1 - m2))
    sparse4 = csr_matrix((n - n1 - n2, m))

    U = scipy.sparse.vstack(
        (
            scipy.sparse.hstack((Ubar1, sparse1)),
            scipy.sparse.hstack((sparse2, U1, sparse3)),
            sparse4,
        )
    )

    return L, U, Q


if __name__ == "__main__":
    toy_stoichemetric_matrix = np.array(
        [
            # 1   2   3   4   5   6   7   8  9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26   27
            [
                1,
                -1,
                0,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # A
            [
                -1,
                1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # B
            [
                0,
                0,
                1,
                -1,
                -1,
                0,
                0,
                0,
                0,
                0,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # C
            [
                0,
                0,
                0,
                0,
                1,
                -1,
                0,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # D
            [
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # E
            [
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                -1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # F
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                0,
                0,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -1,
                0,
            ],  # G
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                0,
                0,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # H
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # I
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # J
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -1,
                0,
                0,
                0,
                1,
                0,
                0,
            ],  # K
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -1,
                1,
                -1,
                0,
                0,
                0,
                0,
                0,
            ],  # L
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                0,
                0,
                0,
                0,
            ],  # M
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -1,
                -1,
                0,
                1,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # N
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -1,
            ],  # O
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                0,
                0,
                0,
            ],  # P
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                -1,
                0,
                0,
                0,
                0,
                1,
                0,
                0,
                0,
            ],  # Q
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],  # R
        ]
    )

    print(len(toy_stoichemetric_matrix[0]))
    # find columns that are not sum 0
    index_of_non_zero_columns = np.where(np.sum(toy_stoichemetric_matrix, axis=0) != 0)[
        0
    ]
    print(index_of_non_zero_columns)
    # remove columns that are sum 0 keep the rest
    toy_stoichemetric_matrix_without_zero_columns = np.delete(
        toy_stoichemetric_matrix, index_of_non_zero_columns, axis=1
    )
    print(len(toy_stoichemetric_matrix_without_zero_columns[0]))
    index_of_non_zero_columns = np.where(
        np.sum(toy_stoichemetric_matrix_without_zero_columns, axis=0) != 0
    )[0]
    # Sum per column
    sum_per_column = toy_stoichemetric_matrix_without_zero_columns.sum(axis=0)
    print("Sum per column:", sum_per_column)

    # Sum per row
    sum_per_row = toy_stoichemetric_matrix_without_zero_columns.sum(axis=1)
    print("Sum per row:", sum_per_row)
    A = csc_matrix(toy_stoichemetric_matrix_without_zero_columns)
    tol = 1e-9
    opt = 2
    if opt == 2:
        calc_left = 0
        calc_right = 1
    do_pivot = 0
    print("d")
    result = sparseNull(A, tol)
    print(result[0].shape)
    indices = np.nonzero(result[0])
    for row, col in zip(*indices):
        print(f"Row: {row + 1}, Col: {col + 1}, Value: {result[0][row, col]}")
    # l, u, q = luq(A, do_pivot, tol)
    # print(l)
    # print(u)
    # print(q)
    #
    # print(len(l), len(l[0]))
    # print((u.shape))
    # print(len(q), len(q[0]))

    # Lbar2 = L(:, I3)*L1;
    # Qbar2 = Q1 * Q(I4,:);
    # L = [Lbar1 Lbar2 L(:, Jl)];
    # Q = [Qbar1;
    # Qbar2;
    # Q(Jq,:)];
    #
    # n1 = length(I);
    # n2 = length(I3);
    # m2 = length(I4);
    # U = [Ubar1 sparse(n1, m - n1);
    # sparse(n2, n1)
    # U1
    # sparse(n2, m - n1 - m2);
    # sparse(n - n1 - n2, m)];

    # I2 = find(max(abs(Ubar2), [], 2) > tol);
    # I5 = find(max(abs(Ubar2), [], 1) > tol);
    #
    # I3 = Jl(I2);
    # I4 = Jq(I5);
    # Jq(I5) = [];
    # Jl(I2) = [];
    # U = [];
    # A = Ubar2(I2, I5);

    # if size(U, 1) == 1 | | size(U, 2) == 1
    #     S = U(1, 1);
    # else
    #     S = diag(U);
    # end
    # I = find(abs(S) > tol);
    # Jl = (1:n)
    # ';
    # Jl(I) = [];
    # Jq = (1:m)
    # ';
    # Jq(I) = [];
    #
    # Ubar1 = U(I, I);
    # Ubar2 = U(Jl, Jq);
    # Qbar1 = Q(I,:);
    # Lbar1 = L(:, I);

    # TRY!!
    # from sympy import Matrix
    # >> > A = [[2, 3, 5], [-4, 2, 3], [0, 0, 0]]
    # >> > A = Matrix(A)
    # >> > A * A.nullspace()[0]
    # Matrix([
    #     [0],
    #     [0],
    #     [0]])
    # >> > A.nullspace()
    # [Matrix([
    #     [-1 / 16],
    #     [-13 / 8],
    #     [1]])]

    # import numpy as np
    # from scipy.linalg import qr
    #
    # def qr_null(A, tol=None):
    #     Q, R, P = qr(A.T, mode='full', pivoting=True)
    #     tol = np.finfo(R.dtype).eps if tol is None else tol
    #     rnk = min(A.shape) - np.abs(np.diag(R))[::-1].searchsorted(tol)
    #     return Q[:, rnk:].conj()

    # have # been # trying # to # find # a # solution # to # the # same # problem.Using # Scipy # 's svds function provides unreliable results for small singular values. Therefore i have been using QR decomposition instead. The sparseqr https://github.com/yig/PySPQR provides a wrapper for Matlabs SuiteSparseQR method, and works reasonably well. Using this the null space can be calculated as:
    # https://stackoverflow.com/questions/33410146/how-can-i-compute-the-null-space-kernel-x-m-x-0-of-a-sparse-matrix-in-pytho
    # from sparseqr import qr
    # Q, _, _, r = qr(M.transpose())
    # N = Q.tocsr()[:, r:]

    # https: // scicomp.stackexchange.com / questions / 10722 / nullspace - algorithm -
    # for -a - sparse - matrix
    # https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#module-scipy.sparse.linalg
    # Singular
    # values
    # problems:
    #
    # svds(A[, k, ncv, tol, which, v0, maxiter, ...])
    #
    # Partial
    # singular
    # value
    # decomposition
    # of
    # a
    # sparse
    # matrix.
    #
    # The
    # svds
    # function
    # supports
    # the
    # following
    # solvers:
    #
    # svds(solver=’arpack’)
    # svds(solver=’lobpcg’)
    # svds(solver=’propack’)
    # Complete or incomplete
    # LU
    # factorizations
    #
    # splu(A[, permc_spec, diag_pivot_thresh, ...])
    #
    # Compute # the # LU # decomposition # of # a # sparse, square # matrix. # # spilu(A[, drop_tol, fill_factor, drop_rule, ...]) # # Compute # an # incomplete # LU # decomposition # for a sparse, square matrix. # # SuperLU() # # LU # factorization # of # a
    # sparse
    # matrix.

    # https: // stackoverflow.com / questions / 5889142 / python - numpy - scipy - finding - the - null - space - of - a - matrix
    # USE Qr 4
    # A faster but less numerically stable method is to use a rank-revealing QR decomposition, such as scipy.linalg.qr with pivoting=True:

    # # LL = np.block([[np.zeros((n-p, p)), np.eye(p)]])
    # print(f"Shape of L in luq: {L.shape}")
    # print(f"Shape of U in luq: {U.shape}")
    # print(f"Shape of Q in luq: {Q.shape}")
    # print(f"Shape of P in luq: {P.shape}")
    # print(f"p in luq: {p}")
    # print(f"n in luq: {n}")
    # print(f"m in luq: {m}")
    # print(f"Shape of [l.T, P.T[n-p:n, :]] in luq: {P[n-p+1:n, :].T.shape}")
    # L = np.block([[L.T, P.T[n-p:n, :]]])
    # U = np.block([[U.T], [np.zeros((p, m))]])
    # print(f"Shape of L in luq: {L.shape}")
    # print(f"Shape of U in luq: {U.shape}")
    # return L, U, Q


def add_loopless(
    cobra_model: "Model", model, zero_cutoff: Optional[float] = None
) -> None:
    """Modify a model so all feasible flux distributions are loopless.

    It adds variables and constraints to a model which will disallow flux
    distributions with loops. The used formulation is described in [1]_.
    This function *will* modify your model.

    In most cases you probably want to use the much faster
    `loopless_solution`. May be used in cases where you want to add complex
    constraints and objecives (for instance quadratic objectives) to the
    model afterwards or use an approximation of Gibbs free energy directions
    in your model.

    Parameters
    ----------
    cobra_model : cobra.Model
        The model to which to add the constraints.
    model : pe.ConcreteModel
    zero_cutoff : positive float, optional
        Cutoff used for null space. Coefficients with an absolute value
        smaller than `zero_cutoff` are considered to be zero. The default
        uses the `model.tolerance` (default None).

    References
    ----------
    .. [1] Elimination of thermodynamically infeasible loops in steady-state
       metabolic models. Schellenberger J, Lewis NE, Palsson BO. Biophys J.
       2011 Feb 2;100(3):544-53. doi: 10.1016/j.bpj.2010.12.3707. Erratum
       in: Biophys J. 2011 Mar 2;100(5):1381.

    """
    time_start = time.time()
    print("Adding loopless constraints", time.time() - time_start)
    zero_cutoff = normalize_cutoff(cobra_model, zero_cutoff)
    print("cutoff", time.time() - time_start)
    pool_reactions_group = cobra_model.groups.get_by_id("Pool reactions")
    all_reactions = [i for i, r in enumerate(cobra_model.reactions)]
    exchange_rxns = [
        i
        for i, r in enumerate(cobra_model.reactions)
        if r.boundary or r in pool_reactions_group.members
    ]
    exchange_rxns = [
        i for i, r in enumerate(cobra_model.reactions) if r.boundary
    ]  # or r  in pool_reactions_group.members]
    internal = [
        i for i, r in enumerate(cobra_model.reactions) if i not in exchange_rxns
    ]

    s_int = create_stoichiometric_matrix(cobra_model)[:, np.array(internal)]
    n_int, ranks = sparseNull(s_int)
    n_int = n_int.T

    print("created nullspace", time.time() - time_start)

    max_bound = max(max(abs(b) for b in r.bounds) for r in cobra_model.reactions)
    print("max bound", time.time() - time_start)

    print("Adding loopless decisions", time.time() - time_start)
    model.internal_reactions = pe.Set(initialize=internal)
    model.indicator = pe.Var(model.internal_reactions, domain=pe.Binary)
    model.delta_g = pe.Var(model.internal_reactions, domain=pe.Reals)

    def on_off_constraint_rule_1(model, rxn):
        return (
            (-max_bound * (1 - model.indicator[rxn - 1]) <= model.V_flux[rxn])
            if rxn - 1 in internal
            else pe.Constraint.Skip
        )
        # return -max_bound <= model.V_flux[rxn+1] - (model.indicator[rxn] * max_bound)

    def on_off_constraint_rule_2(model, rxn):
        return (
            (model.V_flux[rxn] <= max_bound * model.indicator[rxn - 1])
            if rxn - 1 in internal
            else pe.Constraint.Skip
        )

    model.on_off_constraint_1 = pe.Constraint(model.Rxns, rule=on_off_constraint_rule_1)
    model.on_off_constraint_2 = pe.Constraint(model.Rxns, rule=on_off_constraint_rule_2)

    def delta_g_range_rule_1(model, rxn):
        return (
            (
                -max_bound * model.indicator[rxn - 1]
                + 1 * (1 - model.indicator[rxn - 1])
                <= model.delta_g[rxn - 1]
            )
            if rxn - 1 in internal
            else pe.Constraint.Skip
        )

    def delta_g_range_rule_2(model, rxn):
        return (
            (
                model.delta_g[rxn - 1]
                <= (-1 * model.indicator[rxn - 1])
                + max_bound * (1 - model.indicator[rxn - 1])
            )
            if rxn - 1 in internal
            else pe.Constraint.Skip
        )

    model.delta_g_range_1 = pe.Constraint(model.Rxns, rule=delta_g_range_rule_1)
    model.delta_g_range_2 = pe.Constraint(model.Rxns, rule=delta_g_range_rule_2)

    # def nullspace_constraint_rule(model, int_):
    #     return sum(n_int[int_][idx] * model.delta_g[rxn] for idx, rxn in enumerate(model.internal_reactions)) == 0
    # model.Null = pe.RangeSet(1, len(n_int))
    # model.nullspace_constraint = pe.Constraint(model.Null, rule=nullspace_constraint_rule)

    # Assuming 'model' is your Pyomo model and 'n_int' is your nullspace matrix
    for i, row in enumerate(n_int):
        constraint_name = f"nullspace_constraint_{str(i)}"

        # Create a dictionary for the coefficients
        # print("creating coefs ", i, " ", time.time() - time_start)
        coefs = {ridx: row[i] for i, ridx in enumerate(internal) if abs(row[i]) > 1e-5}

        # Define the rule for the constraint
        def constraint_rule(m):
            return sum(m.delta_g[rxn] * coefs[rxn] for rxn in coefs) == 0

        # Add the constraint to the model
        setattr(model, constraint_name, pe.Constraint(rule=constraint_rule))

    # def on_off_constraint_rule(model, rxn):
    #     return  model.V_flux[rxn + 1] - (model.indicator[rxn] * max_bound) <= 0
    #
    # def delta_g_range_rule_1(model, rxn):
    #     return 1 <= model.delta_g[rxn] + ((max_bound + 1) * model.indicator[rxn])
    #
    # def delta_g_range_rule(model, rxn):
    #     return model.delta_g[rxn] + ((max_bound + 1) * model.indicator[rxn]) <= max_bound
    #
    # def nullspace_constraint_rule(model, int_):
    #     return sum(n_int[int_ - 1][idx] * model.delta_g[rxn] for idx, rxn in enumerate(model.internal_reactions)) == 0

    #     # # Add nullspace constraints for G_i
    #     # for i, row in enumerate(n_int):
    #     #     name = f"nullspace_constraint_{str(i)}"
    #     #     nullspace_constraint = prob.Constraint(Zero, lb=0, ub=0, name=name)
    #     #     model.add_cons_vars([nullspace_constraint])
    #     #     coefs = {
    #     #         model.variables[f"delta_g_{model.reactions[ridx].id}"]: row[i] for i, ridx in enumerate(internal) if abs(row[i]) > zero_cutoff
    #     #     }
    #     #     model.constraints[name].set_linear_coefficients(coefs)
    #
    #
    print(
        len(n_int),
        len(n_int[0]),
    )
    # model.pprint()
    #
    # print("Adding loopless constraints", time.time() - time_start)
    # print("first constraint", time.time() - time_start)
    # model.delta_g_range = pe.Constraint(model.internal_reactions, rule=delta_g_range_rule)
    # print("second constraint", time.time() - time_start)
    # model.on_off_constraint_1 = pe.Constraint(model.internal_reactions, rule=on_off_constraint_rule_1)
    # model.delta_g_range_1 = pe.Constraint(model.internal_reactions, rule=delta_g_range_rule_1)
    # model.Null = pe.RangeSet(1, len(n_int))
    # # model.n_int = pe.Var(model.Null)
    # model.nullspace_constraint = pe.Constraint(model.Null, rule=nullspace_constraint_rule)
    print("Loopless constraints added, took: ", time.time() - time_start)

    return model


# def sparse_null_space(A, rcond=None):
#     """
#     Compute the null space of a sparse matrix A.
#     """
#     print("sparse null space")
#     time_start = time.time()
#     u, s, vh = svds(A, k=min(A.shape) - 1)  # Compute the SVD
#     print("svd", time.time() - time_start)
#     tol = np.max(A.shape) * np.spacing(np.max(s)) if rcond is None else rcond
#     print("tol", time.time() - time_start)
#     null_mask = (s <= tol)
#     print("null mask", time.time() - time_start)
#     null_space = np.compress(null_mask, vh, axis=0)
#     print("null space", time.time() - time_start)
#     return null_space.T
#
#
# def sparse_nullspace(A, atol=1e-13, rtol=0):
#     time_start = time.time()
#     A = csc_matrix(A)
#     print("csc", time.time() - time_start)
#     u, s, vh = svds(A)
#     print(u, "\n")
#     print(s, "\n")
#     print(vh, "\n")
#     print("svd", time.time() - time_start)
#     tol = max(atol, rtol * s[0])
#     print("tol", time.time() - time_start)
#     nnz = (s >= tol).sum()
#     print("nnz", time.time() - time_start)
#     ns = vh[nnz:].conj().T
#     return ns
#
# def add_loopless(cobra_model: "Model", model, zero_cutoff: Optional[float] = None) -> None:
#     """Modify a model so all feasible flux distributions are loopless.
#
#     It adds variables and constraints to a model which will disallow flux
#     distributions with loops. The used formulation is described in [1]_.
#     This function *will* modify your model.
#
#     In most cases you probably want to use the much faster
#     `loopless_solution`. May be used in cases where you want to add complex
#     constraints and objecives (for instance quadratic objectives) to the
#     model afterwards or use an approximation of Gibbs free energy directions
#     in your model.
#
#     Parameters
#     ----------
#     cobra_model : cobra.Model
#         The model to which to add the constraints.
#     model : pe.ConcreteModel
#     zero_cutoff : positive float, optional
#         Cutoff used for null space. Coefficients with an absolute value
#         smaller than `zero_cutoff` are considered to be zero. The default
#         uses the `model.tolerance` (default None).
#
#     References
#     ----------
#     .. [1] Elimination of thermodynamically infeasible loops in steady-state
#        metabolic models. Schellenberger J, Lewis NE, Palsson BO. Biophys J.
#        2011 Feb 2;100(3):544-53. doi: 10.1016/j.bpj.2010.12.3707. Erratum
#        in: Biophys J. 2011 Mar 2;100(5):1381.
#
#     """
#     time_start = time.time()
#     print("Adding loopless constraints", time.time() - time_start)
#     zero_cutoff = normalize_cutoff(cobra_model, zero_cutoff)
#     print("cutoff", time.time() - time_start)
#
#     internal = [i for i, r in enumerate(cobra_model.reactions) if not r.boundary]
#     length_internal = len(internal)
#     # internal = internal[:int(length_internal*0.2)]
#     s_int = create_stoichiometric_matrix(cobra_model)[:, np.array(internal)]
#
#     atol= 1e-13
#     rtol = 0.0
#     A = np.atleast_2d(s_int)
#     _, s, vh = np.linalg.svd(A)
#     tol = max(atol, rtol * s[0])
#     nnz = (s >= tol).sum()
#     ns = vh[nnz:].conj().T
#     n_int = ns.T
#     print("nullspace", time.time() - time_start)
#
#
#     # A = s_int
#     # A = np.atleast_2d(A)
#     # # print(np.shape(A))
#     # #
#     # # #
#     # # B = csc_matrix(A)
#     # # u, s, vh = svds(B)
#     # # print("u: \n ", u, "\n")
#     # # print("s: \n ", s, "\n")
#     # # print("vh: \n ", vh, "\n")
#     # #
#     # # print("svds", time.time() - time_start)
#     # # tol = max(atol, rtol * s[0])
#     # # nnz = (s >= tol).sum()
#     # # ns = vh[nnz:].conj().T
#     # # n_int = ns.T
#     # # print(n_int)
#     #
#     #
#     #
#     # np_u, np_s, np_vh = np.linalg.svd(A)
#     # print("np_s: \n ", np_s, "\n")
#     # print("np_vh: \n ", np_vh, "\n")
#     # tol = max(atol, rtol * np_s[0])
#     # nnz = (np_s >= tol).sum()
#     # ns = np_vh[nnz:].conj().T
#     # n_int = ns.T
#     # print("numpy svd", time.time() - time_start)
#     # print(n_int)
#     # print(f"len internal: {len(n_int)}")
#     #
#     #
#     #
#     # sp_u, sp_s, sp_vh = scipy.linalg.svd(A)
#     # print("sp_s: \n ", sp_s, "\n")
#     # print("sp_vh: \n ", sp_vh, "\n")
#     # tol = max(atol, rtol * sp_s[0])
#     # nnz = (sp_s >= tol).sum()
#     # ns = sp_vh[nnz:].conj().T
#     # n_int = ns.T
#     # print("scipy svd", time.time() - time_start)
#     # print(n_int)
#     # print(f"len internal: {len(n_int)}")
#     #
#
#
#     print("created nullspace", time.time() - time_start)
#
#     max_bound = max(max(abs(b) for b in r.bounds) for r in cobra_model.reactions)
#     print("max bound", time.time() - time_start)
#     # # prob = model.problem
#     #
#     # stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(cobra_model)
#     # model = pe.ConcreteModel()
#     # num_metabolites = len(stoichiometric_matrix)  # Example number of metabolites
#     # num_reactions = len(stoichiometric_matrix[0])  # Example number of reactions
#     #
#     # model.Metabolites = pe.RangeSet(1, num_metabolites)
#     # model.Rxns = pe.RangeSet(1, num_reactions)
#     # stoich_dict = {met: {rxn: stoichiometric_matrix[met, rxn] for rxn in range(num_reactions) if stoichiometric_matrix[met, rxn] != 0} for met in
#     #                range(num_metabolites)}
#     # model.V_flux = pe.Var(model.Rxns, domain=pe.Reals)
#     # for rxn in model.Rxns:
#     #     model.V_flux[rxn].setlb(cobra_model.reactions[rxn - 1].lower_bound)
#     #     model.V_flux[rxn].setub(cobra_model.reactions[rxn - 1].upper_bound)
#     #
#     # def sv_zero_rule(model, met):
#     #     return sum(stoich_dict[met - 1][rxn] * model.V_flux[rxn + 1] for rxn in stoich_dict[met - 1]) == 0
#     #
#     # model.SV_constraint = pe.Constraint(model.Metabolites, rule=sv_zero_rule)
#
#
#     print("Adding loopless decisions", time.time() - time_start)
#     model.internal_reactions = pe.Set(initialize=internal)
#     model.indicator = pe.Var(model.internal_reactions, domain=pe.Binary)
#     model.delta_g = pe.Var(model.internal_reactions, domain=pe.Reals)
#
#     def on_off_constraint_rule_1(model, rxn):
#        return -max_bound <= model.V_flux[rxn+1] - (model.indicator[rxn] * max_bound)
#
#     def on_off_constraint_rule(model, rxn):
#         return  model.V_flux[rxn + 1] - (model.indicator[rxn] * max_bound) <= 0
#
#     def delta_g_range_rule_1(model, rxn):
#         return 1 <= model.delta_g[rxn] + ((max_bound + 1) * model.indicator[rxn])
#
#     def delta_g_range_rule(model, rxn):
#         return model.delta_g[rxn] + ((max_bound + 1) * model.indicator[rxn]) <= max_bound
#
#     def nullspace_constraint_rule(model, int_):
#         return sum(n_int[int_ - 1][idx] * model.delta_g[rxn] for idx, rxn in enumerate(model.internal_reactions)) == 0
#
#     print(len(n_int), len(n_int[0]), )
#     #
#     # # Add nullspace constraints for G_i
#     # for i, row in enumerate(n_int):
#     #     name = f"nullspace_constraint_{str(i)}"
#     #     nullspace_constraint = prob.Constraint(Zero, lb=0, ub=0, name=name)
#     #     model.add_cons_vars([nullspace_constraint])
#     #     coefs = {
#     #         model.variables[f"delta_g_{model.reactions[ridx].id}"]: row[i] for i, ridx in enumerate(internal) if abs(row[i]) > zero_cutoff
#     #     }
#     #     model.constraints[name].set_linear_coefficients(coefs)
#
#
#     print("Adding loopless constraints", time.time() - time_start)
#     model.on_off_constraint = pe.Constraint(model.internal_reactions, rule=on_off_constraint_rule)
#     print("first constraint", time.time() - time_start)
#     model.delta_g_range = pe.Constraint(model.internal_reactions, rule=delta_g_range_rule)
#     print("second constraint", time.time() - time_start)
#     model.on_off_constraint_1 = pe.Constraint(model.internal_reactions, rule=on_off_constraint_rule_1)
#     model.delta_g_range_1 = pe.Constraint(model.internal_reactions, rule=delta_g_range_rule_1)
#     model.Null = pe.RangeSet(1, len(n_int))
#     # model.n_int = pe.Var(model.Null)
#     model.nullspace_constraint = pe.Constraint(model.Null, rule=nullspace_constraint_rule)
#     print("Loopless constraints added, took: ", time.time() - time_start)
#
#     return model
#
#     # # Add indicator variables and new constraints
#     # to_add = []
#     # for i in internal:
#     #     rxn = model.reactions[i]
#     #     # indicator variable a_i
#     #     indicator = prob.Variable(f"indicator_{rxn.id}", type="binary")
#     #     # -M*(1 - a_i) <= v_i <= M*a_i
#     #     on_off_constraint = prob.Constraint(
#     #        # -maxbound <= rxn.flux_expression - max_bound * indicator <= 0,
#     #         lb=-max_bound,
#     #         ub=0,
#     #         name=f"on_off_{rxn.id}",
#     #         )
#     #     # -(max_bound + 1) * a_i + 1 <= G_i <= -(max_bound + 1) * a_i + 1000
#     #     delta_g = prob.Variable(f"delta_g_{rxn.id}")
#     #     delta_g_range = prob.Constraint(
#     #         delta_g + (max_bound + 1) * indicator,
#     #         lb=1,
#     #         ub=max_bound,
#     #         name=f"delta_g_range_{rxn.id}",
#     #         )
#     #     to_add.extend([indicator, on_off_constraint, delta_g, delta_g_range])
#     #
#     # model.add_cons_vars(to_add)
