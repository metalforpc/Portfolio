##### Robust Identification and Estimation in Repeated Games
## Antonio Cozzolino, Cristina Gualdani, Lorenzo Magnolfi, Niccolò Lomys

### This is a very simple example of code for the paper  
### "Robust Identification and Inference in Repeated Games" 
### The full code will be available when we circulate the working paper. I kept
### this code simple so that can be run with a simple installation of Python.

### This code illustrate an example of identification routine to identify a repeated 
### Prisoner's Dilemma.

### Import Section, also cython elements to speed up the identification
### function
import cython
import numpy as np
cimport numpy as cnp
cnp.import_array()
DTYPE = np.float64
ctypedef cnp.float64_t DTYPE_t

from scipy.stats import uniform
from scipy.optimize import linprog

class Robust_Identification():
    
    # Constructor of the class Robust Identification
    def __init__(self, ecdf, float β_lb = 0.5, float β_ub = 1.99, float δ_lb = 0.01, float δ_ub = 0.99, float η_lb = 0.01, float η_ub = 0.99):
        """Identification routine for PrisonersDilemma.
            Takes as input an Empirical Distribution of Action Profiles, a 4x1 vector:

            [P(C,C), P(C,D), P(D,C), P(D, D)]"""

        # Global Variables
        self.ecdf = ecdf

        # Bounds on payoff for (C, C)
        self.β_lb = β_lb
        self.β_ub = β_ub

        # Bounds for the discount factor
        self.δ_lb = δ_lb
        self.δ_ub = δ_ub

        # Bounds on the quality of monitoring
        self.η_lb = η_lb
        self.η_ub = η_ub

        # Random Sample from a uniform distribution for numerical integration
        self.Ω = np.random.uniform(-1/2, 1/2, 1000)

    ## Linear Program for solving the max and min
    def continuation_value(self, payoff_matrix, δ, η):
        
        # For CCE Definition
        ε = (2*payoff_matrix[1,1]*η*(δ**2))/(1 - δ)

        # Matrix of Linear Inequality for the Linear Program
        Aub = np.array([[0, 0, (payoff_matrix[0,0] - payoff_matrix[2,0]), (payoff_matrix[1,0] - payoff_matrix[3,0])],
                [(- payoff_matrix[0,0] + payoff_matrix[2,0]), (- payoff_matrix[1,0] + payoff_matrix[3,0]), 0, 0],
                [0, (payoff_matrix[0,1] - payoff_matrix[1,1]), 0, (payoff_matrix[2,1] - payoff_matrix[3,1])],
                [(-payoff_matrix[0,1] + payoff_matrix[1,1]), 0, (-payoff_matrix[2,1] + payoff_matrix[3,1]), 0],
                [-1, 0, 0, 0],
                [0, -1, 0, 0],
                [0, 0, -1, 0],
                [0, 0, 0, -1],
                [-payoff_matrix[0,0], -payoff_matrix[1,0], -payoff_matrix[2,0], -payoff_matrix[3,0]],
                [-payoff_matrix[0,1], -payoff_matrix[1,1], -payoff_matrix[2,1], -payoff_matrix[3,1]]])
        
        # Upper bound 
        bub = np.array([ε, ε, ε, ε, 0, 0, 0, 0, 0, 0])

        # Probability Constraints
        Aeq = np.array([[1, 1, 1, 1]])
        beq = np.array([1])

        # Vector of coefficient for probability distributions
        c = np.array([payoff_matrix[0,0], payoff_matrix[1,0], payoff_matrix[2,0], payoff_matrix[3,0]])

        # Solve the linear programming finding the CCE that attains the maximum
        # continuation value
        res = linprog(c=-c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq)

        M = -res.fun # Maximum Continuation Value
        m = -M  # Minimum Continuation Value

        return M, m
    
        ## Main Identification Routine
    def identify(self, int grid):

        # Parameter space for (β, δ, η)
        cdef cnp.ndarray B = np.linspace(self.β_lb, self.β_ub, grid, dtype=DTYPE)
        cdef cnp.ndarray Δ = np.linspace(self.δ_lb, self.δ_ub, grid, dtype=DTYPE)
        cdef cnp.ndarray H = np.linspace(self.η_lb, self.η_ub, grid, dtype=DTYPE)

        # Instantiate the Identified Set
        ID_set = np.zeros((grid, grid, grid))

        # For each (β, δ, η) in the domain, at this stage we are still developing a method
        # to avoid searching over the entire space. This part is inefficient because 
        # this algorithm spends most of the time in points not identified. 
        for idx_β in range(0, grid):
            for idx_δ in range(0, grid):
                for idx_η in range(0, grid):
                    
                    b = B[idx_β]
                    d = Δ[idx_δ]
                    e = H[idx_η] 

                    # Instantiate payoff matrix
                    payoff_mat = np.array([[b, b], [-1, 2], 
                                            [2, -1], [0, 0]])
                    
                    # Compute the continuation value
                    M, m = self.continuation_value(payoff_mat, d, e)
                        
                    # Sufficient and Necessary conditions for each action profile, numerical integration
                    cc_lb = np.sum(((1-d)*np.minimum(b + self.Ω - 2, -1) + d*m >= 0))/np.shape(self.Ω)[0]
                    cc_ub = np.sum(((1-d)*np.maximum(b + self.Ω - 2, -1) + d*M >= 0))/np.shape(self.Ω)[0]
                    
                    dd_lb = np.sum(((1-d)*np.minimum(2 - b - self.Ω, 1) + d*m >= 0))/np.shape(self.Ω)[0]
                    dd_ub = np.sum(((1-d)*np.maximum(2 - b - self.Ω, 1) + d*M >= 0))/np.shape(self.Ω)[0]

                    cd_lb = np.sum(((1 - d)*np.minimum(b + self.Ω - 2, -1) + d*m >= 0)*((1 - d)*np.minimum(2 - b - self.Ω, 1) + d*m >= 0))/np.shape(self.Ω)[0]
                    cd_ub = np.sum(((1 - d)*np.maximum(b + self.Ω - 2, -1) + d*M >= 0)*((1 - d)*np.maximum(2 - b - self.Ω, 1) + d*M >= 0))/np.shape(self.Ω)[0]

                    dc_lb = cd_lb
                    dc_ub = cd_ub

                    # Identifying conditiond
                    identifying_cond = (self.ecdf[0] >= cc_lb and self.ecdf[0] <= cc_ub) and (self.ecdf[1] >= cd_lb and self.ecdf[1] <= cd_ub) and (self.ecdf[2] >= dc_lb and self.ecdf[2] <= dc_ub) and (self.ecdf[3] >= dd_lb and self.ecdf[3] <= dd_ub)
                    
                    # If identified, put in the set otherwise continue
                    # the identified set can be summarized as an array of 0 and 1 and use
                    # this mask as Carthesian Coordinates.
                    if identifying_cond:
                        ID_set[idx_β, idx_δ, idx_η] = 1
        
        # Generate MashGrid to return
        B, Δ, H = np.meshgrid(B, Δ, H)

        return ID_set, B, Δ, H

    