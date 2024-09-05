import numpy as np
from fedempy.modeler import FedemModeler
from fedempy.fmm_solver import FmmSolver
from fedempy.enums import FmDof, FmDofStat, FmLoadType, FmType, FmVar

import a_make_spool_model


def run_spool_model():
    mySolver = FmmSolver()

    mySolver.solve_all(a_make_spool_model.FMM, True, True)
    status = mySolver.start(a_make_spool_model.FMM, False, True)
    if status < 0:
        print(" *** Something bad happened, check fedem_solver.res", status)

if __name__ == "__main__":
    a_make_spool_model.make_spool_model()
    run_spool_model()