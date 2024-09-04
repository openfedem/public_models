import numpy as np
from fedempy.modeler import FedemModeler
from fedempy.enums import FmDof, FmDofStat, FmLoadType, FmType, FmVar

FMM_TEMPLATE = "0_spool_model.fmm"
FMM = ("1_spool_model.fmm")

def make_spool_model():
    spool_model = FedemModeler(FMM_TEMPLATE)

    # Make triads and beams
    ## Points in 2D defining the spool nodes i.e. triads
    spool_points = np.array([
        [-5,3], [-5,1], [-5,-1], [-4,-2], [-2,-2], [2,-2], [4,0], [4,2], [4,4]], dtype=float) # Define the 2D points of the spool

    ## Appending triads to the model
    triads = [spool_model.make_triad(f"Triad {i}", (x, y, 0.0)) for i, (x, y) in enumerate(spool_points)]

    ## Appending material and beam cross-section properties
    steel_material = spool_model.make_beam_material("Steel", (7850, 2.1e11, 0.3))
    beam_property = spool_model.make_beam_section("Pipe", steel_material, (0.5, 0.45))

    ## Appending beams to the model
    beams = spool_model.make_beam("Spool beams", triads, beam_property)

    ## Defines the 3D point of the hook and appends a triad with mass for that point
    hook_triad = spool_model.make_triad(f"Hook Triad", (-0.43, -0.27 , 5.0))
    spool_model.edit_triad(hook_triad, mass=(500,100,100,100))

    ## Appends springs that represents the slings
    lifting_triads = [triads[i] for i in [1, 4, 7]]
    k = 2.1e11*0.020**2/(4*4)  # E x D^2 / (4 * L)
    xy = [[0.0, k]]
    slings = [spool_model.make_spring("Slings", (lifting_triad, hook_triad),
                                      xy=xy, extrapol_type="FLAT") for lifting_triad in lifting_triads]

    ## Appends a triad for the crane tip as well as a spring representing the crane wire
    crane_tip_triad = spool_model.make_triad(f"Crane tip triad", (-0.43, -0.27, 8.0))
    k = 2.1e11*0.020**2/(4*4)  # E x D^2 / (4 * L)
    xy = [[0.0, 0.0], [0.1, k]]
    cable = spool_model.make_spring("Slings", (hook_triad, crane_tip_triad), xy=xy, extrapol_type="FLAT")

    spool_model.edit_triad(crane_tip_triad,
                           constraints=dict(
                               Tx=FmDofStat.FIXED,
                               Ty=FmDofStat.FIXED,
                               Tz=FmDofStat.PRESCRIBED,
                               Rx=FmDofStat.FIXED,
                               Ry=FmDofStat.FIXED,
                               Rz=FmDofStat.FIXED),
                           motion={"Tz": [[0, 0], [30, -10]]},
                           motion={"Tz": -.1}
                           )

    try:
        print("Saving the model...")
        spool_model.save(FMM)
    except Exception as e:
        print(f"Exception({e})")
        raise e   # "Exception({"Error": "Failed to save current model"})

    spool_model.close()

if __name__ == "__main__":
    make_spool_model()