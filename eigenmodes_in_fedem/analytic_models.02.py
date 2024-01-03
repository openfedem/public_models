"""
https://www.intechopen.com/chapters/72027
http://www.varg.unsw.edu.au/Assets/link%20pdfs/Beam_vibration.pdf
"""
from scipy.optimize import fsolve
#from numpy.linalg import matrix_rank
from scipy.sparse.linalg import svds
from scipy.sparse.linalg import eigsh
import scipy.sparse.linalg as spla
from scipy.linalg import eigh

from scipy.sparse import lil_matrix
import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
import copy

PRINT = True
PLOT = True

pi2 = np.pi ** 2.0


def section_properties(bar=None):
    # DEPRICATED ad replaced with make_cross_Section_properties
    # A, I
    if bar:
        return dict(area=np.pi*bar[0]**2.0, second_moment_of_area=np.pi*bar[0]**4.0/4.0)


def make_cross_section_properties(*vargs):
    # A, I
    if vargs[0] == "bar":
        radius = vargs[1]
        return np.pi*radius**2.0, np.pi*radius**4.0/4.0


def get_requested_modes(modes):
    if isinstance(modes, int):
        requested_modes = range(1, modes+1)
    else:
        requested_modes = modes

    return requested_modes


def beta_for_cantilever(length, num_roots=5, accuracy=1e-6):
    """
    Finds multiple roots of the equation 1 + cos(beta * length) * cosh(beta * length) = 0 for beta.

    Parameters:
    length (float): The value of length in the equation.
    num_roots (int): Number of roots to find.
    accuracy (float): The desired accuracy of the solution.

    Returns:
    list: A list of solutions for beta.
    """
    #print("beta_for_cantilever")
    roots = []

    # Function to find the root near a given guess
    def find_root_near(guess):
        root, = fsolve(lambda beta: 1.0 + np.cos(beta * length) * np.cosh(beta * length), guess, xtol=accuracy)
        return root

    # Initial guess for the first root
    guess = 0.1 / length
    for i in range(num_roots):
        root = find_root_near(guess)
        # Avoid duplicates
        if not any(np.isclose(root, r, atol=accuracy) for r in roots):
            roots.append(root)
        # Update the guess for the next root
        guess = 0.1*(root + np.pi / length)


    #for i in range(num_roots):
    #    n = i + 1
    #    print((2*n+1)*np.pi/length)
    #    roots[i] = (2*n+1)*np.pi/length
    #print(roots)

    roots = [3.516, 22.03, 61.70, 120.9]#, 17.2787]

    return roots


def beta_for_unsupported_beam(length, num_roots=5, accuracy=1e-6):
    """
    cosh(kL)cos(kL)=1
    """
    #print("beta_for_unsupported_beam")
    roots = []

    # Function to find the root near a given guess
    def find_root_near(guess):
        root, = fsolve(lambda beta: 1.0 - np.cos(beta * length) * np.cosh(beta * length), guess, xtol=accuracy)
        return root

    # Initial guess for the first root
    guess = 0.1 / length
    for i in range(num_roots):
        root = find_root_near(guess)
        # Avoid duplicates
        if not any(np.isclose(root, r, atol=accuracy) for r in roots):
            roots.append(root)
        # Update the guess for the next root
        guess = root + np.pi / length
    #print(roots)

    roots = [4.7300/length, 7.8532/length, 10.9956/length, 14.1371/length]

    return roots


def alpha_unsupported_beam(n_modes):
    """
    Calculate the first 'n_modes' alpha values for a free-free beam.

    Args:
    n_modes (int): Number of modes for which to calculate alpha values.

    Returns:
    list: Alpha values for the first 'n_modes' modes.
    """

    def transcendental_eq(alpha, mode):
        # Transcendental equation for free-free beam
        # This is a simplified representation and might need adjustments
        # based on specific beam characteristics.
        return 1 - np.cos(alpha) * np.cosh(alpha)

    alpha_values = []
    for mode in range(1, n_modes + 1):
        # Initial guess for alpha, can be adjusted for better convergence
        alpha_guess = (mode - 0.5) * np.pi
        alpha = fsolve(transcendental_eq, alpha_guess, args=(mode))
        alpha_values.append(alpha[0])

    #print("alpha", alpha_values)

    alpha_values = [4.7300, 7.8532, 10.9956, 14.1371,  17.2787][:n_modes]

    return alpha_values


def unsupported_beam(beam_length, density, elasticity_module, area, second_moment_of_area, modes=5):

    if True:
        # Example: Calculate the first 5 alpha values for a free-free beam
        alpha_values = alpha_unsupported_beam(modes)

        unit_mass = density * area

        eigenvalues = []
        for alpha in alpha_values:
            eigenvalues.append(alpha**2.0 *
                               np.sqrt(elasticity_module * second_moment_of_area / (unit_mass*beam_length**4.0)))

        return eigenvalues

    else:
        beta_values = beta_for_unsupported_beam(beam_length, num_roots=modes, accuracy=ACCURACY)

        unit_mass = density * area

        eigenvalues = []
        for beta in beta_values:
            eigenvalues.append(beta / beam_length ** 2.0 *
                               np.sqrt(elasticity_module * second_moment_of_area / unit_mass))

        return eigenvalues


def supported_beam(beam_length, density, elasticity_module, area, second_moment_of_area, modes=5):
    requested_modes = get_requested_modes(modes)

    unit_mass = density * area

    eigenvalues = []
    for mode_number in requested_modes:
        eigenvalues.append(mode_number ** 2.0 * pi2 / beam_length ** 2.0 *
                           np.sqrt(elasticity_module * second_moment_of_area / unit_mass))

    return eigenvalues


def cantilever_beam(beam_length, density, elasticity_module, area, second_moment_of_area, modes=5):

    beta_values = beta_for_cantilever(beam_length, num_roots=modes, accuracy=ACCURACY)
    unit_mass = density * area

    eigenvalues = []
    for beta in beta_values:
        eigenvalues.append(beta / beam_length ** 2.0 *
                           np.sqrt(elasticity_module * second_moment_of_area / unit_mass))

    return eigenvalues


def make_element_stiffness(element_length, density, elasticity_module, area, second_moment_of_area):
    L = element_length
    L2 = element_length**2.0

    k = (elasticity_module*second_moment_of_area/element_length**3.0)*np.array([
            [0.0]*6,
            [0.0,  12.0,   6.0*L,  0.0, -12.0,    6.0*L],
            [0.0,   6.0*L, 4.0*L2, 0.0,  -6.0*L,  2.0*L2],
            [0.0]*6,
            [0.0, -12.0,  -6.0*L,  0.0,  12.0,  -6.0*L],
            [0.0,   6.0*L, 2.0*L2, 0.0, -6.0*L,  4.0*L2]])

    k[0, 0] = k[3, 3] = elasticity_module * area / L
    k[0, 3] = k[3, 0] = -k[0, 0]

    m = (density*area*element_length/420.0)*np.array([
            [140.0,    0.0,    0.0,   70.0,   0.0,     0.0],
            [  0.0,  156.0,   22.0*L,  0.0,  54.0,   -13.0*L],
            [  0.0,   22.0*L,  4.0*L2, 0.0,  13.0*L,  -3.0*L2],
            [ 70.0,    0.0,    0.0,  140.0,   0.0,     0.0],
            [  0.0,   54.0,   13.0*L,  0.0, 156.0,   -22.0*L],
            [  0.0,  -13.0*L, -3.0*L2, 0.0, -22.0*L,   4.0*L2]])

    return k, m


def append_member_to_system(member_node_names, point_start_, point_end_,
                            material_property, cross_section_property,
                            dofs, system_stiffness, system_mass):

    x1 = np.array(point_start_)
    x2 = np.array(point_end_)
    dx = x2 - x1
    beam_length = np.linalg.norm(dx)
    beam_orientation = np.arctan2(*np.flip(dx))
    density, elasticity_module = material_property
    area, second_moment_of_area = cross_section_property

    element_length = beam_length / (len(member_node_names)-1)

    _element_stiffness, _element_mass = make_element_stiffness(element_length, density, elasticity_module, area, second_moment_of_area)

    transformation_matrix = make_transformation_matrix(beam_orientation)

    element_stiffness = transformation_matrix.T @ _element_stiffness @ transformation_matrix
    element_mass = transformation_matrix.T @ _element_mass @ transformation_matrix
    print("Equations")
    for node_start, node_end in zip(member_node_names[:-1], member_node_names[1:]):
        eqs = np.concatenate([dofs[node_start], dofs[node_end]])-1

        print(node_start, node_end, eqs)

        system_stiffness[np.ix_(eqs, eqs)] += element_stiffness
        system_mass[np.ix_(eqs, eqs)] += element_mass
    return system_stiffness, system_mass


def make_system_nodes_and_dofs(member_name, point_start, point_end, x1, x2, number_of_elements, member_nodes={}, dofs={}):
    dx = (x2-x1)/number_of_elements

    if point_start not in member_nodes:
        member_nodes[point_start] = x1
        dofs[point_start] = np.array([1, 2, 3]) + len(dofs) * 3

    for i in range(1, number_of_elements):
        node_name = f"{member_name}.{i}"
        member_nodes[node_name] = x1 + dx * i
        dofs[node_name] = np.array([1, 2, 3]) + len(dofs) * 3

    if point_end not in member_nodes:
        member_nodes[point_end] = x2
        dofs[point_end] = np.array([1, 2, 3]) + len(dofs) * 3


    return member_nodes, dofs


def append_boundary_conditions(K, M, system_dofs, bcs, inplace=False):
    if inplace:
        model_dofs = system_dofs
    else:
        model_dofs = copy.deepcopy(system_dofs)

    fixed_dofs = []
    for bc in bcs:
        i = bc.rfind(".")
        node, node_dofs = bc[:i], bc[i+1:]
        for node_dof in [int(i) for i in node_dofs]:  #local dof:  1, 2 or 3
            fixed_dofs.append(model_dofs[node][node_dof-1] )  # Minus 1 as dof starts on 1
            model_dofs[node][node_dof-1] = - model_dofs[node][node_dof-1]

    equations = [num-1 for sublist in model_dofs.values() for num in sublist if num > 0]
    _K, _M = K[np.ix_(equations, equations)], M[np.ix_(equations, equations)]

    if PRINT > 2:
        print(equations)
        print(K.shape, _K.shape)
        print(M.shape, _M.shape)
        print(pd.DataFrame(K.toarray()).to_markdown())
        print(pd.DataFrame(_K.toarray()).to_markdown())
    return _K, _M, equations, model_dofs


def make_transformation_matrix(theta, degrees=False):
    if degrees:
        _theta = np.radians(45)  # Example angle in degrees
    else:
        _theta = theta

    c = np.cos(_theta)
    s = np.sin(_theta)
    T_node = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    T = np.block([[T_node, np.zeros((3, 3))], [np.zeros((3, 3)), T_node]])
    return T


def plot(nodes, dofs):
    for member_name, member_nodes in nodes.items():

        xys = np.array([l for l in member_nodes.values()])
        #for node_name, [x, y] in member_nodes[1:].items():

        plt.plot(xys[:, 0], xys[:, 1])

    plt.show()


def make_fea_model(material_properties, cross_section_properties, points, members, boundary_conditions):
    nodes = {}
    dofs = {}

    for member_name, (material_name, cross_section_name, point_start, point_end, number_of_elements) in members.items():
        print(member_name, nodes)
        x1 = np.array(points[point_start])
        x2 = np.array(points[point_end])
        nodes[member_name] = {}
        make_system_nodes_and_dofs(member_name, point_start, point_end, x1, x2, number_of_elements, nodes[member_name], dofs)

    if PRINT:
        print("Nodes:")
        for member_name, member_nodes in nodes.items():
            print(f"  - Member: {member_name}")
            for k, (x, y) in member_nodes.items():
                print(f"    - {k:5}{x:8.3f}{y:8.3f}")

        print("\nDofs, (Degrees of freedom, aka equation numbers) for unsupported system:")
        for node_name, node_dofs in dofs.items():
            print(f"  - {node_name:5}{node_dofs[0]:4}{node_dofs[1]:4}{node_dofs[2]:4}")

    size = max(max(equations) for equations in dofs.values())
    if PRINT:
        print(f"\nSystem size: {size}")

    system_stiffness = lil_matrix((size, size))
    system_mass = lil_matrix((size, size))
    for member_name, (material_name, cross_section_name, point_start, point_end, number_of_elements) in members.items():
        append_member_to_system(
            list(nodes[member_name].keys()), points[point_start], points[point_end],
            material_properties[material_name],
            make_cross_section_properties(*cross_section_properties[cross_section_name]),
            dofs, system_stiffness, system_mass)

    system_stiffness = system_stiffness.tocsr()
    system_mass = system_mass.tocsr()

    models = {}
    if PRINT:
        print("\nModels:")

    for model_name, boundary_condition in boundary_conditions.items():
        models[model_name] = append_boundary_conditions(system_stiffness, system_mass, dofs, boundary_condition)

        if PRINT:
            system_dofs = models[model_name][-1]
            print(f"  - Model dofs for {model_name}")
            for node_name, node_dofs in system_dofs.items():
                print(f"    - {node_name:5}{node_dofs[0]:4}{node_dofs[1]:4}{node_dofs[2]:4}")

    plot(nodes, dofs)

    return models


def solve_modes(K, M, modes, lower_cutoff=1.0e-2):
    eigenvalues, eigenvectors = eigh(K.toarray(), M.toarray())

    #print( type(eigenvalues), type(eigenvalues[0]))

    eigenvalues = np.sqrt(eigenvalues)
    # Eigenfrequencies (square roots of eigenvalues give natural frequencies)
    # Eigenmodes (eigenvectors)

    skip_cutoff = np.argmax(eigenvalues>lower_cutoff)

    return eigenvalues[skip_cutoff:modes + skip_cutoff]



input_data = yaml.safe_load(open("analytic_models.yaml"))
ACCURACY = input_data.get("ACCURACY", 1.0E-8)

number_of_modes = input_data["number_of_modes"]

input_model = input_data["models"]["A"]
model_A = make_fea_model(**input_model)

eigen_frequencies_df = pd.DataFrame(index=range(1, number_of_modes+1))
for name, bc in input_model["boundary_conditions"].items():
        eigen_frequencies_df[f"{name}_fea"] = solve_modes(*model_A[name][:2], number_of_modes)


input_model_B = input_data["models"]["B"]
model_B = make_fea_model(**input_model_B)

#eigen_frequencies_df = pd.DataFrame(index=range(1, number_of_modes+1))
for name, bc in input_model_B["boundary_conditions"].items():
        eigen_frequencies_df[f"{name}_fea"] = solve_modes(*model_B[name][:2], number_of_modes)


"""
import numpy as np
from scipy.sparse.linalg import svds

# Create a sparse matrix
A = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
sparse_A = sps.csr_matrix(A)

# Compute the rank of the sparse matrix
u, s, vt = svds(sparse_A)
rank = len([x for x in s if abs(x) > 1e-12])

print(f"The rank of the sparse matrix is {rank}.")

"""






kwargs = dict(beam_length=1.0, density=8000.0, elasticity_module=2.0e11, **section_properties(bar=[0.005]))
#eigen_frequencies_df = pd.DataFrame(index=range(1, number_of_modes+1))
eigen_frequencies_df["unsupported_beam"] = unsupported_beam(**kwargs, modes=number_of_modes)
eigen_frequencies_df["supported_beam"] = supported_beam(**kwargs, modes=number_of_modes)
eigen_frequencies_df["cantilever_beam"] = cantilever_beam(**kwargs, modes=number_of_modes)

cols = [c for c in eigen_frequencies_df.columns if c.startswith("supported")]
print(eigen_frequencies_df[cols])

cols = [c for c in eigen_frequencies_df.columns if c.startswith("ca")]
print(eigen_frequencies_df[cols])

cols = [c for c in eigen_frequencies_df.columns if c.startswith("unsupported")]
print(eigen_frequencies_df[cols])


