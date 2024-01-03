"""
https://www.intechopen.com/chapters/72027
http://www.varg.unsw.edu.au/Assets/link%20pdfs/Beam_vibration.pdf
"""
from scipy.optimize import fsolve
from scipy.linalg import eigh
import pandas as pd
import numpy as np
import yaml


pi2 = np.pi ** 2.0


def section_properties(bar=None):
    # A, I
    if bar:
        return dict(area=np.pi*bar[0]**2.0, second_moment_of_area=np.pi*bar[0]**4.0/4.0)


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
    print("beta_for_cantilever")
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
    print("beta_for_unsupported_beam")
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
    print(roots)

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
            [ 12.0,   6.0*L, -12.0,    6.0*L],
            [  6.0*L, 4.0*L2, -6.0*L,  2.0*L2],
            [-12.0,  -6.0*L,   12.0,  -6.0*L],
            [  6.0*L, 2.0*L2, -6.0*L,  4.0*L2]])

    m = (density*area*element_length/420.0)*np.array([
            [ 156.0,   22.0*L,   54.0,   -13.0*L],
            [  22.0*L,  4.0*L2,  13.0*L,  -3.0*L2],
            [  54.0,   13.0*L,  156.0,   -22.0*L],
            [ -13.0*L, -3.0*L2, -22.0*L,   4.0*L2]])

    print(L)
    print(pd.DataFrame(k).to_markdown())
    exit()
    return k, m


def make_member(beam_length, density, elasticity_module, area, second_moment_of_area, number_of_elements=10):
    element_length = beam_length / number_of_elements

    element_stiffness, element_mass = make_element_stiffness(element_length, density, elasticity_module, area, second_moment_of_area)

    number_of_nodes = number_of_elements+1
    number_of_dofs = 2 * number_of_nodes

    system_stiffness = np.zeros([number_of_dofs, number_of_dofs])
    system_mass = np.zeros([number_of_dofs, number_of_dofs])
    print(number_of_dofs)
    for i in range(0, number_of_dofs-2, 2):
        print(i, i+4)

        system_stiffness[i:i+4, i:i+4] += element_stiffness
        system_mass[i:i+4, i:i+4] += element_mass

    print(pd.DataFrame(system_stiffness).to_markdown())
    return system_stiffness, system_mass, list(range(number_of_dofs))


def append_boundary_conditions(K, M, bcs):
    number_of_dofs = K.shape[0]  #  eq_nums[-1]
    print(number_of_dofs)

    fixed_dofs = []
    print("equation_numbers")
    for node, dofs in bcs:
        for dof in [int(i) for i in str(dofs)]:
            if node >= 0:
                equation_number = (node-1) * NUMBER_OF_NODAL_DOFS + dof-1
            else:
                equation_number = number_of_dofs + node * NUMBER_OF_NODAL_DOFS + dof-1
            fixed_dofs.append(equation_number)

    """
    np.delete(arr, np.s_[::2], 1)
    array([[ 2,  4],
           [ 6,  8],
           [10, 12]])
    np.delete(arr, [1,3,5], None)
    array([ 1,  3,  5,  7,  8,  9, 10, 11, 12])
    """

    ddd = [item for item in list(range(number_of_dofs)) if item not in fixed_dofs]
    return K[np.ix_(ddd, ddd)], M[np.ix_(ddd, ddd)], ddd


def fea_unsupported_beam(K, M, modes=5, remove_rigid_body_modes=True):
    eigenvalues, eigenvectors = eigh(system_stiffness, system_mass)

    # Eigenfrequencies (square roots of eigenvalues give natural frequencies)
    if remove_rigid_body_modes:
        natural_frequencies = np.sqrt(eigenvalues[NUMBER_OF_NODAL_DOFS:modes+NUMBER_OF_NODAL_DOFS])
    else:
        natural_frequencies = np.sqrt(eigenvalues)

    # Eigenmodes (eigenvectors)
    mode_shapes = eigenvectors

    # Display results
    #print("Natural Frequencies:", natural_frequencies)
    # print("Mode Shapes:\n", mode_shapes)

    return natural_frequencies


def fea_supported_beam(K_, M_, bcs, modes=5):
    K, M, eqs = append_boundary_conditions(K_, M_, bcs)

    print(K.shape, K_.shape)
    eigenvalues, eigenvectors = eigh(K, M)

    # Eigenfrequencies (square roots of eigenvalues give natural frequencies)
    natural_frequencies = np.sqrt(eigenvalues[:modes])

    # Eigenmodes (eigenvectors)
    mode_shapes = eigenvectors

    # Display results
    #print("Natural Frequencies:", natural_frequencies)
    # print("Mode Shapes:\n", mode_shapes)

    return natural_frequencies


def fea_cantilever_beam(K_, M_, bcs, modes=5):
    K, M, eqs = append_boundary_conditions(K_, M_, bcs)

    print(pd.DataFrame(K).to_markdown())


    eigenvalues, eigenvectors = eigh(K, M)

    # Eigenfrequencies (square roots of eigenvalues give natural frequencies)
    natural_frequencies = np.sqrt(eigenvalues[:modes])

    # Eigenmodes (eigenvectors)
    mode_shapes = eigenvectors

    # Display results
    #print("Natural Frequencies:", natural_frequencies)
    # print("Mode Shapes:\n", mode_shapes)

    return natural_frequencies

input_data = yaml.safe_load(open("analytic_models.yaml"))

NUMBER_OF_NODAL_DOFS = 2 #input_data["NUMBER_OF_NODAL_DOFS", 2]
ACCURACY = 1.0E-8 #input_data["ACCURACY", 1.0E-8]

number_of_modes = 4# input_data["number_of_modes"]
number_of_elements = 2 #input_data["number_of_elements"]
#points = input_data["points"]
#members = input_data["members"]




kwargs = dict(beam_length=1.0, density=8000.0, elasticity_module=2.0e11, **section_properties(bar=[0.005]))



eigen_frequencies_df = pd.DataFrame(index=range(1, number_of_modes+1))
eigen_frequencies_df["unsupported_beam"] = unsupported_beam(**kwargs, modes=number_of_modes)
eigen_frequencies_df["supported_beam"] = supported_beam(**kwargs, modes=number_of_modes)
eigen_frequencies_df["cantilever_beam"] = cantilever_beam(**kwargs, modes=number_of_modes)

system_stiffness, system_mass, equation_numbers = make_member(**kwargs, number_of_elements=number_of_elements)

eigen_frequencies_df["unsupported_fea"] = fea_unsupported_beam(system_stiffness, system_mass, number_of_modes)
eigen_frequencies_df["supported_fea"] = fea_supported_beam(system_stiffness, system_mass, [[1, 1], [-1, 1]], number_of_modes)
eigen_frequencies_df["cantilever_fea"] = fea_cantilever_beam(system_stiffness, system_mass, [[1, 12]], number_of_modes)


cols = [c for c in eigen_frequencies_df.columns if c.startswith("supported")]
print(eigen_frequencies_df[cols])

cols = [c for c in eigen_frequencies_df.columns if c.startswith("ca")]
print(eigen_frequencies_df[cols])


cols = [c for c in eigen_frequencies_df.columns if c.startswith("unsupported")]
print(eigen_frequencies_df[cols])


