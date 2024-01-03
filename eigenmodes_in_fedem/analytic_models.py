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

COLORS = "bcgk"
MARKERS = "osx+*"
PRINT = True
PLOT = True

pi2 = np.pi ** 2.0

class BeamTheory:
    def __init__(self, length, density, elasticity_module, area, second_moment_of_area,
                 requested_modes=[0, 1, 2], accuracy=1.0E-8):

        self.length = length
        self.density = density
        self.elasticity_module = elasticity_module
        self.area = area
        self.second_moment_of_area = second_moment_of_area
        requested_modes = input_data.get("requested_modes", 3)
        if isinstance(requested_modes, int):
            self.requested_modes = range(1, requested_modes + 1)
        else:
            self.requested_modes = requested_modes
        self.accuracy = accuracy

    def beta_for_cantilever(self):
        """
        Finds multiple roots of the equation 1 + cos(beta * length) * cosh(beta * length) = 0 for beta.

        Parameters:
        length (float): The value of length in the equation.
        actual_modes (int): Number of roots to find.
        accuracy (float): The desired accuracy of the solution.

        Returns:
        list: A list of solutions for beta.
        """
        # print("beta_for_cantilever")
        roots = []

        # Function to find the root near a given guess
        def find_root_near(guess):
            root, = fsolve(lambda beta: 1.0 + np.cos(beta * self.length) * np.cosh(beta * self.length), guess, xtol=self.accuracy)
            return root

        # Initial guess for the first root
        guess = 0.1 / self.length
        for i in self.requested_modes:
            root = find_root_near(guess)
            # Avoid duplicates
            if not any(np.isclose(root, r, atol=self.accuracy) for r in roots):
                roots.append(root)
            # Update the guess for the next root
            guess = 0.1 * (root + np.pi / self.length)

        # for i in range(actual_modes):
        #    n = i + 1
        #    print((2*n+1)*np.pi/length)
        #    roots[i] = (2*n+1)*np.pi/length
        # print(roots)

        roots = [3.516, 22.03, 61.70, 120.9]  # , 17.2787]

        return roots

    def beta_for_free_beam(self):
        """
        cosh(kL)cos(kL)=1
        """
        # print("beta_for_free_beam")
        roots = []

        # Function to find the root near a given guess
        def find_root_near(guess):
            root, = fsolve(lambda beta: 1.0 - np.cos(beta * self.length) * np.cosh(beta * self.length), guess, xtol=accuracy)
            return root

        # Initial guess for the first root
        guess = 0.1 / self.length
        for i in range(self.requested_modes):
            root = find_root_near(guess)
            # Avoid duplicates
            if not any(np.isclose(root, r, atol=self.accuracy) for r in roots):
                roots.append(root)
            # Update the guess for the next root
            guess = root + np.pi / self.length
        # print(roots)

        roots = [4.7300 / self.length, 7.8532 / self.length, 10.9956 / self.length, 14.1371 / self.length]

        return roots

    def alpha_free_beam(self):
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
        for mode in self.requested_modes:
            # Initial guess for alpha, can be adjusted for better convergence
            alpha_guess = (mode - 0.5) * np.pi
            alpha = fsolve(transcendental_eq, alpha_guess, args=(mode))
            alpha_values.append(alpha[0])

        alpha_values = [4.7300, 7.8532, 10.9956, 14.1371, 17.2787]

        return alpha_values

    def free_beam(self):

        # Example: Calculate the first 5 alpha values for a free-free beam
        alpha_values = self.alpha_free_beam()

        unit_mass = self.density * self.area

        eigenvalues = []
        for alpha in [alpha_values[i-1] for i in self.requested_modes]:
            eigenvalues.append(alpha ** 2.0 *
                               np.sqrt(self.elasticity_module * self.second_moment_of_area / (
                                           unit_mass * self.length ** 4.0)))

        return eigenvalues

    def supported_beam(self):
        unit_mass = self.density * self.area

        eigenvalues = []
        for mode_number in self.requested_modes:
            eigenvalues.append(mode_number ** 2.0 * pi2 / self.length ** 2.0 *
                               np.sqrt(self.elasticity_module * self.second_moment_of_area / unit_mass))

        return eigenvalues

    def cantilever_beam(self):

        beta_values = self.beta_for_cantilever()
        unit_mass = self.density * self.area

        eigenvalues = []
        for beta in [beta_values[i-1] for i in self.requested_modes]:
            eigenvalues.append(beta / self.length ** 2.0 *
                               np.sqrt(self.elasticity_module * self.second_moment_of_area / unit_mass))

        return eigenvalues


class Beam2D:
    def __init__(self, material_properties, cross_section_parameters, points,
                 members, boundary_conditions, **kwargs):

        self.material_properties = material_properties
        self.cross_section_parameters = cross_section_parameters
        self.points = points
        self.members = members
        self.boundary_conditions = boundary_conditions
        make_system = kwargs.get("make_system", False)
        make_plot = kwargs.get("make_plot", False)

        self.dofs = {}
        self.member_nodes = {}
        self.nodes = {k: np.array(v) for k, v in self.points.items()}

        self.stiffness = {}
        self.mass = {}

        if make_system:
            self.make_system()

        if make_plot:
            self.make_plot()

    @staticmethod
    def make_cross_section_properties(*vargs):
        # A, I
        if vargs[0] == "bar":
            radius = vargs[1]
            return np.pi*radius**2.0, np.pi*radius**4.0/4.0

    @staticmethod
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

    @staticmethod
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

    def append_member_to_system(self, member_name, point_start, point_end, material_property, cross_section_property):
        """
        append_member_to_system(
            self.member_nodes[member_name], self.nodes[point_start], self.nodes[point_end],
            material_properties[material],
            make_cross_section_properties(*cross_section_parameters[cross_section]),
            self.dofs, system_stiffness, system_mass)
        """
        member_node_names = self.member_nodes[member_name]
        xy_start = np.array(self.nodes[point_start])
        xy_end = np.array(self.nodes[point_end])

        dx = xy_end - xy_start
        beam_length = np.linalg.norm(dx)
        beam_orientation = np.arctan2(*np.flip(dx))

        density, elasticity_module = material_property
        area, second_moment_of_area = cross_section_property

        element_length = beam_length / (len(member_node_names)-1)

        _element_stiffness, _element_mass = self.make_element_stiffness(element_length, density, elasticity_module, area, second_moment_of_area)

        transformation_matrix = self.make_transformation_matrix(beam_orientation)

        element_stiffness = transformation_matrix.T @ _element_stiffness @ transformation_matrix
        element_mass = transformation_matrix.T @ _element_mass @ transformation_matrix
        for node_start, node_end in zip(member_node_names[:-1], member_node_names[1:]):
            eqs = np.concatenate([self.dofs[node_start], self.dofs[node_end]])-1

            #print(node_start, node_end, eqs)

            self.system_stiffness[np.ix_(eqs, eqs)] += element_stiffness
            self.system_mass[np.ix_(eqs, eqs)] += element_mass

    def make_system_nodes_and_dofs(self, member_name, point_start, point_end, number_of_elements):
        xy_start = self.nodes[point_start]
        xy_end = self.nodes[point_end]
        self.member_nodes[member_name] = []
        member_node_names = self.member_nodes[member_name]

        dx = (xy_end-xy_start)/number_of_elements

        if point_start not in member_node_names:
            member_node_names.append(point_start)
            if point_start not in self.dofs:
                self.nodes[point_start] = xy_start
                self.dofs[point_start] = np.array([1, 2, 3]) + len(self.dofs) * 3

        for i in range(1, number_of_elements):
            node_name = f"{member_name}.{i}"
            member_node_names.append(node_name)
            if node_name not in self.dofs:
                self.nodes[node_name] = xy_start + dx * i
                self.dofs[node_name] = np.array([1, 2, 3]) + len(self.dofs) * 3

        if point_end not in member_node_names:
            member_node_names.append(point_end)
            if point_end not in self.dofs:
                self.nodes[point_end] = xy_end
                self.dofs[point_end] = np.array([1, 2, 3]) + len(self.dofs) * 3

    def make_plot(self):

        fig = plt.figure(figsize=(8, 6))  # Width and height of the figure in inches

        # Adding an axis
        ax = fig.add_subplot(111)
        ax.set_aspect(1)

        for c, (member_name, member_node_names) in enumerate(self.member_nodes.items()):
            xys = np.array([xy for name, xy in self.nodes.items() if name in member_node_names])

            ax.plot(xys[:, 0], xys[:, 1], f"{COLORS[c]}.-", label=member_name, linewidth=3)
            for name in member_node_names:
                # Adding an annotation
                ax.annotate(name, self.nodes[name], textcoords="offset points", xytext=(0, 10), ha='center')

        ax.grid()
        #plt.show()
        plt.savefig("plot.png")

        self.ax = ax

    def make_unconstrained_system_matrices(self):
        for member_name, (_, _, point_start, point_end, number_of_elements) in self.members.items():
            self.make_system_nodes_and_dofs(member_name, point_start, point_end, number_of_elements)

        if PRINT:
            print("Member self.nodes:")
            for member_name, member_i_nodes in self.member_nodes.items():
                print(f"  - Member: {member_name}")
                for k in member_i_nodes:
                    print(f"    - {k:5}")

            print("self.nodes:")
            for node_name, (x, y) in self.nodes.items():
                print(f"    - {node_name:5}{x:8.3f}{y:8.3f}")

            print("\nself.dofs, (Degrees of freedom, aka equation numbers) for free system:")
            for node_name, node_dofs in self.dofs.items():
                print(f"  - {node_name:5}{node_dofs[0]:4}{node_dofs[1]:4}{node_dofs[2]:4}")

        size = max(max(equations) for equations in self.dofs.values())
        if PRINT:
            print(f"\nSystem size: {size}")

        self.system_stiffness = lil_matrix((size, size))
        self.system_mass = lil_matrix((size, size))
        for member_name, (material, cross_section, point_start, point_end, _) in self.members.items():
            self.append_member_to_system(member_name, point_start, point_end,
                self.material_properties[material],
                self.make_cross_section_properties(*self.cross_section_parameters[cross_section]))

        self.system_stiffness = self.system_stiffness.tocsr()
        self.system_mass = self.system_mass.tocsr()
        print(f"Stiffness shape {self.system_stiffness.shape}")
        print(f"Mass shape {self.system_mass.shape}")

    def make_constrained_system_matrices(self, bcs):
        for bc in bcs:
            i = bc.rfind(".")
            node, node_dofs = bc[:i], bc[i+1:]
            for node_dof in [int(i) for i in node_dofs]:
                self.dofs[node][node_dof-1] = - self.dofs[node][node_dof-1]

        self.equations = [num-1 for sublist in self.dofs.values() for num in sublist if num > 0]

        x = 0
        for k, vs in self.dofs.items():
            for i, v in enumerate(vs):
                if v < 0:
                    self.dofs[k][i] = - 1
                    x += 1
                else:
                    self.dofs[k][i] = v - x - 1



class Solvers:
    def __init__(self, dofs, boundary_condition, stiffness, mass, **kwargs):
        self.dofs = copy.deepcopy(dofs)

        self.equations = None
        self.make_equations(boundary_condition)

        self.stiffness = stiffness[np.ix_(self.equations, self.equations)]
        self.mass = mass[np.ix_(self.equations, self.equations)]

        solve_modes = kwargs.get("solve_modes", False)
        add_modes_to_plot = kwargs.get("add_modes_to_plot", False)

        self.eigenvalues = None
        self.eigenvectors = None

        if solve_modes:
            self.calculate_eigenvalues()

        if add_modes_to_plot:
            self.add_modes_to_plot()


    def DEPRICATED_make_models(self):
        if PRINT:
            print("\nModels:")

        for model_name, boundary_condition in self.boundary_conditions.items():
            self.models[model_name] = self.append_boundary_conditions(boundary_condition)

            if PRINT:
                print(f"  - Model self.dofs for {model_name}")
                for node_name, node_dofs in self.dofs.items():
                    print(f"    - {node_name:5}{node_dofs[0]:4}{node_dofs[1]:4}{node_dofs[2]:4}")

        if PRINT > 2:
            print(equations)
            print(K.shape, _K.shape)
            print(M.shape, _M.shape)
            print(pd.DataFrame(K.toarray()).to_markdown())
            print(pd.DataFrame(_K.toarray()).to_markdown())

    def solve_modes(self, requested_modes=None, lower_cutoff=0.01):

        # TODO: replace with spars matrix supported solver
        eigenvalues, eigenvectors = eigh(self.stiffness.toarray(), self.mass.toarray())

        if requested_modes is None:
            requested_modes = range(len(eigenvalues))


        eigenvalues = np.sqrt(eigenvalues)
        skip = np.argmax(eigenvalues > lower_cutoff)
        eigenvalues = eigenvalues[skip:]
        eigenvectors = eigenvectors[:, skip:]

        self.eigenvalues = eigenvalues[np.ix_(requested_modes)]
        self.eigenvectors = eigenvectors[:, np.ix_(requested_modes)]

    def add_modes_to_plot(self):
        for c, (member_name, member_node_names) in enumerate(self.member_nodes.items()):
            gg = np.zeros((len(member_node_names), 2))

            for im, m in enumerate((1, 2)):
                emi = em[:, m] * 0.1
                emi = em[:, m] * 0.1

                for ii, k in enumerate(member_node_names):

                    i, j = self.dofs[k][0:2]
                    x, y = self.nodes[k][:2]

                    dx = emi[i] if i >= 0 else 0.0
                    dy = emi[j] if j >= 0 else 0.0

                    gg[ii, :] = [x + dx, y + dy]

                lab = f"Mode {self.requested_modes[c]}" if im==0 else ""
                ax.plot(gg[:, 0], gg[:, 1], COLORS[c], alpha=0.4, label=lab, marker=MARKERS[m])

        plt.legend()
        plt.savefig("plot.png")


class Solutions(Beam2D):
    def __init__(self, lower_cutoff=1.0e-2, **kwargs):
        super().__init__(**kwargs)
        self.lower_cutoff = lower_cutoff
        calculate_eigenmodes = kwargs.get("calculate_eigenmodes", False)

        self.solutions = {}
        for name, boundary_condition in self.boundary_conditions.items():
            self.solutions[name] = Solvers(self.dofs, boundary_condition, self.system_stiffness, self.system_mass)

        if calculate_eigenmodes:
            self.calculate_eigenmodes()

    def calculate_eigenmodes(self):
        for name, solution in self.solutions.items():
            solution.solve_modes()


    """
    def append_boundary_conditions(self, bcs, inplace=False):
        if inplace:
            dofs = self.dofs
        else:
            dofs = copy.deepcopy(self.dofs)

        fixed_dofs = []
        for bc in bcs:
            i = bc.rfind(".")
            node, node_dofs = bc[:i], bc[i+1:]
            for node_dof in [int(i) for i in node_dofs]:  #local dof:  1, 2 or 3
                fixed_dofs.append(dofs[node][node_dof-1] )  # Minus 1 as dof starts on 1
                dofs[node][node_dof-1] = - dofs[node][node_dof-1]

        equations = [num-1 for sublist in dofs.values() for num in sublist if num > 0]
        _K = self.system_stiffness[np.ix_(equations, equations)]
        _M = self.system_mass[np.ix_(equations, equations)]

        x = 0
        for k, vs in dofs.items():
            for i, v in enumerate(vs):
                if v < 0:
                    dofs[k][i] = - 1
                    x += 1
                else:
                    dofs[k][i] = v - x - 1

    def make_models(self):
        if PRINT:
            print("\nModels:")

        for model_name, boundary_condition in self.boundary_conditions.items():
            self.models[model_name] = self.append_boundary_conditions(boundary_condition)

            if PRINT:
                print(f"  - Model self.dofs for {model_name}")
                for node_name, node_dofs in self.dofs.items():
                    print(f"    - {node_name:5}{node_dofs[0]:4}{node_dofs[1]:4}{node_dofs[2]:4}")

        if PRINT > 2:
            print(equations)
            print(K.shape, _K.shape)
            print(M.shape, _M.shape)
            print(pd.DataFrame(K.toarray()).to_markdown())
            print(pd.DataFrame(_K.toarray()).to_markdown())

    def solve_modes(self):
        eigenvalues, eigenvectors = eigh(self.stiffness.toarray(), self.mass.toarray())

        _eigenvalues = np.sqrt(eigenvalues)
        skip_cutoff = np.argmax(_eigenvalues > lower_cutoff)


        #np.ix_(self.requested_modes)
        #self.eigenvalues = eigenvalues[skip_cutoff:modes + skip_cutoff]
        #self.eigenvectors = eigenvectors[:, skip_cutoff:modes + skip_cutoff]
        self.eigenvalues = eigenvalues[np.ix_(self.requested_modes)]
        self.eigenvectors = eigenvectors[:, np.ix_(self.requested_modes)]

    def append_modes_to_plot(self):

        for c, (member_name, member_node_names) in enumerate(self.member_nodes.items()):
            gg = np.zeros((len(member_node_names), 2))

            for im, m in enumerate((1, 2)):
                emi = em[:, m] * 0.1
                emi = em[:, m] * 0.1

                for ii, k in enumerate(member_node_names):

                    i, j = self.dofs[k][0:2]
                    x, y = self.nodes[k][:2]

                    dx = emi[i] if i >= 0 else 0.0
                    dy = emi[j] if j >= 0 else 0.0

                    gg[ii, :] = [x + dx, y + dy]

                lab = f"Mode {self.requested_modes[c]}" if im==0 else ""
                ax.plot(gg[:, 0], gg[:, 1], COLORS[c], alpha=0.4, label=lab, marker=MARKERS[m])

        plt.legend()
        plt.savefig("plot.png")
    """


if __name__ == "__main__":
    print("Modes\n")
    input_data = yaml.safe_load(open("analytic_models.yaml"))

    # Beam analytics
    beam_analytics = BeamTheory(**input_data["Beam Analytics"])

    eigen_frequencies_df = pd.DataFrame(index=beam_analytics.requested_modes)
    eigen_frequencies_df["free_beam"] = beam_analytics.free_beam()
    eigen_frequencies_df["supported_beam"] = beam_analytics.supported_beam()
    eigen_frequencies_df["cantilever_beam"] = beam_analytics.cantilever_beam()

    # Beam FEA
    model_B = Solutions(**input_data["Beam FEA"]["B"], make_system=True, make_plot=True)#, solve_modes=True, add_modes_to_plot=True)
    model_B.calculate_eigenmodes()
    model_B.add_modes_to_plot()





    #stiffness, mass, nodes, member_nodes, dofs, modes





    #input_model_B = input_data["models"]["B"]
    #model_B, ax, self.nodes, self.member_nodes = make_fea_model(**input_model_B)


    #_, _, equations, model_dofs = model_B["supported_B"]

    #ev, em = solve_modes(*model_B["supported_B"][:2], number_of_modes)
    #print(ev)

    #append_modes_to_plot(ax, self.nodes, equations, model_dofs, em, self.member_nodes)









    # Printing, output and reporting
    cols = [c for c in eigen_frequencies_df.columns if c.startswith("su")]
    print(eigen_frequencies_df[cols])

    cols = [c for c in eigen_frequencies_df.columns if c.startswith("ca")]
    print(eigen_frequencies_df[cols])

    cols = [c for c in eigen_frequencies_df.columns if c.startswith("fr")]
    print(eigen_frequencies_df[cols])





"""
    input_data = yaml.safe_load(open("analytic_models.yaml"))
    accuracy = input_data.get("accuracy", 1.0E-8)

    number_of_modes = input_data["number_of_modes"]

    input_model_B = input_data["models"]["B"]
    model_B, ax, self.nodes, self.member_nodes = make_fea_model(**input_model_B)

    _, _, equations, model_dofs = model_B["supported_B"]

    ev, em = solve_modes(*model_B["supported_B"][:2], number_of_modes)
    print(ev)

    append_modes_to_plot(ax, self.nodes, equations, model_dofs, em, self.member_nodes)

    exit()




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


    
    
    #import numpy as np
    #from scipy.sparse.linalg import svds
    
    # Create a sparse matrix
    #A = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    #sparse_A = sps.csr_matrix(A)
    
    # Compute the rank of the sparse matrix
    #u, s, vt = svds(sparse_A)
    #rank = len([x for x in s if abs(x) > 1e-12])
    
    #print(f"The rank of the sparse matrix is {rank}.")
    
    






    kwargs = dict(beam_length=1.0, density=8000.0, elasticity_module=2.0e11, **section_properties(bar=[0.005]))
    #eigen_frequencies_df = pd.DataFrame(index=range(1, number_of_modes+1))
    eigen_frequencies_df["free_beam"] = free_beam(**kwargs, modes=number_of_modes)
    eigen_frequencies_df["supported_beam"] = supported_beam(**kwargs, modes=number_of_modes)
    eigen_frequencies_df["cantilever_beam"] = cantilever_beam(**kwargs, modes=number_of_modes)

    cols = [c for c in eigen_frequencies_df.columns if c.startswith("supported")]
    print(eigen_frequencies_df[cols])

    cols = [c for c in eigen_frequencies_df.columns if c.startswith("ca")]
    print(eigen_frequencies_df[cols])

    cols = [c for c in eigen_frequencies_df.columns if c.startswith("free")]
    print(eigen_frequencies_df[cols])


"""