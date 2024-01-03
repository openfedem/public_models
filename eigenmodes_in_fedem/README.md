# Eigenmodes in FEDEM

In the context of nonlinear structural dynamics and the use of dynamic super elements (also known as Guyan reduced super element reduction), your description involves several distinct types of modes. Let's first clarify and verify the naming of these modes:

## Terms and Definitions

Free-Free Modes
: These are the vibrational modes of the entire finite element model (FEM) without any constraints or supports. Essentially, they represent the natural vibrations of the structure as if it were floating in space.

Constrained Modes
: These modes are the natural vibrations of the FEM with constraints or supports. The constraints alter the way the structure vibrates by restricting movement in certain directions or points.

Component Modes
: These are the modes of individual components or substructures within the larger FEM. In the context of super element reduction, these modes are considered during the process of reducing the complexity of the model.

Static Modes
: These are not vibrational modes but rather represent static deformations under certain load conditions. They are often included in the reduced model to capture the effects of static displacements.

Modes of the Reduced System (Free-Free and Constrained)
: After applying the Guyan reduction or a similar technique, the reduced system will have its own set of free-free and constrained modes. These are analogous to the modes of the full system but represent the behavior of the simplified, reduced model.


## Introduction

Introduction to Eigenmodes in Nonlinear Structural Dynamics Software Using Dynamic Super Elements

In this chapter, we delve into the intricate world of eigenmodes within the realm of nonlinear structural dynamics, particularly focusing on software that utilizes dynamic super elements, a method akin to Guyan reduced super element reduction. This sophisticated approach allows for a comprehensive understanding and analysis of complex structures through finite element modeling (FEM), providing insights into the vibrational characteristics under various configurations.

We commence with an exploration of the "Free-Free Modes" and "Constrained Modes" of the full FEM. The free-free modes offer a glimpse into the uninhibited vibrational behavior of the structure, assuming no constraints. In contrast, the constrained modes present the vibrational patterns when the structure is subjected to specific boundary conditions or supports, highlighting the impact of such constraints on the structural dynamics.

Next, we incorporate "Component Modes," which examine the vibrations of individual elements or substructures within the overall system. These modes are pivotal in understanding the localized dynamic behavior and are crucial during the reduction process. Alongside these, "Static Modes" are considered, representing the static response of the structure under various loading conditions. These modes are essential to capture the effects of static deformations on the overall dynamic behavior.

Finally, we transition to the "Free-Free and Constrained Modes of the Reduced System." This segment focuses on the behavior of the structure post-reduction, where the complexity of the full FEM is condensed into a more manageable form without significantly compromising the accuracy of dynamic predictions. These modes mirror their counterparts in the full system but are tailored to the characteristics of the reduced model.

In the following sections, we present a series of examples demonstrating these modes in action using our specialized software. These examples are meticulously designed to illustrate the nuances of each mode type and their significance in the comprehensive analysis and design of complex structural systems.

This introduction sets the stage for a detailed exploration of eigenmodes in nonlinear structural dynamics using dynamic super elements, providing a clear framework for understanding the subsequent examples and analyses.

## Classical beam 

\[ \omega \]

\begin{foo} 

\omega

\end{foo}


By the definition of parameter \( \beta \), we can calculate the eigen-frequency \( \omega \):

\[
\beta^4 = \frac{\omega^2 \cdot \overline{m}}{EI_y} \Rightarrow \omega^2 = \frac{EI_y \cdot \beta^4}{\overline{m}} \Rightarrow \omega = \sqrt{\frac{\beta^2 \cdot EI_y}{\overline{m}}}
\]
\(E37\)

Thus, inserting Eq. (36) into Eq. (37), the eigen-frequency \( \omega_n \) directly arises for each \( n \)-value.

\[
\omega_n = \frac{n^2 \cdot \pi^2}{L^2} \cdot \sqrt{\frac{EI_y}{\overline{m}}} \quad n = 1, 2, 3, \ldots
\]
\(E38\)

Therefore, the vibration mode-shape of the examined vertical pendulum arises by Eq. (27)—since previous inserting Eq. (35)—thus:


## Comparrisons

### Classi





### Free-free modes comparison
TODO: Compare the numerically assessed free-free mode eigenvalues and shape deviations from the theoretical shapes. 

Models:
* Supported beam
* Cantilever beam

Solutions (Forced 2D):
1. Exact solution (5 modes)
2. Numerical solution (2, 5 and 10 elements) (Timoshenko beam?)
3. Fedem system beams (2, 5 and 10 elements)
4. Fedem super-element with 10 beam elements (0, 2, 6 component modes)
5. Fedem super-element with 10x2x2 solid elements (0, 2, 6 component modes)


Eigenmode eigenvalues [angular frequency]

| Model                 | Mode 1 | Mode 2 | Mode 3 | Mode 4 | Mode 5 |
|:----------------------|:------:|:------:|:------:|:------:|:------:|
| Classical beam theory |   -    |   -    |   -    |   -    |   -    |
| Numerical 2 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 5 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 10 beams    |   -    |   -    |   -    |   -    |   -    |
| Fedem 2 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 5 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 10 system beams |   -    |   -    |   -    |   -    |   -    |
| Fedem SE, 2 modes     |   -    |   -    |   -    |   -    |   -    |
| Fedem SE, 5 modes     |   -    |   -    |   -    |   -    |   -    |
| Fedem SE, 10 modes    |   -    |   -    |   -    |   -    |   -    |

Eigenmodes, euclidian distance from the exact shape

| Model                 | Mode 1 | Mode 2 | Mode 3 | Mode 4 | Mode 5 |
|:----------------------|:------:|:------:|:------:|:------:|:------:|
| Classical beam theory |   -    |   -    |   -    |   -    |   -    |
| Numerical 2 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 5 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 10 beams    |   -    |   -    |   -    |   -    |   -    |
| Fedem 2 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 5 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 10 system beams |   -    |   -    |   -    |   -    |   -    |
| Fedem SE, 2 modes     |   -    |   -    |   -    |   -    |   -    |
| Fedem SE, 5 modes     |   -    |   -    |   -    |   -    |   -    |
| Fedem SE, 10 modes    |   -    |   -    |   -    |   -    |   -    |

![Classical beams](images/classical_beam.png)
![Numerical beams](images/classical_beam.png)
![Fedem beams](images/classical_beam.png)

### Constrained modes
Eigenmode shapes [angular frequency]

| Model                 | Mode 1 | Mode 2 | Mode 3 | Mode 4 | Mode 5 |
|:----------------------|:------:|:------:|:------:|:------:|:------:|
| Classical beam theory |   -    |   -    |   -    |   -    |   -    |
| Numerical 2 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 5 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 10 beams    |   -    |   -    |   -    |   -    |   -    |
| Fedem 2 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 5 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 10 system beams |   -    |   -    |   -    |   -    |   -    |

Eigenmodes, euclidian distance from the exact shape

| Model                 | Mode 1 | Mode 2 | Mode 3 | Mode 4 | Mode 5 |
|:----------------------|:------:|:------:|:------:|:------:|:------:|
| Classical beam theory |   -    |   -    |   -    |   -    |   -    |
| Numerical 2 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 5 beams     |   -    |   -    |   -    |   -    |   -    |
| Numerical 10 beams    |   -    |   -    |   -    |   -    |   -    |
| Fedem 2 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 5 system beams  |   -    |   -    |   -    |   -    |   -    |
| Fedem 10 system beams |   -    |   -    |   -    |   -    |   -    |


## Internal vibrations
add a Second cantilever towards the end of the primary cantilever. Tune the eigenvalues such that the first and second modes hence are on the first and second 

|    |            0 |        1 |         2 |   3 |   4 |   5 |            6 |        7 |         8 |                9 |              10 |              11 |               12 |           13 |              14 |
|---:|-------------:|---------:|----------:|----:|----:|----:|-------------:|---------:|----------:|-----------------:|----------------:|----------------:|-----------------:|-------------:|----------------:|
|  0 |  3.14159e+07 |     0    |     0     |   0 |   0 |   0 |  0           |     0    |     0     |     -3.14159e+07 |     0           |     0           |      0           |  0           |     0           |
|  1 |  0           |  9424.78 |  2356.19  |   0 |   0 |   0 |  0           |     0    |     0     |      0           | -9424.78        |  2356.19        |      0           |  0           |     0           |
|  2 |  0           |  2356.19 |   785.398 |   0 |   0 |   0 |  0           |     0    |     0     |      0           | -2356.19        |   392.699       |      0           |  0           |     0           |
|  3 |  0           |     0    |     0     |   0 |   0 |   0 |  0           |     0    |     0     |      0           |     0           |     0           |      0           |  0           |     0           |
|  4 |  0           |     0    |     0     |   0 |   0 |   0 |  0           |     0    |     0     |      0           |     0           |     0           |      0           |  0           |     0           |
|  5 |  0           |     0    |     0     |   0 |   0 |   0 |  0           |     0    |     0     |      0           |     0           |     0           |      0           |  0           |     0           |
|  6 |  0           |     0    |     0     |   0 |   0 |   0 |  3.14159e+07 |     0    |     0     |     -3.14159e+07 |     0           |     0           |      0           |  0           |     0           |
|  7 |  0           |     0    |     0     |   0 |   0 |   0 |  0           |  9424.78 | -2356.19  |      0           | -9424.78        | -2356.19        |      0           |  0           |     0           |
|  8 |  0           |     0    |     0     |   0 |   0 |   0 |  0           | -2356.19 |   785.398 |      0           |  2356.19        |   392.699       |      0           |  0           |     0           |
|  9 | -3.14159e+07 |     0    |     0     |   0 |   0 |   0 | -3.14159e+07 |     0    |     0     |      6.29826e+07 |     7.68545e-09 |     0           | -75398.2         | -3.84272e-09 | -9424.78        |
| 10 |  0           | -9424.78 | -2356.19  |   0 |   0 |   0 |  0           | -9424.78 |  2356.19  |      7.68545e-09 |     1.25683e+08 |     0           |     -3.84272e-09 | -6.28319e+07 |     5.77101e-13 |
| 11 |  0           |  2356.19 |   392.699 |   0 |   0 |   0 |  0           | -2356.19 |   392.699 |      0           |     0           |  4712.39        |   9424.78        | -5.77101e-13 |   785.398       |
| 12 |  0           |     0    |     0     |   0 |   0 |   0 |  0           |     0    |     0     | -75398.2         |    -3.84272e-09 |  9424.78        |  75398.2         |  3.84272e-09 |  9424.78        |
| 13 |  0           |     0    |     0     |   0 |   0 |   0 |  0           |     0    |     0     |     -3.84272e-09 |    -6.28319e+07 |    -5.77101e-13 |      3.84272e-09 |  6.28319e+07 |    -5.77101e-13 |
| 14 |  0           |     0    |     0     |   0 |   0 |   0 |  0           |     0    |     0     |  -9424.78        |     5.77101e-13 |   785.398       |   9424.78        | -5.77101e-13 |  1570.8         |
|    |         0 |                1 |               2 |              3 |            4 |        5 |                6 |               7 |              8 |                9 |           10 |              11 |
|---:|----------:|-----------------:|----------------:|---------------:|-------------:|---------:|-----------------:|----------------:|---------------:|-----------------:|-------------:|----------------:|
|  0 |   785.398 |      0           | -2356.19        |  392.699       |  0           |    0     |      0           | -2356.19        |  392.699       |      0           |  0           |     0           |
|  1 |     0     |      6.29826e+07 |     7.68545e-09 |    0           | -3.14159e+07 |    0     |      6.29826e+07 |     7.68545e-09 |    0           | -75398.2         | -3.84272e-09 | -9424.78        |
|  2 | -2356.19  |      7.68545e-09 |     1.25683e+08 |    0           |  0           | 2356.19  |      7.68545e-09 |     1.25683e+08 |    0           |     -3.84272e-09 | -6.28319e+07 |     5.77101e-13 |
|  3 |   392.699 |      0           |     0           | 4712.39        |  0           |  392.699 |      0           |     0           | 4712.39        |   9424.78        | -5.77101e-13 |   785.398       |
|  4 |     0     |     -3.14159e+07 |     0           |    0           |  3.14159e+07 |    0     |     -3.14159e+07 |     0           |    0           |      0           |  0           |     0           |
|  5 |     0     |      0           |  2356.19        |  392.699       |  0           |  785.398 |      0           |  2356.19        |  392.699       |      0           |  0           |     0           |
|  6 |     0     |      6.29826e+07 |     7.68545e-09 |    0           | -3.14159e+07 |    0     |      6.29826e+07 |     7.68545e-09 |    0           | -75398.2         | -3.84272e-09 | -9424.78        |
|  7 | -2356.19  |      7.68545e-09 |     1.25683e+08 |    0           |  0           | 2356.19  |      7.68545e-09 |     1.25683e+08 |    0           |     -3.84272e-09 | -6.28319e+07 |     5.77101e-13 |
|  8 |   392.699 |      0           |     0           | 4712.39        |  0           |  392.699 |      0           |     0           | 4712.39        |   9424.78        | -5.77101e-13 |   785.398       |
|  9 |     0     | -75398.2         |    -3.84272e-09 | 9424.78        |  0           |    0     | -75398.2         |    -3.84272e-09 | 9424.78        |  75398.2         |  3.84272e-09 |  9424.78        |
| 10 |     0     |     -3.84272e-09 |    -6.28319e+07 |   -5.77101e-13 |  0           |    0     |     -3.84272e-09 |    -6.28319e+07 |   -5.77101e-13 |      3.84272e-09 |  6.28319e+07 |    -5.77101e-13 |
| 11 |     0     |  -9424.78        |     5.77101e-13 |  785.398       |  0           |    0     |  -9424.78        |     5.77101e-13 |  785.398       |   9424.78        | -5.77101e-13 |  1570.8         |
