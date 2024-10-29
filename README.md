<header>

# Laplace-2D-Jacobi
This is a report of major exercise of IT4110

</header>

### 2D Electric Potential Solver

Member: (Full Name & Student ID)
- Đỗ Hoàng Minh Hiếu 20225837
- Trịnh Minh Đạt 20225701
- Trần Phúc Sơn 20225666
- Hoàng Trường Giang 20225710
- Ninh Lê Gia Bảo 20225693
- Đặng Huy Hoàng 20225843

Problem & Solution:
- The Laplace equation is a significant equation with applications across various fields in Science and Engineering. These applications are not limited to electrostatics but also extend to fluid dynamics and steady-state heat conduction. It is a second-order partial differential equation (PDE) with an elliptic nature. Solutions to the Laplace equation are often known as harmonic functions, as seen in several fields, and are particularly useful in Science and Engineering, as previously mentioned. The Laplace equation has both analytical and numerical solutions. Numerical solutions for the Laplace equation are found using various methods applicable to many linear PDEs, including the finite difference method (FDM), which is the focus of this study. The approach taken by our group to solve the Laplace equation centers around FDM. Here, we only use the numerical solutions of the Laplace equation, setting up boundary conditions for the equation, which can be derived from an existing equation of the electric potential function using the analytical method. The numerical solutions of the Laplace equation are compared with analytical solutions, showing similar results, allowing the use of common boundary conditions. The finite difference method (FDM) has been implemented through the Jacobi iterative method. This approach is taken to evaluate the efficiency of the Jacobi iterative method in calculating electric potential on a 2D plane while yielding accurate results.

Algorithm:
- In the electric potential equation U=V, the finite difference method is used to divide the shape into a grid of points, with the number of grid points in length being a and in width being b. This gives us a×b points for which the electric potential needs to be calculated.
  <p align="center">
    <img src=https://github.com/user-attachments/assets/3966e88a-b240-4ff7-a099-97c31dbe9f42 alt=celebrate width=300 align=center>
  </p>
- Using the Taylor-Maclaurin expansion, we can derive the following approximation formula:
  <p align="center">
    <img src=https://github.com/user-attachments/assets/dca9ae7a-3d54-44d6-8155-29d7cda6d09d alt=celebrate width=300 align=center>
  </p>
- In practice, when using the finite difference method along the x-axis, the electric potential with a second-order derivative with respect to x approximates the formula:
  <p align="center">
    <img src=https://github.com/user-attachments/assets/db93f133-6387-4aef-81e0-48a7ba4cfc29 alt=celebrate width=300 align=center>
  </p>
- Similarly, for the second-order derivative with respect to y, we derive the following numerical solution formula:
  <p align="center">
    <img src=https://github.com/user-attachments/assets/5ba8b881-175e-4099-a177-6c4f6213a7ed alt=celebrate width=300 align=center>
  </p>
- In the coordinate system, we have:
  <p align="center">
    <img src=https://github.com/user-attachments/assets/d301ac47-80bf-4c87-ab32-637b631341ca alt=celebrate width=300 align=center>
  </p>
- In practice, aside from using the numerical method, we can also use the analytical method to calculate the electric potential. Essentially, the boundaries and potentials calculated at points in both methods are approximately the same.
- In the analytical solution method, the formula u(x,y) is predefined. We know the conditions of the four boundaries of the plane to use in the code.
  <p align="center">
    <img src=https://github.com/user-attachments/assets/6003198c-cad7-431d-bb0f-4984c869f57c alt=celebrate width=300 align=center>
  </p>
- With ω = 1 we have
  1. u(x,b) = sin(x)/sin(a)
  2. u(a,y) = sinh(y)/sinh(b)
  3. u(x,0) = 0
  4. u(0,y) = 0
- We will use these four boundaries to set boundary conditions when creating matrix U in the next code section. In the code, set the numerator of the numerical solution to 1, thus obtaining a Courant coefficient for the numerical solution as follows:
  <p align="center">
    <img src=https://github.com/user-attachments/assets/aba659c4-66b2-4678-a28e-3878df030d58 alt=celebrate width=300 align=center>
  </p>
- Next, use the Jacobi iteration method to calculate the next point. Definition: The Jacobi iteration formula is given by:
  1. Iteration process: Initialize the starting value
  2. Update the value
  3. Check for convergence

Result:
- Input:
  1. Length of Retangle
  2. Width of Retangle
  3. Number of grid points along the x-axis
  4. Number of grid points along the y-axis
  5. Error threshold
  <p align="center">
    <img src=https://github.com/user-attachments/assets/c9e68376-bd2c-433a-83e1-581378883ac2 alt=celebrate width=300 align=center>
  </p>

- Output:
  1. Number of iterations
  2. Final error
  <p align="center">
    <img src=https://github.com/user-attachments/assets/87a6c404-c9f1-43d5-a7ad-1a25880e9c39 alt=celebrate width=300 align=center>
  </p>






