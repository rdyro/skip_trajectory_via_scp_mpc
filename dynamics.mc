load(linearalgebra);
/* x = [v, th, h] */

/* density */
rho: rho0 * exp(-beta * x3);
beta: 1 / 8e3;

/* drag */
Cd0: 0.4;
K: 1.11;
Cd: Cd0 + K * u^2;
r: R + x3;

xn: [x1 + dt * -(Cd * rho * x1^2 * Area) / (2 * m) - g * sin(x2),
x2 + dt * (u * rho * x1 * Area) / (2 * m) + (x1 / r - g / x1) * cos(x2),
x3 + dt * x1 * sin(x2)];

A: jacobian(xn, [x1, x2, x3]);
B: jacobian(xn, [u]);
Acode: fortran(A);
Bcode: fortran(B);
