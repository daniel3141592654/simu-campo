using LinearAlgebra
using Plots
using LaTeXStrings

pythonplot()

x_mesh(start, stop, len) = range(start, stop, len)' .* ones(len)

y_mesh(start, stop, len) = range(start, stop, len) .* ones(len)'

# Constante magnética

K = 1e-7;
I = 5.0; # A

# Radio de la espira

a = 50e-2; # m ???

# TAMAÑO DE TODO

sizeee = 30;

# plano coordenado:

x = x_mesh(-1, 1, sizeee);
y = y_mesh(-1, 1, sizeee);

# Se declara el iterador:

θ = range(0, 2 * pi, 500);
dθ = diff(θ)[end];

ds(i) = [cos(θ[i]); 0; -sin(θ[i])] .* a * dθ

ρ(x, y, i) = [x - a * sin(θ[i]); y; -a * cos(θ[i])]

b_component(x, y, n) = K * I * [cross(ds(i), ρ(x, y, i) / (norm(ρ(x, y, i))^3))[n] for i in eachindex(θ)] |> sum

superBX = reshape([b_component(x[i, j], y[i, j], 1) for j in 1:sizeee for i in 1:sizeee], (sizeee, sizeee));
superBY = reshape([b_component(x[i, j], y[i, j], 2) for j in 1:sizeee for i in 1:sizeee], (sizeee, sizeee));
superBZ = reshape([b_component(x[i, j], y[i, j], 3) for j in 1:sizeee for i in 1:sizeee], (sizeee, sizeee));

Mag_B = reshape([sqrt(superBX[i, j]^2 + superBY[i, j]^2) for j in 1:sizeee for i in 1:sizeee], (sizeee, sizeee));

surf = surface(x, y, Mag_B)
title!(L"Magnitud de $B$ / plano $xy$")

savefig(surf, "surf.pdf")
vector_b = contour(x, y, Mag_B, color=:plasma, clabels=false, fill=true)

quiver!(x, y, quiver=(superBX ./ Mag_B * 5e-2, superBY ./ Mag_B * 5e-2), markersize=1, color=:black)
title!(L"$\mathbf{B}$ vs. plano $xy$")

savefig(vector_b, "campo.pdf")
