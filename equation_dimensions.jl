### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 53a4d7b0-9aa8-11ef-2a21-01318b24dbbe
md"""
# _Report on the properties of the COMBAT model_

The equation of the COMBAT model as presented by Tran *et al.* (2022) are the following:

$$\begin{aligned} \frac{\text {d}B_x}{dt}& = \frac{k_{f}}{Vn_{A}}(n-x+1)AB_{x-1} - k_{r}xB_x - \frac{k_{f}}{Vn_{A}}(n-x)AB_x \\&\quad + k_{r}(x+1)B_{x+1} + \rho _x -r_xB_x \frac{C-\sum _{j=0}^{n}B_j}{C} -d_xB_x \\ \frac{\text {d}A}{dt} &= - \frac{k_{f}}{Vn_{A}}(A\cdot T +\sum _{x=0}^{n-1}(n-x)AB_x) + k_{r}\left(A_T+\sum _{x=1}^{n}xB_x\right) \\ \frac{\text {d}T}{dt} &= - \frac{k_{f}}{Vn_{A}}A\cdot T + k_{r}A_T + \sum _{x=0}^{n}d_x(n-x)B_x \\ \frac{\text {d}A_T}{dt} &= \frac{k_{f}}{Vn_{A}}A\cdot T - k_{r}A_T + \sum _{x=0}^{n}d_xxB_x \\ \rho _x &= 2\sum _{i=x}^{n}f_{i,x}r_i B_i \frac{C-\sum _{j=0}^{n}B_j}{C}, \end{aligned}$$

where the symbols mean the following:

| Symbol | Type | Meaning |
|:---------- | ---------- |:------------:|
| $B_x$    | Model variable | Bacteria with $x$ bound targets ($0\leq x\leq n$)|
| $n$    | Parameter  | Number of targets per bacterium |
| $k_f$  | Parameter | Binding rate of antibiotic to target |
| $k_r$  | Parameter | Unbinding rate of antibiotic to target |
| $V$  | Parameter | Average intracellular volume |
| $n_A$  | Constant | Avodargo number ($\approx 6\cdot 10^{23} \;\mathrm{mol}^{-1}$)|
| $C$  | Parameter | Carrying capacity|
| $A$  | Model variable | Drug concentration|
| $T$  | Model variable | Concentration of extracellular free targets|
| $A_T$  | Model variable | Concentration of extracellular free targets bound to antibotic|
| $\rho_x$  | Function of model variables | Total rate of replication of bacteria with $x$ bound targets|
| $r_x$  | Function of model variables | Replication rate of bacteria with $x$ bound targets|
| $d_x$  | Function of model variables | Death rate of bacteria with $x$ bound targets|
| $f_{i,x}$  | Function of model variables | Hypergeometric distribution function|
"""

# ╔═╡ 4dd3d9bb-d821-4c40-b243-8ebb5b0955ba


# ╔═╡ Cell order:
# ╠═53a4d7b0-9aa8-11ef-2a21-01318b24dbbe
# ╠═4dd3d9bb-d821-4c40-b243-8ebb5b0955ba
