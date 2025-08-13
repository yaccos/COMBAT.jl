### A Pluto.jl notebook ###
# v0.20.15

using Markdown
using InteractiveUtils

# ╔═╡ 2106bf0c-1580-4d22-81fe-6a74f4102d97
import PlutoUI.LocalResource

# ╔═╡ 53a4d7b0-9aa8-11ef-2a21-01318b24dbbe
md"""
# _Report on the properties of the COMBAT model_

## The equations and their symbols
The equations of the COMBAT model as presented by Tran *et al.* (2022) are the following:

$$\begin{aligned} \frac{\text {d}B_x}{\mathrm{d}t}& = \frac{k_{f}}{Vn_{A}}(n-x+1)AB_{x-1} - k_{r}xB_x - \frac{k_{f}}{Vn_{A}}(n-x)AB_x \\&\quad + k_{r}(x+1)B_{x+1} + \rho _x -r_xB_x \frac{C-\sum _{j=0}^{n}B_j}{C} -d_xB_x 
\\ \frac{\text {d}A}{\mathrm{d}t} &= - \frac{k_{f}}{Vn_{A}}\left(A\cdot T +\sum _{x=0}^{n-1}(n-x)AB_x\right) + k_{r}\left(A_T+\sum _{x=1}^{n}xB_x\right) 
\\ \frac{\text {d}T}{\mathrm{d}t} &= - \frac{k_{f}}{Vn_{A}}A\cdot T + k_{r}A_T + \sum _{x=0}^{n}d_x(n-x)B_x 
\\ \frac{\text {d}A_T}{\mathrm{d}t} &= \frac{k_{f}}{Vn_{A}}A\cdot T - k_{r}A_T + \sum _{x=0}^{n}d_xxB_x 
\\ \rho_x &= 2\sum _{i=x}^{n}f_{i,x}r_i B_i \frac{C-\sum _{j=0}^{n}B_j}{C}. \end{aligned}$$

Since $x$ varies from $0$ to $n$ inclusive, these equations constitute a system of of $n+4$ variables where the symbols mean the following:

| Symbol | Type | Meaning |
|:---------- | ---------- |:------------:|
| $B_x$    | Model variable | Bacteria with $x$ bound targets ($0\leq x\leq n$)|
| $n$    | Parameter  | Number of targets per bacterium |
| $k_f$  | Parameter | Binding rate of antibiotic to target |
| $k_r$  | Parameter | Unbinding rate of antibiotic to target |
| $V$  | Parameter | Average intracellular volume |
| $n_A$  | Constant | Avogadro number ($\approx 6\cdot 10^{23} \;\mathrm{mol}^{-1}$)|
| $C$  | Parameter | Carrying capacity|
| $A$  | Model variable | Drug concentration|
| $T$  | Model variable | Concentration of extracellular free targets|
| $A_T$  | Model variable | Concentration of extracellular free targets bound to antibotic|
| $\rho_x$  | Function of model variables | Total rate of replication of bacteria with $x$ bound targets|
| $r_x$  | Function of model variables | Replication rate of bacteria with $x$ bound targets|
| $d_x$  | Function of model variables | Death rate of bacteria with $x$ bound targets|
| $f_{i,x}$  | Function of model variables | Hypergeometric distribution function|

Unless states otherwise, we will only consider the original model which does not have any input file and $A$ will therefore vary according to these equations.

## How the parameters are implemented and interpreted

For initializing the system, the vCOMBAT C implementation (after some modifications) accepts the following parameters with the corresponding units of measurements:

| Symbol in paper | Command line option| Meaning | Unit | Default value |
|:---------- |:---------- | ---------- | ---------- |:------------:|
| $A$    | `-d`| Initial antibiotic concentration | mg/L| $10000$|
| $k_f$    | `-A`| Binding rate of antibiotic to target | L/mol/s| $4450$|
| $k_r$    | `-D`| Dissociation rate of antibiotic from target | 1/s| $0.0023$|
| $B_0$    | `-p`| Initial bacterial population | Cells | $10^6$|
| (Missing)    | `-t`| [total time]:[time between timepoints] | s | $360000:3600$ ($100$ hours with $1$ hour steps)|
| $n$    | `-n`| Number of target molecules | Molecules | $100$|
| $r_T$    | `-r`| Replication threshold | | $\frac{n}{2}$|
| $R_0$    | `-R`| Maximum replication rate | 1/s | $8.34 \cdot 10^{-6}$|
| $k_T$    | `-k`| Killing threshold |  | $\frac{n}{2}$|
| $V$    | `-V`| Intracellular volume | L/Cell | $10^{-15}$|
| $D_0$    | `-K`| Maxium death rate | 1/s | $1.39\cdot 10^{-5}$|
| W    | `-M`| Molecular weight of drug | g/mol | $555.5$|
| $C$    | `-C`| Carrying capacity | Cells | $10^9$|

Furthermore the paper and C implementations operate with the replication and killing rates as:


$r_x = \begin{cases} R_0\left(1-\frac{x}{r_T}\right) & x\leq r_T  \\ 0 & x>r_T \end{cases}$

$d_x = \begin{cases} 0 & x < k_T  \\ D_0 & x\geq k_T \end{cases}$

This is, the graphs look something like this:
"""

# ╔═╡ 9b9b544a-e7b7-4907-a271-d862f99a2de0
LocalResource("replication_killing_rates.png")

# ╔═╡ 6fb53591-d97f-429c-b879-3b7030e65719
md"""

## Internal conversion from macroscopic to molecular units

Internally, the vCOMBAT implementation converts the drug concentration measured in mass per volume to the number of drug molecules per cell before the simulation is run:
$$A_{\text{Molecular}}=\frac{A_{\text{Macroscopic}}\cdot V\cdot n_A}{W}$$
We verify that this actually yields the number of molecules per cell by checking the units:
$$\left[A_{\text{Molecular}}\right]=\left[\frac{A_{\text{Macroscopic}}\cdot V\cdot n_A}{W}\right]=\frac{\left[A_{\text{Macroscopic}}\right]\cdot \left[V\right]\cdot \left[n_A\right]}{\left[W\right]}=\frac{\frac{\mathrm{g}}{\mathrm{L}}\cdot \frac{\mathrm{L}}{\mathrm{Cell}}\cdot \frac{1}{\text{mol}}}{\frac{\text{g}}{\text{mol}}}=\frac{1}{\mathrm{Cell}},$$
which is exactly as intended. By the same logic, the model variables $A_T$ and $T$ are represented in terms of number per cell.

## Dimensional analysis of the equations

We perform dimensional analysis on the equations of the COMBAT model is assessed whether the units used are compatible

### $B_x$
Since the measurement unit of $B_x$ is $\left[B_x\right]=\mathrm{Cells}$, the measurement unit of $\frac{\text {d}B_x}{\mathrm{d}t}$ is $\left[\frac{\text {d}B_x}{\mathrm{d}t}\right]=\frac{\left[B_x\right]}{\left[t\right]}=\frac{\mathrm{Cell}}{s}$. We now prepare to check whether the right side of this equation conforms with these dimensions. First, we assert that the quantity $\frac{C-\sum\limits _{j=0}^{n}B_j}{C}$ is dimensionless:
$$\left[\frac{C-\sum\limits_{j=0}^{n}B_j}{C}\right]=\frac{\left[C\right]-\left[\sum\limits_{j=0}^{n}B_j\right]}{\left[C\right]}=\frac{\mathrm{Cell}-\sum\limits_{j=0}^{n}\left[B_j\right]}{\mathrm{Cell}}=\frac{\mathrm{Cell}-\sum\limits_{j=0}^{n}\mathrm{Cell}}{\mathrm{Cell}}=\frac{\mathrm{Cell}}{\mathrm{Cell}}=1$$

Next, we find the measurement unit of $\rho_x$:

$$\begin{aligned}
\left[\rho_x\right]
=\left[2\right]\left[\sum\limits_{i=x}^{n}f_{i,x}r_i B_i \frac{C-\sum\limits _{j=0}^{n}B_j}{C}\right]=
1\cdot\sum\limits_{i=x}^{n}\left[f_{i,x}\right] \left[r_i\right] \left[B_i\right] \left[\frac{C-\sum\limits _{j=0}^{n}B_j}{C}\right]=\\
\sum\limits_{i=x}^{n}1\cdot \frac{1}{\mathrm{s}}\cdot\mathrm{Cell} \cdot 1=\sum\limits_{i=x}^{n}\frac{\mathrm{Cell}}{\mathrm{s}}=\frac{\mathrm{Cell}}{\mathrm{s}}
\end{aligned}$$

From this, we can proceed with the full expression for $\frac{\text {d}B_x}{\mathrm{d}t}$:

$$\begin{aligned}
\left[\frac{k_{f}}{Vn_{A}}(n-x+1)AB_{x-1} - k_{r}xB_x - \frac{k_{f}}{Vn_{A}}(n-x)AB_x + \right. \\ \left. \quad k_{r}(x+1)B_{x+1} + \rho _x -r_xB_x \frac{C-\sum _{j=0}^{n}B_j}{C} -d_xB_x\right] = \\
\left[\frac{k_{f}}{Vn_{A}}(n-x+1)AB_{x-1}\right] - \left[k_{r}xB_x\right] - \left[\frac{k_{f}}{Vn_{A}}(n-x)AB_x\right] + \\ \quad \left[k_{r}(x+1)B_{x+1}\right] + \left[\rho_x\right] - \left[r_xB_x \frac{C-\sum _{j=0}^{n}B_j}{C}\right] -\left[d_xB_x\right] = \\
\frac{\left[k_{f}\right]}{\left[V\right]\left[n_{A}\right]}\left[n-x+1\right]\left[A\right]\left[B_{x-1}\right] - \left[k_{r}\right]\left[x\right]\left[B_x\right] - \frac{\left[k_{f}\right]}{\left[V\right]\left[n_{A}\right]}\left[n-x\right]\left[A\right]\left[B_{x}\right] + \\ \quad \left[k_{r}\right]\left[x+1\right]\left[B_{x+1}\right] + \left[\rho_x\right] - \left[r_x\right]\left[B_x\right] \left[\frac{C-\sum _{j=0}^{n}B_j}{C}\right] -\left[d_x\right]\left[B_x\right] = \\
\frac{\frac{\mathrm{L}}{\mathrm{mol}\cdot \mathrm{s}}}{\frac{\mathrm{L}}{\mathrm{Cell}}\cdot{\frac{\mathrm{1}}{\mathrm{mol}}}}\cdot 1\cdot\frac{1}{\mathrm{Cell}}\cdot{\mathrm{Cell}} - \frac{1}{\mathrm{s}}\cdot 1\cdot\mathrm{Cell} - \frac{\frac{\mathrm{L}}{\mathrm{mol}\cdot \mathrm{s}}}{\frac{\mathrm{L}}{\mathrm{Cell}}\cdot{\frac{\mathrm{1}}{\mathrm{mol}}}}\cdot 1\cdot\frac{1}{\mathrm{Cell}}\cdot{\mathrm{Cell}} + \\
\quad \frac{1}{\mathrm{s}}\cdot 1\cdot\mathrm{Cell} + \frac{\mathrm{Cell}}{\mathrm{s}} - \frac{1}{\mathrm{s}}\cdot \mathrm{Cell} \cdot 1 - \frac{1}{\mathrm{s}}\cdot \mathrm{Cell} = \\
\frac{\mathrm{Cell}}{\mathrm{s}}\cdot 1 - \frac{\mathrm{Cell}}{\mathrm{s}} - \frac{\mathrm{Cell}}{\mathrm{s}}\cdot 1 + 
\\ \quad \frac{\mathrm{Cell}}{\mathrm{s}} + \frac{\mathrm{Cell}}{\mathrm{s}} - \frac{\mathrm{Cell}}{\mathrm{s}} - \frac{\mathrm{Cell}}{\mathrm{s}}=\\
\frac{\mathrm{Cell}}{\mathrm{s}} - \frac{\mathrm{Cell}}{\mathrm{s}} - \frac{\mathrm{Cell}}{\mathrm{s}} + 
\\ \quad \frac{\mathrm{Cell}}{\mathrm{s}} + \frac{\mathrm{Cell}}{\mathrm{s}} - \frac{\mathrm{Cell}}{\mathrm{s}} - \frac{\mathrm{Cell}}{\mathrm{s}}=\\
\frac{\mathrm{Cell}}{\mathrm{s}}
\end{aligned}$$
Hence, the units of measurement conform for $B_x$.
"""

# ╔═╡ 46f15089-a99c-41ed-8376-644e0c9331ac
md"""

### $A$
We now proceed with the equation for $A$. Since $A$ is internally represented in molecular units, we have $\left[A\right]=\frac{1}{\mathrm{Cell}}$ and hence $\left[\frac{\text {d}A}{\mathrm{d}t}\right]=\frac{\left[A\right]}{\left[t\right]}=\frac{\frac{1}{\mathrm{Cell}}}{s}=\frac{1}{\mathrm{Cell}\cdot\mathrm{s}}$. Before proceeding, we note the interpretation of the terms $\left(n-x\right)B_x$ and $xB_x$:

In the equation for $B_x$, the terms $\frac{k_{f}}{Vn_{A}}\left(n-x+1\right)AB_{x-1}$, $\frac{k_{f}}{Vn_{A}}\left(n-x\right)AB_{x}$, $k_{r}xB_x$, and $k_{r}\left(x+1\right)B_{x+1}$ denote the change in the number of bacteria in compartment $x$ due to binding or unbinding of the antibiotic. The factors $\left(n-x+1\right)$, $\left(n-x\right)$, $x$, and $\left(x+1\right)$ relate to the number of bound molecules or the number of unbound targets in their respective compartments. The way these factors are used in the equation is to denote that the rate of change is proportional to these constants and hence their units of measurements is dimensionsless. In the corresponding terms for $\left(n-x\right)B_x$ and $xB_x$, the quantities relates to how many antibiotic molecules are taken up or relased from compartment $B_x$. There are $B_x$ cells which each have $x$ antibiotic molecules to release and the potential to take up $n-x$ molecules. For this to make sense, we must interpret the factors $\left(n-x\right)$ and $x$ to have the unit of $\frac{1}{\mathrm{Cell}}$.

For the right side of the equation for $\frac{\text{d}A}{\mathrm{d}t}$, we get:


$$\begin{aligned}
\left[- \frac{k_{f}}{Vn_{A}}(A\cdot T +\sum _{x=0}^{n-1}(n-x)AB_x) + k_{r}\left(A_T+\sum _{x=1}^{n}xB_x\right)\right]=\\
- \left[\frac{k_{f}}{Vn_{A}}(A\cdot T +\sum _{x=0}^{n-1}(n-x)AB_x)\right] + \left[k_{r}\left(A_T+\sum _{x=1}^{n}xB_x\right)\right]=\\
- \left[\frac{k_{f}}{Vn_{A}}\right]\left[A\cdot T +\sum _{x=0}^{n-1}(n-x)AB_x\right] + \left[k_r\right]\left[A_T+\sum _{x=1}^{n}xB_x\right]=\\
-\frac{\left[k_{f}\right]}{\left[V\right]\left[n_{A}\right]}\cdot\left(\left[A\cdot T\right]+\left[\sum _{x=0}^{n-1}(n-x)AB_x\right]\right)+\left[k_r\right]\left(\left[A_T\right]+\left[\sum _{x=1}^{n}xB_x\right]\right)=\\
-\frac{\frac{\mathrm{L}}{\mathrm{mol}\cdot \mathrm{s}}}{\frac{\mathrm{L}}{\mathrm{Cell}}\cdot\frac{1}{\mathrm{mol}}}\cdot\left(\left[A\right]\cdot \left[T\right]+\sum _{x=0}^{n-1}\left[(n-x)AB_x\right]\right)+\frac{1}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}+\sum_{x=1}^{n}\left[xB_x\right]\right)=\\
-\frac{\mathrm{Cell}}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}\cdot \frac{1}{\mathrm{Cell}}+\sum _{x=0}^{n-1}\left[n-x\right]\left[A\right]\left[B_x\right]\right)+\frac{1}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}+\sum_{x=1}^{n}\left[x\right]\left[B_x\right]\right)=\\
-\frac{\mathrm{Cell}}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}^2}+\sum _{x=0}^{n-1}\frac{1}{\mathrm{Cell}}\cdot \frac{1}{\mathrm{Cell}}\cdot \mathrm{Cell}\right)+\frac{1}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}+\sum_{x=1}^{n}\frac{1}{\mathrm{Cell}}\cdot\mathrm{Cell}\right)=\\
-\frac{\mathrm{Cell}}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}^2}+\sum _{x=0}^{n-1}\frac{1}{\mathrm{Cell}}\right)+\frac{1}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}+\sum_{x=1}^{n}1\right)=\\
-\frac{\mathrm{Cell}}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}^2}+\frac{1}{\mathrm{Cell}}\right)+\frac{1}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}+1\right)=\\
-\frac{1}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}+1\right)+\frac{1}{\mathrm{s}}\cdot\left(\frac{1}{\mathrm{Cell}}+1\right)
\end{aligned}$$

Here, we arrive at a paradox, the equation tells us to add quantities of unit $1$ (dimensionless) and $\frac{1}{\mathrm{Cell}}$ which is impossible. Hence, the equation must be inconsistent. Also note that the terms originating from $A\cdot T$ and $A_T$ have the expected unit of $\frac{1}{\mathrm{Cell}\cdot\mathrm{s}}$, whereas the terms originating from $\left(n-x\right)AB_x$ and $xB_x$ have the *wrong* unit of $\frac{1}{\mathrm{s}}$. Intuitively, we can reason that the terms $\frac{k_{f}}{Vn_{A}}\cdot (n-x)AB_x$ and $k_rxB_x$ corresponds to the total number of antibiotic molecules being absorbed or released in total for the entire system and not for each individual cell. We could potentially resolve this problem by modifying the equation to:

$$\frac{\text {d}A}{\mathrm{d}t} = - \frac{k_{f}}{Vn_{A}}\left(A\cdot T +\sum _{x=0}^{n-1}\frac{(n-x)AB_x}{\sum\limits_{x=0}^{n}B_x}\right) + k_{r}\left(A_T+\sum _{x=1}^{n}\frac{xB_x}{\sum\limits_{x=0}^{n}B_x}\right),$$

but there is another major disadvantage of this choice of units of measurement: The total number of cells in the system changes when bacteria divide and die even when the total number of antibiotic molecules stay the same.
"""

# ╔═╡ 2b984af8-e36c-4d23-ac30-99cf230855f1
md"""
## Solution to dimension problem
It can be shown that interpreting $V$, as the total volume of which the bacteria are part resolves the dimensionality error. This is $V$ should have the unit of measurement $\left[V\right]=\mathrm{L}$.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "8aa109ae420d50afa1101b40d1430cf3ec96e03e"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═2106bf0c-1580-4d22-81fe-6a74f4102d97
# ╠═53a4d7b0-9aa8-11ef-2a21-01318b24dbbe
# ╠═9b9b544a-e7b7-4907-a271-d862f99a2de0
# ╟─6fb53591-d97f-429c-b879-3b7030e65719
# ╟─46f15089-a99c-41ed-8376-644e0c9331ac
# ╠═2b984af8-e36c-4d23-ac30-99cf230855f1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
