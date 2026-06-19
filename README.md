# Yaskawa HC10 DTP – Jacobian & Singularity Analysis (MATLAB)

A MATLAB implementation of geometric Jacobian computation and kinematic singularity analysis for the **Yaskawa HC10 DTP** collaborative robot, using the **Denavit-Hartenberg (DH)** convention.

---

## Contents

| File | Description |
|------|-------------|
| `dh_transform.m` | Computes a single 4×4 DH transformation matrix |
| `forward_kinematics.m` | Forward kinematics: chains all 6 DH matrices to get TCP pose |
| `singularity_analysis.m` | Main script: reads joint data, computes Jacobian & determinant |

---

## Background & Theory

### The Robot

<div align="center">
<table>
<tr>
<td valign="top">

<img src="HC10_Farbe_gesch.png" height="340" alt="Yaskawa HC10 DTP with colour-coded joint rings">

</td>
<td valign="top">

| i | dᵢ [mm] | ϑᵢ | aᵢ [mm] | αᵢ |
|---|---|---|---|---|
| 1 | 275 | ![v1](https://img.shields.io/badge/%CF%911-1f3a93) | 0 | −π/2 |
| 2 | 0 | ![v2](https://img.shields.io/badge/%CF%912-e67e22) − 90° | 700 | π |
| 3 | 0 | ![v3](https://img.shields.io/badge/%CF%913-27ae60) | 0 | −π/2 |
| 4 | −500 | ![v4](https://img.shields.io/badge/%CF%914-3498db) | 0 | π/2 |
| 5 | −162 | ![v5](https://img.shields.io/badge/%CF%915-8e44ad) | 0 | −π/2 |
| 6 | −170 | ![v6](https://img.shields.io/badge/%CF%916-e74c3c) | 0 | π |

</td>
</tr>
</table>
</div>

The [Yaskawa HC10 DTP](https://www.yaskawa.eu.com) is a 6-DOF collaborative robot. Each joint is marked with a colored ring — the same colors are used in the DH parameter table above.

---

### Denavit-Hartenberg Convention

Denavit and Hartenberg introduced a convention that uses only **four parameters** to describe the relationship between two consecutive coordinate frames (instead of the general six needed for an arbitrary rigid body transformation). For a **revolute joint**, all parameters except θᵢ are constant.

The DH transformation matrix **Aᵢ** is defined as the product of two translations and two rotations:

```
A = Tz(d) · Rz(θ) · Tx(a) · Rx(α)
```

This gives the following 4×4 homogeneous matrix:

```
        | cos(θ)  -sin(θ)·cos(α)   sin(θ)·sin(α)  a·cos(θ) |
A_i =   | sin(θ)   cos(θ)·cos(α)  -cos(θ)·sin(α)  a·sin(θ) |
        |   0         sin(α)           cos(α)           d    |
        |   0           0               0               1    |
```

where:
- **d** – link offset along the previous z-axis [mm]
- **θ** – joint angle about the previous z-axis [rad]
- **a** – link length (perpendicular distance between z-axes) [mm]
- **α** – link twist angle about the x-axis [rad]

---

### DH Parameters – Yaskawa HC10 DTP

The parameters below were determined to match the Yaskawa MotoSim simulation environment using the **absolute coordinate system** (origin at the intersection of axes 1 and 2). Each row is color-coded to match the joint ring colors on the robot image above.

> **Note on Joint 2:** A −90° offset on θ₂ is required to align the DH model with the robot's physical zero position. This was verified against MotoSim for both the zero configuration and the candle position (θ₃ = 90°).

---

### Forward Kinematics

The overall transformation from base frame to TCP is obtained by chaining the six individual DH matrices:

```
T_TCP = A₁ · A₂ · A₃ · A₄ · A₅ · A₆
```

The **translational part** (column 4, rows 1–3) of T_TCP gives the TCP position in millimeters. The **rotational part** (3×3 upper-left block) gives the TCP orientation.

---

### Geometric Jacobian

The geometric Jacobian **J** relates joint velocities **q̇** to Cartesian velocities of the TCP:

```
| ṗₓ |           | q̇₁ |
| ṗᵧ |           | q̇₂ |
| ṗᵤ | = J(q) ·  | q̇₃ |
| φ̇  |           | q̇₄ |
| θ̇  |           | q̇₅ |
| ψ̇  |           | q̇₆ |
```

For the HC10 with 6 revolute joints, **J** is a 6×6 matrix. Each column is computed geometrically using the z-axis and origin of the preceding frame:

```
       | z_{i-1} × (p_TCP − p_{i-1}) |
J_i =  |                             |
       |           z_{i-1}           |
```

where:
- **z_{i-1}** – unit vector of the z-axis of frame i−1 (column 3 of the cumulative transformation matrix)
- **p_{i-1}** – origin of frame i−1 (column 4 of the cumulative transformation matrix)
- **p_TCP** – TCP position (column 4 of T_TCP)

The full Jacobian is assembled as:

```
J_G = [ J₁  J₂  J₃  J₄  J₅  J₆ ]
```

---

### Singularity Detection

A kinematic **singularity** occurs when the robot loses one or more degrees of freedom in a given configuration. This happens when the Jacobian becomes rank-deficient, i.e. when:

```
det(J) = 0
```

The main script computes `det(J)` for every recorded configuration. Values far from zero indicate high manipulability (the robot can move freely in all directions). Values approaching zero indicate proximity to a singularity.

---

## Usage

### 1. Prepare your input file

Create a **CSV file** with N rows and 6 columns — one row per timestep, one column per joint angle in **degrees**. No header row.

```
-43.473, 53.912,  9.587, -0.011, -45.662, -46.557
-43.211, 54.001,  9.612, -0.009, -45.701, -46.489
...
```

> Excel files (`.xlsx`) are also supported. If you export from a Yaskawa controller or MotoPlus application, the data is typically already in the correct 6-column format.

### 2. Set the input filename

Open `singularity_analysis.m` and set the `INPUT_FILE` variable at the top:

```matlab
INPUT_FILE = 'joint_angles.csv';   % or 'MyData.xlsx'
```

### 3. Run the script

```matlab
>> singularity_analysis
```

The script will print progress and save four output CSV files:

| Output file | Contents |
|-------------|----------|
| `determinant.csv` | Jacobian determinant per timestep |
| `position_x.csv` | TCP x-position [mm] |
| `position_y.csv` | TCP y-position [mm] |
| `position_z.csv` | TCP z-position [mm] |

A figure with TCP position plots and the determinant trace is shown automatically.

### 4. Use individual functions directly

```matlab
% Single forward kinematics calculation
angles_deg = [0, 0, 0, 0, 0, 0];
offsets    = [0, -90, 0, 0, 0, 0];   % DH convention offsets
A = deg2rad(angles_deg + offsets);
T = forward_kinematics(A(1), A(2), A(3), A(4), A(5), A(6));

fprintf('TCP position: x=%.1f  y=%.1f  z=%.1f mm\n', T(1,4), T(2,4), T(3,4));
```

```matlab
% Single DH matrix
T = dh_transform(275, deg2rad(0), 0, -pi/2);
```

---

## Requirements

- **MATLAB R2019b** or later (uses `readmatrix` and `writematrix`)
- No additional toolboxes required
- For `.xlsx` input: the `xlsread` fallback is included for older MATLAB versions

---

## Notes on the DH Parameter Derivation

Finding the correct DH parameters for this robot required careful comparison against the Yaskawa MotoSim environment. Key observations:

- **d₁ = 275 mm**: Yaskawa places the first frame at the intersection of axes 1 and 2, not at the base flange. Setting d₁ = 0 would shift all z-values by 275 mm.
- **θ₂ offset (−90°)**: Without this, the zero configuration of the model does not match the physical robot's zero position.
- **Signs of d**: The negative signs on d₄, d₅, d₆ reflect that those frames translate in the negative z-direction of the preceding frame, as verified by the MotoSim coordinate display.
- **α values**: The alternating π and −π/2 values follow from the physical arrangement of the joint axes and were determined by matching both the zero position and the candle position (θ₃ = 90°) in MotoSim.

---

## License

MIT License – feel free to use, modify, and distribute with attribution.

---

## References

- Denavit, J. & Hartenberg, R.S. (1955). *A kinematic notation for lower-pair mechanisms based on matrices.* Journal of Applied Mechanics.
- Siciliano, B. et al. (2009). *Robotics: Modelling, Planning and Control.* Springer.
- Yaskawa Electric Corporation – HC10 DTP product documentation.
