# Cosmic Variance Modeling and Survey Analysis

## Overview
This project implements a cosmic variance estimation pipeline based on the methodology from  
"A Cosmic Variance Cookbook" (Moster et al., 2011).

It computes the cosmic variance (σ_v) and dark matter variance (σ_DM) as a function of:
- Mean redshift (z)
- Stellar mass (log(M/M☉))
- Survey geometry

The pipeline supports multiple astronomical surveys and provides comparative analysis through visualizations and structured datasets.


## Features

- Computes cosmic variance for the following surveys:
  - UDF (3.3′×3.3′)
  - GOODS (10′×16′)
  - GEMS (28′×28′)
  - EGS (10′×70′)
  - COSMOS (10′×70′)

- Implements:
  - Dark matter variance calculation
  - Galaxy bias modeling
  - Survey-dependent scaling

- Generates:
  - Clean CSV datasets for each survey
  - Multi-dimensional parameter grids

- Provides visualizations:
  - Cosmic variance vs redshift (for fixed stellar mass)
  - Cross-survey comparisons
  - σ_DM evolution

---

## Methodology

### Dark Matter Variance
```

σ_DM = σ_a / (z^β + σ_b)

```

### Galaxy Bias Model
```

b(z) = b_0 (z + 1)^b_1 + b_2

```

### Cosmic Variance
```

σ_v = b(z) * σ_DM * sqrt(0.2 / Δz)

```

---

## Project Structure

```

.
├── cosmic_variance_module.py
├── cosmic variance plots.ipynb
├── UDF.csv
├── GOODS.csv
├── GEMS.csv
├── EGS.csv
├── COSMOS.csv
└── README.md

````

---

## Installation

```bash
git clone https://github.com/your-username/cosmic-variance-project.git
cd cosmic-variance-project
pip install -r requirements.txt
````

---

## Usage

### Import the module

```python
from cosmic_variance_module import get_cosmic_variance_array
```

### Run the computation

```python
mean_z = [0.1, 0.3, 0.5]
delta_z = 0.2
mass_array = [8.75, 9.25, 9.75]

cosmic_variance, sigma_DM = get_cosmic_variance_array(
    mean_z, delta_z, mass_array, survey="COSMOS"
)

print(cosmic_variance)
print(sigma_DM)
```

---

## Output

* Cosmic variance matrix:

```
shape = (len(mean_z), len(mass_array))
```

* Dark matter variance (σ_DM) values per redshift

---

## Visualizations

The notebook includes:

* Multi-panel comparison plots
* Survey-wise variance trends
* Mass-dependent variance evolution

---

## Key Insights

* Cosmic variance increases with stellar mass
* Smaller surveys (e.g., UDF) show higher variance
* Larger surveys (e.g., COSMOS) provide more stable estimates
* Dark matter variance decreases with increasing redshift

---

## Tech Stack

* Python
* NumPy
* Pandas
* Matplotlib

---

## Reference

Moster et al., 2011
"A Cosmic Variance Cookbook"

---

## Author

Kaveri Neeli
B.Tech Space Science and Engineering

---

## License

This project is open source and available under the MIT License.

```
```
