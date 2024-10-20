# Project: Derivatives Hedging

## 1. **Objective:**

A project for pricing and Delta hedging of derivatives: Options (Basket, Performance, Asian).

### - **Pricing:**

We model the diffusion dynamics of the underlying asset using a Black-Scholes model:

$$
S_{t,d} = S_{0,d} e^{\left(r - \frac{\sigma_d^2}{2}\right)t + \sigma_d B_{t,d}}, \quad d = 1, \dots, D
$$

The option price is then estimated using **Monte Carlo**.

### - **Hedging:**

We use Delta hedging for our strategies. This means calculating the sensitivities (deltas) of the option price with respect to the underlying asset to construct a hedging portfolio with delta quantities of each underlying asset.

The deltas are calculated as:

$$
\frac{\partial v(t, S_{t_0}, \dots, S_{t_i}, S_t)}{\partial S_t}
$$


## 2. **Development Environment:**

-   **Installation of `PNL` library:**

```bash
mkdir lib
cd lib
cp /relative/path/to/pnl ./
mkdir build
cd build
cmake ..
chmod +x ../split_linker_command.sh
make
make install
```

-  **Installation of `nlohmann-json` library:**

```bash
sudo apt update
sudo apt install nlohmann-json3-dev
```

-   **compilation:**

```bash
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=relative/path/to/lib/pnl/build/ ..
make
```

- **Execution:**

```bash
./price0 data_input.json
./hedge market_file.txt data_input.json
```

-   **vscode configuration:**

    -   Ouvrir le fichier `c_cpp_properties.json`
    -   Ajouter dans `includePath`: **`"relative/path/to/lib/pnl/build/include/"`**

-   **Add executable in `CMakeLists.txt`:**

```Makefile
add_executable(nomExecutable _liste_des_fichiers_cpp)
target_link_libraries(nomExecutable
    ${PNL_LIBRARIES}
    nlohmann_json::nlohmann_json
)

```

-   **Tests:**

```bash
python3  ./CheckArchive.py --check --extract --destdir=tests_rendu/rendu/  --build --pnldir=lib/pnl/build/  tests_rendu/Equipe_2.tar.gz
python3 ./testForPCPD.py --price --exec=tests_rendu/rendu/Equipe_2/build/price0  --datadir=tests_rendu/data/  --outdir=tests_rendu/out
```
