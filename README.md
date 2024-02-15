
# Simulation of Combined Antibiotics in Cross-feeding Communities(SCARCC)

---
## Setting Up Your Environment

To ensure that the Python interpreter can find and import the `scarcc` package, you need to add the source directory to your `PYTHONPATH`. Follow the instructions below specific to your operating system.

### For Linux and macOS

Open a terminal window and navigate to the root of the `scarcc` project. Then, you can temporarily set your `PYTHONPATH` by running:

```bash
export PYTHONPATH="${PYTHONPATH}:/path/to/scarcc/src"
```

To make this change permanent, you can add the export command to your shell's startup script, such as `.bashrc` or `.bash_profile` for Bash, or `.zshrc` for Zsh, replacing `/path/to/scarcc/src` with the actual path to the `scarcc` source directory.

For example, if you're using Bash, you could add the following line to the end of your `.bashrc`:

```bash
echo 'export PYTHONPATH="${PYTHONPATH}:/path/to/scarcc/src"' >> ~/.bashrc
```

Don't forget to replace `/path/to/scarcc/src` with the actual path to the directory. Then, source your `.bashrc` file to apply the changes:

```bash
source ~/.bashrc
```

### For Windows

On Windows, you can set the `PYTHONPATH` environment variable through the System Properties.

1. Search for "Environment Variables" in your Windows search and select "Edit the system environment variables".
2. In the System Properties window, click on the "Environment Variables…" button.
3. In the Environment Variables window, under the "System variables" section, click on the "New…" button to create a new environment variable.
4. Enter `PYTHONPATH` as the variable name and the full path to the `scarcc/src` directory as the variable value.
5. Click OK to close each window.

Alternatively, you can set the `PYTHONPATH` in a Command Prompt or PowerShell session temporarily by running:

```cmd
set PYTHONPATH=%PYTHONPATH%;C:\path\to\scarcc\src
```

Make sure you replace `C:\path\to\scarcc\src` with the actual path to the `scarcc` source directory on your system.

### Verifying the Setup

To verify that your `PYTHONPATH` is set correctly, you can open a Python interpreter and attempt to import the `scarcc` package:

```python
import scarcc
print("Scarcc is successfully imported!")
```

If you do not encounter any import errors, then you have successfully added the `scarcc` source directory to your `PYTHONPATH`.

---


scarccpy
├─ .vscode
│  ├─ PythonImportHelper-v2-Completion.json
│  └─ settings.json
├─ models
│  ├─ iJO1366.xml
│  ├─ iML1515.xml
│  ├─ iML1515_E0.xml
│  ├─ jmc_AM1_KO.txt
│  ├─ jmc_AM1_KO.xml
│  ├─ jmc_AM1_KO_renamed.xml
│  ├─ jmc_EC_iJO_KO.txt
│  ├─ jmc_EC_iJO_KO.xml
│  ├─ jmc_SalMet.txt
│  ├─ jmc_SalMet.xml
│  ├─ STM_v1_0.xml
│  └─ STM_v1_0_S0.xml
├─ src
│  ├─ .vscode
│  │  └─ PythonImportHelper-v2-Completion.json
│  ├─ notebook
│  │  ├─ monoculture_alpha_search.ipynb
│  │  └─ notebook_import.ipynb
│  └─ scarcc
│     ├─ .vs
│     │  └─ slnx.sqlite
│     ├─ .vscode
│     │  └─ PythonImportHelper-v2-Completion.json
│     ├─ data_analysis
│     │  ├─ biomass
│     │  │  ├─ growth_summary.py
│     │  │  ├─ pathway_summary.py
│     │  │  ├─ plot_growth_curve.py
│     │  │  └─ __init__.py
│     │  └─ flux
│     │     ├─ flux_correction.py
│     │     ├─ flux_ratio_analysis.py
│     │     ├─ flux_snapshot.py
│     │     ├─ __init__.py
│     │     └─ __pycache__
│     │        ├─ flux_correction.cpython-311.pyc
│     │        └─ __init__.cpython-311.pyc
│     ├─ generate_alpha_table.py
│     ├─ notebook_import.ipynb
│     ├─ preparation
│     │  ├─ metabolic_model
│     │  │  ├─ basic_model.py
│     │  │  └─ core
│     │  │     ├─ component.py
│     │  │     ├─ medium.py
│     │  │     └─ __init__.py
│     │  ├─ perturbation
│     │  │  ├─ function
│     │  │  │  ├─ iter_species.py
│     │  │  │  ├─ stoichiometry_scaling.py
│     │  │  │  └─ __init__.py
│     │  │  ├─ parameter_derivation
│     │  │  │  ├─ alpha_finder.py
│     │  │  │  ├─ process_double_gene.py
│     │  │  │  ├─ __init__.py
│     │  │  │  └─ __pycache__
│     │  │  │     ├─ basic_model.cpython-311.pyc
│     │  │  │     ├─ component.cpython-311.pyc
│     │  │  │     ├─ component.cpython-312.pyc
│     │  │  │     ├─ medium.cpython-311.pyc
│     │  │  │     ├─ __init__.cpython-311.pyc
│     │  │  │     └─ __init__.cpython-312.pyc
│     │  │  ├─ __init__.py
│     │  │  └─ __pycache__
│     │  │     ├─ iter_species.cpython-311.pyc
│     │  │     ├─ iter_species.cpython-312.pyc
│     │  │     ├─ __init__.cpython-311.pyc
│     │  │     └─ __init__.cpython-312.pyc
│     │  ├─ target_gene
│     │  │  └─ __init__.py
│     │  └─ __init__.py
│     ├─ sim_COMETS.py
│     ├─ sim_COMETS.slurm
│     ├─ sim_engine
│     │  ├─ data_processing
│     │  │  ├─ preprocess_flux.py
│     │  │  └─ __init__.py
│     │  ├─ result_processing.py
│     │  ├─ simulation_builder.py
│     │  ├─ simulation_parameter.py
│     │  ├─ simulation_workflow.py
│     │  ├─ __init__.py
│     │  └─ __pycache__
│     │     ├─ result_processing.cpython-311.pyc
│     │     └─ __init__.cpython-311.pyc
│     ├─ test.py
│     ├─ util
│     │  ├─ README.md
│     │  ├─ util.py
│     │  ├─ __init__.py
│     │  └─ __pycache__
│     │     ├─ util.cpython-311.pyc
│     │     ├─ util.cpython-312.pyc
│     │     ├─ __init__.cpython-311.pyc
│     │     └─ __init__.cpython-312.pyc
│     ├─ __init__.py
│     └─ __pycache__
│        ├─ test.cpython-311.pyc
│        ├─ __init__.cpython-311.pyc
│        └─ __init__.cpython-312.pyc
└─ test
   ├─ e1.py
   ├─ test_basic_model.py
   └─ __init__.py

```