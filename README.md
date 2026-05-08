# Cell vertex model with cell division for vein growth



---

## Description

Source codes for numerical simulations in Sugimura et al. [1].
- `CVN_Vein.cpp`: Main code
- `CVM_Parameters.cpp`: Parameters used in this code. This file is included by CVM_Vein.cpp.
- `CVM_CellDivScenario.cpp`: Control schedule of cell division time. This file is included by CVM_Vein.cpp.
- `CVM_DataOutput.cpp`: Data output
- `lib`: The other libraries for CVM, included by CVM_Vein.cpp.
- `Initialconditions`: Initial cell configuration data
---


## Requirements

* g++ with C++17 support  
* Gnuplot (for visualization, optional)

---

## Usage

1. Clone the repository and inspect the directory structure:  
   ```text
   project/
   ├── CVN_Vein.cpp
   ├── CVM_Parameters.cpp
   ├── CVM_CellDivScenario.cpp
   ├── CVM_DataOutput.cpp
   ├── lib/
   │   ├── geometry.cpp
   │   ├── gnuplot.cpp
   │   ├── head_cell_div.cpp
   │   ├── head_cell_vertex.cpp
   │   ├── IO.cpp
   │   ├── libVein.cpp
   │   ├── ODE.cpp
   │   └── TimeEvolution.cpp
   └── README.md

2. Set the parameters in `CVM_Parameters.cpp` (and `CellDivScenario.cpp` if you want to change the cell division schedule).  `CVM_DataOutput.cpp` describes output. 

3. Compile codes 
```sh
g++ -std=c++17 CVM_Vein.cpp -O3 -o cvm.out
```  

4. Run the simulation:
```sh
mkdir -p results
./cvm.out Initialconditions/<file>
<file>: file in Initialconditions
```
Data and gif animation are generated in results/ .

5. To produce an animated GIF: Edit GNUPLOT_USE to 2 in `CVM_Parameters.cpp` 
   Recompile as above and run the simulation.  

---
## Parameters
 See K. Sugimura et al. [1].


---

## References
 [1] Kaoru Sugimura, Ryu Takayanagi, Toshinori Namba, Zeping Qu, and Shuji Ishihara. Cell size control emerges from the vein-dependent coordinated divisions of distinct cell groups in Drosophila wing. submitted (bioRxiv https://doi.org/10.64898/2025.12.26.696565)



```python

```
