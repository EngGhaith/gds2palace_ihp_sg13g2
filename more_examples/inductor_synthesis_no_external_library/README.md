# Inductor synthesis IHP inductor2/inductor3 using gds2palace

The gds2palace example provided here creates a symmetric inductor2/inductor3 layout in IHP SG13G2 technology fully automatically, based on a target L value at a target frequency. The geometry sweep range where possible implementations are evaluated is defined in the model code: range for width, spacing, number of turns.  

New in this version: The layout is created by code in the model file, there is no need for an external library.

## Principle of operation

This is what this script does for you:

- Step 1: Determine list of possible implementations based on close form equations, by calculating the required diameter and testing if this layout is valid and within a given maximum diameter limit. 
- Step 2: Create GDSII layouts for all these candidates (including ports + ground return) 
- Step 4: Run a fast FEM sweep over all candidates using gds2palace
- Step 5: Evaluate the n candidates with highest Q factor at target frequency and re-tune to target value at target frequency
- Step 6: After m iterations over step 5, select the candidate with the highest Q factor and do a wideband full sweep using gds2palace FEM with full accuracy.
- Step 7: Plot results for L and Q factor of that best candidate
- Step 8: Create a final GDSII file with all extra layout features required for SG13G2 OPDK.

The gds2palace FEM simulation flow runs in non-GUI mode here, so that there is no user action required while the script is processing data. 


## Usage

Acticate the Python venv where you can run gds2palace models. gds2palace must be installed as a Python module: pip install gds2palace. The Palace solver must be available and you must be able to run gds2palace models. If you are not familar with gds2palace, go to the gds2palace documentation [here](https://github.com/VolkerMuehlhaus/gds2palace_ihp_sg13g2).

In the `synthesize_ihp_inductor_v1.py` script, set your target L value and target frequency, and adjust the search range for w,s and number of turns. Then just run the Python script.

### Inductor target and geometry range

'''python
'''


## Change history

07-April-2026: New version with built-in inductor geometry code, no external geometry library required. No limit on number of turns.

![inductor_design](../../doc/png/rfic_inductor_finder3.png)

