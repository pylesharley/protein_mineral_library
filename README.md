# protein_mineral_library

This GitHub repository contains the design models, and the scripts that generated them, associated with the manuscript:

Controlling semiconductor growth with structured de novo protein interfaces
Saragovi A., Pyles H., et al. 2024
DOI: 


Adjust this path to your download location in the following examples

	PATH = /path/to/this/repository

Models:

	Models of all proteins selected for expression and characterization are shown in main text figures and included in PATH/selected_designs/

	Models of all 6730 designs in the yeast displayed library are found in PATH/library_models/

Scripts:

The scripts used to generate the library of potential mineral-binding interfaces are in PATH/library_design/

	Step 1
		01_check_in_dhrs.ipynb 
		
		Takes directory of DHR models (any number of repeats) for example: PATH/library_design/scaffolds/input_dhrs
			need to have interger number of complete repeats! I.E. no trimmed ends. 

		Parses out repeat length from sequence features
			assumes repeat is < 60 residues long by default

		Writes 3 repeat, capped models of DHRs into: PATH/library_design/scaffolds
			also prints tsv with repeat lenght and surface selections for each DHR to PATH/library_design/scaffolds/dhr_surface_raw_selections.tsv

			these selctions need to be double checked, modified as needed, and copied to PATH/library_design/scaffolds/dhr_surface_cooked_selections.tsv


	Step 2
		02_setup_for_sequence_sampling.ipynb
		
		This takes defined amino acid proportions from PATH/library_design/constraints/AA_compositions.tsv and combines it with the input from step 1

		It writes many input files for Rosetta script runs

		AA_compositions.tsv are parsed and written as AA_comp files in PATH/library_design/constraints/

		AA comp constraints are applied to residue selectors defined by PATH/library_design/scaffolds/surface_cooked_selections.tsv

		For each surface + AA comp combination, a repeat and anti-repeat version is produced
			repeat version == strict repeat sym is enforced in surface selection
			antirepeat version == differences at repeat positions are incendized with additional AA comp file (a penalty is applied if same residues is at same position)

		Run the newly generated Rosetta scripts


	Step 3
		03_plot_output.py

		Extracts low energy PDB for each combintorially defined design (Scaffold + AA_comp + Surface + Repeat/Antirepeat) 
			AND the pdb from the lowest energy quartile that has the most sequence differences to the low energy one. 

		Determines score terms and biochemical parameters of sequences from output of step 2

	Step 4
		04_remove_redundant_sequences.py

		Remove sequences that are too similar to other sequences to avoid redundancy in ordered genes

	Step 5
		05_count_surface_compositions.py

		Quantify sequence compositions of interfaces produced by Rosetta in step 3

	Step 6
		06_plot_metrics.py

		Plot scores and sequence compositions of interfaces produced by Rosetta in step 3


The scripts used to generate the designed beta solnoid repeat scaffolds are in PATH/beta_solenoid_design

	Step 1
		01_generate_sequences.ipynb

		This script generates repeat sequences with random length ranges of beta-sheet propensity regions (QxQx) repeats where x positions are hydrophobic residues, connected by random stretches of loop propensity residues (GSTPND)

	Step 2
		Predict the structures of these sequences with RoseTTAFold (M. Baek, et al. 2021)

		Select best scoring structures with 02_select_sequences.ipynb

	Step 3
		03_mutate_sequences.ipynb

		Mutate residues in best scoring sequences to diversify pool

	Step 4
		Repeat steps 2 & 3 until sufficent numbers of well-folding solenoid backbones are obtained
	
	Step 5
		Predict the structures of these sequences with AlphaFold (J. Jumper, et al. 2021)

		Compare structures predicted by RoseTTAFold and AlphaFold and select designs with low RMSD differences with 04_compare_RF_to_AF.ipynb

		Visually inspect models and select ones with flat surfaces and well-packed hydrophobic cores

	Step 6
		06_write_design_xmls.ipynb generates RosettaScript XMLs that redesign the beta-solenoid surfaces to reduce repeativness and promote solubility

		See PATH/beta_solenoid_design/scaffold_surface_design/ for the design scripts used in this study

Scripts used to assign NGS reads to designs in the library and select enriched designs are in PATH/identify_enriched_sequences/

Scripts used to rescaffold Z0-fiber motifs into diffused alpha-beta backbones are in PATH/motif_grafting/

Scripts used to generate oligomers with the Z4 interface are in PATH/oligomer_design
