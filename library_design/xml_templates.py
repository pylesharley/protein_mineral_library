#! /usr/bin/env python

xml_string = '''
<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="BetaNov16" weights="beta_nov16"/>

        <ScoreFunction name="DesignBetaNov16Cst" weights="beta_nov16">
            <Reweight scoretype="atom_pair_constraint" weight="1.0"/>
            <Reweight scoretype="dihedral_constraint" weight="1.0"/>                
            <Reweight scoretype="coordinate_constraint" weight="1.0"/>
            <Reweight scoretype="aa_composition" weight="1.0"/>
            <Reweight scoretype="arg_cation_pi" weight="3" />
            <Reweight scoretype="approximate_buried_unsat_penalty" weight="5.0" />
            <Set approximate_buried_unsat_penalty_assume_const_backbone="true" />
            <Set approximate_buried_unsat_penalty_natural_corrections1="true" />
            <Set approximate_buried_unsat_penalty_hbond_energy_threshold="-0.5" />
            <Set approximate_buried_unsat_penalty_hbond_bonus_cross_chain="0.0" />
            <Set approximate_buried_unsat_penalty_hbond_bonus_ser_to_helix_bb="0.0" />
            <Set approximate_buried_unsat_penalty_burial_probe_radius="3.5"/>
            <Set approximate_buried_unsat_penalty_burial_atomic_depth="3.0"/>
        </ScoreFunction>
    </SCOREFXNS>

    <RESIDUE_SELECTORS>
{0}
        <Chain name="chain_a" chains="A" />
    </RESIDUE_SELECTORS>
    <SIMPLE_METRICS>
        <SasaMetric name="total_sasa" residue_selector="chain_a" sasa_metric_mode="all_sasa" />
        <SasaMetric name="total_polar_sasa" residue_selector="chain_a" sasa_metric_mode="polar_sasa" />
        <SasaMetric name="total_hydrophobic_sasa" residue_selector="chain_a" sasa_metric_mode="hydrophobic_sasa" />
        <SapScoreMetric name="sap_score" />
        <TotalEnergyMetric name="tot_Rscore" scoretype="total_score" scorefxn="BetaNov16" />
        <SelectedResidueCountMetric name="total_no_res_chA" residue_selector="chain_a" />
        <CalculatorMetric name="score_resnorm" equation="TOTAL_SCORE / TOTAL_NUMBER_RESIDUES" >
            <VAR name="TOTAL_SCORE" metric="tot_Rscore" />
            <Var name="TOTAL_NUMBER_RESIDUES" metric="total_no_res_chA"/>
        </CalculatorMetric>
        <CalculatorMetric name="SASA_resnorm" equation="TOTAL_SCORE / TOTAL_NUMBER_RESIDUES" >
            <VAR name="TOTAL_SCORE" metric="total_sasa" />
            <Var name="TOTAL_NUMBER_RESIDUES" metric="total_no_res_chA"/>
        </CalculatorMetric>
        <CalculatorMetric name="polarSASA_resnorm" equation="TOTAL_SCORE / TOTAL_NUMBER_RESIDUES" >
            <VAR name="TOTAL_SCORE" metric="total_polar_sasa" />
            <Var name="TOTAL_NUMBER_RESIDUES" metric="total_no_res_chA"/>
        </CalculatorMetric>
        <CalculatorMetric name="hydrophobicSASA_resnorm" equation="TOTAL_SCORE / TOTAL_NUMBER_RESIDUES" >
            <VAR name="TOTAL_SCORE" metric="total_hydrophobic_sasa" />
            <Var name="TOTAL_NUMBER_RESIDUES" metric="total_no_res_chA"/>
        </CalculatorMetric>
        <CalculatorMetric name="SAP_resnorm" equation="TOTAL_SCORE / TOTAL_NUMBER_RESIDUES" >
            <VAR name="TOTAL_SCORE" metric="sap_score" />
            <Var name="TOTAL_NUMBER_RESIDUES" metric="total_no_res_chA"/>
        </CalculatorMetric>
    </SIMPLE_METRICS>    
    <TASKOPERATIONS>
{1}
        <InitializeFromCommandline name="IFC" />
        <IncludeCurrent name="IC" />
        <LimitAromaChi2 name="aroChi" />
    </TASKOPERATIONS>

    <FILTERS>
        <ScoreType name="DesignBetaNov16CstFilter" scorefxn="DesignBetaNov16Cst" score_type="total_score" threshold="0" confidence="0" />
        <ScoreType name="BetaNov16" scorefxn="BetaNov16" score_type="total_score" threshold="0" confidence="0" />
        <RepeatParameter name="radius" param_type="radius" min="100" numb_repeats="{5}" confidence="0"/>
        <RepeatParameter name="rise" param_type="rise" numb_repeats="{5}" confidence="0"/>
        <RepeatParameter name="omega" param_type="omega" numb_repeats="{5}" confidence="0"/>
    </FILTERS>
    
    <MOVERS>
{2}

        <RandomMutation name="mutate" task_operations="not_ref_res,designable" scorefxn="BetaNov16"/>

        <AddCompositionConstraintMover name="load_aa_comp_cst" filename="%%aa_comp%%" selector="designable" />
{3}

        <FastDesign name="fast_design" scorefxn="DesignBetaNov16Cst" repeats="1" task_operations="designable,not_designable,IFC,IC,aroChi">
            <MoveMap name="fast_design_mm" chi="1" bb="0"/>
        </FastDesign>

        <FastRelax name="repack" scorefxn="DesignBetaNov16Cst" repeats="1" task_operations="not_designable,IFC,IC,aroChi" >
            <MoveMap name="repack_mm" chi="1" bb="0" jump="0"/>
        </FastRelax>

        <FastRelax name="full_repack" scorefxn="DesignBetaNov16Cst" repeats="1" task_operations="IFC,IC,aroChi" >
            <MoveMap name="full_repack_mm" chi="1" bb="0" jump="0"/>
        </FastRelax>

        <ParsedProtocol name="mc_round" mode="sequence" >
            <Add mover_name="mutate" />
            <Add mover_name="NCS" />
            <Add mover_name="repack" /> 
        </ParsedProtocol>

        <GenericSimulatedAnnealer name="monte_carlo_seq_design" mover_name="mc_round" periodic_mover="full_repack" eval_period="11" filter_name="DesignBetaNov16CstFilter" trials="30" sample_type="low" temperature="1" drift="1" recover_low="1" preapply="0" />
        
        <AddSapConstraintMover name="add_sap" />
    </MOVERS>
    
    <PROTOCOLS>
{4}
        <!-- <Add mover_name="add_sap" /> REPLACED WITH SIMPLE METRIC FROM FATIMA -->
        <Add filter_name="BetaNov16" />
        <Add filter_name="radius"/>
        <Add filter_name="rise"/>
        <Add filter_name="omega"/>
        <Add metrics="total_sasa" labels="total_sasa"/>
        <Add metrics="total_polar_sasa" labels="total_polar_sasa"/>
        <Add metrics="total_hydrophobic_sasa" labels="total_hydrophobic_sasa"/>
        <Add metrics="sap_score" labels="sap_score"/>
        <Add metrics="tot_Rscore" labels="tot_Rscore"/>
        <Add metrics="score_resnorm" labels="score_resnorm"/>
        <Add metrics="SASA_resnorm" labels="SASA_resnorm"/>
        <Add metrics="polarSASA_resnorm" labels="polarSASA_resnorm"/>
        <Add metrics="hydrophobicSASA_resnorm" labels="hydrophobicSASA_resnorm"/>
        <Add metrics="SAP_resnorm" labels="SAP_resnorm"/>
    </PROTOCOLS>
    
    <OUTPUT scorefxn="DesignBetaNov16Cst" />
</ROSETTASCRIPTS>'''