import numpy as np
import pandas as pd

INPUT_FOLDER = "/home/jkunz/master_uni_hd/internship_singapore/IF_Tuning"


def pMT_data_preparation():
    """
    - create from .csv file a dict of id - seq

    """
    seq_csv = pd.read_csv(f"{INPUT_FOLDER}/ptm_sequences.csv")
    seq_csv = seq_csv[["Unnamed: 0", "sequence", "activity ", "No.", "mutation"]]

    row_before = 0
    stage_before = "stage 1"
    stages = {}

    for row, val in enumerate(seq_csv["Unnamed: 0"]):
        if val is np.nan:
            continue
        if "stage" in val:
            stage = val.split(",")[0]
            before = (stage_before, row_before)
            stage_before = stage
            row_before = row
            stages[before[0]] = (before[1], row)

    stages[stage_before] = (row_before, len(seq_csv))

    for stage, idx in stages.items():
        seq_csv.loc[idx[0] : idx[1], "stage"] = stage

    seq_csv.drop(columns=["Unnamed: 0"], inplace=True)
    seq_csv["No."] = seq_csv["No."].apply(lambda x: x[1:-1])
    seq_csv.rename(
        columns={"No.": "pos", "mutation": "mut", "activity ": "activity"}, inplace=True
    )

    return seq_csv


def main():
    seq_csv = pMT_data_preparation()

    data = []
    ### EACH STAGE PAIRS ###
    for stage in seq_csv["stage"].unique():
        if stage == "stage 1":
            continue
        stage_df = seq_csv[seq_csv["stage"] == stage]
        pairs = stage_df.merge(stage_df, how="cross", suffixes=("_1", "_2"))

        mask = (pairs["pos_1"] != pairs["pos_2"]) & (
            pairs["activity_1"] > pairs["activity_2"]
        )
        result_df = pairs[mask]
        data.append(
            result_df[
                [
                    "sequence_1",
                    "sequence_2",
                    "stage_1",
                    "stage_2",
                    "mut_1",
                    "mut_2",
                    "activity_1",
                    "activity_2",
                ]
            ].values.tolist()
        )

    ### ACROSS STAGES PAIRS ###
    neg_seq_per_stage = {}
    pos_seq_per_stage = {}
    # first we construct the neg. and pos. sequences for each stage
    stages = seq_csv["stage"].unique()
    stages = stages[stages != "stage 1"]
    for stage in stages:
        if stage == "stage 1":
            continue
        stage_df = seq_csv[seq_csv["stage"] == stage]
        max_act = max(stage_df["activity"])
        pos_seq = []
        neg_seq = []
        for row, val in stage_df.iterrows():
            if val["activity"] < max_act:
                neg_seq.append((
                    val["sequence"],
                    val["pos"],
                    val["mut"],
                    val["activity"],
                ))
            if val["activity"] > 1:
                pos_seq.append((
                    val["sequence"],
                    val["pos"],
                    val["mut"],
                    val["activity"],
                ))
        neg_seq_per_stage[stage] = neg_seq
        pos_seq_per_stage[stage] = pos_seq

    # now we pair all pos. sequences of a stage with all neg. sequences of the lower stages,
    # due to the definition of the activity with the reference sequence,
    # we can be sure that the pos. sequences have higher activity than the neg. sequences
    for id, stage in enumerate(stages):
        if stage == "stage 1" or stage == "stage 2":
            continue
        lower_stages = stages[:id]
        for lower_stage in lower_stages:
            for pos_seq in pos_seq_per_stage[stage]:
                for neg_seq in neg_seq_per_stage[lower_stage]:
                    data.append((
                        pos_seq[0],
                        neg_seq[0],
                        stage,
                        lower_stage,
                        pos_seq[1],
                        neg_seq[1],
                        pos_seq[2],
                        neg_seq[2],
                        pos_seq[3],
                        neg_seq[3],
                    ))

    base_df = pd.DataFrame(
        data,
        columns=[
            "pos_seq",
            "neg_seq",
            "pos_stage",
            "neg_stage",
            "pos_pos",
            "neg_pos",
            "pos_mut",
            "neg_mut",
            "pos_act",
            "neg_act",
        ],
    )
    base_df["PDB_ID"] = "pMT_pdb.pdb"
