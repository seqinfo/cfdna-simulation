import pandas as pd
import numpy as np
from hmmlearn import hmm


# read count HMM
def rc_HMM(df, coverage, fetal_fraction, read_count_params):
    # start probability
    startprob = np.array([0.5, 0.5])

    # transition matrix
    ssr = 10  # stay-switch ratio
    eup_ssr = np.array([ssr, 1])
    tri_ssr = np.array([1, ssr])
    eup_prob = eup_ssr / eup_ssr.sum()
    tri_prob = tri_ssr / tri_ssr.sum()
    transmat = np.array([eup_prob, tri_prob])

    # read params
    rc = read_count_params[
        (read_count_params["coverage"] == coverage)
        & (read_count_params["fetal_fraction"] == fetal_fraction)
    ]

    # mean
    mean_eup = rc.loc[rc["state"] == 1, "mean"].values[0]
    mean_tri = rc.loc[rc["state"] == 2, "mean"].values[0]
    means = np.array([[mean_eup], [mean_tri]])

    # covariance
    var_eup = rc.loc[rc["state"] == 1, "variance"].values[0]
    var_tri = rc.loc[rc["state"] == 2, "variance"].values[0]
    covars = np.array([[[var_eup]], [[var_tri]]])

    # create model
    model = hmm.GaussianHMM(n_components=2, covariance_type="full")
    model.startprob_ = startprob
    model.transmat_ = transmat
    model.means_ = means
    model.covars_ = covars

    # run model
    values = df[["total_count"]].values
    lengths = df.groupby("sample").size().values
    states = model.predict(values, lengths) + 1

    return states


# read count HMM with fixed fetal fraction
def rc_HMM2(df, coverage, read_count_params):
    # fixed fetal fraction
    fetal_fraction = 0.1

    # start probability
    startprob = np.array([0.5, 0.5])

    # transition matrix
    ssr = 10  # stay-switch ratio
    eup_ssr = np.array([ssr, 1])
    tri_ssr = np.array([1, ssr])
    eup_prob = eup_ssr / eup_ssr.sum()
    tri_prob = tri_ssr / tri_ssr.sum()
    transmat = np.array([eup_prob, tri_prob])

    # read params
    # coverage = df.loc[df["chromosome"] != "chr21"]["count"].mean()  # add for experiment
    # coverage = (
    #     500 if coverage // 500 == 0 else coverage // 500 * 500
    # )  # add for experiment
    rc = read_count_params[
        (read_count_params["coverage"] == coverage)
        & (read_count_params["fetal_fraction"] == fetal_fraction)
    ]

    # mean
    mean_eup = rc.loc[rc["state"] == 1, "mean"].values[0]
    mean_tri = rc.loc[rc["state"] == 2, "mean"].values[0]
    means = np.array([[mean_eup], [mean_tri]])

    # covariance
    var_eup = rc.loc[rc["state"] == 1, "variance"].values[0]
    var_tri = rc.loc[rc["state"] == 2, "variance"].values[0]
    covars = np.array([[[var_eup]], [[var_tri]]])

    # create model
    model = hmm.GaussianHMM(n_components=2, covariance_type="full")
    model.startprob_ = startprob
    model.transmat_ = transmat
    model.means_ = means
    model.covars_ = covars

    # run model
    values = df[["total_count"]].values
    lengths = df.groupby("sample").size().values
    # values = df[["count"]].values  # add for experiment
    # lengths = df.groupby("sample").size().values  # add for experiment
    states = model.predict(values, lengths) + 1

    return states


# allelic ratio HMM
def ar_HMM(df, coverage, fetal_fraction, allelic_ratio_params):
    # start probability
    startprob = np.repeat(1 / 7, 7)

    # transition matrix
    ssr = 10  # stay-switch ratio
    eup_ssr = np.array([ssr, 1, ssr, 1, 1, ssr, 1])
    tri_ssr = np.array([1, ssr, 1, ssr, ssr, 1, ssr])
    eup_prob = eup_ssr / eup_ssr.sum()
    tri_prob = tri_ssr / tri_ssr.sum()
    transmat = np.array(
        [
            eup_prob,  # state 1 - euploidy
            tri_prob,  # state 2 - trisomy
            eup_prob,  # state 3 - euploidy
            tri_prob,  # state 4 - paternal trisomy
            tri_prob,  # state 5 - paternal trisomy
            eup_prob,  # state 6 - euploidy
            tri_prob,  # state 7 - trisomy
        ]
    )
    # read params
    ar = allelic_ratio_params[
        (allelic_ratio_params["coverage"] == coverage)
        & (allelic_ratio_params["fetal_fraction"] == fetal_fraction)
    ]

    # mean
    ar1_mean = ar.loc[ar["state"] == 1, "mean"].values[0]
    ar2_mean = ar.loc[ar["state"] == 2, "mean"].values[0]
    ar3_mean = ar.loc[ar["state"] == 3, "mean"].values[0]
    ar4_mean = ar.loc[ar["state"] == 4, "mean"].values[0]
    ar5_mean = ar.loc[ar["state"] == 5, "mean"].values[0]
    ar6_mean = ar.loc[ar["state"] == 6, "mean"].values[0]
    ar7_mean = ar.loc[ar["state"] == 7, "mean"].values[0]
    means = np.array(
        [
            [ar1_mean],  # state 1 - euploidy
            [ar2_mean],  # state 2 - trisomy
            [ar3_mean],  # state 3 - euploidy
            [ar4_mean],  # state 4 - paternal trisomy
            [ar5_mean],  # state 5 - paternal trisomy
            [ar6_mean],  # state 6 - euploidy
            [ar7_mean],  # state 7 - trisomy
        ]
    )

    # covariance
    ar1_var = ar.loc[ar["state"] == 1, "variance"].values[0]
    ar2_var = ar.loc[ar["state"] == 2, "variance"].values[0]
    ar3_var = ar.loc[ar["state"] == 3, "variance"].values[0]
    ar4_var = ar.loc[ar["state"] == 4, "variance"].values[0]
    ar5_var = ar.loc[ar["state"] == 5, "variance"].values[0]
    ar6_var = ar.loc[ar["state"] == 6, "variance"].values[0]
    ar7_var = ar.loc[ar["state"] == 7, "variance"].values[0]
    covars = np.array(
        [
            [[ar1_var]],  # state 1 - euploidy
            [[ar2_var]],  # state 2 - trisomy
            [[ar3_var]],  # state 3 - euploidy
            [[ar4_var]],  # state 4 - paternal trisomy
            [[ar5_var]],  # state 5 - paternal trisomy
            [[ar6_var]],  # state 6 - euploidy
            [[ar7_var]],  # state 7 - trisomy
        ]
    )

    # create model
    model = hmm.GaussianHMM(n_components=7, covariance_type="full")
    model.startprob_ = startprob
    model.transmat_ = transmat
    model.means_ = means
    model.covars_ = covars

    # run model
    df = df.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    values = df[["allelic_ratio"]].values
    lengths = df.groupby("sample").size().values
    states = model.predict(values, lengths) + 1

    return states


# read count and allelic ratio HMM
def rcar_HMM(df, coverage, fetal_fraction, read_count_params, allelic_ratio_params):
    # start probability
    startprob = np.repeat(1 / 7, 7)

    # transition matrix
    ssr = 10  # stay-switch ratio
    eup_ssr = np.array([ssr, 1, ssr, 1, 1, ssr, 1])
    tri_ssr = np.array([1, ssr, 1, ssr, ssr, 1, ssr])
    eup_prob = eup_ssr / eup_ssr.sum()
    tri_prob = tri_ssr / tri_ssr.sum()
    transmat = np.array(
        [
            eup_prob,  # state 1 - euploidy
            tri_prob,  # state 2 - trisomy
            eup_prob,  # state 3 - euploidy
            tri_prob,  # state 4 - paternal trisomy
            tri_prob,  # state 5 - paternal trisomy
            eup_prob,  # state 6 - euploidy
            tri_prob,  # state 7 - trisomy
        ]
    )

    # read params
    ar = allelic_ratio_params[
        (allelic_ratio_params["coverage"] == coverage)
        & (allelic_ratio_params["fetal_fraction"] == fetal_fraction)
    ]
    rc = read_count_params[
        (read_count_params["coverage"] == coverage)
        & (read_count_params["fetal_fraction"] == fetal_fraction)
    ]

    # mean
    ## allelic ratio
    ar1_mean = ar.loc[ar["state"] == 1, "mean"].values[0]
    ar2_mean = ar.loc[ar["state"] == 2, "mean"].values[0]
    ar3_mean = ar.loc[ar["state"] == 3, "mean"].values[0]
    ar4_mean = ar.loc[ar["state"] == 4, "mean"].values[0]
    ar5_mean = ar.loc[ar["state"] == 5, "mean"].values[0]
    ar6_mean = ar.loc[ar["state"] == 6, "mean"].values[0]
    ar7_mean = ar.loc[ar["state"] == 7, "mean"].values[0]

    ## read count
    rc_eup_mean = rc.loc[rc["state"] == 1, "mean"].values[0]
    rc_tri_mean = rc.loc[rc["state"] == 2, "mean"].values[0]

    means = np.array(
        [
            [ar1_mean, rc_eup_mean],  # state 1 - euploidy
            [ar2_mean, rc_tri_mean],  # state 2 - trisomy
            [ar3_mean, rc_eup_mean],  # state 3 - euploidy
            [ar4_mean, rc_tri_mean],  # state 4 - paternal trisomy
            [ar5_mean, rc_tri_mean],  # state 5 - paternal trisomy
            [ar6_mean, rc_eup_mean],  # state 6 - euploidy
            [ar7_mean, rc_tri_mean],  # state 7 - trisomy
        ]
    )

    # covariance
    ## allelic ratio
    ar1_var = ar.loc[ar["state"] == 1, "variance"].values[0]
    ar2_var = ar.loc[ar["state"] == 2, "variance"].values[0]
    ar3_var = ar.loc[ar["state"] == 3, "variance"].values[0]
    ar4_var = ar.loc[ar["state"] == 4, "variance"].values[0]
    ar5_var = ar.loc[ar["state"] == 5, "variance"].values[0]
    ar6_var = ar.loc[ar["state"] == 6, "variance"].values[0]
    ar7_var = ar.loc[ar["state"] == 7, "variance"].values[0]

    ## read count
    rc_eup_var = rc.loc[rc["state"] == 1, "variance"].values[0]
    rc_tri_var = rc.loc[rc["state"] == 2, "variance"].values[0]

    covars = np.array(
        [
            [[ar1_var, 0.0], [0.0, rc_eup_var]],  # state 1 - euploidy
            [[ar2_var, 0.0], [0.0, rc_tri_var]],  # state 2 - trisomy
            [[ar3_var, 0.0], [0.0, rc_eup_var]],  # state 3 - euploidy
            [[ar4_var, 0.0], [0.0, rc_tri_var]],  # state 4 - paternal trisomy
            [[ar5_var, 0.0], [0.0, rc_tri_var]],  # state 5 - paternal trisomy
            [[ar6_var, 0.0], [0.0, rc_eup_var]],  # state 6 - euploidy
            [[ar7_var, 0.0], [0.0, rc_tri_var]],  # state 7 - trisomy
        ]
    )

    # create model
    model = hmm.GaussianHMM(n_components=7, covariance_type="full")
    model.startprob_ = startprob
    model.transmat_ = transmat
    model.means_ = means
    model.covars_ = covars

    # run model
    df = df.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    values = df[["allelic_ratio", "total_count"]].values
    lengths = df.groupby("sample").size().values
    states = model.predict(values, lengths) + 1

    return states
