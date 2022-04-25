from bnspopkne.kne import saee_bns_emgw_with_viewing_angle as saeev


test_inst = saeev(
    m1=2.0,
    m2=1.75,
    EOS="sfho",
    EOS_path="~/Documents/Project1_kNe/kne_modeling/eos_data/",
    kappa_grid_path="~/Documents/Project1_kNe/kne_modeling/korobkin_heating_rates/outputs/thresholded_uncertainties_20_opacity_df_020221.csv",
    hyperparam_file="../data/paper_matern52_hyperparameters.npy"
)

# List out the parameters of the BNS mergers kilonova and binary inspiral
for i in range(12):
    print(getattr(test_inst, f"param{i+1}"))
