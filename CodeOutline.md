Outline of Codes for Chipod/Patches Analysis etc

- _FindPatches_EQ14_Raw.m_ - Identifes overturns in raw Chameleon profiles. Saves a structure 'new_patch_data'
- _run_eq14_for_PATCHES.m_ - Process Chameleon data only over patches (we use eps and chi from this further on). Note we don't use N2 or dTdz from this. This calls _average_data_PATCH_AP.m_, which does the chi and epsilon calculations over the specific patch regions only. 
- _average_data_PATCH_AP.m_ calls _calc_chi.m_ or _calc_chi_AP.m_ ; the difference is that _calc_chi_AP.m_ uses a f max of 7hz when integrating spectra. This is what I use in the chi-pod method, based on looking at spectra for this data.
- _Compute_N2_dTdz_ChamProfiles_V2.m_ - compute N^2 and dT/dz from cham profiles for patches, using a couple different methods. Produces a ‘patches’ structure. It also adds the chameleon chi and epsilon values at the patch locations, for both the patch and binned calculations.
- _MakePlotsGammaAnalysis.m_
- _ComputeChi_Chameleon_Eq14_PATCHES.m_ - This runs the chi-pod method to compute chi and epsilon for just the patches. 
- _CompareProfiles_bin_patch_cham.m_ - This compares epsilon from the chi-pod method to the chameleon epsilon (using shear probes) at patch locations.
