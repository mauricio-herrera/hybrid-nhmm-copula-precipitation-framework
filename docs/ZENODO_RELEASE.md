# GitHub → Zenodo: first software release — completed 2026-10-12

**Published:** GitHub release [v1.0.0](https://github.com/mauricio-herrera/hybrid-nhmm-copula-precipitation-framework/releases/tag/v1.0.0) was archived automatically by Zenodo on 12 October 2026. The Zenodo integration screen shows **DOI [10.5281/zenodo.23322615](https://doi.org/10.5281/zenodo.23322615)** with status **Published**. The Git tag points to commit `9ff10ebc5040aa31b9e7e8c0fed1964ba74d08d6`.

The steps below describe the procedure used, for future reference. **Do not repeat these steps for v1.0.0**; create a new version tag for any future archive.

1. Before publication, verify the input dataset's redistribution/attribution requirements. See [LICENSES.md](../LICENSES.md). Do not assume the MIT software license grants rights over third-party data.
2. Confirm that `CITATION.cff` has the correct author, version `1.0.0`, and **actual** release date. The currently prepared date is 2026-10-12; update it if you publish on a different date.
3. From GitHub choose **Releases → Draft a new release**.
4. Create tag **v1.0.0** from branch **main**. Use the title **Hybrid nHMM–Copula Precipitation Framework v1.0.0**. Do not mark it as a prerelease.
5. Publish the release. Do not separately create a manual Zenodo deposit: with the repository enabled, Zenodo should automatically archive the GitHub release.
6. Go to **Zenodo → GitHub**, check the archiving status, inspect the Zenodo record, and copy both DOI values (version and all-versions).
7. The confirmed DOI has been added to `CITATION.cff` and the README, with a DOI badge. This post-release commit does not retroactively change the archived v1.0.0 snapshot.

Suggested release notes:

> Initial public reproducibility release for the manuscript
> *Characterizing Climate-Driven Shifts in Chilean Rainfall Regimes with a
> Hybrid Hidden Markov–Copula Framework*. Includes processed precipitation
> inputs and an R script for state-conditional vine-copula simulations,
> diagnostics, and 2022 MJJA validation.
> The input contains previously computed rainfall regimes and DTW zones;
> preprocessing and upstream model fitting are not fully reproduced from raw
> atmospheric inputs in this archive. Consult LICENSES.md for CR2MET
> attribution and third-party data terms.

If Zenodo shows an error after the release, open the GitHub integration's **Errors** details and check `CITATION.cff` metadata before issuing another release.
